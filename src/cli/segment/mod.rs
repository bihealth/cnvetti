// Implementation of the "segment" command.

include!(concat!(env!("OUT_DIR"), "/version.rs"));

pub mod options;
pub use self::options::*;

use std::env;
use std::str;

use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::{self, Read as BcfRead};

use shlex;

use slog::Logger;

use cli::shared;

use separator::Separatable;

use rust_segment::{reject_nonaberrant_pvalue, seg_haar};

use cli::shared::process;

/// Epsilon to add normalized metric to prevent -nan through log2()
const PSEUDO_EPSILON: f64 = 1e-12;

/// Process one region.
fn process_region(
    reader: &mut bcf::IndexedReader,
    writer: &mut bcf::Writer,
    logger: &mut Logger,
    options: &Options,
    (chrom, start, end): &(String, u32, u32),
) {
    info!(
        logger,
        "Processing {}:{}-{}",
        chrom,
        (start + 1).separated_string(),
        end.separated_string()
    );

    // Resolve contig name to contig index.
    let rid = reader
        .header()
        .name2rid(chrom.as_bytes())
        .expect("Could not resolve contig name to index");

    // Vectors for the coverage markers.
    let mut ncov: Vec<f64> = Vec::new();
    let mut ncov_log2: Vec<f64> = Vec::new();

    // Key for the metric to segment (do log2 after conversion to f64 here).
    let key = b"NCOV";

    // Whether or not to skip a record.
    fn skip_record(record: &mut bcf::Record, options: &Options) -> bool {
        // TODO: skip on mapability, cohort IQR
        // TODO: record should not be mut, see rust-bio/rust-htslib#78

        let is_missing = record.format(b"COV").float().unwrap()[0][0].is_missing();
        let gap = record
            .info(b"GAP")
            .integer()
            .unwrap_or(Some(&[0_i32]))
            .expect("INFO/GAP empty")[0] != 0;
        let gc: f32 = record
            .info(b"GC")
            .float()
            .expect("no INFO/GC")
            .expect("INFO/GC empty")[0];
        // let few_gc_windows = record.has_filter(b"FEW_GCWINDOWS");

        (gc < options.min_gc) || (gc > options.max_gc) || gap || is_missing // || few_gc_windows
    }

    // First pass, collect normalized coverage and filtered-out mask.
    debug!(logger, "First pass over region (collect coverage)...");
    reader
        .fetch(rid, *start, *end)
        .expect("Could not fetch region from BCF file");
    let mut record = reader.empty_record();
    while reader.read(&mut record).is_ok() {
        writer.translate(&mut record);

        if skip_record(&mut record, options) {
            continue;
        }

        let val = record
            .format(key)
            .float()
            .expect("Record does not have field")[0][0] as f64;
        ncov.push(val);
        ncov_log2.push((val.max(0.0) as f64 + PSEUDO_EPSILON).log2());
    }

    debug!(logger, "Performing segmentation...");

    // Perform segmentation, yielding breakpoints.
    let mut segmentation = seg_haar(
        &ncov_log2,
        None,
        None,
        &[0..(ncov_log2.len())],
        options.haar_seg_breaks_fdr_q,
        options.haar_seg_l_min as u32,
        options.haar_seg_l_max as u32,
    );
    debug!(
        logger,
        "Raw segments: {}",
        segmentation.segments.len().separated_string()
    );
    // Record P-values, to later write out.
    let mut p_values = Vec::new();
    for seg in &segmentation.segments {
        let new_len = p_values.len() + seg.range.len();
        let p_val = if seg.range.len() > 2 {
            seg.p_value_significant_student(options.significant_p_val_thresh) as f32
        } else {
            1.0
        };
        p_values.resize(new_len, p_val);
    }
    // Reject segments not passing p value.
    segmentation =
        reject_nonaberrant_pvalue(&segmentation, &ncov, options.significant_p_val_thresh);
    debug!(
        logger,
        "Segments after selecting aberrant (p): {}",
        segmentation.segments.len().separated_string()
    );

    // Second pass, write segmentation
    debug!(logger, "Second pass over region (write segmentation)...");
    let mut idx = 0; // current index into ncov
    let mut prev_val: Option<(f32, f32)> = None; // (val, log2(val))
    reader
        .fetch(rid, *start, *end)
        .expect("Could not fetch region from BCF file");
    while reader.read(&mut record).is_ok() {
        writer.translate(&mut record);

        let (val, val_log2) = if skip_record(&mut record, options) {
            record.push_filter(writer.header().name_to_id(b"SEG_SKIPPED").unwrap());
            if let Some(prev_val) = prev_val {
                prev_val
            } else {
                (1.0, 0.0)
            }
        } else {
            let val = if segmentation.values[idx] > 0.0 && segmentation.values[idx] < PSEUDO_EPSILON
            {
                (1.0, 0.0)
            } else {
                (
                    (segmentation.values[idx] - PSEUDO_EPSILON) as f32,
                    segmentation.values_log2[idx] as f32,
                )
            };
            prev_val = Some(val);
            idx += 1;
            val
        };
        record
            .push_format_float(b"PVAL", &[p_values[idx]])
            .expect("Could not write FORMAT/PVAL");
        record
            .push_format_float(b"SCOV", &[val])
            .expect("Could not write FORMAT/SCOV");
        record
            .push_format_float(b"SCOV2", &[val_log2])
            .expect("Could not write FORMAT/SCOV2");
        writer.write(&record).expect("Writing the record failed!");
    }
}

/// Main entry point for the "segment" command.
pub fn call(logger: &mut Logger, options: &Options) -> Result<(), String> {
    let options = options.clone();

    // let mut app = CoverageApp::new(logger, options);
    info!(logger, "Running cnvetti segment");
    info!(logger, "Configuration: {:?}", &options);

    debug!(logger, "Opening input file");
    let mut reader = bcf::IndexedReader::from_path(options.input.clone())
        .expect("Could not open input BCF file");
    if options.io_threads > 0 {
        reader
            .set_threads(options.io_threads as usize)
            .expect("Could not set I/O thread count");
    }

    // Open the output file in its own block so we can close before creating the index.
    {
        let mut writer = {
            // Construct extended header.
            let mut header = bcf::Header::from_template(reader.header());
            let lines = vec![
                "##FORMAT=<ID=PVAL,Number=1,Type=Float,Description=\"P value of original segment\">",
                "##FORMAT=<ID=SCOV,Number=1,Type=Float,Description=\"Segmented coverage\">",
                "##FORMAT=<ID=SCOV2,Number=1,Type=Float,Description=\"Segmented coverage \
                 (log2-scaled)\">",
                "##FILTER=<ID=SEG_SKIPPED,Description=\"Window skipped in segmentation\">",
            ];
            for line in lines {
                header.push_record(line.as_bytes());
            }
            header.push_record(format!("##cnvetti_segmentVersion={}", VERSION).as_bytes());
            header.push_record(
                format!(
                    "##cnvetti_segmentCommand={}",
                    env::args()
                        .map(|s| shlex::quote(&s).to_string())
                        .collect::<Vec<String>>()
                        .join(" ")
                ).as_bytes(),
            );

            let uncompressed =
                !options.output.ends_with(".bcf") && !options.output.ends_with(".vcf.gz");
            let vcf = options.output.ends_with(".vcf") || options.output.ends_with(".vcf.gz");

            bcf::Writer::from_path(options.output.clone(), &header, uncompressed, vcf)
                .expect("Could not open output BCF file")
        };
        if options.io_threads > 0 {
            writer
                .set_threads(options.io_threads as usize)
                .expect("Could not set I/O thread count");
        }

        info!(logger, "Processing...");
        process(&mut reader, &mut writer, logger, &options, &process_region);
    }

    // Build index on the output file.
    shared::build_index(logger, &options.output);

    info!(logger, "Done. Have a nice day!");

    Ok(())
}
