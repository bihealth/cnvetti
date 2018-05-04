// Implementation of the "segment" command.

include!(concat!(env!("OUT_DIR"), "/version.rs"));

mod options;

use std::env;
use std::str;

use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::{self, Read as BcfRead};

use shlex;

use slog::Logger;

pub use self::options::*;
use cli::shared;

use separator::Separatable;

use rust_segment::{adjust_breaks, reject_nonaberrant, seg_haar};

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
    let n_samples = reader.header().sample_count() as usize;
    let mut coverage: Vec<Vec<f64>> = Vec::new();
    coverage.resize(n_samples, Vec::new());

    // Key for the metric to segment.
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

        for (i, val) in record
            .format(key)
            .float()
            .expect("Record does not have field")
            .iter()
            .enumerate()
        {
            coverage[i].push(((**val)[0].max(0.0) as f64 + PSEUDO_EPSILON).log2());
        }
    }

    debug!(logger, "Performing segmentation...");
    const MAX_REFINEMENTS: usize = 10; // TODO: get from configuration
    let coverage: Vec<Vec<f64>> = coverage
        .iter()
        .enumerate()
        .map(|(i, ref vals)| {
            debug!(
                logger,
                "Segmenting for {}...",
                str::from_utf8(&reader.header().id_to_sample(bcf::header::Id(i as u32))).unwrap()
            );

            // Perform segmentation, yielding breakpoints.
            let mut haar_res = seg_haar(
                &vals,
                None,
                None,
                &[0..(vals.len())],
                options.haar_seg_breaks_fdr_q,
                options.haar_seg_l_min as u32,
                options.haar_seg_l_max as u32,
            );
            debug!(
                logger,
                "Raw segments: {}",
                haar_res.segments.len().separated_string()
            );
            // Adjust breakpoints.
            for i in 0..MAX_REFINEMENTS {
                let (num_refined, inner_haar_res) = adjust_breaks(&haar_res, &vals);
                debug!(
                    logger,
                    "{} changed after {}-th break adjustments",
                    num_refined,
                    i + 1
                );
                if num_refined > 0 {
                    haar_res = inner_haar_res;
                } else {
                    break; // nothing more to do
                }
            }
            debug!(
                logger,
                "Segments after adjustment: {}",
                haar_res.segments.len().separated_string()
            );
            // Reject non-aberrant segments.
            const ABERRANT_FACTOR: f64 = 3.0;
            haar_res = reject_nonaberrant(&haar_res, &vals, ABERRANT_FACTOR);

            debug!(
                logger,
                "Final segments: {}",
                haar_res.segments.len().separated_string()
            );

            haar_res.seg_values
        })
        .collect();

    // Second pass, write segmentation
    debug!(logger, "Second pass over region (write segmentation)...");
    let mut idx = 0; // current index into coverage[i]
    let mut prev_vals: Option<Vec<f32>> = None; // [log2(val)]
    reader
        .fetch(rid, *start, *end)
        .expect("Could not fetch region from BCF file");
    while reader.read(&mut record).is_ok() {
        writer.translate(&mut record);
        let vals = if skip_record(&mut record, options) {
            record.push_filter(writer.header().name_to_id(b"SEG_SKIPPED").unwrap());
            if let Some(ref prev_vals) = &prev_vals {
                prev_vals.clone()
            } else {
                vec![0_f32; n_samples]
            }
        } else {
            let vals = (0..n_samples)
                .map(|i| {
                    if coverage[i][idx] > 0.0 && coverage[i][idx] < PSEUDO_EPSILON {
                        0.0
                    } else {
                        (coverage[i][idx] - PSEUDO_EPSILON) as f32
                    }
                })
                .collect::<Vec<f32>>();
            prev_vals = Some(vals.clone());
            idx += 1;
            vals
        };
        record
            .push_format_float(b"SCOV", &[2_f32.powf(vals[0])])
            .expect("Could not write FORMAT/SCOV");
        record
            .push_format_float(b"SCOV2", &vals)
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
            let mut header = bcf::Header::with_template(reader.header());
            let lines = vec![
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
