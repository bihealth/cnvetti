// Implementation of the "normalize" command.

include!(concat!(env!("OUT_DIR"), "/version.rs"));

mod options;

use std::env;
use std::fs::File;
use std::io::Read;
use std::str;

use rust_htslib::bcf;

use serde_json;

use shlex;

use slog::Logger;

use cli::coverage::summary::SummarisedMetric;
use cli::shared;
pub use self::options::*;


// TODO: check input file.
// TODO: use index-based readers, is nicer for progress display...
// TODO: could normalize already in coverage step by writing raw coverage to temporary file


/// Perform the actual processing.
fn process(
    reader: &mut bcf::Reader,
    writer: &mut bcf::Writer,
    logger: &mut Logger,
    summaries: Vec<SummarisedMetric>,
    options: &Options,
) {
    let mut record = reader.empty_record();
    let mut prev_rid = Option::None;

    let header = reader.header().clone();

    while reader.read(&mut record).is_ok() {
        if prev_rid.is_some() && prev_rid.unwrap() != record.rid() {
            info!(
                logger,
                "Starting on contig {}",
                str::from_utf8(header.rid2name(record.rid().unwrap())).unwrap();
            );
        }

        // Translate to output header.
        writer.translate(&mut record);

        // Get INFO/GCWINDOWS and maybe add INFO/FEW_GCWINDOWS.
        {
            let gc_windows = record
                .info(b"GCWINDOWS")
                .integer()
                .expect("Could not read INFO/GCWINDOWS")
                .unwrap_or(&[0])
                [0];
            if gc_windows < options.min_gc_window_count {
                record.push_filter(b"FEW_GCWINDOWS").expect(
                    "Could not write INFO/FEW_GCWINDOWS",
                );
            }
        }

        // Get INFO/GC.
        let gc = {
            record
                .info(b"GC")
                .float()
                .expect("Could not read INFO/GC")
                .expect("INFO/GC was empty")
                [0]
        };
        // TODO: gc_step should come frome file header...
        let gc_step: f32 = 0.02;
        let bucket = ((gc as f32 / gc_step).floor() * 100.0 * gc_step) as i32;

        let sample_id = 0;

        // Normalize FORMAT/RC and FORMAT/{COV,SD}.
        match options.count_kind {
            CountKind::Coverage => {
                panic!("Not implemented yet!");
            }
            CountKind::Alignments => {
                let summary = &summaries[sample_id].summaries.get(&bucket).unwrap();
                let median = summary.summary5[2] as f32;
                let nrcs: Vec<f32> = record
                    .format(b"RC")
                    .integer()
                    .expect("Cannot access FORMAT/RC")
                    .iter()
                    .map(|slice| if median < 0.001 {
                        // guard against NaN
                        0.0
                    } else {
                        slice[0] as f32 / median
                    })
                    .collect();
                record.push_format_float(b"NRC", nrcs.as_slice()).expect(
                    "Could not write FORMAT/NRC",
                );
            }
        }

        // Finally, write out record again.
        writer.write(&record).expect("Could not write BCF record!");

        prev_rid = Some(record.rid());
    }
}


/// Main entry point for the "normalize" command.
pub fn call(logger: &mut Logger, options: &Options) -> Result<(), String> {
    let options = options.clone();

    // let mut app = CoverageApp::new(logger, options);
    info!(logger, "Running cnvetti normalize");
    info!(logger, "Configuration: {:?}", &options);

    info!(logger, "Opening input and output files...");
    debug!(logger, "Loading statistics");
    let stats: Vec<SummarisedMetric> = {
        let stats_filename = options.input.clone() + ".stats.txt";
        let mut file = File::open(stats_filename).expect("Could not open stats file for reading");
        let mut json = String::new();
        file.read_to_string(&mut json).expect(
            "Could not read from stats file",
        );
        serde_json::from_str(&json).expect("Could not deserialize stats from JSON")
    };

    debug!(logger, "Opening input file");
    let mut reader =
        bcf::Reader::from_path(options.input.clone()).expect("Could not open input BCF file");
    reader.set_threads(options.io_threads as usize).expect(
        "Could not set I/O thread count",
    );

    // Open the output file in its own block so we can close before creating the index.
    {
        let mut writer = {
            // Construct extended header.
            let mut header = bcf::Header::with_template(reader.header());
            let lines = vec![
                "##FILTER=<ID=FEW_GCWINDOWS,Description=\"Masked because of few windows with \
                    this GC content\">",
                "##FORMAT=<ID=NRC,Number=1,Type=Float,Description=\"Normalized number of \
                    aligning reads\">",
            ];
            for line in lines {
                header.push_record(line.as_bytes());
            }
            header.push_record(format!("##cnvetti_normalizeVersion={}", VERSION).as_bytes());
            header.push_record(
                format!(
                    "##cnvetti_normalizeCommand={}",
                    env::args()
                        .map(|s| shlex::quote(&s).to_string())
                        .collect::<Vec<String>>()
                        .join(" ")
                ).as_bytes(),
            );

            let uncompressed = !options.output.ends_with(".bcf") &&
                !options.output.ends_with(".vcf.gz");
            let vcf = options.output.ends_with(".vcf") || options.output.ends_with(".vcf.gz");

            bcf::Writer::from_path(options.output.clone(), &header, uncompressed, vcf)
                .expect("Could not open output BCF file")
        };
        writer.set_threads(options.io_threads as usize).expect(
            "Could not set I/O thread count",
        );

        info!(logger, "Processing...");
        process(&mut reader, &mut writer, logger, stats, &options);
    }

    // Build index on the output file.
    shared::build_index(logger, &options.output);

    info!(logger, "Done. Have a nice day!");

    Ok(())
}
