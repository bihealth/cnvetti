// Implementation of the "segment" command.

include!(concat!(env!("OUT_DIR"), "/version.rs"));

mod math;
mod options;
pub mod seg_haar; // XXX should not be pub
pub mod seg_utils; // XXX should not be pub

use std::env;
use std::str;

use rust_htslib::bcf::{self, Read as BcfRead};

use shlex;

use slog::Logger;

pub use self::options::*;
use cli::shared;

use separator::Separatable;

/// Generate list of all contigs from BCF header.
fn build_chroms(header: &bcf::header::HeaderView) -> Vec<(String, u32)> {
    let mut result = Vec::new();

    for ref record in header.header_records() {
        if let bcf::HeaderRecord::Contig {
            key,
            key_value_pairs,
        } = record
        {
            assert_eq!(key, "contig");
            let mut name: Option<String> = None;
            let mut length: Option<u32> = None;
            for &(ref key, ref value) in key_value_pairs {
                if key == "ID" {
                    name = Some(value.clone());
                } else if key == "length" {
                    length = Some(value.parse::<u32>().unwrap());
                }
            }
            if let (Some(ref name), Some(length)) = (&name, length) {
                result.push((name.clone(), length));
            } else {
                panic!(
                    "Could not parse both name/length from {:?}/{:?}",
                    name, length
                );
            }
        }
    }

    result
}

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

    // First pass, collect normalized coverage and filtered-out mask.
    debug!(logger, "First pass over region (collect coverage)...");
    reader
        .fetch(rid, *start, *end)
        .expect("Could not fetch region from BCF file");
    let mut record = reader.empty_record();
    while reader.read(&mut record).is_ok() {
    }

    debug!(logger, "Performing segmentation...");

    // Second pass, write segmentation
    debug!(logger, "Second pass over region (write segmentation)...");
    reader
        .fetch(rid, *start, *end)
        .expect("Could not fetch region from BCF file");
    while reader.read(&mut record).is_ok() {
        writer.translate(&mut record);
        writer.write(&record).expect("Writing the record failed!");
    }
}

/// Perform the actual processing.
fn process(
    reader: &mut bcf::IndexedReader,
    writer: &mut bcf::Writer,
    logger: &mut Logger,
    options: &Options,
) {
    let contigs = build_chroms(reader.header());
    for (chrom, length) in contigs.iter() {
        process_region(
            reader,
            writer,
            logger,
            options,
            &(chrom.to_string(), 0, *length),
        );
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
                "##FORMAT=<ID=SRC,Number=1,Type=Float,Description=\"Segmented normalized read \
                 count\">",
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
        process(&mut reader, &mut writer, logger, &options);
    }

    // Build index on the output file.
    shared::build_index(logger, &options.output);

    info!(logger, "Done. Have a nice day!");

    Ok(())
}
