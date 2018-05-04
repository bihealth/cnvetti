include!(concat!(env!("OUT_DIR"), "/version.rs"));

mod options;

use std::env;

use shlex;

use slog::Logger;

pub use self::options::*;
use cli::shared::{self, process};
use separator::Separatable;

use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::{self, Read as BcfRead};

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
            // TODO: remove header fields...
            let lines = vec![
                "##ALT=<ID=DEL,Description=\"Deletion relative to the reference\">",
                "##ALT=<ID=DUP,Description=\"Region of elevated copy number relative to the \
                 reference\">",
                "##FORMAT<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
                "##FORMAT<ID=GP,Number=1,Type=Integer,Description=\"Genotype quality\">",
                "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for \
                 imprecise events\">",
                "##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=\"Copy number genotype \
                 quality for imprecise events\">",
                "##FORMAT=<ID=CNL,Number=G,Type=Float,Description=\"Copy number genotype \
                 likelihood for imprecise events\">",
            ];
            for line in lines {
                header.push_record(line.as_bytes());
            }
            header.push_record(format!("##cnvetti_callVersion={}", VERSION).as_bytes());
            header.push_record(
                format!(
                    "##cnvetti_callCommand={}",
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
