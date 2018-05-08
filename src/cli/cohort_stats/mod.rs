include!(concat!(env!("OUT_DIR"), "/version.rs"));

pub mod options;
pub use self::options::*;

use shlex;

use std::env;

use slog::Logger;

use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::{self, Read as BcfRead};

use cli::shared;
use cli::shared::stats::Stats;

/// Main entry point for the "normalize" command.
pub fn call(logger: &mut Logger, options: &Options) -> Result<(), String> {
    // let mut app = CoverageApp::new(logger, options);
    info!(logger, "Running cnvetti cohort-stats");
    info!(logger, "Configuration: {:?}", &options);

    let mut iqrs = Vec::new();

    // First pass, gather IQRs, collect values.
    {
        let mut reader =
            bcf::Reader::from_path(options.input.clone()).expect("Could not open input BCF file");
        if options.io_threads > 0 {
            reader
                .set_threads(options.io_threads as usize)
                .expect("Could not set I/O thread count");
        }

        let mut record = reader.empty_record();
        while reader.read(&mut record).is_ok() {
            // Compute and write inter quartile range of FORMAT/NCOV.
            let ncovs: Vec<f64> = record
                .format(b"NCOV")
                .float()
                .map_err(|e| format!("Could not access FORMAT/NCOV: {}", e))?
                .iter()
                .map(|a| a[0] as f64)
                .collect();
            let iqr = ncovs.percentile(75.0) - ncovs.percentile(25.0);
            if iqr.is_finite() {
                iqrs.push(iqr as f64);
            }
        }
    }

    // Compute threshold.
    let iqr_thresh = iqrs.percentile(options.percentile_threshold);
    info!(
        logger,
        "Setting IQR threshold at {}% = {} (log2: {}/{})",
        options.percentile_threshold,
        iqr_thresh,
        iqr_thresh.log2(),
        (1.0 + iqr_thresh).log2()
    );

    // Second pass, write out filter tag.
    {
        let mut reader =
            bcf::Reader::from_path(options.input.clone()).expect("Could not open input BCF file");
        if options.io_threads > 0 {
            reader
                .set_threads(options.io_threads as usize)
                .expect("Could not set I/O thread count");
        }

        let mut writer = {
            // Construct extended header.
            let mut header = bcf::Header::from_template(reader.header());
            let lines = vec![
                "##INFO=<ID=NCOV_IQR,Number=1,Type=Float,Description=\"Normalized coverage \
                 cohort IQR\">",
                "##INFO=<ID=NCOV_IQR_PASS,Number=1,Type=Float,Description=\"\
                 Whether or not IQR filter passed\">",
            ];
            for line in lines {
                header.push_record(line.as_bytes());
            }
            header.push_record(format!("##cnvetti_cohortStatsVersion={}", VERSION).as_bytes());
            header.push_record(
                format!(
                    "##cnvetti_cohortStatsCommand={}",
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

        let mut record = reader.empty_record();
        while reader.read(&mut record).is_ok() {
            // Translate to output header.
            writer.translate(&mut record);

            // Compute and write inter quartile range of FORMAT/NCOV.
            let ncovs: Vec<f64> = record
                .format(b"NCOV")
                .float()
                .map_err(|e| format!("Could not access FORMAT/NCOV: {}", e))?
                .iter()
                .map(|a| a[0] as f64)
                .collect();
            let iqr = ncovs.percentile(75.0) - ncovs.percentile(25.0);
            let (iqr, pass) = if iqr.is_finite() {
                (iqr as f32, iqr < iqr_thresh)
            } else {
                (f32::missing(), false)
            };
            record
                .push_info_float(b"NCOV_IQR", &[iqr])
                .map_err(|e| format!("Could not write INFO/NCOV_IQR: {}", e))?;
            record
                .push_info_integer(b"NCOV_IQR_PASS", &[pass as i32])
                .map_err(|e| format!("Could not write INFO/NCOV_IQR_PASS: {}", e))?;

            // Write out record again.
            writer.write(&record).expect("Could not write BCF record!");
        }
    }

    // Build index on the output file.
    shared::build_index(logger, &options.output);

    info!(logger, "Done. Have a nice day!");

    Ok(())
}
