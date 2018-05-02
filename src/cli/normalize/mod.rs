// Implementation of the "normalize" command.

// include!(concat!(env!("OUT_DIR"), "/version.rs"));

mod binned_gc;
mod options;
mod shared;

// use std::env;
// use std::io;
// use std::str;

// use rust_htslib::bcf::record::Numeric;
// use rust_htslib::bcf::{self, Read as BcfRead};

// use regex::Regex;

// use shlex;

use slog::Logger;

use self::binned_gc::normalize_binned_gc;

pub use self::options::*;
// use cli::shared::math;
use cli::shared::build_index;

// TODO: check input file.
// TODO: use index-based readers, is nicer for progress display...
// TODO: could normalize already in coverage step by writing raw coverage to temporary file

/// Main entry point for Lowess-based normalization.
pub fn call_lowess(logger: &mut Logger, _options: &Options) -> Result<(), String> {
    info!(logger, "Normalizing by Lowess via R...");

    // // Read values into vector.
    // let mut gcs: Vec<f64> = Vec::new();
    // let mut maps: Vec<f64> = Vec::new();
    // let mut covs_raw: Vec<f64> = Vec::new();
    // {
    //     debug!(logger, "Opening input BCF file collecting Lowess");
    //     let mut reader =
    //         bcf::Reader::from_path(options.input.clone()).expect("Could not open input BCF file");
    //     if options.io_threads > 0 {
    //         reader
    //             .set_threads(options.io_threads as usize)
    //             .expect("Could not set I/O thread count");
    //     }

    //     let re = Regex::new(&options.contig_regex).unwrap();
    //     let mut record = reader.empty_record();
    //     while reader.read(&mut record).is_ok() {
    //         // Get "is gap" flag.
    //         let is_gap = {
    //             record
    //                 .info(b"GAP")
    //                 .integer()
    //                 .expect("Could not read INFO/GAP")
    //                 .expect("INFO/GAP was empty")[0] != 0
    //         };
    //         // Get "is blacklisted" flag.
    //         let is_blacklist = {
    //             record
    //                 .info(b"BLACKLIST")
    //                 .integer()
    //                 .unwrap_or(Some(&[0]))
    //                 .expect("INFO/BLACKLIST was empty")[0] != 0
    //         };

    //         // Check whether window should be used in Lowess.
    //         let chrom = str::from_utf8(reader.header().rid2name(record.rid().unwrap()))
    //             .unwrap()
    //             .to_string();
    //         let use_chrom = re.is_match(&chrom);

    //         if use_chrom && !is_gap && !is_blacklist {
    //             // Get INFO/GC.
    //             let gc_content = {
    //                 record
    //                     .info(b"GC")
    //                     .float()
    //                     .expect("Could not read INFO/GC")
    //                     .expect("INFO/GC was empty")[0]
    //             };
    //             // TODO: add switch whether mapability is desired, also for when writing out script above.
    //             // Get INFO/MAPABILITY
    //             let mapability = {
    //                 match record.info(b"MAPABILITY").float() {
    //                     Ok(ref mapability) => match mapability {
    //                         Some(ref mapability) => {
    //                             if mapability.len() > 0 {
    //                                 mapability[0]
    //                             } else {
    //                                 0_f32
    //                             }
    //                         }
    //                         None => 0_f32,
    //                     },
    //                     Err(_) => 0_f32,
    //                 }
    //             };

    //             let cov = record
    //                 .format(b"COV")
    //                 .float()
    //                 .expect("FORMAT/COV not found!");
    //             gcs.push(gc_content as f64);
    //             maps.push(mapability as f64);
    //             covs_raw.push(cov[0][0] as f64);
    //         }
    //     }
    // }

    // // TODO: write out points before and after smoothing!
    // let delta = covs_raw.iter().cloned().fold(0. / 0., f64::max)
    //     - covs_raw.iter().cloned().fold(0. / 0., f64::min);
    // let fit_gc = shared::lowess(&gcs, &covs_raw, 0.6, 4, 0.01 * delta).fit;
    // // let covs_gc: Vec<f64> = covs_raw.iter().enumerate().map(|(i, x)| *x / fit_gc[i]).collect();
    // let covs_gc = fit_gc.clone();

    // let delta = fit_gc.iter().cloned().fold(0. / 0., f64::max)
    //     - fit_gc.iter().cloned().fold(0. / 0., f64::min);
    // let fit_map = shared::lowess(&maps, &covs_raw, 0.6, 4, 0.01 * delta).fit;
    // // let covs_map: Vec<f64> = covs_gc.iter().enumerate().map(|(i, x)| *x / fit_map[i]).collect();
    // let covs_map = fit_map.clone();

    // // Block for BCF file reader and writer.
    // {
    //     debug!(
    //         logger,
    //         "Opening input BCF file for generating final BCF file"
    //     );
    //     let mut reader =
    //         bcf::Reader::from_path(options.input.clone()).expect("Could not open input BCF file");
    //     if options.io_threads > 0 {
    //         reader
    //             .set_threads(options.io_threads as usize)
    //             .expect("Could not set I/O thread count");
    //     }

    //     let mut writer = {
    //         // Construct extended header.
    //         let mut header = bcf::Header::with_template(reader.header());
    //         let lines = vec![
    //             "##FORMAT=<ID=NCOV_GC,Number=1,Type=Float,Description=\"Normalized coverage \
    //              (GC only)\">",
    //             "##FORMAT=<ID=NCOV2_GC,Number=1,Type=Float,Description=\"Normalized coverage \
    //              (GC only; log2-scaled)\">",
    //             "##FORMAT=<ID=NCOV,Number=1,Type=Float,Description=\"Normalized coverage\">",
    //             "##FORMAT=<ID=NCOV2,Number=1,Type=Float,Description=\"Normalized coverage \
    //              (log2-scaled)\">",
    //         ];
    //         for line in lines {
    //             header.push_record(line.as_bytes());
    //         }
    //         header.push_record(format!("##cnvetti_normalizeVersion={}", VERSION).as_bytes());
    //         header.push_record(
    //             format!(
    //                 "##cnvetti_normalizeCommand={}",
    //                 env::args()
    //                     .map(|s| shlex::quote(&s).to_string())
    //                     .collect::<Vec<String>>()
    //                     .join(" ")
    //             ).as_bytes(),
    //         );

    //         let uncompressed =
    //             !options.output.ends_with(".bcf") && !options.output.ends_with(".vcf.gz");
    //         let vcf = options.output.ends_with(".vcf") || options.output.ends_with(".vcf.gz");

    //         bcf::Writer::from_path(options.output.clone(), &header, uncompressed, vcf)
    //             .expect("Could not open output BCF file")
    //     };
    //     if options.io_threads > 0 {
    //         writer
    //             .set_threads(options.io_threads as usize)
    //             .expect("Could not set I/O thread count");
    //     }

    //     let re = Regex::new(&options.contig_regex).unwrap();
    //     let mut record = reader.empty_record();
    //     let mut prev_rid = Option::None;

    //     let header = reader.header().clone();

    //     let mut i = 0;
    //     while reader.read(&mut record).is_ok() {
    //         if prev_rid.is_some() && prev_rid.unwrap() != record.rid() {
    //             info!(
    //                 logger,
    //                 "Starting on contig {}",
    //                 str::from_utf8(header.rid2name(record.rid().unwrap())).unwrap()
    //             );
    //         }

    //         // Translate to output header.
    //         writer.translate(&mut record);

    //         // Get "is gap" flag.
    //         let is_gap = {
    //             record
    //                 .info(b"GAP")
    //                 .integer()
    //                 .expect("Could not read INFO/GAP")
    //                 .expect("INFO/GAP was empty")[0] != 0
    //         };
    //         // Get "is blacklisted" flag.
    //         let is_blacklist = {
    //             record
    //                 .info(b"BLACKLIST")
    //                 .integer()
    //                 .unwrap_or(Some(&[0]))
    //                 .expect("INFO/GAP was empty")[0] != 0
    //         };

    //         // Check whether window should be used in Lowess.
    //         let chrom = str::from_utf8(reader.header().rid2name(record.rid().unwrap()))
    //             .unwrap()
    //             .to_string();
    //         let use_chrom = re.is_match(&chrom);

    //         if use_chrom && !is_gap && !is_blacklist {
    //             // TODO: remove restriction to one sample here
    //             let nrcs: Vec<f32> = vec![covs_gc[i] as f32; 1]; // XXX
    //             record
    //                 .push_format_float(b"NCOV_GC", nrcs.as_slice())
    //                 .expect("Could not write FORMAT/NCOV_GC");
    //             // let nrc2 = covs_gc[i].log2();
    //             // let nrcs2: Vec<f32> = vec![
    //             //     if nrc2.is_finite() {
    //             //         nrc2 as f32
    //             //     } else {
    //             //         f32::missing()
    //             //     };
    //             //     1
    //             // ]; // XXX
    //             // record
    //             //     .push_format_float(b"NCOV2_GC", nrcs2.as_slice())
    //             //     .expect("Could not write FORMAT/NCOV2_GC");

    //             // TODO: remove restriction to one sample here
    //             let nrcs: Vec<f32> = vec![covs_map[i] as f32; 1]; // XXX
    //             record
    //                 .push_format_float(b"NCOV", nrcs.as_slice())
    //                 .expect("Could not write FORMAT/NCOV");
    //             // let nrc2 = covs_map[i].log2();
    //             // let nrcs2: Vec<f32> = vec![
    //             //     if nrc2.is_finite() {
    //             //         nrc2 as f32
    //             //     } else {
    //             //         f32::missing()
    //             //     };
    //             //     1
    //             // ]; // XXX
    //             // record
    //             //     .push_format_float(b"NCOV2", nrcs2.as_slice())
    //             //     .expect("Could not write FORMAT/NCOV2");
    //             i += 1;
    //         }

    //         // Finally, write out record again.
    //         writer.write(&record).expect("Could not write BCF record!");

    //         prev_rid = Some(record.rid());
    //     }
    // }

    Ok(())
}

/// Main entry point for the "normalize" command.
pub fn call(logger: &mut Logger, options: &Options) -> Result<(), String> {
    // let mut app = CoverageApp::new(logger, options);
    info!(logger, "Running cnvetti normalize");
    info!(logger, "Configuration: {:?}", &options);

    // TODO: check availability of R.
    match &options.normalization {
        Normalization::BinnedGc => {
            normalize_binned_gc(logger, &options)?;
        }
        Normalization::LowessGc | Normalization::LowessGcMapability => {
            call_lowess(logger, &options)?;
        }
    }

    // Build index on the output file.
    build_index(logger, &options.output);

    info!(logger, "Done. Have a nice day!");

    Ok(())
}
