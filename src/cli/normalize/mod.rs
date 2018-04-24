// Implementation of the "normalize" command.

include!(concat!(env!("OUT_DIR"), "/version.rs"));

mod options;
mod r_scripts;

use std::env;
use std::fs::File;
use std::process::Command;
use std::str;

use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::{self, Read as BcfRead};

use regex::Regex;

use serde_json;

use csv;

use shlex;

use slog::Logger;

use handlebars::Handlebars;

use tempdir::TempDir;

pub use self::options::*;
use cli::coverage::summary::SummarisedMetric;
use cli::shared;

use histogram::Histogram;

// use std::io;
use std::io::prelude::*;

// fn pause() {
//     let mut stdin = io::stdin();
//     let mut stdout = io::stdout();

//     // We want the cursor to stay at the end of the line, so we print without a newline and flush manually.
//     write!(stdout, "Press any key to continue...").unwrap();
//     stdout.flush().unwrap();

//     // Read a single byte and discard
//     let _ = stdin.read(&mut [0u8]).unwrap();
// }

// TODO: check input file.
// TODO: use index-based readers, is nicer for progress display...
// TODO: could normalize already in coverage step by writing raw coverage to temporary file

/// Perform the actual processing for binned GC content.
fn process_binned_gc(
    reader: &mut bcf::Reader,
    writer: &mut bcf::Writer,
    logger: &mut Logger,
    summaries: Vec<SummarisedMetric>,
    options: &Options,
) {
    let mut record = reader.empty_record();
    let mut prev_rid = Option::None;

    let header = reader.header().clone();

    let sample_count = reader.header().sample_count() as usize;
    let mut hists = vec![Histogram::new(); sample_count];

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
                .unwrap_or(&[0])[0];
            if gc_windows < options.min_gc_window_count {
                record.push_filter(
                    header
                        .name_to_id(b"FEW_GCWINDOWS")
                        .expect("FEW_GCWINDOWS not known in header"),
                );
            }
        }

        // Get INFO/GC.
        let gc = {
            record
                .info(b"GC")
                .float()
                .expect("Could not read INFO/GC")
                .expect("INFO/GC was empty")[0]
        };
        let gc_bucket = ((gc as f64 / options.gc_step).floor() * 1000.0 * options.gc_step) as i32;

        let sample_id = 0;

        // Normalize FORMAT/RC and FORMAT/{COV,SD}.
        match options.count_kind {
            CountKind::Coverage => {
                panic!("Not implemented yet!");
            }
            CountKind::Alignments => {
                if let Some(summary) = &summaries[sample_id].summaries.get(&gc_bucket) {
                    let median = summary.summary5[2] as f32;
                    let nrcs: Vec<f32> = record
                        .format(b"RC")
                        .integer()
                        .expect("Cannot access FORMAT/RC")
                        .iter()
                        .map(|slice| {
                            if median < 0.001 {
                                // guard against NaN
                                0.0
                            } else {
                                println!(
                                    "RC={}, median={}, RC/median={}",
                                    slice[0],
                                    median,
                                    slice[0] as f32 / median
                                );
                                slice[0] as f32 / median
                            }
                        })
                        .collect();
                    for (idx, nrc) in nrcs.iter().enumerate() {
                        let val = (nrc * 1000.0) as u64;
                        hists[idx].increment(val + 1).unwrap();
                    }
                    record
                        .push_format_float(b"NRC", nrcs.as_slice())
                        .expect("Could not write FORMAT/NRC");
                } else {
                    record
                        .push_format_float(b"NRC", vec![0.0_f32; sample_count].as_slice())
                        .expect("Could not write FORMAT/NRC");
                }
            }
        }

        // Finally, write out record again.
        writer.write(&record).expect("Could not write BCF record!");

        prev_rid = Some(record.rid());
    }

    info!(logger, "Maximal absolute deviations");
    for (idx, hist) in hists.iter().enumerate() {
        let median = hist.mean().unwrap(); // TODO: should be median
        let mut hist2 = Histogram::new();
        for bucket in hist {
            let val = (bucket.value() as i64 - median as i64).abs() as u64;
            hist2.increment_by(val + 1, bucket.count()).unwrap();
        }

        info!(
            logger,
            "{} => {}",
            idx,
            (hist2.mean().unwrap() - 1) as f64 / 1000.0
        );
    }
}

/// Run normalization with binned GC content.
pub fn call_binned_gc(logger: &mut Logger, options: &Options) -> Result<(), String> {
    info!(logger, "Normalizing by binned GC count...");
    info!(logger, "Opening input and output files...");
    debug!(logger, "Loading statistics");
    let stats: Vec<SummarisedMetric> = {
        let stats_filename = options.input.clone() + ".stats.txt";
        let mut file = File::open(stats_filename).expect("Could not open stats file for reading");
        let mut json = String::new();
        file.read_to_string(&mut json)
            .expect("Could not read from stats file");
        serde_json::from_str(&json).expect("Could not deserialize stats from JSON")
    };

    debug!(logger, "Opening input file");
    let mut reader =
        bcf::Reader::from_path(options.input.clone()).expect("Could not open input BCF file");
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
        process_binned_gc(&mut reader, &mut writer, logger, stats, &options);
    }

    Ok(())
}

/// Main entry point for LOESS-based normalization.
pub fn call_loess(logger: &mut Logger, options: &Options) -> Result<(), String> {
    info!(logger, "Normalizing by LOESS via R...");

    // TODO: Remove limitation to one sample only!

    let tmp_dir = TempDir::new("cnvetti_normalize").expect("Could not create temporary directory");
    // TODO: explicitely drop tmp_dir and look whether tmp_dir is properly cleaned
    let tsv_path = tmp_dir.path().join("input.tsv"); // TODO: rename to input_path
    let output_path = tmp_dir.path().join("output.tsv");
    debug!(logger, "output path is {}", output_path.to_str().unwrap());
    let script_path = tmp_dir.path().join("loess.R");

    info!(logger, "Writing script file");
    debug!(
        logger,
        "Script file path is {}",
        script_path.to_str().unwrap()
    );
    // Block for script file writer.
    {
        let mut script_file =
            File::create(script_path.clone()).expect("Could not create temporary R script file");
        let mut handlebars = Handlebars::new();
        handlebars
            .register_template_string("loess.R", r_scripts::LOESS_R)
            .expect("Could not register loess_map templates");
        handlebars
            .render_to_write(
                "loess.R",
                &json!({
                    "input_file": tsv_path,
                    "output_file": output_path,
                    "loess_mapability": "TRUE",
                    "LOG2_TRANSFORM": "FALSE",
                }),
                &mut script_file,
            )
            .expect("Could not write script file");
    }

    info!(logger, "Writing data file");
    debug!(logger, "TSV file path is {}", tsv_path.to_str().unwrap());
    // Block for BCF reader and data file writer.
    {
        debug!(logger, "Opening input BCF file for writing LOESS input");
        let mut reader =
            bcf::Reader::from_path(options.input.clone()).expect("Could not open input BCF file");
        if options.io_threads > 0 {
            reader
                .set_threads(options.io_threads as usize)
                .expect("Could not set I/O thread count");
        }

        let mut tsv_file = File::create(tsv_path).expect("Could not create temporary TSV file");
        // let samples: Vec<String> = reader
        //     .header()
        //     .samples()
        //     .iter()
        //     .map(|slice| str::from_utf8(slice).unwrap().to_string())
        //     .collect();

        writeln!(tsv_file, "use_loess\tgc_content\tmapability\tcount")
            .expect("Could not write TSV header");

        let re = Regex::new(&options.contig_regex).unwrap();

        let mut record = reader.empty_record();
        while reader.read(&mut record).is_ok() {
            // Get "is gap" flag.
            let is_gap = {
                record
                    .info(b"GAP")
                    .integer()
                    .expect("Could not read INFO/GAP")
                    .expect("INFO/GAP was empty")[0] != 0
            };

            // Check whether window should be used in LOESS.
            let chrom = str::from_utf8(reader.header().rid2name(record.rid().unwrap()))
                .unwrap()
                .to_string();
            let use_chrom = re.is_match(&chrom);

            let use_loess = use_chrom && !is_gap;
            let use_loess = if use_loess { "T" } else { "F" };

            // Get INFO/GC.
            let gc_content = {
                record
                    .info(b"GC")
                    .float()
                    .expect("Could not read INFO/GC")
                    .expect("INFO/GC was empty")[0]
            };
            // TODO: add switch whether mapability is desired, also for when writing out script above.
            // Get INFO/MAPABILITY
            let mapability = {
                match record.info(b"MAPABILITY").float() {
                    Ok(ref mapability) => match mapability {
                        Some(ref mapability) => {
                            if mapability.len() > 0 {
                                mapability[0]
                            } else {
                                0_f32
                            }
                        }
                        None => 0_f32,
                    },
                    Err(_) => 0_f32,
                }
            };

            write!(tsv_file, "{}\t{}\t{}\t", use_loess, gc_content, mapability)
                .expect("Could not write to file");

            match options.count_kind {
                CountKind::Coverage => {
                    panic!("Not implemented yet!");
                }
                CountKind::Alignments => {
                    let rc = record
                        .format(b"RC")
                        .integer()
                        .expect("FORMAT/RC not found!");
                    let vals: Vec<String> = rc.iter().map(|ref rcs| rcs[0].to_string()).collect();
                    writeln!(tsv_file, "{}", vals.join("\t")).expect("Could not write to file");
                }
            }
        }
    }

    // Call out to R.
    Command::new("Rscript")
        .args(&["--vanilla", "--verbose", &script_path.to_str().unwrap()])
        .spawn()
        .expect("Failed to start loess.R script")
        .wait()
        .expect("loess.R script stopped with an error");

    // pause();

    // Block for BCF file reader and writer.
    {
        debug!(
            logger,
            "Opening input BCF file for generating final BCF file"
        );
        let mut reader =
            bcf::Reader::from_path(options.input.clone()).expect("Could not open input BCF file");
        if options.io_threads > 0 {
            reader
                .set_threads(options.io_threads as usize)
                .expect("Could not set I/O thread count");
        }

        let mut writer = {
            // Construct extended header.
            let mut header = bcf::Header::with_template(reader.header());
            let lines = vec![
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
        let mut prev_rid = Option::None;

        let header = reader.header().clone();

        let mut tsv_reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(output_path)
            .expect("Could not open TSV file");

        while reader.read(&mut record).is_ok() {
            if prev_rid.is_some() && prev_rid.unwrap() != record.rid() {
                info!(
                    logger,
                    "Starting on contig {}",
                    str::from_utf8(header.rid2name(record.rid().unwrap())).unwrap()
                );
            }

            let tsv_record: Vec<String> = tsv_reader
                .records()
                .next()
                .expect("Could not read next TSV record")
                .expect("Problem parsing TSV record")
                .iter()
                .map(|x| x.to_string())
                .collect();
            let nrc = tsv_record.last().expect("Empty TSV record");
            // println!("nrc = {}", nrc);
            let nrc = if nrc == "NA" || nrc == "Inf" || nrc == "-Inf" {
                f32::missing()
            } else {
                nrc.parse::<f32>().expect("Could not parse NRC")
            };

            // Translate to output header.
            writer.translate(&mut record);

            // TODO: remove restriction to one sample here
            let nrcs: Vec<f32> = vec![nrc; 1]; // XXX
            record
                .push_format_float(b"NRC", nrcs.as_slice())
                .expect("Could not write FORMAT/NRC");

            // Finally, write out record again.
            writer.write(&record).expect("Could not write BCF record!");

            prev_rid = Some(record.rid());
        }
    }

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
            call_binned_gc(logger, &options)?;
        }
        Normalization::LoessGc | Normalization::LoessGcMapability => {
            call_loess(logger, &options)?;
        }
    }

    // Build index on the output file.
    shared::build_index(logger, &options.output);

    info!(logger, "Done. Have a nice day!");

    Ok(())
}
