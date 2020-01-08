/// Main library module for merging coverage files.
use std::collections::BTreeSet;
use std::env;
use std::str;

extern crate clap;

extern crate chrono;

#[macro_use]
extern crate error_chain;

/// This crate's error-related code, generated by `error-chain`.
mod errors {
    // Create the Error, ErrorKind, ResultExt, and Result types
    error_chain! {}
}

pub use errors::*;

extern crate shlex;

#[macro_use]
extern crate slog;
use slog::Logger;

extern crate rust_htslib;
use rust_htslib::bcf;
use rust_htslib::bcf::record::Numeric;

extern crate lib_shared;
use lib_shared::bcf_utils;

mod options;
pub use options::*;

/// Obtain all field names of the given field type.
fn get_field_names(header: &bcf::header::HeaderView, field_type: &str) -> BTreeSet<String> {
    let mut result: BTreeSet<String> = BTreeSet::new();

    for record in header.header_records() {
        match record {
            bcf::HeaderRecord::Format { key: _, values } => match values.get("Type") {
                Some(this_type) => {
                    if this_type == field_type {
                        result.insert(
                            values
                                .get("ID")
                                .expect("FILTER entry does not have an ID!")
                                .clone(),
                        );
                    }
                }
                _ => (),
            },
            _ => (),
        }
    }

    result
}

/// Build BCF writer.
fn build_writer(
    logger: &mut Logger,
    header: &bcf::header::HeaderView,
    samples: &Vec<&[u8]>,
    options: &MergeCovOptions,
) -> Result<bcf::Writer> {
    debug!(logger, "Opening output file...");

    // Construct extended header.
    let mut header = bcf::Header::from_template_subset(header, &[])
        .chain_err(|| "Problem constructing header from template and no samples")?;
    // TODO: version should come from one central place
    header.push_record(format!("##cnvetti_cmdMergeCovVersion={}", "0.1.0").as_bytes());
    header.push_record(
        format!(
            "##cnvetti_cmdMergeCovCommand={}",
            env::args()
                .map(|s| shlex::quote(&s).to_string())
                .collect::<Vec<String>>()
                .join(" ")
        )
        .as_bytes(),
    );

    for sample in samples {
        header.push_sample(&sample);
    }

    let uncompressed = !options.output.ends_with(".bcf") && !options.output.ends_with(".vcf.gz");
    let vcf = options.output.ends_with(".vcf") || options.output.ends_with(".vcf.gz");

    Ok(
        bcf::Writer::from_path(&options.output, &header, uncompressed, vcf)
            .chain_err(|| format!("Could not open output BCF file {}", options.output))?,
    )
}

/// Merge the coverage BCF files in `reader` to the writer.
pub fn merge_files(
    format_string: &BTreeSet<String>,
    format_float: &BTreeSet<String>,
    reader: &mut bcf::synced::SyncedReader,
    writer: &mut bcf::Writer,
) -> Result<()> {
    while reader
        .read_next()
        .chain_err(|| "Problem reading from input BCF files")?
        != 0
    {
        // TODO: also merge INFO values?
        // Locate first record; will copy INFO from there; FORMAT comes later.
        let mut first: Option<bcf::Record> = None;
        let num_samples = writer.header().sample_count();
        for i in 0..num_samples {
            match (reader.has_line(i), &first) {
                (true, None) => {
                    let mut record = reader.record(i).expect("Could not retrieve record");
                    writer.translate(&mut record);
                    first = Some(record);
                    // break;
                }
                _ => (),
            }
        }

        assert!(first.is_some());
        let record = first.as_mut().unwrap();

        // Push GT, always no-call for coverage records.
        let values = (0..(1 * num_samples))
            .map(|_| bcf::GT_MISSING)
            .collect::<Vec<i32>>();
        record
            .push_format_integer(b"GT", values.as_slice())
            .chain_err(|| "Could not write FORMAT/GT")?;

        // Collect input FORMAT String fields and push to output FORMAT.
        for key in format_string {
            if key == "GT" {
                continue; // already handled above
            }
            let key_b = key.as_bytes();
            let mut values_v: Vec<Vec<u8>> = Vec::new();
            for i in 0..reader.reader_count() {
                if reader.has_line(i) {
                    let mut rec_in = reader
                        .record(i)
                        .expect("We just checked that the record should be there!");
                    let mut tmp_v = rec_in
                        .format(&key_b)
                        .string()
                        .unwrap_or_else(|_| Vec::new())
                        .iter()
                        .map(|v| Vec::from(*v))
                        .collect::<Vec<Vec<u8>>>();
                    if tmp_v.iter().any(|ref x| !x.is_empty()) {
                        values_v.append(&mut tmp_v);
                    }
                } else {
                    let num_samples = reader.header(i).sample_count();
                    for _j in 0..num_samples {
                        values_v.push(b"".to_vec());
                    }
                }
            }
            let values: Vec<&[u8]> = values_v
                .iter()
                .map(|v| v.as_slice())
                .collect::<Vec<&[u8]>>();
            record
                .push_format_string(key_b, values.as_slice())
                .chain_err(|| format!("Could not write FORMAT/{}", key))?;
        }

        // Collect input FORMAT Float fields and push to output FORMAT.
        for key in format_float {
            let key_b = key.as_bytes();
            // First, get dimension of individual arrays.
            let dim = (0..reader.reader_count())
                .filter(|i| reader.has_line(*i))
                .map(|i| {
                    reader
                        .record(i)
                        .expect("Could not get get record")
                        .format(&key_b)
                        .float()
                        .unwrap_or_else(|_| Vec::new())
                        .len()
                })
                .max()
                .expect("Could not compute maximum");
            if dim == 0 {
                continue;
            }
            // Then, build the array.
            let mut values: Vec<f32> = Vec::new();
            for i in 0..reader.reader_count() {
                if !reader.has_line(i) {
                    let header = reader.header(i);
                    let n = dim * header.sample_count() as usize;
                    values.append(&mut (0..n).map(|_| f32::missing()).collect::<Vec<f32>>());
                    continue;
                }
                match reader.record(i).unwrap().format(&key_b).float() {
                    Ok(ref vec) => {
                        for arr in vec {
                            values.append(&mut arr.to_vec());
                        }
                    }
                    Err(_) => {
                        let header = reader.header(i);
                        let n = dim * header.sample_count() as usize;
                        values.append(&mut (0..n).map(|_| f32::missing()).collect::<Vec<f32>>());
                    }
                }
            }
            // Finally, push array into record.
            record
                .push_format_float(key_b, values.as_slice())
                .chain_err(|| format!("Could not write FORMAT/{}", key))?;
        }

        writer
            .write(&record)
            .chain_err(|| "Problem writing BCF record")?;
    }

    Ok(())
}

/// Main entry point for the "cmd build-model-wis" command.
pub fn run(logger: &mut Logger, options: &MergeCovOptions) -> Result<()> {
    info!(logger, "Running: cnvetti cmd merge-cov");
    info!(logger, "Options: {:?}", options);

    // Open reader.
    info!(logger, "Opening input files...");
    let mut reader =
        bcf::synced::SyncedReader::new().chain_err(|| "Could not allocated synced reader")?;
    reader.set_require_index(true);
    reader.set_pairing(bcf::synced::pairing::EXACT);
    for input in &options.input {
        info!(logger, "- {}", input);
        reader
            .add_reader(&input)
            .chain_err(|| format!("Could not open file {} for reading", input))?;
    }
    info!(logger, "=> done");

    // Get all Format fields of type String and Float; "GT" will get special handling.  Also, get
    // all sample names.
    let mut format_string = BTreeSet::new();
    let mut format_float = BTreeSet::new();
    let mut samples_v: Vec<Vec<u8>> = Vec::new();
    for i in 0..reader.reader_count() {
        format_float.append(&mut get_field_names(&reader.header(i), "Float"));
        format_string.append(&mut get_field_names(&reader.header(i), "String"));
        let mut other = reader
            .header(i)
            .samples()
            .iter()
            .map(|s| s.to_vec())
            .collect::<Vec<Vec<u8>>>();
        samples_v.append(&mut other);
    }
    let samples = samples_v
        .iter()
        .map(|v| v.as_slice())
        .collect::<Vec<&[u8]>>();

    // Open output file; construct writer.
    // TODO: do something fancier than taking just the first header
    {
        let mut writer = build_writer(logger, &reader.header(0), &samples, &options)?;
        merge_files(&format_string, &format_float, &mut reader, &mut writer)?;
    }

    // Finally, create index on created output file.
    info!(logger, "Building index for output file...");
    bcf_utils::build_index(logger, &options.output).chain_err(|| "Could not build index")?;
    info!(logger, "All done. Have a nice day!");

    Ok(())
}
