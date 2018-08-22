//! Implementation of EXCAVATOR-style EMRC.
use std::collections::HashMap;
use std::env;

use shlex;

use slog::Logger;

use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::{self, Read};

use super::errors::*;
use options::*;

const EXON_SIZE_BIN_SIZE: i32 = 10;

/// Tuple describing one factor combination.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct EmrcFactor {
    /// The percentage of GC content.
    pub gc: i32,
    /// The rounded mean MQ of the exon.
    pub mapq_mean: i32,
    /// The exon size .
    pub exon_size: i32,
}

impl EmrcFactor {
    pub fn new(gc: i32, mapq_mean: i32, exon_size: i32) -> Self {
        Self {
            gc,
            mapq_mean,
            exon_size,
        }
    }
}

/// Tuple for EMRC correction.
#[derive(Debug, Clone)]
struct EmrcTuple {
    /// Key.
    pub factor: EmrcFactor,
    /// Length-normalized read count.
    pub emrc: f64,
}

/// Load the raw metrics tuples of EMRC values.
fn load_emrc_tuples(_logger: &mut Logger, input: &String) -> Result<Vec<EmrcTuple>> {
    let mut result = Vec::new();
    let mut reader = bcf::Reader::from_path(&input)
        .chain_err(|| format!("Could not open input BCF file {}", input))?;

    let mut record = reader.empty_record();
    loop {
        match reader.read(&mut record) {
            Ok(_) => record.unpack(),
            Err(bcf::ReadError::NoMoreRecord) => break,
            _ => bail!("Error reading BCF record"),
        }

        let gc = (record
            .info(b"GC")
            .float()
            .expect("Could not access INFO/GC")
            .unwrap()[0] * 100.0)
            .round() as i32;
        let mapq_mean = record
            .format(b"MQ")
            .float()
            .expect("Could not access FORMAT/MQ")[0][0]
            .round() as i32;
        let emrc = record
            .format(b"LCV")
            .float()
            .expect("Could not access FORMAT/LCV")[0][0]
            .round() as f64;

        if !emrc.is_finite() {
            continue;
        }

        let pos = record.pos() as i32;
        let end = record
            .info(b"END")
            .integer()
            .expect("Could not access INFO/END")
            .unwrap()[0];
        let exon_size = (end - pos) / EXON_SIZE_BIN_SIZE;

        result.push(EmrcTuple {
            factor: EmrcFactor {
                gc,
                mapq_mean,
                exon_size,
            },
            emrc,
        });
    }

    Ok(result)
}

/// Write out corrected results.
fn write_corrected(
    logger: &mut Logger,
    factors: &HashMap<EmrcFactor, f64>,
    options: &NormalizeOptions,
) -> Result<()> {
    debug!(logger, "Opening input file...");
    let mut reader = bcf::Reader::from_path(&options.input)
        .chain_err(|| format!("Could not open BCF file: {}", options.input))?;
    if options.io_threads > 0 {
        reader
            .set_threads(options.io_threads as usize)
            .chain_err(|| "Could not set threads to BCF reader")?;
    }

    debug!(logger, "Opening output file...");
    let mut writer = {
        // Construct extended header.
        let mut header = bcf::Header::from_template(reader.header());
        // TODO: version should come from one central place
        header.push_record(format!("##cnvetti_cmdNormalizeVersion={}", "0.2.0").as_bytes());
        header.push_record(
            format!(
                "##cnvetti_wiseNormalizeCommand={}",
                env::args()
                    .map(|s| shlex::quote(&s).to_string())
                    .collect::<Vec<String>>()
                    .join(" ")
            ).as_bytes(),
        );

        let uncompressed =
            !options.output.ends_with(".bcf") && !&options.output.ends_with(".vcf.gz");
        let vcf = options.output.ends_with(".vcf") || options.output.ends_with(".vcf.gz");

        bcf::Writer::from_path(options.output.clone(), &header, uncompressed, vcf)
            .chain_err(|| format!("Could not open output BCF file {}", &options.output))?
    };

    let mut record = reader.empty_record();
    loop {
        match reader.read(&mut record) {
            Ok(_) => (),
            Err(bcf::ReadError::NoMoreRecord) => break,
            _ => bail!("Problem reading BCF record"),
        }
        writer.translate(&mut record);

        let pos = record.pos() as i32;
        let end = record
            .info(b"END")
            .integer()
            .expect("Could not access INFO/END")
            .unwrap()[0];
        let exon_size = end - pos;

        // TOOD: look into LCV/EMCR here!
        // let lcv = record
        //     .format(b"LCV")
        //     .float()
        //     .expect("Could not access INFO/LCV")[0][0] as f64;
        let gc = (record
            .info(b"GC")
            .float()
            .expect("Could not access INFO/GC")
            .unwrap()[0] * 100.0)
            .round() as i32;
        let mapq_mean = record
            .format(b"MQ")
            .float()
            .expect("Could not access FORMAT/MQ")[0][0]
            .round() as i32;
        let emrc = record
            .format(b"LCV")
            .float()
            .expect("Could not access FORMAT/LCV")[0][0]
            .round() as f64;

        let cv = *factors
            .get(&EmrcFactor::new(
                gc,
                mapq_mean,
                exon_size / EXON_SIZE_BIN_SIZE,
            ))
            .expect(&format!(
                "Cannot obtain EMRC factor {:?}",
                EmrcFactor::new(gc, mapq_mean, exon_size / EXON_SIZE_BIN_SIZE)
            )) as f64;
        if !emrc.is_finite() {
            record
                .push_format_float(b"CV", &[f32::missing()])
                .chain_err(|| "Problem writing CV to BCF record")?;
        } else {
            record
                .push_format_float(b"CV", &[(cv / emrc) as f32])
                .chain_err(|| "Problem writing CV to BCF record")?;
        }

        writer
            .write(&record)
            .chain_err(|| "Problem writing BCF record to file")?;
    }

    Ok(())
}

/// Perform the exon-mean-read-count correction.
pub fn run_emrc_correction(logger: &mut Logger, options: &NormalizeOptions) -> Result<()> {
    info!(logger, "Loading raw metrics for EMRC...");
    let emrc_tuples =
        load_emrc_tuples(logger, &options.input).chain_err(|| "Could not load EMRC raw metrics")?;

    info!(logger, "Computing normalization factors");
    let factors = {
        let mut sums: HashMap<EmrcFactor, f64> = HashMap::new();
        let mut counts: HashMap<EmrcFactor, usize> = HashMap::new();
        for tuple in &emrc_tuples {
            *sums.entry(tuple.factor).or_insert(0.0) += tuple.emrc;
            *counts.entry(tuple.factor).or_insert(0) += 1;
        }
        sums.iter()
            .map(|(k, sum)| (*k, *sum / (*counts.get(k).unwrap() as f64)))
            .collect::<HashMap<_, _>>()
    };

    write_corrected(logger, &factors, options)
        .chain_err(|| "Problem writing out corrected counts")?;

    Ok(())
}
