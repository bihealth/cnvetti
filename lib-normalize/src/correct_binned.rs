/// Implementation of binning-based GC correction and normalization.
use quantiles::ckms::CKMS;

use std::str;

use slog::Logger;

use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::{self, Read};

use super::errors::*;
use options::*;

use shared::*;

// TODO: the following should be moved to options
/// Allowed error in median computation.
const ERROR: f64 = 1e-6;
/// Minimal number of entries for computing median.
const MIN_COUNT: usize = 10;

/// Compute coverage medians for each GC bin.
///
/// We are computing a percent-wise normalization, such that the result will have 101 entries.
fn compute_format_tag_median(logger: &mut Logger, input: &String, key: &[u8]) -> Result<Vec<f32>> {
    debug!(
        logger,
        "medians of {} in {}",
        str::from_utf8(key).unwrap(),
        &input
    );

    let mut reader = bcf::Reader::from_path(&input)
        .chain_err(|| format!("Could not open input BCF file {}", input))?;

    let mut processors: Vec<CKMS<f32>> = vec![CKMS::<f32>::new(ERROR); 101];

    let mut record = reader.empty_record();
    loop {
        match reader.read(&mut record) {
            Ok(_) => (),
            Err(bcf::ReadError::NoMoreRecord) => break,
            _ => bail!("Error reading BCF record"),
        }

        let gc = record
            .info(b"GC")
            .float()
            .chain_err(|| {
                "Could not access INFO/GC. Did you forget to pass in reference when \
                 computing coverage so GC content was not computed?"
            })?
            .map(|gcs| gcs[0] as f32)
            .unwrap_or(f32::missing());
        if gc.is_finite() {
            let gc = (gc * 100.0).floor() as usize;
            processors[gc].insert(
                record.format(&key).float().chain_err(|| {
                    format!("Could not access FORMAT/{}", str::from_utf8(key).unwrap())
                })?[0][0],
            );
        }
    }

    Ok(processors
        .iter()
        .map(|ckms| {
            if ckms.count() < MIN_COUNT {
                f32::missing()
            } else {
                ckms.query(0.5).unwrap_or((0, f32::missing())).1
            }
        })
        .collect::<Vec<f32>>())
}

pub fn run_binned_correction(logger: &mut Logger, options: &NormalizeOptions) -> Result<()> {
    info!(logger, "Computing GC-wise medians");
    let gc_medians = compute_format_tag_median(logger, &options.input, &b"LCV"[..])?;

    // The processing closure.
    let process = |record: &mut bcf::Record| -> Result<()> {
        // Write normalized CV.
        let lcv = record
            .format(b"LCV")
            .float()
            .chain_err(|| "Problem getting value for LCV from BCF record")?[0][0];

        let gc = record
            .info(b"GC")
            .float()
            .chain_err(|| "Could not access INFO/GC")?
            .map(|gcs| gcs[0] as f32)
            .unwrap_or(f32::missing());

        if gc.is_finite() {
            let gc = (gc * 100.0).floor() as usize;
            let median = gc_medians[gc];
            if !median.is_finite() {
                record
                    .push_info_flag(b"FEW_GCWINDOWS")
                    .chain_err(|| "Could not write INFO/FEW_GCWINDOWS")?;
            } else {
                let norm_lcv = lcv / median as f32;
                record
                    .push_format_float(b"CV", &[norm_lcv])
                    .chain_err(|| "Problem writing CV to BCF record")?;
                if norm_lcv.is_finite() {
                    let cov2 = if norm_lcv == 0.0 {
                        f32::missing()
                    } else {
                        norm_lcv.log2() as f32
                    };
                    record
                        .push_format_float(b"CV2", &[cov2])
                        .chain_err(|| "Could not write FORMAT/CV2")?;
                }
                // Write normalized CV standard deviation.
                let is_ok = record.format(b"LCVSD").float().is_ok();
                if is_ok {
                    let lcvsd = record.format(b"LCVSD").float().unwrap()[0][0];
                    record
                        .push_format_float(b"CVSD", &[lcvsd / median as f32])
                        .chain_err(|| "Problem writing CVSD to BCF record")?;
                }
            }
        }

        Ok(())
    };

    info!(
        logger,
        "Reading input, normalizing (+GC correction), writing output"
    );
    process_bcf(
        logger,
        &options.input,
        &options.output,
        options.io_threads,
        &process,
    )?;
    info!(logger, "Done writing output");

    Ok(())
}
