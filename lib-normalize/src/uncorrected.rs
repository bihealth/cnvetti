//! Implementation normalizations without GC correction.
use std::str;

use slog::Logger;

use rust_htslib::bcf::{self, Read};

use super::errors::*;
use options::*;

use shared::*;

/// Compute summary from all FORMAT values of key `key` in BCF file `input`.
///
/// The input is already assumed to be length-normalized.
///
/// The summary is computed using `summary`.
fn compute_format_tag_sum<S>(
    logger: &mut Logger,
    input: &String,
    key: &[u8],
    summary: S,
) -> Result<f64>
where
    S: Fn(&[f64]) -> f64,
{
    debug!(
        logger,
        "Computing sum of {} in {}",
        str::from_utf8(key).unwrap(),
        &input
    );
    let mut reader = bcf::Reader::from_path(&input)
        .chain_err(|| format!("Could not open input BCF file {}", input))?;

    let mut vals = Vec::new(); // coverage values

    let mut record = reader.empty_record();
    loop {
        match reader.read(&mut record) {
            Ok(_) => (),
            Err(bcf::ReadError::NoMoreRecord) => break,
            _ => bail!("Error reading BCF record"),
        }
        vals.push(
            record
                .format(&key)
                .float()
                .chain_err(|| format!("Could not access FORMAT/{}", str::from_utf8(key).unwrap()))?
                [0][0] as f64,
        );
    }

    Ok(summary(vals.as_slice()))
}

/// Perform a simple normalization of the coverage values calling a function on the vector
/// of all coverages.
pub fn run_simple_normalization<S>(
    logger: &mut Logger,
    options: &NormalizeOptions,
    summary: S,
) -> Result<()>
where
    S: Fn(&[f64]) -> f64,
{
    info!(
        logger,
        "Computing overal normalized coverage sum for length-normalized coverage FORMAT/LCV"
    );
    let lcv_sum = compute_format_tag_sum(logger, &options.input, &b"LCV"[..], summary)?;

    // The processing closure.
    let process = |record: &mut bcf::Record| -> Result<()> {
        // Write normalized CV.
        let lcv = record
            .format(b"LCV")
            .float()
            .chain_err(|| "Problem getting value for LCV from BCF record")?[0][0]
            as f64;
        record
            .push_format_float(b"CV", &[(lcv / lcv_sum) as f32])
            .chain_err(|| "Problem writing CV to BCF record")?;
        record
            .push_format_float(b"CV2", &[(lcv / lcv_sum).log2() as f32])
            .chain_err(|| "Problem writing CV to BCF record")?;

        // Write normalized CV standard deviation.
        let is_ok = record.format(b"LCVSD").float().is_ok();
        if is_ok {
            let lcvsd = record.format(b"LCVSD").float().unwrap()[0][0];
            record
                .push_format_float(b"CVSD", &[(lcvsd as f64 / lcv_sum) as f32])
                .chain_err(|| "Problem writing CVSD to BCF record")?;
        }

        Ok(())
    };

    info!(logger, "Reading input, normalizing, writing to output");
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
