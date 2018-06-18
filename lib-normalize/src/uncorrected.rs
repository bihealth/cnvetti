/// Implementation normalization without GC correction.
use std::str;

use slog::Logger;

use rust_htslib::bcf::{self, Read};

use super::errors::*;
use options::*;

use lib_shared::stats::Stats;

use shared::*;

/// Compute total number of fragments per bp for normalization.
///
/// The input is already assumed to be length-normalized.
fn compute_format_tag_sum(logger: &mut Logger, input: &String, key: &[u8]) -> Result<f64> {
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

    Ok(vals.as_slice().sum())
}

pub fn run_uncorrected_normalization(
    logger: &mut Logger,
    options: &NormalizeOptions,
) -> Result<()> {
    info!(
        logger,
        "Computing overal normalized coverage sum for length-normalized coverage FORMAT/LCV"
    );
    let lcv_sum = compute_format_tag_sum(logger, &options.input, &b"LCV"[..])?;

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
