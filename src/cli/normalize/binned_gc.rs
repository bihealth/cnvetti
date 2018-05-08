// Normalization based on binned GC.

use super::shared::{load_norm_data, write_normalized_data};

use slog::Logger;

use rust_htslib::bcf::record::Numeric;

// use rust_htslib::bcf::{self, Read as BcfRead};

use cli::normalize::options::Options;
use cli::normalize::shared::NormData;

use histogram::Histogram;

/// Compute coverage histograms for each bin.
fn compute_histos(norm_data: &Vec<NormData>, options: &Options) -> Result<Vec<Histogram>, String> {
    let mut result: Vec<Histogram> = Vec::new();

    for ref record in norm_data {
        let bin = record.gc_bin(options.gc_step);
        if record.use_record(0.6) {
            // Resize histogram Vec if necessary.
            if result.len() <= bin {
                result.resize(bin + 1, Histogram::new());
            }
            // Register for coverage (rounded by one decimal place) for GC bin.
            if record.coverage.is_finite() {
                // println!("Incrementing {}", (record.coverage * 10.0 + 1.0).round() as u64);
                result[bin]
                    .increment((record.coverage * 10.0 + 1.0).round() as u64)
                    .map_err(|e| format!("Problem incrementing histogram: {}", e))?;
            }
        }
    }

    Ok(result)
}

/// Normalize coverage based on GC bins.
///
/// This is mostly useful for WGS data as it's fast but needs large number of bins.
pub fn normalize_binned_gc(logger: &mut Logger, options: &Options) -> Result<(), String> {
    info!(logger, "Normalizing by GC bins");

    // Load the data that is to be used for the input of normalization.
    let norm_data = load_norm_data(logger, options)?;

    // Compute histogram for normalization.
    debug!(logger, "Computing histogram for normalization...");
    let histos = compute_histos(&norm_data, options)?;

    // Compute the medians.
    debug!(logger, "Computing medians...");
    let num_bins = (100.0 / options.gc_step).ceil() as usize + 1;
    let medians: Vec<f64> = (0..num_bins)
        .map(|bin| {
            if bin >= histos.len() {
                0.0_f64
            } else {
                histos[bin].percentile(50.0).unwrap_or(0) as f64 / 10.0
            }
        })
        .collect();

    // Compute normalized coverage, normalized means and outlier flags.
    debug!(logger, "Computing normalized coverages...");
    // let mut abs_devs = vec![Histogram::new(); num_bins];
    let normalized: Vec<(f32, f32)> = norm_data
        .iter()
        .map(|record| {
            let bin = record.gc_bin(options.gc_step);
            if medians[bin] > 0.0 {
                // let val = record.coverage / medians[bin];
                // abs_devs[bin]
                //     .increment(scaled_dev(val))
                //     .expect("Could not increment bucket");
                // val
                (
                    (record.coverage / medians[bin]) as f32,
                    (record.coverage_sd / medians[bin]) as f32,
                )
            } else {
                (f32::missing(), f32::missing())
            }
        })
        .collect();

    // Write out the normalized data again.
    write_normalized_data(logger, options, &normalized)?;

    Ok(())
}
