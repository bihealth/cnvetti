// Implementation of the HaarSeg algorithm.

// Based on HaarSeg.c for which the license header is reproduced below

/*
 *    Copyright (C) 2008  Erez Ben-Yaacov
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *    <http://www.gnu.org/licenses/>
 *
 */

use statrs::function::erf;

use slog::Logger;

use cli::segment::math::median_abs_dev;

fn haar_convolution(
    signal: &[f64],
    weight: Option<&[f64]>,
    step_half_size: usize,
    logger: &Logger,
) -> Vec<f64> {
    let signal_size = signal.len();
    if signal_size < 50 {
        debug!(
            logger,
            "signal_size = {}, signal = {:?}", signal_size, signal
        );
    }

    let mut result = vec![0.0; signal_size];

    if step_half_size > signal_size {
        warn!(
            logger,
            "Error? step_half_size = {}, signal_size = {}", step_half_size, signal_size
        );
        return result;
    }

    let mut high_weight_sum: f64 = 0.0;
    let mut high_non_normed: f64 = 0.0;
    let mut low_weight_sum: f64 = 0.0;
    let mut low_non_normed: f64 = 0.0;

    if let Some(weight) = weight {
        // Initialize weight sums.
        for k in 0..step_half_size {
            high_weight_sum += weight[k];
            high_non_normed += weight[k] * signal[k];
        }
        // Circular padding.
        low_weight_sum = high_weight_sum;
        low_non_normed = -high_non_normed;
    }

    for k in 1..signal_size {
        // This loop is the hotspot.
        let high_end = if k + step_half_size >= signal_size + 1 {
            2 * signal_size - k - step_half_size
        } else {
            k + step_half_size - 1
        };
        let low_end = if step_half_size + 1 > k {
            step_half_size + 1 - k
        } else {
            k - step_half_size - 1
        };

        if let Some(weight) = weight {
            low_non_normed += signal[low_end] * weight[low_end] - signal[k - 1] * weight[k - 1];
            high_non_normed += signal[high_end] * weight[high_end] - signal[k - 1] * weight[k - 1];
            low_weight_sum += weight[k - 1] - weight[low_end];
            high_weight_sum += weight[high_end] - weight[k - 1];
            result[k] = (low_non_normed / low_weight_sum + high_non_normed / high_weight_sum)
                * (step_half_size as f64 / 2.0)
                * (step_half_size as f64 / 2.0);
        } else {
            result[k] = result[k - 1] + signal[high_end] + signal[low_end] - 2.0 * signal[k - 1];
        }
    }

    if weight.is_none() {
        let step_norm: f64 = 4.0 * step_half_size as f64 * step_half_size as f64;
        for k in 1..signal_size {
            result[k] /= step_norm;
        }
    }

    result
}

/// Find positions of local peaks in the given `signal` and return them.
fn find_local_peaks(signal: &[f64]) -> Vec<usize> {
    let mut min_suspect: Option<usize> = None;
    let mut max_suspect: Option<usize> = None;
    let mut result: Vec<usize> = Vec::new();

    for k in 1..(signal.len() - 1) {
        // first and last are not considered
        let sig_prev = signal[k - 1];
        let sig_curr = signal[k];
        let sig_next = signal[k + 1];

        if sig_curr > 0.0 {
            // look for local maxima
            if sig_curr > sig_prev && sig_curr > sig_next {
                result.push(k as usize);
            } else if sig_curr > sig_prev && sig_curr == sig_next {
                max_suspect = Some(k);
            } else if sig_curr == sig_prev && sig_curr > sig_next {
                if let Some(pos) = max_suspect {
                    result.push(pos);
                    max_suspect = None;
                }
            } else if sig_curr == sig_prev && sig_curr < sig_next {
                max_suspect = None;
            }
        } else if sig_curr < 0.0 {
            // Look for local maxima.
            if sig_curr < sig_prev && sig_curr < sig_next {
                result.push(k);
            } else if sig_curr < sig_prev && sig_curr == sig_next {
                min_suspect = Some(k);
            } else if sig_curr == sig_prev && sig_curr < sig_next {
                if let Some(pos) = min_suspect {
                    result.push(pos);
                    min_suspect = None;
                }
            } else if sig_curr == sig_prev && sig_curr > sig_next {
                min_suspect = None;
            }
        }
    }

    result
}

/// Perform breakpoint merging as described in Ben-Yaacov et al. (2008).
fn merge_breakpoints(
    prev_breakpoints: &[usize],
    new_breakpoints: &[usize],
    window_size: usize,
) -> Vec<usize> {
    if new_breakpoints.is_empty() {
        prev_breakpoints.to_vec()
    } else {
        let mut result: Vec<usize> = Vec::new();

        let mut it_new = new_breakpoints.iter();
        for it_prev in prev_breakpoints.iter() {
            loop {
                if let Some(val_new) = it_new.next() {
                    if *val_new > *it_prev + window_size {
                        break;
                    } else if *val_new + window_size < *it_prev {
                        result.push(*val_new);
                    }
                } else {
                    break;
                }
            }
            result.push(*it_prev);
        }

        result.extend(it_new);

        result
    }
}

// Perform pulse convolution.
fn pulse_convolution(signal: &[f64], pulse_size: usize) -> Vec<f64> {
    let signal_size = signal.len();
    assert!(
        pulse_size <= signal_size,
        "pulse_size = {}, signal_size = {}",
        pulse_size,
        signal_size
    );
    let pulse_height: f64 = 1.0 / pulse_size as f64;

    // Initialize circular paddig.
    let mut result: Vec<f64> = vec![0.0; signal_size];
    for k in 0..((pulse_size + 1) / 2) {
        result[0] += signal[k];
    }
    for k in 0..(pulse_size / 2) {
        result[0] += signal[k];
    }
    result[0] *= pulse_height;

    let mut n = 1;
    for k in (pulse_size / 2)..(signal_size + (pulse_size / 2) - 1) {
        let tail = if pulse_size <= k {
            k - pulse_size
        } else {
            pulse_size as usize - k as usize - 1
        };
        let head = if k < signal_size {
            k
        } else {
            2 * signal_size - 1 - k
        };
        result[n] = result[n - 1] + (signal[head] - signal[tail]) * pulse_height;
        n += 1;
    }

    result
}

/// Cumulative density function of the normal distribution
fn pnorm(mean: f64, sd: f64) -> f64 {
    0.5 * (1.0 + erf::erf(mean / sd / 2.0_f64.sqrt()))
}

/// Compute FDR threshold.
fn fdr_thresh(x: &[f64], q: f64, stdev: f64, logger: &Logger) -> f64 {
    let x_len = x.len();
    if x_len < 2 {
        return 0.0;
    }

    let mut m = vec![0.0; x_len];
    for i in 0..(m.len()) {
        m[i] = (i + 1) as f64 / x_len as f64;
    }

    let mut x_sorted: Vec<f64> = Vec::with_capacity(x.len());
    x_sorted.extend(x.iter().map(|x| (*x).abs()));
    x_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // The following code means "p <- pnorm(x_sorted);"
    let mut p: Vec<f64> = Vec::with_capacity(x_sorted.len());
    for x in x_sorted.iter() {
        p.push(2.0 * (1.0 - pnorm(*x, stdev)));
    }

    // Get the largets index for which p <= m * q;
    let mut idx: Option<usize> = None;
    for i in 0..(p.len()) {
        if p[i] < m[i] * q {
            idx = Some(i);
        }
    }

    if let Some(idx) = idx {
        x_sorted[idx]
    } else {
        error!(
            logger,
            "No passing p-vaulues: min p={}, min m={}, q={}", p[0], m[0], q
        );
        x_sorted[0] + 1e-16 // 1.e-16 almost 2^-52, like MATLAB "eps"
    }
}

// Signal standard deviation estimates and data for non-stationary variance estimation
struct SigmaEstimates {
    // Estimate of standard deviatiation in signal for peaks.
    peak_est: f64,
    // Estimate of standard deviation in signal for noise, only used when using
    // the raw values for correcting non-stationary variance.
    noise_est: f64,
    // Mask where the raw values are < thresh_nsv (b[n] in the HaarSeg paper).
    var_mask: Option<Vec<f64>>,
}

// Perform estimation of sigmas and compute other values used when raw values are available
// for compensation of non-stationary variance compensation
fn estimate_sigmas(vals: &[f64], raw_vals: Option<&[f64]>, logger: &Logger) -> SigmaEstimates {
    // TODO: should come from options
    let thresh_nsv = 50.0; // using same as in HaarSeg paper
    let diff_vals = haar_convolution(vals, Option::None, 1, logger);

    let (var_mask, peak_est, noise_est) = if let Some(raw_vals) = raw_vals {
        let mut var_mask: Vec<f64> = Vec::with_capacity(vals.len());
        var_mask.extend(
            raw_vals
                .iter()
                .map(|x| if *x < thresh_nsv { 1.0 } else { 0.0 }),
        );

        // Perform pulse convolution.
        let pulse_size = 2;
        let result_conv = pulse_convolution(&var_mask, pulse_size);

        // Estimate peak and noise.
        let mut for_peak: Vec<f64> = Vec::new();
        let mut for_noise: Vec<f64> = Vec::new();
        for x in result_conv {
            if x >= 0.5 {
                for_peak.push(x);
            } else {
                for_noise.push(x);
            }
        }

        (
            Some(var_mask),
            median_abs_dev(for_peak.as_slice()),
            median_abs_dev(for_noise.as_slice()),
        )
    } else {
        (None, median_abs_dev(&diff_vals), 0.0)
    };

    SigmaEstimates {
        peak_est,
        noise_est,
        var_mask,
    }
}

// Perform compensation for non-stationary variance
fn perform_nsv_compensation(
    sigma_est: &mut SigmaEstimates,
    conv_res: &mut [f64],
    step_half_size: usize,
) {
    if let Some(ref var_mask) = sigma_est.var_mask {
        let conv_res_tmp = pulse_convolution(var_mask, 2 * step_half_size);
        let conv_mask: Vec<f64> = conv_res_tmp
            .iter()
            .map(|x| if *x >= 0.5 { 1.0 } else { 0.0 })
            .collect();

        let mut sigma_est_vec: Vec<f64> = Vec::with_capacity(conv_mask.len());
        sigma_est_vec.extend(
            conv_mask
                .iter()
                .map(|x| (1.0 - *x) * sigma_est.peak_est + *x * sigma_est.noise_est),
        );
        for (x, y) in conv_res.iter_mut().zip(sigma_est_vec.iter()) {
            *x /= *y;
        }
        sigma_est.peak_est = 1.0;
    }
}

/// Implementation of Haar segmentation.
pub fn segment_haar_seg(
    vals: &[f64],
    quals: Option<&[f64]>,
    raw_vals: Option<&[f64]>,
    breaks_fdr_q: f64,
    haar_start_level: usize,
    haar_last_level: usize,
    logger: &Logger,
) -> Vec<usize> {
    let haar_end_level = haar_last_level + 1;

    // Short-circuit in the case of empty input.
    if vals.is_empty() {
        return vec![0, 1];
    }

    // Estimate sigma and compute the data necessary for non-stationary variance compenstation.
    let mut sigma_est = estimate_sigmas(vals, raw_vals, logger);

    // Collect breakpoints as positions in `vals`.
    let mut breakpoints: Vec<usize> = Vec::new();
    for level in haar_start_level..haar_end_level {
        let step_half_size = 2_usize.pow(level as u32);

        // Perform Haar convolution and find local peaks.
        let mut conv_res = haar_convolution(vals, quals, step_half_size, logger);
        let local_peak_positions = find_local_peaks(&conv_res);
        debug!(
            logger,
            "Found {} peaks at level {}",
            local_peak_positions.len(),
            level
        );
        if local_peak_positions.len() < 50 {
            debug!(logger, "== {:?}", local_peak_positions);
        }

        // Perform compensation of non-stationary variance.
        perform_nsv_compensation(&mut sigma_est, &mut conv_res, step_half_size);

        // Compute FDR threshold.
        let local_peaks: Vec<f64> = local_peak_positions
            .iter()
            .map(|pos| conv_res[*pos])
            .collect();
        let t = fdr_thresh(&local_peaks, breaks_fdr_q, sigma_est.peak_est, logger);

        // Keep only the peak values where the signal amplitude is large enough.
        let local_peaks_conv_res: Vec<f64> = local_peak_positions
            .iter()
            .map(|pos| conv_res[*pos])
            .collect();
        let addon_peaks: Vec<usize> = local_peaks_conv_res
            .iter()
            .enumerate()
            .filter(|(_pos, val)| val.abs() >= t)
            .map(|(pos, _val)| local_peak_positions[pos])
            .collect();

        if breakpoints.len() < 50 && addon_peaks.len() < 50 {
            debug!(
                logger,
                "breakpoints = {:?}, addon_peaks = {:?}", breakpoints, addon_peaks
            );
        }
        breakpoints = merge_breakpoints(&breakpoints, &addon_peaks, 2_usize.pow(level as u32 - 1));
        if breakpoints.len() < 50 {
            debug!(logger, "==merge_breakpoints==> {:?}", breakpoints);
        }
    }

    // Add first and last element and return result breakpoints vector.
    if *breakpoints.first().unwrap_or(&1_usize) != 0 {
        breakpoints.insert(0, 0);
    }
    if *breakpoints.last().unwrap_or(&1_usize) != vals.len() {
        breakpoints.push(vals.len());
    }
    breakpoints
}

#[cfg(test)]
mod tests {
    use super::*;

    extern crate slog;
    extern crate slog_async;
    extern crate slog_term;
    use slog::Drain;
    use slog::Logger;

    fn _logger() -> Logger {
        let decorator = slog_term::TermDecorator::new().build();
        let drain = slog_term::FullFormat::new(decorator).build().fuse();
        let drain = slog_async::Async::new(drain).build().fuse();
        slog::Logger::root(drain, o!())
    }

    #[test]
    fn call_haar_seg() {
        let logger = _logger();

        let mut vals: Vec<f64> = vec![0.0; 50];
        for i in 20..30 {
            vals[i] = 2.0;
        }
        let breakpoints = segment_haar_seg(&vals, None, None, 1e-7, 1, 5, &logger);
        assert_eq!(breakpoints, vec![0, 20, 30, 50]);
    }
}
