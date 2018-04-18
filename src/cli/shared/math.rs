// Stolen from rust libtest.

use std::cmp::Ordering::{self, Equal, Greater, Less};

fn local_cmp(x: f64, y: f64) -> Ordering {
    // arbitrarily decide that NaNs are larger than everything.
    if y.is_nan() {
        Less
    } else if x.is_nan() {
        Greater
    } else if x < y {
        Less
    } else if x == y {
        Equal
    } else {
        Greater
    }
}

fn local_sort(v: &mut [f64]) {
    v.sort_by(|x: &f64, y: &f64| local_cmp(*x, *y));
}

// Helper function: extract a value representing the `pct` percentile of a sorted sample-set, using
// linear interpolation. If samples are not sorted, return nonsensical value.
fn percentile_of_sorted(sorted_samples: &[f64], pct: f64) -> f64 {
    assert!(!sorted_samples.is_empty());
    if sorted_samples.len() == 1 {
        return sorted_samples[0];
    }
    let zero: f64 = 0.0;
    assert!(zero <= pct);
    let hundred = 100f64;
    assert!(pct <= hundred);
    if pct == hundred {
        return sorted_samples[sorted_samples.len() - 1];
    }
    let length = (sorted_samples.len() - 1) as f64;
    let rank = (pct / hundred) * length;
    let lrank = rank.floor();
    let d = rank - lrank;
    let n = lrank as usize;
    let lo = sorted_samples[n];
    let hi = sorted_samples[n + 1];
    lo + (hi - lo) * d
}

pub fn percentile(vals: &[f64], pct: f64) -> f64 {
    let mut tmp = vals.to_vec();
    local_sort(&mut tmp);
    percentile_of_sorted(&tmp, pct)
}

/// Compute median.
pub fn median(vals: &[f64]) -> f64 {
    percentile(vals, 50.0)
}

/// Compute median absolute deviation.
pub fn median_abs_dev(vals: &[f64]) -> f64 {
    let med = median(vals);
    let abs_devs: Vec<f64> = vals.iter().map(|&v| (med - v).abs()).collect();
    // This constant is derived by smarter statistics brains than me, but it is
    // consistent with how R and other packages treat the MAD.
    let number = 1.4826;
    median(&abs_devs) * number
}
