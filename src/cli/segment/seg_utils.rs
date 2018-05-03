/// Various utility code for segmentation.
use cli::shared::stats::Stats;

use rust_segment::haarseg::{HaarSegResult, HaarSegment};

use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::usize;

/// Element for the max-heap used in `refine_seg_remove_breakpoints`.
#[derive(Debug)]
struct HeapElem {
    /// Value before merging.
    old_val: f64,
    /// Value after changing.
    new_val: f64,
    /// The index of the breakpoint represented.
    breakpoint: usize,
}

impl HeapElem {
    /// Return weight used in sorting heap.
    fn weight(&self) -> f64 {
        (self.old_val as f64) / (self.new_val as f64)
    }
}

impl Ord for HeapElem {
    fn cmp(&self, other: &HeapElem) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl PartialOrd for HeapElem {
    fn partial_cmp(&self, other: &HeapElem) -> Option<Ordering> {
        self.weight().partial_cmp(&other.weight())
    }
}

impl PartialEq for HeapElem {
    fn eq(&self, other: &HeapElem) -> bool {
        (self.old_val, self.new_val) == (other.old_val, other.new_val)
    }
}

impl Eq for HeapElem {}

fn square(x: f64) -> f64 {
    x * x
}

/// Sentinel value for invalid breakpoint idx.
const INVALID: usize = usize::MAX;

/// Return `(old_value, new_value)` for errors after removing breakpoint.
fn compute_error(
    values: &[f64],
    error_squares: &[f64],
    bp_prev: usize,
    bp_this: usize,
    bp_next: usize,
) -> (f64, f64) {
    let left_len = (bp_this - bp_prev) as f64;
    let right_len = (bp_next - bp_this) as f64;

    // Compute old error.
    let old_left = (error_squares[bp_prev..bp_this].iter().sum::<f64>() / left_len).sqrt();
    let old_right = (error_squares[bp_this..bp_next].iter().sum::<f64>() / right_len).sqrt();
    let old_val = (old_left + old_right) / 2.0;

    // Compute potential change for selecting breakpoint at idx for deletion.
    let new_mean = values[bp_prev..bp_next].mean();
    let new_left = (values[bp_prev..bp_this]
        .iter()
        .map(|x| square(x - new_mean))
        .sum::<f64>() / left_len)
        .sqrt();
    let new_right = (values[bp_this..bp_next]
        .iter()
        .map(|x| square(x - new_mean))
        .sum::<f64>() / right_len)
        .sqrt();
    let new_val = (new_left + new_right) / 2.0;

    (old_val, new_val)
}

/// Refine `segmentation` of the values in `intensities` by removing breakpoints as long as the
/// ratios of the standard errors is below `thresh`.
pub fn refine_seg_remove_breakpoints(
    intensities: &[f64],
    segmentation: &HaarSegResult,
    thresh: f64,
) -> HaarSegResult {
    let segments = segmentation.segments.clone();
    let mut seg_values = segmentation.seg_values.clone();

    // Compute error squares, part of the standard error computation.
    let mut error_squares: Vec<f64> = intensities
        .iter()
        .zip(seg_values.iter())
        .map(|(x, y)| { println!("-> sq({} - {}) = {}", x, y, (square(x - y) * 1000.0).round() / 1000.0); square(x - y) })
        .collect::<Vec<f64>>();
    println!("Intensities: {:?}", intensities);
    println!("Error squares: {:?}", error_squares);
    // Build vector of all breakpoint positions.
    let mut breakpoints: Vec<usize> = vec![0];
    for ref seg in segments {
        breakpoints.push(seg.range.end);
    }

    // Mapping from input breakpoint index to output, INVALID if already removed.  This is the
    // simplest way to emulate a heap with "DECREATE-KEY()" using the Rust standard library.
    let mut bp_idx: Vec<usize> = (0..breakpoints.len()).collect();

    println!("Breakpoints: {:?}, bp_idx: {:?}", breakpoints, bp_idx);

    // Max-heap that stores candidate breakpoints to delete.
    let mut heap = BinaryHeap::new();

    // Initialize heap by computing the effect in error change.
    for idx in 1..(bp_idx.len() - 1) {
        let breakpoint = idx;

        let left_range = breakpoints[idx - 1]..breakpoints[idx];
        let right_range = breakpoints[idx]..breakpoints[idx + 1];
        let new_range = breakpoints[idx - 1]..breakpoints[idx + 1];
        println!("{:?} {:?} {:?}", left_range, right_range, new_range);

        let new_len = (right_range.len() + left_range.len()) as f64;

        // Compute old error.
        let old_left = error_squares[left_range.clone()]
            .iter()
            .sum::<f64>();
        let old_right = error_squares[right_range.clone()]
            .iter()
            .sum::<f64>();
        let old_val = ((old_left + old_right) / new_len).sqrt();

        // Compute potential change for selecting breakpoint at idx for deletion.
        let new_mean = intensities[new_range.clone()].mean();
        let new_left = intensities[left_range.clone()]
            .iter()
            .map(|&x| square(x - new_mean))
            .sum::<f64>();
        let new_right = intensities[right_range.clone()]
            .iter()
            .map(|&x| square(x - new_mean))
            .sum::<f64>();
        let new_val = ((new_left + new_right) / new_len).sqrt();

        let elem = HeapElem {
            breakpoint,
            old_val,
            new_val,
        };
        println!("Pushing {:?}", elem);
        heap.push(elem);
    }

    // Consider all elements in heap, if removing the currently best candidate is valid then
    // do it, otherwise skip.
    while !heap.is_empty() {
        // Pop best candidate from heap.
        let elem = heap.pop().unwrap();

        println!("Considering {:?} -> {}", elem, elem.weight());

        // Skip if breakpoint is already removed.
        if bp_idx[elem.breakpoint] == INVALID {
            println!("{} invalid", elem.breakpoint);
            continue;
        }
        // Break if error increase would be too large, no more suitable candidates.
        if elem.weight() > thresh {
            println!("Threshold reached, stopping");
            break;
        }

        // Store indices in breakpoint for the breakpoints left and right of to be removed
        // breakpoint for later
        let mut bp_left = elem.breakpoint - 1;
        while bp_left > 0 && bp_idx[bp_left] == INVALID {
            bp_left -= 1;
        }
        let mut bp_right = elem.breakpoint + 1;
        while bp_right + 1 < bp_idx.len() && bp_idx[bp_right] == INVALID {
            bp_right += 1;
        }

        // Remove breakpoint, shift indices of right breakpoints to the left
        println!("Removing breakpoint {:?}", elem);
        breakpoints.remove(bp_idx[elem.breakpoint]);
        let bp_idx_len = bp_idx.len();
        for idx in &mut bp_idx[(elem.breakpoint)..bp_idx_len] {
            if *idx != INVALID {
                *idx -= 1;
            }
        }
        bp_idx[elem.breakpoint] = INVALID;
        // Update mean and standard error of the new segment
        let new_mean =
            intensities[breakpoints[bp_idx[bp_left]]..breakpoints[bp_idx[bp_right]]].mean();
        for x in &mut seg_values[breakpoints[bp_idx[bp_left]]..breakpoints[bp_idx[bp_right]]] {
            *x = new_mean;
        }
        for i in (breakpoints[bp_idx[bp_left]])..(breakpoints[bp_idx[bp_right]]) {
            error_squares[i] = square(intensities[i] - segmentation.seg_values[i]);
        }

        // Push breakpoints bp_left and bp_right into heap, except if first/last
        if bp_idx[bp_left] > 0 && bp_idx[bp_left] + 1 < breakpoints.len() {
            let (old_val, new_val) = compute_error(
                intensities,
                &error_squares,
                breakpoints[bp_idx[bp_left] - 1],
                breakpoints[bp_idx[bp_left]],
                breakpoints[bp_idx[bp_left] + 1],
            );
            let elem = HeapElem {
                breakpoint: bp_left,
                old_val,
                new_val,
            };
            println!("Pushing {:?}", elem);
            heap.push(elem);
        }
        if bp_idx[bp_right] + 1 < breakpoints.len() {
            let (old_val, new_val) = compute_error(
                intensities,
                &error_squares,
                breakpoints[bp_idx[bp_right] - 1],
                breakpoints[bp_idx[bp_right]],
                breakpoints[bp_idx[bp_right] + 1],
            );
            let elem = HeapElem {
                breakpoint: bp_right,
                old_val,
                new_val,
            };
            println!("Pushing {:?}", elem);
            heap.push(elem);
        }
    }
    println!("Done, building result now!");

    HaarSegResult {
        segments: breakpoints[1..breakpoints.len()]
            .iter()
            .map(|i| HaarSegment {
                range: breakpoints[*i - 1]..breakpoints[*i],
                value: segmentation.seg_values[breakpoints[*i - 1]],
            })
            .collect(),
        seg_values: seg_values,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_segment::haarseg::seg_haar;

    #[test]
    fn test_simple_refinement() {
        // Force bad partitioning by splitting at `chrom_pos`.
        let intensities = vec![0.0, 0.1, 0.2, 0.8, 0.9, 1.1, 1.1, 0.7, 0.0, 0.0];
        let chrom_pos = vec![0..5, 5..9];
        let input = seg_haar(&intensities, None, None, &chrom_pos, 0.0001, 1, 5);
        println!("{:?}", input);

        let output = refine_seg_remove_breakpoints(&intensities, &input, 0.8);
        assert_eq!(1, 2);
    }
}
