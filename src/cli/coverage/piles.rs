// Pile-collecting related.

use std::cmp::max;

use bio::data_structures::interval_tree;

use rust_htslib::bam::{self, Read as BamRead};

use cli::shared::math;

use slog::Logger;

/// Configuration for pile collection for black listing.
#[derive(Clone, Debug)]
struct PileCollectorOptions {
    pile_depth_percentile: f64,
    pile_max_gap: u32,
    pile_mask_window_size: u32,
    min_mapq: u8,
}

/// Collection of read piles for black listing.
#[derive(Debug)]
pub struct PileCollector<'a> {
    bam_reader: &'a mut bam::IndexedReader,
    options: PileCollectorOptions,
}

/// Compute depths of the selected sections in reader.
fn compute_depths(
    bam_reader: &mut bam::IndexedReader,
    options: &PileCollectorOptions,
    contig_length: usize,
) -> Vec<u32> {
    let mut result: Vec<u32> = vec![0; contig_length];

    // Iterate over all pileups
    for pileup in bam_reader.pileup() {
        let pileup = pileup.unwrap();
        let pos = pileup.pos() as usize;

        // Compute depth of "valid" reads (note that only single/first-read coverage)
        // is computed.
        let depth = pileup
            .alignments()
            .filter(|alignment| {
                let record = alignment.record();
                !record.is_duplicate() && !record.is_supplementary() && !record.is_duplicate()
                    && !record.is_quality_check_failed()
                    && (record.mapq() >= options.min_mapq)
            })
            .count();

        result[pos] = depth as u32;
    }

    result
}

/// Collect the piles.
fn collect_piles(depths: &Vec<u32>, options: &PileCollectorOptions) -> Vec<(u32, u32, u32)> {
    // Intervals for the result.
    let mut result = Vec::new();

    // Iterate over pileup depths.
    let mut prev: Option<(u32, u32, u32)> = Option::None; // (start, end, depth)
    for (pos, depth) in depths.iter().enumerate() {
        let pos = pos as u32;
        let depth = *depth;

        // Either extend previous interval or start a new one, adding previous one to `result`.
        if depth > 0 {
            prev = match prev {
                Some((start, end, max_depth)) => {
                    if end + options.pile_max_gap >= pos {
                        // Extend previous interval.
                        Some((start, pos + 1, max(depth, max_depth)))
                    } else {
                        // Store old interval and start new one.
                        if max_depth > 0 {
                            result.push((start, end, max_depth));
                        }
                        Some((pos, pos + 1, depth))
                    }
                }
                None => Some((pos, pos + 1, depth)),
            }
        }
    }

    // Handle the last interval, if any.
    if let Some((start, end, max_depth)) = prev {
        if max_depth > 0 {
            result.push((start, end, max_depth));
        }
    }

    result
}

impl<'a> PileCollector<'a> {
    /// Construct new `PileCollector`.
    pub fn new(
        bam_reader: &'a mut bam::IndexedReader,
        pile_depth_percentile: f64,
        pile_max_gap: u32,
        min_mapq: u8,
        pile_mask_window_size: u32,
    ) -> PileCollector<'a> {
        PileCollector {
            bam_reader,
            options: PileCollectorOptions {
                pile_depth_percentile,
                pile_max_gap,
                pile_mask_window_size,
                min_mapq,
            },
        }
    }

    /// Collect the piles to black list.
    pub fn collect_piles(
        &mut self,
        chrom: &str,
        logger: &Logger,
        depth_threshold: Option<usize>,
    ) -> (usize, u32, interval_tree::IntervalTree<u32, u32>) {
        // Fetch whole contig.
        let tid = self.bam_reader.header().tid(chrom.as_bytes()).unwrap();
        let end = self.bam_reader.header().target_len(tid).unwrap() as usize;
        self.bam_reader
            .fetch(tid, 0, end as u32)
            .expect("Could not fetch contig");

        // Clone options so there is no double-borrow on `self`.
        let options = self.options.clone();

        // Compute coverage depths, properly filtered.
        let depths = compute_depths(&mut self.bam_reader, &self.options, end);
        trace!(logger, "#depths = {}", depths.len());

        // Compute piles.
        let intervals = collect_piles(&depths, &self.options);
        trace!(logger, "#intervals = {}", intervals.len());

        // Collect pile depths for computing percentile.
        let depths: Vec<f64> = intervals
            .iter()
            .map(|(_start, _end, depth)| *depth as f64)
            .collect();

        // Compute pile threshold.
        let depth_threshold = match depth_threshold {
            Some(depth_threshold) => depth_threshold,
            None => {
                let depth_threshold = math::percentile(&depths, options.pile_depth_percentile);
                debug!(
                    logger,
                    "Setting threshold to percentile {} value, is = {}",
                    options.pile_depth_percentile,
                    depth_threshold
                );
                depth_threshold as usize
            }
        };

        if chrom == "1" {
            use std::fs::File;
            use std::io::Write;
            let mut f = File::create("blocked.bed").unwrap();
            for (start, end, depth) in &intervals {
                if (*depth as usize) > depth_threshold {
                    f.write_all(format!("{}\t{}\t{}\t{}\n", chrom, start, end, depth).as_bytes()).unwrap();
                }
            }
        }

        // Build resulting interval tree.
        let mut tree = interval_tree::IntervalTree::new();
        let wlen = self.options.pile_mask_window_size;

        let mut len_sum = 0;
        let mut prev_end = 0;
        for (start, end, depth) in intervals {
            // Compute bucket number of bucket size `wlen`.
            let start = max((start / wlen) * wlen, prev_end);
            let end = ((end + wlen - 1) / wlen) * wlen;

            // Insert one joint interval and count towards sum.
            if (depth as usize) > depth_threshold {
                tree.insert(start..end, depth);
            }
            len_sum += end - start;

            prev_end = end;
        }

        (depth_threshold, len_sum, tree)
    }
}
