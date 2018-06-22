// Pile-collecting related.

use std::cmp::max;
use std::fs;
use std::io::Write;

use bio::data_structures::interval_tree;

use rust_htslib::bam::{self, Read as BamRead};

use slog::Logger;

use lib_shared::stats::Stats;

use super::errors::*;

/// Configuration for pile collection for black listing.
#[derive(Clone, Debug)]
struct PileCollectorOptions {
    pile_depth_percentile: f64,
    pile_max_gap: u32,
    min_mapq: u8,
}

/// Collection of read piles for black listing.
#[derive(Debug)]
pub struct PileCollector<'a> {
    bam_reader: bam::IndexedReader,
    bed_file: Option<&'a mut fs::File>,
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
                !record.is_secondary() && !record.is_duplicate() && !record.is_supplementary()
                    && !record.is_duplicate() && !record.is_quality_check_failed()
                    && (record.mapq() >= options.min_mapq)
                    && (!record.is_paired() || record.is_proper_pair())
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
    let mut prev: Option<(u32, u32, u32)> = Option::None; // (start, end, base count)
    for (pos, depth) in depths.iter().enumerate() {
        let pos = pos as u32;
        let depth = *depth;

        // Either extend previous interval or start a new one, adding previous one to `result`.
        if depth > 0 {
            prev = match prev {
                Some((start, end, num_bases)) => {
                    if end + options.pile_max_gap >= pos {
                        // Extend previous interval.
                        Some((start, pos + 1, max(depth, depth + num_bases)))
                    } else {
                        // Store old interval and start new one.
                        if num_bases > 0 {
                            result.push((start, end, num_bases));
                        }
                        Some((pos, pos + 1, depth))
                    }
                }
                None => Some((pos, pos + 1, depth)),
            }
        }
    }

    // Handle the last interval, if any.
    if let Some((start, end, num_bases)) = prev {
        if num_bases > 0 {
            result.push((start, end, num_bases));
        }
    }

    result
}

impl<'a> PileCollector<'a> {
    /// Construct new `PileCollector`.
    pub fn new(
        bam_path: &str,
        bed_file: Option<&'a mut fs::File>,
        pile_depth_percentile: f64,
        pile_max_gap: u32,
        min_mapq: u8,
    ) -> Result<Self> {
        Ok(PileCollector {
            bam_reader: bam::IndexedReader::from_path(&bam_path).chain_err(|| "Invalid_path")?,
            bed_file,
            options: PileCollectorOptions {
                pile_depth_percentile,
                pile_max_gap,
                min_mapq,
            },
        })
    }

    /// Collect the piles to black list.
    pub fn collect_piles(
        &mut self,
        chrom: &str,
        logger: &Logger,
        num_bases_thresh: Option<usize>,
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

        // Compute piles [(start, end, bases in pile)]
        let intervals = collect_piles(&depths, &self.options);
        trace!(logger, "#intervals = {}", intervals.len());

        // Collect pile base counts for computing percentile.
        let depths: Vec<f64> = intervals
            .iter()
            .map(|(_start, _end, num_bases)| *num_bases as f64)
            .collect();

        // Compute pile threshold.
        let num_bases_thresh = match num_bases_thresh {
            Some(num_bases_thresh) => num_bases_thresh,
            None => {
                let num_bases_thresh = depths.percentile(options.pile_depth_percentile);
                debug!(
                    logger,
                    "Setting threshold to percentile {} value, is = {}",
                    options.pile_depth_percentile,
                    num_bases_thresh
                );
                num_bases_thresh as usize
            }
        };

        // Write out blocked intervals if output path given.
        if let Some(ref mut bed_file) = self.bed_file {
            for (start, end, num_bases) in &intervals {
                if (*num_bases as usize) > num_bases_thresh {
                    bed_file
                        .write_all(
                            format!("{}\t{}\t{}\t{}\n", chrom, start, end, num_bases).as_bytes(),
                        )
                        .expect("Could not write interval to BED file");
                }
            }
        }

        // Build resulting interval tree.
        let mut tree = interval_tree::IntervalTree::new();

        let mut len_sum = 0;
        let mut prev_end = 0;
        for (start, end, num_bases) in intervals {
            // Compute bucket number of bucket size `wlen`.
            let start = max(start, prev_end);
            let end = end;

            // Insert one joint interval and count towards sum.
            if (num_bases as usize) > num_bases_thresh {
                tree.insert(start..end, num_bases);
            }
            len_sum += end - start;

            prev_end = end;
        }

        (num_bases_thresh, len_sum, tree)
    }
}
