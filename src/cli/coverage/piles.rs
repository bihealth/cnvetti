// Pile-collecting related.

use std::cmp::max;

use bio::data_structures::interval_tree;

use rust_htslib::bam::{self, Read as BamRead};


/// Configuration for pile collection for black listing.
#[derive(Clone, Debug)]
struct PileCollectorOptions {
    pile_min_depth: i32,
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


impl<'a> PileCollector<'a> {
    /// Construct new `PileCollector`.
    pub fn new(
        bam_reader: &'a mut bam::IndexedReader,
        pile_min_depth: i32,
        pile_max_gap: u32,
        min_mapq: u8,
        pile_mask_window_size: u32,
    ) -> PileCollector<'a> {
        PileCollector {
            bam_reader,
            options: PileCollectorOptions {
                pile_min_depth,
                pile_max_gap,
                pile_mask_window_size,
                min_mapq,
            },
        }
    }

    /// Collect the piles to black list.
    pub fn collect_piles(&mut self, chrom: &str) -> (u32, interval_tree::IntervalTree<u32, i32>) {
        // Fetch whole contig.
        let tid = self.bam_reader.header().tid(chrom.as_bytes()).unwrap();
        let end = self.bam_reader.header().target_len(tid).unwrap();
        self.bam_reader.fetch(tid, 0, end).expect(
            "Could not fetch contig",
        );

        // Intervals for the result.
        let mut intervals = Vec::new();

        // Clone options so there is no double-borrow on `self`.
        let options = self.options.clone();

        // Previous `(start, end, depth)`.
        let mut prev: Option<(u32, u32, i32)> = Option::None;

        // Iterate over all pileups
        for pileup in self.bam_reader.pileup() {
            let pileup = pileup.unwrap();
            let pos = pileup.pos();

            // Compute depth of "valid" reads (note that only single/first-read coverage)
            // is computed.
            let depth = pileup
                .alignments()
                .filter(|alignment| {
                    let record = alignment.record();
                    !record.is_duplicate() && !record.is_quality_check_failed() &&
                        (!record.is_paired() || record.is_first_in_template()) &&
                        (record.mapq() > options.min_mapq)
                })
                .count() as i32;

            // Either extend previous interval or start a new one, adding previous one to
            // `intervals`.
            prev = match prev {
                Some((start, end, max_depth)) => {
                    if end + options.pile_max_gap >= pos {
                        // Extend previous interval.
                        Some((start, pos, max(depth, max_depth)))
                    } else {
                        // Store old interval and start new one.
                        if max_depth >= options.pile_min_depth {
                            intervals.push((start, end, max_depth));
                        }
                        Some((pos, pos + 1, depth))
                    }
                }
                None => Some((pos, pos + 1, depth)),
            }
        }

        // Handle the last interval, if any.
        if let Some((start, end, max_depth)) = prev {
            if max_depth >= options.pile_min_depth {
                intervals.push((start, end, max_depth));
            }
        }

        // Finally build tree by rounding to `self.pile_mask_window_size`.
        let wlen = self.options.pile_mask_window_size;
        let mut tree = interval_tree::IntervalTree::new();

        let mut len_sum = 0;
        let mut prev_end = 0;
        for (start, end, max_depth) in intervals {
            // Compute bucket number of bucket size `wlen`.
            let start = max((start / wlen) * wlen, prev_end);
            let end = ((end + wlen - 1) / wlen) * wlen;

            // Insert one joint interval and count towards sum.
            tree.insert(start..end, max_depth);
            len_sum += end - start;

            prev_end = end;
        }

        (len_sum, tree)
    }
}
