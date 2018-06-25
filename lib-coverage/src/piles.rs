// Pile-collecting related.

use std::cmp::max;

use rust_htslib::bam::{self, Read as BamRead};

use slog::Logger;

use super::errors::*;

/// Configuration for pile collection for black listing.
#[derive(Clone, Debug)]
struct PileCollectorOptions {
    /// Join piles less than this threshold apart.
    pile_max_gap: u32,
    /// Ignore reads with MAPQ less than this value.
    min_mapq: u8,
    /// Number of threads to use for BAM I/O.
    io_threads: u32,
}

/// Collection of read piles for black listing.
#[derive(Debug)]
pub struct PileCollector {
    /// Indexed BAM reader to use for reading the alignments.
    bam_reader: bam::IndexedReader,
    /// Configuration of the collector.
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
                !record.is_secondary()
                    && !record.is_duplicate()
                    && !record.is_supplementary()
                    && !record.is_duplicate()
                    && !record.is_quality_check_failed()
                    && (record.mapq() >= options.min_mapq)
                    && (!record.is_paired() || record.is_proper_pair())
            })
            .count();

        result[pos] = depth as u32;
    }

    result
}

/// Information for one pile.
#[derive(Debug)]
pub struct PileInfo {
    /// Start position.
    pub start: u32,
    /// End position.
    pub end: u32,
    /// Maximal coverage depth.
    pub size: u32,
}

/// Collect the piles from depth vector.
fn piles_from_depths(depths: &Vec<u32>, options: &PileCollectorOptions) -> Vec<PileInfo> {
    // Intervals for the result.
    let mut result = Vec::new();

    // Iterate over pileup depths.
    let mut prev: Option<PileInfo> = Option::None;
    for (pos, depth) in depths.iter().enumerate() {
        let pos = pos as u32;
        let depth = *depth;

        // Either extend previous interval or start a new one, adding previous one to `result`.
        if depth > 0 {
            prev = match prev {
                Some(PileInfo { start, end, size }) => {
                    if end + options.pile_max_gap >= pos {
                        // Extend previous interval.
                        Some(PileInfo {
                            start,
                            end: pos + 1,
                            size: max(depth, size),
                        })
                    } else {
                        // Store old interval and start new one.
                        if size > 0 {
                            result.push(PileInfo { start, end, size });
                        }
                        Some(PileInfo {
                            start: pos,
                            end: pos + 1,
                            size: depth,
                        })
                    }
                }
                None => Some(PileInfo {
                    start: pos,
                    end: pos + 1,
                    size: depth,
                }),
            }
        }
    }

    // Handle the last interval, if any.
    if let Some(PileInfo { start, end, size }) = prev {
        if size > 0 {
            result.push(PileInfo { start, end, size });
        }
    }

    result
}

impl PileCollector {
    /// Construct new `PileCollector`.
    pub fn new(bam_path: &str, pile_max_gap: u32, min_mapq: u8, io_threads: u32) -> Result<Self> {
        Ok(PileCollector {
            bam_reader: {
                let mut reader =
                    bam::IndexedReader::from_path(&bam_path).chain_err(|| "Invalid_path")?;
                reader
                    .set_threads(io_threads as usize)
                    .chain_err(|| "Could not set threads for reading")?;
                reader
            },
            options: PileCollectorOptions {
                pile_max_gap,
                min_mapq,
                io_threads,
            },
        })
    }

    /// Collect piles from coverage.
    pub fn collect_piles(&mut self, chrom: &str, logger: &Logger) -> Vec<PileInfo> {
        // Fetch whole contig.
        let tid = self.bam_reader.header().tid(chrom.as_bytes()).unwrap();
        let end = self.bam_reader.header().target_len(tid).unwrap() as usize;
        self.bam_reader
            .fetch(tid, 0, end as u32)
            .expect("Could not fetch contig");

        // Compute pile information objects.
        let pile_infos = piles_from_depths(
            &compute_depths(&mut self.bam_reader, &self.options, end),
            &self.options,
        );
        debug!(logger, "#intervals = {}", pile_infos.len());
        // debug!(logger, "{:?}", &pile_infos[1..10]);

        pile_infos
    }
}
