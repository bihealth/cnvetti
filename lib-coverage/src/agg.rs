/// Code for aggregating BAM records in windows and target regions.
use std::cmp::{max, min};
use std::ops::Range;

use options::CoverageOptions;

use bio::data_structures::interval_tree;

use rust_htslib::bam;
use rust_htslib::bcf::record::Numeric;
use rust_htslib::prelude::*;

use lib_shared::regions::GenomeRegions;
use lib_shared::stats::Stats;

/// Struct for the result in one region (window or target).
#[derive(Debug)]
pub struct AggregationStats {
    /// Coverage in the window to use after the "coverage" step.
    pub cov: f32,
    /// Optional standard deviation of coverage in the window.
    pub cov_sd: Option<f32>,

    /// Start of the region on the contig.
    pub start: usize,
    /// End of the region on the contig.
    pub end: usize,
    /// Fraction of region that was removed (e.g., due to piles).
    pub frac_removed: Option<f32>,

    /// Mean MAPQ in region.
    pub mean_mapq: f32,
}

/// Struct with common information for aggregator.
#[derive(Debug)]
pub struct BaseAggregator {
    /// Configuration.
    options: CoverageOptions,
    /// Length of contig process
    contig_length: usize,

    /// Number of processed records.
    num_processed: u32,
    /// Number of skipped records.
    num_skipped: u32,
}

/// Trait for alignment aggregation from BAM files.
pub trait BamRecordAggregator {
    /// Put all `fetch()`ed records from `reader` into the aggregator.
    fn put_fetched_records(&mut self, reader: &mut bam::IndexedReader);

    /// Return statistics for the given window.
    fn get_stats(&self, window_id: usize) -> AggregationStats;

    /// Return number of windows/targets.
    fn num_regions(&self) -> usize;

    /// Number of processed records.
    fn num_processed(&self) -> u32;

    /// Number of skipped records.
    fn num_skipped(&self) -> u32;
}

/// Struct for aggregating fragment counts in a genome-wide fashion.
#[derive(Debug)]
pub struct FragmentsGenomeWideAggregator<'a> {
    /// Common information from all aggregators.
    base: BaseAggregator,

    // TODO: lifetime not necessary?
    /// The `IntervalTree` with black-listed intervals.
    tree: Option<&'a interval_tree::IntervalTree<u32, u32>>,

    /// Per-window fragment counts (each bin is a `u32`).
    counters: Vec<u32>,
    /// Sum of MAPQ values.
    mapq_sums: Vec<u64>,
    /// Window length.
    window_length: usize,
    // TODO: add back writing of output BAM file
    // // Optional BAM record writer.
    // out_bam: Option<&'a mut bam::Writer>,
}

impl<'a> FragmentsGenomeWideAggregator<'a> {
    /// Construct new aggregator with the given BAM `header`.  This information is
    /// necessary to appropriately allocate buffers for all samples in the header.
    pub fn new(
        tree: Option<&'a interval_tree::IntervalTree<u32, u32>>,
        options: CoverageOptions,
        contig_length: usize,
        // out_bam: Option<&'a mut bam::Writer>,
    ) -> FragmentsGenomeWideAggregator<'a> {
        let window_length = options.window_length.unwrap() as usize;
        let num_bins = (contig_length + window_length - 1) / window_length;
        FragmentsGenomeWideAggregator {
            base: BaseAggregator {
                options,
                contig_length,
                num_processed: 0,
                num_skipped: 0,
            },
            tree: tree,
            counters: vec![0; contig_length as usize],
            mapq_sums: vec![0; num_bins],
            window_length: window_length,
            // out_bam,
        }
    }
}

impl<'a> FragmentsGenomeWideAggregator<'a> {
    fn put_bam_record(&mut self, record: &bam::Record) {
        if !self.skip_mapq(record)
            && !self.skip_flags(record)
            && !self.skip_discordant(record)
            && !self.skip_clipping(record)
            && !self.skip_paired_and_all_but_leftmost(record)
        {
            self.base.num_processed += 1;

            let fragment_center = if record.is_paired() {
                record.pos() + record.insert_size() / 2
            } else {
                record.pos()
                    + record
                        .cigar()
                        .end_pos()
                        .expect("Problem interpreting CIGAR string")
            } as u32;

            let pos = fragment_center..(fragment_center + 1);
            let window_length = self.base
                .options
                .window_length
                .expect("Window length must be set here") as usize;
            let bin = fragment_center as usize / window_length;
            if !self.tree.is_some() || self.tree.unwrap().find(pos).next().is_none() {
                self.counters[bin] += 1;
                self.mapq_sums[bin] += record.mapq() as u64;
            } else {
                self.base.num_skipped += 1;
            }
        }
    }

    // Skip `record` based on `MAPQ`?
    fn skip_mapq(&self, record: &bam::Record) -> bool {
        record.mapq() < self.base.options.min_mapq
    }

    // Skip `record` because of flags.
    fn skip_flags(&self, record: &bam::Record) -> bool {
        record.is_secondary()
            || record.is_supplementary()
            || record.is_duplicate()
            || record.is_quality_check_failed()
    }

    /// Skip `record` because of discordant alignment?
    fn skip_discordant(&self, record: &bam::Record) -> bool {
        if !record.is_paired() {
            false // unpaired cannot be discordant
        } else {
            !record.is_proper_pair()
        }
    }

    // Skip `record` because of too much clipping?
    fn skip_clipping(&self, record: &bam::Record) -> bool {
        use rust_htslib::bam::record::Cigar::*;

        let mut num_clipped = 0;
        let mut num_unclipped = 0;
        for c in record.cigar().iter() {
            match c {
                Match(num) | Ins(num) | Equal(num) | Diff(num) => {
                    num_unclipped += num;
                }
                SoftClip(num) | HardClip(num) => {
                    num_clipped += num;
                }
                Del(_) | RefSkip(_) | Pad(_) => (),
            }
        }

        (num_unclipped as f32) / ((num_clipped + num_unclipped) as f32)
            < self.base.options.min_unclipped
    }

    // Skip `record` because it is paired and not the leftmost fragment (or discordant).
    fn skip_paired_and_all_but_leftmost(&self, record: &bam::Record) -> bool {
        record.is_paired() && (record.insert_size() <= 0)
    }

    /// Return number of pile-masked bases in the window.
    fn get_masked_count(&self, window_id: usize) -> u32 {
        let window_length = self.base.options.window_length.unwrap();
        let window = Range {
            start: (window_id * window_length) as u32,
            end: ((window_id + 1) * window_length) as u32,
        };
        match self.tree {
            Some(ref tree) => tree.find(window.clone())
                .map(|entry| {
                    let end = min(window.end, entry.interval().end) as i32;
                    let start = max(window.start, entry.interval().start) as i32;
                    assert!(end >= start, "Should overlap, come from tree query");
                    (end - start) as u32
                })
                .sum::<u32>(),
            None => 0_u32,
        }
    }
}

impl<'a> BamRecordAggregator for FragmentsGenomeWideAggregator<'a> {
    fn put_fetched_records(&mut self, reader: &mut bam::IndexedReader) {
        let mut record = bam::Record::new();
        while reader.read(&mut record).is_ok() {
            self.put_bam_record(&record);
        }
    }

    fn get_stats(&self, window_id: usize) -> AggregationStats {
        let window_length = self.base.options.window_length.unwrap();
        let masked_count = self.get_masked_count(window_id);

        AggregationStats {
            cov: self.counters[window_id as usize] as f32,
            cov_sd: None,

            start: window_id * window_length,
            end: min(self.base.contig_length, (window_id + 1) * window_length),
            frac_removed: Some((masked_count as f32) / (window_length as f32)),
            mean_mapq: if self.counters[window_id as usize] != 0 {
                ((self.mapq_sums[window_id as usize] as f64)
                    / (self.counters[window_id as usize] as f64)) as f32
            } else {
                f32::missing()
            },
        }
    }

    fn num_regions(&self) -> usize {
        (self.base.contig_length + self.window_length - 1) / self.window_length
    }

    fn num_processed(&self) -> u32 {
        self.base.num_processed
    }

    fn num_skipped(&self) -> u32 {
        self.base.num_skipped
    }
}

/// Struct for aggregating fragment counts in a genome-wide fashion.
#[derive(Debug)]
pub struct FragmentsTargetRegionsAggregator {
    /// Common information from all aggregators.
    base: BaseAggregator,

    /// The `IntervalTree` with target region intervals.
    tree: interval_tree::IntervalTree<u32, u32>,
    /// The target region intervals.
    target_regions: GenomeRegions,

    /// Per-target fragment counts (each bin is a `u32`).
    counters: Vec<u32>,
    /// Per-target sum of MAPQ values.
    mapq_sums: Vec<u64>,
}

impl FragmentsTargetRegionsAggregator {
    /// Construct new aggregator.
    pub fn new(
        tree: interval_tree::IntervalTree<u32, u32>,
        target_regions: GenomeRegions,
        options: CoverageOptions,
        contig_length: usize,
    ) -> FragmentsTargetRegionsAggregator {
        let num_targets = target_regions.regions.len();
        FragmentsTargetRegionsAggregator {
            base: BaseAggregator {
                options,
                contig_length,
                num_processed: 0,
                num_skipped: 0,
            },
            tree,
            target_regions,
            counters: vec![0; num_targets],
            mapq_sums: vec![0; num_targets],
        }
    }
}

impl FragmentsTargetRegionsAggregator {
    fn put_bam_record(&mut self, record: &bam::Record) {
        if !self.skip_mapq(record)
            && !self.skip_flags(record)
            && !self.skip_discordant(record)
            && !self.skip_clipping(record)
            && !self.skip_paired_and_all_but_leftmost(record)
        {
            self.base.num_processed += 1;

            // Get center of sequenced fragment.
            let begin_pos = record.pos() as u32;
            let end_pos = if record.is_paired() {
                record.pos() + record.insert_size()
            } else {
                record
                    .cigar()
                    .end_pos()
                    .expect("Problem interpreting CIGAR string")
            } as u32;
            let fragment_center = (begin_pos + end_pos) / 2;

            // Query for overlapping target region.
            let mut best = None;
            for it in self.tree.find(begin_pos..end_pos) {
                let itv = it.interval();
                let idx = *it.data();
                best = match best {
                    None => Some((idx, Self::compute_dist(&**itv, fragment_center))),
                    Some((best_idx, best_dist)) => {
                        let this_dist = Self::compute_dist(&**itv, fragment_center);
                        if this_dist < best_dist {
                            Some((idx, this_dist))
                        } else {
                            Some((best_idx, best_dist))
                        }
                    }
                };
            }
            if let Some(best) = best {
                self.counters[best.0 as usize] += 1;
                self.mapq_sums[best.0 as usize] += record.mapq() as u64;
            } else {
                self.base.num_skipped += 1;
            }
        }
    }

    /// Compute distance between point and interval.
    fn compute_dist(itv: &Range<u32>, pt: u32) -> u32 {
        let start = itv.start;
        let end = itv.end;

        if pt < start {
            start - pt
        } else if pt >= end {
            pt - end + 1
        } else {
            0
        }
    }
    // Skip `record` based on `MAPQ`?
    fn skip_mapq(&self, record: &bam::Record) -> bool {
        record.mapq() < self.base.options.min_mapq
    }

    // Skip `record` because of flags.
    fn skip_flags(&self, record: &bam::Record) -> bool {
        record.is_secondary()
            || record.is_supplementary()
            || record.is_duplicate()
            || record.is_quality_check_failed()
    }

    /// Skip `record` because of discordant alignment?
    fn skip_discordant(&self, record: &bam::Record) -> bool {
        if !record.is_paired() {
            false // unpaired cannot be discordant
        } else {
            !record.is_proper_pair()
        }
    }

    // Skip `record` because of too much clipping?
    fn skip_clipping(&self, record: &bam::Record) -> bool {
        use rust_htslib::bam::record::Cigar::*;

        let mut num_clipped = 0;
        let mut num_unclipped = 0;
        for c in record.cigar().iter() {
            match c {
                Match(num) | Ins(num) | Equal(num) | Diff(num) => {
                    num_unclipped += num;
                }
                SoftClip(num) | HardClip(num) => {
                    num_clipped += num;
                }
                Del(_) | RefSkip(_) | Pad(_) => (),
            }
        }

        (num_unclipped as f32) / ((num_clipped + num_unclipped) as f32)
            < self.base.options.min_unclipped
    }

    // Skip `record` because it is paired and not the leftmost fragment (or discordant).
    fn skip_paired_and_all_but_leftmost(&self, record: &bam::Record) -> bool {
        record.is_paired() && (record.insert_size() <= 0)
    }
}

impl BamRecordAggregator for FragmentsTargetRegionsAggregator {
    fn put_fetched_records(&mut self, reader: &mut bam::IndexedReader) {
        let mut record = bam::Record::new();
        while reader.read(&mut record).is_ok() {
            self.put_bam_record(&record);
        }
    }

    fn get_stats(&self, target_id: usize) -> AggregationStats {
        AggregationStats {
            cov: self.counters[target_id as usize] as f32,
            cov_sd: None,

            start: self.target_regions.regions[target_id].1 as usize,
            end: self.target_regions.regions[target_id].2 as usize,
            frac_removed: None,
            mean_mapq: if self.counters[target_id as usize] != 0 {
                ((self.mapq_sums[target_id as usize] as f64)
                    / (self.counters[target_id as usize] as f64)) as f32
            } else {
                f32::missing()
            },
        }
    }

    fn num_regions(&self) -> usize {
        self.target_regions.regions.len()
    }

    fn num_processed(&self) -> u32 {
        self.base.num_processed
    }

    fn num_skipped(&self) -> u32 {
        self.base.num_skipped
    }
}
// Bin for coverage aggregation.
#[derive(Debug, Clone)]
pub struct CoverageBin {
    /// Mean coverage.
    pub cov_mean: f32,
    /// Coverage standard deviation.
    pub cov_stddev: f32,
    /// Mean MAPQ, weighted by aligned bases.
    pub mapq_mean: f32,
}

impl CoverageBin {
    fn new() -> Self {
        CoverageBin {
            cov_mean: 0_f32,
            cov_stddev: 0_f32,
            mapq_mean: 0_f32,
        }
    }
}

// Struct for aggregating as coverage.
#[derive(Debug)]
pub struct CoverageAggregator {
    /// Common information from all aggregators.
    pub base: BaseAggregator,
    /// Number of bases for each base in the current bin for the one sample in the BAM file.
    pub coverage: Vec<CoverageBin>,

    /// Base-wise coverage information for the current window.
    pub depths: Vec<usize>,
    /// ID of the current window.
    pub window_id: Option<usize>,
    /// Length of the windows.
    pub window_length: usize,
    /// Sum of MAPQ values.
    pub mapqs: Vec<u8>,
}

impl CoverageAggregator {
    pub fn new(options: CoverageOptions, contig_length: usize) -> CoverageAggregator {
        let window_length = options.window_length.unwrap() as usize;
        let num_bins = (contig_length + window_length - 1) / window_length;
        CoverageAggregator {
            base: BaseAggregator {
                options,
                contig_length,
                num_processed: 0,
                num_skipped: 0,
            },
            coverage: vec![CoverageBin::new(); num_bins],
            depths: vec![0; window_length],
            window_id: None,
            mapqs: vec![0; window_length],
            window_length: window_length,
        }
    }

    fn push_window(&mut self, window_id: usize) {
        let depths: Vec<f64> = self.depths.iter().map(|x| *x as f64).collect();
        let mapqs: Vec<f64> = self.mapqs.iter().map(|x| *x as f64).collect();
        self.coverage[window_id] = CoverageBin {
            cov_mean: (&depths).mean() as f32,
            cov_stddev: (&depths).std_dev() as f32,
            mapq_mean: (&mapqs).mean() as f32,
        };
        let window_length = self.base.options.window_length.unwrap() as usize;
        self.depths = vec![0; window_length];
        self.mapqs = vec![0; window_length];
    }
}

impl BamRecordAggregator for CoverageAggregator {
    /// Put all `fetch()`ed records from `reader` into the aggregator.
    fn put_fetched_records(&mut self, reader: &mut bam::IndexedReader) {
        let window_length = self.base.options.window_length.unwrap() as usize;

        // Iterate over all pileups
        let mut prev_window_id = None;
        for pileup in reader.pileup() {
            let pileup = pileup.unwrap();
            let pos = pileup.pos() as usize;

            // On window change, push window to result.
            self.window_id = match (self.window_id, prev_window_id) {
                (Some(window_id), Some(my_prev_window_id)) => {
                    if window_id != my_prev_window_id {
                        self.push_window(my_prev_window_id);
                        prev_window_id = Some(window_id);
                    }
                    Some(pos / window_length)
                }
                (Some(window_id), None) => {
                    prev_window_id = Some(window_id);
                    Some(pos / window_length)
                }
                (None, _) => Some(pos / window_length as usize),
            };

            // Compute depth of "valid" reads (note that only single/first-read coverage)
            // is computed.
            let mapqs = pileup
                .alignments()
                .filter(|alignment| {
                    let record = alignment.record();
                    !record.is_secondary()
                        && !record.is_duplicate()
                        && !record.is_supplementary()
                        && !record.is_duplicate()
                        && !record.is_quality_check_failed()
                        && (record.mapq() >= self.base.options.min_mapq)
                        && (!record.is_paired() || record.is_proper_pair())
                })
                .map(|alignment| alignment.record().mapq())
                .collect::<Vec<u8>>();
            self.depths[pos % window_length] = mapqs.len();
            self.mapqs.extend(&mapqs);
        }

        if let Some(window_id) = self.window_id {
            match prev_window_id {
                Some(prev_window_id) => if prev_window_id != window_id {
                    self.push_window(window_id);
                },
                None => self.push_window(window_id),
            }
        }
    }

    fn get_stats(&self, window_id: usize) -> AggregationStats {
        let window_length = self.base.options.window_length.unwrap();
        AggregationStats {
            cov: self.coverage[window_id as usize].cov_mean,
            cov_sd: Some(self.coverage[window_id as usize].cov_stddev),

            start: window_id * window_length,
            end: min(self.base.contig_length, (window_id + 1) * window_length),
            frac_removed: None,
            mean_mapq: self.coverage[window_id as usize].mapq_mean,
        }
    }

    fn num_regions(&self) -> usize {
        (self.base.contig_length + self.window_length - 1) / self.window_length
    }

    fn num_processed(&self) -> u32 {
        self.base.num_processed
    }

    fn num_skipped(&self) -> u32 {
        self.base.num_skipped
    }
}
