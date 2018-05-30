use std::cmp::{max, min};
use std::collections::HashMap;
use std::ops::Range;

use cli::coverage::options::*;

use bio::data_structures::interval_tree;

use rust_htslib::bam;
use rust_htslib::prelude::*;

use statrs::statistics::Statistics;

// TODO: make this configurable in options.
/// Largest fraction of pile-masked windows before ignoring.
const MAX_MS: f64 = 0.5;

/// Struct with common information for aggregator.
#[derive(Debug)]
pub struct BaseAggregator {
    /// Configuration.
    options: Options,
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

    /// Names of the character fields that will be written out.
    fn character_field_names(&self) -> Vec<String>;

    /// Names of the integer fields that will be written out.
    fn integer_field_names(&self) -> Vec<String>;

    /// Names of the float fields that will be written out.
    fn float_field_names(&self) -> Vec<String>;

    /// Return character values.
    fn character_values(&self, window_id: u32) -> HashMap<String, String>;

    /// Return integer values.
    fn integer_values(&self, window_id: u32) -> HashMap<String, i32>;

    /// Return float values.
    fn float_values(&self, window_id: u32) -> HashMap<String, f32>;

    /// Number of processed records.
    fn num_processed(&self) -> u32;

    /// Number of skipped records.
    fn num_skipped(&self) -> u32;

    /// Whether the window is masked for the sample.
    fn is_masked(&self, window_id: u32) -> bool;
}

/// Struct for aggregating as alignment counts.
#[derive(Debug)]
pub struct CountAlignmentsAggregator<'a> {
    /// Common information from all aggregators.
    base: BaseAggregator,

    // TODO: lifetime not necessary?
    /// The `IntervalTree` with black-listed intervals.
    tree: Option<&'a interval_tree::IntervalTree<u32, u32>>,

    /// Per-sample read counts (each bin is a `u32`).
    counters: Vec<u32>,

    // Optional BAM record writer.
    out_bam: Option<&'a mut bam::Writer>,
}

impl<'a> CountAlignmentsAggregator<'a> {
    /// Construct new aggregator with the given BAM `header`.  This information is
    /// necessary to appropriately allocate buffers for all samples in the header.
    pub fn new(
        tree: Option<&'a interval_tree::IntervalTree<u32, u32>>,
        options: Options,
        contig_length: usize,
        out_bam: Option<&'a mut bam::Writer>,
    ) -> CountAlignmentsAggregator<'a> {
        CountAlignmentsAggregator {
            base: BaseAggregator {
                options,
                contig_length,
                num_processed: 0,
                num_skipped: 0,
            },
            tree: tree,
            counters: vec![0; contig_length as usize],
            out_bam,
        }
    }
}

impl<'a> CountAlignmentsAggregator<'a> {
    fn put_bam_record(&mut self, record: &bam::Record) {
        if !self.skip_mapq(record) && !self.skip_flags(record) && !self.skip_discordant(record)
            && !self.skip_clipping(record) && !self.skip_paired_and_all_but_first(record)
        {
            self.base.num_processed += 1;

            let pos = (record.pos() as u32)..((record.pos() + 1) as u32);
            let window_length = self.base.options.window_length as usize;
            let bin = record.pos() as usize / window_length;
            if self.tree.is_some() && self.tree.unwrap().find(pos).next().is_none() {
                self.out_bam.as_mut().map(|ref mut out_bam| {
                    out_bam.write(record).expect("Could not write BAM record.")
                });
                self.counters[bin] += 1;
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
        record.is_secondary() || record.is_supplementary() || record.is_duplicate()
            || record.is_quality_check_failed()
    }

    /// Skip `record` because of discordant alignment?
    fn skip_discordant(&self, record: &bam::Record) -> bool {
        if !self.base.options.skip_discordant {
            false // skipping based on discordant has to be enabled
        } else if !record.is_paired() {
            false // unpaired cannot be discordant
        } else {
            // record.tid() != record.mtid() || record.is_reverse() == record.is_mate_reverse()
            //     || record.is_unmapped() || record.is_mate_unmapped() ||
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

    // Skip `record` because it is paired and not the first fragment?
    fn skip_paired_and_all_but_first(&self, record: &bam::Record) -> bool {
        record.is_paired() && !record.is_first_in_template()
    }
}

impl<'a> BamRecordAggregator for CountAlignmentsAggregator<'a> {
    /// Put all `fetch()`ed records from `reader` into the aggregator.
    fn put_fetched_records(&mut self, reader: &mut bam::IndexedReader) {
        let mut record = bam::Record::new();
        while reader.read(&mut record).is_ok() {
            self.put_bam_record(&record);
        }
    }

    fn character_field_names(&self) -> Vec<String> {
        vec![String::from("MP")]
    }

    fn integer_field_names(&self) -> Vec<String> {
        Vec::new()
    }

    fn float_field_names(&self) -> Vec<String> {
        vec![String::from("COV")]
    }

    fn integer_values(&self, _window_id: u32) -> HashMap<String, i32> {
        HashMap::new()
    }

    fn character_values(&self, window_id: u32) -> HashMap<String, String> {
        let mut result = HashMap::new();

        let window = Range {
            start: (window_id * self.base.options.window_length as u32) as u32,
            end: ((window_id + 1) * self.base.options.window_length as u32) as u32,
        };

        let ms = match self.tree {
            Some(ref tree) => tree.find(window.clone())
                .map(|entry| {
                    let end = min(window.end, entry.interval().end) as i32;
                    let start = max(window.start, entry.interval().start) as i32;
                    assert!(end >= start, "Should overlap, come from tree query");
                    end - start
                })
                .sum::<i32>(),
            None => 0_i32,
        };

        // Ratio for scaling up to the read count.
        let len = (window.end - window.start) as i32;
        let ratio = (len - ms) as f64 / len as f64;
        result.insert(
            String::from("MP"),
            if ratio < MAX_MS {
                "Y".to_string()
            } else {
                "N".to_string()
            },
        );

        result
    }

    fn float_values(&self, window_id: u32) -> HashMap<String, f32> {
        let window = Range {
            start: (window_id * self.base.options.window_length as u32) as u32,
            end: ((window_id + 1) * self.base.options.window_length as u32) as u32,
        };

        // Pile-masked bases.
        let ms = match self.tree {
            Some(ref tree) => tree.find(window.clone())
                .map(|entry| {
                    let end = min(window.end, entry.interval().end) as i32;
                    let start = max(window.start, entry.interval().start) as i32;
                    assert!(end >= start, "Should overlap, come from tree query");
                    end - start
                })
                .sum::<i32>(),
            None => 0_i32,
        };

        // Ratio for scaling up to the read count.
        let len = (window.end - window.start) as i32;
        let ratio = (len - ms) as f64 / len as f64;

        let mut result = HashMap::new();

        result.insert(
            String::from("COV"),
            if ratio > 0.001 {
                // guard against NaN...
                (self.counters[window_id as usize] as f64 / ratio) as f32
            } else {
                0_f32
            },
        );

        result
    }

    fn num_processed(&self) -> u32 {
        self.base.num_processed
    }

    fn num_skipped(&self) -> u32 {
        self.base.num_skipped
    }

    fn is_masked(&self, window_id: u32) -> bool {
        let window = Range {
            start: (window_id * self.base.options.window_length as u32) as u32,
            end: ((window_id + 1) * self.base.options.window_length as u32) as u32,
        };
        let ms = match self.tree {
            Some(ref tree) => tree.find(window.clone())
                .map(|entry| {
                    let end = min(window.end, entry.interval().end) as i32;
                    let start = max(window.start, entry.interval().start) as i32;
                    assert!(end >= start, "Should overlap, come from tree query");
                    end - start
                })
                .sum::<i32>(),
            None => 0_i32,
        };

        // Ratio for scaling up to the read count.
        let len = (window.end - window.start) as i32;
        let ratio = (len - ms) as f64 / len as f64;

        ratio < MAX_MS
    }
}

// Bin for coverage aggregation.
#[derive(Debug, Clone)]
pub struct CoverageBin {
    /// Mean coverage.
    pub cov_mean: f32,
    /// Coverage standard deviation.
    pub cov_stddev: f32,
}

impl CoverageBin {
    fn new() -> Self {
        CoverageBin {
            cov_mean: 0_f32,
            cov_stddev: 0_f32,
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
}

impl CoverageAggregator {
    pub fn new(options: Options, contig_length: usize) -> CoverageAggregator {
        let window_length = options.window_length as usize;
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
        }
    }

    fn push_window(&mut self, window_id: usize) {
        let depths: Vec<f64> = self.depths.iter().map(|x| *x as f64).collect();
        self.coverage[window_id] = CoverageBin {
            cov_mean: (&depths).mean() as f32,
            cov_stddev: (&depths).std_dev() as f32,
        };
        let window_length = self.base.options.window_length as usize;
        self.depths = vec![0; window_length];
    }
}

impl BamRecordAggregator for CoverageAggregator {
    /// Put all `fetch()`ed records from `reader` into the aggregator.
    fn put_fetched_records(&mut self, reader: &mut bam::IndexedReader) {
        let window_length = self.base.options.window_length as usize;

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
            self.depths[pos % window_length] = pileup
                .alignments()
                .filter(|alignment| {
                    let record = alignment.record();
                    !record.is_secondary() && !record.is_duplicate() && !record.is_supplementary()
                        && !record.is_duplicate()
                        && !record.is_quality_check_failed()
                        && (record.mapq() >= self.base.options.min_mapq)
                        && (!record.is_paired() || record.is_proper_pair())
                })
                .count();
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

    fn character_field_names(&self) -> Vec<String> {
        Vec::new()
    }

    fn integer_field_names(&self) -> Vec<String> {
        Vec::new()
    }

    fn float_field_names(&self) -> Vec<String> {
        vec![String::from("COV"), String::from("COVSD")]
    }

    fn character_values(&self, _window_id: u32) -> HashMap<String, String> {
        HashMap::new()
    }

    fn integer_values(&self, _window_id: u32) -> HashMap<String, i32> {
        HashMap::new()
    }

    fn float_values(&self, window_id: u32) -> HashMap<String, f32> {
        let window_id = window_id as usize;
        let mut result = HashMap::new();

        result.insert(String::from("COV"), self.coverage[window_id].cov_mean);
        result.insert(String::from("COVSD"), self.coverage[window_id].cov_stddev);

        result
    }

    fn num_processed(&self) -> u32 {
        self.base.num_processed
    }

    fn num_skipped(&self) -> u32 {
        self.base.num_skipped
    }

    /// Whether the window is masked for the sample.
    fn is_masked(&self, _window_id: u32) -> bool {
        false // no pile-based masking for coverage
    }
}
