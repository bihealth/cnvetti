use std::cmp::{max, min};
use std::collections::HashMap;
use std::ops::Range;

use cli::coverage::options::*;

use bio::data_structures::interval_tree;

use rust_htslib::bam;

// TODO: make this configurable in options.
/// Largest fraction of pile-masked windows before ignoring.
const MAX_MS: f64 = 0.5;

/// Struct with common information for aggregator.
#[derive(Debug)]
pub struct BaseAggregator {
    /// Names of the samples to aggregate in this sample.
    samples: Vec<String>,
    /// Mapping from read group name to sample.
    read_group_to_sample: HashMap<String, String>,
    /// Mapping from sample name to sample index.
    sample_to_idx: HashMap<String, u32>,

    /// Configuration.
    options: Options,
    /// Length of contig process
    contig_length: u64,

    /// Number of processed records.
    num_processed: u32,
    /// Number of skipped records.
    num_skipped: u32,
}

/// Trait for alignment aggregation from BAM files.
pub trait BamRecordAggregator {
    /// Put the given record into the aggregation.
    fn put_bam_record(&mut self, record: &bam::Record);

    /// Names of the integer fields that will be written out.
    fn integer_field_names(&self) -> Vec<String>;

    /// Names of the float fields that will be written out.
    fn float_field_names(&self) -> Vec<String>;

    /// Return float values.
    fn integer_values(&self, window_id: u32) -> HashMap<String, i32>;

    /// Return float values.
    fn float_values(&self, window_id: u32) -> HashMap<String, f32>;

    /// Number of processed records.
    fn num_processed(&self) -> u32;

    /// Number of skipped records.
    fn num_skipped(&self) -> u32;

    /// Whether the window is masked for the sample.
    fn is_masked(&self, window_id: u32) -> bool;

    /// Coverage/count for the given window.
    fn get_coverage(&self, window_id: u32) -> i32;
}

/// Struct for aggregating as alignment counts.
#[derive(Debug)]
pub struct CountAlignmentsAggregator<'a> {
    /// Common information from all aggregators.
    base: BaseAggregator,

    // TODO: lifetime not necessary?
    /// The `IntervalTree` with black-listed intervals.
    tree: &'a interval_tree::IntervalTree<u32, u32>,

    /// Per-sample read counts (each bin is a `u32`).
    counters: Vec<u32>,
}

impl<'a> CountAlignmentsAggregator<'a> {
    /// Construct new aggregator with the given BAM `header`.  This information is
    /// necessary to appropriately allocate buffers for all samples in the header.
    pub fn new(
        tree: &'a interval_tree::IntervalTree<u32, u32>,
        samples: Vec<String>,
        read_group_to_sample: HashMap<String, String>,
        sample_to_idx: HashMap<String, u32>,
        options: Options,
        contig_length: u64,
    ) -> CountAlignmentsAggregator<'a> {
        CountAlignmentsAggregator {
            base: BaseAggregator {
                samples,
                read_group_to_sample,
                sample_to_idx,
                options,
                contig_length,
                num_processed: 0,
                num_skipped: 0,
            },
            tree: tree,
            counters: vec![0; contig_length as usize],
        }
    }
}

impl<'a> CountAlignmentsAggregator<'a> {
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
    fn put_bam_record(&mut self, record: &bam::Record) {
        use std::str;
        let log =
            (record.tid() == 0) && (record.pos() >= 62_160_000) && (record.pos() <= 62_180_000);
        // let log = false;
        if log {
            println!("Record {}", str::from_utf8(record.qname()).unwrap());
        }

        if !self.skip_mapq(record) && !self.skip_flags(record) && !self.skip_discordant(record)
            && !self.skip_clipping(record) && !self.skip_paired_and_all_but_first(record)
        {
            if log {
                println!("=> considering for counting");
            }
            self.base.num_processed += 1;

            let pos = (record.pos() as u32)..((record.pos() + 1) as u32);
            let window_length = self.base.options.window_length as usize;
            let bin = record.pos() as usize / window_length;
            if self.tree.find(pos).next().is_none() {
                if log {
                    // bin == 94_860/20 {
                    println!("-> Counting {}", str::from_utf8(record.qname()).unwrap());
                }

                self.counters[bin] += 1;
            } else {
                if log {
                    println!("  -> Skipping {}", str::from_utf8(record.qname()).unwrap());
                }
                self.base.num_skipped += 1;
            }
        } else {
            if log {
                println!("=> NOT considering for counting");
            }
        }
    }

    fn integer_field_names(&self) -> Vec<String> {
        vec![
            String::from("RRC"),
            String::from("RC"),
            String::from("MS"),
            String::from("MP"),
        ]
    }

    fn float_field_names(&self) -> Vec<String> {
        Vec::new()
    }

    fn integer_values(&self, window_id: u32) -> HashMap<String, i32> {
        let mut result = HashMap::new();

        let window = Range {
            start: (window_id * self.base.options.window_length as u32) as u32,
            end: ((window_id + 1) * self.base.options.window_length as u32) as u32,
        };

        // Pile-masked bases.
        let ms = self.tree
            .find(window.clone())
            .map(|entry| {
                let end = min(window.end, entry.interval().end) as i32;
                let start = max(window.start, entry.interval().start) as i32;
                assert!(end >= start, "Should overlap, come from tree query");
                end - start
            })
            .sum::<i32>();
        result.insert(String::from("MS"), ms);

        // Raw read counts.
        result.insert(
            String::from("RRC"),
            self.counters[window_id as usize] as i32,
        );

        // Ratio for scaling up to the read count.
        let len = (window.end - window.start) as i32;
        let ratio = (len - ms) as f64 / len as f64;
        result.insert(String::from("MP"), (ratio < MAX_MS) as i32);

        // TODO: maybe better use floats for RC after all?
        result.insert(
            String::from("RC"),
            if ratio > 0.001 {
                // guard against NaN...
                (self.counters[window_id as usize] as f64 / ratio).round() as i32
            } else {
                0
            },
        );

        if window_id == 62_160_000 / 20_000 {
            println!("windows_id = {}, counts = {:?}", window_id, result);
        }

        result
    }

    fn float_values(&self, _window_id: u32) -> HashMap<String, f32> {
        HashMap::new()
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
        let ms = self.tree
            .find(window.clone())
            .map(|entry| {
                let end = min(window.end, entry.interval().end) as i32;
                let start = max(window.start, entry.interval().start) as i32;
                assert!(end >= start, "Should overlap, come from tree query");
                end - start
            })
            .sum::<i32>();

        // Ratio for scaling up to the read count.
        let len = (window.end - window.start) as i32;
        let ratio = (len - ms) as f64 / len as f64;

        ratio < MAX_MS
    }

    fn get_coverage(&self, window_id: u32) -> i32 {
        return self.counters[window_id as usize] as i32;
    }
}

/// Bin for coverage aggregation.
#[derive(Debug)]
pub struct CoverageBin {}

/// Struct for aggregating as coverage.
#[derive(Debug)]
pub struct CoverageAggregator {
    /// Common information from all aggregators.
    pub base: BaseAggregator,
    /// Number of bases for each base in the current bin for the one sample in the BAM file.
    pub coverage: Vec<CoverageBin>,
}

impl CoverageAggregator {
    pub fn new(
        samples: Vec<String>,
        read_group_to_sample: HashMap<String, String>,
        sample_to_idx: HashMap<String, u32>,
        options: Options,
        contig_length: u64,
    ) -> CoverageAggregator {
        CoverageAggregator {
            base: BaseAggregator {
                samples,
                read_group_to_sample,
                sample_to_idx,
                options,
                contig_length,
                num_processed: 0,
                num_skipped: 0,
            },
            coverage: Vec::new(),
        }
    }
}

impl BamRecordAggregator for CoverageAggregator {
    fn put_bam_record(&mut self, _record: &bam::Record) {
        panic!("XXX TODO WRITE ME XXX TODO");
    }

    fn integer_field_names(&self) -> Vec<String> {
        Vec::new()
    }

    fn float_field_names(&self) -> Vec<String> {
        vec![String::from("COV"), String::from("WINSD")]
    }

    fn integer_values(&self, _window_id: u32) -> HashMap<String, i32> {
        panic!("Implement me!");
        // HashMap::new()
    }

    fn float_values(&self, _window_id: u32) -> HashMap<String, f32> {
        panic!("Implement me!");
        // HashMap::new()
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

    fn get_coverage(&self, _window_id: u32) -> i32 {
        panic!("Implement me!");
    }
}
