use std::collections::HashMap;
use std::ops::Range;
use std::cmp::{max, min};

use cli::coverage::options::*;

use bio::data_structures::interval_tree;

use rust_htslib::bam;
use rust_htslib::bcf;


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
}

/// Trait for alignment aggregation from BAM files.
pub trait BamRecordAggregator {
    /// Put the given record into the aggregation.
    fn put_bam_record(&mut self, record: &bam::Record);

    // TODO: more elegant solution for pushing format values.

    /// Initialize arrays for counters in `record`.
    fn prepare_bcf_record(&self, record: &mut bcf::Record, num_sample: usize);

    /// Write out the statistics for the given sample to the given `record`.
    fn fill_bcf_record(&self, record: &mut bcf::Record, sample_id: usize, window_id: u32);
}


/// Struct for aggregating as alignment counts.
#[derive(Debug)]
pub struct CountAlignmentsAggregator<'a> {
    /// Common information from all aggregators.
    base: BaseAggregator,

    // TODO: lifetime not necessary?
    /// The `IntervalTree` with black-listed intervals.
    tree: &'a interval_tree::IntervalTree<u32, i32>,

    /// Per-sample read counts (each bin is a `u32`).
    counters: Vec<u32>,
}


impl<'a> CountAlignmentsAggregator<'a> {
    /// Construct new aggregator with the given BAM `header`.  This information is
    /// necessary to appropriately allocate buffers for all samples in the header.
    pub fn new(
        tree: &'a interval_tree::IntervalTree<u32, i32>,
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
    }

    /// Skip `record` because of discordant alignment?
    fn skip_discordant(&self, record: &bam::Record) -> bool {
        self.base.options.skip_discordant && !record.is_paired() &&
            (record.tid() != record.mtid() || record.is_reverse() == record.is_mate_reverse() ||
                 record.is_unmapped() || record.is_mate_unmapped())
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

        (num_unclipped as f32) / ((num_clipped + num_unclipped) as f32) <
            self.base.options.min_unclipped
    }

    // Skip `record` because it is paired and not the first fragment?
    fn skip_paired_and_all_but_first(&self, record: &bam::Record) -> bool {
        !record.is_paired() && !record.is_first_in_template()
    }
}


impl<'a> BamRecordAggregator for CountAlignmentsAggregator<'a> {
    fn put_bam_record(&mut self, record: &bam::Record) {
        if !self.skip_mapq(record) && !self.skip_flags(record) &&
            !self.skip_discordant(record) && !self.skip_clipping(record) &&
            !self.skip_paired_and_all_but_first(record)
        {
            let pos = (record.pos() as u32)..((record.pos() + 1) as u32);
            let window_length = self.base.options.window_length as usize;
            let bin = record.pos() as usize / window_length;
            if !self.tree.find(pos).next().is_none() {
                self.counters[bin] += 1;
            }
        }
    }

    fn prepare_bcf_record(&self, record: &mut bcf::Record, num_samples: usize) {
        record
            .push_format_integer(b"RC", vec![0 as i32; num_samples].as_slice())
            .expect("Writing FORMAT/RC failed");
        record
            .push_format_integer(b"MS", vec![0 as i32; num_samples].as_slice())
            .expect("Writing FORMAT/MS failed");
    }

    fn fill_bcf_record(&self, record: &mut bcf::Record, sample_id: usize, window_id: u32) {
        let window_id = window_id as usize;

        // Write out read count.
        let rc = [self.counters[window_id] as i32];
        // Scope to only have record mutably once.
        // TODO: is there a more elegant solution?
        {
            let mut rcs = record.format(b"RC").integer_mut().expect(
                "FORMAT/RC not found",
            );
            rcs[0].copy_from_slice(&rc);
            println!("Writing {:?} into RC => {:?}", rc, rcs);
        }

        // Compute number of masked bases in window.
        let window = Range {
            start: (window_id * self.base.options.window_length as usize) as u32,
            end: ((window_id + 1) * self.base.options.window_length as usize) as u32,
        };
        let ms = [
            self.tree
                .find(window.clone())
                .map(|entry| {
                    let end = min(window.end, entry.interval().end) as i32;
                    let start = max(window.start, entry.interval().start) as i32;
                    assert!(end >= start, "Should overlap, come from tree query");
                    end - start
                })
                .sum::<i32>(),
        ];
        // Scope to only have record mutably once.
        // TODO: is there a more elegant solution?
        {
            let mut mss = record.format(b"MS").integer_mut().expect(
                "FORMAT/MS not found",
            );
            mss[0].copy_from_slice(&ms);
            println!("Writing {:?} into MS => {:?}", ms, mss);
        }
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
            },
            coverage: Vec::new(),
        }
    }
}

impl BamRecordAggregator for CoverageAggregator {
    fn put_bam_record(&mut self, record: &bam::Record) {
        panic!("XXX TODO WRITE ME XXX TODO");
    }

    fn prepare_bcf_record(&self, record: &mut bcf::Record, num_sample: usize) {
        panic!("XXX TODO WRITE ME XXX TODO");
    }

    fn fill_bcf_record(&self, record: &mut bcf::Record, sample_id: usize, window_id: u32) {
        panic!("XXX TODO WRITE ME XXX TODO");
    }
}
