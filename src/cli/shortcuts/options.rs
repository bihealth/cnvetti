use clap::ArgMatches;
use tempdir::TempDir;

pub use cli::options::*;

use cli::call::options as call_options;
use cli::coverage::options as coverage_options;
use cli::normalize::options as normalize_options;
use cli::options as cli_options;
use cli::segment::options as segment_options;

/// Options for the "wgs-deep" command.
#[derive(Clone, Debug)]
pub struct WgsDeepOptions {
    /// Path to reference file.
    pub reference: String,
    /// Path to input file.
    pub input: String,
    /// Path to output file.
    pub output: String,
    /// Path to segmentation output file.
    pub output_segmentation: Option<String>,

    // Regular expression of contigs to use for normalization.
    pub contig_regex: String,

    /// Number of threads to use in I/O.
    pub io_threads: u32,
}

// Implementation of constructor.
impl WgsDeepOptions {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Self {
        Self {
            reference: matches.value_of("reference").unwrap().to_string(),
            input: matches.value_of("input").unwrap().to_string(),
            output: matches.value_of("output").unwrap().to_string(),
            output_segmentation: matches
                .value_of("output_segmentation")
                .map(|s| s.to_string()),
            io_threads: matches
                .value_of("io_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
            contig_regex: matches.value_of("contig_regex").unwrap().to_string(),
        }
    }
}

/// Transmogrification into sub-step options.
impl WgsDeepOptions {
    /// Construct configuration for `cnvetti coverage`.
    pub fn to_coverage_options(&self, tmp_dir: &TempDir) -> coverage_options::Options {
        coverage_options::Options {
            // Values taken from command options.
            reference: self.reference.clone(),
            input: self.input.clone(),
            io_threads: self.io_threads,

            // Temporary files.
            output: tmp_dir
                .path()
                .join("coverage.bcf")
                .to_str()
                .unwrap()
                .to_string(),

            // Preconfigured defaults.
            genome_regions: Vec::new(),
            mapability_bed: None,
            blacklist_bed: None,
            output_bam: None,
            output_bed: None,
            window_length: 500, // TODO: make configurable?
            count_kind: cli_options::CountKind::Coverage,
            min_mapq: 0,
            min_unclipped: 0.6,
            skip_discordant: false,
            mask_piles: false,
            pile_depth_percentile: 0.0,
            pile_max_gap: 0,
            gc_step: 0.01,
        }
    }

    /// Construct configuration for `cnvetti normalize`.
    pub fn to_normalize_options(&self, tmp_dir: &TempDir) -> normalize_options::Options {
        let coverage_options = self.to_coverage_options(tmp_dir);

        normalize_options::Options {
            // Values taken from command options.
            io_threads: self.io_threads,
            contig_regex: self.contig_regex.clone(),

            // Temporary files.
            input: tmp_dir
                .path()
                .join("coverage.bcf")
                .to_str()
                .unwrap()
                .to_string(),
            output: tmp_dir
                .path()
                .join("normalized.bcf")
                .to_str()
                .unwrap()
                .to_string(),

            // Preconfigured defaults / values from previous step.
            min_gc_window_count: 10,
            gc_step: coverage_options.gc_step,
            mapability_step: 0.01,
            normalization: normalize_options::Normalization::BinnedGc,
        }
    }

    /// Construct configuration for `cnvetti segment`.
    pub fn to_segment_options(&self, tmp_dir: &TempDir) -> segment_options::Options {
        let coverage_options = self.to_coverage_options(tmp_dir);

        segment_options::Options {
            // Values taken from command options.
            io_threads: self.io_threads,

            // Temporary files.
            input: tmp_dir
                .path()
                .join("normalized.bcf")
                .to_str()
                .unwrap()
                .to_string(),
            output: match self.output_segmentation {
                Some(ref p) => p.clone(), // from command line
                None => tmp_dir
                    .path()
                    .join("segmented.bcf")
                    .to_str()
                    .unwrap()
                    .to_string(),
            },

            // Preconfigured defaults / values from previous step.
            count_kind: coverage_options.count_kind,
            min_mapability: 0.8,
            min_gc: 0.2,
            max_gc: 0.8,
            max_iqr: 0.1, // XXXX
            segmentation: segment_options::Segmentation::HaarSeg,
            haar_seg_l_min: 1,
            haar_seg_l_max: 5,
            haar_seg_breaks_fdr_q: 0.001,
            significant_p_val_thresh: 0.05,
        }
    }

    /// Construct configuration for `cnvetti call`.
    pub fn to_call_options(&self, tmp_dir: &TempDir) -> call_options::Options {
        let seg_opt = self.to_segment_options(tmp_dir);

        call_options::Options {
            // Values taken from command options.
            output: self.output.clone(),
            io_threads: self.io_threads,

            // Temporary files.
            input: match self.output_segmentation {
                Some(ref p) => p.clone(), // from command line
                None => tmp_dir
                    .path()
                    .join("segmented.bcf")
                    .to_str()
                    .unwrap()
                    .to_string(),
            },

            // Preconfigured defaults/values from previous step.
            significant_p_val_thresh: seg_opt.significant_p_val_thresh,
            merge_p_val_thresh: 0.05,
            merge_corr_factor: 2970000000.0, // TODO: can be generated from genome
            max_merges: 10,
        }
    }
}

/// Options for the "wgs-cov-bins" command.
#[derive(Clone, Debug)]
pub struct WgsCovBinsOptions {
    /// Path to reference file.
    pub reference: String,
    /// Path to input file.
    pub input: String,
    /// Path to output file.
    pub output: String,

    // Regular expression of contigs to use for normalization.
    pub contig_regex: String,

    /// Number of threads to use in I/O.
    pub io_threads: u32,
}

// Implementatin of constructor.
impl WgsCovBinsOptions {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Self {
        Self {
            reference: matches.value_of("reference").unwrap().to_string(),
            input: matches.value_of("input").unwrap().to_string(),
            output: matches.value_of("output").unwrap().to_string(),
            io_threads: matches
                .value_of("io_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
            contig_regex: matches.value_of("contig_regex").unwrap().to_string(),
        }
    }
}

/// Transmogrification into sub-step options.
impl WgsCovBinsOptions {
    /// Construct configuration for `cnvetti coverage`.
    pub fn to_coverage_options(&self, tmp_dir: &TempDir) -> coverage_options::Options {
        coverage_options::Options {
            // Values taken from command options.
            reference: self.reference.clone(),
            input: self.input.clone(),
            io_threads: self.io_threads,

            // Temporary files.
            output: tmp_dir
                .path()
                .join("coverage.bcf")
                .to_str()
                .unwrap()
                .to_string(),

            // Preconfigured defaults.
            genome_regions: Vec::new(),
            mapability_bed: None,
            blacklist_bed: None,
            output_bam: None,
            output_bed: None,
            window_length: 500, // TODO: make configurable?
            count_kind: cli_options::CountKind::Coverage,
            min_mapq: 0,
            min_unclipped: 0.6,
            skip_discordant: false,
            mask_piles: false,
            pile_depth_percentile: 0.0,
            pile_max_gap: 0,
            gc_step: 0.01,
        }
    }

    /// Construct configuration for `cnvetti normalize`.
    pub fn to_normalize_options(&self, tmp_dir: &TempDir) -> normalize_options::Options {
        let coverage_options = self.to_coverage_options(tmp_dir);

        normalize_options::Options {
            // Values taken from command options.
            io_threads: self.io_threads,
            contig_regex: self.contig_regex.clone(),

            // Temporary files.
            input: tmp_dir
                .path()
                .join("coverage.bcf")
                .to_str()
                .unwrap()
                .to_string(),
            output: self.output.clone(),

            // Preconfigured defaults / values from previous step.
            min_gc_window_count: 10,
            gc_step: coverage_options.gc_step,
            mapability_step: 0.01,
            // TODO: maybe disable GC correction for ASCAT input?
            normalization: normalize_options::Normalization::BinnedGc,
        }
    }
}
