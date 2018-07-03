/// Types for configuring the "coverage" sub command.
use std::fmt;
use std::str::FromStr;

use clap::{ArgMatches, Values};

/// Enum for selecting count type.
#[derive(Clone, Copy, Debug, PartialEq, EnumString)]
pub enum CountKind {
    Coverage,
    Fragments,
}

impl fmt::Display for CountKind {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

/// Enum for selecting the considered regions.
#[derive(Clone, Copy, Debug, PartialEq, EnumString)]
pub enum ConsideredRegions {
    GenomeWide,
    TargetRegions,
}

impl fmt::Display for ConsideredRegions {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

/// Options for the "coverage" command.
#[derive(Clone, Debug)]
pub struct CoverageOptions {
    // Generic arguments
    /// Number of background threads to use for compression/decompression in I/O.
    pub io_threads: u32,

    // I/O related
    /// Path to sample BAM file.
    pub input: String,
    /// Path to coverage depth BCF file.
    pub output: String,
    /// Path to BED file with masked regions.
    pub output_masked_bed: Option<String>,
    /// Path to reference FASTA file.
    pub reference: Option<String>,

    // Limit processing to region
    /// Optional genome region to limit analysis to.
    pub genome_region: Option<String>,

    // Counting-related
    /// Regular expression to match contigs.
    pub contig_regex: String,
    /// Whether to count coverage/bases or alignments.
    pub count_kind: CountKind,
    /// Optional path to BED file with blacklist.
    pub blacklist_bed: Option<String>,
    /// Whether to count on-target regions or genome-wide.
    pub considered_regions: ConsideredRegions,
    /// Minimal MAPQ of a read to count.
    pub min_mapq: u8,
    /// Minimal ratio of clipped bases a read must have to be considered.
    pub min_unclipped: f32,

    /// Minimal fraction of window that must remain after masking (e.g., for piles).
    pub min_window_remaining: f32,
    /// Minimal raw coverage before ignoring.
    pub min_raw_coverage: usize,

    // Binning and targets
    /// The length of the windows.
    pub window_length: Option<usize>,
    /// Optional path to targets BED file in case of WES.
    pub targets_bed: Option<String>,
    /// Path to WIS or pool-based model BCF file in case of WES with WIS model.
    pub model_bcf: Option<String>,

    // Pile-related
    /// Whether or not to mask based on overlap with large read piles.
    pub mask_piles: bool,
    /// Join piles closer than this number.
    pub pile_max_gap: u32,
    /// Set pile size threshold to get FDR less than or equal to this value using model
    /// assumption of Poisson for noise.
    pub mask_piles_fdr: f64,
}

impl CoverageOptions {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Self {
        Self {
            input: matches.value_of("input").unwrap().to_string(),
            output: matches.value_of("output").unwrap().to_string(),
            output_masked_bed: matches.value_of("output_masked_bed").map(|s| s.to_string()),
            reference: match matches.value_of("reference") {
                Some(x) => Some(x.to_string()),
                None => None,
            },

            genome_region: match matches.value_of("genome_region") {
                Some(x) => Some(x.to_string()),
                None => None,
            },

            contig_regex: matches
                .values_of("contig_regex")
                .unwrap_or(Values::default())
                .collect(),
            count_kind: CountKind::from_str(matches.value_of("count_kind").unwrap())
                .expect("Unknown count kind"),
            blacklist_bed: match matches.value_of("blacklist_bed") {
                Some(x) => Some(x.to_string()),
                None => None,
            },
            considered_regions: ConsideredRegions::from_str(
                matches.value_of("considered_regions").unwrap(),
            ).expect("Unknown considered regions"),
            min_mapq: matches.value_of("min_mapq").unwrap().parse::<u8>().unwrap(),
            min_unclipped: matches
                .value_of("min_unclipped")
                .unwrap()
                .parse::<f32>()
                .unwrap(),
            min_window_remaining: matches
                .value_of("min_window_remaining")
                .unwrap()
                .parse::<f32>()
                .unwrap(),
            min_raw_coverage: matches
                .value_of("min_raw_coverage")
                .unwrap()
                .parse::<usize>()
                .unwrap(),

            window_length: match matches.value_of("window_length") {
                Some(x) => Some(
                    x.to_string()
                        .replace("_", "")
                        .replace(",", "")
                        .parse::<usize>()
                        .unwrap(),
                ),
                None => None,
            },
            targets_bed: match matches.value_of("targets_bed") {
                Some(x) => Some(x.to_string()),
                None => None,
            },
            model_bcf: match matches.value_of("model_bcf") {
                Some(x) => Some(x.to_string()),
                None => None,
            },

            mask_piles: matches.is_present("mask_piles"),
            mask_piles_fdr: matches
                .value_of("mask_piles_fdr")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
            pile_max_gap: matches
                .value_of("pile_max_gap")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
            io_threads: matches
                .value_of("io_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
        }
    }
}
