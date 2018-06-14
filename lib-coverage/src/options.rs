use std::fmt;
/// Types for configuring the "coverage" sub command.
use std::str::FromStr;

use clap::{ArgMatches, Values};

/// Enum for selecting count type.
#[derive(Clone, Debug, PartialEq, EnumString)]
pub enum CountKind {
    Coverage,
    Fragments,
}

impl fmt::Display for CountKind {
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
    /// Minimal MAPQ of a read to count.
    pub min_mapq: u8,
    /// Minimal ratio of clipped bases a read must have to be considered.
    pub min_unclipped: f32,
    /// Whether or not to skip discordant reads (BAM flag).
    pub skip_discordant: bool,

    // Binning and targets
    /// The length of the windows.
    pub window_length: Option<u64>,
    /// Optional path to targets BED file in case of WES.
    pub targets_bed: Option<String>,

    // Pile-related
    /// Whether or not to mask based on overlap with large read piles.
    pub mask_piles: bool,
    /// Read must not be part of a pile in a percentile above this one (ordered by pile size in
    /// number of bases).
    pub pile_size_percentile: f64,
    /// Join piles closer than this number.
    pub pile_max_gap: u32,
}

impl CoverageOptions {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Self {
        CoverageOptions {
            input: matches.value_of("input").unwrap().to_string(),
            output: matches.value_of("output").unwrap().to_string(),
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
            min_mapq: matches.value_of("min_mapq").unwrap().parse::<u8>().unwrap(),
            min_unclipped: matches
                .value_of("min_unclipped")
                .unwrap()
                .parse::<f32>()
                .unwrap(),
            skip_discordant: matches.is_present("skip_discordant"),

            window_length: match matches.value_of("window_length") {
                Some(x) => Some(x.to_string().parse::<u64>().unwrap()),
                None => None,
            },
            targets_bed: match matches.value_of("targets_bed") {
                Some(x) => Some(x.to_string()),
                None => None,
            },

            mask_piles: matches.is_present("mask_piles"),
            pile_size_percentile: matches
                .value_of("pile_size_percentile")
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
