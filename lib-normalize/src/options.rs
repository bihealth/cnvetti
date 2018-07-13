use std::str::FromStr;

use clap::ArgMatches;

/// Define the method for normalization.
#[derive(Clone, EnumString, Debug, PartialEq)]
pub enum Normalization {
    /// Normalize by sum of coverages.
    TotalCoverageSum,
    /// Normalize by median of coverages.
    CoverageMedian,
    /// Perform binning-based GC correction.
    MedianGcBinned,
}

/// Options for "cnvetti cmd coverage".
#[derive(Clone, Debug)]
pub struct NormalizeOptions {
    /// Path to input BCF file.
    pub input: String,
    /// Path to output BCF file.
    pub output: String,
    /// Number of additional threads to use for (de-)compression in I/O.
    pub io_threads: u32,

    /// The normalization method to employ.
    pub normalization: Normalization,
}

impl NormalizeOptions {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Self {
        let normalization = matches.value_of("normalization").unwrap();
        let normalization = Normalization::from_str(&normalization).expect("Unknown normalization");

        Self {
            input: matches.value_of("input").unwrap().to_string(),
            output: matches.value_of("output").unwrap().to_string(),
            io_threads: matches
                .value_of("io_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
            normalization: normalization,
        }
    }
}
