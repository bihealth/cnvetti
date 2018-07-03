//! Options for the segmentation command.

use std::str::FromStr;

use clap::ArgMatches;

/// Define the method for segmentation.
#[derive(Clone, EnumString, Debug, PartialEq)]
pub enum Segmentation {
    /// Segmentation using the Haar-Seg algorithm.
    HaarSeg,
    /// Segmentation using the Circular Binary Segmentation algorithm.
    CircularBinarySegmentation,
    /// Segmentation using the genomic HMM algorithm.
    GenomeHiddenMarkovModel,
    /// Segmentation using the exome HMM algorithm.
    ExomeHiddenMarkovModel,
    /// Segmentation algorithm from WISExome paper.
    WISExome,
}

/// Options for "cnvetti cmd segment".
#[derive(Clone, Debug)]
pub struct SegmentOptions {
    /// Path to input BCF file.
    pub input: String,
    /// Path to output BCF file.
    pub output: String,
    /// Number of additional threads to use for (de-)compression in I/O.
    pub io_threads: u32,

    /// The segmentation method to employ.
    pub segmentation: Segmentation,
}

impl SegmentOptions {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Self {
        let segmentation = matches.value_of("segmentation").unwrap();
        let segmentation = Segmentation::from_str(&segmentation).expect("Unknown segmentation");

        Self {
            input: matches.value_of("input").unwrap().to_string(),
            output: matches.value_of("output").unwrap().to_string(),
            io_threads: matches
                .value_of("io_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
            segmentation: segmentation,
        }
    }
}
