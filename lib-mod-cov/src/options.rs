/// Types for configuring the "coverage" sub command.
use clap::ArgMatches;

/// Options for the "coverage" command.
#[derive(Clone, Debug)]
pub struct ModelBasedCoverageOptions {
    // Generic arguments
    /// Number of background threads to use for compression/decompression in I/O.
    pub io_threads: u32,

    // I/O related
    /// Path to normalized counts BCF file.
    pub input: String,
    /// Path to model input BCF file.
    pub input_model: String,
    /// Path to normalized coverage depth BCF file.
    pub output: String,
}

impl ModelBasedCoverageOptions {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Self {
        Self {
            input: matches.value_of("input").unwrap().to_string(),
            input_model: matches.value_of("input_model").unwrap().to_string(),
            output: matches.value_of("output").unwrap().to_string(),

            io_threads: matches
                .value_of("io_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
        }
    }
}
