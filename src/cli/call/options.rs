use clap::ArgMatches;

pub use cli::options::*;

/// Options for the "coverage" command.
#[derive(Clone, Debug)]
pub struct Options {
    /// Path to input file.
    pub input: String,
    /// Path to output file.
    pub output: String,

    /// Number of threads to use in I/O.
    pub io_threads: u32,

    /// P-value threshold for accepting segment.
    pub significant_p_val_thresh: f64,
    /// P-value threshold for merging segments.
    pub merge_p_val_thresh: f64,
    /// Factor for multiple value correction (assuming human genome and 99% is CN neutral).
    pub merge_corr_factor: f64,
}

impl Options {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Options {
        let segmentation = matches.value_of("segmentation").unwrap();
        let count_kind = matches.value_of("count_kind").unwrap();

        Options {
            input: matches.value_of("input").unwrap().to_string(),
            output: matches.value_of("output").unwrap().to_string(),
            io_threads: matches
                .value_of("io_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
            significant_p_val_thresh: matches
                .value_of("significant_p_val_thresh")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
            merge_p_val_thresh: matches
                .value_of("merge_p_val_thresh")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
            merge_corr_factor: matches
                .value_of("merge_corr_factor")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
        }
    }
}
