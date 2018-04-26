use clap::ArgMatches;

/// Options for the "coverage" command.
#[derive(Clone, Debug)]
pub struct Options {
    pub input: String,
    pub output: String,
    pub io_threads: u32,
    pub percentile_threshold: f64,
}

impl Options {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Options {
        Options {
            input: matches.value_of("input").unwrap().to_string(),
            output: matches.value_of("output").unwrap().to_string(),
            io_threads: matches
                .value_of("io_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
            percentile_threshold: matches
                .value_of("iqr_percentile_threshold")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
        }
    }
}
