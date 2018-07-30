use clap::ArgMatches;

/// Options for "cnvetti cmd ratio"
#[derive(Clone, Debug)]
pub struct RatioOptions {
    /// Path to input file.
    pub input: String,
    /// Path to output file.
    pub output: String,

    /// Name of "numerator" sample.
    pub numerator_sample: String,
    /// Name of "denominator" sample.
    pub denominator_sample: String,

    // Number of additional threads to use for I/O.
    pub io_threads: u32,
}

/// Implementation of constructor.
impl RatioOptions {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Self {
        Self {
            input: matches
                .value_of("input")
                .expect("Problem getting inputargs from command line")
                .to_string(),
            output: matches.value_of("output").unwrap().to_string(),

            numerator_sample: matches
                .value_of("numerator_sample")
                .expect("Problem getting numerator_sample args from command line")
                .to_string(),
            denominator_sample: matches
                .value_of("denominator_sample")
                .expect("Problem getting denominator_sample args from command line")
                .to_string(),

            io_threads: matches
                .value_of("io_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
        }
    }
}
