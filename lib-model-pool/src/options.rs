use clap::ArgMatches;

/// Options for "cnvetti cmd build-model-pool".
#[derive(Clone, Debug)]
pub struct BuildModelPoolOptions {
    /// Path to input file.
    pub input: String,
    /// Path to output file.
    pub output: String,

    // Number of additional threads to use for I/O.
    pub io_threads: u32,
}

/// Implementation of constructor.
impl BuildModelPoolOptions {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Self {
        Self {
            input: matches.value_of("input").unwrap().to_string(),
            output: matches.value_of("output").unwrap().to_string(),
            io_threads: matches
                .value_of("io_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
        }
    }
}
