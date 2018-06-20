use clap::ArgMatches;

/// Options for "cnvetti cmd build-model-wis".
#[derive(Clone, Debug)]
pub struct MergeCovOptions {
    /// Path to input file.
    pub input: Vec<String>,
    /// Path to output file.
    pub output: String,

    // Number of additional threads to use for I/O.
    pub io_threads: u32,
}

/// Implementation of constructor.
impl MergeCovOptions {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Self {
        Self {
            input: matches
                .values_of("input")
                .expect("Problem getting input args from command line")
                .map(|s| s.to_string())
                .collect::<Vec<String>>(),
            output: matches.value_of("output").unwrap().to_string(),
            io_threads: matches
                .value_of("io_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
        }
    }
}
