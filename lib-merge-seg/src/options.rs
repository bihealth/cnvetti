use clap::ArgMatches;

/// Options for "cnvetti cmd merge-seg".
#[derive(Clone, Debug)]
pub struct MergeSegOptions {
    /// Path to input file.
    pub input: Vec<String>,
    /// Path to output file.
    pub output: String,

    // Number of additional threads to use for I/O.
    pub io_threads: u32,

    // Minimal reciprocal overlap as fraction.
    pub reciprocal_overlap: f64,
}

/// Implementation of constructor.
impl MergeSegOptions {
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
            reciprocal_overlap: matches
                .value_of("reciprocal_overlap")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
        }
    }
}
