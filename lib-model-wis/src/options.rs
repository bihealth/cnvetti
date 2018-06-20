use clap::ArgMatches;

/// Options for "cnvetti cmd build-model-wis".
#[derive(Clone, Debug)]
pub struct BuildModelWisOptions {
    /// Path to input file.
    pub input: String,
    /// Path to output file.
    pub output: String,

    // Number of additional threads to use for I/O.
    pub io_threads: u32,
    // Number of compute threads to use.
    pub num_threads: u32,

    // Threshold on z score for filtration.
    pub filter_z_score: f64,
    // Threshold on relative score.
    pub filter_rel: f64,
    // Smallest number of reference targets to accept.
    pub min_ref_targets: usize,
    // Number of reference targets to use for reference.
    pub max_ref_targets: usize,
    // Number of samples a region can be called in before it is flagged as unreliable.
    pub max_samples_reliable: u32,
    // Minimal number of samples that must have coverage here above "min_fragments" from
    // normalize.
    pub min_samples_min_fragments: u32,
}

/// Implementation of constructor.
impl BuildModelWisOptions {
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
            num_threads: matches
                .value_of("num_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
            filter_z_score: matches
                .value_of("filter_z_score")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
            filter_rel: matches
                .value_of("filter_rel")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
            min_ref_targets: matches
                .value_of("min_ref_targets")
                .unwrap()
                .parse::<usize>()
                .unwrap(),
            max_ref_targets: matches
                .value_of("max_ref_targets")
                .unwrap()
                .parse::<usize>()
                .unwrap(),
            max_samples_reliable: matches
                .value_of("max_samples_reliable")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
            min_samples_min_fragments: matches
                .value_of("min_samples_min_fragments")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
        }
    }
}
