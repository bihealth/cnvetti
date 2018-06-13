use clap::ArgMatches;

/// Options for the "wisexome count" command.
#[derive(Clone, Debug)]
pub struct CountOptions {
    /// Path to input file.
    pub input: String,
    /// Path to target BED file.
    pub targets_bed: String,
    /// Path to output file.
    pub output: String,

    // Minimal MAPQ.
    pub min_mapq: u8,
    // Regular expression of contigs to use for normalization.
    pub contig_regex: String,
    // Number of additional threads to use for I/O.
    pub io_threads: u32,
}

// Implementation of constructor.
impl CountOptions {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Self {
        Self {
            input: matches.value_of("input").unwrap().to_string(),
            targets_bed: matches.value_of("targets-bed").unwrap().to_string(),
            output: matches.value_of("output").unwrap().to_string(),
            min_mapq: matches.value_of("min_mapq").unwrap().parse::<u8>().unwrap(),
            io_threads: matches
                .value_of("io_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
            contig_regex: matches.value_of("contig_regex").unwrap().to_string(),
        }
    }
}

/// Options for the "wisexome normalize" command.
#[derive(Clone, Debug)]
pub struct NormalizeOptions {
    /// Path to input file.
    pub input: String,
    /// Path to output file.
    pub output: String,

    // Number of additional threads to use for I/O.
    pub io_threads: u32,
    // Minimal fragment count.
    pub min_fragments: u32,
}

// Implementation of constructor.
impl NormalizeOptions {
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
            min_fragments: matches
                .value_of("min_fragments")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
        }
    }
}

/// Options for the "wisexome build-ref" command.
#[derive(Clone, Debug)]
pub struct BuildRefOptions {
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

// Implementation of constructor.
impl BuildRefOptions {
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

/// Options for the "wisexome call" command.
#[derive(Clone, Debug)]
pub struct CallOptions {
    /// Path to input file.
    pub input: String,
    /// Path to reference input file.
    pub input_ref: String,
    /// Path to output file.
    pub output: String,
    /// Path to per-probe output file.
    pub output_probes: String,

    // Number of additional threads to use for I/O.
    pub io_threads: u32,
    // Number of compute threads to use.
    pub num_threads: u32,
    // Maximum window size for calling.
    pub max_window_size: usize,

    // Threshold on relative coverage deviation.
    pub thresh_rel: f64,
    // Thresehold on z-score.
    pub thresh_z: f64,
}

// Implementation of constructor.
impl CallOptions {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Self {
        Self {
            input: matches.value_of("input").unwrap().to_string(),
            input_ref: matches.value_of("input_ref").unwrap().to_string(),
            output: matches.value_of("output").unwrap().to_string(),
            output_probes: matches.value_of("output_probes").unwrap().to_string(),
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
            max_window_size: matches
                .value_of("max_window_size")
                .unwrap()
                .parse::<usize>()
                .unwrap(),
            thresh_rel: matches
                .value_of("thresh_rel")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
            thresh_z: matches
                .value_of("thresh_z")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
        }
    }
}
