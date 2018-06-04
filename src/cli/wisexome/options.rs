use clap::ArgMatches;

/// Options for the "wisexome count" command.
#[derive(Clone, Debug)]
pub struct CountOptions {
    /// Path to input file.
    pub input: String,
    /// Path to target BED file.
    pub target_bed: String,
    /// Path to output file.
    pub output: String,

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
            target_bed: matches.value_of("targets-bed").unwrap().to_string(),
            output: matches.value_of("output").unwrap().to_string(),
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
        }
    }
}
