/// Options for the visualization in CNVetti.
use clap::ArgMatches;

/// Options for "cnvetti visualize cov-to-igv".
#[derive(Clone, Debug)]
pub struct CovToIgvOptions {
    /// Path to input file.
    pub input: String,
    /// Path to output IGV file with relative coverage.
    pub output_igv_cov: Option<String>,
    /// Path to output IGV file with log2-scaled coverage.
    pub output_igv_cov2: Option<String>,
    /// Path to output IGV file with coverage Z-score.
    pub output_igv_covz: Option<String>,
    /// Path to output IGV file with segmented coverage.
    pub output_igv_seg: Option<String>,
    /// Path to output IGV file with log2-scaled segmented coverage.
    pub output_igv_seg2: Option<String>,

    // Number of additional threads to use for I/O.
    pub io_threads: u32,
}

/// Implementation of constructor.
impl CovToIgvOptions {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Self {
        Self {
            input: matches
                .value_of("input")
                .expect("Problem getting input args from command line")
                .to_string(),
            output_igv_cov: matches.value_of("output_igv_cov").map(|s| s.to_string()),
            output_igv_cov2: matches.value_of("output_igv_cov2").map(|s| s.to_string()),
            output_igv_covz: matches.value_of("output_igv_covz").map(|s| s.to_string()),
            output_igv_seg: matches.value_of("output_igv_seg").map(|s| s.to_string()),
            output_igv_seg2: matches.value_of("output_igv_seg2").map(|s| s.to_string()),
            io_threads: matches
                .value_of("io_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
        }
    }
}
