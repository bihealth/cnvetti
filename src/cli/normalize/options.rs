use clap::ArgMatches;

pub use cli::options::*;

/// Define the normalization count.
#[derive(Clone, Debug, PartialEq)]
pub enum Normalization {
    BinnedGc,
    LowessGc,
    LowessGcMapability,
}

impl Normalization {
    /// Parse `Normalization` from `&str`.
    pub fn from_str(s: &str) -> Option<Normalization> {
        match s {
            "BinnedGc" => Some(Normalization::BinnedGc),
            "LowessGc" => Some(Normalization::LowessGc),
            "LowessGcMapability" => Some(Normalization::LowessGcMapability),
            _ => None,
        }
    }
}

/// Options for the "coverage" command.
#[derive(Clone, Debug)]
pub struct Options {
    pub input: String,
    pub output: String,
    pub min_gc_window_count: i32, // can go away
    pub io_threads: u32,
    pub gc_step: f64,         // should come from stats file
    pub mapability_step: f64, // can go away
    pub normalization: Normalization,
    pub contig_regex: String,
}

impl Options {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Options {
        let normalization = matches.value_of("normalization").unwrap();

        Options {
            input: matches.value_of("input").unwrap().to_string(),
            output: matches.value_of("output").unwrap().to_string(),
            min_gc_window_count: matches
                .value_of("min_gc_window_count")
                .unwrap()
                .parse::<i32>()
                .unwrap(),
            io_threads: matches
                .value_of("io_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
            gc_step: matches.value_of("gc_step").unwrap().parse::<f64>().unwrap(),
            contig_regex: matches.value_of("contig_regex").unwrap().to_string(),
            mapability_step: matches
                .value_of("mapability_step")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
            normalization: Normalization::from_str(normalization).expect("Unknown normalization"),
        }
    }
}
