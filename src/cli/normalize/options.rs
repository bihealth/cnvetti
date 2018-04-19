use clap::ArgMatches;

pub use cli::options::*;

/// Define the normalization count.
#[derive(Clone, Debug, PartialEq)]
pub enum Normalization {
    BinnedGc,
    LoessGc,
    LoessGcMapability,
}

impl Normalization {
    /// Parse `Normalization` from `&str`.
    pub fn from_str(s: &str) -> Option<Normalization> {
        match s {
            "BinnedGc" => Some(Normalization::BinnedGc),
            "LoessGc" => Some(Normalization::LoessGc),
            "LoessGcMapability" => Some(Normalization::LoessGcMapability),
            _ => None,
        }
    }
}

/// Options for the "coverage" command.
#[derive(Clone, Debug)]
pub struct Options {
    pub input: String,
    pub output: String,
    pub count_kind: CountKind,
    pub min_gc_window_count: i32,
    pub io_threads: u32,
    pub gc_step: f64,         // should come from stats file
    pub mapability_step: f64, // can go away
    pub normalization: Normalization,
    pub contig_regex: String,
}

impl Options {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Options {
        // TODO: interpret `--preset`

        let count_kind = matches.value_of("count_kind").unwrap();
        let normalization = matches.value_of("normalization").unwrap();

        let mut options = Options {
            input: matches.value_of("input").unwrap().to_string(),
            output: matches.value_of("output").unwrap().to_string(),
            count_kind: CountKind::from_str(count_kind).expect("Unknown count kind"),
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
        };

        if let Some(preset) = matches.value_of("preset") {
            match OptionsPreset::from_str(preset).expect("Unknown --preset") {
                OptionsPreset::Wgs => {
                    panic!("Wgs preset not implemented yet!");
                }
                OptionsPreset::WesOnTarget => {
                    panic!("WesOnTarget preset not implemented yet!");
                }
                OptionsPreset::WesOffTarget => {
                    // For WES off-target, we are counting reads and require
                    // at least 50 windows per GC content.

                    options.count_kind = CountKind::Alignments;
                    options.min_gc_window_count = 50;
                }
            }
        }

        options
    }
}
