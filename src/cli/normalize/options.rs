use clap::{ArgMatches, Values};

pub use cli::options::*;

/// Options for the "coverage" command.
#[derive(Clone, Debug)]
pub struct Options {
    pub input: Vec<String>,
    pub output: String,
    pub genome_regions: Vec<String>,
    pub count_kind: CountKind,
    pub min_gc_window_count: i32,
}


impl Options {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Options {
        // TODO: interpret `--preset`

        let count_kind = matches.value_of("count_kind").unwrap();

        let mut options = Options {
            input: matches
                .values_of("input")
                .unwrap_or(Values::default())
                .map(|res| res.to_string())
                .collect(),
            output: matches.value_of("output").unwrap().to_string(),
            genome_regions: matches
                .values_of("genome_regions")
                .unwrap_or(Values::default())
                .map(|res| res.to_string())
                .collect(),
            count_kind: CountKind::from_str(count_kind).expect("Unknown count kind"),
            min_gc_window_count: matches
                .value_of("min_gc_window_count")
                .unwrap()
                .parse::<i32>()
                .unwrap(),
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