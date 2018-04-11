use clap::{ArgMatches, Values};


/// Enum for selecting preset.
#[derive(Clone, Debug, PartialEq)]
pub enum OptionsPreset {
    Wgs,
    WesOffTarget,
    WesOnTarget,
}


impl OptionsPreset {
    /// Parse `CountKind` from `&str`.
    fn from_str(s: &str) -> Option<OptionsPreset> {
        match s {
            "Wgs" => Some(OptionsPreset::Wgs),
            "WesOffTarget" => Some(OptionsPreset::WesOffTarget),
            "WesOnTarget" => Some(OptionsPreset::WesOnTarget),
            _ => {
                panic!("Invalid preset {}", s);
            }
        }
    }
}


/// Enum for selecting count type.
#[derive(Clone, Debug, PartialEq)]
pub enum CountKind {
    Coverage,
    Alignments,
}


impl CountKind {
    /// Parse `CountKind` from `&str`.
    pub fn from_str(s: &str) -> Option<CountKind> {
        match s {
            "Coverage" => Some(CountKind::Coverage),
            "Alignments" => Some(CountKind::Alignments),
            _ => {
                panic!("Invalid count type {}", s);
            }
        }
    }
}


/// Options for the "coverage" command.
#[derive(Clone, Debug)]
pub struct Options {
    pub reference: String,
    pub input: Vec<String>,
    pub output: String,
    pub genome_regions: Vec<String>,
    pub mapability_bed: Option<String>,
    pub window_length: u64,
    pub window_overlap: u64,
    pub count_kind: CountKind,
    pub min_mapq: u8,
    pub min_unclipped: f32,
    pub skip_discordant: bool,
    pub mask_piles: bool,
    pub pile_min_depth: i32,
    pub pile_max_gap: u32,
    pub pile_mask_window_size: u32,
}


impl Options {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Options {
        // TODO: interpret `--preset`

        let count_kind = matches.value_of("count_kind").unwrap();

        Options {
            reference: matches.value_of("reference").unwrap().to_string(),
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
            mapability_bed: match matches.value_of("mapability_bed") {
                Some(x) => Some(x.to_string()),
                None => None,
            },
            window_length: matches
                .value_of("window_length")
                .unwrap()
                .parse::<u64>()
                .unwrap(),
            window_overlap: matches
                .value_of("window_overlap")
                .unwrap()
                .parse::<u64>()
                .unwrap(),
            count_kind: CountKind::from_str(count_kind).unwrap(),
            min_mapq: matches.value_of("min_mapq").unwrap().parse::<u8>().unwrap(),
            min_unclipped: matches
                .value_of("min_unclipped")
                .unwrap()
                .parse::<f32>()
                .unwrap(),
            skip_discordant: matches.is_present("skip_discordant"),
            mask_piles: matches.is_present("mask_piles"),
            pile_min_depth: matches
                .value_of("pile_min_depth")
                .unwrap()
                .parse::<i32>()
                .unwrap(),
            pile_max_gap: matches
                .value_of("pile_max_gap")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
            pile_mask_window_size: matches
                .value_of("pile_mask_window_size")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
        }
    }
}
