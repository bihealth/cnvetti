use clap::{ArgMatches, Values};

pub use cli::options::*;

// TODO: add back multi-input file mode.
// TODO: use window overlap

/// Options for the "coverage" command.
#[derive(Clone, Debug)]
pub struct Options {
    /// Path to reference FASTA file.
    pub reference: String,
    /// Path to sample BAM file.
    pub input: String,
    /// Path to coverage depth BCF file.
    pub output: String,

    /// Optional path to BAM file to write counted reads to.
    pub output_bam: Option<String>,
    /// Optional path to BED file to write masked regions to.
    pub output_bed: Option<String>,

    /// Path to tabix-indexed BED file with mapability information in the fourth column.
    pub mapability_bed: Option<String>,
    /// Path to tabix-indexed BED file with black-listed regions.
    pub blacklist_bed: Option<String>,

    /// List of genome regions to limit analysis to.
    pub genome_regions: Vec<String>,

    /// The length of the windows.
    pub window_length: u64,
    /// Whether to count coverage/bases or alignments.
    pub count_kind: CountKind,

    /// Minimal MAPQ of a read to count.
    pub min_mapq: u8,
    /// Minimal ratio of clipped bases a read must have to be considered.
    pub min_unclipped: f32,

    /// Whether or not to skip discordant reads (BAM flag).
    pub skip_discordant: bool,

    /// Whether or not to mask based on overlap with large read piles.
    pub mask_piles: bool,
    /// Read must not be part of a pile in a percentile above this one (ordered by pile size in
    /// number of bases).
    pub pile_depth_percentile: f64,
    /// Join piles closer than this number.
    pub pile_max_gap: u32,

    /// Number of background threads to use for compression/decompression in I/O.
    pub io_threads: u32,

    /// GC bin size when using simple GC normalization.
    pub gc_step: f64,
}

impl Options {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Options {
        // TODO: interpret `--preset`

        let count_kind = matches.value_of("count_kind").unwrap();

        let mut options = Options {
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
            blacklist_bed: match matches.value_of("blacklist_bed") {
                Some(x) => Some(x.to_string()),
                None => None,
            },
            output_bam: match matches.value_of("output_bam") {
                Some(x) => Some(x.to_string()),
                None => None,
            },
            output_bed: match matches.value_of("output_bed") {
                Some(x) => Some(x.to_string()),
                None => None,
            },
            window_length: matches
                .value_of("window_length")
                .unwrap()
                .replace(",", "")
                .replace("_", "")
                .parse::<u64>()
                .unwrap(),
            count_kind: CountKind::from_str(count_kind).expect("Unknown count kind"),
            min_mapq: matches.value_of("min_mapq").unwrap().parse::<u8>().unwrap(),
            min_unclipped: matches
                .value_of("min_unclipped")
                .unwrap()
                .parse::<f32>()
                .unwrap(),
            skip_discordant: matches.is_present("skip_discordant"),
            mask_piles: matches.is_present("mask_piles"),
            pile_depth_percentile: matches
                .value_of("pile_depth_percentile")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
            pile_max_gap: matches
                .value_of("pile_max_gap")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
            io_threads: matches
                .value_of("io_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
            gc_step: matches.value_of("gc_step").unwrap().parse::<f64>().unwrap(),
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
                    // For WES off-target, we count alignments, use a larger window size, and
                    // enable appropriate pile removal and count alignments.
                    options.window_length = 20_000;

                    options.mask_piles = true;
                    options.pile_depth_percentile = 90.0;
                    options.pile_max_gap = 5;
                    options.skip_discordant = true;

                    options.count_kind = CountKind::Alignments;
                }
            }
        }

        options
    }
}
