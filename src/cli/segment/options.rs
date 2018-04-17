use clap::ArgMatches;

pub use cli::options::*;

/// Options for the "coverage" command.
#[derive(Clone, Debug)]
pub struct Options {
    /// Path to input file.
    pub input: String,
    /// Path to output file.
    pub output: String,

    /// Number of threads to use in I/O.
    pub io_threads: u32,

    /// Minimal mapability to require.
    pub min_mapability: f32,
    /// Minimal GC content to require.
    pub min_gc: f32,
    /// Maximal GC content to allow.
    pub max_gc: f32,
    /// Maximal relative coverage/count IQR to allow.
    pub max_iqr: f32,

    /// Segmentation algorithm to use.
    pub segmentation: Segmentation,

    /// HaarSeg: value for L_min.
    pub haar_seg_l_min: i32,
    /// HaarSeg: value for L_max.
    pub haar_seg_l_max: i32,
    /// HaarSeg: value for the FDR "Q" value.
    pub haar_seg_breaks_fdr_q: f64,
}

impl Options {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Options {
        let segmentation = matches.value_of("segmentation").unwrap();

        Options {
            input: matches.value_of("input").unwrap().to_string(),
            output: matches.value_of("output").unwrap().to_string(),
            segmentation: Segmentation::from_str(segmentation)
                .expect("Unknown segmentation algorithm"),
            io_threads: matches
                .value_of("io_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
            min_mapability: matches
                .value_of("min_mapability")
                .unwrap()
                .parse::<f32>()
                .unwrap(),
            min_gc: matches.value_of("min_gc").unwrap().parse::<f32>().unwrap(),
            max_gc: matches.value_of("max_gc").unwrap().parse::<f32>().unwrap(),
            max_iqr: matches.value_of("max_iqr").unwrap().parse::<f32>().unwrap(),
            haar_seg_l_min: matches
                .value_of("haar_seg_l_min")
                .unwrap()
                .parse::<i32>()
                .unwrap(),
            haar_seg_l_max: matches
                .value_of("haar_seg_l_max")
                .unwrap()
                .parse::<i32>()
                .unwrap(),
            haar_seg_breaks_fdr_q: matches
                .value_of("haar_seg_breaks_fdr_q")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
        }
    }
}
