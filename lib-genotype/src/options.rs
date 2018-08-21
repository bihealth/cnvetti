//! Options for the genotyping command.

use std::str::FromStr;

use clap::ArgMatches;

/// Define the method for genotyping.
#[derive(Clone, Copy, EnumString, Debug, PartialEq)]
pub enum GenotypingMethod {
    /// Genotyping using the exome HMM algorithm.
    ExomeHiddenMarkovModel,
    /// Genotyping by segment overlap.
    SegmentOverlap,
}

/// Define the method for segmentation when using `SegmentOverlap` genotyping.
#[derive(Clone, Copy, EnumString, Debug, PartialEq)]
pub enum SegmentationMethod {
    /// HaarSeg-based segmentation.
    HaarSeg,
}

/// Options for "cnvetti cmd genotype".
#[derive(Clone, Debug)]
pub struct GenotypeOptions {
    /// Path to input segmentation or coverage BCF file.
    pub input: String,
    /// Path to input call BCF file.
    pub input_calls: Option<String>,
    /// Path to output BCF file.
    pub output: String,
    /// Number of additional threads to use for (de-)compression in I/O.
    pub io_threads: u32,

    /// The genotyping method to employ.
    pub genotyping: GenotypingMethod,

    /// The segmentation method to employ.
    pub segmentation: SegmentationMethod,
    /// Minimal overlap to call a segment.
    pub overlap: f64,

    // Parameters from XHMM.
    /// Z-score threshold to use.
    pub xhmm_z_score_threshold: f64,
    /// Expected exome-wide CNV rate.
    pub xhmm_cnv_rate: f64,
    /// Expected mean number of targets.
    pub xhmm_mean_target_count: f64,
    /// Expected mean target distance in a CNV.
    pub xhmm_mean_target_dist: f64,

    // Parameters from Haar-Seg.
    /// Value for l_min.
    pub haar_seg_l_min: u32,
    /// Value for l_max.
    pub haar_seg_l_max: u32,
    /// Value for FDR.
    pub haar_seg_fdr: f64,
}

impl GenotypeOptions {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Self {
        println!("{:?}", matches);

        let genotyping = matches.value_of("genotyping").unwrap();
        let genotyping = GenotypingMethod::from_str(&genotyping).expect("Unknown genotyping");

        let segmentation = matches.value_of("segmentation").unwrap();
        let segmentation = SegmentationMethod::from_str(&segmentation).expect("Unknown segmentation");

        Self {
            input: matches.value_of("input").unwrap().to_string(),
            input_calls: matches.value_of("input_calls").map(|s| s.to_string()),
            output: matches.value_of("output").unwrap().to_string(),
            io_threads: matches
                .value_of("io_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),

            genotyping: genotyping,
            segmentation: segmentation,
            overlap: matches.value_of("overlap").unwrap().parse::<f64>().unwrap(),

            xhmm_z_score_threshold: matches
                .value_of("xhmm_z_score_threshold")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
            xhmm_cnv_rate: matches
                .value_of("xhmm_cnv_rate")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
            xhmm_mean_target_count: matches
                .value_of("xhmm_mean_target_count")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
            xhmm_mean_target_dist: matches
                .value_of("xhmm_mean_target_dist")
                .unwrap()
                .parse::<f64>()
                .unwrap(),

            haar_seg_l_min: matches
                .value_of("haar_seg_l_min")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
            haar_seg_l_max: matches
                .value_of("haar_seg_l_max")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
            haar_seg_fdr: matches
                .value_of("haar_seg_fdr")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
        }
    }
}
