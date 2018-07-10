//! Options for the segmentation command.

use std::str::FromStr;

use clap::ArgMatches;

/// Define the method for segmentation.
#[derive(Clone, Copy, EnumString, Debug, PartialEq)]
pub enum Segmentation {
    /// Segmentation using the Haar-Seg algorithm.
    HaarSeg,
    /// Segmentation using the Circular Binary Segmentation algorithm.
    CircularBinarySegmentation,
    /// Segmentation using the genomic HMM algorithm.
    GenomeHiddenMarkovModel,
    /// Segmentation using the exome HMM algorithm.
    ExomeHiddenMarkovModel,
    /// Segmentation algorithm from WISExome paper.
    WISExome,
}

/// Options for "cnvetti cmd segment".
#[derive(Clone, Debug)]
pub struct SegmentOptions {
    /// Path to input BCF file.
    pub input: String,
    /// Path to output BCF file.
    pub output: String,
    /// Number of additional threads to use for (de-)compression in I/O.
    pub io_threads: u32,

    /// The segmentation method to employ.
    pub segmentation: Segmentation,

    /// Parameters for p-value thresholding.
    pub thresh_p_value: f64,

    // Parameters from Haar-Seg.
    /// Value for l_min.
    pub haar_seg_l_min: u32,
    /// Value for l_max.
    pub haar_seg_l_max: u32,
    /// Value for FDR.
    pub haar_seg_fdr: f64,

    // Parameters from WISExome.
    /// Maximal window size in "windowing" step.
    pub wisexome_max_window_size: u32,
    /// Threshold on relative coverage.
    pub wisexome_thresh_rel_cov: f64,
    /// Threshold on Z-score.
    pub wisexome_thresh_z_score: f64,

    // Parameters from XHMM.
    /// Z-score threshold to use.
    pub xhmm_z_score_threshold: f64,
    /// Expected exome-wide CNV rate.
    pub xhmm_cnv_rate: f64,
    /// Expected mean number of targets.
    pub xhmm_mean_target_count: f64,
    /// Expected mean target distance in a CNV.
    pub xhmm_mean_target_dist: f64,
}

impl SegmentOptions {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Self {
        let segmentation = matches.value_of("segmentation").unwrap();
        let segmentation = Segmentation::from_str(&segmentation).expect("Unknown segmentation");

        Self {
            input: matches.value_of("input").unwrap().to_string(),
            output: matches.value_of("output").unwrap().to_string(),
            io_threads: matches
                .value_of("io_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),

            segmentation: segmentation,

            thresh_p_value: matches
                .value_of("thresh_p_value")
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

            wisexome_max_window_size: matches
                .value_of("wisexome_max_window_size")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
            wisexome_thresh_rel_cov: matches
                .value_of("wisexome_thresh_rel_cov")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
            wisexome_thresh_z_score: matches
                .value_of("wisexome_thresh_z_score")
                .unwrap()
                .parse::<f64>()
                .unwrap(),

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
        }
    }
}
