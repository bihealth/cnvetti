//! Implementation of the "cnvetti quick pool-call" command.

use std::str::FromStr;

use clap::ArgMatches;

use slog::Logger;

use tempdir::TempDir;

use lib_coverage::{self, CoverageOptions};
use lib_genotype::{self, GenotypeOptions, GenotypingMethod};
use lib_mod_cov::{self, ModelBasedCoverageOptions};
use lib_normalize::{self, NormalizeOptions};
use lib_segment::{self, SegmentOptions, Segmentation};
use lib_visualize::{self, CovToIgvOptions};

use super::errors::*;

/// Options for "cnvetti cmd coverage".
#[derive(Clone, Debug)]
pub struct QuickPoolCallOptions {
    /// Path to input BCF file.
    pub input: String,
    /// Path to indexed BCF file with targets.
    pub input_model: String,

    /// Path to output call BCF file.
    pub output: String,
    /// Path to output per-target BCF file.
    pub output_targets: Option<String>,

    /// Path to IGV file with CV information.
    pub output_igv_cov: Option<String>,
    /// Path to IGV file with CV2 information.
    pub output_igv_cov2: Option<String>,
    /// Path to IGV file with CVZ information.
    pub output_igv_covz: Option<String>,
    /// Path to IGV file with SG information.
    pub output_igv_seg: Option<String>,
    /// Path to IGV file with SG2 information.
    pub output_igv_seg2: Option<String>,

    /// The segmentation method to employ.
    pub segmentation: Segmentation,

    /// Parameter for p-value thresholding.
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
    /// Threshold on Z-score.
    pub xhmm_z_score_threshold: f64,
    /// Expected exome-wide CNV rate.
    pub xhmm_cnv_rate: f64,
    /// Expected mean number of targets.
    pub xhmm_mean_target_count: f64,
    /// Expected mean target distance in a CNV.
    pub xhmm_mean_target_dist: f64,
}

/// Conversion into CoverageOptions.
impl QuickPoolCallOptions {
    fn into_coverage_options(&self, input: &String, output: &String) -> CoverageOptions {
        CoverageOptions {
            io_threads: 0,
            input: input.clone(),
            output: output.clone(),
            output_masked_bed: None,

            reference: None,
            genome_region: None,

            contig_regex: "^(chr)?\\d\\d?$".to_string(),
            count_kind: lib_coverage::CountKind::Fragments,
            blacklist_bed: None,
            considered_regions: lib_coverage::ConsideredRegions::TargetRegions,
            min_mapq: 0,
            min_unclipped: 0.6,

            min_window_remaining: 0.5,
            min_raw_coverage: 10,

            window_length: None,
            targets_bed: None,
            model_bcf: Some(self.input_model.clone()),

            mask_piles: false,
            mask_piles_fdr: 0.0,
            pile_max_gap: 0,
        }
    }

    fn into_normalize_options(&self, input: &String, output: &String) -> NormalizeOptions {
        NormalizeOptions {
            input: input.clone(),
            output: output.clone(),

            io_threads: 0,
            normalization: lib_normalize::Normalization::CoverageMedian,
        }
    }

    fn into_build_model_based_coverage_options(
        &self,
        input: &String,
        input_model: &String,
        output: &String,
    ) -> ModelBasedCoverageOptions {
        ModelBasedCoverageOptions {
            input: input.clone(),
            input_wis_model: None,
            input_pool_model: Some(input_model.clone()),

            output: output.clone(),

            io_threads: 0,
        }
    }

    fn into_segment_options(&self, input: &String, output: &String) -> SegmentOptions {
        SegmentOptions {
            input: input.clone(),
            output: output.clone(),
            segmentation: self.segmentation,

            io_threads: 0,

            thresh_p_value: self.thresh_p_value,

            haar_seg_l_min: self.haar_seg_l_min,
            haar_seg_l_max: self.haar_seg_l_max,
            haar_seg_fdr: self.haar_seg_fdr,

            wisexome_max_window_size: self.wisexome_max_window_size,
            wisexome_thresh_rel_cov: self.wisexome_thresh_rel_cov,
            wisexome_thresh_z_score: self.wisexome_thresh_z_score,

            xhmm_z_score_threshold: self.xhmm_z_score_threshold,
            xhmm_cnv_rate: self.xhmm_cnv_rate,
            xhmm_mean_target_count: self.xhmm_mean_target_count,
            xhmm_mean_target_dist: self.xhmm_mean_target_dist,
        }
    }

    fn into_cov_to_igv_options(&self, input: &String) -> CovToIgvOptions {
        CovToIgvOptions {
            input: input.clone(),
            output_igv_cov: self.output_igv_cov.clone(),
            output_igv_cov2: self.output_igv_cov2.clone(),
            output_igv_covz: self.output_igv_covz.clone(),
            output_igv_seg: self.output_igv_seg.clone(),
            output_igv_seg2: self.output_igv_seg2.clone(),
            io_threads: 0,
        }
    }

    fn into_genotype_options(&self, input: &String) -> GenotypeOptions {
        GenotypeOptions {
            input: input.clone(),
            input_calls: None,
            output: self.output.clone(),
            io_threads: 0,

            genotyping: GenotypingMethod::ExomeHiddenMarkovModel,

            xhmm_z_score_threshold: self.xhmm_z_score_threshold,
            xhmm_cnv_rate: self.xhmm_cnv_rate,
            xhmm_mean_target_count: self.xhmm_mean_target_count,
            xhmm_mean_target_dist: self.xhmm_mean_target_dist,
        }
    }
}

impl QuickPoolCallOptions {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Self {
        let segmentation = matches.value_of("segmentation").unwrap();
        let segmentation = Segmentation::from_str(&segmentation).expect("Unknown segmentation");

        Self {
            input: matches
                .value_of("input")
                .expect("Problem getting input args from command line")
                .to_string(),
            input_model: matches
                .value_of("input_model")
                .expect("Problem getting input_model args from command line")
                .to_string(),

            output: matches.value_of("output").unwrap().to_string(),
            output_targets: matches.value_of("output_targets").map(|s| s.to_string()),

            output_igv_cov: matches.value_of("output_igv_cov").map(|s| s.to_string()),
            output_igv_cov2: matches.value_of("output_igv_cov2").map(|s| s.to_string()),
            output_igv_covz: matches.value_of("output_igv_covz").map(|s| s.to_string()),
            output_igv_seg: matches.value_of("output_igv_seg").map(|s| s.to_string()),
            output_igv_seg2: matches.value_of("output_igv_seg2").map(|s| s.to_string()),

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

pub fn run(logger: &mut Logger, options: &QuickPoolCallOptions) -> Result<()> {
    info!(logger, "Running: cnvetti quick pool-call");
    info!(logger, "Options: {:?}", options);

    let tmp_dir = TempDir::new("cnvetti_quick_pool_call")
        .chain_err(|| "Could not create temporary directory.")?;

    // Parallel coverage computation and normalization.
    info!(
        logger,
        "Running cnvetti cmd coverage && cnvetti cmd normalize for sample."
    );
    let cov_out = tmp_dir
        .path()
        .join(format!("coverage.bcf"))
        .to_str()
        .unwrap()
        .to_string();
    lib_coverage::run(
        &mut logger.new(o!("step" => "coverage")),
        &options.into_coverage_options(&options.input, &cov_out),
    ).chain_err(|| format!("Problem computing coverage on {}", &options.input))?;
    info!(logger, " => done");

    let norm_out = tmp_dir
        .path()
        .join(format!("normalized.bcf"))
        .to_str()
        .unwrap()
        .to_string();
    lib_normalize::run(logger, &options.into_normalize_options(&cov_out, &norm_out))
        .chain_err(|| format!("Problem normalizing on {}", &cov_out))?;
    info!(logger, " => done");

    // Compute coverage relative to pool-based model.
    info!(logger, "Compute model-based coverage normalization");
    let mod_cov_out = tmp_dir
        .path()
        .join(format!("output_mod_cov.bcf"))
        .to_str()
        .unwrap()
        .to_string();
    lib_mod_cov::run(
        &mut logger.new(o!("step" => "mod-coverage")),
        &options.into_build_model_based_coverage_options(
            &norm_out,
            &options.input_model,
            &mod_cov_out,
        ),
    ).chain_err(|| "Problem with merging coverage file")?;
    info!(logger, " => done");

    // Compute segmentation using the normalized WIS file.
    info!(logger, "Computing segmentation");
    let output_targets = if let Some(ref output_targets) = options.output_targets {
        output_targets.clone()
    } else {
        tmp_dir
            .path()
            .join(format!("output_segments.bcf"))
            .to_str()
            .unwrap()
            .to_string()
    };
    lib_segment::run(
        logger,
        &options.into_segment_options(&mod_cov_out, &output_targets),
    ).chain_err(|| format!("Problem segmenting on {}", &mod_cov_out))?;
    info!(logger, " => done");

    // Generate IGV output files for coverage.
    info!(logger, "Generate IGV output files");
    lib_visualize::cov_to_igv::run(
        &mut logger.new(o!("step" => "visualize:cov-to-igv")),
        &options.into_cov_to_igv_options(&output_targets),
    ).chain_err(|| "Problem with merging coverage file")?;
    info!(logger, " => done");

    // Generate Call output file.
    info!(logger, "Generate genotype output files");
    lib_genotype::run(
        &mut logger.new(o!("step" => "mod-genotype")),
        &options.into_genotype_options(&output_targets),
    ).chain_err(|| "Problem genotyping file")?;
    info!(logger, " => done");

    warn!(logger, "Generating visualization for genotype not implemented yet!");

    info!(logger, "All done. Have a nice day!");
    Ok(())
}
