/// Implementation of the "cnvetti quick wis-build-model" command.
use std::env;

use clap::ArgMatches;

use rayon::prelude::*;

use slog::Logger;

use tempdir::TempDir;

use lib_coverage::{self, CoverageOptions};
use lib_merge_cov::{self, MergeCovOptions};
use lib_model_wis::{self, BuildModelWisOptions};
use lib_normalize::{self, NormalizeOptions};

use super::errors::*;

/// Options for "cnvetti cmd coverage".
#[derive(Clone, Debug)]
pub struct QuickWisBuildModelOptions {
    /// Path to input BCF file.
    pub input: Vec<String>,
    /// Path to tabix-indexed targets BED file.
    pub targets_bed: String,

    /// Path to output model BCF file.
    pub output: String,
    /// Path to output joint coverage BCF file.
    pub output_cov: Option<String>,

    /// Number of threads to use for parallel coverage computation, normalization, and
    /// model building.
    pub num_threads: u32,
}

/// Conversion into CoverageOptions.
impl QuickWisBuildModelOptions {
    fn into_coverage_options(&self, input: &String, output: &String) -> CoverageOptions {
        CoverageOptions {
            io_threads: 0,
            input: input.clone(),
            output: output.clone(),

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
            targets_bed: Some(self.targets_bed.clone()),
            wis_model_bcf: None,

            mask_piles: false,
            pile_size_percentile: 0.0,
            pile_max_gap: 0,
        }
    }

    fn into_normalize_options(&self, input: &String, output: &String) -> NormalizeOptions {
        NormalizeOptions {
            input: input.clone(),
            output: output.clone(),

            io_threads: 0,
            normalization: lib_normalize::Normalization::TotalCovSum,
        }
    }

    fn into_merge_cov_options(&self, input: &Vec<String>, output: &String) -> MergeCovOptions {
        MergeCovOptions {
            input: input.clone(),
            output: output.clone(),
            io_threads: self.num_threads,
        }
    }

    fn into_build_model_wis_options(
        &self,
        input: &String,
        output: &String,
    ) -> BuildModelWisOptions {
        BuildModelWisOptions {
            input: input.clone(),
            output: output.clone(),
            io_threads: 0,
            num_threads: self.num_threads,
            filter_z_score: 5.64,
            filter_rel: 0.35,
            min_ref_targets: 10,
            max_ref_targets: 100,
            max_samples_reliable: 8,
            min_samples_min_fragments: 10,
        }
    }
}

impl QuickWisBuildModelOptions {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Self {
        Self {
            input: matches
                .values_of("input")
                .expect("Problem getting input args from command line")
                .map(|s| s.to_string())
                .collect::<Vec<String>>(),
            targets_bed: matches.value_of("targets_bed").unwrap().to_string(),

            output: matches.value_of("output").unwrap().to_string(),
            output_cov: match matches.value_of("output_cov") {
                Some(output_cov) => Some(output_cov.to_string()),
                None => None,
            },

            num_threads: matches
                .value_of("num_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),
        }
    }
}

pub fn run(logger: &mut Logger, options: &QuickWisBuildModelOptions) -> Result<()> {
    // Set number of threads to use by rayon.
    env::set_var("RAYON_NUM_THREADS", format!("{}", options.num_threads));

    info!(logger, "Running: cnvetti quick wis-build-model");
    info!(logger, "Options: {:?}", options);

    let tmp_dir = TempDir::new("cnvetti_quick_wis_build_model")
        .chain_err(|| "Could not create temporary directory.")?;

    // Parallel coverage computation and normalization.
    info!(
        logger,
        "Running cnvetti cmd coverage && cnvetti cmd normalize for all samples."
    );
    let norm_out = options
        .input
        .par_iter()
        .enumerate()
        .map(|(idx, input)| -> Result<String> {
            let mut logger =
                logger.new(o!("step" => "cov-norm", "input-file" => format!("{}", idx)));
            info!(logger, "This thread is processing {}", input);

            let cov_out = tmp_dir
                .path()
                .join(format!("coverage.{}.bcf", idx))
                .to_str()
                .unwrap()
                .to_string();
            lib_coverage::run(
                &mut logger,
                &options.into_coverage_options(&input, &cov_out),
            ).chain_err(|| format!("Problem computing coverage on {}", &input))?;

            let norm_out = tmp_dir
                .path()
                .join(format!("normalized.{}.bcf", idx))
                .to_str()
                .unwrap()
                .to_string();
            lib_normalize::run(
                &mut logger,
                &options.into_normalize_options(&cov_out, &norm_out),
            ).chain_err(|| format!("Problem normalizing on {}", &cov_out))?;

            Ok(norm_out)
        })
        .collect::<Result<Vec<String>>>()
        .chain_err(|| "Problem with coverage-normalize step")?;

    // Merging of the resulting files, maybe writing to user output.
    info!(logger, "Merge result of normalized file");
    let merge_out = if let Some(ref output_cov) = options.output_cov {
        output_cov.clone()
    } else {
        tmp_dir
            .path()
            .join(format!("normalized.merged.bcf"))
            .to_str()
            .unwrap()
            .to_string()
    };
    lib_merge_cov::run(
        &mut logger.new(o!("step" => "merge-cov")),
        &options.into_merge_cov_options(&norm_out, &merge_out),
    ).chain_err(|| "Problem with merging coverage file")?;

    // Build within-sample model.
    info!(logger, "Build model and write out");
    lib_model_wis::run(
        &mut logger.new(o!("step" => "build-wis-model")),
        &options.into_build_model_wis_options(&merge_out, &options.output),
    ).chain_err(|| "Problem building the model")?;

    info!(logger, "All done. Have a nice day!");
    Ok(())
}
