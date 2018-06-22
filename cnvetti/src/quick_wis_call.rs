/// Implementation of the "cnvetti quick wis-call" command.
use clap::ArgMatches;

use slog::Logger;

use tempdir::TempDir;

use lib_coverage::{self, CoverageOptions};
use lib_mod_cov::{self, ModelBasedCoverageOptions};
use lib_normalize::{self, NormalizeOptions};

use super::errors::*;

/// Options for "cnvetti cmd coverage".
#[derive(Clone, Debug)]
pub struct QuickWisCallOptions {
    /// Path to input BCF file.
    pub input: String,
    /// Path to indexed BCF file with targets.
    pub input_model: String,

    /// Path to output call BCF file.
    pub output: String,
    /// Path to output per-target BCF file.
    pub output_targets: Option<String>,
}

/// Conversion into CoverageOptions.
impl QuickWisCallOptions {
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
            targets_bed: None,
            wis_model_bcf: Some(self.input_model.clone()),

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

    fn into_build_model_based_coverage_options(
        &self,
        input: &String,
        input_model: &String,
        output: &String,
    ) -> ModelBasedCoverageOptions {
        ModelBasedCoverageOptions {
            input: input.clone(),
            input_model: input_model.clone(),

            output: output.clone(),

            io_threads: 0,
        }
    }
}

impl QuickWisCallOptions {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Self {
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
            output_targets: match matches.value_of("output_targets") {
                Some(output_targets) => Some(output_targets.to_string()),
                None => None,
            },
        }
    }
}

pub fn run(logger: &mut Logger, options: &QuickWisCallOptions) -> Result<()> {
    info!(logger, "Running: cnvetti quick wis-call");
    info!(logger, "Options: {:?}", options);

    let tmp_dir = TempDir::new("cnvetti_quick_wis_call")
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

    // Compute coverage relative to WIS model.
    info!(logger, "Compute model-based coverage normalization");
    let output_targets = if let Some(ref output_targets) = options.output_targets {
        output_targets.clone()
    } else {
        tmp_dir
            .path()
            .join(format!("output_targets.bcf"))
            .to_str()
            .unwrap()
            .to_string()
    };
    lib_mod_cov::run(
        &mut logger.new(o!("step" => "mod-coverage")),
        &options.into_build_model_based_coverage_options(
            &norm_out,
            &options.input_model,
            &output_targets,
        ),
    ).chain_err(|| "Problem with merging coverage file")?;
    info!(logger, " => done");

    warn!(logger, "Actual calling step has not been implemented yet!");

    info!(logger, "All done. Have a nice day!");
    Ok(())
}
