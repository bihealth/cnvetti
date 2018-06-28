/// Main library module for building coverage files from a model and BAM files.
extern crate shlex;

extern crate clap;

#[macro_use]
extern crate error_chain;

#[macro_use]
extern crate slog;
use slog::Logger;

extern crate rust_htslib;

extern crate lib_shared;
use lib_shared::bcf_utils;

/// This crate's error-related code, generated by `error-chain`.
mod errors {
    // Create the Error, ErrorKind, ResultExt, and Result types
    error_chain!{}
}

pub use errors::*;

mod options;
pub use options::*;
mod cov_pool;
mod cov_wis;

/// Main entry point for the `coverage` command.
pub fn run(logger: &mut Logger, options: &ModelBasedCoverageOptions) -> Result<()> {
    info!(logger, "Running: cnvetti cmd mod-coverage");
    info!(logger, "Options: {:?}", options);

    if let Some(input_wis_model) = &options.input_wis_model {
        // Load normalized counts and input header.
        info!(logger, "Loading normalized counts...");
        let normalized_counts = cov_wis::load_counts(&options.input)?;
        info!(logger, " => done");

        // Extract the target-wise information.
        info!(logger, "Extracting target-wise information...");
        let target_infos = cov_wis::load_target_infos(&input_wis_model, &normalized_counts)?;
        info!(logger, " => done");

        // Write output file.
        info!(logger, "Writing model-based coverage file...");
        cov_wis::write_mod_cov(&target_infos, &options.input, &options.output)?;
        info!(logger, " => done");
    } else if let Some(input_pool_model) = &options.input_pool_model {
        // Write output file.
        cov_pool::run(logger, &input_pool_model, &options.input, &options.output)?;
    } else {
        bail!("Neither WIS or pool-based model option found");
    }

    // Finally, create index on created output file.
    info!(logger, "Building index for output file...");
    bcf_utils::build_index(logger, &options.output).chain_err(|| "Could not build index")?;
    info!(logger, "All done. Have a nice day!");

    Ok(())
}
