// Implementation of the "normalize" command.

mod binned_gc;
mod options;
mod shared;

use slog::Logger;

use self::binned_gc::normalize_binned_gc;

pub use self::options::*;
use cli::shared::build_index;

// TODO: check input file.
// TODO: use index-based readers, is nicer for progress display...
// TODO: could normalize already in coverage step by writing raw coverage to temporary file

/// Main entry point for Lowess-based normalization.
pub fn call_lowess(logger: &mut Logger, _options: &Options) -> Result<(), String> {
    info!(logger, "Normalizing by LOWESS fit...");

    Ok(())
}

/// Main entry point for the "normalize" command.
pub fn call(logger: &mut Logger, options: &Options) -> Result<(), String> {
    // let mut app = CoverageApp::new(logger, options);
    info!(logger, "Running cnvetti normalize");
    info!(logger, "Configuration: {:?}", &options);

    // TODO: check availability of R.
    match &options.normalization {
        Normalization::BinnedGc => {
            normalize_binned_gc(logger, &options)?;
        }
        Normalization::LowessGc | Normalization::LowessGcMapability => {
            call_lowess(logger, &options)?;
        }
    }

    // Build index on the output file.
    build_index(logger, &options.output);

    info!(logger, "Done. Have a nice day!");

    Ok(())
}
