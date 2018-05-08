// Shortcuts, preconfigured with best practice parameters.

// TODO: get rid of the logger object...

use slog::Logger;
use tempdir::TempDir;

pub mod options;
pub use self::options::*;

use cli::call::call as call_call;
use cli::coverage::call as call_coverage;
use cli::normalize::call as call_normalize;
use cli::segment::call as call_segment;

/// Main entry point for the "segment" command.
pub fn call_quick_wgs_deep(logger: &mut Logger, options: &WgsDeepOptions) -> Result<(), String> {
    let options = options.clone();

    info!(logger, "Running cnvetti quick wgs-deep");
    info!(logger, "Configuration: {:?}", &options);

    // Create the temporary directory that we will use for calling the wrapped commands.
    let tmp_dir = TempDir::new("cnvetti_quick_wgs-deep")
        .map_err(|e| format!("Could not create temporary directory: {:?}", e))?;

    // Perform the individual steps.
    call_coverage(logger, &options.to_coverage_options(&tmp_dir))?;
    call_normalize(logger, &options.to_normalize_options(&tmp_dir))?;
    call_segment(logger, &options.to_segment_options(&tmp_dir))?;
    call_call(logger, &options.to_call_options(&tmp_dir))?;

    Ok(())
}
