mod options;

use slog::Logger;

pub use self::options::*;


/// Main entry point for the "normalize" command.
pub fn call(logger: &mut Logger, options: &Options) -> Result<(), String> {
    // let mut app = CoverageApp::new(logger, options);
    info!(logger, "Running normalization...");

    Ok(())
}
