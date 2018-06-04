// Implementation of the WISExome method.

use slog::Logger;

// TODO: get rid of the logger object...

pub mod options;
pub use self::options::*;

/// Main entry point for the "wisexome count" command.
pub fn call_wise_count(logger: &mut Logger, options: &CountOptions) -> Result<(), String> {
    Ok(())
}

/// Main entry point for the "wisexome normalize" command.
pub fn call_wise_normalize(logger: &mut Logger, options: &NormalizeOptions) -> Result<(), String> {
    Ok(())
}

/// Main entry point for the "wisexome build_ref" command.
pub fn call_wise_build_ref(logger: &mut Logger, options: &BuildRefOptions) -> Result<(), String> {
    Ok(())
}
