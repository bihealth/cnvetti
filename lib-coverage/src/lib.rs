#[macro_use]
extern crate error_chain;

#[macro_use]
extern crate clap;

#[macro_use]
extern crate slog;

use slog::Logger;

extern crate strum;
#[macro_use]
extern crate strum_macros;

mod options;

pub use options::*;

mod errors {
    // Create the Error, ErrorKind, ResultExt, and Result types
    error_chain!{}
}

pub use errors::*;

pub fn run(logger: &mut Logger, options: &CoverageOptions) -> Result<()> {
    info!(logger, "Running: cnvetti cmd coverage");
    info!(logger, "Options: {:?}", options);

    Ok(())
}
