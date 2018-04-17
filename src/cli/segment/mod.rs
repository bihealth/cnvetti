// Implementation of the "segment" command.

include!(concat!(env!("OUT_DIR"), "/version.rs"));

mod math;
mod options;
pub mod seg_haar; // XXX should not be pub
pub mod seg_utils; // XXX should not be pub

// use std::env;
// use std::fs::File;
// use std::io::Read as IoRead;
use std::str;

// use rust_htslib::bcf::{self, Read as BcfRead};

// use serde_json;

// use shlex;

use slog::Logger;

pub use self::options::*;
// use cli::shared;

/// Main entry point for the "segment" command.
pub fn call(logger: &mut Logger, options: &Options) -> Result<(), String> {
    let options = options.clone();

    // let mut app = CoverageApp::new(logger, options);
    info!(logger, "Running cnvetti segment");
    info!(logger, "Configuration: {:?}", &options);

    Ok(())
}
