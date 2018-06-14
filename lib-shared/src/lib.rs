/// Module with shared code.
extern crate bio;
extern crate rust_htslib;

#[macro_use]
extern crate error_chain;

#[macro_use]
extern crate slog;

extern crate regex;

pub mod bam_utils;
pub mod bcf_utils;
pub mod regions;
pub mod stats;
