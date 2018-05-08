#![feature(iterator_step_by)]

extern crate bio;
extern crate chrono;
extern crate clap;
extern crate histogram;
extern crate regex;
extern crate rust_htslib;
extern crate separator;
extern crate shlex;
#[macro_use]
extern crate slog;
extern crate statrs;
#[macro_use]
extern crate quick_error;
extern crate tempdir;

extern crate rust_segment;

pub mod cli;
