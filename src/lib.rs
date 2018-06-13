#![feature(iterator_step_by)]
#![feature(iterator_flatten)]

extern crate bio;
extern crate chrono;
extern crate clap;
extern crate histogram;
extern crate ordered_float;
extern crate rayon;
extern crate regex;
extern crate rust_htslib;
extern crate separator;
extern crate shlex;
#[macro_use]
extern crate slog;
extern crate pdqselect;
extern crate statrs;
extern crate strum;
#[macro_use]
extern crate strum_macros;
#[macro_use]
extern crate quick_error;
extern crate tempdir;

extern crate rust_segment;

pub mod cli;
