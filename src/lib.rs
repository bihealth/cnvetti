#![feature(iterator_step_by)]

#[macro_use]
extern crate serde_derive;
extern crate serde_json;

extern crate bio;
extern crate chrono;
extern crate clap;
extern crate histogram;
extern crate rust_htslib;
extern crate separator;
extern crate shlex;
#[macro_use]
extern crate slog;

pub mod cli;
