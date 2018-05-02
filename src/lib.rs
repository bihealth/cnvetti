#![feature(iterator_step_by)]

#[macro_use]
extern crate serde_derive;
// #[macro_use]
extern crate serde_json;

#[macro_use]
extern crate quick_error;

extern crate bio;
extern crate chrono;
extern crate clap;
extern crate handlebars;
extern crate histogram;
extern crate regex;
extern crate rust_htslib;
extern crate separator;
extern crate shlex;
#[macro_use]
extern crate slog;
extern crate statrs;

pub mod cli;
