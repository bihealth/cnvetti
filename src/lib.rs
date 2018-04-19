#![feature(iterator_step_by)]

extern crate serde;
#[macro_use]
extern crate serde_derive;
extern crate serde_json;

extern crate bio;
extern crate chrono;
extern crate clap;
extern crate histogram;
#[macro_use] extern crate itertools;
extern crate regex;
extern crate rust_htslib;
extern crate separator;
extern crate shlex;
#[macro_use]
extern crate slog;
extern crate statrs;

pub mod cli;
