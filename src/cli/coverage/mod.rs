// Implementation of the "coverage" command.

use clap::{ArgMatches, Values};
use slog::Logger;

/// Options for the "coverage" command.
#[derive(Debug)]
pub struct Options {
    input: String,
    output: String,
    genome_regions: Vec<String>,
    mapability_bed: String,
}

/// Build options from ArgMatches.
pub fn build_options(matches: &ArgMatches) -> Options {
    Options {
        input: matches.value_of("input").unwrap().to_string(),
        output: matches.value_of("output").unwrap().to_string(),
        genome_regions: matches
            .values_of("genome-regions")
            .unwrap_or(Values::default())
            .map(|res| res.to_string())
            .collect(),
        mapability_bed: matches.value_of("mapability-bed").unwrap_or("").to_string(),
    }
}

/// Main entry point for the "coverage" command.
pub fn call(logger: Logger, options: Options) -> Result<(), String> {
    println!("{:?}", options);

    Ok(())
}
