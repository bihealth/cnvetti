//! Options for the genotyping command.

use std::str::FromStr;

use clap::ArgMatches;

/// Define the method for genotyping.
#[derive(Clone, Copy, EnumString, Debug, PartialEq)]
pub enum GenotypingMethod {
    /// Genotyping using the exome HMM algorithm.
    ExomeHiddenMarkovModel,
}

/// Options for "cnvetti cmd genotype".
#[derive(Clone, Debug)]
pub struct GenotypeOptions {
    /// Path to input segmentation or coverage BCF file.
    pub input: String,
    /// Path to input call BCF file.
    pub input_calls: Option<String>,
    /// Path to output BCF file.
    pub output: String,
    /// Number of additional threads to use for (de-)compression in I/O.
    pub io_threads: u32,

    /// The genotyping method to employ.
    pub genotyping: GenotypingMethod,

    // Parameters from XHMM.
    /// Z-score threshold to use.
    pub xhmm_z_score_threshold: f64,
    /// Expected exome-wide CNV rate.
    pub xhmm_cnv_rate: f64,
    /// Expected mean number of targets.
    pub xhmm_mean_target_count: f64,
    /// Expected mean target distance in a CNV.
    pub xhmm_mean_target_dist: f64,
}

impl GenotypeOptions {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Self {
        println!("{:?}", matches);
        let genotyping = matches.value_of("genotyping").unwrap();
        let genotyping =
            GenotypingMethod::from_str(&genotyping).expect("Unknown genotyping");

        Self {
            input: matches.value_of("input").unwrap().to_string(),
            input_calls: matches.value_of("input_calls").map(|s| s.to_string()),
            output: matches.value_of("output").unwrap().to_string(),
            io_threads: matches
                .value_of("io_threads")
                .unwrap()
                .parse::<u32>()
                .unwrap(),

            genotyping: genotyping,

            xhmm_z_score_threshold: matches
                .value_of("xhmm_z_score_threshold")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
            xhmm_cnv_rate: matches
                .value_of("xhmm_cnv_rate")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
            xhmm_mean_target_count: matches
                .value_of("xhmm_mean_target_count")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
            xhmm_mean_target_dist: matches
                .value_of("xhmm_mean_target_dist")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
        }
    }
}
