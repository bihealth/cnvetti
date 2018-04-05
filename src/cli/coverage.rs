// Implementation of the "coverage" command.

include!(concat!(env!("OUT_DIR"), "/version.rs"));

use std::collections::HashMap;
use std::env;
use std::fs::File;

use bio::io::bed;
use bio::io::fasta;

use chrono;

use clap::{ArgMatches, Values};

use rust_htslib::bcf;
use rust_htslib::bam::{self, Read};

use shlex;

use slog::Logger;

/// Options for the "coverage" command.
#[derive(Debug)]
pub struct Options {
    reference: String,
    input: Vec<String>,
    output: String,
    genome_regions: Vec<String>,
    mapability_bed: Option<String>,
}

/// Build options from ArgMatches.
pub fn build_options(matches: &ArgMatches) -> Options {
    Options {
        reference: matches.value_of("reference").unwrap().to_string(),
        input: matches
            .values_of("input")
            .unwrap_or(Values::default())
            .map(|res| res.to_string())
            .collect(),
        output: matches.value_of("output").unwrap().to_string(),
        genome_regions: matches
            .values_of("genome-regions")
            .unwrap_or(Values::default())
            .map(|res| res.to_string())
            .collect(),
        mapability_bed: match matches.value_of("mapability-bed") {
            Some(x) => Some(x.to_string()),
            None => None,
        },
    }
}

/// Structure holding the input files readers and related data structures.
struct CoverageInput {
    ref_reader: fasta::IndexedReader<File>,
    bed_reader: Option<bed::Reader<File>>,
    bam_readers: Vec<bam::IndexedReader>,
    samples: Vec<String>,
    sample_to_index: HashMap<String, u32>,
    read_group_to_sample: HashMap<String, String>,
}

/// Parse @RG lane into triple (id, sm).
fn parse_line_rg(line: String) -> Option<(String, String)> {
    let line_split = line.split("\t");
    let mut id: Option<String> = None;
    let mut sm: Option<String> = None;
    for s in line_split {
        let token: Vec<&str> = s.split(":").collect();
        if token.len() >= 2 {
            match token[0] {
                "ID" => {
                    id = Some(token[1].to_string());
                }
                "SM" => {
                    sm = Some(token[1].to_string());
                }
                _ => (),
            }
        }
    }

    match (id, sm) {
        (Some(id), Some(sm)) => Some((id, sm)),
        _ => None,
    }
}

/// Open input file handlers.
fn open_input(logger: &Logger, options: &Options) -> CoverageInput {
    // Open BAM files.
    let bam_readers: Vec<bam::IndexedReader> = options
        .input
        .iter()
        .map(|path| match bam::IndexedReader::from_path(path) {
            Ok(reader) => reader,
            Err(error) => {
                panic!("Could create BAM reader for path {}! {:?}", path, error);
            }
        })
        .collect();
    // Build mapping from sample name (from read group) to numeric ID.
    let mut samples = Vec::new();
    let mut read_group_to_sample: HashMap<String, String> = HashMap::new();
    let mut sample_to_index: HashMap<String, u32> = HashMap::new();
    let mut idx: u32 = 0;
    for reader in &bam_readers {
        let text = String::from_utf8(Vec::from(reader.header().as_bytes())).unwrap();

        for line in text.lines() {
            if line.starts_with("@RG") {
                match parse_line_rg(line.to_string()) {
                    Some((id, sm)) => {
                        debug!(logger, "RG '{}' => '{}'", &id, &sm);
                        debug!(logger, "SM '{}' => '{}'", &sm, &idx);
                        samples.push(sm.clone());
                        read_group_to_sample.insert(id.clone(), sm.clone());
                        sample_to_index.insert(sm.clone(), idx.clone());
                        idx += 1;
                    }
                    None => (),
                }
            }
        }
    }

    CoverageInput {
        // Create FASTA reader for reference.
        ref_reader: match fasta::IndexedReader::from_file(&options.reference) {
            Ok(reader) => reader,
            Err(error) => {
                panic!("Could not open indexed FASTA File! {:?}", error);
            }
        },
        // Create reader for mapability BED file, if any.
        bed_reader: match &options.mapability_bed {
            &Some(ref path) => match bed::Reader::from_file(&path) {
                Ok(reader) => Some(reader),
                Err(error) => {
                    panic!("Could create BED reader {:?}", error);
                }
            },
            &None => None,
        },
        // Create readers for BAM files.
        bam_readers: bam_readers,
        // The samples in the order defined in input files.
        samples: samples,
        // Mapping from sample name to index.
        sample_to_index: sample_to_index,
        // Mapping from read group name to sample name.
        read_group_to_sample: read_group_to_sample,
    }
}

/// Structure holding the output files readers and related data structures.
struct CoverageOutput {
    header: bcf::Header,
    bcf_writer: bcf::Writer,
}

/// Create BCF header from input files.
fn build_header(input: CoverageInput) -> bcf::Header {
    let mut header = bcf::Header::new();

    // TODO(holtgrewe): compare contig names in BAM and FAI file?
    // TODO(holtgrewe): what about the VCF version?

    // Put overall meta information into the BCF header.
    let now = chrono::Utc::now();
    header.push_record(format!("##fileDate={}", now.format("%Y%m%d").to_string()).as_bytes());

    // Put creating tool version an dcall into file.
    header.push_record(format!("##cnvetti_coverageVersion={}", VERSION).as_bytes());
    header.push_record(
        format!(
            "##cnvetti_coverageCommand={}",
            env::args()
                .map(|s| shlex::quote(&s).to_string())
                .collect::<Vec<String>>()
                .join(" ")
        ).as_bytes(),
    );

    // Put contig information into BCF header.
    for seq in input.ref_reader.index.sequences() {
        header.push_record(format!("##contig=<ID={},length={}>", seq.name, seq.len).as_bytes());
    }

    // Add samples to BCF header.
    for sample in input.samples {
        header.push_sample(sample.as_bytes());
    }

    // Write out FILTER, INFO, FORMAT, and ALT fields.
    // TODO: some of these could go away eventually
    let lines = vec![
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Window end\">",
        "##INFO=<ID=GC,Number=1,Type=Float,Description=\"Reference GC content in percent\">",
        "##INFO=<ID=MAPABILITY,Number=1,Type=Float,Description=\"Mean mapability in the window\">",
        "##INFO=<ID=GAP,Number=0,Type=Flag,Description=\"Window overlaps with N in reference (gap)\">",
        "##INFO=<ID=GCWINDOWS,Number=1,Type=Integer,Description=\"Number of windows with same GC content\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=COV,Number=1,Type=Float,Description=\"Average coverage with non-q0 reads\">",
        "##FORMAT=<ID=COV0,Number=1,Type=Float,Description=\"Average coverage with q0+non-q0 reads\">",
        "##FORMAT=<ID=WINSD,Number=1,Type=Float,Description=\"Per-window coverage SD (non-q0 reads)\">",
        "##FORMAT=<ID=WINSD0,Number=1,Type=Float,Description=\"Per-window coverage SD (q0+non-q0 reads)\">",
        "##FORMAT=<ID=RC,Number=1,Type=Float,Description=\"Number of aligning non-q0 reads\">",
        "##FORMAT=<ID=RC0,Number=1,Type=Float,Description=\"Number of aligning q0+non-q0 reads\">",
        "##ALT=<ID=COUNT,Description=\"Record describes a window for read counting\">",
    ];
    for line in lines {
        header.push_record(line.as_bytes());
    }

    header
}

/// Open output file handlers.
fn open_output(logger: &Logger, input: CoverageInput, options: Options) -> CoverageOutput {
    let uncompressed = !options.output.ends_with(".bcf") && !options.output.ends_with(".vcf.gz");
    let vcf = options.output.ends_with(".vcf") || options.output.ends_with(".vcf.gz");

    let header = build_header(input);
    info!(logger, "Writing BCF header {:?}", &header);

    CoverageOutput {
        bcf_writer: match bcf::Writer::from_path(&options.output, &header, uncompressed, vcf) {
            Ok(writer) => writer,
            Err(error) => {
                panic!(
                    "Could not open BCF file for output {}. {:?}",
                    options.output, error
                );
            }
        },
        header: header,
    }
}

/// Main entry point for the "coverage" command.
pub fn call(logger: Logger, options: Options) -> Result<(), String> {
    info!(logger, "Running cnvetti coverage");
    info!(logger, "Configuration: {:?}", options);

    info!(logger, "Opening input files");
    let input = open_input(&logger, &options);

    info!(logger, "Parsing genomic regions");
    info!(logger, "Opening output files");
    let output = open_output(&logger, input, options);

    info!(logger, "Processing");
    info!(logger, "Writing statistics");
    info!(logger, "Done. Have a nice day!");

    Ok(())
}
