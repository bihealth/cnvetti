use std::path::Path;
use std::str;

use regex::Regex;

use super::regions::GenomeRegions;
use rust_htslib::bam::{self, Read};

mod errors {
    // Create the Error, ErrorKind, ResultExt, and Result types
    error_chain! {}
}

pub use self::errors::*;

/// Generate list of all contigs from BAM header.
pub fn build_chroms_bam(header: &bam::HeaderView, regex: Option<String>) -> Result<GenomeRegions> {
    let mut pairs: Vec<(String, usize)> = Vec::new();

    let re = if let Some(regex) = regex {
        Some(Regex::new(&regex).unwrap())
    } else {
        None
    };

    for (no, name) in header.target_names().iter().enumerate() {
        let name = str::from_utf8(name)
            .chain_err(|| "Could not decode contig name")?
            .to_string();
        let len = header.target_len(no as u32).unwrap() as usize;
        let do_push = match &re {
            Some(re) => re.is_match(&name),
            None => true,
        };
        if do_push {
            pairs.push((name, len));
        }
    }

    Ok(GenomeRegions::from_name_and_length(&pairs))
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

/// Construct sample list from BAM header.
fn samples_from_header(header: &Vec<u8>) -> Vec<String> {
    let text = String::from_utf8(header.to_vec()).unwrap();
    let mut samples = Vec::new();

    for line in text.lines() {
        if line.starts_with("@RG") {
            match parse_line_rg(line.to_string()) {
                Some((_id, sm)) => {
                    if !samples.contains(&sm) {
                        samples.push(sm.clone());
                    }
                }
                None => (),
            }
        }
    }

    samples
}

/// Construct sample list from BAM file.
pub fn samples_from_file<P: AsRef<Path>>(path: P) -> Result<Vec<String>> {
    let path = path.as_ref().to_str().unwrap();
    let reader = bam::Reader::from_path(path)
        .chain_err(|| format!("Could not open BAM file for reading"))?;

    Ok(samples_from_header(&Vec::from(reader.header().as_bytes())))
}
