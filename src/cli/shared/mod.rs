pub mod regions;
pub mod stats;

use std::path::Path;

use slog::Logger;

use std::str;

use rust_htslib::bam::{self, Read as BamRead};
use rust_htslib::bcf::{self, Read as BcfRead};

/// Build index file for the VCF/BCF file at `path`.
pub fn build_index(logger: &mut Logger, path: &String) {
    if path.ends_with(".vcf.gz") {
        info!(logger, "Writing .tbi index file...");
        // TODO: the unsafe stuff should go into rust_htslib
        use rust_htslib;
        use std::ffi;
        let res = unsafe {
            rust_htslib::htslib::tbx_index_build(
                ffi::CString::new(path.as_bytes()).unwrap().as_ptr(),
                0,
                &rust_htslib::htslib::tbx_conf_vcf,
            )
        };
        if res != 0 {
            panic!("Could not create .tbi index");
        }
    } else if path.ends_with(".bcf") {
        info!(logger, "Writing .csi index file...");
        use rust_htslib;
        use std::ffi;
        let res = unsafe {
            rust_htslib::htslib::bcf_index_build(
                ffi::CString::new(path.as_bytes()).unwrap().as_ptr(),
                14,
            )
        };
        if res != 0 {
            panic!("Could not create .tbi index");
        }
    } else {
        info!(
            logger,
            "Not building index, output file does not having ending .bcf or .vcf.gz"
        );
    }
}

/// Tokenize genome region strings with 1-based positions.
pub fn tokenize_genome_regions(regions: &Vec<String>) -> Vec<(String, u64, u64)> {
    regions
        .iter()
        .map(|region| {
            let region_split = region.split(":").collect::<Vec<&str>>();
            if region_split.len() != 2 {
                panic!("Invalid region: {}", region);
            }

            let number_split = region_split[1].split("-").collect::<Vec<&str>>();
            if number_split.len() != 2 {
                panic!("Invalid region: {}", region);
            }

            let begin = number_split[0]
                .to_string()
                .replace(",", "")
                .parse::<u64>()
                .unwrap();
            let end = number_split[1]
                .to_string()
                .replace(",", "")
                .parse::<u64>()
                .unwrap();

            (region_split[0].to_string(), begin, end)
        })
        .collect()
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

/// Generate list of all contigs from BAM header.
pub fn build_chroms_bam(header: &bam::HeaderView) -> Vec<(String, u32)> {
    let mut result = Vec::new();

    for (no, name) in header.target_names().iter().enumerate() {
        let len = header.target_len(no as u32).unwrap();
        result.push((
            str::from_utf8(name)
                .expect("Could not decode contig name")
                .to_string(),
            len,
        ));
    }

    result
}

/// Construct sample list from BAM header.
pub fn sample_from_header(header: &Vec<u8>) -> Vec<String> {
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
pub fn sample_from_bam<P: AsRef<Path>>(path: P) -> Result<Vec<String>, String> {
    let path = path.as_ref().to_str().unwrap();
    let reader = bam::Reader::from_path(path)
        .map_err(|e| format!("Could not open BAM file for reading {}: {}", &path, e))?;

    Ok(sample_from_header(&Vec::from(reader.header().as_bytes())))
}

/// Generate list of all contigs from BCF header.
pub fn build_chroms(header: &bcf::header::HeaderView) -> Vec<(String, u32)> {
    let mut result = Vec::new();

    for ref record in header.header_records() {
        if let bcf::HeaderRecord::Contig { key, values } = record {
            assert_eq!(key, "contig");
            let mut name: Option<String> = None;
            let mut length: Option<u32> = None;
            for (ref key, ref value) in values {
                if key.as_str() == "ID" {
                    name = Some(value.to_string());
                } else if key.as_str() == "length" {
                    length = Some(value.parse::<u32>().unwrap());
                }
            }
            if let (Some(ref name), Some(length)) = (&name, length) {
                result.push((name.clone(), length));
            } else {
                panic!(
                    "Could not parse both name/length from {:?}/{:?}",
                    name, length
                );
            }
        }
    }

    result
}

pub fn process<Opt>(
    reader: &mut bcf::IndexedReader,
    writer: &mut bcf::Writer,
    logger: &mut Logger,
    options: &Opt,
    process_region: &Fn(
        &mut bcf::IndexedReader,
        &mut bcf::Writer,
        &mut Logger,
        &Opt,
        &(String, u32, u32),
    ) -> (),
) {
    let contigs = build_chroms(reader.header());
    for (chrom, length) in contigs.iter() {
        process_region(
            reader,
            writer,
            logger,
            options,
            &(chrom.to_string(), 0, *length),
        );
    }
}
