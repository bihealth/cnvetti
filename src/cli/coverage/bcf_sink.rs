// Code for writing

include!(concat!(env!("OUT_DIR"), "/version.rs"));

use chrono;

use std::env;

use shlex;

use std::path::Path;

use rust_htslib::bam::{self, Read};
use rust_htslib::bcf;

use bio::io::fasta;

use slog::Logger;

pub struct CoverageBcfSink {
    pub header: bcf::Header,
    pub writer: bcf::Writer,
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

pub fn sample_from_bam<P: AsRef<Path>>(path: P) -> Result<Vec<String>, String> {
    let path = path.as_ref().to_str().unwrap();
    let reader = bam::Reader::from_path(path)
        .map_err(|e| format!("Could not open BAM file for reading {}: {}", &path, e))?;

    let text = String::from_utf8(Vec::from(reader.header().as_bytes())).unwrap();
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

    Ok(samples)
}

impl CoverageBcfSink {
    pub fn from_path<P: AsRef<Path>>(
        path: P,
        samples: &Vec<String>,
        fasta_index: &fasta::Index,
        _logger: &mut Logger,
    ) -> Result<Self, String> {
        let path = path.as_ref().to_str().unwrap();

        let uncompressed = !path.ends_with(".bcf") && !path.ends_with(".vcf.gz");
        let vcf = path.ends_with(".vcf") || path.ends_with(".vcf.gz");

        let header = Self::build_header(samples, fasta_index);
        let writer = bcf::Writer::from_path(&path, &header, uncompressed, vcf)
            .map_err(|e| format!("Could not open BCF writer: {}", e))?;

        Ok(CoverageBcfSink { header, writer })
    }

    fn build_header(samples: &Vec<String>, fasta_index: &fasta::Index) -> bcf::Header {
        let mut header = bcf::Header::new();
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

        // Add samples to BCF header.
        for sample in samples {
            header.push_sample(sample.as_bytes());
        }

        // Put contig information into BCF header.
        for seq in fasta_index.sequences() {
            header.push_record(format!("##contig=<ID={},length={}>", seq.name, seq.len).as_bytes());
        }

        // TODO: GAP should be a Flag
        // TODO: Change description based on configured count.
        let lines = vec![
            // Misc fields
            "##ALT=<ID=COUNT,Description=\"Record describes a window for read counting\">",
            // TODO: FEW_GCWINDOWS should go into normalize but then cannot be found?
            "##FILTER=<ID=FEW_GCWINDOWS,Description=\"Masked because of few windows with \
             this GC content\">",
            // INFO fields describing the window
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Window end\">",
            "##INFO=<ID=GC,Number=1,Type=Float,Description=\"Reference GC content in percent\">",
            "##INFO=<ID=MAPABILITY,Number=1,Type=Float,Description=\"Mean mapability in the \
             window\">",
            "##INFO=<ID=BLACKLIST,Number=1,Type=Float,Description=\"Interval overlaps with \
             blacklist site\">",
            "##INFO=<ID=GAP,Number=1,Type=Integer,Description=\"Window overlaps with N in \
             reference (gap)\">",
            // Generic FORMAT fields
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##FORMAT=<ID=COV,Number=1,Type=Float,Description=\"Average coverage\">",
            // FORMAT fields used for coverage
            "##FORMAT=<ID=WINSD,Number=1,Type=Float,Description=\"Per-window coverage SD)\">",
            // FORMAT fields used for read counts
            "##FORMAT=<ID=MP,Description=\"Masked for sample because too much masked because of \
             piles\",Type=Integer,Number=1>",
            "##FORMAT=<ID=RCOV,Number=1,Type=Integer,Description=\"Coverage in number of \
             aligning reads before scaling for read piles\">",
            "##FORMAT=<ID=MS,Number=1,Type=Integer,Description=\"Number of bases in window \
             masked because of read piles\">",
        ];
        for line in lines {
            header.push_record(line.as_bytes());
        }

        header
    }
}
