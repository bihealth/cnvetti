// Code for writing

include!(concat!(env!("OUT_DIR"), "/version.rs"));

use chrono;

use std::env;
use std::path::Path;

use shlex;

use rust_htslib::bcf;

use bio::io::fasta;

use slog::Logger;

pub struct CoverageBcfSink {
    pub header: bcf::Header,
    pub writer: bcf::Writer,
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

        // Put creating tool version and call into file.
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

        // TODO: Change description based on configured count.
        let lines = vec![
            // Misc fields
            "##ALT=<ID=COUNT,Description=\"Record describes a window for read counting\">",
            // INFO fields describing the window
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Window end\">",
            "##INFO=<ID=GC,Number=1,Type=Float,Description=\"Reference GC content in percent\">",
            "##INFO=<ID=MAPABILITY,Number=1,Type=Float,Description=\"Mean mapability in the \
             window\">",
            "##INFO=<ID=BLACKLIST,Number=0,Type=Flag,Description=\"Interval overlaps with \
             blacklist site\">",
            "##INFO=<ID=GAP,Number=0,Type=Flag,Description=\"Window overlaps with N in \
             reference (gap)\">",
            // Generic FORMAT fields
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##FORMAT=<ID=COV,Number=1,Type=Float,Description=\"Mean coverage\">",
            // FORMAT fields used when calling by coverage (deep WGS)
            "##FORMAT=<ID=COVSD,Number=1,Type=Float,Description=\"Per-window coverage SD)\">",
            // FORMAT fields used for read counts (off-target WES, shallow WGS)
            // TODO: use Character here with 'Y'/'N'
            "##FORMAT=<ID=MP,Description=\"Masked for sample because too much masked because of \
             piles (Y/N)\",Type=Character,Number=1>",
        ];
        for line in lines {
            header.push_record(line.as_bytes());
        }

        header
    }
}
