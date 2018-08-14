//! Shared code for lib-segment.

use chrono;
use rust_htslib::bcf::{self, header};
use shlex;
use slog::Logger;
use std::env;

use lib_shared::bcf_utils;

use super::errors::*;
use rust_segment::shared::{CopyState, Segment, Segmentation};

/// Create BCF file with segment.
pub fn open_segment_file(
    path: &String,
    old_header: &bcf::header::HeaderView,
) -> Result<bcf::Writer> {
    let mut header = bcf::Header::new();

    // Put overall meta information into the BCF header.
    let now = chrono::Utc::now();
    header.push_record(format!("##fileDate={}", now.format("%Y%m%d").to_string()).as_bytes());

    // Put creating tool version and call into file.
    header.push_record(format!("##cnvetti_cmdSegmentVersion={}", "0.1.0").as_bytes());
    header.push_record(
        format!(
            "##cnvetti_cmdSegmentCommand={}",
            env::args()
                .map(|s| shlex::quote(&s).to_string())
                .collect::<Vec<String>>()
                .join(" ")
        ).as_bytes(),
    );

    // Copy over samples to BCF header.
    for sample in old_header.samples() {
        header.push_sample(sample);
    }

    // Put contig information into BCF header.
    let regions = bcf_utils::extract_chroms(&old_header);
    for (name, _, length) in &regions.regions {
        header.push_record(format!("##contig=<ID={},length={}>", name, length).as_bytes());
    }

    // Push the relevant header records.
    // TODO: later decide about commented-out lines
    let lines = vec![
        // Define ALT column <SEG>
        "##ALT=<ID=SEG,Description=\"Record describes a segment\">",
        // INFO fields describing the segment
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Window end\">",
        // FILTER fields
        // Generic FORMAT fields
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=CV,Number=1,Type=Float,Description=\"Normalized mean coverage value\">",
        "##FORMAT=<ID=CV2,Number=1,Type=Float,Description=\"Log2-scaled mean coverage value\">",
        "##FORMAT=<ID=CVSD,Number=1,Type=Float,Description=\"Normalized std dev coverage value\">",
        "##FORMAT=<ID=CV2SD,Number=1,Type=Float,Description=\"Log2-scaled std dev coverage
         value\">",
    ];
    for line in lines {
        header.push_record(line.as_bytes());
    }

    let uncompressed = !path.ends_with(".bcf") && !path.ends_with(".vcf.gz");
    let vcf = path.ends_with(".vcf") || path.ends_with(".vcf.gz");

    bcf::Writer::from_path(path.clone(), &header, uncompressed, vcf)
        .chain_err(|| "Could not open output BCF file")
}

/// Write out `Segment` to BCF file.
pub fn write_segment(
    writer: &mut bcf::Writer,
    chrom: &String,
    segment: &Segment,
    borders: &Vec<usize>,
) -> Result<()> {
    let mut record = writer.empty_record();
    let rid = writer
        .header()
        .name2rid(chrom.as_bytes())
        .chain_err(|| format!("Could not find contig {}", chrom))?;
    let sample = String::from_utf8(record.header().samples()[0].to_vec())
        .expect("Could not decode sample name from UTF-8");

    // Columns: CHROM, POS, ID, REF, ALT, (FILTER)
    let alleles_v = vec![Vec::from("N"), Vec::from("<SEG>")];
    let alleles = alleles_v
        .iter()
        .map(|x| x.as_slice())
        .collect::<Vec<&[u8]>>();

    record.set_rid(&Some(rid));
    record.set_pos(borders[segment.range.start] as i32);
    record
        .set_id(
            format!(
                "{}_{}_{}_{}",
                &sample, chrom, borders[segment.range.start], borders[segment.range.end]
            ).as_bytes(),
        )
        .chain_err(|| "Could not update ID")?;
    record
        .set_alleles(&alleles)
        .chain_err(|| "Could nto udpate alleles")?;

    // Columns: INFO
    record
        .push_info_integer(b"END", &[borders[segment.range.end] as i32])
        .chain_err(|| "Could not write INFO/END")?;

    // Columns: FORMAT/{GT,CV,CV2}
    record
        .push_format_integer(b"GT", &[1])
        .chain_err(|| "Could not write FORMAT/GT: {}")?;
    record
        .push_format_float(b"CV", &[segment.mean as f32])
        .chain_err(|| "Could not write FORMAT/CV: {}")?;
    record
        .push_format_float(b"CVSD", &[segment.std_dev as f32])
        .chain_err(|| "Could not write FORMAT/CVSD: {}")?;
    record
        .push_format_float(b"CV2", &[segment.mean_log2 as f32])
        .chain_err(|| "Could not write FORMAT/CV2: {}")?;
    record
        .push_format_float(b"CV2SD", &[segment.std_dev_log2 as f32])
        .chain_err(|| "Could not write FORMAT/CV2SD: {}")?;

    writer
        .write(&record)
        .chain_err(|| "Problem writing segment BCF record to file")
}
