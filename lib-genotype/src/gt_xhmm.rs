//! Genotyping using the XHMM approach.

use std::env;

use bio::stats::hmm::{viterbi, State};
use bio::stats::probs::PHREDProb;
use chrono;
use rust_htslib::bcf::{self, Read};
// use rust_segment::{reject_nonaberrant_pvalue, seg_haar};
// use separator::Separatable;
use shlex;
use slog::Logger;

use super::errors::*;
use options::*;

use lib_shared::bcf_utils;
use lib_shared::regions::GenomeRegions;
use rust_segment::shared::{CopyState, Segment, Segmentation};

/// Write out segmentation.
fn write_result(
    filtered_segmentation: &Vec<(bool, CnvGenotypeInfo)>,
    writer: &mut bcf::Writer,
) -> Result<()> {
    bail!("Implmenet me!")
}

/// Filter segmentation.
fn filter_segmentation(
    segmentation: &Vec<CnvGenotypeInfo>,
    options: &GenotypeOptions,
) -> Result<Vec<(bool, CnvGenotypeInfo)>> {
    bail!("Implement me!")
}

/// Struct with quality values and other metrics from XHMM describing an CNV.
struct XhmmCnvQuals {
    /// Probability for "exact deletion".
    qExactDel: PHREDProb,
    /// Probability for "some deletion".
    qSomeDel: PHREDProb,
    /// Probability for "no deletion".
    qNoDel: PHREDProb,
    /// Probability for "left deletion breakpoint".
    qLeftDelBp: PHREDProb,
    /// Probability for "right deletion breakpoint".
    qRightDelBp: PHREDProb,

    /// Probability for "exact duplication".
    qExactDup: PHREDProb,
    /// Probability for "some duplication".
    qSomeDup: PHREDProb,
    /// Probability for "no duplication".
    qNoDup: PHREDProb,
    /// Probability for "left duplication breakpoint".
    qLeftDupBp: PHREDProb,
    /// Probability for "right duplication breakpoint".
    qRightDupBp: PHREDProb,

    /// Probability for "is not".
    qNotDiploid: PHREDProb,
    /// Probability for "is diploid".
    qDiploid: PHREDProb,

    /// Mean normalized read depth throughout the CNV.
    meanReadDepth: f64,
    /// Mean Z-score of read depth through the CNV.
    meanZScore: f64,
    /// Number of target region in the cnv.
    numTargetRegions: u32,
}

/// Description of a CNV as called by XHMM.
pub struct CnvGenotypeInfo {
    /// Segment with region information etc.
    segment: Segment,
    /// Qualities associated with the CNV.
    quals: XhmmCnvQuals,
}

/// Compute metrics associated with each segment.
fn compute_seg_metrics(
    segmentation: &Segmentation,
    cvs: &Vec<f64>,
    covzs: &Vec<f64>,
) -> Result<Vec<CnvGenotypeInfo>> {
    bail!("Implement me!")
}

/// Read segmentation and region-wise coverage from input file(s).
fn read_seg_and_cov(
    reader: &mut bcf::IndexedReader,
    reader_calls: Option<&mut bcf::IndexedReader>,
) -> Result<(Segmentation, Vec<f64>, Vec<f64>)> {
    bail!("Implement me!")
}

/// Build header for the output BCF file.
///
/// This defines all values used throughout the whole window/target specific BAM files,
/// regardless whether they are actually used in the file.
///
/// Note that we use shared FORMAT tags for coverage and fragment count, such that we get
/// unified processing of copy number data in BCF files.
fn build_header(samples: &Vec<String>, contigs: &GenomeRegions) -> bcf::Header {
    let mut header = bcf::Header::new();

    // Put overall meta information into the BCF header.
    let now = chrono::Utc::now();
    header.push_record(format!("##fileDate={}", now.format("%Y%m%d").to_string()).as_bytes());

    // Put creating tool version and call into file.
    header.push_record(format!("##cnvetti_cmdGenotypeVersion={}", "0.1.0").as_bytes());
    header.push_record(
        format!(
            "##cnvetti_cmdGentypeCommand={}",
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
    for (name, _, length) in &contigs.regions {
        header.push_record(format!("##contig=<ID={},length={}>", name, length).as_bytes());
    }

    // Push the relevant header records.
    let lines = vec![
        // Define ALT column <DEL>/<DUP>/<CNV>
        "##ALT=<ID=DEL,Description=\"Record describes a deletion, with respect to the
         reference\">",
        "##ALT=<ID=DUP,Description=\"Record describes a duplication, with respect to the
         reference\">",
        "##ALT=<ID=CNV,Description=\"Record describes a copy number variant with respect
         to the reference, results from merging two overlapping DEL/DUP regions or from
         an inconclusive \"missing genotype\" call, e.g., using Exome HMM method\">",
        // INFO fields describing the window
        "##INFO=<ID=CENTER,Number=1,Type=Integer,Description=\"Mid-point of the CNV, to have a
         single number, e.g., for plotting.\">",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Window end\">",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV\">",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV in bp\">",
        // FILTER fields
        "##FILTER=<ID=LowQual,Description=\"Low-quality call or genotype.\">",
        // Generic FORMAT fields
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"The genotype in the sample\">",
        "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype-wise filter, \
         semicolon-separated.\">",
        // Coverage- and quality-related FORMAT fields.
        "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Predicted copy number;
         0=no-call, 1=del, 2=diploid, 3=duplication\">",
        "##FORMAT=<ID=QED,Number=1,Type=Float,Description=\"Phred-scaled probability for \
         'exact del/dup'\">",
        "##FORMAT=<ID=QSD,Number=1,Type=Float,Description=\"Phred-scaled probability for \
         'some del/dup'\">",
        "##FORMAT=<ID=QND,Number=1,Type=Float,Description=\"Phred-scaled probability for \
         'no del/dup'\">",
        "##FORMAT=<ID=QLD,Number=1,Type=Float,Description=\"Phred-scaled probability for \
         'left del/dup breakpoint'\">",
        "##FORMAT=<ID=QRD,Number=1,Type=Float,Description=\"Phred-scaled probability for \
         'right del/dup breakpoint'\">",
        "##FORMAT=<ID=QED,Number=1,Type=Float,Description=\"Phred-scaled probability for \
         'exact del/dup'\">",
        "##FORMAT=<ID=CV,Number=1,Type=Float,Description=\"Mean coverage over the CNV region\">",
        "##FORMAT=<ID=CV2,Number=1,Type=Float,Description=\"Mean log2-scaled coverage over
         the CNV region\">",
        "##FORMAT=<ID=CVZ,Number=1,Type=Float,Description=\"Mean Z-score over the CNV region\">",
    ];
    for line in lines {
        header.push_record(line.as_bytes());
    }

    header
}

/// Build bcf::Writer with appropriate header.
fn build_bcf_writer(path: &String, reader: &bcf::IndexedReader) -> Result<bcf::Writer> {
    let uncompressed = !path.ends_with(".bcf") && !path.ends_with(".vcf.gz");
    let vcf = path.ends_with(".vcf") || path.ends_with(".vcf.gz");

    let chroms = bcf_utils::extract_chroms(&reader.header());
    let samples = reader
        .header()
        .samples()
        .iter()
        .map(|s| {
            String::from_utf8(s.to_vec()).expect(&format!("Could not decode sample name: {:?}", s))
        })
        .collect::<Vec<String>>();

    let header = build_header(&samples, &chroms);
    bcf::Writer::from_path(&path, &header, uncompressed, vcf)
        .chain_err(|| "Could not open BCF file for writing")
}

/// Perform genotyping using the XHMM algorithm.
pub fn run_genotyping(logger: &mut Logger, options: &GenotypeOptions) -> Result<()> {
    info!(logger, "Computing genotyping using XHMM");

    debug!(logger, "Opening input file(s)");
    let mut reader = bcf::IndexedReader::from_path(options.input.clone())
        .chain_err(|| "Could not open input BCF file")?;
    let mut reader_calls: Option<bcf::IndexedReader> = options
        .input_calls
        .as_ref()
        .map(|ref path| {
            bcf::IndexedReader::from_path(path.clone())
                .chain_err(|| "Could not open input BCF file")
        })
        .map_or(Ok(None), |r| r.map(Some))?;

    debug!(logger, "Opening output file");
    let mut writer = build_bcf_writer(&options.output, &reader)?;
    if options.io_threads > 0 {
        writer
            .set_threads(options.io_threads as usize)
            .chain_err(|| "Could not set I/O thread count")?;
    }

    debug!(logger, "Performing the actual work.");
    let chroms = bcf_utils::extract_chroms(&reader.header());
    for (chrom, _, chrom_len) in &chroms.regions {
        debug!(logger, "Processing chrom {}", &chrom);

        trace!(logger, "Jumping in BCF file(s)");
        let rid = reader
            .header()
            .name2rid(chrom.as_bytes())
            .chain_err(|| format!("[Regions BCF] Could not translate header name {}", chrom))?;
        reader
            .fetch(rid, 0, *chrom_len as u32)
            .chain_err(|| format!("[Regions BCF] Could not fetch chromosome {}", chrom))?;
        if let Some(ref mut reader_calls) = reader_calls.as_mut() {
            let rid = reader_calls
                .header()
                .name2rid(chrom.as_bytes())
                .chain_err(|| format!("[Calls BCF] Could not translate header name {}", chrom))?;
            reader_calls
                .fetch(rid, 0, *chrom_len as u32)
                .chain_err(|| format!("[Calls BCF] Could not fetch chromosome {}", chrom))?;
        }

        trace!(logger, "Loading coverage and segmentation");
        let (segmentation, covs, covzs) = read_seg_and_cov(&mut reader, reader_calls.as_mut())?;

        trace!(logger, "Computing quality metrics of the segmentation");
        let gt_infos = compute_seg_metrics(&segmentation, &covs, &covzs)?;

        trace!(logger, "Post-filtering segments");
        let gt_infos = filter_segmentation(&gt_infos, &options)?;

        trace!(logger, "Write out calls BCF file");
        write_result(&gt_infos, &mut writer)?;
    }

    info!(logger, "=> OK");
    Ok(())
}
