//! Genotyping using the segmentation overlap approach.

use std::env;
use std::ops::Range;

use bio::data_structures::interval_tree::IntervalTree;
use chrono;
use rust_htslib::bcf::{self, Read};
use rust_segment::seg_haar;
use separator::Separatable;
use shlex;
use slog::Logger;

use super::errors::*;
use options::*;

use lib_shared::bcf_utils;
use lib_shared::regions::GenomeRegions;
// use lib_shared::stats::Stats;
use rust_segment::shared::{CopyState, Segment, Segmentation};

use ranges::RangeSet;

/// Compute metrics associated with each segment.
fn compute_seg_metrics(
    sites: &Vec<Range<usize>>,
    segmentation: &Segmentation,
    ranges: &Vec<Range<usize>>,
    options: &GenotypeOptions,
) -> Result<Segmentation> {
    let tree = sites
        .iter()
        .enumerate()
        .map(|(i, site)| (site, i))
        .collect::<IntervalTree<usize, usize>>();
    let mut uncovered = sites
        .iter()
        .map(|site| RangeSet::from_range(site.clone()))
        .collect::<Vec<RangeSet<_>>>();

    for segment in &segmentation.segments {
        let start = ranges[segment.range.start].start;
        let end = ranges[segment.range.end - 1].end;
        for entry in tree.find(start..end) {
            let idx = *entry.data();
            let res = uncovered[idx].sub(&(start..end));
            uncovered[idx] = res;
        }
    }

    let range_tree = ranges
        .iter()
        .enumerate()
        .map(|(i, range)| (range, i))
        .collect::<IntervalTree<usize, usize>>();

    let mut segments = Vec::new();
    let mut cn_states = Vec::new();
    for (i, site) in sites.iter().enumerate() {
        let site_len = site.len() as f64;
        let uncov_len = uncovered[i]
            .entries
            .iter()
            .map(|ref e| e.end - e.start)
            .sum::<usize>() as f64;
        let overlap = (site_len - uncov_len) / site_len;
        // println!(
        //     "site = {:?}, uncovered = {:?}, overlap = {}",
        //     &site, &uncovered, overlap
        // );
        let genotyped = overlap >= options.overlap;

        let start = *range_tree.find(site).map(|it| it.data()).min().unwrap();
        let end = *range_tree.find(site).map(|it| it.data()).max().unwrap();

        let mut segment = Segment {
            range: start..end,
            copy_state: Some(CopyState::Neutral),
            mean: 0.0,
            std_dev: 0.0,
            mean_log2: 0.0,
            std_dev_log2: 0.0,
        };
        segment.update_stats(&segmentation.values);
        if genotyped {
            // TODO: This really is chromosome dependent
            if 2_f64.powf(segment.mean_log2).round() > 2.0 {
                segment.copy_state = Some(CopyState::Duplication);
            } else if 2_f64.powf(segment.mean_log2).round() < 2.0 {
                segment.copy_state = Some(CopyState::Deletion);
            }
        }

        cn_states.push(segment.copy_state.unwrap());
        segments.push(segment);
    }

    Ok(Segmentation {
        segments: segments,
        values: segmentation.values.clone(),
        values_log2: segmentation.values_log2.clone(),
        cn_states: Some(cn_states),
    })
}

/// Write out segmentation.
fn write_result(
    logger: &mut Logger,
    segmentation: &Segmentation,
    chrom: &String,
    // ranges: &Vec<Range<usize>>,
    reader: &mut bcf::IndexedReader,
    writer: &mut bcf::Writer,
) -> Result<()> {
    let rid = writer
        .header()
        .name2rid(chrom.as_bytes())
        .chain_err(|| format!("Could not find contig {}", chrom))?;

    let mut in_record = reader.empty_record();
    let mut written = 0;
    for entry in &segmentation.segments {
        match reader.read(&mut in_record) {
            Ok(_) => in_record.unpack(),
            Err(bcf::ReadError::NoMoreRecord) => bail!("Fewer calls than genotyped?"),
            _ => bail!("Could not read BCF record"),
        }

        let alleles_v = vec![Vec::from("N"), Vec::from("<CNV>")];
        let alleles = alleles_v
            .iter()
            .map(|x| x.as_slice())
            .collect::<Vec<&[u8]>>();

        let mut record = writer.empty_record();
        record.set_rid(&Some(rid));
        record.set_pos(in_record.pos() as i32);
        record
            .set_id(&in_record.id())
            .chain_err(|| "Could not set ID")?;
        record
            .set_alleles(&alleles)
            .chain_err(|| "Could not update alleles")?;
        let end = in_record
            .info(b"END")
            .integer()
            .chain_err(|| "Could not acess INFO/END")?
            .expect("Could not extract INFO/END")[0] as i32;
        record
            .push_info_integer(b"END", &[end])
            .chain_err(|| "Could not write INFO/END")?;

        // Columns: FORMAT/GT
        // TODO: what to push as GT here?
        let gt = if entry.copy_state.unwrap() == CopyState::Neutral {
            0
        } else {
            1
        };
        record
            .push_format_integer(b"GT", &[(gt + 1) << 1])
            .chain_err(|| "Could not write FORMAT/GT")?;
        record
            .push_format_integer(b"CN", &[2_f64.powf(entry.mean_log2).round() as i32])
            .chain_err(|| "Could not write FORMAT/CN")?;
        record
            .push_format_float(b"CV", &[entry.mean as f32])
            .chain_err(|| "Could not write FORMAT/CV")?;
        record
            .push_format_float(b"CV2", &[entry.mean_log2 as f32])
            .chain_err(|| "Could not write FORMAT/CV2")?;

        writer
            .write(&record)
            .chain_err(|| "Could not write BCF record")?;
        written += 1;
    }

    debug!(logger, "Wrote {} record", written);

    Ok(())
}

/// Whether or not to skip the record.
fn skip_record(record: &mut bcf::Record) -> bool {
    if !record.format(b"CV").float().is_ok() || !record.format(b"CV2").float().is_ok() {
        return true;
    }

    let cv = record
        .format(b"CV")
        .float()
        .expect("Could not access FORMAT/CV")[0][0] as f64;
    let cv2 = record
        .format(b"CV2")
        .float()
        .expect("Could not access FORMAT/CV2")[0][0] as f64;

    !cv.is_finite() || !cv2.is_finite()
}

/// Load site list as vector of ranges.
fn read_site_list(reader: &mut bcf::IndexedReader) -> Result<Vec<Range<usize>>> {
    // Read segments from segmentation, also build interval tree.
    let mut ranges = Vec::new();
    let mut record = reader.empty_record();
    loop {
        match reader.read(&mut record) {
            Ok(_) => record.unpack(),
            Err(bcf::ReadError::NoMoreRecord) => break,
            _ => bail!("Could not read BCF record"),
        }

        let pos = record.pos() as usize;
        let end = record
            .info(b"END")
            .integer()
            .chain_err(|| "Could not acess INFO/END")?
            .expect("Could not extract INFO/END")[0] as usize;

        ranges.push(pos..end);
    }

    Ok(ranges)
}

/// Load coverages and window starts/ends.
fn read_coverages(reader: &mut bcf::IndexedReader) -> Result<(Vec<f64>, Vec<Range<usize>>)> {
    // Read segments from segmentation, also build interval tree.
    let mut cov2s = Vec::new();
    let mut ranges = Vec::new();
    let mut record = reader.empty_record();
    loop {
        match reader.read(&mut record) {
            Ok(_) => record.unpack(),
            Err(bcf::ReadError::NoMoreRecord) => break,
            _ => bail!("Could not read BCF record"),
        }

        if skip_record(&mut record) {
            continue;
        }

        let pos = record.pos() as usize;
        let end = record
            .info(b"END")
            .integer()
            .chain_err(|| "Could not acess INFO/END")?
            .expect("Could not extract INFO/END")[0] as usize;
        ranges.push(pos..end);

        let cov2 = record
            .format(b"CV2")
            .float()
            .chain_err(|| "Could not access FORMAT/CV2")?[0][0] as f64;
        cov2s.push(cov2);
    }

    Ok((cov2s, ranges))
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
    header.push_record(format!("##cnvetti_cmdGenotypeVersion={}", "0.2.0").as_bytes());
    header.push_record(
        format!(
            "##cnvetti_cmdGenotypeCommand={}",
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
        "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Length of the SV in bp\">",
        // FILTER fields
        "##FILTER=<ID=LowQual,Description=\"Low-quality call or genotype.\">",
        // Generic FORMAT fields
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"The genotype in the sample\">",
        "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype-wise filter, \
         semicolon-separated.\">",
        // Coverage- and quality-related FORMAT fields.
        "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Predicted copy number, \
         2 is diploid\">",
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

/// Perform genotyping using the SegOverlap algorithm.
pub fn run_genotyping(logger: &mut Logger, options: &GenotypeOptions) -> Result<()> {
    info!(logger, "Computing genotyping using SegOverlap");

    debug!(logger, "Opening input file(s)");
    let mut reader = bcf::IndexedReader::from_path(options.input.clone())
        .chain_err(|| "Could not open input BCF file")?;
    let mut reader_calls = if let Some(ref path) = options.input_calls.as_ref() {
        bcf::IndexedReader::from_path(path.clone()).chain_err(|| "Could not open input BCF file")
    } else {
        bail!("Calls BCF has to be given when using SegOverlap genotyping")
    }?;

    debug!(logger, "Opening output file");
    let mut writer = build_bcf_writer(&options.output, &reader)?;
    if options.io_threads > 0 {
        writer
            .set_threads(options.io_threads as usize)
            .chain_err(|| "Could not set I/O thread count")?;
    }

    debug!(logger, "Performing the actual work.");
    // TODO: check compatibility between BCF contig headers
    let chroms = bcf_utils::extract_chroms(&reader.header());
    for (chrom, _, chrom_len) in &chroms.regions {
        debug!(logger, "Processing chrom {}", &chrom);

        trace!(logger, "Loading coverage BCF");
        let rid = reader
            .header()
            .name2rid(chrom.as_bytes())
            .chain_err(|| format!("[Regions BCF] Could not translate header name {}", chrom))?;
        reader
            .fetch(rid, 0, *chrom_len as u32)
            .chain_err(|| format!("[Regions BCF] Could not fetch chromosome {}", chrom))?;
        let (cov2s, ranges) = read_coverages(&mut reader)?;

        trace!(logger, "Loading site list BCF");
        let rid_calls = reader_calls
            .header()
            .name2rid(chrom.as_bytes())
            .chain_err(|| format!("[Calls BCF] Could not translate header name {}", chrom))?;
        reader_calls
            .fetch(rid_calls, 0, *chrom_len as u32)
            .chain_err(|| format!("[Calls BCF] Could not fetch chromosome {}", chrom))?;
        let sites = read_site_list(&mut reader_calls)?;

        trace!(logger, "Re-segmenting");
        trace!(logger, "# samples: {:?}", cov2s.len());
        let segmentation = seg_haar(
            &cov2s,
            None,
            None,
            &[0..(cov2s.len())],
            options.haar_seg_fdr,
            options.haar_seg_l_min,
            options.haar_seg_l_max,
        );
        debug!(
            logger,
            "Raw segments: {}",
            segmentation.segments.len().separated_string()
        );

        trace!(logger, "Perform genotyping by segment interval overlapping");
        let gt_infos = compute_seg_metrics(&sites, &segmentation, &ranges, &options)?;

        trace!(logger, "Write out calls BCF file");
        reader_calls
            .fetch(rid_calls, 0, *chrom_len as u32)
            .chain_err(|| format!("[Calls BCF] Could not fetch chromosome {}", chrom))?;
        write_result(logger, &gt_infos, &chrom, &mut reader_calls, &mut writer)?;
    }

    info!(logger, "=> OK");
    Ok(())
}
