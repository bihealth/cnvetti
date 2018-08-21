//! Genotyping using the segmentation overlap approach.

use std::env;
use std::ops::Range;

use shlex;
use slog::Logger;
use chrono;
use rust_htslib::bcf::{self, Read};

use super::errors::*;
use options::*;

use lib_segment::seg_haar;
use lib_shared::bcf_utils;
use lib_shared::regions::GenomeRegions;
use lib_shared::stats::Stats;
use rust_segment::shared::{CopyState, Segment, Segmentation};

/// Write out segmentation.
fn write_result(
    filtered_segmentation: &Vec<(bool, CnvGenotypeInfo)>,
    chrom: &String,
    ranges: &Vec<Range<usize>>,
    writer: &mut bcf::Writer,
) -> Result<()> {
    let rid = writer
        .header()
        .name2rid(chrom.as_bytes())
        .chain_err(|| format!("Could not find contig {}", chrom))?;

    for (is_ok, ref cnv_info) in filtered_segmentation {
        assert!(cnv_info.quals.copy_state != CopyState::Neutral);
        let mut record = writer.empty_record();

        // Columns: CHROM, POS, ID, REF, ALT, (FILTER)
        let pos = ranges[cnv_info.segment.range.start].start;
        let end = ranges[cnv_info.segment.range.end].end;
        let alleles_v = vec![
            Vec::from("N"),
            if cnv_info.quals.copy_state == CopyState::Deletion {
                Vec::from("<DEL>")
            } else {
                Vec::from("<DUP>")
            },
        ];
        let alleles = alleles_v
            .iter()
            .map(|x| x.as_slice())
            .collect::<Vec<&[u8]>>();

        let (sv_type, sv_len) = if cnv_info.quals.copy_state == CopyState::Deletion {
            ("DEL", -((end - pos) as i32))
        } else {
            ("DUP", (end - pos) as i32)
        };
        let num_targets = cnv_info.segment.range.len() as i32;

        record.set_rid(&Some(rid));
        record.set_pos(pos as i32);
        record
            .set_id(format!("{}_{}:{}-{}", sv_type, chrom, pos + 1, end).as_bytes())
            .chain_err(|| "Could not update ID")?;
        record
            .set_alleles(&alleles)
            .chain_err(|| "Could not update alleles")?;

        // Columns: INFO
        record
            .push_info_integer(b"END", &[end as i32])
            .chain_err(|| "Could not write INFO/END")?;
        record
            .push_info_string(b"SVTYPE", &[sv_type.as_bytes()])
            .chain_err(|| "Could not write INFO/SVTYPE")?;
        record
            .push_info_integer(b"SVLEN", &[sv_len])
            .chain_err(|| "Could not write INFO/SVLEN")?;
        record
            .push_info_integer(b"TARGETS", &[num_targets])
            .chain_err(|| "Could not write INFO/TARGETS")?;
        record
            .push_info_integer(b"CENTER", &[((pos + end) / 2) as i32])
            .chain_err(|| "Could not write INFO/END")?;

        // Columns: FORMAT/GT
        // TODO: what to push as GT here?
        record
            .push_format_integer(b"GT", &[bcf::GT_MISSING])
            .chain_err(|| "Could not write FORMAT/GT")?;
        if !is_ok {
            record
                .push_format_string(b"FT", &[b"LowQual"])
                .chain_err(|| "Could not find FORMAT/FILTER LowQual")?;
        }
        record
            .push_format_integer(b"CN", &[cs_to_idx(cnv_info.quals.copy_state) as i32 + 1])
            .chain_err(|| "Could not write FORMAT/CN")?;

        record
            .push_format_float(b"EQ1", &[*cnv_info.quals.q_exact_del as f32])
            .chain_err(|| "Could not write FORMAT/EQ1")?;
        record
            .push_format_float(b"SQ1", &[*cnv_info.quals.q_some_del as f32])
            .chain_err(|| "Could not write FORMAT/SQ1")?;
        record
            .push_format_float(b"NQ1", &[*cnv_info.quals.q_no_del as f32])
            .chain_err(|| "Could not write FORMAT/NQ1")?;
        record
            .push_format_float(b"LQ1", &[*cnv_info.quals.q_left_del_bp as f32])
            .chain_err(|| "Could not write FORMAT/LQ1")?;
        record
            .push_format_float(b"RQ1", &[*cnv_info.quals.q_right_del_bp as f32])
            .chain_err(|| "Could not write FORMAT/RQ1")?;

        record
            .push_format_float(b"EQ3", &[*cnv_info.quals.q_exact_dup as f32])
            .chain_err(|| "Could not write FORMAT/EQ3")?;
        record
            .push_format_float(b"SQ3", &[*cnv_info.quals.q_some_dup as f32])
            .chain_err(|| "Could not write FORMAT/SQ3")?;
        record
            .push_format_float(b"NQ3", &[*cnv_info.quals.q_no_dup as f32])
            .chain_err(|| "Could not write FORMAT/NQ3")?;
        record
            .push_format_float(b"LQ3", &[*cnv_info.quals.q_left_dup_bp as f32])
            .chain_err(|| "Could not write FORMAT/LQ3")?;
        record
            .push_format_float(b"RQ3", &[*cnv_info.quals.q_right_dup_bp as f32])
            .chain_err(|| "Could not write FORMAT/RQ3")?;

        record
            .push_format_float(b"NDQ", &[*cnv_info.quals.q_not_diploid as f32])
            .chain_err(|| "Could not write FORMAT/NDQ")?;
        record
            .push_format_float(b"DQ", &[*cnv_info.quals.q_diploid as f32])
            .chain_err(|| "Could not write FORMAT/DQ")?;

        record
            .push_format_float(b"CV", &[cnv_info.quals.mean_cov as f32])
            .chain_err(|| "Could not write FORMAT/CV")?;
        record
            .push_format_float(b"CV2", &[cnv_info.quals.mean_cov.log2() as f32])
            .chain_err(|| "Could not write FORMAT/CV2")?;
        record
            .push_format_float(b"CVZ", &[cnv_info.quals.mean_z_score as f32])
            .chain_err(|| "Could not write FORMAT/CVZ")?;

        writer
            .write(&record)
            .chain_err(|| "Could not write BCF record")?;
    }

    Ok(())
}

/// Filter segmentation.
fn filter_segmentation(
    segmentation: &Vec<CnvGenotypeInfo>,
    _options: &GenotypeOptions,
) -> Result<Vec<(bool, CnvGenotypeInfo)>> {
    // TODO: actually implement!
    let result = segmentation
        .iter()
        .map(|x| (true, x.clone()))
        .collect::<Vec<(bool, CnvGenotypeInfo)>>();
    Ok(result)
}

/// Whether or not to skip the record.
fn skip_record(record: &mut bcf::Record) -> bool {
    if !record.format(b"CV").float().is_ok()
        || !record.format(b"CV2").float().is_ok()
        || !record.format(b"CVZ").float().is_ok()
    {
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
    let cvz = record
        .format(b"CVZ")
        .float()
        .expect("Could not access FORMAT/CVZ")[0][0] as f64;

    !cv.is_finite() || !cv2.is_finite() || !cvz.is_finite()
}

/// Read segmentation and region-wise coverage from input file(s).
fn read_seg_and_cov(
    logger: &mut Logger,
    reader: &mut bcf::IndexedReader,
    _reader_calls: Option<&mut bcf::IndexedReader>,
    options: &GenotypeOptions,
) -> Result<(Segmentation, Vec<Range<usize>>, Vec<f64>)> {
    let mut ranges = Vec::new();
    let mut covs = Vec::new();
    let mut cov2s = Vec::new();
    let mut covzs = Vec::new();
    let mut cn_states = Vec::new();

    let mut record = reader.empty_record();

    // Load values from file.
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

        covs.push(
            record
                .format(b"CV")
                .float()
                .chain_err(|| "Could not access FORMAT/CV")?[0][0] as f64,
        );
        cov2s.push(
            record
                .format(b"CV2")
                .float()
                .chain_err(|| "Could not access FORMAT/CV2")?[0][0] as f64,
        );

        let cvz = record
            .format(b"CVZ")
            .float()
            .chain_err(|| "Could not access FORMAT/CVZ")?[0][0] as f64;
        const Z_SCORE_LIMIT_FACTOR: f64 = 5.0;
        let z_score_limit = options.xhmm_z_score_threshold * Z_SCORE_LIMIT_FACTOR;
        let cvz = if cvz < -z_score_limit {
            -z_score_limit
        } else if cvz > z_score_limit {
            z_score_limit
        } else {
            cvz
        };

        covzs.push(cvz);

        let cn_state = record
            .format(b"SGS")
            .integer()
            .chain_err(|| "Could not access FORMAT/SGS")?[0][0];
        cn_states.push(match cn_state {
            0 => CopyState::Deletion,
            1 => CopyState::Neutral,
            2 => CopyState::Duplication,
            _ => bail!("Could not decode copy state"),
        });
    }

    // Extract segmentation information.
    let mut curr = None;
    let mut begin = 0;
    let mut end = 0;
    let mut segments = Vec::new();
    for i in 0..(covs.len()) {
        if let Some(state) = curr {
            end = i;
            if state != cn_states[i] {
                let mut segment = Segment {
                    range: begin..end,
                    copy_state: Some(state),
                    mean: 0.0,
                    std_dev: 0.0,
                    mean_log2: 0.0,
                    std_dev_log2: 0.0,
                };
                segment.update_stats(&covs);
                segments.push(segment);
                curr = Some(cn_states[i]);
                begin = i;
            }
        } else {
            curr = Some(cn_states[i]);
            begin = i;
            end = i;
        }
    }
    if let Some(state) = curr {
        if begin != end {
            let mut segment = Segment {
                range: begin..end,
                copy_state: Some(state),
                mean: 0.0,
                std_dev: 0.0,
                mean_log2: 0.0,
                std_dev_log2: 0.0,
            };
            segment.update_stats(&covs);
            segments.push(segment);
        }
    }

    trace!(logger, "# segments = {}", segments.len());
    trace!(logger, "# values = {}", covs.len());
    trace!(logger, "# values_log2 = {}", cov2s.len());

    // Create Segmentation struct and return.
    let segmentation = Segmentation {
        segments: segments,
        values: covs,
        values_log2: cov2s,
        cn_states: Some(cn_states),
    };

    Ok((segmentation, ranges, covzs))
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
        "##INFO=<ID=TARGETS,Number=.,Type=Integer,Description=\"Number of target regions
         used for calling\">",
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
    // TODO: check compatibility between BCF contig headers
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
        let (segmentation, ranges, covzs) =
            read_seg_and_cov(logger, &mut reader, reader_calls.as_mut(), &options)?;

        trace!(logger, "Computing quality metrics of the segmentation");
        let gt_infos = compute_seg_metrics(logger, &segmentation, &ranges, &covzs, &options)?;

        trace!(logger, "Post-filtering segments");
        let gt_infos = filter_segmentation(&gt_infos, &options)?;

        trace!(logger, "Write out calls BCF file");
        write_result(&gt_infos, &chrom, &ranges, &mut writer)?;
    }

    info!(logger, "=> OK");
    Ok(())
}
