// Implementation of the WISExome method.
include!(concat!(env!("OUT_DIR"), "/version.rs"));

use chrono;
use regex::Regex;
use shlex;
use slog::Logger;
use std::cell::RefCell;
use std::cmp::min;
use std::collections::HashMap;
use std::env;
use std::f32;
use std::ops::Range;
use std::str;
use std::string::ToString;
use std::sync::atomic::{AtomicUsize, Ordering};

// TODO: get rid of the logger object...

pub mod options;
pub use self::options::*;

use bio::data_structures::interval_tree::IntervalTree;

use ordered_float::OrderedFloat;
use pdqselect;

use rust_htslib::bam::{self, Read as BamRead};
use rust_htslib::bcf::{self, Read as BcfRead};
use rust_htslib::tbx::{self, Read as TbxRead};

use cli::shared;
use cli::shared::regions::GenomeRegions;
use cli::shared::stats::Stats;
use cli::shared::{build_chroms_bam, sample_from_header};

use separator::Separatable;

use rayon::prelude::*;

/// Construct BCF Header from input BAM File.
fn wise_count_build_header(bam_reader: &bam::IndexedReader) -> bcf::Header {
    let mut header = bcf::Header::new();
    // TODO(holtgrewe): what about the VCF version?

    // Put overall meta information into the BCF header.
    let now = chrono::Utc::now();
    header.push_record(format!("##fileDate={}", now.format("%Y%m%d").to_string()).as_bytes());

    // Put creating tool version and call into file.
    header.push_record(format!("##cnvetti_wisexomeCountVersion={}", VERSION).as_bytes());
    header.push_record(
        format!(
            "##cnvetti_coverageCommand={}",
            env::args()
                .map(|s| shlex::quote(&s).to_string())
                .collect::<Vec<String>>()
                .join(" ")
        ).as_bytes(),
    );

    let chrom_lens = build_chroms_bam(&bam_reader.header());
    for (chrom, len) in chrom_lens {
        header.push_record(format!("##contig=<ID={},length={}>", &chrom, len).as_bytes());
    }

    let lines = vec![
        // Misc fields
        "##ALT=<ID=TARGET,Description=\"Record describes coverage of a target region\">",
        // INFO fields describing the window
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Window end\">",
        // FORMAT fields
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=FC,Number=1,Type=Float,Description=\"(Unnormalized) fragment count\">",
        // FILTER fields
        "##FILTER=<ID=LOW_COV,Description=\"Coverage is too low\">",
    ];
    for line in lines {
        header.push_record(line.as_bytes());
    }

    // Add samples to BCF header.
    let samples = sample_from_header(&Vec::from(bam_reader.header().as_bytes()));
    for sample in samples {
        header.push_sample(sample.as_bytes());
    }

    header
}

/// Load genome regions from file.
fn load_regions(
    logger: &mut Logger,
    tbx_reader: &mut tbx::Reader,
    chrom: &str,
    len: u32,
) -> Result<GenomeRegions, String> {
    let tid = match tbx_reader.tid(&chrom) {
        Ok(tid) => tid,
        Err(e) => {
            debug!(logger, "Could not map contig to ID: {}: {}", &chrom, e);
            return Ok(GenomeRegions::new());
        }
    };
    tbx_reader
        .fetch(tid, 0, len)
        .map_err(|e| format!("Could not fetch region {}:1-{}: {}", &chrom, len, e))?;
    let mut result = GenomeRegions::new();

    for buf in tbx_reader.records() {
        let s = String::from_utf8(buf.unwrap()).unwrap();
        let arr: Vec<&str> = s.split('\t').collect();
        if arr.len() < 3 {
            return Err(format!(
                "Targets BED file had too few columns: {} (< 3)",
                arr.len()
            ));
        }
        let begin = arr[1].parse::<usize>().unwrap();
        let end = arr[2].parse::<usize>().unwrap();
        result.regions.push((chrom.to_string(), begin, end));
    }

    Ok(result)
}

/// Compute distance between point and interval.
fn compute_dist(itv: &Range<i32>, pt: i32) -> i32 {
    let start = itv.start;
    let end = itv.end;

    if pt < start {
        start - pt
    } else if pt >= end {
        pt - end + 1
    } else {
        0
    }
}

/// Count reads/pairs for each target.
fn do_counting(
    logger: &mut Logger,
    regions: &GenomeRegions,
    bam_reader: &mut bam::IndexedReader,
    chrom: &str,
    len: u32,
    options: &CountOptions,
) -> Result<Vec<u32>, String> {
    let mut result = vec![0; regions.regions.len()];

    // Build interval tree.
    let mut intervals: IntervalTree<i32, usize> = IntervalTree::new();
    for (i, (_, start, end)) in regions.regions.iter().enumerate() {
        let start = *start as i32;
        let end = *end as i32;
        intervals.insert(start..end, i);
    }

    // Stream over BAM file and assign counts to intervals.
    let tid = match bam_reader.header().tid(chrom.as_bytes()) {
        Some(tid) => tid,
        None => {
            info!(logger, "Could not jump to chrom {}", chrom);
            return Ok(result);
        }
    };
    bam_reader
        .fetch(tid, 0, len)
        .map_err(|e| format!("Could not fetch {}:1-{}: {}", chrom, len, e))?;
    let mut record = bam::Record::new();
    loop {
        // First, short-circuit in case of errors, end of file, or when we want to skip this
        // particular record.s
        match bam_reader.read(&mut record) {
            Err(bam::ReadError::NoMoreRecord) => break,
            Err(err) => return Err(format!("Problem reading error: {}", err)),
            _ => {
                if record.is_paired() && !record.is_proper_pair() {
                    continue; // skip; paired but not properly
                } else if record.is_paired() && record.insert_size() <= 0 {
                    continue; // skip; paired but not left mate
                } else if record.is_secondary()
                    || record.is_quality_check_failed()
                    || record.is_supplementary()
                    || record.is_duplicate()
                    || record.mapq() < options.min_mapq
                {
                    continue; // skip; QC-related failure
                }
            }
        };

        // Next, get center of fragment.
        let begin_pos = record.pos();
        let end_pos = if record.is_paired() {
            begin_pos + record.insert_size()
        } else {
            record
                .cigar()
                .end_pos()
                .map_err(|e| format!("Problem decoding CIGAR: {}", e))?
        };
        let center_pos = (end_pos + begin_pos) / 2;

        // Query for overlapping regions.
        let mut best = None;
        for it in intervals.find(begin_pos..end_pos) {
            let itv = it.interval();
            let idx = *it.data();
            best = match best {
                None => Some((idx, compute_dist(&**itv, center_pos))),
                Some((best_idx, best_dist)) => {
                    let this_dist = compute_dist(&**itv, center_pos);
                    if this_dist < best_dist {
                        Some((idx, this_dist))
                    } else {
                        Some((best_idx, best_dist))
                    }
                }
            };
        }
        if let Some(best) = best {
            result[best.0] += 1;
        }
    }

    Ok(result)
}

/// Perform counting for "wisexome count" command.
fn wise_count_count(
    logger: &mut Logger,
    tbx_reader: &mut tbx::Reader,
    bam_reader: &mut bam::IndexedReader,
    writer: &mut bcf::Writer,
    chrom: &str,
    len: u32,
    options: &CountOptions,
) -> Result<(), String> {
    debug!(logger, "Loading targets...");
    let targets = load_regions(logger, tbx_reader, &chrom, len)?;
    debug!(
        logger,
        "=> Loaded {} targets",
        targets.regions.len().separated_string()
    );

    debug!(logger, "Assigning counts to targets...");
    let counts = do_counting(logger, &targets, bam_reader, &chrom, len, options)?;
    debug!(
        logger,
        "=> assigned {} fragments",
        counts.iter().sum::<u32>().separated_string()
    );

    debug!(logger, "Writing target-wise counts...");
    for i in 0..targets.regions.len() {
        let mut record = writer.empty_record();

        let rid = writer
            .header()
            .name2rid(chrom.as_bytes())
            .map_err(|e| format!("Could not find contig {}: {}", chrom, e))?;
        record.set_rid(&Some(rid));

        record.set_pos(targets.regions[i].1 as i32);
        record
            .set_id(
                format!(
                    "{}:{}-{}",
                    chrom,
                    targets.regions[i].1 as i32 + 1,
                    targets.regions[i].2 as i32
                ).as_bytes(),
            )
            .map_err(|e| format!("Could not write ID: {}", e))?;

        let alleles_v = vec![Vec::from("N"), Vec::from("<TARGET>")];
        let alleles = alleles_v
            .iter()
            .map(|x| x.as_slice())
            .collect::<Vec<&[u8]>>();
        record
            .set_alleles(&alleles)
            .map_err(|e| format!("Could not update alleles: {}", e))?;

        record
            .push_info_integer(b"END", &[targets.regions[i].2 as i32])
            .map_err(|e| format!("Could not write INFO/END: {}", e))?;

        record
            .push_format_integer(b"GT", &[bcf::GT_MISSING])
            .map_err(|e| format!("Could not write FORMAT/GT: {}", e))?;

        record
            .push_format_float(b"FC", &[counts[i] as f32])
            .map_err(|e| format!("Could not write FORMAT/RP: {}", e))?;

        writer
            .write(&record)
            .map_err(|e| format!("Could not write BCF record: {}", e))?;
    }
    debug!(
        logger,
        "=> Wrote {} targets",
        targets.regions.len().separated_string()
    );

    Ok(())
}

/// Main entry point for the "wisexome count" command.
pub fn call_wise_count(logger: &mut Logger, options: &CountOptions) -> Result<(), String> {
    let options = options.clone();

    info!(logger, "Running cnvetti wisexome count");
    info!(logger, "Configuration: {:?}", &options);

    debug!(logger, "Opening input targets BED file");
    let mut tbx_reader = tbx::Reader::from_path(&options.targets_bed).map_err(|x| x.to_string())?;
    if options.io_threads > 0 {
        tbx_reader
            .set_threads(options.io_threads as usize)
            .expect("Could not set I/O thread count");
    }

    debug!(logger, "Opening input BAM file");
    let mut bam_reader = bam::IndexedReader::from_path(&options.input).map_err(|e| e.to_string())?;
    // TODO: set threads?

    // Generate chromsomes and lengths from input BAM file.
    let chrom_lens = build_chroms_bam(&bam_reader.header());
    // Construct BCF writer (and header) from input BAM file.
    let header = wise_count_build_header(&bam_reader);
    let uncompressed = !options.output.ends_with(".bcf") && !options.output.ends_with(".vcf.gz");
    let vcf = options.output.ends_with(".vcf") || options.output.ends_with(".vcf.gz");
    {
        let mut writer = bcf::Writer::from_path(&options.output, &header, uncompressed, vcf)
            .map_err(|e| format!("Could not open BCF writer: {}", e))?;

        let re = Regex::new(&options.contig_regex).unwrap();

        for (chrom, len) in chrom_lens {
            if !re.is_match(&chrom) {
                continue;
            }
            info!(
                logger,
                "Processing chrom {} ({} bp)",
                &chrom,
                len.separated_string()
            );
            wise_count_count(
                logger,
                &mut tbx_reader,
                &mut bam_reader,
                &mut writer,
                &chrom,
                len,
                &options,
            )?;
        }
    }

    // Build index on the output file.
    shared::build_index(logger, &options.output);

    info!(logger, "Done. Have a nice day!");

    Ok(())
}

/// Compute total number of fragments per bp for normalization.
fn wise_normalize_total_normalized_cov(
    logger: &mut Logger,
    options: &NormalizeOptions,
) -> Result<f64, String> {
    info!(
        logger,
        "Computing normalized fragment count per target bp..."
    );
    let mut reader = bcf::Reader::from_path(&options.input)
        .map_err(|e| format!("Could not open BCF file: {}", e))?;

    let mut vals = Vec::new();

    let mut record = reader.empty_record();
    loop {
        match reader.read(&mut record) {
            Ok(_) => (),
            Err(bcf::ReadError::NoMoreRecord) => break,
            _ => return Err("Could not read record".to_string()),
        }
        let begin = record.pos();
        let end = record
            .info(b"END")
            .integer()
            .map_err(|e| format!("Could not get target END: {}", e))?
            .expect("Could not access END")[0];

        let frag_count = record
            .format(b"FC")
            .float()
            .map_err(|e| format!("Could not access field FC: {}", e))?[0][0]
            as f64;

        let len = (end - begin as i32) as f64;

        vals.push(frag_count / len);
    }

    Ok(vals.as_slice().sum())
}

/// Write out normalized number of fragments per bp.
fn wise_write_normalized_cov(
    logger: &mut Logger,
    options: &NormalizeOptions,
    total_cov: f64,
    writer: &mut bcf::Writer,
) -> Result<(), String> {
    info!(logger, "Computing normalized coverage per target region...");
    let mut reader = bcf::Reader::from_path(&options.input)
        .map_err(|e| format!("Could not open BCF file: {}", e))?;

    let mut record = reader.empty_record();
    loop {
        match reader.read(&mut record) {
            Ok(_) => (),
            Err(bcf::ReadError::NoMoreRecord) => break,
            _ => return Err("Could not read record".to_string()),
        }
        let begin = record.pos();
        let end = record
            .info(b"END")
            .integer()
            .map_err(|e| format!("Could not get target END: {}", e))?
            .expect("Could not access END")[0];

        let frag_count = record
            .format(b"FC")
            .float()
            .map_err(|e| format!("Could not access field FC: {}", e))?[0][0]
            as f64;

        let len = (end - begin as i32) as f64;

        writer.translate(&mut record);

        if frag_count < options.min_fragments as f64 {
            record.push_filter(
                writer
                    .header()
                    .name_to_id(b"LOW_COV")
                    .expect("FILTER 'LOW_COV' unknown"),
            );
        }

        let norm1 = (frag_count / len) as f32;
        let norm2 = (frag_count / len / total_cov) as f32;
        record
            .push_format_float(b"LFC", &[norm1 as f32])
            .map_err(|e| format!("Could not write record: {}", e))?;
        record
            .push_format_float(b"NFC", &[norm2 as f32])
            .map_err(|e| format!("Could not write record: {}", e))?;

        writer
            .write(&record)
            .map_err(|e| format!("Could not write to BCF file: {}", e))?;
    }

    Ok(())
}

/// Main entry point for the "wisexome normalize" command.
pub fn call_wise_normalize(logger: &mut Logger, options: &NormalizeOptions) -> Result<(), String> {
    let options = options.clone();

    info!(logger, "Running cnvetti wisexome normalize");
    info!(logger, "Configuration: {:?}", &options);

    info!(logger, "Opening input file...");
    let mut reader = bcf::Reader::from_path(&options.input)
        .map_err(|e| format!("Could not open BCF file: {}", e))?;
    if options.io_threads > 0 {
        reader
            .set_threads(options.io_threads as usize)
            .map_err(|e| format!("Could not set threads to BCF reader: {}", e))?;
    }

    // Compute total normalized coverage.
    let total_cov = wise_normalize_total_normalized_cov(logger, &options)?;
    info!(logger, "Total normalized coverage is: {}", total_cov);

    {
        info!(logger, "Opening output file...");
        let mut writer = {
            // Construct extended header.
            let mut header = bcf::Header::from_template(reader.header());
            let lines = vec![
                "##FORMAT=<ID=LFC,Number=1,Type=Float,Description=\"Length-normalized  \
                 fragment count\">",
                "##FORMAT=<ID=NFC,Number=1,Type=Float,Description=\"Normalized fragment \
                 count\">",
            ];
            for line in lines {
                header.push_record(line.as_bytes());
            }
            header.push_record(format!("##cnvetti_wiseNormalizeVersion={}", VERSION).as_bytes());
            header.push_record(
                format!(
                    "##cnvetti_wiseNormalizeCommand={}",
                    env::args()
                        .map(|s| shlex::quote(&s).to_string())
                        .collect::<Vec<String>>()
                        .join(" ")
                ).as_bytes(),
            );

            let uncompressed =
                !options.output.ends_with(".bcf") && !options.output.ends_with(".vcf.gz");
            let vcf = options.output.ends_with(".vcf") || options.output.ends_with(".vcf.gz");

            bcf::Writer::from_path(options.output.clone(), &header, uncompressed, vcf)
                .map_err(|e| format!("Could not open output BCF file {}", e))?
        };

        // Compute and write out normalized coverages.
        wise_write_normalized_cov(logger, &options, total_cov, &mut writer)?;
    }

    // Build index on the output file.
    shared::build_index(logger, &options.output);

    info!(logger, "Done. Have a nice day!");

    Ok(())
}

/// Construct BCF Header from input BCF header.
fn call_wise_build_ref_header(header_in: &bcf::header::HeaderView) -> bcf::Header {
    let mut header = bcf::Header::new();
    // TODO(holtgrewe): what about the VCF version?

    // Put overall meta information into the BCF header.
    let now = chrono::Utc::now();
    header.push_record(format!("##fileDate={}", now.format("%Y%m%d").to_string()).as_bytes());

    // Put creating tool version and call into file.
    header.push_record(format!("##cnvetti_wisexomeBuildRefVersion={}", VERSION).as_bytes());
    header.push_record(
        format!(
            "##cnvetti_wisexomeBuildRefCommand={}",
            env::args()
                .map(|s| shlex::quote(&s).to_string())
                .collect::<Vec<String>>()
                .join(" ")
        ).as_bytes(),
    );

    for record in header_in.header_records() {
        match record {
            bcf::header::HeaderRecord::Contig { key: _, values } => {
                header.push_record(
                    format!(
                        "##contig=<ID={},length={}",
                        values.get("ID").unwrap(),
                        values.get("length").unwrap()
                    ).as_bytes(),
                );
            }
            _ => {}
        };
    }

    let lines = vec![
        // Misc fields
        "##ALT=<ID=TARGET,Description=\"Record describes coverage of a target region\">",
        // FILTER tags
        "##FILTER=<ID=FEW_REF,Description=\"Few reference targets available\">",
        "##FILTER=<ID=UNRELIABLE,Description=\"Unreliable; called in too many training samples\">",
        // INFO fields describing the reference.
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Window end\">",
        "##INFO=<ID=REF_TARGETS,Number=.,Type=String,Description=\"The reference targets\">",
        "##INFO=<ID=NUM_CALLS,Number=1,Type=Integer, Description=\"Number of training samples \
         this target was called in; used for UNRELIABLE flag\">",
    ];
    for line in lines {
        header.push_record(line.as_bytes());
    }

    header
}

/// Main entry point for the "wisexome build-ref" command.
pub fn call_wise_build_ref(logger: &mut Logger, options: &BuildRefOptions) -> Result<(), String> {
    let options = options.clone();

    info!(logger, "Running cnvetti wisexome build-ref");
    info!(logger, "Configuration: {:?}", &options);

    // Set number of threads to use by rayon.
    env::set_var("RAYON_NUM_THREADS", format!("{}", options.num_threads));

    info!(logger, "Opening input file...");
    let mut reader = bcf::Reader::from_path(&options.input)
        .map_err(|e| format!("Could not open input BCF file: {}", e))?;
    if options.io_threads > 0 {
        reader
            .set_threads(options.io_threads as usize)
            .map_err(|e| format!("Could not set threads to BCF reader: {}", e))?;
    }

    info!(logger, "Reading normalized fragment counts...");
    let mut target_regions: Vec<(String, u32, u32)> = Vec::new();
    let mut tids: Vec<u32> = Vec::new();
    let mut nfcs: Vec<Vec<f32>> = Vec::new();
    let mut record = reader.empty_record();
    loop {
        match reader.read(&mut record) {
            Ok(_) => (),
            Err(bcf::ReadError::NoMoreRecord) => break,
            _ => return Err("Could not read record".to_string()),
        }
        tids.push(record.rid().unwrap());
        target_regions.push((
            String::from_utf8(reader.header().rid2name(record.rid().unwrap()).to_vec())
                .map_err(|e| format!("Could not parse from UTF-8: {}", e))?,
            record.pos(),
            record.info(b"END").integer().unwrap().unwrap()[0] as u32,
        ));
        let mut tmp = Vec::new();
        for nfc in record
            .format(b"NFC")
            .float()
            .map_err(|e| format!("Could not read FORMAT/NFC: {}", e))?
        {
            tmp.push(nfc[0]);
        }
        nfcs.push(tmp);
    }
    info!(
        logger,
        " => number of targets {}",
        target_regions.len().separated_string()
    );

    let distances = if options.num_threads == 0 {
        info!(logger, "Computing distances (single-threaded)...");
        let mut distances: Vec<Vec<(f32, usize)>> = Vec::new();
        let mut dist: Vec<(f32, usize)> = Vec::with_capacity(nfcs.len());
        for (i, row) in nfcs.iter().enumerate() {
            dist.clear();
            if i % 1000 == 0 {
                debug!(logger, "Done with {} / {}", i, nfcs.len());
            }
            for (j, row2) in nfcs.iter().enumerate() {
                if tids[i] == tids[j] {
                    dist.push((f32::INFINITY, j as usize));
                } else {
                    let mut y = 0_f32;
                    for i in 0..row.len() {
                        let x = row[i] - row2[i];
                        y += x * x;
                    }
                    let y = y.sqrt();
                    dist.push((y, j as usize));
                }
            }
            let size = min(dist.len(), options.max_ref_targets);
            pdqselect::select_by(&mut dist, size, |a, b| a.0.partial_cmp(&b.0).unwrap());
            let mut part = dist[0..size].to_vec();
            part.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            distances.push(part);
        }
        distances
    } else {
        info!(
            logger,
            "Computing distances (with {} threads)...", options.num_threads
        );
        let progress = AtomicUsize::new(0);
        nfcs.par_iter()
            .enumerate()
            .map(|(i, row)| {
                thread_local!{static DIST_RC: RefCell<Vec<(f32, usize)>> = RefCell::new(Vec::new());}

                DIST_RC.with(|dist_rc| {
                    let mut dist = dist_rc.borrow_mut();
                    dist.clear();
                    for (j, row2) in nfcs.iter().enumerate() {
                        if tids[i] == tids[j] {
                            dist.push((f32::INFINITY, j as usize));
                        } else {
                            let mut y = 0_f32;
                            for i in 0..row.len() {
                                let x = row[i] - row2[i];
                                y += x * x;
                            }
                            let y = y.sqrt();
                            dist.push((y, j as usize));
                        }
                    }

                    let size = min(dist.len(), options.max_ref_targets);
                    pdqselect::select_by(&mut dist, size, |a, b| a.0.partial_cmp(&b.0).unwrap());

                    let mut part = dist[0..size].to_vec();
                    part.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

                    loop {
                        let val = progress.load(Ordering::Relaxed);
                        if progress
                            .compare_exchange(val, val + 1, Ordering::Acquire, Ordering::Relaxed)
                            .is_ok()
                        {
                            if val % 10000 == 0 {
                                info!(
                                    logger,
                                    "Processed {} targets so far",
                                    val.separated_string()
                                );
                            }
                            break;
                        }
                    }

                    part
                })
            })
            .collect::<Vec<Vec<(f32, usize)>>>()
    };
    info!(logger, " => done");

    info!(logger, "Pruning distances...");
    // Compute distance error statistics.
    let error_stats: Vec<f64> = distances
        .iter()
        .map(|xs| xs[0].0 as f64)
        .collect::<Vec<f64>>();
    let mean = error_stats.as_slice().mean();
    let std_dev = error_stats.as_slice().std_dev();
    let threshold = (mean + options.filter_z_score * std_dev) as f32;

    // Actually collect the reliable regions and flag for reliability.
    let mut ref_targets: Vec<Vec<usize>> = Vec::new();

    for row in distances.iter() {
        ref_targets.push(
            row.iter()
                .filter(|(dist, _)| *dist <= threshold)
                .map(|(_, idx)| *idx)
                .collect::<Vec<usize>>(),
        );
    }
    info!(logger, " => done");

    // For each probe, count samples that have a call in it based on z-score.
    info!(logger, "Count per-probe calls for UNRELIABLE flag...");
    let num_calls: Vec<u32> = {
        let mut num_calls = vec![0; target_regions.len()];

        for (target_idx, matched_targets) in ref_targets.iter().enumerate() {
            for sample_idx in 0..(reader.header().sample_count() as usize) {
                // Get empirical reference distribution for target `target_idx` for sample `sample_idx`.
                let ref_dist = matched_targets
                    .iter()
                    .map(|idx| nfcs[*idx][sample_idx] as f64)
                    .collect::<Vec<f64>>();
                if ref_dist.len() < options.min_ref_targets {
                    continue;
                }
                let ref_dist = ref_dist.as_slice();
                // Compute z-score
                let x = nfcs[target_idx][sample_idx] as f64;
                let z_score = (x - ref_dist.mean()).abs() / ref_dist.std_dev();
                // Compute relative value.
                let rel = x / ref_dist.mean();
                // Increment per-probe counter for target.
                if z_score > options.filter_z_score && (rel - 1.0).abs() > options.filter_rel {
                    num_calls[target_idx] += 1;
                }
            }
        }
        info!(logger, " => done");

        num_calls
    };

    // Actualy collect the reliable regions and flag for reliability.
    info!(logger, "Writing output BCF file...");
    let header = call_wise_build_ref_header(reader.header());
    let uncompressed = !options.output.ends_with(".bcf") && !options.output.ends_with(".vcf.gz");
    let vcf = options.output.ends_with(".vcf") || options.output.ends_with(".vcf.gz");
    {
        let mut writer = bcf::Writer::from_path(&options.output, &header, uncompressed, vcf)
            .map_err(|e| format!("Could not open BCF writer: {}", e))?;

        for (i, tgt) in target_regions.iter().enumerate() {
            let mut record = writer.empty_record();

            record.set_rid(&Some(tids[i]));
            record.set_pos(tgt.1 as i32);
            record
                .set_id(format!("{}:{}-{}", tgt.0, tgt.1 + 1, tgt.2).as_bytes())
                .map_err(|e| format!("Could not write ID: {}", e))?;

            let header = record.header().clone();
            if ref_targets[i].len() < options.min_ref_targets {
                record.push_filter(
                    header
                        .name_to_id(b"FEW_REF")
                        .expect("FILTER 'FEW_REF' unknown"),
                );
            }
            if num_calls[i] > options.max_samples_reliable {
                record.push_filter(
                    header
                        .name_to_id(b"UNRELIABLE")
                        .expect("Filter 'UNRELIABLE' unknown"),
                );
            }

            let alleles_v = vec![Vec::from("N"), Vec::from("<TARGET>")];
            let alleles = alleles_v
                .iter()
                .map(|x| x.as_slice())
                .collect::<Vec<&[u8]>>();
            record
                .set_alleles(&alleles)
                .map_err(|e| format!("Could not update alleles: {}", e))?;
            record
                .push_info_integer(b"END", &[tgt.2 as i32])
                .map_err(|e| format!("Could not write INFO/END: {}", e))?;
            record
                .push_info_integer(b"NUM_CALLS", &[num_calls[i] as i32])
                .map_err(|e| format!("Could not write INFO/NUM_CALLS: {}", e))?;

            let ref_targets = ref_targets[i]
                .iter()
                .map(|idx| {
                    format!(
                        "{}:{}-{}",
                        target_regions[*idx].0,
                        target_regions[*idx].1 + 1,
                        target_regions[*idx].2
                    )
                })
                .collect::<Vec<String>>();
            let ref_targets_b = ref_targets
                .iter()
                .map(|s| s.as_bytes())
                .collect::<Vec<&[u8]>>();
            record
                .push_info_string(b"REF_TARGETS", ref_targets_b.as_slice())
                .map_err(|e| format!("Could not write INFO/REF_TARGETS: {}", e))?;

            writer
                .write(&record)
                .map_err(|e| format!("Problem writing to output: {}", e))?;
        }
    }
    // Build index on the output file.
    shared::build_index(logger, &options.output);
    info!(logger, " => done");

    Ok(())
}

// Calling state.
#[derive(Debug, Display, Clone, Copy, EnumString, PartialEq)]
enum CallState {
    // Called as reference.
    Ref,
    // Called as duplicate.
    Dup,
    // Called as deletion.
    Del,
}

// Information used in the calling.
//
// Created only for non-filtered probes.
#[derive(Debug, Clone)]
struct ProbeInfo {
    // Reference chromosome.
    chrom: String,
    // Start position.
    start: u32,
    // End position.
    end: u32,
    // Name of the probe.
    name: String,
    // Normalized coverage at the probe.
    norm_cov: f64,
    // Reference mean.
    ref_mean: f64,
    // Reference std deviation.
    ref_std_dev: f64,
    // Coverage relative to ref_norm_covs; called "effect size" by Straver et al.
    rel: f64,
    // Window-wide "effect size", starts out as window coverage.
    window_rel: f64,
    // Z score of norm_cov with respect to ref_norm_covs.
    z_score: f64,
    // Call state.
    call_state: CallState,
}

impl ProbeInfo {
    /// Update call state based on probe effect size and z-score.
    fn update_call_state_probe(&mut self, thresh_rel: f64, thresh_z: f64) {
        self.call_state = if (self.rel - 1.0).abs() < thresh_rel {
            CallState::Ref
        } else if self.z_score <= -thresh_z {
            CallState::Del
        } else if self.z_score >= thresh_z {
            CallState::Dup
        } else {
            CallState::Ref
        };
    }

    /// Update call state based on window effect size, can only switch to reference if threshold
    /// not reached.
    fn update_call_state_window(&mut self, thresh_rel: f64) {
        if (self.window_rel - 1.0).abs() < thresh_rel {
            self.call_state = CallState::Ref;
        }
    }
}

/// Information of a CNV.
struct CnvInfo {
    /// Chromosome of the CNV.
    chrom: String,
    /// Start of the CNV.
    start: u32,
    /// End of the CNV.
    end: u32,
    /// Relative count over the window.
    rel_count: f64,
    /// Number of targets.
    num_targets: u32,
    /// Del/dup
    call_state: CallState,
}

/// Construct BCF Header from input BCF header.
fn call_wise_build_call_header(header_in: &bcf::header::HeaderView) -> bcf::Header {
    let mut header = bcf::Header::new();
    // TODO(holtgrewe): what about the VCF version?

    // Put overall meta information into the BCF header.
    let now = chrono::Utc::now();
    header.push_record(format!("##fileDate={}", now.format("%Y%m%d").to_string()).as_bytes());

    // Put creating tool version and call into file.
    header.push_record(format!("##cnvetti_wisexomeCallVersion={}", VERSION).as_bytes());
    header.push_record(
        format!(
            "##cnvetti_wisexomeCallCommand={}",
            env::args()
                .map(|s| shlex::quote(&s).to_string())
                .collect::<Vec<String>>()
                .join(" ")
        ).as_bytes(),
    );

    for record in header_in.header_records() {
        match record {
            bcf::header::HeaderRecord::Contig { key: _, values } => {
                header.push_record(
                    format!(
                        "##contig=<ID={},length={}",
                        values.get("ID").unwrap(),
                        values.get("length").unwrap()
                    ).as_bytes(),
                );
            }
            _ => {}
        };
    }

    let lines = vec![
        // Misc fields
        "##ALT=<ID=DUP,Description=\"Record describing duplication event\">",
        "##ALT=<ID=DEL,Description=\"Record describing deletion event\">",
        // FILTER tags
        // INFO fields describing the reference.
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"SV end\">",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">",
        "##INFO=<ID=REGION,Number=1,Type=Integer,Description=\"Affected genomics region\">",
        "##INFO=<ID=NUM_TARGETS,Number=1,Type=Integer,Description=\"The number of targets in \
         CNV (similar to aCGH probes)\">",
        // FORMAT fields
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=RC,Number=1,Type=Float,Description=\"Relative coverage\">",
    ];
    for line in lines {
        header.push_record(line.as_bytes());
    }

    for sample in &header_in.samples() {
        header.push_sample(sample);
    }

    header
}

/// Main entry point for the "wisexome call" command.
pub fn call_wise_call(logger: &mut Logger, options: &CallOptions) -> Result<(), String> {
    info!(logger, "Running cnvetti wisexome call");
    info!(logger, "Configuration: {:?}", &options);

    // Set number of threads to use by rayon.
    env::set_var("RAYON_NUM_THREADS", format!("{}", options.num_threads));

    info!(logger, "Loading normalized counts...");
    let reader_header;
    let normalized_counts = {
        let mut reader = bcf::Reader::from_path(&options.input)
            .map_err(|e| format!("Could not open input BCF file: {}", e))?;
        reader_header = reader.header().clone();
        let mut normalized_counts = HashMap::new();
        let mut record = reader.empty_record();
        loop {
            match reader.read(&mut record) {
                Ok(_) => {}
                Err(bcf::ReadError::NoMoreRecord) => break,
                _ => return Err("Could not read record".to_string()),
            }

            let name = String::from_utf8(record.id())
                .map_err(|e| format!("Problem converting from UTF-8: {}", e))?;

            normalized_counts.insert(
                name,
                record
                    .format(b"NFC")
                    .float()
                    .map_err(|e| format!("Could not access field FC: {}", e))?[0][0]
                    as f64,
            );
        }
        normalized_counts
    };
    info!(logger, " => done");

    info!(logger, "Extracting target-wise information...");
    let mut probe_infos = {
        let mut reader = bcf::Reader::from_path(&options.input_ref)
            .map_err(|e| format!("Could not open reference BCF file: {}", e))?;
        let mut record = reader.empty_record();
        let mut probe_infos: Vec<ProbeInfo> = Vec::new();
        loop {
            match reader.read(&mut record) {
                Ok(_) => {
                    // Skip any filtered record from reference.
                    if record
                        .filters()
                        .map(|id| {
                            String::from_utf8(record.header().id_to_name(id))
                                .expect(&format!("Could not transflate from UTF-8"))
                        })
                        .any(|s| s != "PASS")
                    {
                        continue;
                    }
                }
                Err(bcf::ReadError::NoMoreRecord) => break,
                _ => return Err("Could not read record".to_string()),
            }

            let name = String::from_utf8(record.id())
                .map_err(|e| format!("Cannot decode from UTF-8: {}", e))?;

            // Get normalized count at this probe and the empirical distribution from the reference.
            let norm_cov = *normalized_counts
                .get(&name)
                .expect(&format!("Could not get norm count for probe: {}", name));
            let empty: Vec<&[u8]> = Vec::new();
            let ref_norm_covs = record
                .info(b"REF_TARGETS")
                .string()
                .map_err(|e| format!("Problem accessing REF_TARGETS: {}", e))?
                .unwrap_or(empty)
                .iter()
                .map(|v| {
                    let s: String =
                        String::from_utf8(v.to_vec()).expect("Decoding from UTF-8 failed");
                    let xs: Vec<String> =
                        s.split(",").map(|s| s.to_string()).collect::<Vec<String>>();
                    xs
                })
                .flatten()
                .map(|s| {
                    *normalized_counts
                        .get(&s)
                        .expect(format!("Could not find: {}", s).as_ref())
                })
                .collect::<Vec<f64>>();

            // Compute initial set of statistics.
            if !ref_norm_covs.is_empty() {
                let ref_mean = ref_norm_covs.as_slice().mean();
                let ref_std_dev = ref_norm_covs.as_slice().std_dev();
                let (z_score, rel) = if ref_std_dev == 0.0 {
                    (0.0, 0.0)
                } else {
                    ((norm_cov - ref_mean) / ref_std_dev, norm_cov / ref_mean)
                };

                let name = String::from_utf8(record.id())
                    .map_err(|e| format!("Cannot decode from UTF-8: {}", e))?;

                probe_infos.push(ProbeInfo {
                    chrom: String::from_utf8(
                        reader.header().rid2name(record.rid().unwrap()).to_vec(),
                    ).unwrap(),
                    start: record.pos(),
                    end: record
                        .info(b"END")
                        .integer()
                        .map_err(|e| format!("Could not get target END: {}", e))?
                        .expect("Could not access END")[0] as u32,
                    name: name,
                    norm_cov,
                    ref_mean,
                    ref_std_dev,
                    rel,
                    window_rel: rel,
                    z_score,
                    // Note that the call state will be updated by window later below.
                    call_state: CallState::Ref,
                });
            }
        }
        probe_infos
    };
    info!(logger, " => done");

    info!(logger, "Performing calling and segmentation...");
    info!(logger, "-> windowing & score updates");
    // Perform windowing and z-score updates.
    //
    // First compute z scores of varying size windows.
    let n = probe_infos.len();
    let m = options.max_window_size / 2 + 1;
    let mut window_z_scores: Vec<Vec<f64>> = vec![vec![0.0; m]; n];
    for i in 0..m {
        debug!(logger, "..> window size = {}", 2 * i + 1);
        let delta = 2 * i;
        for j in 0..n {
            let left = if j < delta { 0 } else { j - delta };
            let right = if j + delta + 1 >= n { n } else { j + delta + 1 };
            window_z_scores[j][i] = probe_infos[left..right]
                .iter()
                .map(|i| i.z_score)
                .sum::<f64>() / ((right - left) as f64).sqrt();
        }
    }
    // Update z score of probes.
    for (i, zs) in window_z_scores.iter().enumerate() {
        probe_infos[i].z_score = *zs.iter().max_by_key(|x| OrderedFloat(x.abs())).unwrap();
        probe_infos[i].update_call_state_probe(options.thresh_rel, options.thresh_z);
    }

    // Compute window effect size and update call state.
    //
    // First, collect windows, compute window and update effect sizes.
    info!(logger, "-> effect sizes");
    let mut curr_window: Option<(usize, usize, CallState)> = None;
    for i in 0..probe_infos.len() {
        let info = probe_infos[i].clone();
        curr_window = match (&info.call_state, &curr_window) {
            (CallState::Ref, None) => None,
            (_, None) => Some((i, i + 1, info.call_state)),
            (CallState::Ref, Some((start, end, _))) => {
                let window_rel = probe_infos[*start..*end]
                    .iter()
                    .map(|info| info.rel)
                    .collect::<Vec<f64>>()
                    .as_slice()
                    .median();
                for i in *start..*end {
                    probe_infos[i].window_rel = window_rel;
                }
                None
            }
            (_, Some((start, end, call_state))) => {
                if info.call_state == *call_state {
                    Some((*start, i + 1, *call_state))
                } else {
                    // TODO: extract these two statements into one function?
                    let window_rel = probe_infos[*start..*end]
                        .iter()
                        .map(|info| info.rel)
                        .collect::<Vec<f64>>()
                        .as_slice()
                        .median();
                    for i in *start..*end {
                        probe_infos[i].window_rel = window_rel;
                    }
                    Some((i, i + 1, info.call_state))
                }
            }
        }
    }
    if let Some((start, end, _)) = curr_window {
        let window_rel = probe_infos[start..end]
            .iter()
            .map(|info| info.rel)
            .collect::<Vec<f64>>()
            .as_slice()
            .median();
        for i in start..end {
            probe_infos[i].window_rel = window_rel;
        }
    }
    // Then, update call state based on window's effect size.
    for info in &mut probe_infos {
        info.update_call_state_window(options.thresh_rel);
    }

    // Perform fine-tuning on all windows.
    //
    // First, compute window ranges.
    info!(logger, "-> fine-tuning");
    let mut windows: Vec<(usize, usize, CallState)> = Vec::new();
    let mut curr_window: Option<(usize, usize, CallState)> = None;
    for i in 0..probe_infos.len() {
        let info = probe_infos[i].clone();
        curr_window = match (&info.call_state, &curr_window) {
            (CallState::Ref, None) => None,
            (_, None) => Some((i, i + 1, info.call_state)),
            (CallState::Ref, Some(triple)) => {
                windows.push(*triple);
                None
            }
            (_, Some((start, end, call_state))) => {
                if info.call_state == *call_state {
                    windows.push((*start, *end, *call_state));
                    Some((*start, i + 1, *call_state))
                } else {
                    Some((i, i + 1, info.call_state))
                }
            }
        }
    }
    if let Some(triple) = curr_window {
        windows.push(triple);
    }
    // Then, perform the actual fine-tuning.
    const DELTA: usize = 8;
    for (start, end, call_state) in windows {
        let start = if start > DELTA { start - DELTA } else { 0 };
        let end = if end + DELTA > probe_infos.len() {
            probe_infos.len()
        } else {
            end + DELTA
        };

        let mut best: Option<(usize, usize, CallState, f64)> = None;
        for i in start..end {
            for j in i..end {
                if i < j {
                    // For determining the window, we are using the *mean* and not the *median*.
                    let window_rel = probe_infos[i..j]
                        .iter()
                        .map(|info| info.rel)
                        .collect::<Vec<f64>>()
                        .as_slice()
                        .mean();

                    if (1.0 - window_rel).abs() < options.thresh_rel {
                        continue; // below threshold, ignore
                    }

                    best = match best {
                        None => Some((start, end, call_state, window_rel)),
                        Some((s, e, c, r)) => {
                            if window_rel > r {
                                Some((start, end, call_state, window_rel))
                            } else {
                                Some((s, e, c, r))
                            }
                        }
                    }
                }
            }
        }

        if let Some((start, end, call_state, _)) = best {
            // Compute true window_rel with median.
            let window_rel = probe_infos[start..end]
                .iter()
                .map(|info| info.rel)
                .collect::<Vec<f64>>()
                .as_slice()
                .median();
            for info in &mut probe_infos[start..end] {
                info.window_rel = window_rel;
                info.call_state = call_state;
                info.update_call_state_window(options.thresh_rel);
            }
        }
    }
    info!(logger, " => done");

    // Compute empirical p-values of for each window.
    info!(logger, "(Not computing empirical p values yet...)");
    // TODO

    info!(logger, "Collecting CNVs...");
    let mut cnvs: Vec<CnvInfo> = Vec::new();
    let mut current: Option<CnvInfo> = None;
    for ref probe_info in &probe_infos {
        if probe_info.call_state == CallState::Ref {
            if let Some(inner_current) = current {
                cnvs.push(inner_current);
                current = None;
            }
            continue;
        }
        current = match current {
            None => Some(CnvInfo {
                chrom: probe_info.chrom.clone(),
                start: probe_info.start,
                end: probe_info.end,
                rel_count: probe_info.window_rel,
                call_state: probe_info.call_state,
                num_targets: 0,
            }),
            Some(current) => {
                if probe_info.chrom != current.chrom || probe_info.call_state != current.call_state
                {
                    cnvs.push(current);
                    Some(CnvInfo {
                        chrom: probe_info.chrom.clone(),
                        start: probe_info.start,
                        end: probe_info.end,
                        rel_count: probe_info.window_rel,
                        call_state: probe_info.call_state,
                        num_targets: 0,
                    })
                } else {
                    Some(CnvInfo {
                        end: probe_info.end,
                        num_targets: current.num_targets + 1,
                        ..current
                    })
                }
            }
        }
    }
    if let Some(cnv_info) = current {
        cnvs.push(cnv_info);
    }

    info!(logger, "Writing output call BCF file...");
    let uncompressed = !options.output.ends_with(".bcf") && !options.output.ends_with(".vcf.gz");
    let vcf = options.output.ends_with(".vcf") || options.output.ends_with(".vcf.gz");
    {
        let mut writer = bcf::Writer::from_path(
            &options.output,
            &call_wise_build_call_header(&reader_header),
            uncompressed,
            vcf,
        ).map_err(|e| format!("Could not open BCF writer: {}", e))?;

        for ref cnv_info in cnvs {
            let mut record = writer.empty_record();

            let rid = writer
                .header()
                .name2rid(cnv_info.chrom.as_bytes())
                .map_err(|e| format!("Could not find contig {}: {}", cnv_info.chrom, e))?;
            record.set_rid(&Some(rid));

            record.set_pos(cnv_info.start as i32);

            let alleles_v = vec![
                Vec::from("N"),
                if cnv_info.call_state == CallState::Del {
                    Vec::from("<DEL>")
                } else {
                    Vec::from("<DUP>")
                },
            ];
            let alleles = alleles_v
                .iter()
                .map(|x| x.as_slice())
                .collect::<Vec<&[u8]>>();
            record
                .set_alleles(&alleles)
                .map_err(|e| format!("Could not update alleles: {}", e))?;

            record
                .push_info_integer(b"END", &[cnv_info.end as i32])
                .map_err(|e| format!("Could not write INFO/END: {}", e))?;

            let svlen = if cnv_info.call_state == CallState::Del {
                (cnv_info.start as i32) - (cnv_info.end as i32)
            } else {
                (cnv_info.end as i32) - (cnv_info.start as i32)
            };
            record
                .push_info_integer(b"SVLEN", &[svlen])
                .map_err(|e| format!("Could not write INFO/END: {}", e))?;
            record
                .push_info_integer(b"NUM_TARGETS", &[cnv_info.num_targets as i32])
                .map_err(|e| format!("Could not write INFO/NUM_TARGETS: {}", e))?;

            let region = format!("{}:{}-{}", cnv_info.chrom, cnv_info.start + 1, cnv_info.end);
            record
                .push_info_string(b"REGION", &[region.as_bytes()])
                .map_err(|e| format!("Could not write INFO/NUM_TARGETS: {}", e))?;

            record
                .push_info_string(
                    b"SVTYPE",
                    &[if cnv_info.call_state == CallState::Del {
                        b"DEL"
                    } else {
                        b"DUP"
                    }],
                )
                .map_err(|e| format!("Could not write INFO/END: {}", e))?;

            record
                .push_format_integer(b"GT", &[1])
                .map_err(|e| format!("Could not write FORMAT/GT: {}", e))?;
            record
                .push_format_float(b"RC", &[cnv_info.rel_count as f32])
                .map_err(|e| format!("Could not write FORMAT/RC: {}", e))?;

            writer
                .write(&record)
                .map_err(|e| format!("Problem writing record: {}", e))?;
        }
    }

    // Build index on the output file.
    shared::build_index(logger, &options.output);

    if !options.output_probes.is_empty() {
        info!(logger, "Writing per-probe output BCF file...");
        let mut out_infos: HashMap<String, (f64, f64, f64, CallState)> = HashMap::new();
        for info in &probe_infos {
            out_infos.insert(
                info.name.clone(),
                (info.z_score, info.rel, info.window_rel, info.call_state),
            );
        }

        let mut reader = bcf::Reader::from_path(&options.input)
            .map_err(|e| format!("Could not open BCF file: {}", e))?;
        let uncompressed =
            !options.output_probes.ends_with(".bcf") && !options.output_probes.ends_with(".vcf.gz");
        let vcf =
            options.output_probes.ends_with(".vcf") || options.output_probes.ends_with(".vcf.gz");
        {
            let mut writer = {
                // Construct extended header.
                let mut header = bcf::Header::from_template(reader.header());
                let lines = vec![
                    "##FILTER=<ID=FEW_REF,Description=\"Few reference targets available\">",
                    "##FORMAT=<ID=TZ,Number=1,Type=Float,Description=\"Target z-score\">",
                    "##FORMAT=<ID=TE,Number=1,Type=Float,Description=\"Target effect size\">",
                    "##FORMAT=<ID=WE,Number=1,Type=Float,Description=\"Window-wise effect size\">",
                    "##FORMAT=<ID=CALL,Number=1,Type=String,Description=\"Call state\">",
                ];
                for line in lines {
                    header.push_record(line.as_bytes());
                }
                header.push_record(format!("##cnvetti_wiseCallVersion={}", VERSION).as_bytes());
                header.push_record(
                    format!(
                        "##cnvetti_wiseCallCommand={}",
                        env::args()
                            .map(|s| shlex::quote(&s).to_string())
                            .collect::<Vec<String>>()
                            .join(" ")
                    ).as_bytes(),
                );

                bcf::Writer::from_path(&options.output_probes, &header, uncompressed, vcf)
                    .map_err(|e| format!("Could not open BCF writer: {}", e))?
            };

            let mut record = reader.empty_record();
            loop {
                match reader.read(&mut record) {
                    Ok(_) => (),
                    Err(bcf::ReadError::NoMoreRecord) => break,
                    _ => return Err("Could not read record".to_string()),
                }

                writer.translate(&mut record);

                let name = String::from_utf8(record.id())
                    .map_err(|e| format!("Cannot decode from UTF-8: {}", e))?;

                if let Some((z_score, rel, window_rel, call_state)) = out_infos.get(&name) {
                    record
                        .push_format_float(b"TZ", &[*z_score as f32])
                        .map_err(|e| format!("Problem setting FORMAT/TZ: {}", e))?;
                    record
                        .push_format_float(b"TE", &[*rel as f32])
                        .map_err(|e| format!("Problem setting FORMAT/TE: {}", e))?;
                    record
                        .push_format_float(b"WE", &[*window_rel as f32])
                        .map_err(|e| format!("Problem setting FORMAT/WE: {}", e))?;
                    record
                        .push_format_string(b"CALL", &[(*call_state).to_string().as_bytes()])
                        .map_err(|e| format!("Problem setting FORMAT/CALL: {}", e))?;
                } else {
                    record.push_filter(
                        writer
                            .header()
                            .name_to_id(b"FEW_REF")
                            .expect("FILTER 'FEW_REF' unknown in header"),
                    );
                }

                writer
                    .write(&record)
                    .map_err(|e| format!("Problem writing record: {}", e))?;
            }
        }

        // Build index on the per-probe output file.
        shared::build_index(logger, &options.output_probes);
    }

    Ok(())
}
