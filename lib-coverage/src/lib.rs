/// Main library module for building initial coverage BCF files from BAM files.
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::Write;

extern crate bio;
use bio::data_structures::interval_tree::IntervalTree;

extern crate chrono;

#[macro_use]
extern crate error_chain;

extern crate clap;

extern crate separator;
use separator::Separatable;

#[macro_use]
extern crate slog;
use slog::Logger;

extern crate strum;
#[macro_use]
extern crate strum_macros;

extern crate rust_htslib;
use rust_htslib::bam::{self, Read as BamRead};
use rust_htslib::bcf::{self, Read as BcfRead};
use rust_htslib::tbx::{self, Read as TbxRead};

extern crate shlex;

extern crate statrs;
use statrs::distribution::{Discrete, Poisson, Univariate};

extern crate lib_shared;
use lib_shared::bam_utils;
use lib_shared::bcf_utils;

mod agg;
use agg::*;
mod options;
pub use options::*;
mod piles;
mod reference;
use reference::ReferenceStats;

use lib_shared::regions::GenomeRegions;

/// This crate's error-related code, generated by `error-chain`.
mod errors {
    // Create the Error, ErrorKind, ResultExt, and Result types
    error_chain!{}
}

pub use errors::*;

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
    header.push_record(format!("##cnvetti_cmdCoverageVersion={}", "0.1.0").as_bytes());
    header.push_record(
        format!(
            "##cnvetti_cmdCoverageCommand={}",
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
    // TODO: later decide about commented-out lines
    let lines = vec![
        // Define ALT column <WINDOW>/<TARGET>
        "##ALT=<ID=WINDOW,Description=\"Record describes a window for read or coverage \
         counting\">",
        "##ALT=<ID=TARGET,Description=\"Record describes a target for fragment counting\">",
        // INFO fields describing the window
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Window end\">",
        "##INFO=<ID=GC,Number=1,Type=Float,Description=\"Reference GC fraction, if reference \
         FASTA file was given\">",
        "##INFO=<ID=MEAN_MAPQ,Number=1,Type=Float,Description=\"Mean MAPQ value across samples \
         for approximating mapability\">",
        "##INFO=<ID=FEW_GCWINDOWS,Number=0,Type=Flag,Description=\"There are too few windows \
         with this GC content for robust statistics, ignoring (only in window-based GC \
         correction)\">",
        "##INFO=<ID=REF_MEAN,Number=1,Type=Flag,Description=\"Reference mean coverage\">",
        "##INFO=<ID=REF_STD_DEV,Number=1,Type=Flag,Description=\"Reference mean standard \
         deviation\">",
        "##INFO=<ID=REF_IQR,Number=1,Type=Flag,Description=\"Inter-quartile range of reference\">",
        "##INFO=<ID=REF_5SUMMARY,Number=5,Type=Flag,Description=\"Five-number summary of \
         reference coverage\">",
        // FILTER fields
        "##FILTER=<ID=PILE,Description=\"Window masked in sample because piles make up a
         too large fraction of it.\">",
        "##FILTER=<ID=LOWCOV,Description=\"Window has too low raw coverage (in number of fragment;
         only used for targeted sequencing CNVs)\">",
        "##FILTER=<ID=NO_REF,Description=\"No reference targets in WIS method for this \
         target\">",
        "##FILTER=<ID=SKIPPED_SEG,Description=\"Skipped region on segmentation.\">",
        // "##INFO=<ID=BLACKLIST,Number=0,Type=Flag,Description=\"Interval overlaps with \
        //  blacklist site\">",
        "##INFO=<ID=GAP,Number=0,Type=Flag,Description=\"Window overlaps with N in \
         reference (gap)\">",
        // Generic FORMAT fields
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=MQ,Number=1,Type=Float,Description=\"Mean read MAPQ from region\">",
        // Note that the meaning of the term "coverage" differs considerably between
        // the different approaches. When processing targeted data, it is the length
        // normalized fraction of reads in this sample in the current target. In
        // contrast, for WGS or off-target reads, it is a normalized coverag where
        // distribution of all values on the autosomes is expected to center around
        // the center of diploid positions without any alterations.
        "##FORMAT=<ID=RCV,Number=1,Type=Float,Description=\"Raw coverage value\">",
        "##FORMAT=<ID=RCVSD,Number=1,Type=Float,Description=\"Raw coverage standard deviation\">",
        "##FORMAT=<ID=LCV,Number=1,Type=Float,Description=\"Length-normalized coverage \
         value\">",
        "##FORMAT=<ID=LCVSD,Number=1,Type=Float,Description=\"Length-normalized coverage \
         standard deviation\">",
        "##FORMAT=<ID=NCV,Number=1,Type=Float,Description=\"Further normalized coverage \
         value\">",
        "##FORMAT=<ID=CV,Number=1,Type=Float,Description=\"Finalized coverage value as input \
         to discovery step (or similar steps such as segmentation)\">",
        "##FORMAT=<ID=CVSD,Number=1,Type=Float,Description=\"Finalized coverage standard deviation \
         (or similar steps such as segmentation); only when considering base-wise coverage\">",
        "##FORMAT=<ID=CV2,Number=1,Type=Float,Description=\"Log2-scaled, finalized coverage \
         value as input to discovery step (or similar steps such as segmentation)\">",
        "##FORMAT=<ID=CVZ,Number=1,Type=Float,Description=\"Finalized coverage value as Z-score
         of reference targets distribution.\">",
        "##FORMAT=<ID=FRM,Number=1,Type=Float,Description=\"Fraction of window that was masked
         in the sample.\">",
        "##FORMAT=<ID=SG,Number=1,Type=Float,Description=\"Coverage value in current segment.\">",
        "##FORMAT=<ID=SG2,Number=1,Type=Float,Description=\"Log2-scale coverage value in current \
         segment.\">",
        "##FORMAT=<ID=SGP,Number=1,Type=Float,Description=\"P-value of original segment.\">",
        "##FORMAT=<ID=SGZ,Number=1,Type=Float,Description=\"Coverage Z-score value in current \
         segment.\">",
        "##FORMAT=<ID=SGS,Number=1,Type=Integer,Description=\"Copy number state of current
         segment from HMM (0: del, 1: neutral, 2: dup).\">",
        // TODO: we might need to add a "TNCV" here for "total-normalized coverage" as generate
        //       when normalizing counts for the WIS approach.
        // FORMAT fields used when calling by coverage (deep WGS only)
        "##FORMAT=<ID=CVSD,Number=1,Type=Float,Description=\"Per-window coverage SD\">",
        // FORMAT fields used for read counts (off-target WES, shallow WGS)
        "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Filter values for the genotype\">",
    ];
    for line in lines {
        header.push_record(line.as_bytes());
    }

    header
}

/// Build bcf::Writer with appropriate header.
fn build_bcf_writer(
    path: &String,
    samples: &Vec<String>,
    ref_contigs: &GenomeRegions,
) -> Result<bcf::Writer> {
    let uncompressed = !path.ends_with(".bcf") && !path.ends_with(".vcf.gz");
    let vcf = path.ends_with(".vcf") || path.ends_with(".vcf.gz");

    let header = build_header(samples, ref_contigs);
    bcf::Writer::from_path(&path, &header, uncompressed, vcf)
        .chain_err(|| "Could not open BCF file for writing")
}

/// Load genome regions from tabix-indexed BED file.
///
/// This can be used for loading blacklists or target regions.
fn load_bed_regions_from_tabix(
    logger: &mut Logger,
    tbx_reader: &mut tbx::Reader,
    chrom: &str,
    len: usize,
) -> Result<GenomeRegions> {
    let tid = match tbx_reader.tid(&chrom) {
        Ok(tid) => tid,
        Err(e) => {
            debug!(logger, "Could not map contig to ID: {}: {}", &chrom, e);
            return Ok(GenomeRegions::new());
        }
    };
    tbx_reader
        .fetch(tid, 0, len as u32)
        .chain_err(|| format!("Could not fetch region {}:1-{}", &chrom, len))?;
    let mut result = GenomeRegions::new();

    for buf in tbx_reader.records() {
        let s = String::from_utf8(buf.unwrap()).unwrap();
        let arr: Vec<&str> = s.split('\t').collect();
        if arr.len() < 3 {
            bail!("BED file line had too few columns! {} (<3)");
        }
        let begin = arr[1].parse::<usize>().unwrap();
        let end = arr[2].parse::<usize>().unwrap();
        result.regions.push((chrom.to_string(), begin, end));
    }

    Ok(result)
}

/// Load regions from WIS or pool-based model BCF file.
fn load_regions_from_model_bcf(
    logger: &mut Logger,
    filename: &String,
    chrom: &str,
    len: usize,
) -> Result<GenomeRegions> {
    let mut reader =
        bcf::IndexedReader::from_path(&filename).chain_err(|| "Could not open BCF file")?;

    let tid = match reader.header().name2rid(chrom.as_bytes()) {
        Ok(tid) => tid,
        Err(e) => {
            debug!(logger, "Could not map contig to ID: {}: {}", &chrom, e);
            return Ok(GenomeRegions::new());
        }
    };

    reader
        .fetch(tid, 0, len as u32)
        .chain_err(|| format!("Could not fetch region {}:1-{}", &chrom, len))?;

    let mut result = GenomeRegions::new();
    let mut record = reader.empty_record();

    loop {
        match reader.read(&mut record) {
            Ok(_) => (),
            Err(bcf::ReadError::NoMoreRecord) => break,
            _ => bail!("Could not read BCF record"),
        }

        let end = record.info(b"END").integer().unwrap().unwrap()[0] as usize;
        result
            .regions
            .push((chrom.to_string(), record.pos() as usize, end));
    }

    Ok(result)
}

/// Compute FDR threshold given absolute frequency histogram.
///
/// Use same logic as in `peakCutoff` of `Bioconductor/chipseq`.
///
/// ## Parameters
///
/// - `y` -- absolute frequencies
fn compute_threshold(logger: &mut Logger, y: &Vec<usize>, fdr_cutoff: f64) -> Result<u32> {
    let k_max: u32 = 20;
    let k_lim: u32 = k_max + 1;
    let k: Vec<u32> = (2..k_lim).collect::<Vec<u32>>();

    if y.len() < 3 {
        bail!("No coverage above 1?");
    } else if y[0] != 0 {
        bail!("Logic error, y[0] must be == 0");
    }

    // Get local copy of sufficiently sized `y`.
    let mut y = y.clone();
    let new_len = std::cmp::max(k_lim as usize, y.len());
    y.resize(new_len, 0);

    // Get empirical estimate of `\lambda`.
    let y_1 = y[1] as f64;
    let y_2 = y[2] as f64;
    let lambda = 2.0 * y_2 / y_1;
    debug!(logger, " -> lambda = {}", lambda);

    // Estimate factor `n`.
    let poiss = Poisson::new(lambda).chain_err(|| {
        format!(
            "Could not construct Poisson distribution for lambda = {}",
            lambda
        )
    })?;
    let n = (y_1.ln() - poiss.ln_pmf(1)).exp();
    debug!(logger, " -> n = {}", lambda);

    // Compute expected false discoveries.
    let exp_fd = (0..k_lim)
        .map(|x| {
            if x == 0 {
                0.0
            } else {
                n * (1.0 - poiss.cdf(x as f64))
            }
        })
        .collect::<Vec<f64>>();
    debug!(logger, " -> exp_fd = {:?}", &exp_fd);

    // Compute observed discoveries.
    let mut obs_d = vec![0; k_lim as usize]; // UNUSED: `obs_d[0]`
    for i in 1..(k_max - 1) {
        // println!("i = {}", i);
        for j in (i + 1)..k_max {
            // println!("  -> j = {}", j);
            obs_d[i as usize] += y[j as usize];
        }
    }
    debug!(logger, " -> obs_d = {:?}", &obs_d);

    // Compute ratio of expected FPs over observed counts.
    let ratios = exp_fd
        .iter()
        .zip(obs_d.iter())
        .map(|(e, o)| {
            let o = *o as f64;
            if o != 0.0 {
                e / o as f64
            } else {
                0.0
            }
        })
        .collect::<Vec<f64>>();
    debug!(logger, " -> ratios = {:?}", &ratios);

    // Select smallest possible threshold with requested FDR.
    if let Some(chosen) = ratios
        .iter()
        .enumerate()
        .filter(|(i, r)| (*i > 0) && (**r < fdr_cutoff))
        .map(|(i, _)| i)
        .peekable()
        .peek()
    {
        let chosen = *chosen;
        Ok((k[chosen - 1] as f64
            + (ratios[chosen - 1] - fdr_cutoff) / (ratios[chosen - 1] - ratios[chosen]).ceil())
            as u32)
    } else {
        bail!("Could not find appropriate peak threshold");
    }
}

/// Collect piles.
///
/// This will compute the piles genom-wide, compute the depth threshold base on a Poisson
/// noise assumption and an FDR value given by the user on the command line, and then return
/// the piles with depth above this threshold.
fn gather_piles(
    logger: &mut Logger,
    options: &CoverageOptions,
    regions: &GenomeRegions,
    mut out_bed: Option<&mut File>,
) -> Result<HashMap<String, IntervalTree<u32, u32>>> {
    info!(logger, "Collecting piles...");
    let mut collector = piles::PileCollector::new(
        &options.input,
        options.pile_max_gap,
        options.min_mapq,
        options.io_threads,
    )?;
    let mut chrom_pile_infos = HashMap::new();
    for (chrom, _, _) in &regions.regions {
        debug!(logger, " -> for contig {}", chrom);
        chrom_pile_infos.insert(chrom.clone(), collector.collect_piles(&chrom, logger));
    }
    info!(logger, " => done");

    info!(logger, "Computing depth frequencies");
    let mut abs_freq: Vec<usize> = Vec::new();
    for (chrom, _, _) in &regions.regions {
        for pile in chrom_pile_infos.get(&chrom.to_string()).unwrap() {
            if abs_freq.len() <= (pile.size as usize) {
                abs_freq.resize(pile.size as usize + 1, 0);
            }
            abs_freq[pile.size as usize] += 1;
        }
    }
    info!(logger, " => done");

    info!(logger, "Computing depth threshold");
    let thresh = compute_threshold(logger, &abs_freq, options.mask_piles_fdr)?;
    info!(logger, "threshold = {}", thresh);

    let mut result = HashMap::new();
    for (chrom, pile_infos) in &chrom_pile_infos {
        let mut interval_tree = IntervalTree::new();
        let mut len_sum = 0;
        for info in pile_infos.iter() {
            if info.size >= thresh {
                interval_tree.insert(info.start..info.end, info.size);
                len_sum += info.end - info.start;
                if let Some(ref mut file) = out_bed {
                    file.write_all(format!("{}\t{}\t{}\n", chrom, info.start, info.end).as_bytes())
                        .chain_err(|| "Problem writing to BED file")?;
                }
            }
        }
        info!(
            logger,
            "Black-listing {} bp on contig {}",
            len_sum.separated_string(),
            chrom
        );
        result.insert(chrom.clone(), interval_tree);
    }
    info!(logger, " => done");

    Ok(result)
}

/// Process one region.
fn process_region(
    logger: &mut Logger,
    options: &CoverageOptions,
    (chrom, start, end): &(String, usize, usize),
    contig_length: usize,
    piles_tree: Option<&IntervalTree<u32, u32>>,
    bcf_writer: &mut bcf::Writer,
) -> Result<()> {
    info!(
        logger,
        "Processing region {}:{}-{}",
        chrom,
        (start + 1).separated_string(),
        end.separated_string()
    );

    // Load statistics for the given contig.
    // TODO: -> function
    let ref_stats = match (&options.reference, options.window_length) {
        (Some(reference), Some(window_length)) => {
            info!(logger, "Loading reference statistics from FASTA file");
            Some(ReferenceStats::from_path(
                reference,
                chrom,
                window_length as usize,
                logger,
            )?)
        }
        _ => {
            debug!(logger, "NOT loading reference statistics from FASTA file");
            None
        }
    };

    let target = if let Some(targets_bed) = &options.targets_bed {
        info!(logger, "Loading target regions from BED file");
        let mut tbx_reader =
            tbx::Reader::from_path(targets_bed).chain_err(|| "Could not open targets BED file")?;
        let regions = load_bed_regions_from_tabix(logger, &mut tbx_reader, &chrom, contig_length)?;
        // Build interval tree.
        let mut intervals: IntervalTree<u32, u32> = IntervalTree::new();
        for (i, (_, start, end)) in regions.regions.iter().enumerate() {
            let start = *start as u32;
            let end = *end as u32;
            intervals.insert(start..end, i as u32);
        }
        Some((regions, intervals))
    } else if let Some(model_bcf) = &options.model_bcf {
        info!(
            logger,
            "Loading target regions from WIS model BCF file file"
        );
        let regions = load_regions_from_model_bcf(logger, &model_bcf, &chrom, contig_length)?;
        // TODO: extract this into function
        // Build interval tree.
        let mut intervals: IntervalTree<u32, u32> = IntervalTree::new();
        for (i, (_, start, end)) in regions.regions.iter().enumerate() {
            let start = *start as u32;
            let end = *end as u32;
            intervals.insert(start..end, i as u32);
        }
        Some((regions, intervals))
    } else {
        info!(logger, "Not loading any targets");
        None
    };

    info!(logger, "Loading black-listed regions from BED file");
    let blacklist_tree = if let Some(blacklist_bed) = &options.blacklist_bed {
        let mut tbx_reader = tbx::Reader::from_path(blacklist_bed)
            .chain_err(|| "Could not open blacklist BED file")?;
        let regions = load_bed_regions_from_tabix(logger, &mut tbx_reader, &chrom, contig_length)?;
        // TODO: extract this into function
        // Build interval tree.
        let mut intervals: IntervalTree<usize, usize> = IntervalTree::new();
        for (i, (_, start, end)) in regions.regions.iter().enumerate() {
            let start = *start as usize;
            let end = *end as usize;
            intervals.insert(start..end, i);
        }
        Some(intervals)
    } else {
        None
    };

    // Construct aggregator for the records.
    info!(logger, "Constructing aggregator...");
    let mut aggregator: Box<BamRecordAggregator> =
        match (options.count_kind, options.considered_regions) {
            // TODO: next step is to integrate the target fragment counting here
            (CountKind::Fragments, ConsideredRegions::GenomeWide) => Box::new(
                FragmentsGenomeWideAggregator::new(piles_tree, options.clone(), *end),
            ),
            (CountKind::Fragments, ConsideredRegions::TargetRegions) => {
                let (target_regions, target_tree) = target.unwrap();
                Box::new(FragmentsTargetRegionsAggregator::new(
                    target_tree.clone(),
                    target_regions.clone(),
                    options.clone(),
                    contig_length,
                ))
            }
            (CountKind::Coverage, ConsideredRegions::GenomeWide) => {
                Box::new(CoverageAggregator::new(options.clone(), *end))
            }
            (_, _) => bail!("Invalid combination of coverage/on-target regions"),
        };

    // TODO: 2 blocks -> function
    // Jump to region with BAM reader.
    let mut bam_reader =
        bam::IndexedReader::from_path(&options.input).chain_err(|| "Could not open BAM file")?;
    if options.io_threads > 0 {
        bam_reader
            .set_threads(options.io_threads as usize)
            .chain_err(|| "Could not set threads on BAM reader")?;
    }
    let tid: u32 = bam_reader.header().tid(chrom.as_bytes()).unwrap();
    bam_reader
        .fetch(tid, *start as u32, *end as u32)
        .chain_err(|| format!("Could not seek to region {}:{}-{}", chrom, start, end))?;

    // Main loop for region: pass all BAM records in region through aggregator.
    info!(logger, "Computing coverage...");
    aggregator.put_fetched_records(&mut bam_reader);
    debug!(
        logger,
        "Processed {}, skipped {} records ({:.2}% were processed)",
        aggregator.num_processed().separated_string(),
        aggregator.num_skipped().separated_string(),
        100.0 * (aggregator.num_processed() - aggregator.num_skipped()) as f64
            / aggregator.num_processed() as f64,
    );

    // TODO: full block -> function
    // Create the BCF records for this region.
    info!(logger, "Writing BCF with coverage information...");
    for region_id in 0..aggregator.num_regions() {
        let stats = aggregator.get_stats(region_id);
        let mut record = bcf_writer.empty_record();
        let rid = bcf_writer
            .header()
            .name2rid(chrom.as_bytes())
            .chain_err(|| format!("Could not find contig {}", chrom))?;

        // Columns: CHROM, POS, ID, REF, ALT, (FILTER)
        let pos = stats.start;
        let window_end = stats.end;
        let alleles_v = vec![
            Vec::from("N"),
            if options.considered_regions == ConsideredRegions::TargetRegions {
                Vec::from("<TARGET>")
            } else {
                Vec::from("<WINDOW>")
            },
        ];
        let alleles = alleles_v
            .iter()
            .map(|x| x.as_slice())
            .collect::<Vec<&[u8]>>();

        record.set_rid(&Some(rid));
        record.set_pos(pos as i32);
        record
            .set_id(format!("{}:{}-{}", chrom, pos + 1, window_end).as_bytes())
            .chain_err(|| "Could not update ID")?;
        record
            .set_alleles(&alleles)
            .chain_err(|| "Could not update alleles")?;

        // Columns: INFO
        record
            .push_info_integer(b"END", &[window_end as i32])
            .chain_err(|| "Could not write INFO/END")?;
        if let Some(ref_stats) = &ref_stats {
            if ref_stats.has_gap[region_id] {
                record
                    .push_info_flag(b"GAP")
                    .chain_err(|| "Could not write INFO/GAP")?;
            }
            record
                .push_info_float(b"GC", &[ref_stats.gc_content[region_id] as f32])
                .chain_err(|| "Could not write INFO/GC: {}")?;
        }
        if let Some(ref blacklist) = &blacklist_tree {
            if blacklist.find(pos..window_end).peekable().peek().is_some() {
                record
                    .push_info_flag(b"BLACKLIST")
                    .chain_err(|| "Could not write INFO/BLACKLIST: {}")?;
            }
        }

        // Vector for collecting the filter values.
        let mut fts: Vec<String> = Vec::new();

        // Columns: FORMAT/GT
        record
            .push_format_integer(b"GT", &[bcf::GT_MISSING])
            .chain_err(|| "Could not write FORMAT/GT: {}")?;

        // Columns: FORMAT/CV etc.
        record
            .push_format_float(b"RCV", &[stats.cov])
            .chain_err(|| "Could not write FORMAT/RCV")?;
        if options.considered_regions == ConsideredRegions::TargetRegions
            && stats.cov < options.min_raw_coverage as f32
        {
            fts.push("LOWCOV".to_string());
        }
        if let Some(cov_sd) = stats.cov_sd {
            record
                .push_format_float(b"RCVSD", &[cov_sd])
                .chain_err(|| "Could not write FORMAT/RCVSD")?;
        }
        let length = (stats.end - stats.start) as f32;
        record
            .push_format_float(b"LCV", &[stats.cov / length])
            .chain_err(|| "Could not write FORMAT/CV")?;
        if let Some(cov_sd) = stats.cov_sd {
            record
                .push_format_float(b"LCVSD", &[cov_sd / length])
                .chain_err(|| "Could not write FORMAT/CVSD")?;
        }

        if let Some(frac_removed) = stats.frac_removed {
            record
                .push_format_float(b"FRM", &[frac_removed])
                .chain_err(|| "Could not write FORMAT/FRM")?;
            if frac_removed > options.min_window_remaining {
                fts.push("PILE".to_string());
            }
        }

        record
            .push_format_float(b"MQ", &[stats.mean_mapq])
            .chain_err(|| "Could not write FORMAT/MQ")?;

        if !fts.is_empty() {
            let value = fts.join(";");
            record
                .push_format_char(b"FT", value.as_bytes())
                .chain_err(|| "Could not write FORMAT/PILE")?;
        }

        bcf_writer
            .write(&record)
            .chain_err(|| "Could not write BCF record")?;
    }

    Ok(())
}

/// Main entry point for the `coverage` command.
pub fn run(logger: &mut Logger, options: &CoverageOptions) -> Result<()> {
    info!(logger, "Running: cnvetti cmd coverage");
    info!(logger, "Options: {:?}", options);

    // Warn about unusual parametrizations.
    if options.window_length.is_some() && options.reference.is_none() {
        warn!(
            logger,
            "You are counting per-window coverage but have specified no reference FASTA file.  \
             This is required for GC correction and usually desired.  Please reconsider your \
             choice!"
        );
    }

    // Bail out on coverage on target regions (invalid combination).
    if options.considered_regions == ConsideredRegions::TargetRegions
        && options.count_kind == CountKind::Coverage
    {
        bail!(
            "Cannot combine target regions with considering coverage. Either consider \
             genome-wide coverage or consider fragment counts in target regions."
        );
    }

    // Get list of regions to process.
    let ref regions = match &options.genome_region {
        Some(regions) => GenomeRegions::from_string_list(&vec![regions.clone()])
            .chain_err(|| "Problem parsing region")?,
        None => {
            let bam_reader = bam::IndexedReader::from_path(&options.input)
                .chain_err(|| "Could not open BAI-indexed BAM file")?;
            bam_utils::build_chroms_bam(bam_reader.header(), Some(options.contig_regex.clone()))
                .chain_err(|| "Problem getting contig-length list from BAM header")?
        }
    };
    // Print list of regions that will be processed.
    info!(logger, "Will process {} regions", regions.regions.len());
    debug!(
        logger,
        "Regions: {:?}",
        regions
            .regions
            .iter()
            .map(|(name, _, len)| format!("{}:1-{}", name, len.separated_string()))
            .collect::<Vec<String>>()
    );

    let mut out_masked = if let Some(ref output_masked_bed) = options.output_masked_bed {
        info!(logger, "Creating BED file with masked regions");
        Some(File::create(output_masked_bed).chain_err(|| "Could not create BED file")?)
    } else {
        info!(logger, "NOT creating BED file with masked regions");
        None
    };

    // Create output file writer and kick off processing.  This is done in its own block such
    // that the file is definitely closed when building the index below.
    {
        let contigs = {
            let bam_reader = bam::IndexedReader::from_path(&options.input)
                .chain_err(|| "Could not open BAI-indexed BAM file")?;
            bam_utils::build_chroms_bam(bam_reader.header(), Some(options.contig_regex.clone()))
                .chain_err(|| "Problem getting contig-length list from BAM header")?
        };
        let mut contig_lengths: HashMap<String, usize> = HashMap::new();
        for (chrom, _, end) in &contigs.regions {
            contig_lengths.insert(chrom.clone(), *end);
        }

        let samples = bam_utils::samples_from_file(&options.input)
            .chain_err(|| "Could not get sample names from BAM file")?;
        let mut writer = build_bcf_writer(&options.output, &samples, &contigs)?;

        // TODO: add back piles BAM output file?

        info!(logger, "Collecting read piles and computing FDR...");
        let piles_by_chrom = if options.mask_piles {
            Some(gather_piles(
                logger,
                &options,
                &regions,
                out_masked.as_mut(),
            )?)
        } else {
            None
        };
        info!(logger, " => done");

        // Process each region.
        for region in &regions.regions {
            process_region(
                logger,
                &options,
                &region,
                *contig_lengths
                    .get(&*region.0)
                    .expect("Could not resolve region"),
                piles_by_chrom
                    .as_ref()
                    .map(|m| m.get(&*region.0).expect("Unknown contig")),
                &mut writer,
            )?;
        }
    }

    // Finally, create index on created output file.
    info!(logger, "Building index for output file...");
    bcf_utils::build_index(logger, &options.output).chain_err(|| "Could not build index")?;
    info!(logger, "All done. Have a nice day!");

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_threshold() {
        let y = vec![0, 10, 8, 7, 6, 5, 3];
        // assert_eq!(7, compute_threshold(&y, 0.001).expect("No threshold!"));
    }
}
