/// Main library module for building the WIS model.
use std::cell::RefCell;
use std::cmp::min;
use std::env;
use std::f32;
use std::sync::atomic::{AtomicUsize, Ordering};

extern crate clap;

extern crate chrono;

#[macro_use]
extern crate error_chain;

extern crate pdqselect;

#[macro_use]
extern crate slog;
use slog::Logger;

extern crate rust_htslib;
use rust_htslib::bcf::{self, Read};

extern crate shlex;

extern crate rayon;
use rayon::prelude::*;

extern crate separator;
use separator::Separatable;

extern crate lib_shared;
use lib_shared::bcf_utils;
use lib_shared::regions::GenomeRegions;
use lib_shared::stats::Stats;

mod options;
pub use options::*;

/// This crate's error-related code, generated by `error-chain`.
mod errors {
    // Create the Error, ErrorKind, ResultExt, and Result types
    error_chain!{}
}

pub use errors::*;

/// Structure for representing the data extracted from input BCF file.
struct ModelInput {
    /// Contigs from input BCF file.
    pub contigs: GenomeRegions,
    /// Number of samples.
    pub num_samples: usize,
    /// Number of target regions.
    pub num_regions: usize,
    /// The genomic regions for which coverage is computed.
    pub regions: GenomeRegions, // TODO: rename to targets?
    /// For each `GenomicRegions` entry, the numeric reference/contig ID.
    pub tids: Vec<u32>, // TODO: rename
    /// For each region, the normalized coverage/fragment count per sample.
    pub ncovs: Vec<Vec<f32>>,
}

/// Load the normalized fragment counts, return region-wise vector of counts.
fn read_coverage_bcf(logger: &mut Logger, options: &BuildModelWisOptions) -> Result<ModelInput> {
    info!(logger, "Opening input file...");
    let mut reader =
        bcf::Reader::from_path(&options.input).chain_err(|| "Could not open input BCF file.")?;
    if options.io_threads > 0 {
        reader
            .set_threads(options.io_threads as usize)
            .chain_err(|| "Could not set threads to BCF reader.")?;
    }

    let mut regions = GenomeRegions::new();
    let mut tids = Vec::new();
    let mut ncovs = Vec::new();

    let mut record = reader.empty_record();
    loop {
        match reader.read(&mut record) {
            Ok(_) => (),
            Err(bcf::ReadError::NoMoreRecord) => break,
            _ => bail!("Could not read record"),
        }

        regions.regions.push((
            String::from_utf8(reader.header().rid2name(record.rid().unwrap()).to_vec())
                .chain_err(|| "Could not parse from UTF-8")?,
            record.pos() as usize,
            record.info(b"END").integer().unwrap().unwrap()[0] as usize,
        ));

        tids.push(record.rid().unwrap());

        let mut tmp = Vec::new();
        for ncov in record
            .format(b"CV")
            .float()
            .chain_err(|| "Could not read FORMAT/CV.")?
        {
            tmp.push(ncov[0]);
        }
        ncovs.push(tmp);
    }

    info!(
        logger,
        " => number of targets in BCF file: {}",
        regions.regions.len().separated_string()
    );

    let contigs = bcf_utils::extract_chroms(&reader.header());
    let num_samples = reader.header().sample_count() as usize;
    let num_regions = regions.regions.len() as usize;

    Ok(ModelInput {
        contigs,
        num_samples,
        num_regions,
        regions,
        tids,
        ncovs,
    })
}

/// Compute the distance matrix (in parallel).
fn compute_distances(
    logger: &mut Logger,
    tids: &Vec<u32>,
    ncovs: &Vec<Vec<f32>>,
    options: &BuildModelWisOptions,
) -> Vec<Vec<(f32, usize)>> {
    // Set number of threads to use by rayon.
    env::set_var("RAYON_NUM_THREADS", format!("{}", options.num_threads));

    info!(
        logger,
        "Computing distances (with {} threads)...", options.num_threads
    );

    let progress = AtomicUsize::new(0);

    let result = ncovs
        .par_iter()
        .enumerate()
        .map(|(i, row)| {
            thread_local!{static DIST_RC: RefCell<Vec<(f32, usize)>> = RefCell::new(Vec::new());}

            DIST_RC.with(|dist_rc| {
                let mut dist = dist_rc.borrow_mut();
                dist.clear();
                for (j, row2) in ncovs.iter().enumerate() {
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
        .collect::<Vec<Vec<(f32, usize)>>>();

    info!(logger, " => done");
    result
}

/// Prune distances, return top partner targets for each target.
fn prune_distances(
    logger: &mut Logger,
    distances: &Vec<Vec<(f32, usize)>>,
    options: &BuildModelWisOptions,
) -> (Vec<Vec<usize>>, Vec<Vec<f32>>) {
    info!(logger, "Pruning distances...");
    // Compute distance error statistics.
    let error_stats: Vec<f64> = distances
        .iter()
        .map(|xs| xs[0].0 as f64)
        .collect::<Vec<f64>>();
    if error_stats.len() < 2 {
        return (Vec::new(), Vec::new());
    }
    let median = error_stats.as_slice().median();
    let mad = error_stats.as_slice().median_abs_dev();
    let threshold = (median + options.filter_z_score * mad) as f32;

    // Actually collect the reliable regions and flag for reliability.
    let mut ref_targets: Vec<Vec<usize>> = Vec::new();
    let mut ref_distances: Vec<Vec<f32>> = Vec::new();

    for row in distances.iter() {
        ref_targets.push(
            row.iter()
                .filter(|(dist, _)| *dist <= threshold)
                .map(|(_, idx)| *idx)
                .collect::<Vec<usize>>(),
        );
        ref_distances.push(
            row.iter()
                .filter(|(dist, _)| *dist <= threshold)
                .map(|(dist, _)| *dist)
                .collect::<Vec<f32>>(),
        );
    }
    info!(logger, " => done");

    (ref_targets, ref_distances)
}

/// Compute per-probe calls for determining UNRELIABLE flag.
fn count_per_probe_calls(
    logger: &mut Logger,
    num_regions: usize,
    num_samples: usize,
    ref_targets: &Vec<Vec<usize>>,
    ncovs: &Vec<Vec<f32>>,
    options: &BuildModelWisOptions,
) -> Vec<u32> {
    info!(logger, "Count per-probe calls for UNRELIABLE flag...");
    let mut num_calls = vec![0; num_regions];

    for (target_idx, matched_targets) in ref_targets.iter().enumerate() {
        for sample_idx in 0..num_samples {
            // Get empirical reference distribution for target `target_idx` for sample `sample_idx`.
            let ref_dist = matched_targets
                .iter()
                .map(|idx| ncovs[*idx][sample_idx] as f64)
                .collect::<Vec<f64>>();
            if ref_dist.len() < options.min_ref_targets {
                continue;
            }
            let ref_dist = ref_dist.as_slice();
            // Compute z-score
            let x = ncovs[target_idx][sample_idx] as f64;
            let z_score = (x - ref_dist.mean()) / ref_dist.std_dev();
            // Compute relative value.
            let rel = x / ref_dist.mean();
            // Increment per-probe counter for target.
            if z_score.abs() > options.filter_z_score && (rel - 1.0).abs() > options.filter_rel {
                num_calls[target_idx] += 1;
            }
        }
    }
    info!(logger, " => done");

    num_calls
}

/// Build header for the output BCF file.
///
/// This defines all values used throughout the whole window/target specific BAM files,
/// regardless whether they are actually used in the file.
///
/// Note that we use shared FORMAT tags for coverage and fragment count, such that we get
/// unified processing of copy number data in BCF files.
fn build_header(contigs: &GenomeRegions) -> bcf::Header {
    let mut header = bcf::Header::new();

    // Put overall meta information into the BCF header.
    let now = chrono::Utc::now();
    header.push_record(format!("##fileDate={}", now.format("%Y%m%d").to_string()).as_bytes());

    // Put creating tool version and call into file.
    header.push_record(format!("##cnvetti_cmdBuildModelWisVersion={}", "0.1.0").as_bytes());
    header.push_record(
        format!(
            "##cnvetti_cmdBuildModelWisCommand={}",
            env::args()
                .map(|s| shlex::quote(&s).to_string())
                .collect::<Vec<String>>()
                .join(" ")
        ).as_bytes(),
    );

    // Note: no samples

    // Put contig information into BCF header.
    for (name, _, length) in &contigs.regions {
        header.push_record(format!("##contig=<ID={},length={}>", name, length).as_bytes());
    }

    // Push the relevant header records.
    let lines = vec![
        // Define ALT column <TARGET>
        "##ALT=<ID=TARGET,Description=\"Record describes a target for fragment counting\">",
        // INFO fields describing the window
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Window end\">",
        "##INFO=<ID=NUM_CALLS,Number=1,Type=Integer,Description=\"Number of samples with \
         calls\">",
        "##INFO=<ID=REF_TARGETS,Number=.,Type=String,Description=\"List of target regions that \
         are used as within-sample reference.\">",
        "##INFO=<ID=REF_TARGETS2,Number=.,Type=String,Description=\"List of target regions that \
         are used as within-sample reference plus distances.\">",
        // FILTER fields
        "##FILTER=<ID=FEW_REF,Description=\"Target is masked because there were not enough
         reference targets.\">",
        "##FILTER=<ID=UNRELIABLE,Description=\"Target is masked because it was called in
         too many reference panel samples.\">",
    ];
    for line in lines {
        header.push_record(line.as_bytes());
    }

    header
}

/// Build bcf::Writer with appropriate header.
fn build_bcf_writer(path: &String, ref_contigs: &GenomeRegions) -> Result<bcf::Writer> {
    let uncompressed = !path.ends_with(".bcf") && !path.ends_with(".vcf.gz");
    let vcf = path.ends_with(".vcf") || path.ends_with(".vcf.gz");

    let header = build_header(ref_contigs);
    bcf::Writer::from_path(&path, &header, uncompressed, vcf)
        .chain_err(|| "Could not open BCF file for writing")
}

fn write_output(
    logger: &mut Logger,
    model_input: &ModelInput,
    num_calls: &Vec<u32>,
    ref_targets: &Vec<Vec<usize>>,
    ref_distances: &Vec<Vec<f32>>,
    options: &BuildModelWisOptions,
) -> Result<()> {
    // Actually collect the reliable regions and flag for reliability.
    info!(logger, "Writing output BCF file...");
    {
        let mut writer = build_bcf_writer(&options.output, &model_input.contigs)
            .chain_err(|| "Could not open output BCF file for writing")?;

        for (i, tgt) in model_input.regions.regions.iter().enumerate() {
            let mut record = writer.empty_record();

            record.set_rid(&Some(model_input.tids[i]));
            record.set_pos(tgt.1 as i32);
            record
                .set_id(format!("{}:{}-{}", tgt.0, tgt.1 + 1, tgt.2).as_bytes())
                .chain_err(|| "Could not write ID")?;

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
                .chain_err(|| "Could not update alleles")?;
            record
                .push_info_integer(b"END", &[tgt.2 as i32])
                .chain_err(|| "Could not write INFO/END")?;
            record
                .push_info_integer(b"NUM_CALLS", &[num_calls[i] as i32])
                .chain_err(|| "Could not write INFO/NUM_CALLS")?;

            let ref_targets2 = ref_targets[i]
                .iter()
                .map(|idx| {
                    format!(
                        "{}:{}-{}",
                        model_input.regions.regions[*idx].0,
                        model_input.regions.regions[*idx].1 + 1,
                        model_input.regions.regions[*idx].2
                    )
                })
                .collect::<Vec<String>>();
            let ref_targets_b = ref_targets2
                .iter()
                .map(|s| s.as_bytes())
                .collect::<Vec<&[u8]>>();
            record
                .push_info_string(b"REF_TARGETS", ref_targets_b.as_slice())
                .chain_err(|| "Could not write INFO/REF_TARGETS")?;

            let ref_targets2 = ref_targets[i]
                .iter()
                .enumerate()
                .map(|(j, idx)| {
                    format!(
                        "{}:{}-{}({})",
                        model_input.regions.regions[*idx].0,
                        model_input.regions.regions[*idx].1 + 1,
                        model_input.regions.regions[*idx].2,
                        ref_distances[i][j]
                    )
                })
                .collect::<Vec<String>>();
            let ref_targets_b = ref_targets2
                .iter()
                .map(|s| s.as_bytes())
                .collect::<Vec<&[u8]>>();
            record
                .push_info_string(b"REF_TARGETS2", ref_targets_b.as_slice())
                .chain_err(|| "Could not write INFO/REF_TARGETS2")?;

            writer
                .write(&record)
                .chain_err(|| "Problem writing to output")?;
        }
    }

    Ok(())
}

/// Main entry point for the "cmd build-model-wis" command.
pub fn run(logger: &mut Logger, options: &BuildModelWisOptions) -> Result<()> {
    info!(logger, "Running: cnvetti cmd build-model-wis");
    info!(logger, "Options: {:?}", options);

    let model_input = read_coverage_bcf(logger, &options)?;
    let distances = compute_distances(logger, &model_input.tids, &model_input.ncovs, &options);
    let (ref_targets, ref_distances) = prune_distances(logger, &distances, &options);
    let num_calls = count_per_probe_calls(
        logger,
        model_input.num_regions,
        model_input.num_samples,
        &ref_targets,
        &model_input.ncovs,
        options,
    );

    // Write the output file.
    write_output(
        logger,
        &model_input,
        &num_calls,
        &ref_targets,
        &ref_distances,
        options,
    )?;

    // Finally, create index on created output file.
    info!(logger, "Building index for output file...");
    bcf_utils::build_index(logger, &options.output).chain_err(|| "Could not build index")?;
    info!(logger, "All done. Have a nice day!");

    Ok(())
}
