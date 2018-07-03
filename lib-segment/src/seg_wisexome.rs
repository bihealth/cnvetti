//! Implementation of segmentation algorithm from WISExome.

use std::env;

use rust_htslib::bcf::{self, Read};
use shlex;
use slog::Logger;

use super::errors::*;
use options::*;

use lib_shared::bcf_utils;
use lib_shared::stats::Stats;

use ordered_float::OrderedFloat;

/// Epsilon to add normalized metric to prevent -nan through log2()
const PSEUDO_EPSILON: f64 = 1e-12;

/// Whether or not to skip the record.
fn skip_record(record: &mut bcf::Record) -> bool {
    if !record.format(b"CV").float().is_ok() || !record.format(b"CVZ").float().is_ok() {
        return true;
    }

    let cv = record
        .format(b"CV")
        .float()
        .expect("Could not access FORMAT/CV")[0][0] as f64;
    let cvz = record
        .format(b"CVZ")
        .float()
        .expect("Could not access FORMAT/CVZ")[0][0] as f64;
    let is_gap = record.info(b"GAP").flag().unwrap_or(false);

    !cv.is_finite() || !cvz.is_finite() || is_gap
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

    // Coverage relative to reference targets; called "effect size" by Straver et al.
    cov_rel: f64,
    // Window-wide "effect size", starts out as window coverage.
    window_rel: f64,
    // Z score of cov_rel with respect to reference targets.
    z_score: f64,

    // Call state.
    call_state: CallState,
}

impl ProbeInfo {
    /// Update call state based on probe effect size and z-score.
    fn update_call_state_probe(&mut self, thresh_rel: f64, thresh_z: f64) {
        self.call_state = if (self.cov_rel - 1.0).abs() < thresh_rel {
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

/// Segmentation using the algorithm from the WISExome paper.
fn seg_wisexome(
    logger: &mut Logger,
    probe_infos: &mut Vec<ProbeInfo>,
    options: &SegmentOptions,
) -> Vec<ProbeInfo> {
    // Perform windowing and z-score updates.
    //
    // First compute z scores of varying size windows.
    let n = probe_infos.len();
    let m = (options.wisexome_max_window_size / 2 + 1) as usize;
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
        probe_infos[i].update_call_state_probe(
            options.wisexome_thresh_rel_cov,
            options.wisexome_thresh_z_score,
        );
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
                    .map(|info| info.cov_rel)
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
                        .map(|info| info.cov_rel)
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
            .map(|info| info.cov_rel)
            .collect::<Vec<f64>>()
            .as_slice()
            .median();
        for i in start..end {
            probe_infos[i].window_rel = window_rel;
        }
    }
    // Then, update call state based on window's effect size.
    for info in probe_infos.iter_mut() {
        info.update_call_state_window(options.wisexome_thresh_rel_cov);
    }
    return probe_infos.clone();

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
                        .map(|info| info.cov_rel)
                        .collect::<Vec<f64>>()
                        .as_slice()
                        .mean();

                    if (1.0 - window_rel).abs() < options.wisexome_thresh_rel_cov {
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
                .map(|info| info.cov_rel)
                .collect::<Vec<f64>>()
                .as_slice()
                .median();
            for info in &mut probe_infos[start..end] {
                info.window_rel = window_rel;
                info.call_state = call_state;
                info.update_call_state_window(options.wisexome_thresh_rel_cov);
            }
        }
    }

    probe_infos.clone()
}

/// Perform segmentation using the algorithm from the WISExome paper.
pub fn run_segmentation(logger: &mut Logger, options: &SegmentOptions) -> Result<()> {
    info!(logger, "Computing segmentation using Haar-Seg");

    debug!(logger, "Opening input file");
    let mut reader = bcf::IndexedReader::from_path(options.input.clone())
        .chain_err(|| "Could not open input BCF file")?;

    debug!(logger, "Opening output file");
    let mut writer = {
        // Construct extended header.
        let mut header = bcf::Header::from_template(reader.header());
        // let lines = vec![];
        // for line in lines {
        //     header.push_record(line.as_bytes());
        // }
        header.push_record(format!("##cnvetti_segmentVersion={}", "0.1.0").as_bytes());
        header.push_record(
            format!(
                "##cnvetti_segmentCommand={}",
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
            .chain_err(|| "Could not open output BCF file")?
    };
    if options.io_threads > 0 {
        writer
            .set_threads(options.io_threads as usize)
            .chain_err(|| "Could not set I/O thread count")?;
    }

    debug!(logger, "Performing the actual work.");
    let chroms = bcf_utils::extract_chroms(&reader.header());
    for (chrom, _, chrom_len) in &chroms.regions {
        debug!(logger, "Processing chrom {}", &chrom);

        trace!(logger, "Jumping in coverage BCF file");
        let rid = reader
            .header()
            .name2rid(chrom.as_bytes())
            .chain_err(|| format!("Could not translate header name {}", chrom))?;
        reader
            .fetch(rid, 0, *chrom_len as u32)
            .chain_err(|| format!("Could not fetch chromosome {}", chrom))?;

        trace!(logger, "Loading probe infos");
        let mut probe_infos: Vec<ProbeInfo> = Vec::new();
        let mut record = reader.empty_record();
        loop {
            match reader.read(&mut record) {
                Ok(_) => (),
                Err(bcf::ReadError::NoMoreRecord) => break,
                _ => bail!("Error reading BCF record"),
            }

            if !skip_record(&mut record) {
                let cv = record
                    .format(b"CV")
                    .float()
                    .expect("Could not access FORMAT/CV")[0][0] as f64;
                let cvz = record
                    .format(b"CVZ")
                    .float()
                    .expect("Could not access FORMAT/CVZ")[0][0] as f64;

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
                    cov_rel: cv,
                    window_rel: cv,
                    z_score: cvz,
                    call_state: CallState::Ref,
                });
            }
        }
        trace!(logger, " => done");

        trace!(logger, "Performing segmentation...");
        let segmentation = seg_wisexome(logger, &mut probe_infos, options);
        trace!(logger, " => done");

        // Compute empirical p-values of for each window.
        info!(logger, "(Not computing empirical p values yet...)");
        // TODO

        trace!(logger, "Write out segmentation");
        let mut idx = 0; // current index into ncov
        let mut prev_val: Option<(f32, f32, f32)> = None; // (val, log2(val), z_score)
        reader
            .fetch(rid, 0, *chrom_len as u32)
            .chain_err(|| format!("Could not fetch chromosome {}", chrom))?;
        while reader.read(&mut record).is_ok() {
            writer.translate(&mut record);
            let (val, val_log2, z_score) = if !skip_record(&mut record) {
                record.push_filter(writer.header().name_to_id(b"SKIPPED_SEG").unwrap());
                if let Some(prev_val) = prev_val {
                    prev_val
                } else {
                    (1.0, 0.0, 1.0)
                }
            } else {
                let val = if segmentation[idx].window_rel > 0.0
                    && segmentation[idx].window_rel < PSEUDO_EPSILON
                {
                    (1.0, 0.0, 1.0)
                } else {
                    (
                        (segmentation[idx].window_rel - PSEUDO_EPSILON) as f32,
                        segmentation[idx].window_rel.log2() as f32,
                        segmentation[idx].z_score as f32,
                    )
                };
                prev_val = Some(val);
                idx += 1;
                val
            };
            record
                .push_format_float(b"SG", &[val])
                .chain_err(|| "Could not write FORMAT/SG")?;
            record
                .push_format_float(b"SG2", &[val_log2])
                .chain_err(|| "Could not write FORMAT/SG2")?;
            record
                .push_format_float(b"SGZ", &[z_score])
                .chain_err(|| "Could not write FORMAT/SGZ")?;
            writer.write(&record).expect("Writing the record failed!");
        }
    }

    info!(logger, "=> OK");
    Ok(())
}
