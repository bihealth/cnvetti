//! Genotyping using the XHMM approach.

// TODO: we're using chromosome-wide targets here instead of exome-wide targets as in XHMM paper.
// TODO: can we really do this?

use std::env;
use std::ops::Range;

use bio::stats::hmm::{backward, forward, Model, State};
use bio::stats::probs::{LogProb, PHREDProb};
use chrono;
use ndarray::prelude::*;
use rust_htslib::bcf::{self, Read};
// use rust_segment::{reject_nonaberrant_pvalue, seg_haar};
// use separator::Separatable;
use shlex;
use slog::Logger;

use super::errors::*;
use options::*;

use lib_segment::seg_xhmm::xhmm;
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

        let (sv_type, sv_len) = if cnv_info.quals.copy_state == CopyState::Deletion {
            ("DEL", -(cnv_info.segment.range.len() as i32))
        } else {
            ("DUP", cnv_info.segment.range.len() as i32)
        };

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
            .push_info_integer(b"CENTER", &[((pos + end) / 2) as i32])
            .chain_err(|| "Could not write INFO/END")?;
        record
            .push_info_string(b"SVTYPE", &[sv_type.as_bytes()])
            .chain_err(|| "Could not write INFO/SVTYPE")?;
        record
            .push_info_integer(b"SVLEN", &[sv_len])
            .chain_err(|| "Could not write INFO/SVTYPE")?;

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
            .push_format_float(b"EQ", &[*cnv_info.quals.q_exact_cnv as f32])
            .chain_err(|| "Could not write FORMAT/EQ")?;
        record
            .push_format_float(b"SQ", &[*cnv_info.quals.q_some_cnv as f32])
            .chain_err(|| "Could not write FORMAT/SQ")?;
        record
            .push_format_float(b"NQ", &[*cnv_info.quals.q_no_cnv as f32])
            .chain_err(|| "Could not write FORMAT/NQ")?;
        record
            .push_format_float(b"LQ", &[*cnv_info.quals.q_left_cnv_bp as f32])
            .chain_err(|| "Could not write FORMAT/LQ")?;
        record
            .push_format_float(b"RQ", &[*cnv_info.quals.q_right_cnv_bp as f32])
            .chain_err(|| "Could not write FORMAT/RQ")?;
        record
            .push_format_float(b"NDQ", &[*cnv_info.quals.q_right_cnv_bp as f32])
            .chain_err(|| "Could not write FORMAT/NDQ")?;
        record
            .push_format_float(b"DQ", &[*cnv_info.quals.q_right_cnv_bp as f32])
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

/// Struct with quality values and other metrics from XHMM describing an CNV.
#[derive(Clone, Debug)]
struct XhmmCnvQuals {
    // The copy state assumed below for "CNV/Cnv".
    copy_state: CopyState,

    /// Probability for "exact CNV".
    q_exact_cnv: PHREDProb,
    /// Probability for "some CNV".
    q_some_cnv: PHREDProb,
    /// Probability for "no CNV".
    q_no_cnv: PHREDProb,
    /// Probability for "left CNV breakpoint".
    q_left_cnv_bp: PHREDProb,
    /// Probability for "right CNV breakpoint".
    q_right_cnv_bp: PHREDProb,

    /// Probability for "is not".
    q_not_diploid: PHREDProb,
    /// Probability for "is diploid".
    q_diploid: PHREDProb,

    /// Mean normalized coverage throughout the CNV.
    mean_cov: f64,
    /// Mean Z-score of read depth through the CNV.
    mean_z_score: f64,
    /// Number of target region in the cnv.
    num_target_regions: usize,
}

/// Description of a CNV as called by XHMM.
#[derive(Clone, Debug)]
pub struct CnvGenotypeInfo {
    /// Segment with region information etc.
    segment: Segment,
    /// Qualities associated with the CNV.
    quals: XhmmCnvQuals,
}

use std::fmt::Debug;

/// Restricted version of the forward algorithm, sets probability to `0.0` for `impossible_state`,
/// states in `impossible_range`.
pub fn restricted_forward<O, M: Model<O>>(
    hmm: &M,
    observations: &[O],
    impossible_state: &State,
    impossible_range: &Range<usize>,
) -> (Array2<LogProb>, LogProb)
where
    O: Debug,
{
    // The matrix with probabilities.
    let mut vals = Array2::<LogProb>::zeros((observations.len(), hmm.num_states()));

    // Compute matrix.
    for (i, o) in observations.iter().enumerate() {
        if i == 0 {
            // Initial column.
            for s in hmm.states() {
                if (i >= impossible_range.start)
                    && (i < impossible_range.end)
                    && (s == *impossible_state)
                {
                    vals[[0, *s]] = LogProb::ln_zero();
                } else {
                    vals[[0, *s]] = hmm.initial_prob(s) + hmm.observation_prob(s, o);
                }
            }
        } else {
            // Subsequent columns.
            for j in hmm.states() {
                if (i >= impossible_range.start)
                    && (i < impossible_range.end)
                    && (j == *impossible_state)
                {
                    vals[[i, *j]] = LogProb::ln_zero();
                // println!(
                //     "i = {:?}, *j= {:?}, vals[[i, *j]] = {:?}",
                //     i,
                //     *j,
                //     &vals[[i, *j]]
                // );
                } else {
                    let xs = hmm.states()
                        .map(|k| {
                            // println!(
                            //     "i = {:?}, *j = {:?}, *k = {:?}, vals[[i - 1, *k]] = {:?}, \
                            //      hmm.transition_prob_idx(k, j, i) = {:?}, \
                            //      hmm.observation_prob(j, o) = {:?}, o = {:?}",
                            //     &i,
                            //     &*j,
                            //     &*k,
                            //     &vals[[i - 1, *k]],
                            //     &hmm.transition_prob_idx(k, j, i),
                            //     &hmm.observation_prob(j, o),
                            //     &o
                            // );
                            vals[[i - 1, *k]]
                                + hmm.transition_prob_idx(k, j, i)
                                + hmm.observation_prob(j, o)
                        })
                        .collect::<Vec<LogProb>>();
                    vals[[i, *j]] = LogProb::ln_sum_exp(&xs);
                    // println!(
                    //     "i = {:?}, *j= {:?}, xs = {:?}, vals[[i, *j]] = {:?}",
                    //     i,
                    //     *j,
                    //     &xs,
                    //     &vals[[i, *j]]
                    // );
                }
            }
        }
    }

    // Compute final probability.
    let prob = LogProb::ln_sum_exp(vals.row(observations.len() - 1).into_slice().unwrap());

    (vals, prob)
}

/// Restricted version of the backward algorithm, sets probability to `0.0` for `impossible_state`,
/// states in `impossible_range`.
pub fn restricted_backward<O, M: Model<O>>(
    hmm: &M,
    observations: &[O],
    impossible_state: &State,
    impossible_range: &Range<usize>,
) -> (Array2<LogProb>, LogProb) {
    // The matrix with probabilities.
    let mut vals = Array2::<LogProb>::zeros((observations.len(), hmm.num_states()));

    // Compute matrix.
    let n = observations.len();
    for (i, o) in observations.iter().rev().enumerate() {
        if i == 0 {
            for j in hmm.states() {
                if (i >= impossible_range.start)
                    && (i < impossible_range.end)
                    && (j == *impossible_state)
                {
                    vals[[0, *j]] = LogProb::ln_zero();
                } else {
                    let maybe_initial = if i == observations.len() - 1 {
                        hmm.initial_prob(j)
                    } else {
                        LogProb::ln_one()
                    };
                    vals[[0, *j]] = LogProb::ln_one() + hmm.observation_prob(j, o) + maybe_initial;
                }
            }
        } else {
            // Previous columns.
            for j in hmm.states() {
                if (i >= impossible_range.start)
                    && (i < impossible_range.end)
                    && (j == *impossible_state)
                {
                    vals[[i, *j]] = LogProb::ln_zero();
                } else {
                    let maybe_initial = if i == observations.len() - 1 {
                        hmm.initial_prob(j)
                    } else {
                        LogProb::ln_one()
                    };
                    let xs = hmm.states()
                        .map(|k| {
                            vals[[i - 1, *k]]
                                + hmm.transition_prob_idx(j, k, n - i)
                                + hmm.observation_prob(j, o)
                                + maybe_initial
                        })
                        .collect::<Vec<LogProb>>();
                    vals[[i, *j]] = LogProb::ln_sum_exp(&xs);
                }
            }
        }
    }

    // Compute final probability.
    let prob = LogProb::ln_sum_exp(vals.row(observations.len() - 1).into_slice().unwrap());

    (vals, prob)
}

/// Convert copy state to idx.
fn cs_to_idx(state: CopyState) -> usize {
    match state {
        CopyState::Deletion => 0,
        CopyState::Neutral => 1,
        CopyState::Duplication => 2,
    }
}

/// Compute metrics associated with each segment.
fn compute_seg_metrics(
    logger: &Logger,
    segmentation: &Segmentation,
    ranges: &Vec<Range<usize>>,
    covzs: &Vec<f64>,
    options: &GenotypeOptions,
) -> Result<Vec<CnvGenotypeInfo>> {
    let mut result = Vec::new();

    let pos = ranges.iter().map(|r| r.start).collect::<Vec<usize>>();

    // Construct XHMM model for computing qualities.
    let model = xhmm::ExomeModel::new(
        options.xhmm_cnv_rate,
        options.xhmm_mean_target_dist,
        options.xhmm_mean_target_count,
        options.xhmm_z_score_threshold,
        pos,
    );

    // Compute f/b from XHMM paper.
    let (f, _) = forward(&model, &covzs);
    let (b, _) = backward(&model, &covzs);

    // Compute data likelihood ``Pr(y_{1:E})`.
    let data_likelihood = LogProb::ln_sum_exp(&[
        f[[0, 0]] + b[[0, 0]],
        f[[0, 1]] + b[[0, 1]],
        f[[0, 2]] + b[[0, 2]],
    ]);

    // Probability that copy state is equal to a given value within the given range.
    let pr_cs_equals = |range: Range<usize>, cs: CopyState| -> LogProb {
        let state = State(cs_to_idx(cs));
        let total_trans_prob: LogProb = ((range.start + 1)..range.end)
            .map(|i| model.transition_prob_idx(state, state, i))
            .sum();
        let total_emiss_prob: LogProb = ((range.start + 1)..range.end)
            .map(|i| model.observation_prob(state, &covzs[i]))
            .sum();

        total_trans_prob
            + total_emiss_prob
            + f[[range.start, cs_to_idx(cs)]]
            + b[[range.end - 1, cs_to_idx(cs)]] - data_likelihood
    };

    // Probability that copy state is one of the given values.
    let pr_cs_one_of = |range: Range<usize>,
                        b_restr: &Array2<LogProb>,
                        f_restr: &Array2<LogProb>,
                        cs1: CopyState,
                        cs2: CopyState|
     -> LogProb {
        let t = range.end - 1;
        let s1 = cs_to_idx(cs1);
        let s2 = cs_to_idx(cs2);

        LogProb::ln_sum_exp(&[
            f_restr[[t, s1]] + b_restr[[t, s1]],
            f_restr[[t, s2]] + b_restr[[t, s2]],
        ]) - data_likelihood
    };

    for (i, ref segment) in segmentation.segments.iter().enumerate() {
        debug!(
            logger,
            "Segment {} of {}",
            i + 1,
            segmentation.segments.len()
        );
        let copy_state = segment
            .copy_state
            .expect("Segment must have copy state set!");
        if copy_state == CopyState::Neutral {
            continue;
        }

        let state = State(cs_to_idx(copy_state));
        let start = segment.range.start;
        let end = segment.range.end;
        // println!("=> segment = {:?}", &segment);

        // Compute restricted forward and backward probability (\in {1,2})
        let (b_restr12, _) = restricted_backward(&model, &covzs, &State(2), &segment.range);
        let (f_restr12, _) = restricted_forward(&model, &covzs, &State(2), &segment.range);

        // Compute restricted forward and backward probability (\in {2,3})
        let (b_restr23, _) = restricted_backward(&model, &covzs, &State(0), &segment.range);
        let (f_restr23, _) = restricted_forward(&model, &covzs, &State(0), &segment.range);

        let q_exact_cnv = PHREDProb::from(pr_cs_equals(segment.range.clone(), copy_state));
        let (q_some_cnv, q_no_cnv) = if state == State(0) {
            // State is deletion, "some/no CNV" = "some/no DEL"
            let q_some_del_a = pr_cs_one_of(
                segment.range.clone(),
                &b_restr12,
                &f_restr12,
                CopyState::Deletion,
                CopyState::Neutral,
            );
            let q_some_del_b = pr_cs_equals(segment.range.clone(), CopyState::Neutral);
            // println!(
            //     "q_some_del_a={:?}, q_some_del_b={:?}",
            //     &q_some_del_a, &q_some_del_b
            // );
            let q_no_del = PHREDProb::from(pr_cs_one_of(
                segment.range.clone(),
                &b_restr23,
                &f_restr23,
                CopyState::Neutral,
                CopyState::Duplication,
            ));
            (
                PHREDProb::from(q_some_del_a.ln_sub_exp(q_some_del_b)),
                q_no_del,
            )
        } else {
            // State is duplication, "some/no CNV" = "some/no DUP"
            let q_some_dup_a = pr_cs_one_of(
                segment.range.clone(),
                &b_restr23,
                &f_restr23,
                CopyState::Neutral,
                CopyState::Duplication,
            );
            let q_some_dup_b = pr_cs_equals(segment.range.clone(), CopyState::Neutral);
            // if segment.range.len() < 2 {
            //     println!("f_restr23={:?}", &f_restr23);
            //     println!("b_restr23={:?}", &b_restr23);
            //     println!(
            //         "q_some_dup_a={:?}, q_some_dup_b={:?}",
            //         &q_some_dup_a, &q_some_dup_b
            //     );
            // }
            let q_no_dup = PHREDProb::from(pr_cs_one_of(
                segment.range.clone(),
                &b_restr12,
                &f_restr12,
                CopyState::Deletion,
                CopyState::Neutral,
            ));
            (
                PHREDProb::from(q_some_dup_a.ln_sub_exp(q_some_dup_b)),
                q_no_dup,
            )
        };

        let q_left_cnv_bp = PHREDProb::from(
            f[[segment.range.start - 1, cs_to_idx(CopyState::Neutral)]]
                + b[[segment.range.end, cs_to_idx(copy_state)]] - data_likelihood,
        );
        let q_right_cnv_bp = PHREDProb::from(
            f[[segment.range.end - 1, cs_to_idx(copy_state)]]
                + b[[segment.range.end, cs_to_idx(CopyState::Neutral)]]
                - data_likelihood,
        );

        let pr_diploid = pr_cs_equals(segment.range.clone(), CopyState::Neutral);
        let q_diploid = PHREDProb::from(pr_diploid);
        let q_not_diploid = PHREDProb::from(pr_diploid.ln_one_minus_exp());

        let mean_cov = segmentation.values[segment.range.clone()].mean();
        let mean_z_score = covzs[segment.range.clone()].mean();
        let num_target_regions = segment.range.len();

        let segment = (*segment).clone();
        let quals = XhmmCnvQuals {
            copy_state,
            q_exact_cnv,
            q_some_cnv,
            q_no_cnv,
            q_left_cnv_bp,
            q_right_cnv_bp,
            q_not_diploid,
            q_diploid,
            mean_cov,
            mean_z_score,
            num_target_regions,
        };
        result.push(CnvGenotypeInfo { segment, quals });
    }

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
    reader_calls: Option<&mut bcf::IndexedReader>,
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
            Ok(_) => (),
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
        "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Length of the SV in bp\">",
        // FILTER fields
        "##FILTER=<ID=LowQual,Description=\"Low-quality call or genotype.\">",
        // Generic FORMAT fields
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"The genotype in the sample\">",
        "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype-wise filter, \
         semicolon-separated.\">",
        // Coverage- and quality-related FORMAT fields.
        "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Predicted copy number;
         0=no-call, 1=del, 2=diploid, 3=duplication\">",
        "##FORMAT=<ID=EQ,Number=1,Type=Float,Description=\"Phred-scaled probability for \
         'exact del/dup'\">",
        "##FORMAT=<ID=SQ,Number=1,Type=Float,Description=\"Phred-scaled probability for \
         'some del/dup'\">",
        "##FORMAT=<ID=NQ,Number=1,Type=Float,Description=\"Phred-scaled probability for \
         'no del/dup'\">",
        "##FORMAT=<ID=LQ,Number=1,Type=Float,Description=\"Phred-scaled probability for \
         'left del/dup breakpoint'\">",
        "##FORMAT=<ID=RQ,Number=1,Type=Float,Description=\"Phred-scaled probability for \
         'right del/dup breakpoint'\">",
        "##FORMAT=<ID=NDQ,Number=1,Type=Float,Description=\"Phred-scaled probability for \
         'not diploid'\">",
        "##FORMAT=<ID=DQ,Number=1,Type=Float,Description=\"Phred-scaled probability for \
         'diploid'\">",
        // TODO: Artificially add RD and ORD as in XHMM?
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
        let (segmentation, ranges, covzs) =
            read_seg_and_cov(logger, &mut reader, reader_calls.as_mut(), &options)?;

        trace!(logger, "Computing quality metrics of the segmentation");
        let gt_infos = compute_seg_metrics(logger, &segmentation, &ranges, &covzs, &options)?;

        trace!(logger, "Post-filtering segments");
        let gt_infos = filter_segmentation(&gt_infos, &options)?;

        trace!(logger, "Write out calls BCF file");
        write_result(&gt_infos, &chrom, &ranges, &mut writer)?;
        break;
    }

    info!(logger, "=> OK");
    Ok(())
}
