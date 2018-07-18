//! Segmentation using the XHMM approach.

use std::env;

use bio::stats::hmm::{viterbi, State};
use rust_htslib::bcf::{self, Read};
// use rust_segment::{reject_nonaberrant_pvalue, seg_haar};
// use separator::Separatable;
use shlex;
use slog::Logger;

use super::errors::*;
use options::*;

use lib_shared::bcf_utils;
use rust_segment::shared::{CopyState, Segment, Segmentation};

/// Code for implementing the XHMM-based segmentation.
mod xhmm {
    use bio::stats::hmm::{Model, State, StateIter, StateTransitionIter};
    use bio::stats::{LogProb, Prob};
    use statrs::distribution::{Continuous, Normal};

    /// Integer used for marking "deleted" state.
    pub const DEL: usize = 0;
    /// Integer used for marking "neutral" state.
    pub const NEUTRAL: usize = 1;
    //// Integer used for marking "duplicated" state.
    pub const DUP: usize = 2;

    /// Implementation of XHMM variation of the Hidden Markov Model.
    #[derive(Debug)]
    pub struct ExomeModel {
        // Parameters of the Exome Hidden Markov Model.
        /// Expected CNV rate.
        cnv_rate: f64,
        /// Expected mean target distance inside a CNV.
        mean_target_dist: f64,
        /// Expected mean target count per CNV.
        mean_target_count: f64,
        /// Z-score threshold for emission probability.
        z_score_threshold: f64,

        /// Positions of targets on the genome.
        obs_pos: Vec<usize>,

        /// The distributions for the states `DEL`, `NEUTRAL`, `DUP`.
        dists: Vec<Normal>,
    }

    impl ExomeModel {
        pub fn new(
            cnv_rate: f64,
            mean_target_dist: f64,
            mean_target_count: f64,
            z_score_threshold: f64,
            obs_pos: Vec<usize>,
        ) -> Self {
            Self {
                cnv_rate,
                mean_target_dist,
                mean_target_count,
                obs_pos,
                z_score_threshold,
                dists: vec![
                    Normal::new(-z_score_threshold, 1.0).unwrap(),
                    Normal::new(0.0, 1.0).unwrap(),
                    Normal::new(z_score_threshold, 1.0).unwrap(),
                ],
            }
        }

        /// Return distance of target at `to_idx` and the previous one.
        fn distance(&self, to_idx: usize) -> f64 {
            (self.obs_pos[to_idx] - self.obs_pos[to_idx - 1]) as f64
        }
    }

    impl Model<f64> for ExomeModel {
        fn num_states(&self) -> usize {
            3
        }

        /// Return iterator over the states of an HMM.
        fn states(&self) -> StateIter {
            StateIter::new(self.num_states())
        }

        fn transitions(&self) -> StateTransitionIter {
            StateTransitionIter::new(self.num_states())
        }

        fn initial_prob(&self, state: State) -> LogProb {
            match state {
                State(NEUTRAL) => LogProb::from(Prob::from(1.0 - 2.0 * self.cnv_rate)),
                _ => LogProb::from(Prob::from(self.cnv_rate)),
            }
        }

        fn transition_prob(&self, _from: State, _to: State) -> LogProb {
            panic!("Should not be called, expecte calls to transition_prob_idx instead");
        }

        fn transition_prob_idx(&self, from: State, to: State, to_idx: usize) -> LogProb {
            let f = (-self.distance(to_idx) / self.mean_target_dist).exp();
            let p = self.cnv_rate;
            let q = 1.0 / self.mean_target_count;

            let res = match (from, to) {
                (State(DEL), State(DEL)) => f * (1.0 - q) + (1.0 - f) * p,
                (State(DEL), State(NEUTRAL)) => f * q + (1.0 - f) * (1.0 - 2.0 * p),
                (State(DEL), State(DUP)) => (1.0 - f) * p,
                (State(NEUTRAL), State(DEL)) => p,
                (State(NEUTRAL), State(NEUTRAL)) => 1.0 - 2.0 * p,
                (State(NEUTRAL), State(DUP)) => p,
                (State(DUP), State(DEL)) => (1.0 - f) * p,
                (State(DUP), State(NEUTRAL)) => f * q + (1.0 - f) * (1.0 - 2.0 * p),
                (State(DUP), State(DUP)) => f * (1.0 - q) + (1.0 - f) * p,
                _ => panic!("Invalid state combination"),
            };
            LogProb::from(Prob::from(res))
        }

        fn observation_prob(&self, state: State, observation: &f64) -> LogProb {
            LogProb::from(Prob::from(self.dists[*state].pdf(*observation)))
        }
    }

}

/// Perfor XHMM-based segmentation.
fn xhmm_seg(
    logger: &mut Logger,
    pos: Vec<usize>,
    cvs: &Vec<f64>,
    cv2s: &Vec<f64>,
    cvzs: &Vec<f64>,
    options: &SegmentOptions,
) -> Segmentation {
    info!(logger, "Constructing XHMM...");
    let model = xhmm::ExomeModel::new(
        options.xhmm_cnv_rate,
        options.xhmm_mean_target_dist,
        options.xhmm_mean_target_count,
        options.xhmm_z_score_threshold,
        pos,
    );

    info!(logger, "Computing Viterbi path...");
    let viterbi_path = viterbi(&model, cvzs).0;

    let mut segments: Vec<Segment> = Vec::new();
    let mut cn_states: Vec<CopyState> = Vec::new();
    let mut start: usize = 0;
    for nxt in 1..viterbi_path.len() {
        let cur = nxt - 1;
        if viterbi_path[cur] != viterbi_path[nxt] {
            segments.push(Segment {
                range: start..nxt,
                copy_state: match viterbi_path[cur] {
                    State(xhmm::DEL) => Some(CopyState::Deletion),
                    State(xhmm::NEUTRAL) => Some(CopyState::Neutral),
                    State(xhmm::DUP) => Some(CopyState::Duplication),
                    _ => panic!("Unknown state"),
                },
                mean: 0.0,
                std_dev: 0.0,
                mean_log2: 0.0,
                std_dev_log2: 0.0,
            });
            start = nxt;
        }

        cn_states.push(match viterbi_path[cur] {
            State(xhmm::DEL) => CopyState::Deletion,
            State(xhmm::NEUTRAL) => CopyState::Neutral,
            State(xhmm::DUP) => CopyState::Duplication,
            _ => panic!("Invalid state!"),
        });
    }
    if start < viterbi_path.len() {
        segments.push(Segment {
            range: start..(viterbi_path.len()),
            copy_state: match viterbi_path[start] {
                State(xhmm::DEL) => Some(CopyState::Deletion),
                State(xhmm::NEUTRAL) => Some(CopyState::Neutral),
                State(xhmm::DUP) => Some(CopyState::Duplication),
                _ => panic!("Unknown state"),
            },
            mean: 0.0,
            std_dev: 0.0,
            mean_log2: 0.0,
            std_dev_log2: 0.0,
        });
        cn_states.push(match viterbi_path[start] {
            State(xhmm::DEL) => CopyState::Deletion,
            State(xhmm::NEUTRAL) => CopyState::Neutral,
            State(xhmm::DUP) => CopyState::Duplication,
            _ => panic!("Invalid state!"),
        });
    }
    for ref mut seg in &mut segments {
        seg.update_stats(&cvs);
    }

    info!(logger, " => Done with Viterbi");

    let mut cvs = Vec::new();
    let mut cv2s = Vec::new();
    for ref seg in &segments {
        let old_len = cvs.len();
        cvs.resize(old_len + seg.range.len(), seg.mean);
        cv2s.resize(old_len + seg.range.len(), seg.mean);
    }

    Segmentation {
        segments,
        values: cvs,
        values_log2: cv2s,
        cn_states: Some(cn_states),
    }
}

// TODO: remove rundant code, cmp. seg_haar.rs

/// Epsilon to add normalized metric to prevent -nan through log2()
const PSEUDO_EPSILON: f64 = 1e-12;

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
    let is_gap = record.info(b"GAP").flag().unwrap_or(false);

    !cv.is_finite() || !cv2.is_finite() || is_gap
}

/// Perform segmentation using the Haar-Seg algorithm.
pub fn run_segmentation(logger: &mut Logger, options: &SegmentOptions) -> Result<()> {
    info!(logger, "Computing segmentation using XHMM");

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

        trace!(logger, "Load coverage Z scores");
        let mut cvzs: Vec<f64> = Vec::new();
        let mut cvs: Vec<f64> = Vec::new();
        let mut cv2s: Vec<f64> = Vec::new();
        let mut pos: Vec<usize> = Vec::new();
        let mut record = reader.empty_record();
        loop {
            match reader.read(&mut record) {
                Ok(_) => (),
                Err(bcf::ReadError::NoMoreRecord) => break,
                _ => bail!("Error reading BCF record"),
            }

            if !skip_record(&mut record) {
                let end = record
                    .info(b"END")
                    .integer()
                    .expect("Could not access FORMAT/END")
                    .unwrap()[0] as usize;
                let center = (record.pos() as usize + end) / 2;

                let cv = record
                    .format(b"CV")
                    .float()
                    .expect("Could not access FORMAT/CV")[0][0] as f64;
                let cvz = record
                    .format(b"CVZ")
                    .float()
                    .expect("Could not access FORMAT/CVZ")[0][0] as f64;
                let cv2 = record
                    .format(b"CV2")
                    .float()
                    .expect("Could not access FORMAT/CV2")[0][0] as f64;

                // Limit outliers. TODO: think of something smarter.
                const Z_SCORE_LIMIT_FACTOR: f64 = 5.0;
                let z_score_limit = options.xhmm_z_score_threshold * Z_SCORE_LIMIT_FACTOR;
                let cvz = if cvz < -z_score_limit {
                    -z_score_limit
                } else if cvz > z_score_limit {
                    z_score_limit
                } else {
                    cvz
                };

                pos.push(center);
                cvs.push(cv);
                cvzs.push(cvz);
                cv2s.push(cv2);
            }
        }

        info!(logger, "Computing segmentation");
        let segmentation = xhmm_seg(logger, pos, &cvs, &cv2s, &cvzs, options);

        trace!(logger, "Write out segmentation");
        let mut idx = 0; // current index into ncov
        let mut prev_val: Option<(f32, f32, f32, i32)> = None; // (val, log2(val), p_value, cn_state)
        reader
            .fetch(rid, 0, *chrom_len as u32)
            .chain_err(|| format!("Could not fetch chromosome {}", chrom))?;
        assert_eq!(
            segmentation.values.len(),
            segmentation.cn_states.as_ref().unwrap().len()
        );
        while reader.read(&mut record).is_ok() {
            writer.translate(&mut record);
            let (val, val_log2, p_value, cn_state) = if skip_record(&mut record) {
                record.push_filter(writer.header().name_to_id(b"SKIPPED_SEG").unwrap());
                if let Some(prev_val) = prev_val {
                    prev_val
                } else {
                    (1.0, 0.0, 1.0, 1)
                }
            } else {
                let val = if segmentation.values[idx] > 0.0
                    && segmentation.values[idx] < PSEUDO_EPSILON
                {
                    (1.0, 0.0, 1.0, 1)
                } else {
                    (
                        (segmentation.values[idx] - PSEUDO_EPSILON) as f32,
                        segmentation.values_log2[idx] as f32,
                        1.0, //p_values[idx],  // TODO: compute empirical p value
                        match segmentation.cn_states.as_ref().unwrap()[idx] {
                            CopyState::Deletion => 0,
                            CopyState::Neutral => 1,
                            CopyState::Duplication => 2,
                        },
                    )
                };
                prev_val = Some(val);
                idx += 1;
                val
            };
            record
                .push_format_float(b"SGP", &[p_value])
                .chain_err(|| "Could not write FORMAT/SGP")?;
            record
                .push_format_float(b"SG", &[val])
                .chain_err(|| "Could not write FORMAT/SG")?;
            record
                .push_format_float(b"SG2", &[val_log2])
                .chain_err(|| "Could not write FORMAT/SG2")?;
            record
                .push_format_integer(b"SGS", &[cn_state])
                .chain_err(|| "Could not write FORMAT/SGS")?;
            writer.write(&record).expect("Writing the record failed!");
        }
    }

    info!(logger, "=> OK");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::xhmm::*;
    use bio::stats::hmm::{viterbi, State};
    use bio::stats::{LogProb, Prob};

    #[test]
    fn test_xhmm_simple() {
        let obs: Vec<f64> = vec![
            -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 3.0, 3.0, 3.0,
            3.0, 3.0,
        ];
        let obs_pos: Vec<usize> = vec![
            1_100_000, 1_110_000, 1_120_000, 1_130_000, 1_140_000, 1_150_000, 2_100_000, 2_110_000,
            2_120_000, 2_130_000, 2_140_000, 2_150_000, 3_100_000, 3_110_000, 3_120_000, 3_130_000,
            3_140_000, 3_150_000,
        ];
        let model: ExomeModel = ExomeModel::new(1e-06, 50_000.0, 6.0, 2.0, obs_pos);
        let (path, log_prob) = viterbi(&model, &obs);
        let prob = Prob::from(log_prob);

        let expected = vec![0, 1, 2]
            .iter()
            .map(|i| State(*i))
            .collect::<Vec<State>>();
        assert_eq!(expected, path);
        assert_eq!(0.0_f64, *prob);
    }
}
