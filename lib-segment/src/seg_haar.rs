//! Integration of Haar-Seg.

use std::env;

use rust_htslib::bcf::{self, Read};
use rust_segment::{reject_nonaberrant_pvalue, seg_haar};
use separator::Separatable;
use shlex;
use slog::Logger;

use super::errors::*;
use options::*;

use lib_shared::bcf_utils;

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

        trace!(logger, "Load log2-scaled coverage");
        let mut cvs: Vec<f64> = Vec::new();
        let mut cv2s: Vec<f64> = Vec::new();
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
                let cv2 = record
                    .format(b"CV2")
                    .float()
                    .expect("Could not access FORMAT/CV2")[0][0] as f64;

                cvs.push(cv);
                cv2s.push(cv2);
            }
        }

        trace!(logger, "Perform segmentation, yield breakpoints");
        let segmentation = seg_haar(
            &cv2s,
            None,
            None,
            &[0..(cv2s.len())],
            options.haar_seg_fdr,
            options.haar_seg_l_min,
            options.haar_seg_l_max,
        );
        debug!(
            logger,
            "Raw segments: {}",
            segmentation.segments.len().separated_string()
        );

        trace!(logger, "Compute P-values");
        let mut p_values = Vec::new();
        for seg in &segmentation.segments {
            let new_len = p_values.len() + seg.range.len();
            let p_val = if seg.range.len() > 2 {
                seg.p_value_significant_student(options.thresh_p_value) as f32
            } else {
                1.0
            };
            p_values.resize(new_len, p_val);
        }

        trace!(logger, "Reject segments not passing P-value filter");
        let segmentation = reject_nonaberrant_pvalue(&segmentation, &cvs, options.thresh_p_value);
        debug!(
            logger,
            "Segments after selecting aberrant (p): {}",
            segmentation.segments.len().separated_string()
        );

        trace!(logger, "Write out segmentation");
        let mut idx = 0; // current index into ncov
        let mut prev_val: Option<(f32, f32, f32)> = None; // (val, log2(val), p_value)
        reader
            .fetch(rid, 0, *chrom_len as u32)
            .chain_err(|| format!("Could not fetch chromosome {}", chrom))?;
        while reader.read(&mut record).is_ok() {
            writer.translate(&mut record);
            let (val, val_log2, p_value) = if !skip_record(&mut record) {
                record.push_filter(writer.header().name_to_id(b"SKIPPED_SEG").unwrap());
                if let Some(prev_val) = prev_val {
                    prev_val
                } else {
                    (1.0, 0.0, 1.0)
                }
            } else {
                let val = if segmentation.values[idx] > 0.0
                    && segmentation.values[idx] < PSEUDO_EPSILON
                {
                    (1.0, 0.0, 1.0)
                } else {
                    (
                        (segmentation.values[idx] - PSEUDO_EPSILON) as f32,
                        segmentation.values_log2[idx] as f32,
                        p_values[idx],
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
            writer.write(&record).expect("Writing the record failed!");
        }
    }

    info!(logger, "=> OK");
    Ok(())
}
