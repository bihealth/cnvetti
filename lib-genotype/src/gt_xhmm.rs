//! Genotyping using the XHMM approach.

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

/// Perform genotyping using the XHMM algorithm.
pub fn run_genotyping(logger: &mut Logger, options: &GenotypeOptions) -> Result<()> {
    info!(logger, "Computing genotyping using XHMM");

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
        header.push_record(format!("##cnvetti_genotypeVersion={}", "0.1.0").as_bytes());
        header.push_record(
            format!(
                "##cnvetti_genotypeCommand={}",
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

        // trace!(logger, "Jumping in coverage BCF file");
        // let rid = reader
        //     .header()
        //     .name2rid(chrom.as_bytes())
        //     .chain_err(|| format!("Could not translate header name {}", chrom))?;
        // reader
        //     .fetch(rid, 0, *chrom_len as u32)
        //     .chain_err(|| format!("Could not fetch chromosome {}", chrom))?;

        // trace!(logger, "Load coverage Z scores");
        // let mut cvzs: Vec<f64> = Vec::new();
        // let mut cvs: Vec<f64> = Vec::new();
        // let mut cv2s: Vec<f64> = Vec::new();
        // let mut pos: Vec<usize> = Vec::new();
        // let mut record = reader.empty_record();
        // loop {
        //     match reader.read(&mut record) {
        //         Ok(_) => (),
        //         Err(bcf::ReadError::NoMoreRecord) => break,
        //         _ => bail!("Error reading BCF record"),
        //     }

        //     if !skip_record(&mut record) {
        //         let end = record
        //             .info(b"END")
        //             .integer()
        //             .expect("Could not access FORMAT/END")
        //             .unwrap()[0] as usize;
        //         let center = (record.pos() as usize + end) / 2;

        //         let cv = record
        //             .format(b"CV")
        //             .float()
        //             .expect("Could not access FORMAT/CV")[0][0] as f64;
        //         let cvz = record
        //             .format(b"CVZ")
        //             .float()
        //             .expect("Could not access FORMAT/CVZ")[0][0] as f64;
        //         let cv2 = record
        //             .format(b"CV2")
        //             .float()
        //             .expect("Could not access FORMAT/CV2")[0][0] as f64;

        //         // Limit outliers. TODO: think of something smarter.
        //         const Z_SCORE_LIMIT_FACTOR: f64 = 5.0;
        //         let z_score_limit = options.xhmm_z_score_threshold * Z_SCORE_LIMIT_FACTOR;
        //         let cvz = if cvz < -z_score_limit {
        //             -z_score_limit
        //         } else if cvz > z_score_limit {
        //             z_score_limit
        //         } else {
        //             cvz
        //         };

        //         pos.push(center);
        //         cvs.push(cv);
        //         cvzs.push(cvz);
        //         cv2s.push(cv2);
        //     }
        // }

        info!(logger, "Computing genotype");
        // let segmentation = xhmm_seg(logger, pos, &cvs, &cv2s, &cvzs, options);

        trace!(logger, "Write out genotypes");
        // let mut idx = 0; // current index into ncov
        // let mut prev_val: Option<(f32, f32, f32, i32)> = None; // (val, log2(val), p_value, cn_state)
        // reader
        //     .fetch(rid, 0, *chrom_len as u32)
        //     .chain_err(|| format!("Could not fetch chromosome {}", chrom))?;
        // assert_eq!(
        //     segmentation.values.len(),
        //     segmentation.cn_states.as_ref().unwrap().len()
        // );
        // while reader.read(&mut record).is_ok() {
        //     writer.translate(&mut record);
        //     let (val, val_log2, p_value, cn_state) = if skip_record(&mut record) {
        //         record.push_filter(writer.header().name_to_id(b"SKIPPED_SEG").unwrap());
        //         if let Some(prev_val) = prev_val {
        //             prev_val
        //         } else {
        //             (1.0, 0.0, 1.0, 1)
        //         }
        //     } else {
        //         let val = if segmentation.values[idx] > 0.0
        //             && segmentation.values[idx] < PSEUDO_EPSILON
        //         {
        //             (1.0, 0.0, 1.0, 1)
        //         } else {
        //             (
        //                 (segmentation.values[idx] - PSEUDO_EPSILON) as f32,
        //                 segmentation.values_log2[idx] as f32,
        //                 1.0, //p_values[idx],  // TODO: compute empirical p value
        //                 match segmentation.cn_states.as_ref().unwrap()[idx] {
        //                     CopyState::Deletion => 0,
        //                     CopyState::Neutral => 1,
        //                     CopyState::Duplication => 2,
        //                 },
        //             )
        //         };
        //         prev_val = Some(val);
        //         idx += 1;
        //         val
        //     };
        //     record
        //         .push_format_float(b"SGP", &[p_value])
        //         .chain_err(|| "Could not write FORMAT/SGP")?;
        //     record
        //         .push_format_float(b"SG", &[val])
        //         .chain_err(|| "Could not write FORMAT/SG")?;
        //     record
        //         .push_format_float(b"SG2", &[val_log2])
        //         .chain_err(|| "Could not write FORMAT/SG2")?;
        //     record
        //         .push_format_integer(b"SGS", &[cn_state])
        //         .chain_err(|| "Could not write FORMAT/SGS")?;
        //     writer.write(&record).expect("Writing the record failed!");
        // }
    }

    info!(logger, "=> OK");
    Ok(())
}
