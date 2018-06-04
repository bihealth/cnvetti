include!(concat!(env!("OUT_DIR"), "/version.rs"));

pub mod options;
pub use self::options::*;

use std::env;
use std::mem;
use std::ops::Range;

use shlex;

use slog::Logger;

use cli::shared::{self, process};
use separator::Separatable;

use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::{self, Read as BcfRead};

use statrs::distribution::{StudentsT, Univariate};

/// Collect information about one segment.
#[derive(Debug, Clone)]
struct SegmentInfo {
    segment: Range<usize>,
    length: usize,
    cov_mean: f32,
    // Error square sum
    cov_err: f64,
}

impl SegmentInfo {
    /// Compute standard deviation from length and `cov_err` error sums.
    fn cov_sd(&self) -> f64 {
        ((self.cov_err as f64) / ((self.length - 1) as f64)).sqrt()
    }

    /// Merge `other` into `self` and return merge result.
    fn merge(&self, other: &Self) -> Self {
        let a = (self.length as f32) / ((self.length + other.length) as f32);
        let b = 1.0 - a;

        SegmentInfo {
            segment: Range {
                start: self.segment.start,
                end: other.segment.end,
            },
            length: self.length + other.length,
            cov_mean: self.cov_mean * a + other.cov_mean * b,
            cov_err: self.cov_err + other.cov_err,
        }
    }

    /// Compute whether merging should be performed with given `threshold` on P-value.
    ///
    /// Uses Welch's test.
    fn can_merge(&self, other: &Self, threshold: f64) -> bool {
        use statrs::distribution::StudentsT;

        let x_1 = self.cov_mean as f64;
        let x_2 = other.cov_mean as f64;
        let s_1 = self.cov_err.sqrt();
        let s_2 = other.cov_err.sqrt();
        let n_1 = self.length as f64;
        let n_2 = other.length as f64;

        // Helper for computing `x^2`.
        fn sq(x: f64) -> f64 {
            x * x
        }
        /// Helper for computing `x^4`.
        fn q(x: f64) -> f64 {
            x * x * x * x
        }

        // Compute test statistic and degrees of freedom for t Distribution.
        let t = (x_1 - x_2) / (sq(s_1) / n_1 + sq(s_2) / n_2).sqrt();
        let nu = sq(sq(s_1) / n_1 + sq(s_2) / n_2)
            / ((q(s_1) / sq(n_1) / (n_1 - 1.0)) + q(s_2) / sq(n_2) / (n_2 - 1.0));

        let dist = StudentsT::new(0.0, 1.0, nu as f64).unwrap();
        (1.0 - dist.cdf(t)) < threshold
    }

    /// Compute P-value of segment.
    fn p_value(&self, options: &Options) -> f64 {
        // TODO: robust handling of X/Y chromosome
        // TODO: is this correct?
        // Note that coverages are already normalized.

        // Compute test statistics and obtain distribution.
        let t = (1.0 - self.cov_mean as f64) / self.cov_sd() * (self.length as f64).sqrt();
        let dist = StudentsT::new(0.0, 1.0, (self.length - 1) as f64)
            .expect("Could not find t distribution");
        // Compute p value.
        let p = 1.0 - dist.cdf(t);
        // println!("Uncorrected P-value: {}", p);
        // Return P-value corrected for multiple testing.
        p * options.significant_p_val_thresh / (self.length as f64)
    }

    /// Return `SegmentInfo` with mean set to copy number neutral state.
    fn with_reset_mean(&self) -> Self {
        SegmentInfo {
            cov_mean: 1.0, // TODO: properly reset on X and Y chrom
            ..self.clone()
        }
    }
}

/// Collect information about all segments from `reader`.
fn collect_segments(reader: &mut bcf::IndexedReader) -> Vec<SegmentInfo> {
    let mut result: Vec<SegmentInfo> = Vec::new();

    // Temporary segment.
    #[derive(Debug, Clone)]
    struct TempSegmentInfo {
        // This segment's range.
        segment: Range<usize>,
        // Number of skipped reference so far.
        skipped_ref: usize,
        // Current mean (all the same).
        cov_mean: f32,
        // Current sum of variance squares.
        cov_err: f64,
    };
    impl TempSegmentInfo {
        fn length(&self) -> usize {
            self.segment.len() - self.skipped_ref
        }
        fn cov_sd(&self) -> f64 {
            ((self.cov_err as f64) / ((self.length() - 1) as f64)).sqrt()
        }
    }
    let mut tmp: Option<TempSegmentInfo> = None;

    let mut record = reader.empty_record();
    let mut prev_rid = None;
    while reader.read(&mut record).is_ok() {
        if let Some(prev_rid) = prev_rid {
            if record.rid() != prev_rid {
                break; // contig change!
            }
        }

        // Extract the interesting information about the segment from BCF record.
        let start = record.pos() as usize;
        let end = record
            .info(b"END")
            .integer()
            .expect("Could not read FORMAT/END")
            .expect("Could not access FORMAT/END")[0] as usize;
        // TODO: we must use SCOV here, should we scale NCOVSD?
        let cov_mean = record
            .format(b"SCOV")
            .float()
            .expect("Could not access FORMAT/SCOV")[0][0];
        let cov_sd = record
            .format(b"NCOVSD")
            .float()
            .expect("Could not access FORMAT/NCOVSD")[0][0];

        // Possibly ignore this window.
        if cov_mean.is_missing()
            || !cov_mean.is_finite()
            || cov_sd.is_missing()
            || !cov_sd.is_finite()
        {
            continue;
        }
        let cov_sd = cov_sd as f64;
        // println!("{} {} {}", cov_sd, cov_sd * cov_sd, cov_sd * cov_sd * ((end - start) as f64));
        // println!("{:?}", tmp);

        // Update temporary segment, possibly saving starting new one and pushing current.
        tmp = Some(match tmp {
            Some(tmp) => {
                println!("tmp = {:?}, sd = {}", tmp, tmp.cov_sd());
                if cov_mean == tmp.cov_mean {
                    // Extend current segment.
                    TempSegmentInfo {
                        segment: Range {
                            start: tmp.segment.start,
                            end: end,
                        },
                        skipped_ref: tmp.skipped_ref + (start - tmp.segment.end),
                        cov_mean: tmp.cov_mean,
                        cov_err: tmp.cov_err + cov_sd * cov_sd * ((end - start) as f64),
                    }
                } else {
                    // Save current segment and start new.
                    let length = tmp.segment.len() - tmp.skipped_ref;
                    result.push(SegmentInfo {
                        segment: tmp.segment.clone(),
                        length: length,
                        cov_mean: tmp.cov_mean,
                        cov_err: tmp.cov_err,
                    });
                    TempSegmentInfo {
                        segment: Range { start, end },
                        skipped_ref: 0,
                        cov_mean: cov_mean,
                        cov_err: tmp.cov_err,
                    }
                }
            }
            None => {
                // Start first segment.
                TempSegmentInfo {
                    segment: Range { start, end },
                    skipped_ref: 0,
                    cov_mean: cov_mean,
                    cov_err: cov_sd * cov_sd * ((end - start) as f64),
                }
            }
        });
        // println!("tmp = {:?}", tmp);

        prev_rid = Some(record.rid());
    }

    if let Some(tmp) = tmp {
        // Save last segment, if any.
        let length = tmp.segment.len() - tmp.skipped_ref;
        result.push(SegmentInfo {
            segment: tmp.segment.clone(),
            length: length,
            cov_mean: tmp.cov_mean,
            cov_err: tmp.cov_err,
        });
    }

    result
}

/// Set segment means to copy number neutral state if not significantly deviated...
///
/// Merging sements is handled in `merge_segments`.
fn validate_segments(segments: &[SegmentInfo], options: &Options) -> Vec<SegmentInfo> {
    segments
        .iter()
        .map(|seg| {
            // println!("seg = {:?}", seg);
            if seg.p_value(options) < options.significant_p_val_thresh {
                seg.with_reset_mean()
            } else {
                seg.clone()
            }
        })
        .collect::<Vec<SegmentInfo>>()
}

/// Merge segments from lef to right using the strategy also used in CNVnator.
fn merge_segments(segments: &[SegmentInfo], options: &Options) -> Vec<SegmentInfo> {
    if segments.len() <= 2 {
        return Vec::from(segments);
    }

    let mut result: Vec<SegmentInfo> = Vec::new();
    result.push(segments[0].clone());

    for i in 1..segments.len() {
        let prev = result.last().unwrap().clone();
        let current = segments[i].clone();

        if current.cov_mean == prev.cov_mean || prev.can_merge(&current, options.merge_p_val_thresh)
        {
            result.pop();
            result.push(prev.merge(&current));
        } else {
            result.push(current);
        }
    }

    result
}

/// Process one region.
fn process_region(
    reader: &mut bcf::IndexedReader,
    writer: &mut bcf::Writer,
    logger: &mut Logger,
    options: &Options,
    (chrom, start, end): &(String, u32, u32),
) {
    info!(
        logger,
        "Processing {}:{}-{}",
        chrom,
        (start + 1).separated_string(),
        end.separated_string()
    );

    // Resolve contig name to contig index.
    let rid = reader
        .header()
        .name2rid(chrom.as_bytes())
        .expect("Could not resolve contig name to index");

    // Go over the region, collecting segment mean coverage, as well as standard deviations.
    reader
        .fetch(rid, *start, *end)
        .expect("Could not fetch region from BCF file");
    println!("Collecting...");
    let segment_infos = collect_segments(reader);
    for j in 0..10 {
        println!("{} -- {:?}", j, segment_infos[j]);
    }

    // Perform t-test on segments and reset their mean to neutral copy number state if the test
    // fails (p value threshold taken from `options`).
    println!("Validating...");
    let mut segment_infos = validate_segments(&segment_infos, options);
    for j in 0..10 {
        println!("{} -- {:?}", j, segment_infos[j]);
    }
    info!(
        logger,
        "Collected {} segments",
        segment_infos.len().separated_string()
    );

    // Merge segments until the number of segments does not fall any more, up the number of
    // times as configured in `options`.
    println!("Merging...");
    for _i in 0..options.max_merges {
        let mut merged_infos = merge_segments(&segment_infos, options);
        if merged_infos.len() != segment_infos.len() {
            mem::swap(&mut merged_infos, &mut segment_infos);
            for j in 0..10 {
                if j < segment_infos.len() {
                    println!("{} -- {:?}", j, segment_infos[j]);
                }
            }
        } else {
            break;
        }
    }
    info!(
        logger,
        "Merging => {} segments",
        segment_infos.len().separated_string()
    );

    // TODO: "knit" long CNV calls with short non-CNV segments

    // TODO: also do last step in CNVnator?

    // Finally, construct CNV call records and write to output.
    for info in segment_infos.iter().filter(|r| r.cov_mean != 1.0) {
        let mut record = writer.empty_record();

        let alleles_v = vec![
            Vec::from("N"),
            if info.cov_mean < 1.0 {
                Vec::from("<DEL>")
            } else {
                Vec::from("<DUP>")
            },
        ];
        let alleles = alleles_v
            .iter()
            .map(|x| x.as_slice())
            .collect::<Vec<&[u8]>>();

        // REF..QUAL
        record.inner_mut().rid = rid as i32;
        record.set_pos(info.segment.start as i32);
        record
            .set_id(format!("CNV_{}_{}_{}", chrom, info.segment.start, info.segment.end).as_bytes())
            .expect("Cannot set ID");
        record
            .set_alleles(&alleles)
            .map_err(|e| format!("Could not update alleles: {}", e))
            .unwrap();

        // INFO/*
        record
            .push_info_integer(b"END", &[info.segment.end as i32])
            .expect("Cannot set INFO/END");
        record
            .push_info_integer(b"CNVLEN", &[info.segment.len() as i32])
            .expect("Cannot set INFO/CNVLEN");

        // FORMAT/GT
        record
            .push_format_integer(b"GT", &[0, 1])
            .map_err(|e| format!("Could not write FORMAT/GT: {}", e))
            .unwrap();
        record
            .push_format_integer(b"CN", &[(info.cov_mean * 2.0).round() as i32])
            .map_err(|e| format!("Could not write FORMAT/CN: {}", e))
            .unwrap();
        record
            .push_format_float(b"CNR", &[(info.cov_mean * 2.0) as f32])
            .map_err(|e| format!("Could not write FORMAT/CNR: {}", e))
            .unwrap();

        writer.write(&record).expect("Writing to BCF file failed!");
    }
}

/// Main entry point for the "segment" command.
pub fn call(logger: &mut Logger, options: &Options) -> Result<(), String> {
    let options = options.clone();

    // let mut app = CoverageApp::new(logger, options);
    info!(logger, "Running cnvetti segment");
    info!(logger, "Configuration: {:?}", &options);

    debug!(logger, "Opening input file");
    let mut reader = bcf::IndexedReader::from_path(options.input.clone())
        .expect("Could not open input BCF file");
    if options.io_threads > 0 {
        reader
            .set_threads(options.io_threads as usize)
            .expect("Could not set I/O thread count");
    }

    // Open the output file in its own block so we can close before creating the index.
    {
        let mut writer = {
            // Construct extended header.
            let mut header = bcf::Header::from_template(reader.header());
            // TODO: remove header fields...
            let lines = vec![
                "##ALT=<ID=DEL,Description=\"Deletion relative to the reference\">",
                "##ALT=<ID=DUP,Description=\"Region of elevated copy number relative to the \
                 reference\">",
                "##INFO=<ID=CNVLEN,Number=1,Type=Integer,Description=\"Length of CNV\">",
                "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
                "##FORMAT=<ID=GP,Number=1,Type=Integer,Description=\"Genotype quality\">",
                "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for \
                 imprecise events\">",
                "##FORMAT=<ID=CNR,Number=1,Type=Integer,Description=\"Real-valued copy number \
                 genotype for imprecise events\">",
                "##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=\"Copy number genotype \
                 quality for imprecise events\">",
                "##FORMAT=<ID=CNL,Number=G,Type=Float,Description=\"Copy number genotype \
                 likelihood for imprecise events\">",
            ];
            for line in lines {
                header.push_record(line.as_bytes());
            }
            header.push_record(format!("##cnvetti_callVersion={}", VERSION).as_bytes());
            header.push_record(
                format!(
                    "##cnvetti_callCommand={}",
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
                .expect("Could not open output BCF file")
        };
        if options.io_threads > 0 {
            writer
                .set_threads(options.io_threads as usize)
                .expect("Could not set I/O thread count");
        }

        info!(logger, "Processing...");
        process(&mut reader, &mut writer, logger, &options, &process_region);
    }

    // Build index on the output file.
    shared::build_index(logger, &options.output);

    info!(logger, "Done. Have a nice day!");

    Ok(())
}
