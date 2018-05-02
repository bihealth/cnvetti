/// Re-usable code for normalization.

include!(concat!(env!("OUT_DIR"), "/version.rs"));

use regex::Regex;

use shlex;

use slog::Logger;

use std::env;
use std::str;

use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::{self, Read as BcfRead};

use cli::normalize::options::Options;

/// Epsilon to use to prevent log2(0)
const EPSILON: f32 = 1e-7;

/// Information needed for normalization.
#[derive(Debug, Clone, Copy)]
pub struct NormData {
    pub use_chrom: bool,
    pub is_gap: bool,
    pub is_blacklisted: bool,
    pub gc_content: f64,
    pub mapability: f64,
    pub coverage: f64,
}

impl NormData {
    /// Construct new, empty `NormData`.
    pub fn new(
        use_chrom: bool,
        is_gap: bool,
        is_blacklisted: bool,
        gc_content: f64,
        mapability: f64,
        coverage: f64,
    ) -> NormData {
        NormData {
            use_chrom,
            is_gap,
            is_blacklisted,
            gc_content,
            mapability,
            coverage,
        }
    }

    /// Whether to use this record, also based on minimal mapability.
    pub fn use_record(&self, min_mapability: f64) -> bool {
        self.use_chrom && !self.is_gap && !self.is_blacklisted
            && (self.mapability >= min_mapability)
    }

    /// Return GC bin ID.
    pub fn gc_bin(&self, gc_step: f64) -> usize {
        if self.gc_content.is_finite() {
            (self.gc_content / gc_step).round() as usize
        } else {
            0_usize
        }
    }
}

/// Load `NormData` given the specific configuration in `options`.
pub fn load_norm_data(logger: &Logger, options: &Options) -> Result<Vec<NormData>, String> {
    let mut res: Vec<NormData> = Vec::new();

    // Regular expression for matching contigs.
    let re = Regex::new(&options.contig_regex).unwrap();

    {
        debug!(
            logger,
            "Reading in BCF for GAP/GC/MAP/COV from {}", options.input
        );
        let mut reader = bcf::Reader::from_path(&options.input)
            .map_err(|e| format!("Could not open BCF file: {}", e))?;
        if options.io_threads > 0 {
            reader
                .set_threads(options.io_threads as usize)
                .map_err(|e| format!("Could not set threads to BCF reader: {}", e))?;
        }

        // The BCF record to use for the buffer.
        let mut record = reader.empty_record();
        // Previous reference ID, for printing messages to user.
        let mut prev_rid = Option::None;
        // We need a copy of the lightweight header.
        let header = reader.header().clone();

        while reader.read(&mut record).is_ok() {
            if !prev_rid.is_some() || prev_rid.unwrap() != record.rid() {
                info!(
                    logger,
                    "Starting on contig {}",
                    str::from_utf8(header.rid2name(record.rid().unwrap())).unwrap()
                );
            }

            // Get "is gap" flag.
            let is_gap = {
                record
                    .info(b"GAP")
                    .integer()
                    .map_err(|e| format!("Could not load gaps: {}", e))?
                    .ok_or("INFO/GAP was empty")?[0] != 0
            };
            // Get "is blacklisted" flag.
            let is_blacklisted = {
                record
                    .info(b"BLACKLIST")
                    .integer()
                    .unwrap_or(Some(&[0]))
                    .ok_or("INFO/BLACKLIST was empty")?[0] != 0
            };
            // Get INFO/GC.
            let gc_content = {
                record
                    .info(b"GC")
                    .float()
                    .map_err(|e| format!("Could not read INFO/GC: {}", e))
                    .map_err(|e| format!("INFO/GC was empty: {}", e))?
                    .ok_or("Could not access INFO/GC")?[0]
            };
            // TODO: handle missing value in a better way?
            // Get INFO/MAPABILITY.
            let mapability = {
                match record.info(b"MAPABILITY").float() {
                    Ok(ref mapability) => match mapability {
                        Some(ref mapability) => {
                            if mapability.len() > 0 {
                                mapability[0]
                            } else {
                                0_f32
                            }
                        }
                        None => 0_f32,
                    },
                    Err(_) => 0_f32,
                }
            };
            // Check whether window should be used for GC correction.
            let chrom = str::from_utf8(reader.header().rid2name(record.rid().unwrap()))
                .map_err(|e| format!("Cannot convert chrom name from UTF-8: {}", e))?
                .to_string();
            let use_chrom = re.is_match(&chrom);
            // Get the coverage.
            let coverage = record
                .format(b"COV")
                .float()
                .map_err(|e| format!("FORMAT/COV not found: {}", e))?[0][0];

            res.push(NormData::new(
                use_chrom,
                is_gap,
                is_blacklisted,
                gc_content.into(),
                mapability.into(),
                coverage.into(),
            ));
            // Update book-keeping.
            prev_rid = Some(record.rid());
        }
    }

    Ok(res)
}

fn chi_square_distance(a: f32, b: f32) -> f32 {
    let mu = (a + b) / 2.0;
    ((a - mu) * (a - mu) + (b - mu) * (b - mu)) / mu
}

/// Write out normalized coverage data again.
pub fn write_normalized_data(
    logger: &Logger,
    options: &Options,
    normalized_coverage: &[f32],
) -> Result<(), String> {
    // Block for BCF file reader and writer.
    {
        debug!(
            logger,
            "Opening input BCF file for generating final BCF file"
        );

        // Construct reader...
        let mut reader = bcf::Reader::from_path(options.input.clone())
            .map_err(|e| format!("Could not open input BCF file: {}", e))?;
        if options.io_threads > 0 {
            reader
                .set_threads(options.io_threads as usize)
                .map_err(|e| format!("Could not set I/O thread count {}", e))?;
        }

        // Construct writer...
        let mut writer = {
            // Construct extended header.
            let mut header = bcf::Header::with_template(reader.header());
            let lines = vec![
                "##FORMAT=<ID=OUTLIER,Number=1,Type=Integer,Description=\"Is outlier\">",
                "##FORMAT=<ID=NCOV,Number=1,Type=Float,Description=\"Normalized coverage\">",
                "##FORMAT=<ID=NCOV2,Number=1,Type=Float,Description=\"Normalized coverage \
                 (log2-scaled)\">",
            ];
            for line in lines {
                header.push_record(line.as_bytes());
            }
            header.push_record(format!("##cnvetti_normalizeVersion={}", VERSION).as_bytes());
            header.push_record(
                format!(
                    "##cnvetti_normalizeCommand={}",
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
        if options.io_threads > 0 {
            writer
                .set_threads(options.io_threads as usize)
                .map_err(|e| format!("Could not set I/O thread count: {}", e))?;
        }

        // Buffer to use for reading records.
        let mut record = reader.empty_record();
        // Previous reference ID, for printing messages to user.
        let mut prev_rid = Option::None;
        // We need a copy of the lightweight header.
        let header = reader.header().clone();
        // Bookkeeping, index in `normalized_coverage`.
        let mut i = 0;

        // Update all records from input.
        while reader.read(&mut record).is_ok() {
            if !prev_rid.is_some() || prev_rid.unwrap() != record.rid() {
                info!(
                    logger,
                    "Starting on contig {}",
                    str::from_utf8(header.rid2name(record.rid().unwrap())).unwrap()
                );
            }

            // Translate to output header.
            writer.translate(&mut record);

            // Write normalized coverage into record.
            let nrcs: Vec<f32> = vec![normalized_coverage[i] as f32; 1];
            record
                .push_format_float(b"NCOV", nrcs.as_slice())
                .map_err(|e| format!("Could not write FORMAT/NCOV: {}", e))?;

            let nrc2 = (normalized_coverage[i] + EPSILON).log2();
            let nrcs2: Vec<f32> = vec![
                if nrc2.is_finite() {
                    nrc2 as f32
                } else {
                    f32::missing()
                };
                1
            ];
            record
                .push_format_float(b"NCOV2", nrcs2.as_slice())
                .expect("Could not write FORMAT/NCOV2");

            // Compute flag for single point outlier.
            let is_outlier = if i > 0 && i + 1 < normalized_coverage.len() {
                const CHI_SQUAR_THRESH: f32 = 6.635;
                chi_square_distance(normalized_coverage[i - 1], normalized_coverage[i])
                    > CHI_SQUAR_THRESH
                    && chi_square_distance(normalized_coverage[i + 1], normalized_coverage[i])
                        > CHI_SQUAR_THRESH
            } else {
                false
            };
            let vals: Vec<i32> = vec![is_outlier as i32; 1];
            record
                .push_format_integer(b"OUTLIER", vals.as_slice())
                .map_err(|e| format!("Could not write FORMAT/OUTLIER: {}", e))?;

            // Finally, write out record again.
            writer.write(&record).expect("Could not write BCF record!");

            // Update book-keeping.
            prev_rid = Some(record.rid());
            i += 1;
        }
    }

    Ok(())
}
