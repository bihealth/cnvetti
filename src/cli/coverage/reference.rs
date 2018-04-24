/// Analysis of the reference sequence for GC content and gaps.
use std::path::Path;

use slog::Logger;

use bio::io::fasta;
use bio::utils::Text;

use separator::Separatable;

use rust_htslib::prelude::*;

/// Some statistics on the reference.
pub struct ReferenceStats {
    /// Reference sequence length.
    pub length: usize,
    /// GC content of windows.
    pub gc_content: Vec<f32>,
    /// Whether or not the window contains a gap (`N`).
    pub has_gap: Vec<bool>,
}

impl ReferenceStats {
    /// Create a new statistics from a given path.
    pub fn from_path<P: AsRef<Path>>(
        path: P,
        chrom: &str,
        window_length: usize,
        logger: &Logger,
    ) -> Result<Self, ReferenceStatsError> {
        match path.as_ref().to_str() {
            Some(p) if path.as_ref().exists() => {
                Ok(try!(Self::build(p.as_ref(), chrom, window_length, logger)))
            }
            _ => Err(ReferenceStatsError::InvalidPath),
        }
    }

    /// Internal builder function.
    fn build(
        path: &str,
        chrom: &str,
        window_length: usize,
        logger: &Logger,
    ) -> Result<Self, ReferenceStatsError> {
        info!(logger, "Loading GC content and gap (is-N) status...");

        let seq = Self::load_seq(path, chrom, logger)?;
        let length = seq.len();

        let (gc_content, has_gap) = Self::look_at_chars(&seq, window_length, logger)?;

        Ok(ReferenceStats {
            length,
            gc_content,
            has_gap,
        })
    }

    /// Load sequence from `path`.
    fn load_seq(path: &str, chrom: &str, logger: &Logger) -> Result<Text, ReferenceStatsError> {
        debug!(logger, "Loading reference sequence {}: {}...", path, chrom);

        trace!(logger, "Opening indexed FASTA reader");
        let mut ref_reader = match fasta::IndexedReader::from_file(&path) {
            Ok(reader) => reader,
            Err(_error) => {
                return Err(ReferenceStatsError::InvalidPath);
            }
        };

        trace!(logger, "Reading reference {}", chrom);
        let mut seq = Text::new();
        match ref_reader.read_all(chrom, &mut seq) {
            Ok(_) => (),
            Err(_e) => {
                return Err(ReferenceStatsError::LoadingFailed);
            }
        }
        let end = seq.len();
        debug!(
            logger,
            "Reference seq {} has {} characters",
            chrom,
            end.separated_string()
        );

        Ok(seq)
    }

    /// Perform analysis for GC ratio and "has N" flags.
    fn look_at_chars(
        seq: &Text,
        window_length: usize,
        logger: &Logger,
    ) -> Result<(Vec<f32>, Vec<bool>), ReferenceStatsError> {
        trace!(logger, "Analyzing sequence composition...");

        let num_buckets = ((seq.len() + window_length - 1) / window_length) as usize;

        // Count, GC characters, non-N chracters, and record "is N" flags.
        let mut gc_count = vec![0 as i32; num_buckets];
        let mut non_n_count = vec![0 as i32; num_buckets];
        let mut has_gap = vec![false; num_buckets];

        // Count GC chars and establish gap status.
        debug!(logger, "Counting GC and N characters...");
        for (i, c) in seq.iter().enumerate() {
            let bucket = i / window_length;

            match *c as char {
                'g' | 'c' | 'G' | 'C' => {
                    gc_count[bucket] += 1;
                    non_n_count[bucket] += 1;
                }
                'a' | 't' | 'A' | 'T' => {
                    non_n_count[bucket] += 1;
                }
                'n' | 'N' => {
                    has_gap[bucket] = true;
                }
                _ => (),
            }
        }

        // Compute GC content from counts.
        debug!(logger, "Computing GC content from GC counts...");
        let mut gc_ratio = vec![0 as f32; num_buckets];
        for (i, count) in gc_count.iter().enumerate() {
            gc_ratio[i] = if non_n_count[i] == 0 {
                f32::missing()
            } else {
                (*count as f32) / (non_n_count[i] as f32)
            };
        }

        Ok((gc_ratio, has_gap))
    }
}

quick_error! {
    /// Error for giving invalid (e.g., non-existing) path to reference FASTA file.
    #[derive(Debug, Clone)]
    pub enum ReferenceStatsError {
        InvalidPath {
            description("invalid path")
        }
        LoadingFailed {
            description("loading failed")
        }
    }
}
