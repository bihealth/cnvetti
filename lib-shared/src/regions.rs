/// Helper code for handling genomic regions.
use std::path::Path;

use bio::io::fasta;

mod errors {
    // Create the Error, ErrorKind, ResultExt, and Result types
    error_chain! {}
}

pub use self::errors::*;

/// Representation of a list of genomic regions.
#[derive(Debug, Clone)]
pub struct GenomeRegions {
    /// The region specification as (chr, start, end).
    pub regions: Vec<(String, usize, usize)>,
}

impl GenomeRegions {
    /// Create empty one.
    pub fn new() -> Self {
        GenomeRegions {
            regions: Vec::new(),
        }
    }

    /// Parse `GenomeRegions` from vector of `String`s if given, otherwise load from FAI index.
    pub fn from_string_list_or_path<P: AsRef<Path>>(
        strings: &Vec<String>,
        path: P,
    ) -> Result<Self> {
        Ok(if !strings.is_empty() {
            try!(Self::from_string_list(strings))
        } else {
            try!(Self::from_path(path))
        })
    }

    /// Load from FAI file.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let ref_reader = fasta::IndexedReader::from_file(&path)
            .chain_err(|| "Loading FAI for genome regions failed")?;

        Ok(GenomeRegions {
            regions: ref_reader
                .index
                .sequences()
                .iter()
                .map(|ref seq| (seq.name.clone(), 1, seq.len as usize))
                .collect(),
        })
    }

    /// Parse from list of strings.
    pub fn from_string_list(strings: &Vec<String>) -> Result<Self> {
        let regions: Result<Vec<(String, usize, usize)>> = strings
            .iter()
            .map(|region| {
                let region_split = region.split(":").collect::<Vec<&str>>();
                if region_split.len() != 2 {
                    bail!("Could not parse genome region (not exactly one ':')");
                }

                let number_split = region_split[1].split("-").collect::<Vec<&str>>();
                if number_split.len() != 2 {
                    bail!("Could not parse genome region (not exatly one '-')");
                }

                let start = number_split[0]
                    .to_string()
                    .replace(",", "")
                    .replace("_", "")
                    .parse::<usize>()
                    .unwrap();
                let end = number_split[1]
                    .to_string()
                    .replace(",", "")
                    .replace("_", "")
                    .parse::<usize>()
                    .unwrap();

                Ok((region_split[0].to_string(), start - 1, end))
            })
            .collect();
        Ok(GenomeRegions { regions: regions? })
    }

    /// Construct from list of names and lengths.
    pub fn from_name_and_length(pairs: &Vec<(String, usize)>) -> Self {
        GenomeRegions {
            regions: pairs
                .iter()
                .map(|(chr, len)| (chr.clone(), 0, *len))
                .collect::<Vec<(String, usize, usize)>>(),
        }
    }
}
