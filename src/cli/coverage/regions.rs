/// Helper code for handling genomic regions.
use std::path::Path;

/// Representation of a list of genomic regions.
pub struct GenomeRegions {
    /// The region specification as (chr, start, end).
    regions: Vec<(String, usize, usize)>,
}

impl GenomeRegions {
    /// Parse `GenomeRegions` from vector of `String`s if given, otherwise load from FAI index.
    pub fn from_list_or_path<P: AsRef<Path>>(
        strings: Vec<String>,
        path: P,
    ) -> Result<Self, GenomeRegionsError> {
        Ok(if !strings.is_empty() {
            try!(Self::from_list(&strings))
        } else {
            try!(Self::from_path(&path))
        })
    }

    /// Load from FAI file.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, GenomeRegionsError> {
        Err(GenomeRegionsError::InvalidRegion)
    }

    /// Parse from list of strings.
    pub fn from_list(strings: &Vec<String>) -> Result<Self, GenomeRegionsError> {
        Ok(GenomeRegions {
            regions: strings
                .iter()
                .map(|region| {
                    let region_split = region.split(":").collect::<Vec<&str>>();
                    if region_split.len() != 2 {
                        return Err(GenomeRegionsError::InvalidRegion);
                    }

                    let number_split = region_split[1].split("-").collect::<Vec<&str>>();
                    if number_split.len() != 2 {
                        return Err(GenomeRegionsError::InvalidRegion);
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

                    Ok((region_split[0].to_string(), start, end))
                })
                .collect()?,
        })
    }
}

quick_error! {
    /// Error construction `GenomeRegions`.
    #[derive(Debug, Clone)]
    pub enum GenomeRegionsError {
        InvalidPath {
            description("invalid path")
        }
        LoadingFailed {
            description("loading failed")
        }
        InvalidRegion {
            description("invalid region")
        }
    }
}
