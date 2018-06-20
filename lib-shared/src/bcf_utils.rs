use slog::Logger;

use rust_htslib;
use rust_htslib::bcf;
use std::ffi;

mod errors {
    // Create the Error, ErrorKind, ResultExt, and Result types
    error_chain!{}
}

pub use self::errors::*;

use regions::GenomeRegions;

/// Generate list of all contigs from BCF header.
pub fn extract_chroms(header: &bcf::header::HeaderView) -> GenomeRegions {
    let mut result = GenomeRegions::new();

    for ref record in header.header_records() {
        if let bcf::HeaderRecord::Contig { key, values } = record {
            assert_eq!(key, "contig");
            let mut name: Option<String> = None;
            let mut length: Option<u32> = None;
            for (ref key, ref value) in values {
                if key.as_str() == "ID" {
                    name = Some(value.to_string());
                } else if key.as_str() == "length" {
                    length = Some(value.parse::<u32>().unwrap());
                }
            }
            if let (Some(ref name), Some(length)) = (&name, length) {
                result.regions.push((name.clone(), 0, length as usize));
            } else {
                panic!(
                    "Could not parse both name/length from {:?}/{:?}",
                    name, length
                );
            }
        }
    }

    result
}

/// Build index file for the VCF/BCF file at `path`.
pub fn build_index(logger: &mut Logger, path: &String) -> Result<()> {
    if path.ends_with(".vcf.gz") {
        info!(logger, "Writing .tbi index file...");
        // TODO: the unsafe stuff should go into rust_htslib
        use rust_htslib;
        use std::ffi;
        let res = unsafe {
            rust_htslib::htslib::tbx_index_build(
                ffi::CString::new(path.as_bytes()).unwrap().as_ptr(),
                0,
                &rust_htslib::htslib::tbx_conf_vcf,
            )
        };
        if res != 0 {
            bail!("Could not create .tbi index");
        }
    } else if path.ends_with(".bcf") {
        info!(logger, "Writing .csi index file...");
        let res = unsafe {
            rust_htslib::htslib::bcf_index_build(
                ffi::CString::new(path.as_bytes()).unwrap().as_ptr(),
                14,
            )
        };
        if res != 0 {
            bail!("Could not create .csi index");
        }
    } else {
        info!(
            logger,
            "Not building index, output file does not having ending .bcf or .vcf.gz"
        );
    }

    Ok(())
}
