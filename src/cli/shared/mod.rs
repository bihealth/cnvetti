pub mod lowess;
pub mod math;

pub use self::lowess::{lowess, LowessResults};
use slog::Logger;

/// Build index file for the VCF/BCF file at `path`.
pub fn build_index(logger: &mut Logger, path: &String) {
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
            panic!("Could not create .tbi index");
        }
    } else if path.ends_with(".bcf") {
        info!(logger, "Writing .csi index file...");
        use rust_htslib;
        use std::ffi;
        let res = unsafe {
            rust_htslib::htslib::bcf_index_build(
                ffi::CString::new(path.as_bytes()).unwrap().as_ptr(),
                14,
            )
        };
        if res != 0 {
            panic!("Could not create .tbi index");
        }
    } else {
        info!(
            logger,
            "Not building index, output file does not having ending .bcf or .vcf.gz"
        );
    }
}

/// Tokenize genome region strings with 1-based positions.
pub fn tokenize_genome_regions(regions: &Vec<String>) -> Vec<(String, u64, u64)> {
    regions
        .iter()
        .map(|region| {
            let region_split = region.split(":").collect::<Vec<&str>>();
            if region_split.len() != 2 {
                panic!("Invalid region: {}", region);
            }

            let number_split = region_split[1].split("-").collect::<Vec<&str>>();
            if number_split.len() != 2 {
                panic!("Invalid region: {}", region);
            }

            let begin = number_split[0]
                .to_string()
                .replace(",", "")
                .parse::<u64>()
                .unwrap();
            let end = number_split[1]
                .to_string()
                .replace(",", "")
                .parse::<u64>()
                .unwrap();

            (region_split[0].to_string(), begin, end)
        })
        .collect()
}
