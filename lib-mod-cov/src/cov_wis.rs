//! Code for computing model-based coverage for the within-sample model.

use std::collections::HashMap;
use std::env;

extern crate shlex;

extern crate clap;

extern crate rust_htslib;
use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::{self, Read};

extern crate lib_shared;
use lib_shared::stats::Stats;

use super::errors::*;

/// Load the normalized counts from input BCF file.
pub fn load_counts(input: &String) -> Result<HashMap<String, f64>> {
    let mut reader = bcf::Reader::from_path(&input).chain_err(|| "Could not open input BCF file")?;
    let mut normalized_counts = HashMap::new();
    let mut record = reader.empty_record();
    loop {
        match reader.read(&mut record) {
            Ok(_) => record.unpack(),
            Err(bcf::ReadError::NoMoreRecord) => break,
            _ => bail!("Could not read record"),
        }

        let name = String::from_utf8(record.id()).chain_err(|| "Problem converting from UTF-8")?;

        normalized_counts.insert(
            name,
            record
                .format(b"CV")
                .float()
                .chain_err(|| "Could not access field FORMAT/CV")?[0][0] as f64,
        );
    }

    Ok(normalized_counts)
}

/// Information to be written to output.
///
/// Only to be created for non-filtered probes.
#[derive(Debug, Clone)]
pub struct ProbeInfo {
    // Normalized coverage at the probe.
    norm_cov: f64,

    // Reference mean.
    ref_median: f64,
    // Reference std deviation.
    ref_median_abs_dev: f64,
}

impl ProbeInfo {
    /// Return relative coverage.
    fn cov_rel(&self) -> f64 {
        self.norm_cov / self.ref_median
    }

    /// Return coverage as Z-score of reference.
    fn cov_z_score(&self) -> f64 {
        (self.norm_cov - self.ref_median) / self.ref_median_abs_dev
    }
}

/// Load reference and compute per-target coverage information to write out to result file.
pub fn load_target_infos(
    input_ref: &String,
    normalized_counts: &HashMap<String, f64>,
) -> Result<HashMap<String, ProbeInfo>> {
    let mut reader =
        bcf::Reader::from_path(&input_ref).chain_err(|| "Could not open reference BCF file")?;

    let mut record = reader.empty_record();
    let mut probe_infos = HashMap::new();

    loop {
        match reader.read(&mut record) {
            Ok(_) => {
                // Skip any filtered record from reference.
                if record
                    .filters()
                    .map(|id| {
                        String::from_utf8(record.header().id_to_name(id))
                            .expect(&format!("Could not transflate from UTF-8"))
                    })
                    .any(|s| s != "PASS")
                {
                    continue;
                }
                else
                {
                    record.unpack();
                }
            }
            Err(bcf::ReadError::NoMoreRecord) => break,
            _ => bail!("Could not read record"),
        }

        let name = String::from_utf8(record.id()).chain_err(|| "Cannot decode from UTF-8")?;

        // Get normalized count at this probe and the empirical distribution from the reference.
        let norm_cov = *normalized_counts
            .get(&name)
            .expect(&format!("Could not get norm count for probe: {}", name));
        let empty: Vec<&[u8]> = Vec::new();
        // TODO: the following can be simplified with flatten.
        let mut names = Vec::new();
        let id = String::from_utf8(record.id()).chain_err(|| "Cannot decode from UTF-8")?;
        let tmp = record
            .info(b"REF_TARGETS")
            .string()
            .chain_err(|| "Problem accessing REF_TARGETS")?
            .unwrap_or(empty);
        for vec in &tmp {
            let s: String = String::from_utf8(vec.to_vec()).expect("Decoding from UTF-8 failed");
            let mut xs: Vec<String> = s.split(",").map(|s| s.to_string()).collect::<Vec<String>>();
            names.append(&mut xs);
        }
        let ref_norm_covs = names
            .iter()
            .map(|s| {
                *normalized_counts
                    .get(s)
                    .expect(format!("Could not find: {}", s).as_ref())
            })
            .collect::<Vec<f64>>();

        // Compute initial set of statistics.
        if !ref_norm_covs.is_empty() {
            let ref_median = ref_norm_covs.as_slice().median();
            let ref_median_abs_dev = ref_norm_covs.as_slice().median_abs_dev();

            probe_infos.insert(
                id,
                ProbeInfo {
                    norm_cov,
                    ref_median,
                    ref_median_abs_dev,
                },
            );
        }
    }

    Ok(probe_infos)
}

/// Read through input file again, and write output file.
pub fn write_mod_cov(
    probe_infos: &HashMap<String, ProbeInfo>,
    input: &String,
    output: &String,
) -> Result<()> {
    // Open input file.
    let mut reader = bcf::Reader::from_path(&input).chain_err(|| "Could not open input BCF file")?;

    // Build output header.
    let mut header = bcf::Header::from_template(reader.header());
    header.push_record(format!("##cnvetti_cmdModCoverageVersion={}", "0.1.0").as_bytes());
    header.push_record(
        format!(
            "##cnvetti_cmdModCoverageCommand={}",
            env::args()
                .map(|s| shlex::quote(&s).to_string())
                .collect::<Vec<String>>()
                .join(" ")
        ).as_bytes(),
    );

    // Open output file.
    let uncompressed = !output.ends_with(".bcf") && !output.ends_with(".vcf.gz");
    let vcf = output.ends_with(".vcf") || output.ends_with(".vcf.gz");
    let mut writer = bcf::Writer::from_path(&output, &header, uncompressed, vcf)
        .chain_err(|| "Could not open output BCF file")?;

    // Read records from input, modify, and write out again.
    let mut record = reader.empty_record();
    loop {
        match reader.read(&mut record) {
            Ok(_) => record.unpack(),
            Err(bcf::ReadError::NoMoreRecord) => break,
            _ => bail!("Could not read record"),
        }
        writer.translate(&mut record);

        let target_id = String::from_utf8(record.id())
            .map_err(|e| format!("Cannot decode ID from UTF-8: {}", e))?;

        if let Some(probe_info) = probe_infos.get(&target_id) {
            record
                .push_info_float(b"REF_MEAN", &[probe_info.ref_median as f32])
                .chain_err(|| "Could not write INFO/REF_MEAN")?;
            record
                .push_info_float(b"REF_STD_DEV", &[probe_info.ref_median_abs_dev as f32])
                .chain_err(|| "Could not write INFO/REF_STD_DEV")?;

            record
                .push_format_float(b"CV", &[probe_info.cov_rel() as f32])
                .chain_err(|| "Could not write FORMAT/CV")?;
            record
                .push_format_float(b"CVZ", &[probe_info.cov_z_score() as f32])
                .chain_err(|| "Could not write FORMAT/CVZ")?;
            if probe_info.cov_rel().is_finite() {
                let cov2 = if probe_info.cov_rel() == 0.0 {
                    f32::missing()
                } else {
                    probe_info.cov_rel().log2() as f32
                };
                record
                    .push_format_float(b"CV2", &[cov2])
                    .chain_err(|| "Could not write FORMAT/CV2")?;
            }
        } else {
            record.push_filter(
                writer
                    .header()
                    .name_to_id(b"NO_REF")
                    .expect("FILTER 'NO_REF' unknown"),
            );
        }

        writer
            .write(&record)
            .chain_err(|| "Problem writing BCF record.")?;
    }

    Ok(())
}
