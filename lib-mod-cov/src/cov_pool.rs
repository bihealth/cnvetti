//! Code for computing model-based coverage for the pool-based model.

use std::collections::HashMap;
use std::env;

extern crate shlex;

extern crate clap;

extern crate rust_htslib;
use rust_htslib::bcf::record::Numeric;
use rust_htslib::bcf::{self, Read};

use slog::Logger;

extern crate lib_shared;
use lib_shared::bcf_utils;

use super::errors::*;

/// Summary statistics of one region.
#[derive(Debug, Default)]
pub struct RegionInfo {
    /// Coverage mean in region.
    pub mean: f64,
    /// Standard deviation of coverage in region.
    pub std_dev: f64,
    /// Five-number summary of coverage in region.
    pub summary5: [f64; 5],
}

/// Load region information by region name (targeted sequencing) for contig.
fn load_region_infos(
    path_model: &String,
    contig: &String,
    contig_length: u32,
) -> Result<HashMap<String, RegionInfo>> {
    let mut result = HashMap::new();

    let mut reader = bcf::IndexedReader::from_path(&path_model)
        .chain_err(|| "Could not open reference BCF file")?;
    let rid = reader
        .header()
        .name2rid(contig.as_bytes())
        .chain_err(|| format!("Could not translate header name {}", contig))?;
    reader
        .fetch(rid, 0, contig_length)
        .chain_err(|| "Could not jump to region")?;

    let mut record = reader.empty_record();
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

        let mean = record.info(b"CV_MEAN").float().unwrap().unwrap()[0];
        let std_dev = record.info(b"CV_STD_DEV").float().unwrap().unwrap()[0];
        let summary5 = record
            .info(b"CV_5SUMMARY")
            .float()
            .unwrap()
            .unwrap()
            .iter()
            .map(|x| *x as f64)
            .collect::<Vec<f64>>();

        let mut region_info: RegionInfo = Default::default();
        region_info.mean = mean.into();
        region_info.std_dev = std_dev.into();
        region_info.summary5.copy_from_slice(&summary5);

        result.insert(name, region_info);
    }

    Ok(result)
}

/// Process one contig.
fn process_contig(
    logger: &mut Logger,
    contig: &String,
    contig_length: u32,
    path_model: &String,
    reader: &mut bcf::IndexedReader,
    writer: &mut bcf::Writer,
) -> Result<()> {
    debug!(logger, "Load region summary from model BCF file");
    let region_infos = load_region_infos(path_model, contig, contig_length)?;
    debug!(logger, "=> done");

    debug!(logger, "Jumping to region in coverage file");
    let rid = reader
        .header()
        .name2rid(contig.as_bytes())
        .chain_err(|| format!("Could not translate header name {}", contig))?;
    reader
        .fetch(rid, 0, contig_length)
        .chain_err(|| "Could not jump to region")?;
    debug!(logger, "=> done");

    let mut record = reader.empty_record();
    loop {
        match reader.read(&mut record) {
            Ok(_) => record.unpack(),
            Err(bcf::ReadError::NoMoreRecord) => break,
            _ => bail!("Could not read record"),
        }
        writer.translate(&mut record);

        let region_id = String::from_utf8(record.id())
            .map_err(|e| format!("Cannot decode ID from UTF-8: {}", e))?;

        if let Some(region_info) = region_infos.get(&region_id) {
            record
                .push_info_float(b"REF_MEAN", &[region_info.mean as f32])
                .chain_err(|| "Could not write INFO/REF_MEAN")?;
            record
                .push_info_float(b"REF_STD_DEV", &[region_info.std_dev as f32])
                .chain_err(|| "Could not write INFO/REF_STD_DEV")?;
            record
                .push_info_float(
                    b"REF_IQR",
                    &[(region_info.summary5[3] - region_info.summary5[1]) as f32],
                )
                .chain_err(|| "Could not write INFO/REF_IQR")?;
            record
                .push_info_float(
                    b"REF_5SUMMARY",
                    &[
                        region_info.summary5[0] as f32,
                        region_info.summary5[1] as f32,
                        region_info.summary5[2] as f32,
                        region_info.summary5[3] as f32,
                        region_info.summary5[4] as f32,
                    ],
                )
                .chain_err(|| "Could not write INFO/REF_5SUMMARY")?;

            let cov: f32 = record.format(b"CV").float().unwrap()[0][0];
            let cov_rel = if region_info.mean == 0.0 {
                f32::missing()
            } else {
                cov / region_info.mean as f32
            };
            let cov_z_score = if region_info.std_dev == 0.0 || !cov_rel.is_finite() {
                f32::missing()
            } else {
                (cov - region_info.mean as f32) / region_info.std_dev as f32
            };

            record
                .push_format_float(b"CV", &[cov_rel])
                .chain_err(|| "Could not write FORMAT/CV")?;
            record
                .push_format_float(b"CVZ", &[cov_z_score])
                .chain_err(|| "Could not write FORMAT/CVZ")?;
            if cov_rel.is_finite() {
                let cov2 = if cov_rel == 0.0 || !cov_rel.is_finite() {
                    f32::missing()
                } else {
                    cov_rel.log2() as f32
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

/// Perform pool-based calling.
pub fn run(
    logger: &mut Logger,
    path_model: &String,
    path_input: &String,
    path_output: &String,
) -> Result<()> {
    // Open input BCF file for reading.
    info!(
        logger,
        "Opening input file {} and extracting contig", path_input
    );
    let mut reader =
        bcf::IndexedReader::from_path(path_input).chain_err(|| "Could not open BCF input file")?;

    // Get contig names.
    let contigs = bcf_utils::extract_chroms(&reader.header());

    // Open output BCF file for writing.
    info!(logger, "Opening output file {}", path_output);

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

    // Actually open output file.
    let uncompressed = !path_output.ends_with(".bcf") && !path_output.ends_with(".vcf.gz");
    let vcf = path_output.ends_with(".vcf") || path_output.ends_with(".vcf.gz");
    let mut writer = bcf::Writer::from_path(&path_output, &header, uncompressed, vcf)
        .chain_err(|| "Could not open output BCF file")?;

    // Process each contig individually.
    for (contig, _, length) in &contigs.regions {
        info!(logger, "Processing contig {}", &contig);
        process_contig(
            logger,
            contig,
            *length as u32,
            path_model,
            &mut reader,
            &mut writer,
        )?;
    }

    Ok(())
}
