/// Shared code between the different normalization variants.
use std::env;

use shlex;

use slog::Logger;

use rust_htslib::bcf::{self, Read};

use super::errors::*;

/// Read records from input file, write to output file, processing each record to throug the
/// functor.
pub fn process_bcf(
    logger: &mut Logger,
    input: &String,
    output: &String,
    io_threads: u32,
    f: &Fn(&mut bcf::Record) -> Result<()>,
) -> Result<()> {
    debug!(logger, "Opening input file...");
    let mut reader =
        bcf::Reader::from_path(&input).chain_err(|| format!("Could not open BCF file: {}", input))?;
    if io_threads > 0 {
        reader
            .set_threads(io_threads as usize)
            .chain_err(|| "Could not set threads to BCF reader")?;
    }

    debug!(logger, "Opening output file...");
    let mut writer = {
        // Construct extended header.
        let mut header = bcf::Header::from_template(reader.header());
        // TODO: version should come from one central place
        header.push_record(format!("##cnvetti_cmdNormalizeVersion={}", "0.2.0").as_bytes());
        header.push_record(
            format!(
                "##cnvetti_wiseNormalizeCommand={}",
                env::args()
                    .map(|s| shlex::quote(&s).to_string())
                    .collect::<Vec<String>>()
                    .join(" ")
            ).as_bytes(),
        );

        let uncompressed = !output.ends_with(".bcf") && !output.ends_with(".vcf.gz");
        let vcf = output.ends_with(".vcf") || output.ends_with(".vcf.gz");

        bcf::Writer::from_path(output.clone(), &header, uncompressed, vcf)
            .chain_err(|| format!("Could not open output BCF file {}", output))?
    };

    let mut record = reader.empty_record();
    loop {
        match reader.read(&mut record) {
            Ok(_) => (),
            Err(bcf::ReadError::NoMoreRecord) => break,
            _ => bail!("Problem reading BCF record"),
        }
        writer.translate(&mut record);
        f(&mut record)?;
        writer
            .write(&record)
            .chain_err(|| "Problem writing BCF record to file")?;
    }

    Ok(())
}
