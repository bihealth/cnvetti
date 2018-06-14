use std::env;

extern crate bio;

extern crate chrono;

#[macro_use]
extern crate error_chain;

extern crate clap;

extern crate separator;
use separator::Separatable;

#[macro_use]
extern crate slog;
use slog::Logger;

extern crate strum;
#[macro_use]
extern crate strum_macros;

extern crate rust_htslib;
use rust_htslib::bam::{self, Read};
use rust_htslib::bcf;

extern crate shlex;

extern crate lib_shared;
use lib_shared::bam_utils;
use lib_shared::bcf_utils;
use lib_shared::regions::GenomeRegions;

mod options;
pub use options::*;

mod errors {
    // Create the Error, ErrorKind, ResultExt, and Result types
    error_chain!{}
}

pub use errors::*;

/// Build header for the output BCF file.
///
/// This defines all values used throughout the whole window/target specific BAM files,
/// regardless whether they are actually used in the file.
///
/// Note that we use shared FORMAT tags for coverage and fragment count, such that we get
/// unified processing of copy number data in BCF files.
fn build_header(samples: &Vec<String>, contigs: &GenomeRegions, is_targeted: bool) -> bcf::Header {
    let mut header = bcf::Header::new();
    // TODO(holtgrewe): what about the VCF version?

    // Put overall meta information into the BCF header.
    let now = chrono::Utc::now();
    header.push_record(format!("##fileDate={}", now.format("%Y%m%d").to_string()).as_bytes());

    // Put creating tool version and call into file.
    header.push_record(format!("##cnvetti_cmdCoverageVersion={}", "0.1.0").as_bytes());
    header.push_record(
        format!(
            "##cnvetti_cmdCoverageCommand={}",
            env::args()
                .map(|s| shlex::quote(&s).to_string())
                .collect::<Vec<String>>()
                .join(" ")
        ).as_bytes(),
    );

    // Add samples to BCF header.
    for sample in samples {
        header.push_sample(sample.as_bytes());
    }

    // Put contig information into BCF header.
    for (name, _, length) in &contigs.regions {
        header.push_record(format!("##contig=<ID={},length={}>", name, length).as_bytes());
    }

    // Push the relevant header records.
    // TODO: later decide about commented-out lines
    let lines = vec![
        // Define ALT column <WINDOW>/<TARGET>
        "##ALT=<ID=WINDOW,Description=\"Record describes a window for read or coverage \
         counting\">",
        "##ALT=<ID=TARGET,Description=\"Record describes a target for fragment counting\">",
        // INFO fields describing the window
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Window end\">",
        "##INFO=<ID=GC,Number=1,Type=Float,Description=\"Reference GC fraction, if reference \
         FASTA file was given\">",
        "##INFO=<ID=MAPQ_MEAN,Number=1,Type=Float,Description=\"Mean MAPQ value for \
         approximating mapability\">",
        // "##INFO=<ID=BLACKLIST,Number=0,Type=Flag,Description=\"Interval overlaps with \
        //  blacklist site\">",
        "##INFO=<ID=GAP,Number=0,Type=Flag,Description=\"Window overlaps with N in \
         reference (gap)\">",
        // Generic FORMAT fields
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        // Note that the meaning of the term "coverage" differs considerably between
        // the different approaches. When processing targeted data, it is the length
        // normalized fraction of reads in this sample in the current target. In
        // contrast, for WGS or off-target reads, it is a normalized coverag where
        // distribution of all values on the autosomes is expected to center around
        // the center of diploid positions without any alterations.
        "##FORMAT=<ID=RCV,Number=1,Type=Float,Description=\"Raw coverage value\">",
        "##FORMAT=<ID=LCV,Number=1,Type=Float,Description=\"Length-normalized coverage \
         value\">",
        "##FORMAT=<ID=NCV,Number=1,Type=Float,Description=\"Further normalized coverage \
         value\">",
        "##FORMAT=<ID=CV,Number=1,Type=Float,Description=\"Finalized coverage value as input \
         to discovery step (or similar steps such as segmentation)\">",
        // FORMAT fields used when calling by coverage (deep WGS only)
        "##FORMAT=<ID=CVSD,Number=1,Type=Float,Description=\"Per-window coverage SD\">",
        // FORMAT fields used for read counts (off-target WES, shallow WGS)
        // TODO: use Character here with 'Y'/'N'
        "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Filter values for the genotype\">",
    ];
    for line in lines {
        header.push_record(line.as_bytes());
    }

    header
}

/// Build bcf::Writer with appropriate header.
fn build_bcf_writer(
    path: &String,
    samples: &Vec<String>,
    ref_contigs: &GenomeRegions,
    is_targeted: bool,
) -> Result<bcf::Writer> {
    let uncompressed = !path.ends_with(".bcf") && !path.ends_with(".vcf.gz");
    let vcf = path.ends_with(".vcf") || path.ends_with(".vcf.gz");

    let header = build_header(samples, ref_contigs, is_targeted);
    bcf::Writer::from_path(&path, &header, uncompressed, vcf)
        .chain_err(|| "Could not open BCF file for writing")
}

/// Main entry point for the `coverage` command.
pub fn run(logger: &mut Logger, options: &CoverageOptions) -> Result<()> {
    info!(logger, "Running: cnvetti cmd coverage");
    info!(logger, "Options: {:?}", options);

    // Warn about unusual parametrizations.
    if options.window_length.is_some() && options.reference.is_none() {
        warn!(
            logger,
            "You are counting per-window coverage but have specified no reference FASTA file.  \
             This is required for GC correction and usually desired.  Please reconsider your \
             choice!"
        );
    }

    // Opening BAM file.
    info!(logger, "Opening BAI-indexed BAM file...");
    let bam_reader = bam::IndexedReader::from_path(&options.input)
        .chain_err(|| "Could not open BAI-indexed BAM file")?;

    // Get list of regions to process.
    let ref regions = match &options.genome_region {
        Some(regions) => GenomeRegions::from_string_list(&vec![regions.clone()])
            .chain_err(|| "Problem parsing region")?,
        None => {
            bam_utils::build_chroms_bam(bam_reader.header(), Some(options.contig_regex.clone()))
                .chain_err(|| "Problem getting contig-length list from BAM header")?
        }
    };
    // Print list of regions that will be processed.
    info!(logger, "Will process {} regions", regions.regions.len());
    debug!(
        logger,
        "Regions: {:?}",
        regions
            .regions
            .iter()
            .map(|(name, _, len)| format!("{}:1-{}", name, len.separated_string()))
            .collect::<Vec<String>>()
    );

    // Create output file writer and kick off processing.  This is done in its own block such
    // that the file is definitely closed when building the index below.
    {
        let contigs =
            bam_utils::build_chroms_bam(bam_reader.header(), Some(options.contig_regex.clone()))
                .chain_err(|| "Problem getting contig-length list from BAM header")?;

        let samples = bam_utils::samples_from_file(&options.input)
            .chain_err(|| "Could not get sample names from BAM file")?;
        let mut writer = build_bcf_writer(
            &options.output,
            &samples,
            &contigs,
            options.targets_bed.is_some(),
        )?;
    }

    // Finally, create index on created output file.
    info!(logger, "Building index for output file...");
    bcf_utils::build_index(logger, &options.output).chain_err(|| "Could not build index")?;

    Ok(())
}
