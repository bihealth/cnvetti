// Implementation of the "coverage" command.

pub mod options;
pub use self::options::*;

use std::cmp::{max, min};
use std::fs::File;

use bio::data_structures::interval_tree;
use bio::io::fasta;

use rust_htslib::bam::{self, Read as BamRead};
use rust_htslib::bcf;
use rust_htslib::tbx::{self, Read as TbxRead};

use separator::Separatable;

use slog::Logger;

mod agg;
mod bcf_sink;
mod piles;
mod reference;
mod regions;

use self::agg::*;
use self::piles::PileCollector;

use self::bcf_sink::{sample_from_bam, CoverageBcfSink};
use self::reference::ReferenceStats;
use self::regions::GenomeRegions;

use cli::shared;

// TODO: remove restriction to same-length windows.
// TODO: implement overlapping windows

/// Load mapability if any.
fn load_mapability(
    logger: &mut Logger,
    options: &Options,
    chrom: &str,
    end: usize,
) -> Result<Option<Vec<f64>>, String> {
    // Open tabix reader when mapability BED file is given, otherwise return early with `None`.
    if options.mapability_bed.is_none() {
        return Ok(None);
    }
    let mut tbx_reader = tbx::Reader::from_path(options.mapability_bed.as_ref().unwrap())
        .map_err(|x| x.to_string())?;
    if options.io_threads > 0 {
        tbx_reader
            .set_threads(options.io_threads as usize)
            .expect("Could not set I/O thread count");
    }

    // Compute number of buckets and allocate array.
    let window_length = options.window_length as usize;
    let num_buckets = ((end + window_length - 1) / window_length) as usize;
    let mut result = vec![0_f64; num_buckets];

    // Get numeric index of chrom and fetch region.
    let tid = match tbx_reader.tid(chrom) {
        Ok(tid) => tid,
        Err(_) => {
            info!(logger, "No mapability information for {}", chrom);
            return Ok(Some(result));
        }
    };
    tbx_reader
        .fetch(tid, 0, end as u32)
        .map_err(|e| format!("Could not fetch chrom {}: {}", chrom, e))?;

    // TODO: can we make this loop tighter without unwrapping all the time?
    for buf in tbx_reader.records() {
        // TODO: mapability is shifted by 1/2*k to the right?
        // Parse begin/end/mapability from BED file.
        let s = String::from_utf8(buf.unwrap()).unwrap();
        let arr: Vec<&str> = s.split('\t').collect();
        if arr.len() != 4 {
            panic!("Mapability BED file had {} instead of 4 columns", arr.len());
        }
        let begin = arr[1].parse::<usize>().unwrap();
        let end = arr[2].parse::<usize>().unwrap();
        let mapability = arr[3].parse::<f64>().unwrap();

        // Modify buckets in `result`.
        let mut window_id = begin as usize / window_length;
        let window_begin = |window_id: usize| window_id * window_length;
        let window_end = |window_id: usize| window_begin(window_id) + window_length;
        while window_begin(window_id) < end {
            let len = min(window_end(window_id), end) - max(window_begin(window_id), begin);
            result[window_id as usize] += (len as f64 / window_length as f64) * mapability;
            window_id += 1;
        }
    }

    Ok(Some(result))
}

/// Load blacklist if any.
fn load_blacklist(
    logger: &mut Logger,
    options: &Options,
    chrom: &str,
    end: usize,
) -> Result<Option<interval_tree::IntervalTree<usize, bool>>, String> {
    // Open tabix reader when blacklist BED file is given, otherwise return early with `None`.
    if options.blacklist_bed.is_none() {
        return Ok(None);
    }
    let mut tbx_reader =
        tbx::Reader::from_path(options.blacklist_bed.as_ref().unwrap()).map_err(|x| x.to_string())?;
    if options.io_threads > 0 {
        tbx_reader
            .set_threads(options.io_threads as usize)
            .expect("Could not set I/O thread count");
    }

    // Compute number of buckets and allocate array.
    let mut result = interval_tree::IntervalTree::new();

    // Get numeric index of chrom and fetch region.
    let tid = match tbx_reader.tid(chrom) {
        Ok(tid) => tid,
        Err(_) => {
            info!(logger, "No blacklist information for {}", chrom);
            return Ok(Some(result));
        }
    };
    tbx_reader
        .fetch(tid, 0, end as u32)
        .map_err(|e| format!("Could not fetch chrom {}: {}", chrom, e))?;

    // TODO: can we make this loop tighter without unwrapping all the time?
    for buf in tbx_reader.records() {
        // TODO: mapability is shifted by 1/2*k to the right?
        // Parse begin/end from BED file.
        let s = String::from_utf8(buf.unwrap()).unwrap();
        let arr: Vec<&str> = s.split('\t').collect();
        if arr.len() != 3 {
            panic!("Blacklist BED file had {} instead of 3 columns", arr.len());
        }
        let begin = arr[1].parse::<usize>().unwrap();
        let end = arr[2].parse::<usize>().unwrap();
        result.insert(begin..end, true);
    }

    Ok(Some(result))
}

/// Collect pile information if any.
fn collect_pile_info(
    logger: &mut Logger,
    options: &Options,
    chrom: &str,
    out_bed: Option<&mut File>,
    depth_threshold: Option<usize>,
) -> Result<(usize, interval_tree::IntervalTree<u32, u32>), String> {
    debug!(logger, "Collecting piles on reference {}...", chrom);
    let mut collector = PileCollector::new(
        &options.input,
        out_bed,
        options.pile_depth_percentile,
        options.pile_max_gap,
        options.min_mapq,
    ).map_err(|e| e.to_string())?;
    let threshold_empty = depth_threshold.is_none();
    let (depth_threshold, len_sum, tree) = collector.collect_piles(chrom, logger, depth_threshold);

    debug!(
        logger,
        "Black-listed {} bp on ref seq {}",
        len_sum.separated_string(),
        chrom
    );
    if threshold_empty {
        debug!(
            logger,
            "Setting threshold to {}",
            depth_threshold.separated_string()
        );
    }

    Ok((depth_threshold, tree))
}

/// Process one region.
fn process_region(
    logger: &mut Logger,
    options: &Options,
    (chrom, start, end): &(String, usize, usize),
    out_bcf: &mut CoverageBcfSink,
    out_bed: Option<&mut File>,
    out_bam: Option<&mut bam::Writer>,
    pile_threshold: Option<usize>,
) -> Result<Option<usize>, String> {
    info!(logger, "Processing region {}:{}-{}", chrom, start, end);

    // Load statistics for the given contig.
    let ref_stats = ReferenceStats::from_path(
        &options.reference,
        chrom,
        options.window_length as usize,
        logger,
    ).map_err(|e| e.to_string())?;

    // Collect pile information, if requested and counting alignments.
    let (depth_threshold, tree) =
        if options.mask_piles && options.count_kind == CountKind::Alignments {
            debug!(logger, "Collecting piles...");
            let (depth_threshold, tree) =
                collect_pile_info(logger, options, chrom, out_bed, pile_threshold)?;
            (Some(depth_threshold), Some(tree))
        } else {
            debug!(logger, "Not collecting piles...");
            (None, None)
        };

    // Collect mapability information, if any.
    debug!(logger, "Collecting mapability...");
    let mapability = load_mapability(logger, options, chrom, *end)?;

    // Collect blacklist information, if any.
    debug!(logger, "Collecting blacklist...");
    let blacklist = load_blacklist(logger, options, chrom, *end)?;

    // Construct aggregator for the records.
    let mut aggregator: Box<BamRecordAggregator> = match options.count_kind {
        CountKind::Alignments => Box::new(CountAlignmentsAggregator::new(
            tree.as_ref(),
            options.clone(),
            *end,
            out_bam,
        )),
        CountKind::Coverage => Box::new(CoverageAggregator::new(options.clone(), *end)),
    };

    // Jump to region with BAM reader.
    let mut bam_reader = bam::IndexedReader::from_path(&options.input).map_err(|e| e.to_string())?;
    if options.io_threads > 0 {
        bam_reader
            .set_threads(options.io_threads as usize)
            .map_err(|e| e.to_string())?;
    }
    let tid: u32 = bam_reader.header().tid(chrom.as_bytes()).unwrap();
    bam_reader
        .fetch(tid, (start - 1) as u32, *end as u32)
        .map_err(|e| {
            format!(
                "Could not seek to region {}:{}-{}: {}",
                chrom, start, end, e
            )
        })?;

    // Main loop for region: pass all BAM records in region through aggregator.
    info!(logger, "Computing coverage...");
    aggregator.put_fetched_records(&mut bam_reader);
    debug!(
        logger,
        "Processed {}, skipped {} records ({:.2}% are off-target)",
        aggregator.num_processed().separated_string(),
        aggregator.num_skipped().separated_string(),
        100.0 * (aggregator.num_processed() - aggregator.num_skipped()) as f64
            / aggregator.num_processed() as f64,
    );

    // Create the BCF records for this region.
    info!(logger, "Writing BCF with coverage information...");

    // Get `u32` 0-base coordinates.
    for (wid, gc) in ref_stats.gc_content.iter().enumerate() {
        let mut record = out_bcf.writer.empty_record();
        let rid = out_bcf
            .writer
            .header()
            .name2rid(chrom.as_bytes())
            .map_err(|e| format!("Could not find contig {}: {}", chrom, e))?;

        // Columns: CHROM, POS, ID, REF, ALT, (FILTER)
        let window_length = options.window_length as usize;
        let pos = wid * window_length;
        let window_end = min(*end as usize, (wid + 1) * window_length);
        let alleles_v = vec![Vec::from("N"), Vec::from("<COUNT>")];
        let alleles = alleles_v
            .iter()
            .map(|x| x.as_slice())
            .collect::<Vec<&[u8]>>();

        record.inner_mut().rid = rid as i32;
        record.set_pos(pos as i32);
        record
            .set_id(format!("WIN_{}_{}_{}", chrom, pos + 1, window_end).as_bytes())
            .map_err(|e| format!("Could not update ID: {}", e))?;
        record
            .set_alleles(&alleles)
            .map_err(|e| format!("Could not update alleles: {}", e))?;

        // Columns: INFO
        record
            .push_info_integer(b"END", &[window_end as i32])
            .map_err(|e| format!("Could not write INFO/END: {}", e))?;
        if ref_stats.has_gap[wid] {
            record
                .push_info_flag(b"GAP")
                .map_err(|e| format!("Could not write INFO/GAP: {}", e))?;
        }
        record
            .push_info_float(b"GC", &[*gc as f32])
            .map_err(|e| format!("Could not write INFO/GC: {}", e))?;
        if let Some(ref mapability) = &mapability {
            record
                .push_info_float(b"MAPABILITY", &[mapability[wid] as f32])
                .map_err(|e| format!("Could not write INFO/MAPABILITY: {}", e))?;
        }
        if let Some(ref blacklist) = &blacklist {
            if blacklist.find(pos..window_end).peekable().peek().is_some() {
                record
                    .push_info_flag(b"BLACKLIST")
                    .map_err(|e| format!("Could not write INFO/BLACKLIST: {}", e))?;
            }
        }

        // Columns: FORMAT/GT
        record
            .push_format_integer(b"GT", &[bcf::GT_MISSING])
            .map_err(|e| format!("Could not write FORMAT/GT: {}", e))?;

        // Columns: FORMAT/COV etc.
        for field in aggregator.character_field_names() {
            let value: String = aggregator
                .character_values(wid as u32)
                .get(&field)
                .unwrap()
                .to_string();
            record
                .push_format_char(field.as_bytes(), value.as_bytes())
                .map_err(|e| format!("Could not write FORMAT/{}: {}", field, e))?;
        }
        for field in aggregator.integer_field_names() {
            let value = *aggregator.integer_values(wid as u32).get(&field).unwrap();
            record
                .push_format_integer(field.as_bytes(), &[value])
                .map_err(|e| format!("Could not write FORMAT/{}: {}", field, e))?;
        }
        for field in aggregator.float_field_names() {
            let value = *aggregator.float_values(wid as u32).get(&field).unwrap();
            record
                .push_format_float(field.as_bytes(), &[value])
                .map_err(|e| format!("Could not write FORMAT/{}: {}", field, e))?;
        }

        out_bcf
            .writer
            .write(&record)
            .map_err(|e| format!("Could not write BCF record: {}", e))?;
    }

    Ok(depth_threshold)
}

/// Main entry point for the "coverage" command.
pub fn call(logger: &mut Logger, options: &Options) -> Result<(), String> {
    info!(logger, "Configuration: {:?}", options);

    // Get list of regions to process.
    let regions = GenomeRegions::from_list_or_path(&options.genome_regions, &options.reference)
        .map_err(|e| e.to_string())?;

    // Open output file (in its own block so we can create the index below).
    {
        let index = fasta::Index::from_file(&format!("{}.fai", options.reference))
            .map_err(|e| format!("Could not read fasta index: {}", e.to_string()))?;
        let mut out_bcf = CoverageBcfSink::from_path(
            &options.output,
            &sample_from_bam(&options.input)?,
            &index,
            logger,
        )?;
        // Output file for writing BED file with blocked regions.
        let mut out_bed = match options.output_bed {
            Some(ref path) => Some(File::create(&path).map_err(|e| e.to_string())?),
            None => None,
        };
        // Output file for writing BAM file with passing alignments.
        let mut out_bam = match options.output_bam {
            Some(ref path) => {
                let bam_reader =
                    bam::IndexedReader::from_path(&options.input).map_err(|e| e.to_string())?;
                let header = bam::Header::from_template(bam_reader.header());
                Some(bam::Writer::from_path(path, &header).map_err(|e| e.to_string())?)
            }
            None => None,
        };

        // The threshold for pileup masking is estimated from chr1.
        let mut pileup_threshold = None;

        // Process each region.
        for region in regions.regions {
            pileup_threshold = process_region(
                logger,
                &options,
                &region,
                &mut out_bcf,
                out_bed.as_mut(),
                out_bam.as_mut(),
                pileup_threshold,
            ).map_err(|e| e.to_string())?;
        }
    }

    // Build index on the output file.
    shared::build_index(logger, &options.output);

    Ok(())
}
