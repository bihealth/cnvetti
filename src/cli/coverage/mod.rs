// Implementation of the "coverage" command.

include!(concat!(env!("OUT_DIR"), "/version.rs"));

use std::cmp::{max, min};
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::Write;
use std::iter;

use bio::data_structures::interval_tree;
use bio::io::fasta;
use bio::utils::Text;

use chrono;

use regex::Regex;

use rust_htslib::bam::{self, Read as BamRead};
use rust_htslib::bcf;
use rust_htslib::tbx::{self, Read as TbxRead};

use separator::Separatable;

use serde_json;

use shlex;

use slog::Logger;

mod agg;
mod options;
mod piles;
pub mod summary; // TODO: move into mod shared?

use self::agg::*;
pub use self::options::*;
use self::piles::PileCollector;
use self::summary::CoverageSummarizer;

use cli::shared;

// TODO: remove restriction to same-length windows.
// TODO: implement overlapping windows
// TODO: filter to autosomal + sex chromosomes, ignore the rest

/// Structure holding the input files readers and related data structures.
struct CoverageInput {
    /// Reader for FAI-indexed FASTA file.
    ref_reader: fasta::IndexedReader<File>,
    /// An optional BED reader for the --mapability-bed file.
    tbx_reader: Option<tbx::Reader>,
    /// All indexed BAM readers.
    bam_readers: Vec<bam::IndexedReader>,
    /// Tokenized 1-based genomic region as `(chrom, begin, end)`.
    genome_regions: Vec<(String, u64, u64)>,
    /// Samples from all BAM files.
    samples: Vec<String>,
    /// Mapping from sample name to 0-based numeric index from all BAM files.
    sample_to_index: HashMap<String, u32>,
    /// Mapping from read group to sample name.
    read_group_to_sample: HashMap<String, String>,
}

/// Parse @RG lane into triple (id, sm).
fn parse_line_rg(line: String) -> Option<(String, String)> {
    let line_split = line.split("\t");
    let mut id: Option<String> = None;
    let mut sm: Option<String> = None;
    for s in line_split {
        let token: Vec<&str> = s.split(":").collect();
        if token.len() >= 2 {
            match token[0] {
                "ID" => {
                    id = Some(token[1].to_string());
                }
                "SM" => {
                    sm = Some(token[1].to_string());
                }
                _ => (),
            }
        }
    }

    match (id, sm) {
        (Some(id), Some(sm)) => Some((id, sm)),
        _ => None,
    }
}

/// Tokenize genome region strings with 1-based positions.
fn tokenize_genome_regions(regions: &Vec<String>) -> Vec<(String, u64, u64)> {
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

impl CoverageInput {
    /// Open input file handlers.
    fn open_files(options: &Options, logger: &Logger) -> CoverageInput {
        // Open BAM files.
        let bam_readers: Vec<bam::IndexedReader> = options
            .input
            .iter()
            .map(|path| match bam::IndexedReader::from_path(path) {
                Ok(mut reader) => {
                    if options.io_threads > 0 {
                        reader
                            .set_threads(options.io_threads as usize)
                            .expect("Could not set I/O thread count");
                    }
                    reader
                }
                Err(error) => {
                    panic!("Could create BAM reader for path {}! {:?}", path, error);
                }
            })
            .collect();
        // Build mapping from sample name (from read group) to numeric ID.
        let mut samples = Vec::new();
        let mut read_group_to_sample: HashMap<String, String> = HashMap::new();
        let mut sample_to_index: HashMap<String, u32> = HashMap::new();
        let mut idx: u32 = 0;
        for reader in &bam_readers {
            let text = String::from_utf8(Vec::from(reader.header().as_bytes())).unwrap();
            // The one sample for the current reader.
            let mut seen_sm = Option::None;
            // All samples seen so far.
            let mut seen_before = samples.clone();

            for line in text.lines() {
                if line.starts_with("@RG") {
                    match parse_line_rg(line.to_string()) {
                        Some((id, sm)) => {
                            debug!(logger, "RG '{}' => '{}'", &id, &sm);
                            debug!(logger, "SM '{}' => '{}'", &sm, &idx);

                            // Protect against two different samples in the same file.
                            seen_sm = match seen_sm {
                                Some(seen_sm) => {
                                    if seen_sm != sm {
                                        panic!("Seen more than one sample in same BAM file");
                                    }
                                    Some(seen_sm)
                                }
                                None => Some(sm.clone()),
                            };

                            // Protect against same sample in two files.
                            if seen_before.contains(&sm) {
                                panic!("Seen sample in more than one file");
                            }

                            samples.push(sm.clone());
                            read_group_to_sample.insert(id.clone(), sm.clone());
                            sample_to_index.insert(sm.clone(), idx.clone());
                            idx += 1;
                        }
                        None => (),
                    }
                }
            }
        }

        CoverageInput {
            // Create FASTA reader for reference.
            ref_reader: match fasta::IndexedReader::from_file(&options.reference) {
                Ok(reader) => reader,
                Err(error) => {
                    panic!("Could not open indexed FASTA File! {:?}", error);
                }
            },
            // Create reader for mapability BED file, if any.
            tbx_reader: match &options.mapability_bed {
                &Some(ref path) => match tbx::Reader::from_path(&path) {
                    Ok(mut reader) => {
                        if options.io_threads > 0 {
                            reader
                                .set_threads(options.io_threads as usize)
                                .expect("Could not set I/O thread count");
                        }
                        Some(reader)
                    }
                    Err(error) => {
                        panic!("Could create BED reader {:?}", error);
                    }
                },
                &None => None,
            },
            // Parse out the genomic regions.
            genome_regions: tokenize_genome_regions(&options.genome_regions),
            // Create readers for BAM files.
            bam_readers: bam_readers,
            // The samples in the order defined in input files.
            samples: samples,
            // Mapping from sample name to index.
            sample_to_index: sample_to_index,
            // Mapping from read group name to sample name.
            read_group_to_sample: read_group_to_sample,
        }
    }
}

/// Structure holding the output files readers and related data structures.
#[derive(Debug)]
struct CoverageOutput {
    header: bcf::Header,
    bcf_writer: bcf::Writer,
}

impl CoverageOutput {
    /// Open output file handles.
    pub fn open_files(
        options: &Options,
        input: &CoverageInput,
        logger: &mut Logger,
    ) -> CoverageOutput {
        let uncompressed =
            !options.output.ends_with(".bcf") && !options.output.ends_with(".vcf.gz");
        let vcf = options.output.ends_with(".vcf") || options.output.ends_with(".vcf.gz");

        let header = CoverageOutput::build_header(&input);
        info!(logger, "Writing BCF header {:?}", &header);

        CoverageOutput {
            bcf_writer: match bcf::Writer::from_path(&options.output, &header, uncompressed, vcf) {
                Ok(mut writer) => {
                    if options.io_threads > 0 {
                        writer
                            .set_threads(options.io_threads as usize)
                            .expect("Could not set I/O thread count");
                    }
                    writer
                }
                Err(error) => {
                    panic!(
                        "Could not open BCF file for output {}. {:?}",
                        options.output, error
                    );
                }
            },
            header: header,
        }
    }

    /// Create BCF header from input files.
    fn build_header(input: &CoverageInput) -> bcf::Header {
        let mut header = bcf::Header::new();

        // TODO(holtgrewe): compare contig names in BAM and FAI file?
        // TODO(holtgrewe): what about the VCF version?

        // Put overall meta information into the BCF header.
        let now = chrono::Utc::now();
        header.push_record(format!("##fileDate={}", now.format("%Y%m%d").to_string()).as_bytes());

        // Put creating tool version an dcall into file.
        header.push_record(format!("##cnvetti_coverageVersion={}", VERSION).as_bytes());
        header.push_record(
            format!(
                "##cnvetti_coverageCommand={}",
                env::args()
                    .map(|s| shlex::quote(&s).to_string())
                    .collect::<Vec<String>>()
                    .join(" ")
            ).as_bytes(),
        );

        // Put contig information into BCF header.
        for seq in input.ref_reader.index.sequences() {
            header.push_record(format!("##contig=<ID={},length={}>", seq.name, seq.len).as_bytes());
        }

        // Add samples to BCF header.
        for sample in &input.samples {
            header.push_sample(sample.as_bytes());
        }

        // Write out FILTER, INFO, FORMAT, and ALT fields.
        // TODO: some of these could go away eventually
        let lines = vec![
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Window end\">",
            "##INFO=<ID=GC,Number=1,Type=Float,Description=\"Reference GC content in percent\">",
            "##INFO=<ID=MAPABILITY,Number=1,Type=Float,Description=\"Mean mapability in the \
             window\">",
            "##INFO=<ID=GAP,Number=1,Type=Integer,Description=\"Window overlaps with N in \
            reference (gap)\">",
            "##INFO=<ID=GCWINDOWS,Number=1,Type=Integer,Description=\"Number of windows with same \
             GC content\">",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##FORMAT=<ID=COV,Number=1,Type=Float,Description=\"Average coverage\">",
            "##FORMAT=<ID=WINSD,Number=1,Type=Float,Description=\"Per-window coverage SD)\">",
            "##FORMAT=<ID=MP,Description=\"Masked for sample because too much masked because of \
             piles\",Type=Integer,Number=1>",
            "##FORMAT=<ID=RC,Number=1,Type=Integer,Description=\"Number of aligning reads, \
             scaled in the case of pile-masking\">",
            "##FORMAT=<ID=RRC,Number=1,Type=Integer,Description=\"Raw number of aligning reads\">",
            "##FORMAT=<ID=MS,Number=1,Type=Integer,Description=\"Number of bases in window \
             masked because of read piles\">",
            "##ALT=<ID=COUNT,Description=\"Record describes a window for read counting\">",
            // TODO: FEW_GCWINDOWS should go into normalize but then cannot be found?
            "##FILTER=<ID=FEW_GCWINDOWS,Description=\"Masked because of few windows with \
             this GC content\">",
        ];
        for line in lines {
            header.push_record(line.as_bytes());
        }

        header
    }
}

/// Implementation of the `cnvetti coverage` command.
pub struct CoverageApp<'a> {
    /// Logging struct
    logger: &'a mut Logger,

    /// Configuration of the application.
    options: Options,

    /// Input files.
    input: CoverageInput,

    /// Summariser.
    summariser: CoverageSummarizer,

    /// Depth threshold.
    depth_threshold: Option<usize>,
}

impl<'a> CoverageApp<'a> {
    /// Construct and initialize new `CoverageApp` instance.
    pub fn new(logger: &'a mut Logger, options: &Options) -> CoverageApp<'a> {
        info!(logger, "Running cnvetti coverage");

        let options = options.clone();
        info!(logger, "Configuration: {:?}", &options);

        info!(logger, "Opening input files and parsing genomic regions...");
        let input = CoverageInput::open_files(&options, logger);

        // Computation of summaries.
        CoverageApp {
            summariser: CoverageSummarizer::new(
                (options.gc_step * 1000.0) as i32,
                &input.samples,
                &options.count_kind.to_string(),
            ),
            logger,
            options,
            input,
            depth_threshold: None,
        }
    }

    /// Kick-off the actual processing.
    pub fn run(&mut self) -> Result<(), String> {
        info!(self.logger, "Processing");
        let regions = if self.input.genome_regions.is_empty() {
            self.gen_all_regions()
        } else {
            self.input.genome_regions.clone()
        };

        // Open output files and start processing in its own scope.
        {
            // The output file has to be opened here so it can go out of scope and be closed
            // before we build the index.
            info!(self.logger, "Opening output files...");
            let output = CoverageOutput::open_files(&self.options, &self.input, self.logger);

            self.process_regions(output, &regions);
        }

        // Build index on the output file.
        shared::build_index(&mut self.logger, &self.options.output);

        info!(self.logger, "Writing statistics files...");
        // Write out statistics JSON file.
        {
            let path_stats = self.options.output.clone() + ".stats.txt"; // TODO => .json
            let mut file_stats =
                File::create(path_stats).expect("Could not open stats file for writing");
            let summaries = self.summariser.summaries();
            file_stats
                .write_all(serde_json::to_string(&summaries).unwrap().as_bytes())
                .expect("Could not write statistics to output file");
        }

        info!(self.logger, "Done. Have a nice day!");

        Ok(())
    }

    /// Process all given regions.
    fn process_regions(&mut self, mut output: CoverageOutput, regions: &Vec<(String, u64, u64)>) {
        for &(ref chrom, start, end) in regions {
            self.process_region(&mut output, &chrom, start, end);
        }
    }

    /// Process one region.
    fn process_region(&mut self, output: &mut CoverageOutput, chrom: &str, start: u64, end: u64) {
        info!(
            self.logger,
            "Processing region {}:{}-{}",
            chrom,
            start.separated_string(),
            end.separated_string()
        );

        // Determin whether or not to collect statistics.
        let re = Regex::new(&self.options.contig_regex).unwrap();
        let collect_stats = re.is_match(&chrom);
        debug!(
            self.logger,
            "{} statistics for from {}",
            if collect_stats {
                "Collect"
            } else {
                "Do not collect"
            },
            chrom
        );

        // Get `u32` 0-based coordinates.
        let begin = (start - 1) as u32;
        let end = end as u32;

        let (chrom_len, gc_content, has_gap) = self.analyze_reference(chrom).unwrap();
        let mapability = self.maybe_load_mapability(chrom, chrom_len);
        let piles = self.maybe_collect_piles(chrom);

        info!(self.logger, "Computing coverage...");
        // Seek to region in readers.
        for bam_reader in &mut self.input.bam_readers {
            let tid = bam_reader
                .header()
                .tid(chrom.as_bytes())
                .expect("Could not resolve contig name to integer");
            bam_reader
                .fetch(tid, begin, end)
                .expect("Could not fetch region!");
        }

        // Initialize the BAM aggregation.
        // TODO: this could be made prettier with factory pattern (or the rust equivalent).
        let mut aggregators: Vec<Box<BamRecordAggregator>> = self.input
            .bam_readers
            .iter()
            .enumerate()
            .map(|(i, _)| {
                let result: Box<BamRecordAggregator> = match self.options.count_kind {
                    CountKind::Coverage => Box::new(CoverageAggregator::new(
                        self.input.samples.clone(),
                        self.input.read_group_to_sample.clone(),
                        self.input.sample_to_index.clone(),
                        self.options.clone(),
                        chrom_len,
                    )),
                    CountKind::Alignments => Box::new(CountAlignmentsAggregator::new(
                        &piles[i],
                        self.input.samples.clone(),
                        self.input.read_group_to_sample.clone(),
                        self.input.sample_to_index.clone(),
                        self.options.clone(),
                        chrom_len,
                    )),
                };

                result
            })
            .collect();

        // Count reads / compute coverage.
        let mut record = bam::Record::new();
        for (i, ref mut bam_reader) in self.input.bam_readers.iter_mut().enumerate() {
            // Read all records in region for reader `i`.
            while bam_reader.read(&mut record).is_ok() {
                aggregators[i].put_bam_record(&mut record);
            }
            debug!(
                self.logger,
                "For sample {} -- processed {}, skipped {} records ({:.2}% are off-target)",
                self.input.samples[i],
                aggregators[i].num_processed().separated_string(),
                aggregators[i].num_skipped().separated_string(),
                100.0 * (aggregators[i].num_processed() - aggregators[i].num_skipped()) as f64
                    / aggregators[i].num_processed() as f64,
            );
        }

        // Compute histogram of GC contents.
        let gc_step = self.options.gc_step;
        let mut gc_histo: HashMap<i32, i32> = HashMap::new();
        for gc in gc_content.iter() {
            let rounded_gc: i32 = (*gc as f64 / gc_step) as i32;
            let count = gc_histo.entry(rounded_gc).or_insert(0);
            *count += 1;
        }

        // Write out records to BCF file.
        let num_windows = gc_content.len();
        info!(
            self.logger,
            "Writing {} windows to BCF file...",
            num_windows.separated_string()
        );
        let rid = output
            .bcf_writer
            .header()
            .name2rid(chrom.as_bytes())
            .unwrap() as i32;
        for (wid, gc) in gc_content.iter().enumerate() {
            let rounded_gc = (*gc as f64 / gc_step) as i32;

            let mut record = output.bcf_writer.empty_record();
            // Columns: CHROM, POS, ID, REF, ALT, FILTER
            record.inner_mut().rid = rid;
            let pos = wid * self.options.window_length as usize;
            record.set_pos(pos as i32);
            let window_end = min(
                end as usize,
                (wid + 1) * self.options.window_length as usize,
            ) as i32;
            record
                .set_id(format!("WIN_{}_{}_{}", chrom, pos + 1, window_end).as_bytes())
                .expect("Could not update ID");
            let alleles = &[Vec::from("N"), Vec::from("<COUNT>")];
            record
                .set_alleles(alleles)
                .expect("Could not update alleles!");

            // INFO fields
            record
                .push_info_integer(b"END", &[window_end])
                .expect("Could not write INFO/END");
            record
                .push_info_float(b"GC", &[*gc])
                .expect("Could not write INFO/GC");
            if let Some(ref mapability) = &mapability {
                record
                    .push_info_float(b"MAPABILITY", &[mapability[wid] as f32])
                    .expect("Could not write INFO/MAPABILITY");
            }
            record
                .push_info_integer(b"GAP", &[has_gap[wid] as i32])
                .expect("Could not write INFO/GAP");
            record
                .push_info_integer(b"GCWINDOWS", &[*gc_histo.get(&rounded_gc).unwrap_or(&0)])
                .expect("Could not write INFO/GCWINDOWS");

            // FORMAT/GT
            record
                .push_format_integer(
                    b"GT",
                    iter::repeat(bcf::GT_MISSING)
                        .take(self.input.sample_to_index.len())
                        .collect::<Vec<i32>>()
                        .as_slice(),
                )
                .expect("Could not push genotypes");

            // The remaining FORMAT fields are done by `BamRecordAggregator`s.
            for field in aggregators[0].integer_field_names() {
                let values: Vec<i32> = aggregators
                    .iter()
                    .map(|ref agg| *agg.integer_values(wid as u32).get(&field).unwrap())
                    .collect();
                record
                    .push_format_integer(field.as_bytes(), values.as_slice())
                    .expect("Could not write FORMAT field");
            }
            for field in aggregators[0].float_field_names() {
                let values: Vec<f32> = aggregators
                    .iter()
                    .map(|ref agg| *agg.float_values(wid as u32).get(&field).unwrap())
                    .collect();
                record
                    .push_format_float(field.as_bytes(), values.as_slice())
                    .expect("Could not write FORMAT field");
            }

            // Update summary statistics.
            if collect_stats {
                for (sample_id, agg) in aggregators.iter().enumerate() {
                    if !agg.is_masked(wid as u32) {
                        self.summariser.increment(
                            &self.input.samples[sample_id],
                            *gc,
                            agg.get_coverage(wid as u32),
                        );
                    }
                }
            }

            // Actually write the record.
            output
                .bcf_writer
                .write(&record)
                .expect("Could not write record!");
        }
    }

    /// Build Vec of chromosomes as regions.
    fn gen_all_regions(&self) -> Vec<(String, u64, u64)> {
        self.input
            .ref_reader
            .index
            .sequences()
            .iter()
            .map(|seq| (seq.name.clone(), 1, seq.len))
            .collect()
    }

    /// Load mapability if configured to do so.
    fn maybe_load_mapability(&mut self, chrom: &str, chrom_len: u64) -> Option<Vec<f64>> {
        if self.input.tbx_reader.is_some() {
            let tid = self.input.tbx_reader.as_ref().unwrap().tid(chrom);
            match tid {
                Ok(_) => {
                    info!(self.logger, "Loading mapability...");
                    Some(
                        self.load_mapability(&chrom, chrom_len)
                            .expect("loading mapability failed"),
                    )
                }
                _ => {
                    info!(self.logger, "No mapability for chrom {}", chrom);
                    None
                }
            }
        } else {
            info!(self.logger, "Mapability is not considered.");
            None
        }
    }

    /// Load mapability of `chrom` from BED file.
    fn load_mapability(&mut self, chrom: &str, chrom_len: u64) -> Result<Vec<f64>, String> {
        // Compute number of buckets and allocate array.
        let window_length = self.options.window_length as u64;
        let num_buckets = ((chrom_len + window_length - 1) / window_length) as usize;
        let mut result = vec![0_f64; num_buckets];

        // Get numeric index of chrom and fetch region.
        let tid = self.input
            .tbx_reader
            .as_ref()
            .map(|reader| reader.tid(chrom).expect("sequence not found"))
            .unwrap();
        match self.input
            .tbx_reader
            .as_mut()
            .map(|reader| reader.fetch(tid, 0, chrom_len as u32))
            .unwrap()
        {
            Err(x) => panic!("Could not fetch in mapability BED file {:?}", x),
            _ => (),
        };

        // TODO: can we make this loop tighter without unwrapping all the time?
        for buf in self.input.tbx_reader.as_mut().unwrap().records() {
            // TODO: mapability is shifted by 1/2*k to the right?
            // Parse begin/end/mapability from BED file.
            let s = String::from_utf8(buf.unwrap()).unwrap();
            let arr: Vec<&str> = s.split('\t').collect();
            if arr.len() != 4 {
                panic!("Mapability BED file had {} instead of 4 columns", arr.len());
            }
            let begin = arr[1].parse::<u64>().unwrap();
            let end = arr[2].parse::<u64>().unwrap();
            let mapability = arr[3].parse::<f64>().unwrap();

            // Modify buckets in `result`.
            let mut window_id = begin as u64 / window_length;
            let window_begin = |window_id: u64| window_id * window_length;
            let window_end = |window_id: u64| window_begin(window_id) + window_length;
            while window_begin(window_id) < end {
                let len = min(window_end(window_id), end) - max(window_begin(window_id), begin);
                result[window_id as usize] += (len as f64 / window_length as f64) * mapability;
                window_id += 1;
            }
        }

        Ok(result)
    }

    /// Load piles for black-listing or an empty `IntervalTree` for each sample.
    fn maybe_collect_piles(&mut self, chrom: &str) -> Vec<interval_tree::IntervalTree<u32, u32>> {
        // Only collect piles when masking by pile and counting alignments.
        if self.options.mask_piles && self.options.count_kind == CountKind::Alignments {
            info!(self.logger, "Computing piles for black listing");
            let options = self.options.clone();
            let samples = self.input.samples.clone();

            let mut result = Vec::new();
            for (i, bam_reader) in self.input.bam_readers.iter_mut().enumerate() {
                let mut collector = PileCollector::new(
                    bam_reader,
                    options.pile_depth_percentile,
                    options.pile_max_gap,
                    options.min_mapq,
                    options.pile_mask_window_size,
                );
                let (depth_threshold, len_sum, tree) =
                    collector.collect_piles(chrom, self.logger, self.depth_threshold);
                debug!(
                    self.logger,
                    "Black listed {} bp for {}",
                    len_sum.separated_string(),
                    samples[i]
                );
                if let None = self.depth_threshold {
                    self.depth_threshold = Some(depth_threshold);
                }

                result.push(tree);
            }

            result
        } else {
            info!(self.logger, "Piles for black listing are not considered");
            (0..(self.input.bam_readers.len()))
                .map(|_| interval_tree::IntervalTree::new())
                .collect()
        }
    }

    /// Analyze reference for GC content, and is-gap (=N) status.
    fn analyze_reference(&mut self, chrom: &str) -> Result<(u64, Vec<f32>, Vec<bool>), String> {
        info!(self.logger, "Loading GC content and gap (is-N) status...");
        let window_length = self.options.window_length as usize;

        // Read reference sequence.
        debug!(self.logger, "Loading reference sequence {}...", chrom);
        let mut seq = Text::new();
        match self.input.ref_reader.read_all(chrom, &mut seq) {
            Ok(_) => (),
            Err(e) => {
                return Err(format!("Could not read reference. {:?}", e));
            }
        }
        let end = seq.len();
        debug!(
            self.logger,
            "Reference seq {} has {} characters",
            chrom,
            end.separated_string()
        );

        let num_buckets = ((end + window_length - 1) / window_length) as usize;

        let mut gc_count = vec![0 as i32; num_buckets];
        let mut has_gap = vec![false; num_buckets];

        // Count GC chars and establish gap status.
        debug!(self.logger, "Counting GC and N characters...");
        for (i, c) in seq.iter().enumerate() {
            let bucket = i / window_length;

            // TODO: Ns should not be counted as non-GC, are not sequence.
            match *c as char {
                'G' | 'C' => {
                    gc_count[bucket] += 1;
                }
                'N' => {
                    has_gap[bucket] = true;
                }
                _ => (),
            }
        }

        // Compute GC content from count.
        debug!(self.logger, "Computing GC content from GC counts...");
        let mut result = vec![0 as f32; num_buckets];
        for (i, count) in gc_count.iter().enumerate() {
            let bucket_end = (i + 1) * window_length;
            let bucket_len = if end > bucket_end {
                window_length
            } else {
                end - i * window_length
            };
            result[i] = (*count as f32) / (bucket_len as f32);
        }

        Ok((end as u64, result, has_gap))
    }
}

/// Main entry point for the "coverage" command.
pub fn call(logger: &mut Logger, options: &Options) -> Result<(), String> {
    let mut app = CoverageApp::new(logger, options);
    app.run()
}
