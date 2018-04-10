// Implementation of the "coverage" command.

include!(concat!(env!("OUT_DIR"), "/version.rs"));

use std::collections::HashMap;
use std::cmp::{max, min};
use std::env;
use std::iter;
use std::fs::File;

use bio::io::fasta;
use bio::utils::Text;

use chrono;

use clap::{ArgMatches, Values};

use rust_htslib::bcf;
use rust_htslib::tbx::{self, Read as TbxRead};
use rust_htslib::bam::{self, Read as BamRead};

use separator::Separatable;

use shlex;

use slog::Logger;

// TODO: check that no two BAM files have overlapping sample names.
// TODO: remove restriction to same-length windows.

/// Enum for selecting preset.
#[derive(Clone, Debug)]
enum OptionsPreset {
    Wgs,
    WesOffTarget,
    WesOnTarget,
}


impl OptionsPreset {
    /// Parse `CountKind` from `&str`.
    fn from_str(s: &str) -> Option<OptionsPreset> {
        match s {
            "Wgs" => Some(OptionsPreset::Wgs),
            "WesOffTarget" => Some(OptionsPreset::WesOffTarget),
            "WesOnTarget" => Some(OptionsPreset::WesOnTarget),
            _ => {
                panic!("Invalid preset {}", s);
            }
        }
    }
}


/// Enum for selecting count type.
#[derive(Clone, Debug)]
enum CountKind {
    Coverage,
    Alignments,
}


impl CountKind {
    /// Parse `CountKind` from `&str`.
    fn from_str(s: &str) -> Option<CountKind> {
        match s {
            "Coverage" => Some(CountKind::Coverage),
            "Alignments" => Some(CountKind::Alignments),
            _ => {
                panic!("Invalid count type {}", s);
            }
        }
    }
}


/// Options for the "coverage" command.
#[derive(Clone, Debug)]
pub struct Options {
    reference: String,
    input: Vec<String>,
    output: String,
    genome_regions: Vec<String>,
    mapability_bed: Option<String>,
    window_length: u64,
    window_overlap: u64,
    count_kind: CountKind,
    min_mapq: u8,
    min_unclipped: f32,
    skip_discordant: bool,
}

impl Options {
    /// Build options from ArgMatches.
    pub fn new(matches: &ArgMatches) -> Options {
        // TODO: interpret `--preset`

        let count_kind = matches.value_of("count_kind").unwrap();

        Options {
            reference: matches.value_of("reference").unwrap().to_string(),
            input: matches
                .values_of("input")
                .unwrap_or(Values::default())
                .map(|res| res.to_string())
                .collect(),
            output: matches.value_of("output").unwrap().to_string(),
            genome_regions: matches
                .values_of("genome_regions")
                .unwrap_or(Values::default())
                .map(|res| res.to_string())
                .collect(),
            mapability_bed: match matches.value_of("mapability_bed") {
                Some(x) => Some(x.to_string()),
                None => None,
            },
            window_length: matches
                .value_of("window_length")
                .unwrap()
                .parse::<u64>()
                .unwrap(),
            window_overlap: matches
                .value_of("window_overlap")
                .unwrap()
                .parse::<u64>()
                .unwrap(),
            count_kind: CountKind::from_str(count_kind).unwrap(),
            min_mapq: matches.value_of("min_mapq").unwrap().parse::<u8>().unwrap(),
            min_unclipped: matches
                .value_of("min_unclipped")
                .unwrap()
                .parse::<f32>()
                .unwrap(),
            skip_discordant: matches.is_present("skip_discordant"),
        }
    }
}


/// Struct with common information for aggregator.
struct BaseAggregator {
    /// Names of the samples to aggregate in this sample.
    samples: Vec<String>,
    /// Mapping from read group name to sample.
    read_group_to_sample: HashMap<String, String>,
    /// Mapping from sample name to sample index.
    sample_to_idx: HashMap<String, u32>,

    /// Configuration.
    options: Options,
    /// Length of contig process
    contig_length: u64,
}


impl BaseAggregator {
    fn new(
        samples: &Vec<String>,
        read_group_to_sample: &HashMap<String, String>,
        sample_to_idx: &HashMap<String, u32>,
        options: &Options,
        contig_length: u64,
    ) -> BaseAggregator {
        BaseAggregator {
            samples: samples.clone(),
            read_group_to_sample: read_group_to_sample.clone(),
            sample_to_idx: sample_to_idx.clone(),
            options: options.clone(),
            contig_length,
        }
    }
}


/// Trait for alignment aggregation from BAM files.
trait BamRecordAggregator {
    /// Put the given record into the aggregation.
    fn put_bam_record(&mut self, record: &bam::Record);

    /// Write out the counts to the given `record`.
    fn fill_bcf_record(&self, record: &mut bcf::Record, window_id: u32);
}


/// Struct for aggregating as alignment counts.
struct CountAlignmentsAggregator {
    /// Common information from all aggregators.
    base: BaseAggregator,

    /// Per-sample read counts (each bin is a `u32`).
    counters: Vec<Vec<u32>>,
}


impl CountAlignmentsAggregator {
    /// Construct new aggregator with the given BAM `header`.  This information is
    /// necessary to appropriately allocate buffers for all samples in the header.
    fn new(
        samples: &Vec<String>,
        read_group_to_sample: &HashMap<String, String>,
        sample_to_idx: &HashMap<String, u32>,
        options: &Options,
        contig_length: u64,
    ) -> CountAlignmentsAggregator {
        CountAlignmentsAggregator {
            base: BaseAggregator::new(
                samples,
                read_group_to_sample,
                sample_to_idx,
                options,
                contig_length,
            ),
            counters: Vec::new(),
        }
    }
}


impl BamRecordAggregator for CountAlignmentsAggregator {
    fn put_bam_record(&mut self, record: &bam::Record) {}

    fn fill_bcf_record(&self, record: &mut bcf::Record, window_id: u32) {}
}


/// Bin for coverage aggregation.
struct CoverageBin {}


/// Struct for aggregating as coverage.
struct CoverageAggregator {
    /// Common information from all aggregators.
    base: BaseAggregator,
    /// Number of bases for each base in the current bin for each sample.
    coverage: Vec<Vec<CoverageBin>>,
}


impl CoverageAggregator {
    fn new(
        samples: &Vec<String>,
        read_group_to_sample: &HashMap<String, String>,
        sample_to_idx: &HashMap<String, u32>,
        options: &Options,
        contig_length: u64,
    ) -> CoverageAggregator {
        CoverageAggregator {
            base: BaseAggregator::new(
                samples,
                read_group_to_sample,
                sample_to_idx,
                options,
                contig_length,
            ),
            coverage: Vec::new(),
        }
    }
}

impl BamRecordAggregator for CoverageAggregator {
    fn put_bam_record(&mut self, record: &bam::Record) {}

    fn fill_bcf_record(&self, record: &mut bcf::Record, window_id: u32) {}
}


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
                Ok(reader) => reader,
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

            for line in text.lines() {
                if line.starts_with("@RG") {
                    match parse_line_rg(line.to_string()) {
                        Some((id, sm)) => {
                            debug!(logger, "RG '{}' => '{}'", &id, &sm);
                            debug!(logger, "SM '{}' => '{}'", &sm, &idx);
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
                &Some(ref path) => {
                    match tbx::Reader::from_path(&path) {
                        Ok(reader) => Some(reader),
                        Err(error) => {
                            panic!("Could create BED reader {:?}", error);
                        }
                    }
                }
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
        let uncompressed = !options.output.ends_with(".bcf") &&
            !options.output.ends_with(".vcf.gz");
        let vcf = options.output.ends_with(".vcf") || options.output.ends_with(".vcf.gz");

        let header = CoverageOutput::build_header(&input);
        info!(logger, "Writing BCF header {:?}", &header);

        CoverageOutput {
            bcf_writer: match bcf::Writer::from_path(&options.output, &header, uncompressed, vcf) {
                Ok(writer) => writer,
                Err(error) => {
                    panic!(
                        "Could not open BCF file for output {}. {:?}",
                        options.output,
                        error
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
        header.push_record(
            format!("##fileDate={}", now.format("%Y%m%d").to_string()).as_bytes(),
        );

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
            header.push_record(
                format!("##contig=<ID={},length={}>", seq.name, seq.len).as_bytes(),
            );
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
            "##INFO=<ID=GAP,Number=0,Type=Flag,Description=\"Window overlaps with N in reference \
            (gap)\">",
            "##INFO=<ID=GCWINDOWS,Number=1,Type=Integer,Description=\"Number of windows with same \
            GC content\">",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "##FORMAT=<ID=COV,Number=1,Type=Float,Description=\"Average coverage with non-q0 \
            reads\">",
            "##FORMAT=<ID=COV0,Number=1,Type=Float,Description=\"Average coverage with \
            q0+non-q0 reads\">",
            "##FORMAT=<ID=WINSD,Number=1,Type=Float,Description=\"Per-window coverage SD (non-q0 \
            reads)\">",
            "##FORMAT=<ID=WINSD0,Number=1,Type=Float,Description=\"Per-window coverage SD \
            (q0+non-q0 reads)\">",
            "##FORMAT=<ID=RC,Number=1,Type=Float,Description=\"Number of aligning non-q0 \
            reads\">",
            "##FORMAT=<ID=RC0,Number=1,Type=Float,Description=\"Number of aligning q0+non-q0 \
            reads\">",
            "##ALT=<ID=COUNT,Description=\"Record describes a window for read counting\">",
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
    /// Output files.
    output: CoverageOutput,
}

impl<'a> CoverageApp<'a> {
    /// Construct and initialize new `CoverageApp` instance.
    pub fn new(logger: &'a mut Logger, options: &Options) -> CoverageApp<'a> {
        info!(logger, "Running cnvetti coverage");

        let options = options.clone();
        info!(logger, "Configuration: {:?}", &options);

        info!(logger, "Opening input files and parsing genomic regions...");
        let input = CoverageInput::open_files(&options, logger);

        info!(logger, "Opening output files...");
        let output = CoverageOutput::open_files(&options, &input, logger);

        CoverageApp {
            logger,
            options,
            input,
            output,
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
        self.process_regions(&regions);

        info!(self.logger, "Writing statistics");
        info!(self.logger, "Done. Have a nice day!");

        Ok(())
    }

    /// Process all given regions.
    fn process_regions(&mut self, regions: &Vec<(String, u64, u64)>) {
        for &(ref chrom, start, end) in regions {
            self.process_region(&chrom, start, end);
        }
    }

    /// Process one region.
    fn process_region(&mut self, chrom: &str, start: u64, end: u64) {
        info!(
            self.logger,
            "Processing region {}:{}-{}",
            chrom,
            start.separated_string(),
            end.separated_string()
        );

        info!(self.logger, "Loading GC content and gap (is-N) status...");
        let begin = start - 1; // 1-based to 0-based
        let (chrom_len, gc_content, has_gap) = self.analyze_reference(chrom).unwrap();

        let mapability = if self.input.tbx_reader.is_some() {
            info!(self.logger, "Loading mapability...");
            Some(self.load_mapability(&chrom, chrom_len).expect(
                "loading mapability failed",
            ))
        } else {
            info!(self.logger, "Mapability is not considered.");
            None
        };

        info!(self.logger, "Processing alignments...");
        // Seek to region in readers.
        for bam_reader in &mut self.input.bam_readers {
            let tid = bam_reader.header().tid(chrom.as_bytes()).expect(
                "Could not resolve contig name to integer",
            );
            bam_reader.fetch(tid, start as u32, end as u32).expect(
                "Could not fetch region!",
            );
        }

        // Initialize the BAM aggregation.
        // TODO: this could be made prettier with factory pattern (or the rust equivalent).
        let mut aggregators: Vec<Box<BamRecordAggregator>> = self.input
            .bam_readers
            .iter()
            .map(|bam_reader| {
                let mut samples = Vec::new();
                let text = String::from_utf8(Vec::from(bam_reader.header().as_bytes())).unwrap();
                for line in text.lines() {
                    if line.starts_with("@RG") {
                        match parse_line_rg(line.to_string()) {
                            Some((id, sm)) => {
                                samples.push(sm.clone());
                            }
                            None => (),
                        }
                    }
                }

                let result: Box<BamRecordAggregator> = match self.options.count_kind {
                    CountKind::Coverage => Box::new(CoverageAggregator::new(
                        &samples,
                        &self.input.read_group_to_sample,
                        &self.input.sample_to_index,
                        &self.options,
                        chrom_len,
                    )),
                    CountKind::Alignments => Box::new(CountAlignmentsAggregator::new(
                        &samples,
                        &self.input.read_group_to_sample,
                        &self.input.sample_to_index,
                        &self.options,
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
        }

        // Compute histogram of GC contents.
        // TODO: allow rounding via command line.
        let gc_step: f32 = 0.02;
        let mut gc_histo: HashMap<i32, i32> = HashMap::new();
        for gc in gc_content.iter() {
            let rounded_gc: i32 = (*gc / gc_step) as i32;
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
        let rid = self.output
            .bcf_writer
            .header()
            .name2rid(chrom.as_bytes())
            .unwrap() as i32;
        for (wid, gc) in gc_content.iter().enumerate() {
            let rounded_gc = (*gc / gc_step) as i32;

            let mut record = self.output.bcf_writer.empty_record();
            // Columns: CHROM, POS, ID, REF, ALT, FILTER
            record.inner_mut().rid = rid;
            let pos = wid * self.options.window_length as usize;
            record.set_pos(pos as i32);
            let window_end = min(
                end as usize,
                (wid + 1) * self.options.window_length as usize,
            ) as i32;
            record.update_id(
                format!("WIN_{}_{}_{}", chrom, pos + 1, window_end).as_bytes(),
            );
            record.update_alleles_str(b"N,<COUNT>");

            // INFO fields
            record.push_info_integer(b"END", &[window_end]).expect(
                "Could not write INFO/END",
            );
            record.push_info_float(b"GC", &[*gc]).expect(
                "Could not write INFO/GC",
            );
            record.push_info_integer(b"GAP", &[has_gap[wid] as i32]).expect(
                "Could not write INFO/GAP",
            );
            record
                .push_info_integer(b"GCWINDOWS", &[*gc_histo.get(&rounded_gc).unwrap_or(&0)])
                .expect("Could not write INFO/GCWINDOWS");

            // FORMAT/GT
            record
                .push_format_genotypes(
                    iter::repeat(bcf::gt_missing)
                        .take(self.input.sample_to_index.len())
                        .collect::<Vec<i32>>()
                        .as_slice(),
                )
                .expect("Could not push genotypes");

            // The remaining FORMAT fields are done by `BamRecordAggregator`s.
            for agg in &aggregators {
                agg.fill_bcf_record(&mut record, wid as u32);
            }

            // Actually write the record.
            self.output.bcf_writer.write(&record).expect(
                "Could not write record!",
            );
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

    /// Load mapability of `chrom` from BED file.
    fn load_mapability(&mut self, chrom: &str, chrom_len: u64) -> Result<Vec<f32>, String> {
        // Compute number of buckets and allocate array.
        let window_length = self.options.window_length as u64;
        let num_buckets = ((chrom_len + window_length - 1) / window_length) as usize;
        let mut result = vec![0 as f32; num_buckets];

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
            .unwrap() {
            Err(x) => panic!("Could not fetch in mapability BED file {:?}", x),
            _ => (),
        };

        // Load line by line.
        // let mut buf = Vec::new();
        // // TODO: can we make this loop tighter without unwrapping all the time?
        // while self.input
        //     .tbx_reader
        //     .as_mut()
        //     .map(|reader| reader.read(&mut buf).is_ok())
        //     .unwrap()
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
            let mapability = arr[3].parse::<f32>().unwrap();

            // Modify buckets in `result`.
            let mut window_id = begin as u64 / window_length;
            let window_begin = |window_id: u64| window_id * window_length;
            let window_end = |window_id: u64| window_begin(window_id) + window_length;
            while window_begin(window_id) < end {
                let len = min(window_end(window_id), end) - max(window_begin(window_id), begin);
                window_id += 1;
                result[window_id as usize] += (len as f32 / window_length as f32) * mapability;
            }
        }

        Ok(result)
    }

    /// Analyze reference for GC content, and is-gap (=N) status.
    fn analyze_reference(&mut self, chrom: &str) -> Result<(u64, Vec<f32>, Vec<bool>), String> {
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
