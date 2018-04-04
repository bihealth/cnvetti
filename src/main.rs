extern crate cnvetti;

include!(concat!(env!("OUT_DIR"), "/version.rs"));

extern crate clap;

#[macro_use]
extern crate slog;
extern crate slog_async;
extern crate slog_term;

use slog::Drain;

use std::sync::{atomic, Arc};
use std::sync::atomic::Ordering;
use std::result;

use std::process;
use clap::{App, Arg, ArgMatches, SubCommand};

use cnvetti::cli::coverage;

/// Custom Drain logic
struct RuntimeLevelFilter<D> {
    drain: D,
    log_level: Arc<atomic::AtomicIsize>,
}

impl<D> Drain for RuntimeLevelFilter<D>
where
    D: Drain,
{
    type Ok = Option<D::Ok>;
    type Err = Option<D::Err>;

    fn log(
        &self,
        record: &slog::Record,
        values: &slog::OwnedKVList,
    ) -> result::Result<Self::Ok, Self::Err> {
        let current_level = match self.log_level.load(Ordering::Relaxed) {
            0 => slog::Level::Error,
            1 => slog::Level::Warning,
            2 => slog::Level::Info,
            _ => slog::Level::Trace,
        };

        if record.level().is_at_least(current_level) {
            self.drain.log(record, values).map(Some).map_err(Some)
        } else {
            Ok(None)
        }
    }
}

fn run(matches: ArgMatches) -> Result<(), String> {
    // Logging setup ------------------------------------------------------------------------------

    // Atomic variable controlling logging level
    let log_level = Arc::new(atomic::AtomicIsize::new(1));

    // Perform slog setup
    let decorator = slog_term::TermDecorator::new().build();
    let drain = slog_term::FullFormat::new(decorator).build();
    let drain = RuntimeLevelFilter {
        drain: drain,
        log_level: log_level.clone(),
    }.fuse();
    let drain = slog_async::Async::new(drain).build().fuse();

    let logger = slog::Logger::root(drain, o!());

    // Switch log level
    if matches.is_present("quiet") {
        log_level.store(0, Ordering::Relaxed);
    } else {
        log_level.store(
            1 + matches.occurrences_of("verbose") as isize,
            Ordering::Relaxed,
        );
    }

    match matches.subcommand() {
        ("coverage", Some(m)) => {
            cnvetti::cli::coverage::call(logger, cnvetti::cli::coverage::build_options(&m))
        }
        _ => Err("Invalid command".to_string()),
    }
}

fn main() {
    let matches = App::new("cnvetti")
        .version(VERSION)
        .author("Manuel Holtgrewe")
        .about("CNV calling for the masses")
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .multiple(true)
                .help("verbosity level"),
        )
        .arg(
            Arg::with_name("quiet")
                .short("q")
                .help("Only show error messages"),
        )
        .subcommand(
            SubCommand::with_name("coverage")
                .about("Collect per-sample coverage information")
                .arg(
                    Arg::with_name("reference")
                        .long("reference")
                        .short("r")
                        .required(true)
                        .value_name("REFERENCE.fa")
                        .help("Path to FAI-indexed FASTA file"),
                )
                .arg(
                    Arg::with_name("input")
                        .long("input")
                        .short("i")
                        .required(true)
                        .value_name("INPUT.bam")
                        .help("Path to BAI-indexed BAM file"),
                )
                .arg(
                    Arg::with_name("output")
                        .long("output")
                        .short("o")
                        .required(true)
                        .value_name("OUTPUT.{bcf,vcf}")
                        .help("Path to output VCF/BCF file"),
                )
                .arg(
                    Arg::with_name("genome-region")
                        .long("genome-region")
                        .value_name("REGION")
                        .multiple(true)
                        .help("Genome region(s) (e.g., chr1:1,000,000-1,100,000) to process"),
                )
                .arg(
                    Arg::with_name("mapability-bed")
                        .long("mapability-bed")
                        .value_name("MAPABILITY.bed")
                        .help("Path to mapability BED file"),
                ),
        )
        .get_matches();

    if let Err(e) = run(matches) {
        println!("Application error: {}", e);
        process::exit(1);
    }
}
