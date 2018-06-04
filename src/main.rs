extern crate cnvetti;

include!(concat!(env!("OUT_DIR"), "/version.rs"));

#[macro_use]
extern crate clap;

#[macro_use]
extern crate slog;
extern crate slog_async;
extern crate slog_term;

use slog::Drain;

use std::result;
use std::sync::atomic::Ordering;
use std::sync::{atomic, Arc};

use clap::{App, ArgMatches};
use std::process;

use cnvetti::cli::call;
use cnvetti::cli::cohort_stats;
use cnvetti::cli::coverage;
use cnvetti::cli::normalize;
use cnvetti::cli::segment;
use cnvetti::cli::shortcuts;
use cnvetti::cli::wisexome;

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
            0 => slog::Level::Warning,
            1 => slog::Level::Info,
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

    let mut logger = slog::Logger::root(drain, o!());

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
        ("coverage", Some(m)) => coverage::call(&mut logger, &coverage::Options::new(&m)),
        ("normalize", Some(m)) => normalize::call(&mut logger, &normalize::Options::new(&m)),
        ("segment", Some(m)) => segment::call(&mut logger, &segment::Options::new(&m)),
        ("cohort-stats", Some(m)) => {
            cohort_stats::call(&mut logger, &cohort_stats::Options::new(&m))
        }
        ("call", Some(m)) => call::call(&mut logger, &call::Options::new(&m)),
        ("quick", Some(m)) => match m.subcommand() {
            ("wgs-deep", Some(m)) => {
                shortcuts::call_quick_wgs_deep(&mut logger, &shortcuts::WgsDeepOptions::new(&m))
            }
            ("wgs-cov-bins", Some(m)) => shortcuts::call_quick_wgs_cov_bins(
                &mut logger,
                &shortcuts::WgsCovBinsOptions::new(&m),
            ),
            _ => Err("Invalid command".to_string()),
        },
        ("wisexome", Some(m)) => match m.subcommand() {
            ("count", Some(m)) => {
                wisexome::call_wise_count(&mut logger, &wisexome::CountOptions::new(&m))
            }
            ("normalize", Some(m)) => {
                wisexome::call_wise_normalize(&mut logger, &wisexome::NormalizeOptions::new(&m))
            }
            ("build-ref", Some(m)) => {
                wisexome::call_wise_build_ref(&mut logger, &wisexome::BuildRefOptions::new(&m))
            }
            _ => Err("Invalid command".to_string()),
        },
        _ => Err("Invalid command".to_string()),
    }
}

fn main() {
    let yaml = load_yaml!("cli.yaml");
    // TODO: the VERSION from below is ignored :(
    let matches = App::from_yaml(yaml).version(VERSION).get_matches();

    if let Err(e) = run(matches) {
        println!("Application error: {}", e);
        process::exit(1);
    }
}
