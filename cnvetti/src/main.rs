// `error_chain!` can recurse deeply.
#![recursion_limit = "1024"]

// We are using `error-chain`.
#[macro_use]
extern crate error_chain;

// We are using the `clap` crate for command line argument parsing.
#[macro_use]
extern crate clap;

extern crate rayon;

// We are using the `slog` crate for logging.
#[macro_use]
extern crate slog;
extern crate slog_async;
extern crate slog_term;

extern crate tempdir;

use slog::Drain;

use std::result;
use std::sync::atomic::Ordering;
use std::sync::{atomic, Arc};

use clap::{App, ArgMatches};

mod errors {
    // Create the Error, ErrorKind, ResultExt, and Result types
    error_chain!{}
}

pub use errors::*;

extern crate lib_coverage;
extern crate lib_merge_cov;
extern crate lib_model_wis;
extern crate lib_normalize;

mod quick_wis_build_model;

/// Custom `slog` Drain logic
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

fn run(matches: ArgMatches) -> Result<()> {
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
    };

    // Dispatch commands from command line.
    match matches.subcommand() {
        ("cmd", Some(m)) => match m.subcommand() {
            ("coverage", Some(m)) => {
                lib_coverage::run(&mut logger, &lib_coverage::CoverageOptions::new(&m))
                    .chain_err(|| "Could not execute 'cmd coverage'")?
            }
            ("normalize", Some(m)) => {
                lib_normalize::run(&mut logger, &lib_normalize::NormalizeOptions::new(&m))
                    .chain_err(|| "Could not execute 'cmd normalize'")?
            }
            ("filter", Some(_m)) => bail!("cmd filter not implemented!"),
            ("merge-cov", Some(m)) => {
                lib_merge_cov::run(&mut logger, &lib_merge_cov::MergeCovOptions::new(&m))
                    .chain_err(|| "Could not execute 'cmd merge-cov'")?
            }
            ("de-bias", Some(_m)) => bail!("cmd de-bias not implemented!"),
            ("build-model-pool", Some(_m)) => bail!("cmd build-model-pool not implemented!"),
            ("build-model-wis", Some(m)) => {
                lib_model_wis::run(&mut logger, &lib_model_wis::BuildModelWisOptions::new(&m))
                    .chain_err(|| "Could not execute 'cmd build-model-wis'")?
            }
            ("mod-coverage", Some(_m)) => bail!("cmd mod-coverage not implemented!"),
            ("discover", Some(_m)) => bail!("cmd discover not implemented!"),
            ("genotype", Some(_m)) => bail!("cmd genotype not implemented!"),
            _ => bail!("Invalid command: {}", m.subcommand().0),
        },
        ("quick", Some(m)) => match m.subcommand() {
            ("wis-build-model", Some(m)) => {
                quick_wis_build_model::run(
                    &mut logger,
                    &quick_wis_build_model::QuickWisBuildModelOptions::new(&m),
                ).chain_err(|| "Could not execute 'cmd quick wis-build-models")?
            }
            _ => bail!("Invalid command: {}", m.subcommand().0),
        },
        _ => bail!("Invalid command: {}", matches.subcommand().0),
    }

    Ok(())
}

fn main() {
    let yaml = load_yaml!("cli.yaml");
    let matches = App::from_yaml(yaml).get_matches();

    if let Err(ref e) = run(matches) {
        eprintln!("error: {}", e);

        for e in e.iter().skip(1) {
            eprintln!("caused by: {}", e);
        }

        // The backtrace is not always generated. Try to run this example
        // with `RUST_BACKTRACE=1`.
        if let Some(backtrace) = e.backtrace() {
            eprintln!("backtrace: {:?}", backtrace);
        }

        ::std::process::exit(1);
    }
}
