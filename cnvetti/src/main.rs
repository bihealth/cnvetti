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
extern crate lib_genotype;
extern crate lib_merge_cov;
extern crate lib_mod_cov;
extern crate lib_model_pool;
extern crate lib_model_wis;
extern crate lib_normalize;
extern crate lib_segment;
extern crate lib_visualize;

mod quick_pool_build_model;
mod quick_pool_call;
mod quick_wis_build_model;
mod quick_wis_call;

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
        // cnvetti cmd <coverage|normalize|...>
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
            ("segment", Some(m)) => {
                lib_segment::run(&mut logger, &lib_segment::SegmentOptions::new(&m))
                    .chain_err(|| "Could not execute 'cmd segment'")?
            }
            ("de-bias", Some(_m)) => bail!("cmd de-bias not implemented!"),
            ("build-model-pool", Some(m)) => {
                lib_model_pool::run(&mut logger, &lib_model_pool::BuildModelPoolOptions::new(&m))
                    .chain_err(|| "Could not execute 'cmd build-model-pool'")?
            }
            ("build-model-wis", Some(m)) => {
                lib_model_wis::run(&mut logger, &lib_model_wis::BuildModelWisOptions::new(&m))
                    .chain_err(|| "Could not execute 'cmd build-model-wis'")?
            }
            ("mod-coverage", Some(m)) => lib_mod_cov::run(
                &mut logger,
                &lib_mod_cov::ModelBasedCoverageOptions::new(&m),
            ).chain_err(|| "Could not execute 'cmd mod-coverage'")?,
            ("discover", Some(_m)) => bail!("cmd discover not implemented!"),
            ("genotype", Some(_m)) => {
                lib_genotype::run(&mut logger, &lib_genotype::GenotypeOptions::new(&m))
                    .chain_err(|| "Could not execute 'cmd genotype'")?
            }
            _ => bail!("Invalid command: {}", m.subcommand().0),
        },
        // cnvetti quick <wis-build-model|wis-call|pool-build-model|pool-call>
        ("quick", Some(m)) => match m.subcommand() {
            ("pool-build-model", Some(m)) => {
                quick_pool_build_model::run(
                    &mut logger,
                    &quick_pool_build_model::QuickPoolBuildModelOptions::new(&m),
                ).chain_err(|| "Could not execute 'cmd quick pool-build-models")?
            }
            ("pool-call", Some(m)) => {
                quick_pool_call::run(&mut logger, &quick_pool_call::QuickPoolCallOptions::new(&m))
                    .chain_err(|| "Could not execute 'cmd quick pool-call")?
            }
            ("wis-build-model", Some(m)) => {
                quick_wis_build_model::run(
                    &mut logger,
                    &quick_wis_build_model::QuickWisBuildModelOptions::new(&m),
                ).chain_err(|| "Could not execute 'cmd quick wis-build-models")?
            }
            ("wis-call", Some(m)) => {
                quick_wis_call::run(&mut logger, &quick_wis_call::QuickWisCallOptions::new(&m))
                    .chain_err(|| "Could not execute 'cmd quick wis-call")?
            }
            _ => bail!("Invalid command: {}", m.subcommand().0),
        },
        // cnvetti visualize <cov-to-igv>
        ("visualize", Some(m)) => match m.subcommand() {
            ("cov-to-igv", Some(m)) => lib_visualize::cov_to_igv::run(
                &mut logger,
                &lib_visualize::CovToIgvOptions::new(&m),
            ).chain_err(|| "Could not execute 'visualize cov-to-igv'")?,
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
