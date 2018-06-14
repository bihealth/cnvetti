// `error_chain!` can recurse deeply.
#![recursion_limit = "1024"]

// We are using `error-chain`.
#[macro_use]
extern crate error_chain;

// We are using the `clap` crate for command line argument parsing.
#[macro_use]
extern crate clap;

// We are using the `slog` crate for logging.
#[macro_use]
extern crate slog;
extern crate slog_async;
extern crate slog_term;

use std::result;
use std::sync::atomic::Ordering;
use std::sync::{atomic, Arc};

use clap::{App, ArgMatches};

use slog::Drain;

mod errors {
    // Create the Error, ErrorKind, ResultExt, and Result types
    error_chain!{}
}

pub use errors::*;

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

extern crate lib_coverage;

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
