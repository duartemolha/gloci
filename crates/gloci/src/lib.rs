//! gloci command-line framework.
pub mod cli;
pub mod commands;
pub mod logging;

use clap::{CommandFactory, Parser};
use cli::{Cli, Commands};
use gloci_io::{IoError, UsageError, install_cleanup_handlers, shell_quote_args};
use std::ffi::OsString;
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::process::ExitCode;

/// Per-invocation context handed to every command.
pub struct Ctx {
    pub threads: NonZeroUsize,
    pub tmp_dir: PathBuf,
    pub command_line: String,
    pub version: &'static str,
}

/// The spill directory: the command's `--tmp-dir` (which clap also fills from
/// `GLOCI_TMPDIR`), else a non-empty `GLOCI_TMPDIR`, else the system temp dir.
fn resolve_tmp_dir(arg: Option<PathBuf>) -> PathBuf {
    let non_empty = |p: &PathBuf| !p.as_os_str().is_empty();
    arg.filter(non_empty)
        .or_else(|| {
            std::env::var_os("GLOCI_TMPDIR")
                .map(PathBuf::from)
                .filter(non_empty)
        })
        .unwrap_or_else(std::env::temp_dir)
}

/// The clap command tree (used by xtask for man pages and completions).
pub fn command() -> clap::Command {
    Cli::command()
}

/// What an error means for the exit status.
enum Outcome {
    Usage(String),
    BrokenPipe,
    Failure(String),
}

fn classify(err: &anyhow::Error) -> Outcome {
    for cause in err.chain() {
        if let Some(u) = cause.downcast_ref::<UsageError>() {
            return Outcome::Usage(u.0.clone());
        }
        if let Some(IoError::Usage(u)) = cause.downcast_ref::<IoError>() {
            return Outcome::Usage(u.0.clone());
        }
    }
    for cause in err.chain() {
        if let Some(e) = cause.downcast_ref::<std::io::Error>()
            && e.kind() == std::io::ErrorKind::BrokenPipe
        {
            return Outcome::BrokenPipe;
        }
    }
    // IoError messages already embed their source and hint; `{:#}` would repeat it.
    if err.downcast_ref::<IoError>().is_some() {
        Outcome::Failure(err.to_string())
    } else {
        Outcome::Failure(format!("{err:#}"))
    }
}

/// Parse arguments, run the command and map the result to an exit code.
pub fn run(args: Vec<OsString>) -> ExitCode {
    let cli = match Cli::try_parse_from(&args) {
        Ok(c) => c,
        Err(e) => {
            let code = e.exit_code();
            let _ = e.print();
            return ExitCode::from(code as u8);
        }
    };
    logging::init(cli.quiet, cli.verbose);
    install_cleanup_handlers();

    let (requested, tmp_arg) = match &cli.command {
        Commands::Index(a) => (a.threads.threads, None),
        Commands::Sort(a) => (a.threads.threads, a.tmp_dir.tmp_dir.clone()),
        Commands::View(a) => (a.threads.threads, None),
    };
    let threads = requested
        .unwrap_or_else(|| std::thread::available_parallelism().unwrap_or(NonZeroUsize::MIN));
    // An already-built global pool is fine.
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.get())
        .build_global();

    let ctx = Ctx {
        threads,
        tmp_dir: resolve_tmp_dir(tmp_arg),
        command_line: shell_quote_args(&args),
        version: env!("CARGO_PKG_VERSION"),
    };

    let result = match &cli.command {
        Commands::Index(a) => commands::index::run(&ctx, a),
        Commands::Sort(a) => commands::sort::run(&ctx, a),
        Commands::View(a) => commands::view::run(&ctx, a),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => match classify(&e) {
            Outcome::Usage(msg) => {
                eprintln!("error: {msg}");
                ExitCode::from(2)
            }
            Outcome::BrokenPipe => ExitCode::SUCCESS,
            Outcome::Failure(msg) => {
                eprintln!("error: {msg}");
                ExitCode::from(1)
            }
        },
    }
}
