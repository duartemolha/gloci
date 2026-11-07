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

/// The arguments without thread options (`-t N`, `-tN`, `--threads N`,
/// `--threads=N`), so the recorded `##glociCommand` does not depend on them.
/// Everything after `--` is kept verbatim.
fn strip_thread_options(args: &[OsString]) -> Vec<OsString> {
    let mut out = Vec::with_capacity(args.len());
    let mut it = args.iter().enumerate();
    while let Some((i, a)) = it.next() {
        // argv[0] is the program name
        if i == 0 {
            out.push(a.clone());
            continue;
        }
        let Some(s) = a.to_str() else {
            out.push(a.clone());
            continue;
        };
        match s {
            "--" => {
                out.extend(args[i..].iter().cloned());
                break;
            }
            "-t" | "--threads" => {
                it.next();
            }
            _ if s.starts_with("--threads=") => {}
            _ if s.len() > 2
                && s.starts_with("-t")
                && s[2..].bytes().all(|b| b.is_ascii_digit()) => {}
            _ => out.push(a.clone()),
        }
    }
    out
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
        Commands::Check(a) => (a.threads.threads, None),
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
        command_line: shell_quote_args(&strip_thread_options(&args)),
        version: env!("CARGO_PKG_VERSION"),
    };

    let mut all_ok = true;
    let result = match &cli.command {
        Commands::Index(a) => commands::index::run(&ctx, a),
        Commands::Sort(a) => commands::sort::run(&ctx, a),
        Commands::View(a) => commands::view::run(&ctx, a),
        Commands::Check(a) => commands::check::run(&ctx, a).map(|ok| all_ok = ok),
    };
    match result {
        Ok(()) if all_ok => ExitCode::SUCCESS,
        Ok(()) => ExitCode::from(1),
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Error, ErrorKind};

    fn strip(a: &[&str]) -> Vec<String> {
        let v: Vec<OsString> = a.iter().map(OsString::from).collect();
        strip_thread_options(&v)
            .into_iter()
            .map(|s| s.into_string().unwrap())
            .collect()
    }

    #[test]
    fn thread_options_are_stripped_from_command_line() {
        assert_eq!(
            strip(&["gloci", "sort", "-t", "4", "a.bed"]),
            ["gloci", "sort", "a.bed"]
        );
        assert_eq!(
            strip(&["gloci", "sort", "-t4", "a.bed"]),
            ["gloci", "sort", "a.bed"]
        );
        assert_eq!(
            strip(&["gloci", "sort", "--threads", "4", "a.bed"]),
            ["gloci", "sort", "a.bed"]
        );
        assert_eq!(
            strip(&["gloci", "sort", "--threads=4", "a.bed"]),
            ["gloci", "sort", "a.bed"]
        );
        // others untouched
        assert_eq!(
            strip(&[
                "gloci",
                "sort",
                "-o",
                "t.vcf",
                "--tmp-dir",
                "t",
                "-tx",
                "-q"
            ]),
            [
                "gloci",
                "sort",
                "-o",
                "t.vcf",
                "--tmp-dir",
                "t",
                "-tx",
                "-q"
            ]
        );
        // after `--` everything is verbatim
        assert_eq!(
            strip(&["gloci", "view", "-t", "2", "--", "-t", "x"]),
            ["gloci", "view", "--", "-t", "x"]
        );
    }

    fn pipe() -> Error {
        Error::from(ErrorKind::BrokenPipe)
    }

    /// An error whose `source()` is a `UsageError`.
    #[derive(Debug)]
    struct HasUsage(UsageError);
    impl std::fmt::Display for HasUsage {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "has usage")
        }
    }
    impl std::error::Error for HasUsage {
        fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
            Some(&self.0)
        }
    }

    #[test]
    fn usage_error_beats_broken_pipe_in_same_chain() {
        // chain: io::Error(BrokenPipe) -> UsageError (a custom io::Error
        // forwards `source()` of its payload), so both are in one chain.
        let io = Error::new(
            ErrorKind::BrokenPipe,
            HasUsage(UsageError("bad flag".into())),
        );
        let e = anyhow::Error::new(io);
        let kinds: Vec<_> = e
            .chain()
            .map(|c| c.downcast_ref::<UsageError>().is_some())
            .collect();
        assert!(kinds.contains(&true), "fixture must hold both: {kinds:?}");
        assert!(matches!(classify(&e), Outcome::Usage(m) if m == "bad flag"));
        let wrapped = anyhow::Error::new(IoError::Usage(UsageError("u".into()))).context("ctx");
        assert!(matches!(classify(&wrapped), Outcome::Usage(_)));
    }

    #[test]
    fn write_broken_pipe_is_quiet_success() {
        let e = anyhow::Error::new(IoError::Write {
            path: "out".into(),
            source: pipe(),
        });
        assert!(matches!(classify(&e), Outcome::BrokenPipe));
    }

    #[test]
    fn bare_broken_pipe_is_quiet_success() {
        assert!(matches!(
            classify(&anyhow::Error::new(pipe())),
            Outcome::BrokenPipe
        ));
        let other = anyhow::Error::new(Error::from(ErrorKind::NotFound));
        assert!(matches!(classify(&other), Outcome::Failure(_)));
    }
}
