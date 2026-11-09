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

/// The arguments without the thread option, so the recorded `##glociCommand`
/// does not depend on it. Uses clap's option metadata, so every spelling clap
/// accepts is handled: `-t N`, `-tN`, `-t=N`, clusters such as `-qt2` or
/// `-vt=2`, `--threads N` and `--threads=N`. Values of other options (even if
/// they look like `-t`) and everything after `--` are kept verbatim.
fn strip_thread_options(args: &[OsString]) -> Vec<OsString> {
    let mut root = command();
    root.build();
    let takes_value = |a: &clap::Arg| a.get_action().takes_values();
    let is_threads = |a: &clap::Arg| a.get_long() == Some("threads");
    let mut cmd: &clap::Command = &root;
    let mut out = Vec::with_capacity(args.len());
    let mut i = 0;
    if let Some(first) = args.first() {
        out.push(first.clone());
        i = 1;
    }
    while i < args.len() {
        let a = &args[i];
        i += 1;
        let Some(s) = a.to_str() else {
            out.push(a.clone());
            continue;
        };
        if s == "--" {
            out.push(a.clone());
            out.extend(args[i..].iter().cloned());
            break;
        }
        if let Some(long) = s.strip_prefix("--") {
            let (name, has_eq) = match long.split_once('=') {
                Some((n, _)) => (n, true),
                None => (long, false),
            };
            let arg = cmd.get_arguments().find(|x| x.get_long() == Some(name));
            let separate = !has_eq && arg.is_some_and(takes_value);
            if arg.is_some_and(is_threads) {
                if separate {
                    i += 1;
                }
                continue;
            }
            out.push(a.clone());
            if separate && i < args.len() {
                out.push(args[i].clone());
                i += 1;
            }
            continue;
        }
        if s.len() > 1 && s.starts_with('-') {
            let chars: Vec<char> = s[1..].chars().collect();
            let mut kept = String::from("-");
            let mut next_is_value = false;
            let mut drop_next = false;
            let mut k = 0;
            while k < chars.len() {
                let c = chars[k];
                k += 1;
                let arg = cmd.get_arguments().find(|x| x.get_short() == Some(c));
                let Some(arg) = arg else {
                    // unknown to clap: keep the rest verbatim
                    kept.push(c);
                    kept.extend(&chars[k..]);
                    break;
                };
                if takes_value(arg) {
                    let rest: String = chars[k..].iter().collect();
                    if is_threads(arg) {
                        drop_next = rest.is_empty();
                    } else {
                        kept.push(c);
                        kept.push_str(&rest);
                        next_is_value = rest.is_empty();
                    }
                    break;
                }
                kept.push(c);
            }
            if kept.len() > 1 {
                out.push(OsString::from(kept));
            }
            if (drop_next || next_is_value) && i < args.len() {
                if next_is_value {
                    out.push(args[i].clone());
                }
                i += 1;
            }
            continue;
        }
        // a bare word before the subcommand is the subcommand
        if let Some(sub) = cmd.find_subcommand(s) {
            cmd = sub;
        }
        out.push(a.clone());
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
        Commands::Genomes(_) => (None, None),
        Commands::Index(a) => (a.threads.threads, None),
        Commands::Sort(a) => (a.threads.threads, a.tmp_dir.tmp_dir.clone()),
        Commands::View(a) => (a.threads.threads, None),
        Commands::Check(a) => (a.threads.threads, None),
        Commands::Intersect(a) => (a.threads.threads, a.sort.tmp_dir.clone()),
        Commands::Merge(a) => (a.threads.threads, a.sort.tmp_dir.clone()),
        Commands::Complement(a) => (a.threads.threads, a.sort.tmp_dir.clone()),
        Commands::Partition(a) => (a.threads.threads, a.sort.tmp_dir.clone()),
        Commands::Slop(a) => (a.threads.threads, None),
        Commands::Flank(a) => (a.threads.threads, None),
        Commands::Windows(a) => (a.threads.threads, None),
        Commands::Convert(a) => (a.threads.threads, None),
        Commands::Subtract(a) => (a.threads.threads, a.sort.tmp_dir.clone()),
        Commands::Closest(a) => (a.threads.threads, a.sort.tmp_dir.clone()),
        Commands::Map(a) => (a.threads.threads, a.sort.tmp_dir.clone()),
        Commands::Coverage(a) => (a.threads.threads, a.sort.tmp_dir.clone()),
        Commands::Groupby(a) => (a.threads.threads, None),
        Commands::Vcf(a) => match &a.command {
            cli::VcfCommand::Filter(f) => (f.threads.threads, None),
        },
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
        Commands::Genomes(a) => commands::genomes::run(a),
        Commands::Index(a) => commands::index::run(&ctx, a),
        Commands::Sort(a) => commands::sort::run(&ctx, a),
        Commands::View(a) => commands::view::run(&ctx, a),
        Commands::Check(a) => commands::check::run(&ctx, a).map(|ok| all_ok = ok),
        Commands::Intersect(a) => commands::intersect::run(&ctx, a),
        Commands::Merge(a) => commands::merge::run(&ctx, a),
        Commands::Complement(a) => commands::complement::run(&ctx, a),
        Commands::Partition(a) => commands::partition::run(&ctx, a),
        Commands::Slop(a) => commands::slop::run(&ctx, a),
        Commands::Flank(a) => commands::flank::run(&ctx, a),
        Commands::Windows(a) => commands::windows::run(&ctx, a),
        Commands::Convert(a) => commands::convert::run(&ctx, a),
        Commands::Subtract(a) => commands::subtract::run(&ctx, a),
        Commands::Closest(a) => commands::closest::run(&ctx, a),
        Commands::Map(a) => commands::map::run(&ctx, a),
        Commands::Coverage(a) => commands::coverage::run(&ctx, a),
        Commands::Groupby(a) => commands::groupby::run(&ctx, a),
        Commands::Vcf(a) => match &a.command {
            cli::VcfCommand::Filter(f) => commands::vcf::run_filter(&ctx, f),
        },
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
        let want = ["gloci", "sort", "a.bed"];
        for t in [
            &["-t", "2"][..],
            &["-t2"],
            &["-t=2"],
            &["--threads", "2"],
            &["--threads=2"],
        ] {
            let mut v = vec!["gloci", "sort", "a.bed"];
            v.extend_from_slice(t);
            assert_eq!(strip(&v), want, "{t:?}");
            let mut v = vec!["gloci", "sort"];
            v.extend_from_slice(t);
            v.push("a.bed");
            assert_eq!(strip(&v), want, "{t:?} first");
        }
        // clusters, keeping the other flags
        assert_eq!(
            strip(&["gloci", "sort", "a.bed", "-qt2"]),
            ["gloci", "sort", "a.bed", "-q"]
        );
        assert_eq!(
            strip(&["gloci", "sort", "a.bed", "-qt", "2"]),
            ["gloci", "sort", "a.bed", "-q"]
        );
        assert_eq!(
            strip(&["gloci", "sort", "a.bed", "-vt=2"]),
            ["gloci", "sort", "a.bed", "-v"]
        );
        assert_eq!(
            strip(&["gloci", "sort", "a.bed", "-qvt2"]),
            ["gloci", "sort", "a.bed", "-qv"]
        );
        assert_eq!(
            strip(&["gloci", "-qv", "sort", "a.bed", "-t2"]),
            ["gloci", "-qv", "sort", "a.bed"]
        );
        // values that merely look like the option are kept
        assert_eq!(
            strip(&["gloci", "sort", "a.bed", "-o", "-t", "-t", "3"]),
            ["gloci", "sort", "a.bed", "-o", "-t"]
        );
        assert_eq!(
            strip(&["gloci", "sort", "a.bed", "--tmp-dir", "-t2", "-t2"]),
            ["gloci", "sort", "a.bed", "--tmp-dir", "-t2"]
        );
        assert_eq!(
            strip(&["gloci", "sort", "-o", "t.vcf", "--tmp-dir=t", "-q"]),
            ["gloci", "sort", "-o", "t.vcf", "--tmp-dir=t", "-q"]
        );
        // after `--` everything is verbatim
        assert_eq!(
            strip(&["gloci", "view", "-t", "2", "--", "-t2", "x"]),
            ["gloci", "view", "--", "-t2", "x"]
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
