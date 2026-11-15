//! gloci command-line framework.
pub mod cli;
pub mod commands;
pub mod embed;
pub mod logging;

use cli::{Cli, Commands};
use embed::Outcome;
use gloci_io::{install_cleanup_handlers, shell_quote_args};
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

/// The stable core: the commands whose options and output change only with a version bump and a
/// `CHANGELOG.md` entry. Every other command is experimental and says so at the top of its help
/// ([`cli::EXPERIMENTAL`]). To promote a command, add it here and take the marker off its variant
/// in `cli/mod.rs`; `docs/book/src/stability.md` lists what promotion asks for.
pub const CORE_COMMANDS: &[&str] = &[
    "genomes",
    "index",
    "sort",
    "view",
    "intersect",
    "merge",
    "subtract",
    "closest",
    "coverage",
];

/// Is `name` (a top-level command) in the stable core?
pub fn is_core(name: &str) -> bool {
    CORE_COMMANDS.contains(&name)
}

/// The clap command tree (used by xtask for man pages and completions).
pub fn command() -> clap::Command {
    cli::command()
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

/// An error whose message was already printed: the exit status is 1 and nothing more is said.
/// `run` maps the other errors to an exit status by what they are ([`embed::classify`]): usage
/// errors exit 2, a broken pipe exits 0, and the rest print their message and exit 1.
#[derive(Debug)]
pub struct Reported;

impl std::fmt::Display for Reported {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("already reported")
    }
}

impl std::error::Error for Reported {}

/// Parse arguments, run the command and map the result to an exit code.
pub fn run(args: Vec<OsString>) -> ExitCode {
    let (code, message) = run_reporting(args);
    if let Some(m) = message {
        eprintln!("error: {m}");
    }
    code
}

/// [`run`], except that the final `error:` message is returned instead of printed.
///
/// The CLI is [`embed::parse`] (here clap's own error, so help, versions and colours print
/// as clap does), then the logger and the signal and panic handlers, then [`embed::execute`]
/// and [`embed::classify`].
#[doc(hidden)]
pub fn run_reporting(args: Vec<OsString>) -> (ExitCode, Option<String>) {
    let cli = match cli::try_parse_from(&args) {
        Ok(c) => c,
        Err(e) => {
            let code = e.exit_code();
            let _ = e.print();
            return (ExitCode::from(code as u8), None);
        }
    };
    logging::init(cli.quiet, cli.verbose);
    install_cleanup_handlers();
    let opts = embed::ExecOptions {
        command_line: shell_quote_args(&strip_thread_options(&args)),
        tmp_dir: None,
    };
    match embed::execute(&cli, &opts) {
        Ok(true) => (ExitCode::SUCCESS, None),
        Ok(false) => (ExitCode::from(1), None),
        Err(e) => match embed::classify(&e) {
            Outcome::Usage(msg) => (ExitCode::from(2), Some(msg)),
            Outcome::BrokenPipe => (ExitCode::SUCCESS, None),
            Outcome::Cancelled => (ExitCode::from(1), Some("interrupted".into())),
            Outcome::Reported => (ExitCode::from(1), None),
            Outcome::Failure(msg) => (ExitCode::from(1), Some(msg)),
        },
    }
}

/// Run the parsed command: its result, and false when `check` found a problem.
pub(crate) fn dispatch(cli: &Cli, ctx: &Ctx) -> (anyhow::Result<()>, bool) {
    let mut all_ok = true;
    let result = match &cli.command {
        Commands::Genomes(a) => commands::genomes::run(a),
        Commands::Index(a) => commands::index::run(ctx, a),
        Commands::Sort(a) => commands::sort::run(ctx, a),
        Commands::View(a) => commands::view::run(ctx, a),
        Commands::Cat(a) => commands::cat::run(ctx, a),
        Commands::Check(a) => commands::check::run(ctx, a).map(|ok| all_ok = ok),
        Commands::Intersect(a) => commands::intersect::run(ctx, a),
        Commands::Merge(a) => commands::merge::run(ctx, a),
        Commands::Cluster(a) => commands::cluster::run(ctx, a),
        Commands::Spacing(a) => commands::spacing::run(ctx, a),
        Commands::Complement(a) => commands::complement::run(ctx, a),
        Commands::Partition(a) => commands::partition::run(ctx, a),
        Commands::Genomecov(a) => commands::genomecov::run(ctx, a),
        Commands::Slop(a) => commands::slop::run(ctx, a),
        Commands::Flank(a) => commands::flank::run(ctx, a),
        Commands::Shift(a) => commands::shift::run(ctx, a),
        Commands::Windows(a) => commands::windows::run(ctx, a),
        Commands::Convert(a) => commands::convert::run(ctx, a),
        Commands::Getfasta(a) => commands::getfasta::run(ctx, a),
        Commands::Nuc(a) => commands::nuc::run(ctx, a),
        Commands::Subtract(a) => commands::subtract::run(ctx, a),
        Commands::Closest(a) => commands::closest::run(ctx, a),
        Commands::Map(a) => commands::map::run(ctx, a),
        Commands::Coverage(a) => commands::coverage::run(ctx, a),
        Commands::Groupby(a) => commands::groupby::run(ctx, a),
        Commands::Maskfasta(a) => commands::maskfasta::run(ctx, a),
        Commands::Shuffle(a) => commands::shuffle::run(ctx, a),
        Commands::Sample(a) => commands::sample::run(ctx, a),
        Commands::Random(a) => commands::random::run(ctx, a),
        Commands::Stats(a) => match &a.command {
            cli::StatsCommand::Jaccard(j) => commands::stats::run_jaccard(ctx, j),
            cli::StatsCommand::Fisher(f) => commands::stats::run_fisher(ctx, f),
            cli::StatsCommand::Reldist(r) => commands::stats::run_reldist(ctx, r),
        },
        Commands::PairIntersect(a) => commands::pairs::run_intersect(ctx, a),
        Commands::PairPair(a) => commands::pairs::run_pair(ctx, a),
        Commands::Dict(a) => commands::dict::run(ctx, a),
        Commands::Flags(a) => commands::flags::run(a),
        Commands::Idxstats(a) => commands::idxstats::run(ctx, a),
        Commands::Flagstat(a) => commands::flagstat::run(ctx, a),
        Commands::Depth(a) => commands::depth::run(ctx, a),
        Commands::Bedcov(a) => commands::bedcov::run(ctx, a),
        Commands::Reheader(a) => commands::reheader::run(ctx, a),
        Commands::Igv(a) => commands::igv::run(ctx, a),
        Commands::Links(a) => commands::links::run(ctx, a),
        Commands::Tag(a) => commands::tag::run(ctx, a),
        Commands::Faidx(a) => commands::faidx::run(ctx, a),
        Commands::Head(a) => commands::head::run(ctx, a),
        Commands::Samples(a) => commands::samples::run(ctx, a),
        Commands::Vcf(a) => match &a.command {
            cli::VcfCommand::Filter(f) => commands::vcf::run_filter(ctx, f),
            cli::VcfCommand::Extract(e) => commands::vcf::run_extract(ctx, e),
            cli::VcfCommand::Stats(e) => commands::vcf::run_stats(ctx, e),
            cli::VcfCommand::Concat(c) => commands::vcf::run_concat(ctx, c),
            cli::VcfCommand::Split(c) => commands::vcf_split::run(ctx, c),
            cli::VcfCommand::Norm(c) => commands::vcf_norm::run(ctx, c),
            cli::VcfCommand::Annotate(c) => commands::vcf_annotate::run(ctx, c),
            cli::VcfCommand::Merge(c) => commands::vcf_merge::run(ctx, c),
            cli::VcfCommand::Compare(c) => commands::vcf::run_compare(ctx, c),
            cli::VcfCommand::Popgen(c) => commands::vcf::run_popgen(ctx, c),
        },
    };
    (result, all_ok)
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
