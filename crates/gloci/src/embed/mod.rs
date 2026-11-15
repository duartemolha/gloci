//! Running gloci commands in-process, for a host such as the Python library (D83). The host
//! builds an argv (from keyword arguments, [`kwargs_to_argv`]), [`parse`]s it with the CLI's
//! own clap definition, and [`execute`]s it through the same dispatch as the command line, so
//! the defaults, validation, messages and output bytes are the CLI's.
//!
//! A step's whole argv, with its inputs where its binding puts them and the terminal's options
//! (`-o`, `-t`) before any `--`, is [`step_argv`]; nothing may be appended to it.
//!
//! What a host must do, since `execute` does none of it:
//! - **Cancellation (D91).** Make one [`RunToken`] per run and enter it on the thread that
//!   calls `execute` ([`RunGuard::enter`]); `execute` installs it on every thread of the run
//!   (the command's thread and its rayon pool). [`RunToken::cancel`] then makes the run fail
//!   with a cancellation, removing partial output files. Judge the error under the same token:
//!   [`classify`] on that thread, or [`classify_under`] elsewhere.
//! - **Outputs.** A token made with [`RunToken::recording`] records every output file the
//!   run commits ([`RunToken::outputs`]): the exact list of what it wrote, which no
//!   directory listing can give while other writers share the directory. [`StageSet`] makes
//!   such a token ([`StageSet::outputs`]).
//! - **In-memory endpoints (D84).** Inputs and outputs are names from
//!   [`gloci_io::memio::register_input`] and `register_output`. Every name registered must be
//!   discarded ([`gloci_io::memio::discard`]) once the run is over, whatever its result, or its
//!   source or sink stays alive: hold each in a [`gloci_io::memio::Registered`], which discards
//!   it on drop. Classify the run's error (and read [`error_site`]) before dropping them, since
//!   messages show the registered display names.
//! - **Messages.** Warnings and progress go through the `log` crate; the host installs its
//!   logger and sets the maximum level (`execute` installs neither). Render each record with
//!   [`crate::logging::render`] inside `Log::log`, which shows in-memory names by their
//!   display names.
//! - **Handlers.** `execute` never prints (the output goes to `-o`, which a host always gives),
//!   never installs signal or panic handlers and never exits the process. A panic in a command
//!   goes on in the caller, which should catch it (`catch_unwind`).
//!
//! A pipeline's steps run at once, chained by bounded in-memory pipes ([`pipe`]), in a
//! [`StageSet`] (D85): it does the token, the endpoint discarding, the classifying and the
//! panic catching above for each step, and reports one error for the set.

mod args;
mod bindings;
mod pipe;
mod signature;
mod spec;
mod stages;

pub use args::{
    KwValue, RESERVED_ALWAYS, RESERVED_ON_STEP_OF, RESERVED_ON_STEPS, Role,
    drop_needless_separator, input_arg, kwargs_to_argv, kwargs_to_argv_as, reserved,
    rewrite_flags_as_kwargs, step_argv,
};
pub use bindings::{
    BINDINGS, Bind, Binding, CLI_ONLY, MODULE_FUNCTIONS, Output, RECORDS, SOURCES, Shape, UNBOUND,
    binding, is_cli_only, is_unbound, method_name, module_function,
};
pub use gloci_io::cancel::{RunGuard, RunToken};
pub use pipe::{PIPE_CHUNK, PIPE_CHUNKS, PipeFail, PipeReader, PipeWriter, pipe};
pub use signature::{
    Kw, OUT_HELP, Signature, Skipped, function_inputs, function_signature, step_signature,
};
pub use spec::{ArgSpec, CommandSpec, FILE_VALUE_NAMES, command_specs, find_spec};
pub use stages::{HostStage, StageError, StageFn, StageSet, StageSpec};

use crate::cli::{self, Cli, Commands};
use crate::{Ctx, Reported};
use gloci_core::VcfSpan;
use gloci_io::cancel::{self, Watch};
use gloci_io::memio::scrub;
use gloci_io::{IoError, Loc, OrderError, UsageError};
use std::ffi::OsString;
use std::num::NonZeroUsize;
use std::path::PathBuf;

/// Parse an argv (`gloci` first) with the CLI's definition. The error is clap's message as
/// plain text: no colours, no `error: ` and no "For more information" line.
pub fn parse(args: &[OsString]) -> Result<Cli, UsageError> {
    crate::cli::try_parse_from(args).map_err(|e| UsageError(clap_message(&e)))
}

fn clap_message(e: &clap::Error) -> String {
    let text = e.render().to_string();
    let text = text.strip_prefix("error: ").unwrap_or(&text);
    let kept: Vec<&str> = text
        .lines()
        .filter(|l| !l.starts_with("For more information"))
        .collect();
    kept.join("\n").trim_end().to_string()
}

/// How a host reads a stage's output back as records (D89): what the output is, and the
/// footprint of its VCF records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReadBack {
    pub shape: Shape,
    /// The stage's `--vcf-span` when it has one, else the default, [`VcfSpan::Allele`]: the
    /// records read back have the footprints the command used.
    pub vcf_span: VcfSpan,
}

/// How the output of the stage `argv` (`gloci` first, as [`step_argv`] builds it) reads back:
/// its binding's [`Output`] with the flags given ([`Output::with`]), and its `--vcf-span`.
/// `groupby` writes records only when its keys (`-g`) start with columns 1, 2 and 3 (BED's
/// chrom, start and end), else a table. A command that is no step (a module source,
/// `windows -g`) writes records.
pub fn read_back(argv: &[OsString]) -> Result<ReadBack, UsageError> {
    let matches = crate::command()
        .try_get_matches_from(argv)
        .map_err(|e| UsageError(clap_message(&e)))?;
    let mut path = Vec::new();
    let mut m = &matches;
    while let Some((name, sub)) = m.subcommand() {
        path.push(name);
        m = sub;
    }
    let output = binding(&path).map_or(RECORDS, |b| b.output);
    let given = |id: &str| matches!(m.try_get_one::<bool>(id), Ok(Some(true)));
    let mut shape = output
        .with
        .iter()
        .find(|(id, _)| given(id))
        .map_or(output.shape, |(_, s)| *s);
    if path == ["groupby"] && shape == Shape::Records && !keyed_on_bed(m) {
        shape = Shape::Table;
    }
    let span = m
        .try_get_raw("vcf_span")
        .ok()
        .flatten()
        .and_then(|mut v| v.next());
    let vcf_span = if span.is_some_and(|v| v == "ref") {
        VcfSpan::Ref
    } else {
        VcfSpan::Allele
    };
    Ok(ReadBack { shape, vcf_span })
}

/// Whether `groupby`'s keys (`-g`, comma-separated, possibly repeated) start with columns 1, 2
/// and 3 in order, so that its rows start with a BED interval.
fn keyed_on_bed(m: &clap::ArgMatches) -> bool {
    let Ok(Some(values)) = m.try_get_raw("group") else {
        return false;
    };
    let keys: Vec<String> = values
        .flat_map(|v| {
            v.to_string_lossy()
                .split(',')
                .map(|k| k.trim().to_string())
                .collect::<Vec<_>>()
        })
        .collect();
    keys.len() >= 3 && keys[..3] == ["1", "2", "3"]
}

/// How a run is told apart from the command line.
#[derive(Debug, Clone, Default)]
pub struct ExecOptions {
    /// The command line recorded in output headers (`##glociCommand`, `@PG`), with no `-t`;
    /// in-memory names in it are shown by their display names.
    pub command_line: String,
    /// The spill directory when the command has no `--tmp-dir`; else `GLOCI_TMPDIR`, else the
    /// system temp directory.
    pub tmp_dir: Option<PathBuf>,
}

/// Run a parsed command: `Ok(false)` when it ran but found a problem (`check`), `Ok(true)`
/// otherwise.
///
/// The command runs on a worker of a rayon pool of its own, sized by its `-t` (else every
/// core), whose threads have the stack CRAM decoding needs ([`gloci_io::WORKER_STACK`]); the
/// process-global pool is never touched. The calling thread's run token ([`RunGuard`]) is
/// installed on every thread of the pool. See the module docs for what the host does.
pub fn execute(cli: &Cli, o: &ExecOptions) -> anyhow::Result<bool> {
    let (requested, tmp_arg) = thread_and_tmp(&cli.command);
    let threads = requested
        .unwrap_or_else(|| std::thread::available_parallelism().unwrap_or(NonZeroUsize::MIN));
    let ctx = Ctx {
        threads,
        tmp_dir: resolve_tmp_dir(tmp_arg.or_else(|| o.tmp_dir.clone())),
        command_line: scrub(&o.command_line),
        version: env!("CARGO_PKG_VERSION"),
    };
    let (result, all_ok) = in_run_pool(threads, &Watch::current(), || crate::dispatch(cli, &ctx))?;
    result.map(|()| all_ok)
}

/// Run `f` on a worker of a fresh rayon pool of `threads` threads with
/// [`gloci_io::WORKER_STACK`], each carrying `watch`'s token for its whole life: rayon work
/// that `f` starts runs on this pool, never on the global one. A panic in `f` goes on in the
/// caller.
fn in_run_pool<T: Send>(
    threads: NonZeroUsize,
    watch: &Watch,
    f: impl FnOnce() -> T + Send,
) -> anyhow::Result<T> {
    let w = watch.clone();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads.get())
        .stack_size(gloci_io::WORKER_STACK)
        .thread_name(|_| "gloci".into())
        .start_handler(move |_| w.install_for_thread())
        .build()
        .map_err(|e| {
            anyhow::anyhow!(
                "cannot start the command's threads ({e}); check the system's thread limits"
            )
        })?;
    Ok(pool.install(f))
}

/// The spill directory: `arg` (the command's `--tmp-dir`, which clap also fills from
/// `GLOCI_TMPDIR`, else the host's default), else a non-empty `GLOCI_TMPDIR`, else the system
/// temp dir.
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

/// The command's `-t` and `--tmp-dir`, for the commands that have them.
pub fn thread_and_tmp(command: &Commands) -> (Option<NonZeroUsize>, Option<PathBuf>) {
    match command {
        Commands::Genomes(_) => (None, None),
        Commands::Index(a) => (a.threads.threads, None),
        Commands::Sort(a) => (a.threads.threads, a.tmp_dir.tmp_dir.clone()),
        Commands::View(a) => (a.threads.threads, None),
        Commands::Cat(a) => (a.threads.threads, None),
        Commands::Check(a) => (a.threads.threads, None),
        Commands::Intersect(a) => (a.threads.threads, a.sort.tmp_dir.clone()),
        Commands::Merge(a) => (a.threads.threads, a.sort.tmp_dir.clone()),
        Commands::Cluster(a) => (a.threads.threads, a.sort.tmp_dir.clone()),
        Commands::Spacing(a) => (a.threads.threads, a.sort.tmp_dir.clone()),
        Commands::Complement(a) => (a.threads.threads, a.sort.tmp_dir.clone()),
        Commands::Partition(a) => (a.threads.threads, a.sort.tmp_dir.clone()),
        Commands::Genomecov(a) => (a.threads.threads, a.sort.tmp_dir.clone()),
        Commands::Slop(a) => (a.threads.threads, None),
        Commands::Flank(a) => (a.threads.threads, None),
        Commands::Shift(a) => (a.threads.threads, None),
        Commands::Windows(a) => (a.threads.threads, None),
        Commands::Convert(a) => (a.threads.threads, None),
        Commands::Getfasta(a) => (a.threads.threads, None),
        Commands::Nuc(a) => (a.threads.threads, None),
        Commands::Subtract(a) => (a.threads.threads, a.sort.tmp_dir.clone()),
        Commands::Closest(a) => (a.threads.threads, a.sort.tmp_dir.clone()),
        Commands::Map(a) => (a.threads.threads, a.sort.tmp_dir.clone()),
        Commands::Coverage(a) => (a.threads.threads, a.sort.tmp_dir.clone()),
        Commands::Groupby(a) => (a.threads.threads, None),
        Commands::Maskfasta(a) => (a.threads.threads, None),
        Commands::Stats(a) => match &a.command {
            cli::StatsCommand::Jaccard(j) => (j.threads.threads, j.sort.tmp_dir.clone()),
            cli::StatsCommand::Fisher(f) => (f.threads.threads, f.sort.tmp_dir.clone()),
            cli::StatsCommand::Reldist(r) => (r.threads.threads, r.sort.tmp_dir.clone()),
        },
        Commands::PairIntersect(a) => (a.threads.threads, None),
        Commands::PairPair(a) => (a.threads.threads, None),
        Commands::Random(a) => (a.threads.threads, None),
        Commands::Sample(a) => (a.threads.threads, None),
        Commands::Shuffle(a) => (a.threads.threads, None),
        Commands::Faidx(a) => (a.threads.threads, None),
        Commands::Dict(a) => (a.threads.threads, None),
        Commands::Flags(_) => (None, None),
        Commands::Idxstats(a) => (a.threads.threads, None),
        Commands::Flagstat(a) => (a.threads.threads, None),
        Commands::Head(a) => (a.threads.threads, None),
        Commands::Samples(_) => (None, None),
        Commands::Depth(a) => (a.threads.threads, None),
        Commands::Bedcov(a) => (a.threads.threads, None),
        Commands::Reheader(a) => (a.threads.threads, None),
        Commands::Igv(a) => (a.threads.threads, None),
        Commands::Links(a) => (a.threads.threads, None),
        Commands::Vcf(a) => match &a.command {
            cli::VcfCommand::Filter(f) => (f.threads.threads, None),
            cli::VcfCommand::Extract(e) => (e.threads.threads, e.tmp_dir.tmp_dir.clone()),
            cli::VcfCommand::Stats(e) => (e.threads.threads, None),
            cli::VcfCommand::Concat(c) => (c.threads.threads, None),
            cli::VcfCommand::Split(c) => (c.threads.threads, None),
            cli::VcfCommand::Norm(c) => (c.threads.threads, None),
            cli::VcfCommand::Annotate(c) => (c.threads.threads, None),
            cli::VcfCommand::Merge(c) => (c.threads.threads, None),
            cli::VcfCommand::Compare(c) => (c.threads.threads, None),
            cli::VcfCommand::Popgen(c) => (c.threads.threads, None),
        },
    }
}

/// What a failed run's error is, for the exit status (2, 0, 1) and the message, which never
/// holds `error: ` and shows in-memory names by their display names.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Outcome {
    /// A usage problem (exit 2).
    Usage(String),
    /// The output's reader went away: a quiet success (exit 0).
    BrokenPipe,
    /// The run was cancelled (D91), whatever error the cancellation surfaced as (exit 1,
    /// `interrupted`).
    Cancelled,
    /// The message was already logged ([`Reported`]; exit 1).
    Reported,
    /// Any other failure (exit 1).
    Failure(String),
}

/// [`classify_under`] the token of the calling thread (the run's, when called where the run's
/// [`RunGuard`] is entered).
pub fn classify(err: &anyhow::Error) -> Outcome {
    classify_under(err, cancel::current().as_ref())
}

/// What `err`, the error of a run whose token is `token`, is. A failed run whose token was
/// cancelled was cancelled, whichever error it surfaced as (a write that failed, an open, a
/// lost worker).
pub fn classify_under(err: &anyhow::Error, token: Option<&RunToken>) -> Outcome {
    if token.is_some_and(RunToken::is_cancelled) || err.chain().any(cancel::is_cancellation) {
        return Outcome::Cancelled;
    }
    if err.chain().any(|c| c.downcast_ref::<Reported>().is_some()) {
        return Outcome::Reported;
    }
    for cause in err.chain() {
        if let Some(u) = cause.downcast_ref::<UsageError>() {
            return Outcome::Usage(scrub(&u.0));
        }
        if let Some(IoError::Usage(u)) = cause.downcast_ref::<IoError>() {
            return Outcome::Usage(scrub(&u.0));
        }
        if let Some(e @ IoError::Gate { .. }) = cause.downcast_ref::<IoError>() {
            return Outcome::Usage(scrub(&e.to_string()));
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
    let msg = if err.downcast_ref::<IoError>().is_some() {
        err.to_string()
    } else {
        format!("{err:#}")
    };
    Outcome::Failure(scrub(&msg))
}

/// Which gloci-io error an [`ErrorSite`] comes from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiteKind {
    Open,
    Create,
    Read,
    Write,
    /// A malformed record or header.
    Parse,
    /// The format could not be told, or the input is empty.
    Detect,
    /// A BGZF file without its end-of-file marker.
    Truncated,
    /// A record out of sort order ([`ErrorSite::order`] has both records).
    Unsorted,
}

/// The structured fields of a failed run's data error, for a host's exception attributes
/// (D90): the first gloci-io error in the chain that names an input or output. The path, the
/// quoted records and the hint show in-memory names by their display names, like messages.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ErrorSite {
    pub kind: SiteKind,
    pub path: String,
    /// The line or record of a parse error, or of the later record of an unsorted pair.
    pub loc: Option<Loc>,
    /// Both records of an unsorted pair, and the hint.
    pub order: Option<OrderError>,
}

/// The [`ErrorSite`] of `err`, if it has one.
pub fn error_site(err: &anyhow::Error) -> Option<ErrorSite> {
    err.chain().find_map(|c| {
        if let Some(e) = c.downcast_ref::<IoError>() {
            return site_of_io(e);
        }
        c.downcast_ref::<OrderError>().map(site_of_order)
    })
}

fn site_of_io(e: &IoError) -> Option<ErrorSite> {
    let site = |kind, path: &str, loc| ErrorSite {
        kind,
        path: scrub(path),
        loc,
        order: None,
    };
    Some(match e {
        IoError::Open { path, .. } => site(SiteKind::Open, path, None),
        IoError::Create { path, .. } => site(SiteKind::Create, path, None),
        IoError::Read { path, .. } => site(SiteKind::Read, path, None),
        IoError::Write { path, .. } => site(SiteKind::Write, path, None),
        IoError::Parse { path, loc, .. } => site(SiteKind::Parse, path, Some(*loc)),
        IoError::Detect { path, .. } => site(SiteKind::Detect, path, None),
        IoError::Truncated { path } => site(SiteKind::Truncated, path, None),
        IoError::Order(o) => site_of_order(o),
        IoError::Gate { .. } | IoError::Usage(_) | IoError::Internal(_) | IoError::Cancelled => {
            return None;
        }
    })
}

fn site_of_order(o: &OrderError) -> ErrorSite {
    let order = OrderError {
        path: scrub(&o.path),
        prev_loc: o.prev_loc,
        prev: scrub(&o.prev),
        loc: o.loc,
        rec: scrub(&o.rec),
        hint: scrub(&o.hint),
    };
    ErrorSite {
        kind: SiteKind::Unsorted,
        path: order.path.clone(),
        loc: Some(o.loc),
        order: Some(order),
    }
}

#[cfg(test)]
mod tests;
