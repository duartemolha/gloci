//! The stage runtime (D85): the steps of a pipeline run at once, one thread each, joined by
//! bounded pipes ([`super::pipe`]), and the set reports one result.
//!
//! What the set reports ([`StageSet::wait`]):
//! - **Success**, `Ok(all_ok)` (the AND of the stages' results), when every stage succeeded or
//!   ended quietly. A broken pipe is quiet: the stage's reader went away, because the consumer
//!   stopped early or the next stage failed.
//! - **Cancelled**, when [`StageSet::cancel`] was called and any stage failed: the earliest
//!   failed stage, with [`Outcome::Cancelled`].
//! - Otherwise **the earliest primary error**. An error is derived when it is a broken pipe or
//!   comes from reading a pipe whose writer failed (`upstream step failed: ...`), and primary
//!   otherwise. So a data error in stage 0 beats the derived error it causes in stage 1, and a
//!   usage error in stage 1 beats the broken pipe it causes in stage 0. (With only derived
//!   errors, which a working pipeline never has, the earliest one.)
//!
//! A failure never cancels the other stages: it travels through the pipes. A stage that fails
//! fails the pipes it writes ([`StageSpec::fail_downstream`]), and drops its endpoints, so its
//! readers stop at once and its writers' next writes fail with a broken pipe.

use super::pipe::{PipeFail, mentions_upstream_failure, upstream_failure};
use super::{ErrorSite, ExecOptions, Outcome, RunGuard, RunToken, classify_under, error_site};
use crossbeam_channel::{Receiver, RecvTimeoutError, Sender, unbounded};
use gloci_io::IoError;
use gloci_io::memio::Registered;
use std::any::Any;
use std::ffi::OsString;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::PathBuf;
use std::thread::JoinHandle;
use std::time::{Duration, Instant};

/// One step of a pipeline: a command and the in-memory endpoints it reads and writes.
#[derive(Debug, Default)]
pub struct StageSpec {
    /// The whole argv, `gloci` first ([`super::step_argv`]).
    pub argv: Vec<OsString>,
    /// [`ExecOptions::command_line`].
    pub command_line: String,
    /// [`ExecOptions::tmp_dir`].
    pub tmp_dir: Option<PathBuf>,
    /// The registrations of the names in `argv`: pipe ends, and the host's sources and sinks.
    /// The stage's thread drops them when the stage ends, whatever its result (after judging
    /// its error, whose message shows their display names), so a stage that fails before it
    /// opens one still lets its neighbours go on.
    pub endpoints: Vec<Registered>,
    /// The fail handles of the pipes this stage writes ([`super::PipeWriter::fail_handle`]).
    /// When the stage fails, other than with a broken pipe, each is failed with its message
    /// before the endpoints are dropped, so the next stage gets `upstream step failed: {msg}`
    /// rather than a clean end of a cut-off stream. The handles also hold the pipes' EOF back
    /// until the stage has ended, so the failure always arrives first. The set binds its
    /// run token to these pipes, so a cancel stops a read or write parked on one from any
    /// thread.
    ///
    /// Holding EOF back assumes each stage writes at most one pipe (gloci commands have one
    /// `-o`, so a pipeline is a tree): a stage that wrote two pipes, read one after the other by
    /// one consumer, would never end the first while it waited to write the second.
    pub fail_downstream: Vec<PipeFail>,
}

/// The stage whose error the set reports, judged on its thread before its endpoints were
/// dropped.
#[derive(Debug)]
pub struct StageError {
    /// The stage's index in the set.
    pub stage: usize,
    pub error: anyhow::Error,
    /// What the error is ([`classify_under`] the set's token). It is [`Outcome::Cancelled`]
    /// whenever the set was cancelled, whatever `error` is.
    pub outcome: Outcome,
    /// [`error_site`]; `None` for a cancelled set.
    pub site: Option<ErrorSite>,
}

impl std::fmt::Display for StageError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&message(&self.outcome))
    }
}

impl std::error::Error for StageError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(self.error.as_ref())
    }
}

/// The message of a judged error, as a failed pipe's reader is told it.
fn message(outcome: &Outcome) -> String {
    match outcome {
        Outcome::Usage(m) | Outcome::Failure(m) => m.clone(),
        Outcome::Cancelled => "interrupted".into(),
        Outcome::BrokenPipe => "broken pipe".into(),
        // The command logged the message itself.
        Outcome::Reported => "the step's error was logged".into(),
    }
}

/// A stage thread's stack: the size of a main thread's, where the CLI parses its argv. It is
/// set rather than left to the default, which `RUST_MIN_STACK` can make too small for clap
/// (64 KiB overflows in a debug build); [`super::execute`] runs the command itself on its
/// pool's [`gloci_io::WORKER_STACK`] threads.
const STAGE_STACK: usize = 8 << 20;

/// A stage's work: a host stage's code ([`HostStage`]), or a test's ([`StageSet::spawn_fns`]).
pub type StageFn = Box<dyn FnOnce() -> anyhow::Result<bool> + Send>;

/// A stage that runs the host's own code rather than a command, under the set's token, panic
/// handling and error priority: the host's reader of the pipeline's output (the records it
/// hands to Python). It owns what it reads, such as the reading end of the last
/// command's pipe, and drops it when it ends; a command stage that writes that pipe keeps its
/// fail handle, so a failure upstream reaches it as `upstream step failed: ...`, a derived
/// error, and its own failures are primary.
pub struct HostStage {
    pub run: StageFn,
}

struct Stage {
    run: StageFn,
    endpoints: Vec<Registered>,
    fail_downstream: Vec<PipeFail>,
}

/// How a stage ended.
enum Ended {
    Ok(bool),
    Err {
        // Boxed: an error is rare and large (clippy::large_enum_variant).
        err: Box<StageError>,
        /// A broken pipe: the stage's reader went away.
        quiet: bool,
        /// Caused by another stage's failure: a broken pipe or a failed pipe read.
        derived: bool,
    },
}

/// The stages of one pipeline, running. Each runs on its own thread (`gloci-stage-{i}`) under
/// the set's one [`RunToken`].
///
/// Dropping the set cancels it if any stage thread is still running, then joins every stage
/// thread; a set whose stages have all ended is only joined.
///
/// Both [`wait`](Self::wait) and drop block until stage threads end, so a host must call them
/// with any lock its endpoints need released: a Python-backed source or sink takes the GIL to
/// read or write, and to be dropped. A cancel cannot interrupt a host read or write already in
/// progress, so a drop or wait holding the GIL while a stage is inside one never returns.
pub struct StageSet {
    token: RunToken,
    handles: Vec<Option<JoinHandle<()>>>,
    rx: Receiver<(usize, Ended)>,
    ended: Vec<Option<Ended>>,
    reported: bool,
}

impl std::fmt::Debug for StageSet {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("StageSet")
            .field("stages", &self.handles.len())
            .field("ended", &self.ended.iter().filter(|e| e.is_some()).count())
            .field("cancelled", &self.token.is_cancelled())
            .finish()
    }
}

impl StageSet {
    /// Start every stage: [`super::parse`] its argv and [`super::execute`] it, on its thread.
    /// The set's token records the output files the stages commit ([`StageSet::outputs`]).
    pub fn spawn(specs: Vec<StageSpec>) -> StageSet {
        StageSet::spawn_with_hosts(specs, Vec::new())
    }

    /// [`spawn`](Self::spawn) the command stages `specs`, and the host stages `hosts` after
    /// them (numbered from `specs.len()`), under one token.
    pub fn spawn_with_hosts(specs: Vec<StageSpec>, hosts: Vec<HostStage>) -> StageSet {
        let token = RunToken::recording();
        let stages: Vec<Stage> = specs
            .into_iter()
            .map(|s| {
                for f in &s.fail_downstream {
                    f.bind(&token);
                }
                let StageSpec {
                    argv,
                    command_line,
                    tmp_dir,
                    endpoints,
                    fail_downstream,
                } = s;
                let run: StageFn = Box::new(move || {
                    let cli = super::parse(&argv)?;
                    super::execute(
                        &cli,
                        &ExecOptions {
                            command_line,
                            tmp_dir,
                        },
                    )
                });
                Stage {
                    run,
                    endpoints,
                    fail_downstream,
                }
            })
            .chain(hosts.into_iter().map(|h| Stage {
                run: h.run,
                endpoints: Vec::new(),
                fail_downstream: Vec::new(),
            }))
            .collect();
        StageSet::start(token, stages)
    }

    /// Start stages that run closures instead of commands, under the same token, panic
    /// handling and error priority; each closure owns the endpoints it uses. For tests.
    #[doc(hidden)]
    pub fn spawn_fns(fns: Vec<StageFn>) -> StageSet {
        let stages = fns
            .into_iter()
            .map(|run| Stage {
                run,
                endpoints: Vec::new(),
                fail_downstream: Vec::new(),
            })
            .collect();
        StageSet::start(RunToken::new(), stages)
    }

    fn start(token: RunToken, stages: Vec<Stage>) -> StageSet {
        let (tx, rx) = unbounded();
        let mut handles = Vec::with_capacity(stages.len());
        let mut ended = Vec::with_capacity(stages.len());
        for (i, stage) in stages.into_iter().enumerate() {
            let fails = stage.fail_downstream.clone();
            let (tx, t) = (tx.clone(), token.clone());
            let spawned = std::thread::Builder::new()
                .name(format!("gloci-stage-{i}"))
                .stack_size(STAGE_STACK)
                .spawn(move || run_stage(i, stage, &t, &tx));
            match spawned {
                Ok(h) => {
                    handles.push(Some(h));
                    ended.push(None);
                }
                Err(e) => {
                    // The stage (and its endpoints) went with the closure; say why downstream.
                    let error = anyhow::anyhow!(
                        "cannot start a thread for step {i} ({e}); check the system's thread limits"
                    );
                    handles.push(None);
                    ended.push(Some(judge(i, Err(error), &token, &fails)));
                }
            }
        }
        StageSet {
            token,
            handles,
            rx,
            ended,
            reported: false,
        }
    }

    /// The set's result once every stage thread has ended (they are joined by then), or `None`
    /// if `timeout` passes first (`Duration::MAX` waits without a deadline). A host waits in
    /// short slices to look at its own signals in between. The result is returned once; later
    /// calls return `None` at once.
    pub fn wait(&mut self, timeout: Duration) -> Option<Result<bool, StageError>> {
        if self.reported {
            return None;
        }
        // A timeout too long to add to now (Duration::MAX) means no deadline.
        let deadline = Instant::now().checked_add(timeout);
        while self.ended.iter().any(Option::is_none) {
            let next = match deadline {
                Some(d) => self
                    .rx
                    .recv_timeout(d.saturating_duration_since(Instant::now())),
                None => self.rx.recv().map_err(|_| RecvTimeoutError::Disconnected),
            };
            match next {
                Ok((i, e)) => self.ended[i] = Some(e),
                Err(RecvTimeoutError::Timeout) => return None,
                Err(RecvTimeoutError::Disconnected) => {
                    // Every thread is gone; one that sent nothing died outside its guard.
                    for (i, slot) in self.ended.iter_mut().enumerate() {
                        if slot.is_none() {
                            let lost = IoError::Internal(format!(
                                "the thread of step {i} ended without a result"
                            ));
                            *slot = Some(judge(i, Err(lost.into()), &self.token, &[]));
                        }
                    }
                }
            }
        }
        self.join();
        self.reported = true;
        let ended = self.ended.drain(..).flatten().collect();
        Some(verdict(ended, &self.token).map_err(|e| *e))
    }

    /// Ask every stage to stop (D91): reads, writes and long loops fail at their next look at
    /// the token, and the set reports [`Outcome::Cancelled`].
    pub fn cancel(&self) {
        self.token.cancel();
    }

    /// The output files the stages have committed so far, each as its command named it, in the
    /// order they were committed: every file a run renamed into place
    /// ([`gloci_io::AtomicFile::commit`], [`gloci_io::StagedSet::publish`]), and no other.
    /// Once [`StageSet::wait`] has returned the result, this is everything the run wrote.
    pub fn outputs(&self) -> Vec<PathBuf> {
        self.token.outputs()
    }

    fn join(&mut self) {
        for h in self.handles.iter_mut().filter_map(Option::take) {
            let _ = h.join();
        }
    }
}

impl Drop for StageSet {
    fn drop(&mut self) {
        if self.handles.iter().flatten().any(|h| !h.is_finished()) {
            self.token.cancel();
        }
        self.join();
    }
}

/// A stage's thread: run it under the set's token, judge its result while its endpoints are
/// registered, fail its downstream pipes if it failed, drop its endpoints, and report.
fn run_stage(i: usize, stage: Stage, token: &RunToken, tx: &Sender<(usize, Ended)>) {
    let Stage {
        run,
        endpoints,
        fail_downstream,
    } = stage;
    let mut ended = {
        let _guard = RunGuard::enter(token);
        let result = catch_unwind(AssertUnwindSafe(run)).unwrap_or_else(|p| Err(panicked(&*p)));
        judge(i, result, token, &fail_downstream)
    };
    drop(fail_downstream);
    // A source or sink may run host code when dropped; a panic there is this stage's error.
    if let Err(p) = catch_unwind(AssertUnwindSafe(move || drop(endpoints)))
        && matches!(ended, Ended::Ok(_))
    {
        ended = judge(i, Err(panicked(&*p)), token, &[]);
    }
    let _ = tx.send((i, ended));
}

/// A panic as an error: [`IoError::Internal`], whose message says it is a bug.
fn panicked(payload: &(dyn Any + Send)) -> anyhow::Error {
    let what = payload
        .downcast_ref::<&str>()
        .map(|s| s.to_string())
        .or_else(|| payload.downcast_ref::<String>().cloned())
        .unwrap_or_else(|| "with no message".into());
    IoError::Internal(format!("the step panicked ({what})")).into()
}

/// Judge stage `i`'s result under the set's token and, when it failed other than with a broken
/// pipe, fail the pipes it writes: with the message the failed pipe it read carried (so a
/// chain of stages passes the first failure on), else with its own.
fn judge(i: usize, result: anyhow::Result<bool>, token: &RunToken, fails: &[PipeFail]) -> Ended {
    let error = match result {
        Ok(ok) => return Ended::Ok(ok),
        Err(e) => e,
    };
    let outcome = classify_under(&error, Some(token));
    let site = error_site(&error);
    let quiet = outcome == Outcome::BrokenPipe;
    let upstream = upstream_failure(&error);
    let derived = quiet || upstream.is_some() || mentions_upstream_failure(&format!("{error:#}"));
    if !quiet {
        let msg = upstream.unwrap_or_else(|| message(&outcome));
        for f in fails {
            f.fail(&msg);
        }
    }
    Ended::Err {
        err: Box::new(StageError {
            stage: i,
            error,
            outcome,
            site,
        }),
        quiet,
        derived,
    }
}

/// The set's result from its stages', in stage order (see the module docs).
fn verdict(ended: Vec<Ended>, token: &RunToken) -> Result<bool, Box<StageError>> {
    let mut all_ok = true;
    let mut failed = Vec::new();
    for e in ended {
        match e {
            Ended::Ok(ok) => all_ok &= ok,
            Ended::Err { quiet: true, .. } => {}
            Ended::Err { err, derived, .. } => failed.push((err, derived)),
        }
    }
    if failed.is_empty() {
        return Ok(all_ok);
    }
    if token.is_cancelled() {
        let (mut err, _) = failed.swap_remove(0);
        err.outcome = Outcome::Cancelled;
        err.site = None;
        return Err(err);
    }
    let pick = failed.iter().position(|(_, derived)| !derived).unwrap_or(0);
    Err(failed.swap_remove(pick).0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::embed::pipe::pipe;
    use std::io::{self, Read};

    fn failure(msg: &str) -> anyhow::Error {
        anyhow::anyhow!("{msg}")
    }

    fn broken_pipe() -> anyhow::Error {
        IoError::Write {
            path: "<out>".into(),
            source: io::Error::from(io::ErrorKind::BrokenPipe),
        }
        .into()
    }

    /// The error a command returns after reading a pipe failed with `msg`.
    fn read_of_failed_pipe(msg: &str) -> anyhow::Error {
        let (w, mut r) = pipe(1, 4);
        w.fail(msg);
        let source = r.read(&mut [0u8; 4]).unwrap_err();
        anyhow::Error::new(IoError::Read {
            path: "<pipe>".into(),
            source,
        })
        .context("while merging")
    }

    fn judged(results: Vec<anyhow::Result<bool>>, token: &RunToken) -> Vec<Ended> {
        results
            .into_iter()
            .enumerate()
            .map(|(i, r)| judge(i, r, token, &[]))
            .collect()
    }

    fn picked(
        results: Vec<anyhow::Result<bool>>,
        token: &RunToken,
    ) -> Result<bool, Box<StageError>> {
        verdict(judged(results, token), token)
    }

    #[test]
    fn success_is_every_stage_ok_or_quiet() {
        let t = RunToken::new();
        assert!(picked(vec![Ok(true), Ok(true)], &t).unwrap());
        // `check` finding a problem in one stage.
        assert!(!picked(vec![Ok(true), Ok(false)], &t).unwrap());
        assert!(picked(vec![Err(broken_pipe()), Err(broken_pipe()), Ok(true)], &t).unwrap());
    }

    #[test]
    fn the_earliest_primary_error_wins() {
        let t = RunToken::new();
        let e = picked(
            vec![
                Err(broken_pipe()),
                Err(read_of_failed_pipe("x")),
                Err(failure("stage 2")),
                Err(failure("stage 3")),
            ],
            &t,
        )
        .unwrap_err();
        assert_eq!(e.stage, 2);
        assert_eq!(e.outcome, Outcome::Failure("stage 2".into()));
        // A failed read that only kept the message is derived too.
        let text_only = failure("error reading <p>: upstream step failed: boom");
        let e = picked(vec![Err(text_only), Err(failure("real"))], &t).unwrap_err();
        assert_eq!(e.stage, 1);
        // Only derived errors: the earliest is reported rather than a success.
        let e = picked(vec![Ok(true), Err(read_of_failed_pipe("x"))], &t).unwrap_err();
        assert_eq!(e.stage, 1);
    }

    #[test]
    fn a_cancelled_set_reports_cancelled_from_the_earliest_failed_stage() {
        let t = RunToken::new();
        // Judged before the cancel: a data error with a site.
        let site_err: anyhow::Error = IoError::Parse {
            path: "a.bed".into(),
            loc: gloci_io::Loc::Line(2),
            msg: "bad".into(),
        }
        .into();
        let mut ended = judged(vec![Ok(true), Err(site_err)], &t);
        t.cancel();
        ended.extend(judged(vec![Err(failure("later"))], &t));
        let e = verdict(ended, &t).unwrap_err();
        assert_eq!(
            (e.stage, &e.outcome, &e.site),
            (1, &Outcome::Cancelled, &None)
        );
        // A cancel that came after every stage succeeded changes nothing.
        assert!(picked(vec![Ok(true), Ok(true)], &t).unwrap());
    }

    #[test]
    fn a_failed_stage_fails_its_downstream_pipes_with_the_first_failure() {
        let t = RunToken::new();
        // Its own error, when it is primary.
        let (w, mut r) = pipe(1, 4);
        judge(0, Err(failure("bad line 2")), &t, &[w.fail_handle()]);
        let e = r.read(&mut [0u8; 4]).unwrap_err();
        assert_eq!(e.to_string(), "upstream step failed: bad line 2");
        // The message it read, when it failed reading a failed pipe.
        let (w, mut r) = pipe(1, 4);
        judge(
            1,
            Err(read_of_failed_pipe("bad line 2")),
            &t,
            &[w.fail_handle()],
        );
        let e = r.read(&mut [0u8; 4]).unwrap_err();
        assert_eq!(e.to_string(), "upstream step failed: bad line 2");
        // Not when it ended with a broken pipe: its reader is gone anyway.
        let (w, mut r) = pipe(1, 4);
        let f = w.fail_handle();
        judge(0, Err(broken_pipe()), &t, std::slice::from_ref(&f));
        drop((w, f));
        assert_eq!(r.read(&mut [0u8; 4]).unwrap(), 0);
    }
}
