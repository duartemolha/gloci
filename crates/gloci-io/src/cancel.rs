//! Cooperative cancellation (D91). A host that embeds gloci (the Python library) makes one
//! [`RunToken`] per run, installs it on the threads that run work ([`RunGuard`]), and calls
//! [`RunToken::cancel`] when the user interrupts. Everything a run builds (readers, outputs,
//! sorters, worker pools) captures the thread's token when it is constructed, as a [`Watch`];
//! a cancelled watch makes reads and writes fail with [`IoError::Cancelled`]-classified errors
//! and long loops stop at their next poll, after which the usual error path removes any
//! partial output *file*. Runs never share a token, so cancelling one leaves the others alone,
//! and a token that is never cancelled costs one relaxed load per poll. The CLI installs no
//! token: its run gets one that is never cancelled, so its watch never fires.
//!
//! **The run's outputs.** A token made with [`RunToken::recording`] also records each output
//! file the run commits: an [`crate::AtomicFile`] renamed into place, or every file of a
//! [`crate::StagedSet`] once it is published. The file captures the token when it is created,
//! like a reader, so the record reaches every thread of the run the way the cancel flag does,
//! and a host learns exactly what its run wrote ([`RunToken::outputs`]), whatever else
//! changes on disk meanwhile. Without a recording token (the CLI) nothing is recorded.
//!
//! **Standard input.** A run reads standard input for at most one of its inputs: a second
//! would find it empty (or with the bytes the first left) and the run would succeed with a
//! wrong answer. The token holds which input takes it. A host gives the inputs the command line
//! names `-` before the run starts ([`RunToken::give_stdin`]), so that two of them are refused
//! by the names the user wrote; and every reader of `-` claims it as it opens
//! ([`claim_stdin`], through [`crate::open_stdin`]), so that a second reader is refused even
//! when no host named it. Runs never share a token, so one run's use of standard input never
//! reaches the next. Without a token (code that runs outside any run) nothing is checked.

use crate::error::IoError;
use std::cell::RefCell;
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, PoisonError};

/// One run's token: its cancel flag, and the output files it committed when it records them. Clones
/// share both.
#[derive(Clone, Debug, Default)]
pub struct RunToken(Arc<RunState>);

#[derive(Debug, Default)]
struct RunState {
    cancelled: AtomicBool,
    /// The committed output files, in the order they were committed; `None` when not recorded.
    outputs: Option<Mutex<Vec<PathBuf>>>,
    stdin: Mutex<StdinUse>,
}

/// Which input of a run takes standard input.
#[derive(Debug, Default)]
struct StdinUse {
    /// The input the command line names `-`, as the user wrote it ([`RunToken::give_stdin`]).
    given: Option<String>,
    /// Whether a reader has opened standard input ([`claim_stdin`]).
    opened: bool,
}

/// The usage error for standard input given to two inputs, `first` and `second`.
fn stdin_twice(first: &str, second: &str) -> IoError {
    IoError::Usage(crate::UsageError(format!(
        "standard input is given for both {first} and {second}; give one of them as a file"
    )))
}

impl RunToken {
    pub fn new() -> RunToken {
        RunToken::default()
    }

    /// A token that also records the output files the run commits ([`RunToken::outputs`]).
    pub fn recording() -> RunToken {
        RunToken(Arc::new(RunState {
            outputs: Some(Mutex::new(Vec::new())),
            ..RunState::default()
        }))
    }

    /// A [`Watch`] of this token, for code that runs the run's work on other threads.
    pub fn watch(&self) -> Watch {
        Watch(Some(self.clone()))
    }

    /// Record that the command line gives `-` (standard input) to `input`, named as the user
    /// wrote it (`-a`, `-b (value 2)`, `A`). The second input given `-` in the run is a usage
    /// error that names both. Call it for each such input before the run reads any: the reader
    /// that then opens `-` takes what was given ([`claim_stdin`]).
    pub fn give_stdin(&self, input: &str) -> Result<(), IoError> {
        let mut s = self.0.stdin.lock().unwrap_or_else(PoisonError::into_inner);
        if let Some(first) = &s.given {
            return Err(stdin_twice(first, input));
        }
        s.given = Some(input.to_string());
        Ok(())
    }

    fn claim_stdin(&self) -> Result<(), IoError> {
        let mut s = self.0.stdin.lock().unwrap_or_else(PoisonError::into_inner);
        if s.opened {
            let msg = match &s.given {
                Some(first) => format!(
                    "standard input is given for {first} and read for another input too; give \
one of them as a file"
                ),
                None => "standard input is read for two inputs; give one of them as a file".into(),
            };
            return Err(IoError::Usage(crate::UsageError(msg)));
        }
        s.opened = true;
        Ok(())
    }

    /// Ask the run to stop.
    pub fn cancel(&self) {
        self.0.cancelled.store(true, Ordering::SeqCst);
    }

    #[inline]
    pub fn is_cancelled(&self) -> bool {
        self.0.cancelled.load(Ordering::Relaxed)
    }

    /// The output files the run has committed so far, in the order it committed them, each as
    /// the command named it; empty for a token that does not record.
    pub fn outputs(&self) -> Vec<PathBuf> {
        self.0.outputs.as_ref().map_or_else(Vec::new, |o| {
            o.lock().unwrap_or_else(PoisonError::into_inner).clone()
        })
    }

    fn record(&self, path: &Path) {
        if let Some(o) = &self.0.outputs {
            o.lock()
                .unwrap_or_else(PoisonError::into_inner)
                .push(path.to_path_buf());
        }
    }
}

thread_local! {
    static CURRENT: RefCell<Option<RunToken>> = const { RefCell::new(None) };
}

/// The token installed on this thread, if any.
pub fn current() -> Option<RunToken> {
    CURRENT.with(|c| c.borrow().clone())
}

/// Take standard input for the reader that is opening `-`, under this thread's token: the first
/// reader of a run takes it, and a second is a usage error, whichever inputs they are. Every
/// reader of `-` calls it before it reads ([`crate::open_stdin`]). Without a token nothing is
/// checked.
pub fn claim_stdin() -> Result<(), IoError> {
    match current() {
        Some(t) => t.claim_stdin(),
        None => Ok(()),
    }
}

/// Whether this thread's token is cancelled.
pub fn is_cancelled() -> bool {
    CURRENT.with(|c| c.borrow().as_ref().is_some_and(RunToken::is_cancelled))
}

/// Installs a token on the current thread; the previous one is restored on drop. A thread a
/// run spawns enters the run's token at its start.
#[must_use = "the token is uninstalled when the guard drops"]
pub struct RunGuard {
    prev: Option<RunToken>,
    // Restores thread-local state: must drop on the thread that made it.
    _not_send: std::marker::PhantomData<*const ()>,
}

impl RunGuard {
    pub fn enter(token: &RunToken) -> RunGuard {
        RunGuard::install(Some(token.clone()))
    }

    fn install(t: Option<RunToken>) -> RunGuard {
        RunGuard {
            prev: CURRENT.with(|c| c.replace(t)),
            _not_send: std::marker::PhantomData,
        }
    }
}

impl Drop for RunGuard {
    fn drop(&mut self) {
        let prev = self.prev.take();
        CURRENT.with(|c| *c.borrow_mut() = prev);
    }
}

/// A token captured when something is built, to poll later from any thread. Empty (never
/// fires) when no run installed one.
#[derive(Clone, Debug, Default)]
pub struct Watch(Option<RunToken>);

impl Watch {
    /// Capture this thread's token.
    pub fn current() -> Watch {
        Watch(current())
    }

    #[inline]
    pub fn is_cancelled(&self) -> bool {
        self.0.as_ref().is_some_and(RunToken::is_cancelled)
    }

    /// `Err(IoError::Cancelled)` once cancelled: the poll for long loops.
    #[inline]
    pub fn check(&self) -> Result<(), IoError> {
        if self.is_cancelled() {
            Err(IoError::Cancelled)
        } else {
            Ok(())
        }
    }

    /// The `io::Error` a read or write fails with once cancelled.
    #[inline]
    pub fn check_io(&self) -> io::Result<()> {
        if self.is_cancelled() {
            Err(io_error())
        } else {
            Ok(())
        }
    }

    /// Record `path` as an output file the run committed, when the captured token records
    /// ([`RunToken::recording`]); nothing otherwise, as in the CLI.
    pub fn record_output(&self, path: &Path) {
        if let Some(t) = &self.0 {
            t.record(path);
        }
    }

    /// Install the captured token on this thread (for a thread a run spawns).
    pub fn enter(&self) -> RunGuard {
        RunGuard::install(self.0.clone())
    }

    /// Install the captured token on this thread for the rest of its life, replacing any: for
    /// the start handler of a thread pool a run builds, whose workers outlive any guard.
    pub fn install_for_thread(&self) {
        let t = self.0.clone();
        let prev = CURRENT.with(|c| c.replace(t));
        drop(prev);
    }
}

/// `Err(IoError::Cancelled)` when this thread's token is cancelled. Loops should capture a
/// [`Watch`] once instead.
pub fn check() -> Result<(), IoError> {
    if is_cancelled() {
        Err(IoError::Cancelled)
    } else {
        Ok(())
    }
}

/// The `io::Error` a cancelled read or write fails with.
pub fn io_error() -> io::Error {
    io::Error::other(Cancelled)
}

/// Whether `e` is the error of a cancelled read or write.
pub fn is_cancel_error(e: &io::Error) -> bool {
    e.get_ref().is_some_and(|i| i.is::<Cancelled>())
}

/// Whether `e`, or anything in its `source()` chain, is a cancellation: an
/// [`IoError`] for which [`IoError::is_cancelled`] holds, or an `io::Error` carrying the cancel
/// marker. (`io::Error::source` skips a custom inner error, so that is looked at directly.)
pub fn is_cancellation(e: &(dyn std::error::Error + 'static)) -> bool {
    let mut cur: Option<&(dyn std::error::Error + 'static)> = Some(e);
    while let Some(err) = cur {
        if let Some(io) = err.downcast_ref::<IoError>()
            && io.is_cancelled()
        {
            return true;
        }
        if let Some(io) = err.downcast_ref::<io::Error>()
            && is_cancel_error(io)
        {
            return true;
        }
        cur = err.source();
    }
    false
}

#[derive(Debug)]
struct Cancelled;

impl std::fmt::Display for Cancelled {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("interrupted")
    }
}

impl std::error::Error for Cancelled {}

/// A reader that fails once the run that built it is cancelled.
pub struct CancelRead<R> {
    inner: R,
    watch: Watch,
}

impl<R> CancelRead<R> {
    pub fn new(inner: R) -> Self {
        CancelRead {
            inner,
            watch: Watch::current(),
        }
    }
}

impl<R: Read> Read for CancelRead<R> {
    #[inline]
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.watch.check_io()?;
        self.inner.read(buf)
    }
}

impl<R: io::Seek> io::Seek for CancelRead<R> {
    fn seek(&mut self, pos: io::SeekFrom) -> io::Result<u64> {
        self.watch.check_io()?;
        self.inner.seek(pos)
    }
}
