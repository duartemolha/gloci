//! Cooperative cancellation (D91). A host that embeds gloci (the Python library) makes one
//! [`CancelToken`] per run, installs it on the threads that run work ([`RunGuard`]), and calls
//! [`CancelToken::cancel`] when the user interrupts. Everything a run builds (readers, outputs,
//! sorters, worker pools) captures the thread's token when it is constructed, as a [`Watch`];
//! a cancelled watch makes reads and writes fail with [`IoError::Cancelled`]-classified errors
//! and long loops stop at their next poll, after which the usual error path removes any
//! partial output *file*. Runs never share a token, so cancelling one leaves the others alone,
//! and a token that is never cancelled costs one relaxed load per poll. The CLI installs no
//! token: its watch is empty and never fires.

use crate::error::IoError;
use std::cell::RefCell;
use std::io::{self, Read};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

/// The cancel flag of one run. Clones share it.
#[derive(Clone, Debug, Default)]
pub struct CancelToken(Arc<AtomicBool>);

impl CancelToken {
    pub fn new() -> CancelToken {
        CancelToken::default()
    }

    /// Ask the run to stop.
    pub fn cancel(&self) {
        self.0.store(true, Ordering::SeqCst);
    }

    #[inline]
    pub fn is_cancelled(&self) -> bool {
        self.0.load(Ordering::Relaxed)
    }
}

thread_local! {
    static CURRENT: RefCell<Option<CancelToken>> = const { RefCell::new(None) };
}

/// The token installed on this thread, if any.
pub fn current() -> Option<CancelToken> {
    CURRENT.with(|c| c.borrow().clone())
}

/// Whether this thread's token is cancelled.
pub fn is_cancelled() -> bool {
    CURRENT.with(|c| c.borrow().as_ref().is_some_and(CancelToken::is_cancelled))
}

/// Installs a token on the current thread; the previous one is restored on drop. A thread a
/// run spawns enters the run's token at its start.
#[must_use = "the token is uninstalled when the guard drops"]
pub struct RunGuard {
    prev: Option<CancelToken>,
    // Restores thread-local state: must drop on the thread that made it.
    _not_send: std::marker::PhantomData<*const ()>,
}

impl RunGuard {
    pub fn enter(token: &CancelToken) -> RunGuard {
        RunGuard::install(Some(token.clone()))
    }

    fn install(t: Option<CancelToken>) -> RunGuard {
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
pub struct Watch(Option<CancelToken>);

impl Watch {
    /// Capture this thread's token.
    pub fn current() -> Watch {
        Watch(current())
    }

    #[inline]
    pub fn is_cancelled(&self) -> bool {
        self.0.as_ref().is_some_and(CancelToken::is_cancelled)
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

    /// Install the captured token on this thread (for a thread a run spawns).
    pub fn enter(&self) -> RunGuard {
        RunGuard::install(self.0.clone())
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
