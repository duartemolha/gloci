//! Cooperative cancellation (D91). A host that embeds gloci (the Python library) calls
//! [`request_cancel`] when the user interrupts a run; every read through [`open_input`]
//! (crate::open_input) and every write to an [`Output`](crate::Output) then fails with
//! [`IoError::Cancelled`], long loops poll [`check`], and the usual error path removes the
//! partial output. The CLI never cancels: it has its own signal handlers.

use crate::error::IoError;
use std::io::{self, Read};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

static CANCEL: AtomicBool = AtomicBool::new(false);
static GUARDS: AtomicUsize = AtomicUsize::new(0);

/// Marks a run. While at least one guard is alive a cancellation stays requested; when the
/// last one drops the flag is cleared, so the next run starts clean.
#[must_use = "the cancel flag is cleared when the last RunGuard drops"]
pub struct RunGuard(());

impl RunGuard {
    pub fn enter() -> RunGuard {
        GUARDS.fetch_add(1, Ordering::SeqCst);
        RunGuard(())
    }
}

impl Drop for RunGuard {
    fn drop(&mut self) {
        if GUARDS.fetch_sub(1, Ordering::SeqCst) == 1 {
            CANCEL.store(false, Ordering::SeqCst);
        }
    }
}

/// Ask every running gloci read, write and long loop to stop.
pub fn request_cancel() {
    CANCEL.store(true, Ordering::SeqCst);
}

#[inline]
pub fn is_cancelled() -> bool {
    CANCEL.load(Ordering::Relaxed)
}

/// `Err(IoError::Cancelled)` once a cancellation was requested: the poll for long loops.
#[inline]
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

#[derive(Debug)]
struct Cancelled;

impl std::fmt::Display for Cancelled {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("interrupted")
    }
}

impl std::error::Error for Cancelled {}

/// A reader that fails once a cancellation was requested.
pub struct CancelRead<R>(pub R);

impl<R: Read> Read for CancelRead<R> {
    #[inline]
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if is_cancelled() {
            return Err(io_error());
        }
        self.0.read(buf)
    }
}

/// A writer that fails once a cancellation was requested.
pub struct CancelWrite<W>(pub W);

impl<W: io::Write> io::Write for CancelWrite<W> {
    #[inline]
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if is_cancelled() {
            return Err(io_error());
        }
        self.0.write(buf)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.0.flush()
    }
}
