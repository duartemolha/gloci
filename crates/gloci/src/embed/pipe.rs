//! Bounded in-memory pipes between the steps of a pipeline (D85): the writer is the output of
//! one step and the reader the input of the next, each registered as an in-memory endpoint
//! ([`gloci_io::memio::Registered`]). The pipe holds at most `chunks` chunks of `chunk_bytes`,
//! plus the one the writer is filling, so a fast step waits for a slow one instead of growing
//! memory.
//!
//! - **EOF.** Dropping the writer sends what it holds, then ends the stream; when fail handles
//!   ([`PipeWriter::fail_handle`]) exist, the stream ends once they are dropped too.
//! - **Failure.** [`PipeWriter::fail`] or [`PipeFail::fail`] makes every later read fail with
//!   `upstream step failed: {msg}` instead of reading on or ending cleanly.
//! - **Broken pipe.** Once the reader is dropped, writes fail with
//!   [`std::io::ErrorKind::BrokenPipe`].
//! - **Cancellation (D91).** A read or write parked on an empty or full pipe fails with the
//!   cancellation error once the run's token is cancelled. The token is the calling thread's
//!   ([`gloci_io::cancel::RunGuard`]), or the one captured when the pipe was made, or the one
//!   [`crate::embed::StageSet::spawn`] binds through the pipe's fail handle.

use crossbeam_channel::{Receiver, RecvTimeoutError, SendTimeoutError, Sender, bounded};
use gloci_io::cancel::{self, RunToken};
use std::io::{self, Read, Write};
use std::sync::{Arc, Mutex, OnceLock};
use std::time::Duration;

/// How many chunks a pipeline's pipe holds.
pub const PIPE_CHUNKS: usize = 4;
/// The size of a pipeline pipe's chunk: the writer sends when it has this much.
pub const PIPE_CHUNK: usize = 1 << 20;

/// How often a parked read or write checks the run token for a cancel.
const POLL: Duration = Duration::from_millis(20);

const UPSTREAM: &str = "upstream step failed";

/// What the two ends share.
#[derive(Debug)]
struct State {
    /// The first failure reported; later ones are ignored.
    failed: OnceLock<String>,
    /// The run's token: captured when the pipe is made, replaced by the stage set's.
    token: Mutex<Option<RunToken>>,
}

impl State {
    fn fail(&self, msg: &str) {
        let _ = self.failed.set(msg.to_string());
    }

    /// Whether the calling thread's run, or the pipe's, was cancelled. Only asked while
    /// parked, so the lock is taken once per [`POLL`] at most.
    fn cancelled(&self) -> bool {
        cancel::is_cancelled()
            || self
                .token
                .lock()
                .unwrap_or_else(|p| p.into_inner())
                .as_ref()
                .is_some_and(RunToken::is_cancelled)
    }
}

/// The error a read of a failed pipe returns.
#[derive(Debug)]
struct UpstreamFailed(String);

impl std::fmt::Display for UpstreamFailed {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{UPSTREAM}: {}", self.0)
    }
}

impl std::error::Error for UpstreamFailed {}

/// The message a failed pipe carried, when `err` comes from reading one. (`io::Error::source`
/// skips a custom payload, so each `io::Error` is looked into directly.)
pub(crate) fn upstream_failure(err: &anyhow::Error) -> Option<String> {
    err.chain().find_map(|c| {
        c.downcast_ref::<UpstreamFailed>()
            .or_else(|| {
                c.downcast_ref::<io::Error>()
                    .and_then(io::Error::get_ref)
                    .and_then(|i| i.downcast_ref::<UpstreamFailed>())
            })
            .map(|u| u.0.clone())
    })
}

/// Whether `text`, an error's message, says it read a failed pipe: for errors that kept the
/// message but not the error itself.
pub(crate) fn mentions_upstream_failure(text: &str) -> bool {
    text.contains(&format!("{UPSTREAM}: "))
}

/// A bounded pipe of `chunks` chunks of `chunk_bytes` (each at least 1); a pipeline uses
/// [`PIPE_CHUNKS`] and [`PIPE_CHUNK`].
pub fn pipe(chunks: usize, chunk_bytes: usize) -> (PipeWriter, PipeReader) {
    let (tx, rx) = bounded(chunks.max(1));
    // Drained chunks go back to the writer, so a long stream reuses a few buffers.
    let (spare_tx, spare_rx) = bounded(2);
    let state = Arc::new(State {
        failed: OnceLock::new(),
        token: Mutex::new(cancel::current()),
    });
    let w = PipeWriter {
        tx,
        spare: spare_rx,
        buf: Vec::new(),
        chunk: chunk_bytes.max(1),
        state: state.clone(),
        broken: false,
    };
    let r = PipeReader {
        rx,
        spare: spare_tx,
        cur: Vec::new(),
        pos: 0,
        state,
        eof: false,
    };
    (w, r)
}

fn broken() -> io::Error {
    io::Error::from(io::ErrorKind::BrokenPipe)
}

/// The writing end. It buffers a chunk before sending it; [`Write::flush`] sends a partial one.
#[derive(Debug)]
pub struct PipeWriter {
    tx: Sender<Vec<u8>>,
    spare: Receiver<Vec<u8>>,
    buf: Vec<u8>,
    chunk: usize,
    state: Arc<State>,
    /// The reader is gone.
    broken: bool,
}

impl PipeWriter {
    /// A handle that fails this pipe, and keeps its stream from ending until it is dropped. It
    /// works after the writer itself has been moved into a command's output and dropped, so a
    /// stage that fails can say so to the next one. Clones share the pipe.
    pub fn fail_handle(&self) -> PipeFail {
        PipeFail {
            state: self.state.clone(),
            _hold: self.tx.clone(),
        }
    }

    /// End the stream with a failure: the reader gets `upstream step failed: {msg}` instead of
    /// what is left or EOF.
    pub fn fail(mut self, msg: &str) {
        self.state.fail(msg);
        self.buf.clear();
    }

    /// Send the buffer, waiting while the pipe is full.
    fn send_buf(&mut self) -> io::Result<()> {
        if self.buf.is_empty() {
            return Ok(());
        }
        if self.broken {
            return Err(broken());
        }
        let mut chunk = std::mem::take(&mut self.buf);
        loop {
            match self.tx.send_timeout(chunk, POLL) {
                Ok(()) => return Ok(()),
                Err(SendTimeoutError::Timeout(back)) => {
                    if self.state.cancelled() {
                        return Err(cancel::io_error());
                    }
                    chunk = back;
                }
                Err(SendTimeoutError::Disconnected(_)) => {
                    self.broken = true;
                    return Err(broken());
                }
            }
        }
    }
}

impl Write for PipeWriter {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        if data.is_empty() {
            return Ok(0);
        }
        if self.broken {
            return Err(broken());
        }
        if self.buf.capacity() == 0 {
            self.buf = self
                .spare
                .try_recv()
                .unwrap_or_else(|_| Vec::with_capacity(self.chunk));
        }
        let n = data.len().min(self.chunk - self.buf.len());
        self.buf.extend_from_slice(&data[..n]);
        if self.buf.len() == self.chunk {
            self.send_buf()?;
        }
        Ok(n)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.send_buf()
    }
}

impl Drop for PipeWriter {
    /// Send what is left (not after a failure or a cancel), then let the stream end.
    fn drop(&mut self) {
        if !self.buf.is_empty()
            && !self.broken
            && self.state.failed.get().is_none()
            && !self.state.cancelled()
        {
            let _ = self.send_buf();
        }
    }
}

/// Fails a pipe from outside its writer ([`PipeWriter::fail_handle`]). While a handle lives,
/// the reader does not see EOF: drop every handle once the writing step has ended (the stage
/// runtime does). That is deadlock-free when each step writes at most one pipe, as in a
/// pipeline of gloci commands (one `-o` each): the reader waits only for its own writer's
/// step, which waits for nothing but its inputs and that reader.
#[derive(Debug, Clone)]
pub struct PipeFail {
    state: Arc<State>,
    /// A sender held only so that the stream does not end before the handle says how.
    _hold: Sender<Vec<u8>>,
}

impl PipeFail {
    /// The reader's later reads fail with `upstream step failed: {msg}`. The first failure
    /// reported wins. Never blocks.
    pub fn fail(&self, msg: &str) {
        self.state.fail(msg);
    }

    /// Let the token of the stage set that owns this pipe stop a parked read or write, from
    /// any thread.
    pub(crate) fn bind(&self, token: &RunToken) {
        *self.state.token.lock().unwrap_or_else(|p| p.into_inner()) = Some(token.clone());
    }
}

/// The reading end.
#[derive(Debug)]
pub struct PipeReader {
    rx: Receiver<Vec<u8>>,
    spare: Sender<Vec<u8>>,
    cur: Vec<u8>,
    pos: usize,
    state: Arc<State>,
    eof: bool,
}

impl Read for PipeReader {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        loop {
            // A failure is seen at once, before any data still in the pipe.
            if let Some(msg) = self.state.failed.get() {
                return Err(io::Error::other(UpstreamFailed(msg.clone())));
            }
            if self.pos < self.cur.len() {
                let n = out.len().min(self.cur.len() - self.pos);
                out[..n].copy_from_slice(&self.cur[self.pos..self.pos + n]);
                self.pos += n;
                return Ok(n);
            }
            if self.eof {
                return Ok(0);
            }
            match self.rx.recv_timeout(POLL) {
                Ok(next) => {
                    let mut done = std::mem::replace(&mut self.cur, next);
                    self.pos = 0;
                    done.clear();
                    if done.capacity() > 0 {
                        // Full, or the writer is gone: the buffer is simply freed.
                        let _ = self.spare.try_send(done);
                    }
                }
                Err(RecvTimeoutError::Timeout) => {
                    if self.state.cancelled() {
                        return Err(cancel::io_error());
                    }
                }
                // Look at `failed` once more: a failure comes before the end.
                Err(RecvTimeoutError::Disconnected) => self.eof = true,
            }
        }
    }
}
