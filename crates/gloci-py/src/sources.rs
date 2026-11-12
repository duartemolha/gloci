//! Python-object sources (D88): file objects, iterables of lines or tuples, and Arrow streams
//! ([`crate::frames`]), each read by one stage as an in-memory input (D84).
//!
//! A source has two halves:
//! - a [`SourceReader`], the input the stage reads: chunks from a bounded channel, so reading
//!   never needs the GIL;
//! - a [`Pump`], kept by the driver, which turns the Python object's items into bytes and sends
//!   them. It runs in the driver loop with the GIL held, at most [`TURN_BYTES`] a turn and only
//!   while the channel has room, so the driver never waits holding the GIL and a source is never
//!   in memory whole (R4).
//!
//! **Errors are delivered in the data's order.** When a source raises (its iterator, or its
//! `read()`), or holds what no line can (a tab in a value, a null coordinate), the bytes before
//! it are sent, then the reader's next read fails. The stage fails there as it would at a bad
//! line of a file, so an error earlier in the data is still the one reported. When the source's
//! failure is the run's error, its exception is raised: for an exception the source raised, the
//! same object (D90). A `BaseException` that is not an `Exception` (`KeyboardInterrupt`,
//! `SystemExit`) cancels the run at once instead.

use crate::errors;
use crate::frames::Frame;
use crate::records::Record;
use crossbeam_channel::{Receiver, RecvTimeoutError, Sender, TrySendError, bounded};
use gloci::embed::StageError;
use gloci_io::cancel::{self, Watch};
use pyo3::exceptions::{PyException, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{
    PyByteArray, PyBytes, PyDict, PyFloat, PyInt, PyIterator, PyList, PyString, PyTuple,
};
use std::io::{self, Read};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, OnceLock};
use std::time::Duration;

/// The bytes a pump sends at once.
const CHUNK: usize = 256 << 10;
/// The chunks a source's channel holds: with the chunk being read and the one the pump is
/// making, about 1.5 MiB per source.
const CHUNKS: usize = 4;
/// The most a pump reads from its Python object in one turn of the driver loop (R4).
pub(crate) const TURN_BYTES: usize = 1 << 20;
/// How often a read parked on an empty channel looks at the run's cancel token.
const POLL: Duration = Duration::from_millis(20);

/// What kind of Python object a source is (`classify_source` in `gloci._pipeline`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    /// An object with `read()`: bytes, or text encoded as UTF-8.
    File,
    /// An iterable of lines (`str` or `bytes`) or tuples of values.
    Lines,
    /// An object with `__arrow_c_stream__`.
    Frame,
}

/// A plan's source: `{"kind": "file"|"lines"|"frame", "obj": ..., "display": "<list>"}`.
pub(crate) struct SourceSpec {
    kind: Kind,
    obj: Py<PyAny>,
    display: String,
}

pub(crate) fn read_spec(token: &str, spec: &Bound<'_, PyAny>) -> PyResult<SourceSpec> {
    let invalid = |msg: String| PyValueError::new_err(format!("invalid gloci plan: {msg}"));
    let spec = spec
        .cast::<PyDict>()
        .map_err(|_| invalid(format!("source {token:?} is not a dict")))?;
    let item = |key: &str| {
        spec.get_item(key)?
            .ok_or_else(|| invalid(format!("source {token:?} has no \"{key}\"")))
    };
    let kind = match item("kind")?.extract::<String>()?.as_str() {
        "file" => Kind::File,
        "lines" => Kind::Lines,
        "frame" => Kind::Frame,
        other => return Err(invalid(format!("source {token:?} has the kind {other:?}"))),
    };
    let display: String = item("display")?.extract()?;
    if display.contains('\0') {
        return Err(invalid(format!(
            "source {token:?} has a NUL in its display name"
        )));
    }
    Ok(SourceSpec {
        kind,
        obj: item("obj")?.unbind(),
        display,
    })
}

/// What the two halves share.
#[derive(Debug, Default)]
struct Shared {
    /// The source ended cleanly: the reader's end of the channel is its EOF.
    finished: AtomicBool,
    /// The source failed after the bytes sent: the reader's end of the channel is this error.
    failed: OnceLock<String>,
    /// A read returned the failure.
    delivered: AtomicBool,
}

/// The error a stage's read gets when its source failed.
#[derive(Debug)]
struct SourceFailed(String);

impl std::fmt::Display for SourceFailed {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for SourceFailed {}

/// The input a stage reads.
pub(crate) struct SourceReader {
    rx: Receiver<Vec<u8>>,
    /// Read chunks go back to the pump, so a long source reuses a few buffers.
    spare: Sender<Vec<u8>>,
    cur: Vec<u8>,
    pos: usize,
    shared: Arc<Shared>,
    /// The run's cancel token, as the reading thread has it (taken at the first read).
    watch: Option<Watch>,
    eof: bool,
}

impl Read for SourceReader {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        loop {
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
                    if done.capacity() > 0 {
                        done.clear();
                        let _ = self.spare.try_send(done);
                    }
                }
                Err(RecvTimeoutError::Timeout) => {
                    let watch = self.watch.get_or_insert_with(Watch::current);
                    if watch.is_cancelled() || cancel::is_cancelled() {
                        return Err(cancel::io_error());
                    }
                }
                // Every chunk sent has been read; how the source ended decides the rest.
                Err(RecvTimeoutError::Disconnected) => {
                    if let Some(msg) = self.shared.failed.get() {
                        self.shared.delivered.store(true, Ordering::SeqCst);
                        return Err(io::Error::other(SourceFailed(msg.clone())));
                    }
                    if self.shared.finished.load(Ordering::SeqCst) {
                        self.eof = true;
                        continue;
                    }
                    // The driver let the source go without ending it: the run was cancelled.
                    return Err(cancel::io_error());
                }
            }
        }
    }
}

/// What a pump turn left.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Turn {
    /// The channel has room: the turn's budget was spent.
    Room,
    /// The channel is full.
    Full,
    /// The source has nothing to give now (a non-blocking `read()` returned None).
    Starved,
    /// Nothing more to send: the source ended or failed, or its reader is gone.
    Done,
}

enum End {
    Finished,
    Failed,
}

/// What filling a chunk from the Python object gave.
pub(crate) enum Fill {
    /// The chunk is full, and there may be more.
    More,
    Starved,
    /// The source is exhausted.
    End,
    /// The source failed: the exception to raise if this failure is the run's error.
    Failed(PyErr),
    /// A `BaseException` that is not an `Exception`: the run is cancelled at once.
    Interrupt(PyErr),
}

/// Where the bytes come from.
enum Feed {
    /// `read1()` (one read of the stream, so a pipe's data comes as it arrives), else `read()`.
    File {
        read: Py<PyAny>,
    },
    Lines {
        iter: Py<PyIterator>,
    },
    Frame(Box<Frame>),
}

/// The driver's half of a source.
pub(crate) struct Pump {
    /// The stage that reads it.
    stage: usize,
    display: String,
    feed: Feed,
    tx: Option<Sender<Vec<u8>>>,
    spare: Receiver<Vec<u8>>,
    shared: Arc<Shared>,
    /// A chunk made but not sent yet: the channel was full.
    pending: Option<Vec<u8>>,
    /// How the source ended, said once `pending` is sent.
    end: Option<End>,
    /// The items (lines, tuples) or rows read so far.
    items: u64,
    /// The failure's message, as the reader gets it.
    message: String,
    /// What to raise when the source's failure is the run's error.
    error: Option<PyErr>,
}

/// Open the source `spec`, read by stage `stage`; `k` numbers the plan's sources.
pub(crate) fn open(
    py: Python<'_>,
    spec: SourceSpec,
    k: usize,
    stage: usize,
) -> PyResult<(Pump, SourceReader, String)> {
    let obj = spec.obj.bind(py);
    let feed = match spec.kind {
        Kind::File => {
            let read = if obj.hasattr("read1")? {
                obj.getattr("read1")?
            } else {
                obj.getattr("read")?
            };
            Feed::File {
                read: read.unbind(),
            }
        }
        Kind::Lines => Feed::Lines {
            iter: obj.try_iter()?.unbind(),
        },
        Kind::Frame => Feed::Frame(Box::new(Frame::open(py, obj, &spec.display)?)),
    };
    let (tx, rx) = bounded(CHUNKS);
    let (spare_tx, spare_rx) = bounded(2);
    let shared = Arc::new(Shared::default());
    let reader = SourceReader {
        rx,
        spare: spare_tx,
        cur: Vec::new(),
        pos: 0,
        shared: shared.clone(),
        watch: None,
        eof: false,
    };
    let pump = Pump {
        stage,
        message: format!("the Python source {} failed (source {k})", spec.display),
        display: spec.display.clone(),
        feed,
        tx: Some(tx),
        spare: spare_rx,
        shared,
        pending: None,
        end: None,
        items: 0,
        error: None,
    };
    Ok((pump, reader, spec.display))
}

impl Pump {
    /// One turn: send what was made, then read on from the Python object, until the channel
    /// is full, [`TURN_BYTES`] were sent, or the source ends. An `Err` is a `BaseException`
    /// that is not an `Exception`, which must cancel the run.
    pub(crate) fn pump(&mut self, py: Python<'_>) -> PyResult<Turn> {
        let mut sent = 0;
        loop {
            let Some(tx) = &self.tx else {
                return Ok(Turn::Done);
            };
            if let Some(chunk) = self.pending.take() {
                let n = chunk.len();
                match tx.try_send(chunk) {
                    Ok(()) => sent += n,
                    Err(TrySendError::Full(chunk)) => {
                        self.pending = Some(chunk);
                        return Ok(Turn::Full);
                    }
                    // The stage is done with it (it ended, or needs no more): read no further.
                    Err(TrySendError::Disconnected(_)) => {
                        self.tx = None;
                        self.end = None;
                        return Ok(Turn::Done);
                    }
                }
            }
            if let Some(end) = self.end.take() {
                match end {
                    End::Finished => self.shared.finished.store(true, Ordering::SeqCst),
                    End::Failed => {
                        let _ = self.shared.failed.set(self.message.clone());
                    }
                }
                // The reader sees the end once it has read every chunk.
                self.tx = None;
                return Ok(Turn::Done);
            }
            if sent >= TURN_BYTES {
                return Ok(Turn::Room);
            }
            let mut buf = self
                .spare
                .try_recv()
                .unwrap_or_else(|_| Vec::with_capacity(CHUNK + (CHUNK >> 4)));
            let fill = match &mut self.feed {
                Feed::File { read } => fill_file(py, read.bind(py), &mut buf, &self.display),
                Feed::Lines { iter } => {
                    fill_lines(py, iter.bind(py), &mut buf, &mut self.items, &self.display)
                }
                Feed::Frame(f) => f.fill(py, &mut buf, CHUNK, &mut self.items, &self.display),
            };
            if !buf.is_empty() {
                self.pending = Some(buf);
            }
            match fill {
                Fill::More => {}
                Fill::Starved if self.pending.is_none() => return Ok(Turn::Starved),
                Fill::Starved => {}
                Fill::End => self.end = Some(End::Finished),
                Fill::Failed(e) => {
                    self.error = Some(e);
                    self.end = Some(End::Failed);
                }
                Fill::Interrupt(e) => return Err(e),
            }
        }
    }

    /// The sending end, for the driver to wait on for room (only while it pumps).
    pub(crate) fn sender(&self) -> Option<Sender<Vec<u8>>> {
        self.tx.clone()
    }

    /// Let the source go without ending it, as the run is being cancelled: its reader fails
    /// with the cancellation.
    pub(crate) fn abandon(&mut self) {
        self.tx = None;
        self.pending = None;
        self.end = None;
    }

    /// The exception to raise for the run's error `e`, when `e` is this source's failure as
    /// its stage read it.
    pub(crate) fn claim(&mut self, e: &StageError) -> Option<PyErr> {
        if e.stage != self.stage || !self.shared.delivered.load(Ordering::SeqCst) {
            return None;
        }
        let ours = e.error.chain().any(|c| {
            c.downcast_ref::<SourceFailed>()
                .or_else(|| {
                    c.downcast_ref::<io::Error>()
                        .and_then(io::Error::get_ref)
                        .and_then(|i| i.downcast_ref::<SourceFailed>())
                })
                .is_some_and(|f| f.0 == self.message)
        }) || format!("{:#}", e.error).contains(&self.message);
        if ours { self.error.take() } else { None }
    }
}

/// An exception a source raised: kept to raise later, unless it must stop the run now.
fn raised(py: Python<'_>, e: PyErr) -> Fill {
    if e.is_instance_of::<PyException>(py) {
        Fill::Failed(e)
    } else {
        Fill::Interrupt(e)
    }
}

/// Read up to [`CHUNK`] from a file object: `read1(n)` or `read(n)` returns bytes, or str for
/// a text file (encoded as UTF-8), empty at its end, or None when a non-blocking stream has
/// nothing now.
fn fill_file(py: Python<'_>, read: &Bound<'_, PyAny>, buf: &mut Vec<u8>, display: &str) -> Fill {
    while buf.len() < CHUNK {
        let got = match read.call1((CHUNK - buf.len(),)) {
            Ok(got) => got,
            Err(e) => return raised(py, e),
        };
        let before = buf.len();
        if got.is_none() {
            return Fill::Starved;
        } else if let Ok(b) = got.cast::<PyBytes>() {
            buf.extend_from_slice(b.as_bytes());
        } else if let Ok(s) = got.cast::<PyString>() {
            match s.to_str() {
                Ok(s) => buf.extend_from_slice(s.as_bytes()),
                Err(e) => return raised(py, e),
            }
        } else if let Ok(b) = got.cast::<PyByteArray>() {
            buf.extend_from_slice(&b.to_vec());
        } else {
            let what = type_name(&got);
            return Fill::Failed(PyTypeError::new_err(format!(
                "{display}: read() returned {what}; a file object's read() returns bytes or str"
            )));
        }
        if buf.len() == before {
            return Fill::End;
        }
    }
    Fill::More
}

fn type_name(obj: &Bound<'_, PyAny>) -> String {
    obj.get_type()
        .name()
        .map_or_else(|_| "an object".into(), |n| n.to_string())
}

/// Why an item cannot be a line.
enum Bad {
    /// A data problem, said after `{display} item {n}: `.
    Data(String),
    /// An exception from the object (its `__str__`, or a str that is not valid UTF-8).
    Raised(PyErr),
}

impl From<PyErr> for Bad {
    fn from(e: PyErr) -> Bad {
        Bad::Raised(e)
    }
}

/// Read items from an iterable until the chunk is full: each a line.
fn fill_lines(
    py: Python<'_>,
    iter: &Bound<'_, PyIterator>,
    buf: &mut Vec<u8>,
    items: &mut u64,
    display: &str,
) -> Fill {
    let mut iter = iter.clone();
    while buf.len() < CHUNK {
        let item = match iter.next() {
            None => return Fill::End,
            Some(Err(e)) => return raised(py, e),
            Some(Ok(item)) => item,
        };
        *items += 1;
        let mark = buf.len();
        if let Err(bad) = write_item(&item, buf) {
            buf.truncate(mark);
            return match bad {
                Bad::Raised(e) => raised(py, e),
                Bad::Data(msg) => {
                    let n = *items;
                    Fill::Failed(errors::data_at(
                        py,
                        format!("{display} item {}: {msg}", gloci_core::thousands(n)),
                        display,
                        n,
                    ))
                }
            };
        }
    }
    Fill::More
}

/// One item as a line: a `gloci.Record`'s text, a `str` or `bytes` line (its line break added
/// if it has none), or a tuple (or list) of values, tab-joined.
fn write_item(item: &Bound<'_, PyAny>, buf: &mut Vec<u8>) -> Result<(), Bad> {
    if let Ok(r) = item.cast::<Record>() {
        // A record of an output (D89): its text, as the output had it (a FASTA record's
        // holds its sequence lines).
        buf.extend_from_slice(r.get().line_str(item.py()).to_str()?.as_bytes());
        buf.push(b'\n');
        return Ok(());
    }
    if let Ok(s) = item.cast::<PyString>() {
        return write_line(s.to_str()?.as_bytes(), buf);
    }
    if let Ok(b) = item.cast::<PyBytes>() {
        return write_line(b.as_bytes(), buf);
    }
    if let Ok(t) = item.cast::<PyTuple>() {
        return write_row(t.iter(), t.len(), "tuple", buf);
    }
    if let Ok(l) = item.cast::<PyList>() {
        return write_row(l.iter(), l.len(), "list", buf);
    }
    Err(Bad::Data(format!(
        "{} is not a line or a row; an iterable source yields lines (str or bytes) or tuples \
of values",
        a_type(item)
    )))
}

/// "an int", "a dict".
fn a_type(obj: &Bound<'_, PyAny>) -> String {
    let name = type_name(obj);
    let article = if name.starts_with(['a', 'e', 'i', 'o', 'u']) {
        "an"
    } else {
        "a"
    };
    format!("{article} {name}")
}

fn write_line(line: &[u8], buf: &mut Vec<u8>) -> Result<(), Bad> {
    let body = line.strip_suffix(b"\n").unwrap_or(line);
    if memchr::memchr(b'\n', body).is_some() {
        return Err(Bad::Data(
            "the line holds a line break inside it; give one line per item".into(),
        ));
    }
    buf.extend_from_slice(body);
    buf.push(b'\n');
    Ok(())
}

/// A row's values, tab-joined (D88): `None` is `.`, anything else its `str()` (computed here
/// for `str`, `int` and `float`, as Python writes them).
fn write_row<'py>(
    values: impl Iterator<Item = Bound<'py, PyAny>>,
    len: usize,
    what: &str,
    buf: &mut Vec<u8>,
) -> Result<(), Bad> {
    if len == 0 {
        return Err(Bad::Data(format!("the {what} is empty")));
    }
    for (i, v) in values.enumerate() {
        if i > 0 {
            buf.push(b'\t');
        }
        if v.is_none() {
            buf.push(b'.');
            continue;
        }
        let at = buf.len();
        if v.is_exact_instance_of::<PyInt>() {
            if let Ok(n) = v.extract::<i64>() {
                buf.extend_from_slice(itoa::Buffer::new().format(n).as_bytes());
                continue;
            }
        } else if v.is_exact_instance_of::<PyFloat>() {
            push_float(v.extract::<f64>()?, buf);
            continue;
        }
        if let Ok(s) = v.cast_exact::<PyString>() {
            buf.extend_from_slice(s.to_str()?.as_bytes());
        } else {
            buf.extend_from_slice(v.str()?.to_str()?.as_bytes());
        }
        if let Some(what) = breaks(&buf[at..]) {
            let shown = v.repr().and_then(|r| r.to_str().map(str::to_string))?;
            return Err(Bad::Data(format!(
                "value {} ({shown}) holds {what}, which would split the line; replace or \
remove it",
                i + 1
            )));
        }
    }
    buf.push(b'\n');
    Ok(())
}

/// What in `value` would split a line: a tab or a line break.
pub(crate) fn breaks(value: &[u8]) -> Option<&'static str> {
    match memchr::memchr3(b'\t', b'\n', b'\r', value).map(|i| value[i]) {
        Some(b'\t') => Some("a tab"),
        Some(_) => Some("a line break"),
        None => None,
    }
}

/// `x` as Python's `str(x)` writes it: the shortest digits that read back as `x`, in fixed
/// notation when the exponent is from -4 to 15 (with `.0` for a whole number), else as
/// `1.5e+16`; `inf`, `-inf` and `nan`.
pub(crate) fn push_float(x: f64, buf: &mut Vec<u8>) {
    if x.is_nan() {
        buf.extend_from_slice(b"nan");
        return;
    }
    if x.is_infinite() {
        buf.extend_from_slice(if x > 0.0 { b"inf" } else { b"-inf" });
        return;
    }
    // Rust's `{:e}` gives the shortest round-trip digits too: "-1.2345e-7", "0e0".
    let sci = format!("{x:e}");
    let (mantissa, exp) = sci.split_once('e').expect("{:e} has an exponent");
    let exp: i32 = exp.parse().expect("{:e} has an integer exponent");
    let (negative, mantissa) = match mantissa.strip_prefix('-') {
        Some(m) => (true, m),
        None => (false, mantissa),
    };
    let digits: Vec<u8> = mantissa.bytes().filter(|&c| c != b'.').collect();
    if negative {
        buf.push(b'-');
    }
    if (-4..16).contains(&exp) {
        if exp >= 0 {
            let whole = exp as usize + 1;
            if digits.len() <= whole {
                buf.extend_from_slice(&digits);
                buf.resize(buf.len() + whole - digits.len(), b'0');
                buf.extend_from_slice(b".0");
            } else {
                buf.extend_from_slice(&digits[..whole]);
                buf.push(b'.');
                buf.extend_from_slice(&digits[whole..]);
            }
        } else {
            buf.extend_from_slice(b"0.");
            buf.resize(buf.len() + (-exp - 1) as usize, b'0');
            buf.extend_from_slice(&digits);
        }
    } else {
        buf.push(digits[0]);
        if digits.len() > 1 {
            buf.push(b'.');
            buf.extend_from_slice(&digits[1..]);
        }
        buf.push(b'e');
        buf.push(if exp < 0 { b'-' } else { b'+' });
        let e = exp.unsigned_abs();
        if e < 10 {
            buf.push(b'0');
        }
        buf.extend_from_slice(itoa::Buffer::new().format(e).as_bytes());
    }
}
