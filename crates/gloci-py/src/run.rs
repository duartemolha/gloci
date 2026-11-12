//! Running a pipeline's plan (D85, D91): its stages on threads joined by bounded pipes
//! ([`StageSet`]), driven from the calling Python thread.
//!
//! A plan is `{"stages": [{"argv", "command_line"}], "sources": {token: source},
//! "terminal": "save"|"bytes"|"count", "names"}`. Each stage's argv is its whole command line,
//! `gloci` first, with these tokens where paths go, which become in-memory names (D84):
//! - `"\0pipe:<j>"`: the pipe stage `j` writes (its `-o`), read by exactly one later stage; a
//!   chain's steps, and a nested pipeline (a Pipeline given as a source), are such stages;
//! - `"\0out"`: the last stage's output, collected for the terminal (`bytes`, `count`);
//! - `"\0src:<k>"`: the Python-object source of `"sources"` with that token
//!   (`{"kind", "obj", "display"}`, [`sources`]), read by exactly one stage.
//!
//! **The driver loop** holds the GIL only to pump the Python-object sources (at most
//! [`sources::TURN_BYTES`] each a turn, and only into room their channel has), forward log
//! records and check signals. It releases the GIL while it waits, at most [`TURN`] at a time,
//! for room in a source's channel or for the stages to end. A signal handler's exception
//! (`KeyboardInterrupt`) cancels the run, which removes its partial output files, and is raised
//! once every stage has ended. The stage set is dropped with the GIL released (its drop joins
//! the stage threads, R2).

use crate::sources::{self, Pump, SourceSpec, Turn};
use crate::{errors, logs};
use crossbeam_channel::{Select, Sender};
use gloci::embed::{
    self, CommandSpec, PIPE_CHUNK, PIPE_CHUNKS, PipeFail, StageError, StageSet, StageSpec,
};
use gloci_io::memio::Registered;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict};
use std::collections::BTreeMap;
use std::ffi::OsString;
use std::io::{self, Write};
use std::sync::{Arc, Mutex, MutexGuard, PoisonError};
use std::time::Duration;

/// The longest the driver waits for the stages with the GIL released before it looks at
/// signals and logs again.
pub(crate) const TURN: Duration = Duration::from_millis(50);

/// The longest the driver waits before it asks a starved source (a non-blocking file object
/// that had nothing) again.
const STARVED: Duration = Duration::from_millis(5);

const OUT: &str = "\0out";
const PIPE: &str = "\0pipe:";
const SRC: &str = "\0src:";

/// What the caller gets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Terminal {
    /// `save(path)`: the last stage writes the file; the result is `None`.
    Save,
    /// `to_bytes()`: the last stage's output.
    Bytes,
    /// `count()`: the number of records in it.
    Count,
}

struct PlanStage {
    argv: Vec<String>,
    command_line: String,
}

fn invalid(msg: impl Into<String>) -> PyErr {
    PyValueError::new_err(format!("invalid gloci plan: {}", msg.into()))
}

/// A plan's parts: its stages, its terminal, and its Python-object sources by token.
type Plan = (Vec<PlanStage>, Terminal, BTreeMap<String, SourceSpec>);

fn read_plan(plan: &Bound<'_, PyDict>) -> PyResult<Plan> {
    let item = |key: &str| {
        plan.get_item(key)?
            .ok_or_else(|| invalid(format!("no \"{key}\"")))
    };
    let mut stages = Vec::new();
    for s in item("stages")?.try_iter()? {
        let s = s?;
        stages.push(PlanStage {
            argv: s.get_item("argv")?.extract()?,
            command_line: s.get_item("command_line")?.extract()?,
        });
    }
    if stages.is_empty() {
        return Err(invalid("no stages"));
    }
    let terminal = match item("terminal")?.extract::<String>()?.as_str() {
        "save" => Terminal::Save,
        "bytes" => Terminal::Bytes,
        "count" => Terminal::Count,
        other => return Err(invalid(format!("unknown terminal {other:?}"))),
    };
    let mut sources = BTreeMap::new();
    if let Some(given) = plan.get_item("sources")?
        && !given.is_none()
    {
        let given = given
            .cast::<PyDict>()
            .map_err(|_| invalid("\"sources\" is not a dict"))?;
        for (token, spec) in given.iter() {
            let token: String = token.extract()?;
            if source_index(&token).is_none() {
                return Err(invalid(format!("bad source token {token:?}")));
            }
            let spec = sources::read_spec(&token, &spec)?;
            sources.insert(token, spec);
        }
    }
    Ok((stages, terminal, sources))
}

/// `k` of a source token `"\0src:<k>"`.
fn source_index(token: &str) -> Option<usize> {
    let k = token.strip_prefix(SRC)?;
    (!k.is_empty() && k.bytes().all(|c| c.is_ascii_digit()))
        .then(|| k.parse().ok())
        .flatten()
}

/// The command a stage runs, from its argv (`gloci vcf filter ...`).
pub(crate) fn spec_of(argv: &[String]) -> Option<&'static CommandSpec> {
    let words = argv.get(1..).unwrap_or_default();
    (1..=words.len().min(2))
        .rev()
        .find_map(|n| embed::find_spec(&words[..n]))
}

/// A sink the terminal reads once the run is over.
#[derive(Default)]
struct Shared<T>(Arc<Mutex<T>>);

impl<T> Clone for Shared<T> {
    fn clone(&self) -> Self {
        Shared(self.0.clone())
    }
}

impl<T> Shared<T> {
    fn lock(&self) -> MutexGuard<'_, T> {
        self.0.lock().unwrap_or_else(PoisonError::into_inner)
    }
}

impl Write for Shared<Vec<u8>> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.lock().extend_from_slice(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// Counts the records of a text output as it is written.
/// - FASTA (the first non-blank line starts with `>`): its `>` lines, one per sequence.
/// - Anything else: its lines, but blank ones, `#` header and comment lines, `track` and
///   `browser` lines, and SAM's `@XY` header lines before the first record.
///
/// Known gap, for P2b Task 8: a summary table whose header row has no `#` (`stats jaccard`,
/// `stats fisher`, `stats reldist`, ...) counts that row as a record. Task 8 makes `count()`
/// count what iteration reads, through its re-read of the output
/// (`test_count_summary_tables`, a strict xfail until then).
#[derive(Default)]
struct RecordCount {
    records: u64,
    /// The first bytes of the line being written, enough to tell a header line.
    head: [u8; 8],
    head_len: usize,
    /// Bytes of the current line have been written.
    in_line: bool,
    /// Whether the output is FASTA, once its first non-blank line has ended.
    fasta: Option<bool>,
}

impl RecordCount {
    fn take(&mut self, bytes: &[u8]) {
        if bytes.is_empty() {
            return;
        }
        self.in_line = true;
        let n = (self.head.len() - self.head_len).min(bytes.len());
        self.head[self.head_len..self.head_len + n].copy_from_slice(&bytes[..n]);
        self.head_len += n;
    }

    fn end_line(&mut self) {
        let head = &self.head[..self.head_len];
        let blank = !self.in_line || head == b"\r";
        if !blank {
            let fasta = *self.fasta.get_or_insert(head[0] == b'>');
            let word = |w: &[u8]| {
                head.starts_with(w)
                    && matches!(head.get(w.len()), None | Some(b' ' | b'\t' | b'\r'))
            };
            let sam_header = self.records == 0
                && head.len() >= 3
                && head[0] == b'@'
                && head[1..3].iter().all(u8::is_ascii_alphabetic)
                && matches!(head.get(3), None | Some(b'\t' | b'\r'));
            let record = if fasta {
                head[0] == b'>'
            } else {
                !(head[0] == b'#' || word(b"track") || word(b"browser") || sam_header)
            };
            if record {
                self.records += 1;
            }
        }
        self.head_len = 0;
        self.in_line = false;
    }

    /// The count, with a last line that has no line break.
    fn total(&mut self) -> u64 {
        if self.in_line {
            self.end_line();
        }
        self.records
    }
}

impl Write for Shared<RecordCount> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let mut c = self.lock();
        let mut rest = buf;
        while let Some(i) = memchr::memchr(b'\n', rest) {
            c.take(&rest[..i]);
            c.end_line();
            rest = &rest[i + 1..];
        }
        c.take(rest);
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// Where the last stage's output goes.
enum Sink {
    None,
    Bytes(Shared<Vec<u8>>),
    Count(Shared<RecordCount>),
}

/// The stage specs of a plan, with every token replaced by an in-memory name registered for
/// the stage that uses it; the terminal's sink; and the pumps of the Python-object sources.
fn build(
    py: Python<'_>,
    stages: &[PlanStage],
    terminal: Terminal,
    mut sources: BTreeMap<String, SourceSpec>,
) -> PyResult<(Vec<StageSpec>, Sink, Vec<Pump>)> {
    let n = stages.len();
    let mut endpoints: Vec<Vec<Registered>> = (0..n).map(|_| Vec::new()).collect();
    let mut fails: Vec<Vec<PipeFail>> = (0..n).map(|_| Vec::new()).collect();
    // Pipe j: the name of its writing end, and its reading end until a stage takes it.
    let mut writers: BTreeMap<usize, String> = BTreeMap::new();
    let mut readers: BTreeMap<usize, Option<Registered>> = BTreeMap::new();
    let mut sink = Sink::None;
    let mut argvs: Vec<Vec<OsString>> = Vec::with_capacity(n);
    let mut pumps = Vec::new();

    // Make every pipe its writer names, first, so that readers can take theirs in any order.
    for (j, s) in stages.iter().enumerate() {
        let writes = format!("{PIPE}{j}");
        let uses = s.argv.iter().filter(|a| **a == writes).count();
        if uses == 0 {
            continue;
        }
        if uses > 1 || j + 1 == n {
            return Err(invalid(format!("stage {j} writes its pipe {uses} times")));
        }
        let display = format!(
            "<{} output>",
            spec_of(&s.argv).map_or_else(|| "step".into(), embed::method_name)
        );
        let (w, r) = embed::pipe(PIPE_CHUNKS, PIPE_CHUNK);
        // Ruling R1: the one fail handle goes to the writing stage, and no copy is kept.
        fails[j].push(w.fail_handle());
        let w = Registered::output(&display, Box::new(w));
        writers.insert(j, w.name().to_string());
        endpoints[j].push(w);
        readers.insert(j, Some(Registered::input(&display, Box::new(r))));
    }

    for (i, s) in stages.iter().enumerate() {
        let mut argv = Vec::with_capacity(s.argv.len());
        for a in &s.argv {
            let name = if let Some(j) = a.strip_prefix(PIPE) {
                let j: usize = j
                    .parse()
                    .map_err(|_| invalid(format!("bad pipe token {a:?}")))?;
                if j == i {
                    // This stage writes the pipe, registered above.
                    writers[&j].clone()
                } else {
                    let r = readers
                        .get_mut(&j)
                        .ok_or_else(|| invalid(format!("stage {i} reads pipe {j}, never written")))?
                        .take()
                        .ok_or_else(|| invalid(format!("pipe {j} is read twice")))?;
                    if j > i {
                        return Err(invalid(format!(
                            "stage {i} reads pipe {j} of a later stage"
                        )));
                    }
                    let name = r.name().to_string();
                    endpoints[i].push(r);
                    name
                }
            } else if a == OUT {
                if i + 1 != n || terminal == Terminal::Save || !matches!(sink, Sink::None) {
                    return Err(invalid(format!("stage {i} writes \\0out")));
                }
                let (reg, made) = match terminal {
                    Terminal::Count => {
                        let c = Shared::<RecordCount>::default();
                        (
                            Registered::output("<count()>", Box::new(c.clone())),
                            Sink::Count(c),
                        )
                    }
                    _ => {
                        let b = Shared::<Vec<u8>>::default();
                        (
                            Registered::output("<to_bytes()>", Box::new(b.clone())),
                            Sink::Bytes(b),
                        )
                    }
                };
                sink = made;
                let name = reg.name().to_string();
                endpoints[i].push(reg);
                name
            } else if a.starts_with(SRC) {
                let k =
                    source_index(a).ok_or_else(|| invalid(format!("bad source token {a:?}")))?;
                let spec = sources.remove(a).ok_or_else(|| {
                    invalid(format!(
                        "stage {i} reads source {a:?}, which is not in the plan's sources or is \
read twice"
                    ))
                })?;
                let (pump, reader, display) = sources::open(py, spec, k, i)?;
                let reg = Registered::input(&display, Box::new(reader));
                let name = reg.name().to_string();
                endpoints[i].push(reg);
                pumps.push(pump);
                name
            } else {
                a.clone()
            };
            argv.push(OsString::from(name));
        }
        argvs.push(argv);
    }
    if let Some((j, _)) = readers.iter().find(|(_, r)| r.is_some()) {
        return Err(invalid(format!("the pipe of stage {j} is never read")));
    }
    if terminal != Terminal::Save && matches!(sink, Sink::None) {
        return Err(invalid("the last stage does not write \\0out"));
    }
    if let Some(token) = sources.keys().next() {
        return Err(invalid(format!("source {token:?} is never read")));
    }

    let specs = argvs
        .into_iter()
        .zip(stages)
        .zip(endpoints.into_iter().zip(fails))
        .map(|((argv, s), (endpoints, fail_downstream))| StageSpec {
            argv,
            command_line: s.command_line.clone(),
            tmp_dir: None,
            endpoints,
            fail_downstream,
        })
        .collect();
    Ok((specs, sink, pumps))
}

/// How a run ended, as the driver saw it.
pub(crate) struct Driven {
    result: Result<bool, StageError>,
    /// An exception to raise instead of the result: a signal handler's (`KeyboardInterrupt`)
    /// or a log handler's, which cancelled the run.
    pending: Option<PyErr>,
    /// The last error record forwarded.
    logged: Option<String>,
}

/// A running set, dropped with the GIL released however the driver leaves (R2): its drop
/// cancels any stage still running and joins every stage thread, and a stage may need the GIL
/// to end (a Python-backed source or sink), even when the driver is unwinding from a panic.
pub(crate) struct Running<'py> {
    py: Python<'py>,
    set: Option<StageSet>,
}

impl<'py> Running<'py> {
    pub(crate) fn new(py: Python<'py>, set: StageSet) -> Running<'py> {
        Running { py, set: Some(set) }
    }
}

impl Drop for Running<'_> {
    fn drop(&mut self) {
        if let Some(set) = self.set.take() {
            self.py.detach(move || drop(set));
        }
    }
}

/// Wait, with the GIL released, until one of `senders` has room (or its reader is gone), at
/// most `timeout`.
fn wait_for_room(senders: &[Sender<Vec<u8>>], timeout: Duration) {
    let mut sel = Select::new();
    for s in senders {
        sel.send(s);
    }
    let _ = sel.ready_timeout(timeout);
}

/// Drive a running set to its end (D91): pump the sources, forward logs and check signals with
/// the GIL held, and wait with it released, [`TURN`] at a time: for room in a source's channel
/// while any source is being pumped, else for the stages.
pub(crate) fn drive(py: Python<'_>, set: StageSet, pumps: &mut [Pump]) -> Driven {
    let mut running = Running::new(py, set);
    let set = running.set.as_mut().expect("the set is taken only on drop");
    let mut pending: Option<PyErr> = None;
    let mut logged = None;
    let mut stop = |set: &StageSet, pumps: &mut [Pump], e: PyErr| {
        if pending.is_none() {
            set.cancel();
            // The sources' readers fail with the cancellation rather than see an end.
            pumps.iter_mut().for_each(Pump::abandon);
            pending = Some(e);
        }
    };
    let result = loop {
        if let Err(e) = logs::forward(py, &mut logged) {
            stop(set, pumps, e);
        }
        if let Err(e) = py.check_signals() {
            stop(set, pumps, e);
        }
        let mut room: Vec<Sender<Vec<u8>>> = Vec::new();
        let mut starved = false;
        for i in 0..pumps.len() {
            match pumps[i].pump(py) {
                Ok(Turn::Done) => {}
                Ok(Turn::Room | Turn::Full) => room.extend(pumps[i].sender()),
                Ok(Turn::Starved) => starved = true,
                Err(e) => stop(set, pumps, e),
            }
        }
        let r = if room.is_empty() && !starved {
            py.detach(|| set.wait(TURN))
        } else {
            py.detach(|| {
                wait_for_room(&room, if starved { STARVED } else { TURN });
                drop(room);
                set.wait(Duration::ZERO)
            })
        };
        if let Some(r) = r {
            break r;
        }
    };
    if let Err(e) = logs::forward(py, &mut logged) {
        stop(set, pumps, e);
    }
    // Every stage has ended, so dropping the set only joins.
    drop(running);
    Driven {
        result,
        pending,
        logged,
    }
}

impl Driven {
    /// The run's result as Python sees it: the pending exception; else, when the run failed on
    /// a Python-object source's failure as its stage read it, what the source raised (or the
    /// `DataError` for what it held); else the run's error as a gloci exception (`spec` gives
    /// each stage's command, for usage messages).
    pub(crate) fn settle(
        self,
        py: Python<'_>,
        spec: impl Fn(usize) -> Option<&'static CommandSpec>,
        pumps: &mut [Pump],
    ) -> PyResult<bool> {
        if let Some(e) = self.pending {
            return Err(e);
        }
        self.result.map_err(|e| {
            pumps
                .iter_mut()
                .find_map(|p| p.claim(&e))
                .unwrap_or_else(|| {
                    errors::from_stage(py, &e, spec(e.stage), self.logged.as_deref())
                })
        })
    }
}

/// Run a plan and return the terminal's result.
pub(crate) fn run_plan(py: Python<'_>, plan: &Bound<'_, PyDict>) -> PyResult<Py<PyAny>> {
    let (stages, terminal, sources) = read_plan(plan)?;
    logs::sync_level(py)?;
    let (specs, sink, mut pumps) = build(py, &stages, terminal, sources)?;
    let set = StageSet::spawn(specs);
    let driven = drive(py, set, &mut pumps);
    driven.settle(py, |i| spec_of(&stages[i].argv), &mut pumps)?;
    Ok(match sink {
        Sink::None => py.None(),
        Sink::Bytes(b) => {
            let bytes = std::mem::take(&mut *b.lock());
            PyBytes::new(py, &bytes).into_any().unbind()
        }
        Sink::Count(c) => c.lock().total().into_pyobject(py)?.into_any().unbind(),
    })
}
