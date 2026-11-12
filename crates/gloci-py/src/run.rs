//! Running a pipeline's plan (D85, D91): its stages on threads joined by bounded pipes
//! ([`StageSet`]), driven from the calling Python thread.
//!
//! A plan is `{"stages": [{"argv", "command_line"}], "sources": {token: obj},
//! "terminal": "save"|"bytes"|"count", "names"}`. Each stage's argv is its whole command line,
//! `gloci` first, with these tokens where paths go, which become in-memory names (D84):
//! - `"\0pipe:<j>"`: the pipe stage `j` writes (its `-o`), read by exactly one later stage;
//! - `"\0out"`: the last stage's output, collected for the terminal (`bytes`, `count`);
//! - `"\0src:<k>"`: the Python-object source `k` of `"sources"` (not taken yet).
//!
//! **The driver loop** holds the GIL only to forward log records and check signals, and
//! releases it while it waits for the stages, at most [`TURN`] at a time. A signal handler's
//! exception (`KeyboardInterrupt`) cancels the run, which removes its partial output files, and
//! is raised once every stage has ended. The stage set is dropped with the GIL released (its
//! drop joins the stage threads, R2).

use crate::{errors, logs};
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

fn read_plan(plan: &Bound<'_, PyDict>) -> PyResult<(Vec<PlanStage>, Terminal)> {
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
    if let Some(sources) = plan.get_item("sources")?
        && !sources.is_none()
        && sources.len()? > 0
    {
        return Err(invalid(
            "Python-object sources are not supported yet; pass paths",
        ));
    }
    Ok((stages, terminal))
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

/// Counts the records of a text output as it is written: its lines, but blank ones, `#`
/// header and comment lines, `track` and `browser` lines, and SAM's `@XY` header lines before
/// the first record.
#[derive(Default)]
struct RecordCount {
    records: u64,
    /// The first bytes of the line being written, enough to tell a header line.
    head: [u8; 8],
    head_len: usize,
    /// Bytes of the current line have been written.
    in_line: bool,
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
        let word = |w: &[u8]| {
            head.starts_with(w) && matches!(head.get(w.len()), None | Some(b' ' | b'\t' | b'\r'))
        };
        let sam_header = self.records == 0
            && head.len() >= 3
            && head[0] == b'@'
            && head[1..3].iter().all(u8::is_ascii_alphabetic)
            && matches!(head.get(3), None | Some(b'\t' | b'\r'));
        let blank = !self.in_line || head == b"\r";
        if !(blank || head[0] == b'#' || word(b"track") || word(b"browser") || sam_header) {
            self.records += 1;
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
/// the stage that uses it, and the terminal's sink.
fn build(stages: &[PlanStage], terminal: Terminal) -> PyResult<(Vec<StageSpec>, Sink)> {
    let n = stages.len();
    let mut endpoints: Vec<Vec<Registered>> = (0..n).map(|_| Vec::new()).collect();
    let mut fails: Vec<Vec<PipeFail>> = (0..n).map(|_| Vec::new()).collect();
    // Pipe j: the name of its writing end, and its reading end until a stage takes it.
    let mut writers: BTreeMap<usize, String> = BTreeMap::new();
    let mut readers: BTreeMap<usize, Option<Registered>> = BTreeMap::new();
    let mut sink = Sink::None;
    let mut argvs: Vec<Vec<OsString>> = Vec::with_capacity(n);

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
                return Err(invalid(format!(
                    "source {a:?} is not in the plan's sources"
                )));
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
    Ok((specs, sink))
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

/// Drive a running set to its end (D91): forward logs and check signals with the GIL held,
/// and wait for the stages with it released, [`TURN`] at a time.
pub(crate) fn drive(py: Python<'_>, set: StageSet) -> Driven {
    let mut running = Running::new(py, set);
    let set = running.set.as_mut().expect("the set is taken only on drop");
    let mut pending: Option<PyErr> = None;
    let mut logged = None;
    let mut stop = |set: &StageSet, e: PyErr| {
        if pending.is_none() {
            set.cancel();
            pending = Some(e);
        }
    };
    let result = loop {
        if let Err(e) = logs::forward(py, &mut logged) {
            stop(set, e);
        }
        if let Err(e) = py.check_signals() {
            stop(set, e);
        }
        if let Some(r) = py.detach(|| set.wait(TURN)) {
            break r;
        }
    };
    if let Err(e) = logs::forward(py, &mut logged) {
        stop(set, e);
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
    /// The run's result as Python sees it: the pending exception, else the run's error as a
    /// gloci exception (`spec` gives each stage's command, for usage messages).
    pub(crate) fn settle(
        self,
        py: Python<'_>,
        spec: impl Fn(usize) -> Option<&'static CommandSpec>,
    ) -> PyResult<bool> {
        if let Some(e) = self.pending {
            return Err(e);
        }
        self.result
            .map_err(|e| errors::from_stage(py, &e, spec(e.stage), self.logged.as_deref()))
    }
}

/// Run a plan and return the terminal's result.
pub(crate) fn run_plan(py: Python<'_>, plan: &Bound<'_, PyDict>) -> PyResult<Py<PyAny>> {
    let (stages, terminal) = read_plan(plan)?;
    logs::sync_level(py)?;
    let (specs, sink) = build(&stages, terminal)?;
    let set = StageSet::spawn(specs);
    drive(py, set).settle(py, |i| spec_of(&stages[i].argv))?;
    Ok(match sink {
        Sink::None => py.None(),
        Sink::Bytes(b) => {
            let bytes = std::mem::take(&mut *b.lock());
            PyBytes::new(py, &bytes).into_any().unbind()
        }
        Sink::Count(c) => c.lock().total().into_pyobject(py)?.into_any().unbind(),
    })
}
