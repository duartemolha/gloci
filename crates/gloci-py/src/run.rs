//! Running a pipeline's plan (D85, D91): its stages on threads joined by bounded pipes
//! ([`StageSet`]), driven from the calling Python thread.
//!
//! A plan is `{"stages": [{"argv", "command_line"}], "sources": {token: source},
//! "terminal": "save"|"bytes"|"count"|"iter"|"arrow", "names"}`. Each stage's argv is its whole
//! command line, `gloci` first, with these tokens where paths go, which become in-memory names
//! (D84):
//! - `"\0pipe:<j>"`: the pipe stage `j` writes (its `-o`), read by exactly one later stage; a
//!   chain's steps, and a nested pipeline (a Pipeline given as a source), are such stages;
//! - `"\0out"`: the last stage's output, for the terminal: its bytes (`bytes`), or its records
//!   (`count`, `iter`, `arrow`), which a host stage after the commands reads back
//!   ([`records::reread`]) as the step's [`embed::read_back`] says;
//! - `"\0src:<k>"`: the Python-object source of `"sources"` with that token
//!   (`{"kind", "obj", "display"}`, [`sources`]), read by exactly one stage.
//!
//! **The driver loop** ([`Driver`]) holds the GIL only to pump the Python-object sources (at
//! most [`sources::TURN_BYTES`] each a turn, and only into room their channel has), forward log
//! records and check signals. It releases the GIL while it waits, at most [`TURN`] at a time,
//! for room in a source's channel, for a batch of records (iteration), or for the stages to
//! end. A signal handler's exception (`KeyboardInterrupt`) cancels the run, which removes its
//! partial output files, and is raised once every stage has ended. The stage set is dropped
//! with the GIL released (its drop joins the stage threads).

use crate::arrow::{ArrowSink, ArrowStream, Table};
use crate::records::{self, Batch, CountSink, IterSink, RecordIterator, RowSink};
use crate::sources::{self, Pump, SourceSpec, Turn};
use crate::{errors, logs};
use crossbeam_channel::{Receiver, Select, Sender, bounded};
use gloci::embed::{
    self, CommandSpec, HostStage, PIPE_CHUNK, PIPE_CHUNKS, PipeFail, PipeReader, StageError,
    StageSet, StageSpec,
};
use gloci_io::memio::Registered;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict};
use std::collections::BTreeMap;
use std::ffi::OsString;
use std::io::{self, Write};
use std::path::PathBuf;
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
    /// `save(path)`, and the module functions: the stages write files; the result is the list
    /// of the output files the run committed, as the commands named them.
    Save,
    /// `to_bytes()`: the last stage's output.
    Bytes,
    /// `count()`: the number of records it reads back as (counted in Rust).
    Count,
    /// Iteration: its records, handed to Python as they come ([`RecordIterator`]).
    Iter,
    /// `to_arrow_stream()` and the frame terminals: its records as an Arrow table.
    Arrow,
}

impl Terminal {
    /// Whether the terminal reads the output back as records.
    fn reads_back(self) -> bool {
        matches!(self, Terminal::Count | Terminal::Iter | Terminal::Arrow)
    }
}

struct PlanStage {
    argv: Vec<String>,
    command_line: String,
}

fn invalid(msg: impl Into<String>) -> PyErr {
    PyValueError::new_err(format!("invalid gloci plan: {}", msg.into()))
}

/// A plan's parts: its stages, its terminal, its Python-object sources by token, and the
/// column names of an Arrow table (`names=`).
type Plan = (
    Vec<PlanStage>,
    Terminal,
    BTreeMap<String, SourceSpec>,
    Option<Vec<String>>,
);

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
        "iter" => Terminal::Iter,
        "arrow" => Terminal::Arrow,
        other => return Err(invalid(format!("unknown terminal {other:?}"))),
    };
    let names: Option<Vec<String>> = match plan.get_item("names")? {
        Some(n) if !n.is_none() => Some(n.extract()?),
        _ => None,
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
    Ok((stages, terminal, sources, names))
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

/// Where the last stage's output goes.
enum Sink {
    None,
    Bytes(Shared<Vec<u8>>),
    /// A pipe a host stage reads the records from: its reading end, and its name.
    Records(PipeReader, String),
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
        // The one fail handle goes to the writing stage, and no copy is kept.
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
                let reg = if terminal.reads_back() {
                    // A pipe to the host stage that reads the records back; the one fail
                    // handle goes to this stage.
                    let display = format!(
                        "<{} output>",
                        spec_of(&s.argv).map_or_else(|| "step".into(), embed::method_name)
                    );
                    let (w, r) = embed::pipe(PIPE_CHUNKS, PIPE_CHUNK);
                    fails[i].push(w.fail_handle());
                    sink = Sink::Records(r, display.clone());
                    Registered::output(&display, Box::new(w))
                } else {
                    let b = Shared::<Vec<u8>>::default();
                    sink = Sink::Bytes(b.clone());
                    Registered::output("<to_bytes()>", Box::new(b))
                };
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

/// A running set and its driver loop (D91), from its start to its result: [`run_plan`]
/// drives it to the end at once, a [`RecordIterator`] a turn at a time as records are asked
/// for.
///
/// However it is left (the end, an error, an iterator dropped, a panic unwinding), the set is
/// dropped with the GIL released: its drop cancels any stage still running and joins
/// every stage thread, and a stage may need the GIL to end.
pub(crate) struct Driver {
    set: Option<StageSet>,
    pumps: Vec<Pump>,
    /// An exception to raise instead of the result: a signal handler's (`KeyboardInterrupt`)
    /// or a log handler's, which cancelled the run.
    pending: Option<PyErr>,
    /// The last error record forwarded.
    logged: Option<String>,
    /// The set's result, once every stage has ended.
    result: Option<Result<bool, StageError>>,
    /// The output files the run committed ([`StageSet::outputs`]).
    outputs: Vec<PathBuf>,
}

/// Wait, with the GIL released, until one of `senders` has room (or its reader is gone), or
/// `ready` has something (or its sender is gone), at most `timeout`.
fn wait_for(senders: &[Sender<Vec<u8>>], ready: Option<&Receiver<Batch>>, timeout: Duration) {
    let mut sel = Select::new();
    for s in senders {
        sel.send(s);
    }
    if let Some(r) = ready {
        sel.recv(r);
    }
    let _ = sel.ready_timeout(timeout);
}

impl Driver {
    pub(crate) fn new(set: StageSet, pumps: Vec<Pump>) -> Driver {
        Driver {
            set: Some(set),
            pumps,
            pending: None,
            logged: None,
            result: None,
            outputs: Vec::new(),
        }
    }

    /// Every stage has ended.
    pub(crate) fn ended(&self) -> bool {
        self.result.is_some()
    }

    /// A signal or log handler's exception cancelled the run.
    pub(crate) fn interrupted(&self) -> bool {
        self.pending.is_some()
    }

    fn stop(&mut self, e: PyErr) {
        if self.pending.is_none() {
            if let Some(set) = &self.set {
                set.cancel();
            }
            // The sources' readers fail with the cancellation rather than see an end.
            self.pumps.iter_mut().for_each(Pump::abandon);
            self.pending = Some(e);
        }
    }

    /// One turn of the loop: forward logs, check signals and pump the sources with the GIL
    /// held, then wait with it released, at most [`TURN`]: for room in a source's channel while
    /// any is being pumped, or for a batch on `ready`, else for the stages. True once every
    /// stage has ended.
    pub(crate) fn turn(&mut self, py: Python<'_>, ready: Option<&Receiver<Batch>>) -> bool {
        if self.result.is_some() {
            return true;
        }
        if let Err(e) = logs::forward(py, &mut self.logged) {
            self.stop(e);
        }
        if let Err(e) = py.check_signals() {
            self.stop(e);
        }
        let mut room: Vec<Sender<Vec<u8>>> = Vec::new();
        let mut starved = false;
        for i in 0..self.pumps.len() {
            match self.pumps[i].pump(py) {
                Ok(Turn::Done) => {}
                Ok(Turn::Room | Turn::Full) => room.extend(self.pumps[i].sender()),
                Ok(Turn::Starved) => starved = true,
                Err(e) => self.stop(e),
            }
        }
        let set = self
            .set
            .as_mut()
            .expect("the set is taken only once it has ended");
        let r = if room.is_empty() && !starved && ready.is_none() {
            py.detach(|| set.wait(TURN))
        } else {
            py.detach(|| {
                wait_for(&room, ready, if starved { STARVED } else { TURN });
                drop(room);
                set.wait(Duration::ZERO)
            })
        };
        let Some(r) = r else {
            return false;
        };
        self.outputs = set.outputs();
        self.result = Some(r);
        if let Err(e) = logs::forward(py, &mut self.logged) {
            self.stop(e);
        }
        // Every stage has ended, so dropping the set only joins.
        if let Some(set) = self.set.take() {
            py.detach(move || drop(set));
        }
        true
    }

    /// Drive the run to its end.
    pub(crate) fn finish(&mut self, py: Python<'_>) {
        while !self.turn(py, None) {}
    }

    /// The run's result as Python sees it, once it has ended: the pending exception; else, when
    /// the run failed on a Python-object source's failure as its stage read it, what the source
    /// raised (or the `DataError` for what it held); else the run's error as a gloci exception
    /// (`spec` gives each stage's command, for usage messages).
    pub(crate) fn settle(
        &mut self,
        py: Python<'_>,
        spec: impl Fn(usize) -> Option<&'static CommandSpec>,
    ) -> PyResult<bool> {
        let result = self.result.take().expect("settled once the run has ended");
        if let Some(e) = self.pending.take() {
            return Err(e);
        }
        let logged = self.logged.as_deref();
        result.map_err(|e| {
            self.pumps
                .iter_mut()
                .find_map(|p| p.claim(&e))
                .unwrap_or_else(|| errors::from_stage(py, &e, spec(e.stage), logged))
        })
    }
}

impl Drop for Driver {
    fn drop(&mut self) {
        let Some(set) = self.set.take() else {
            return;
        };
        set.cancel();
        self.pumps.iter_mut().for_each(Pump::abandon);
        // With the GIL released. `try_attach` runs the closure unless the interpreter is shutting
        // down; then the closure, and the set in it, is dropped here (the pumps are abandoned,
        // so no stage needs the GIL to end).
        let _ = Python::try_attach(move |py| py.detach(move || drop(set)));
    }
}

/// The spec of stage `i`'s command; none for a host stage.
fn stage_spec(stages: &[PlanStage], i: usize) -> Option<&'static CommandSpec> {
    stages.get(i).and_then(|s| spec_of(&s.argv))
}

/// The host stage that reads the output `r` back into `sink`.
fn host<S: RowSink + Send + 'static>(
    r: PipeReader,
    display: String,
    back: embed::ReadBack,
    mut sink: S,
    done: impl FnOnce(S) + Send + 'static,
) -> HostStage {
    HostStage {
        run: Box::new(move || {
            records::reread(Box::new(r), &display, back, &mut sink)?;
            done(sink);
            Ok(true)
        }),
    }
}

/// Run a plan and return the terminal's result: the files written (`save`), the bytes, the
/// count, a [`RecordIterator`] (which drives the run as it is read) or an [`ArrowStream`].
pub(crate) fn run_plan(py: Python<'_>, plan: &Bound<'_, PyDict>) -> PyResult<Py<PyAny>> {
    let (stages, terminal, sources, names) = read_plan(plan)?;
    logs::sync_level(py)?;
    let (specs, sink, pumps) = build(py, &stages, terminal, sources)?;
    let last: Vec<OsString> = specs.last().map(|s| s.argv.clone()).unwrap_or_default();
    let reread = |r: PipeReader, display: String| -> PyResult<_> {
        let back = embed::read_back(&last).map_err(|e| errors::internal(py, e.0))?;
        Ok((r, display, back))
    };
    let count = Arc::new(Mutex::new(0u64));
    let table: Arc<Mutex<Option<Table>>> = Arc::default();
    let mut batches = None;
    let mut hosts = Vec::new();
    let sink = match sink {
        Sink::Records(r, display) => {
            let (r, display, back) = reread(r, display)?;
            hosts.push(match terminal {
                Terminal::Iter => {
                    // One batch waiting, one being filled: about 8 MiB besides the one read.
                    let (tx, rx) = bounded(1);
                    batches = Some(rx);
                    host(r, display, back, IterSink::new(tx), drop)
                }
                Terminal::Arrow => host(r, display, back, ArrowSink::new(table.clone()), drop),
                _ => {
                    let count = count.clone();
                    host(r, display, back, CountSink::default(), move |c| {
                        *count.lock().unwrap_or_else(PoisonError::into_inner) = c.0;
                    })
                }
            });
            None
        }
        other => Some(other),
    };
    let set = StageSet::spawn_with_hosts(specs, hosts);
    let mut driver = Driver::new(set, pumps);
    if let Some(rx) = batches {
        let argvs: Vec<Vec<String>> = stages.iter().map(|s| s.argv.clone()).collect();
        let it = RecordIterator::new(rx, driver, &argvs);
        return Ok(Bound::new(py, it)?.into_any().unbind());
    }
    driver.finish(py);
    let outputs = std::mem::take(&mut driver.outputs);
    driver.settle(py, |i| stage_spec(&stages, i))?;
    drop(driver);
    Ok(match (sink, terminal) {
        (Some(Sink::None), _) => {
            let paths: Vec<std::ffi::OsString> =
                outputs.into_iter().map(|p| p.into_os_string()).collect();
            paths.into_pyobject(py)?.into_any().unbind()
        }
        (Some(Sink::Bytes(b)), _) => {
            let bytes = std::mem::take(&mut *b.lock());
            PyBytes::new(py, &bytes).into_any().unbind()
        }
        (_, Terminal::Arrow) => {
            let mut built = table
                .lock()
                .unwrap_or_else(PoisonError::into_inner)
                .take()
                .ok_or_else(|| errors::internal(py, "the Arrow table was not built"))?;
            if let Some(names) = names {
                built = built.renamed(&names).map_err(|e| errors::usage(py, e))?;
            }
            Bound::new(py, ArrowStream::new(built))?.into_any().unbind()
        }
        _ => {
            let n = *count.lock().unwrap_or_else(PoisonError::into_inner);
            n.into_pyobject(py)?.into_any().unbind()
        }
    })
}
