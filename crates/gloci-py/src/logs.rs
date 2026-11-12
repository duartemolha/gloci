//! Log forwarding (D91): gloci's warnings and progress go to the Python logger `gloci`.
//!
//! The `log` crate's logger is process-wide, and its records are made on the run's threads,
//! which never take the GIL for them. Each record is rendered at once (while the run's
//! in-memory names are registered, so it shows their display names) and queued; a run's driver
//! loop forwards the queue to `logging.getLogger("gloci")` on the calling thread, at each turn
//! and once more at the end. Runs share the queue, so with runs on several Python threads at
//! once, a record may be forwarded by another run's driver; it still reaches the same logger.

use log::{Level, LevelFilter, Log, Metadata, Record};
use pyo3::prelude::*;
use std::sync::{Mutex, MutexGuard, PoisonError};

/// The name of the Python logger.
const LOGGER: &str = "gloci";

static QUEUE: Mutex<Vec<(Level, String)>> = Mutex::new(Vec::new());

fn queue() -> MutexGuard<'static, Vec<(Level, String)>> {
    QUEUE.lock().unwrap_or_else(PoisonError::into_inner)
}

struct Forward;

impl Log for Forward {
    fn enabled(&self, m: &Metadata) -> bool {
        m.level() <= log::max_level()
    }

    fn log(&self, r: &Record) {
        if self.enabled(r.metadata()) {
            let text = gloci::logging::render(r);
            queue().push((r.level(), text));
        }
    }

    fn flush(&self) {}
}

static FORWARD: Forward = Forward;

/// Install the forwarding logger, with nothing enabled until a run sets the level.
pub(crate) fn install() {
    // The extension has its own copy of `log`, so no other logger can be set in it.
    let _ = log::set_logger(&FORWARD);
    log::set_max_level(LevelFilter::Off);
}

fn logger(py: Python<'_>) -> PyResult<Bound<'_, PyAny>> {
    py.import("logging")?.call_method1("getLogger", (LOGGER,))
}

/// The `log` level that makes the records the Python level lets through.
fn filter_for(python_level: i64) -> LevelFilter {
    match python_level {
        ..=5 => LevelFilter::Trace,
        6..=10 => LevelFilter::Debug,
        11..=20 => LevelFilter::Info,
        21..=30 => LevelFilter::Warn,
        31..=40 => LevelFilter::Error,
        _ => LevelFilter::Off,
    }
}

fn python_level(level: Level) -> i64 {
    match level {
        Level::Error => 40,
        Level::Warn => 30,
        Level::Info => 20,
        Level::Debug => 10,
        Level::Trace => 5,
    }
}

/// Make only the records the logger `gloci` would handle: its effective level, read at the
/// start of each run (the CLI's `-q` and `-v` are this level in Python, D86).
pub(crate) fn sync_level(py: Python<'_>) -> PyResult<()> {
    let level: i64 = logger(py)?.call_method0("getEffectiveLevel")?.extract()?;
    log::set_max_level(filter_for(level));
    Ok(())
}

/// Forward the queued records to the logger `gloci`; `last_error` keeps the text of the last
/// error record (the message of a failure a command logged rather than returned). An exception
/// from a handler is returned after every record has been offered.
pub(crate) fn forward(py: Python<'_>, last_error: &mut Option<String>) -> PyResult<()> {
    let records = std::mem::take(&mut *queue());
    if records.is_empty() {
        return Ok(());
    }
    let logger = logger(py)?;
    let mut failed = None;
    for (level, text) in records {
        if level == Level::Error {
            *last_error = Some(text.clone());
        }
        if let Err(e) = logger.call_method1("log", (python_level(level), text)) {
            failed.get_or_insert(e);
        }
    }
    failed.map_or(Ok(()), Err)
}
