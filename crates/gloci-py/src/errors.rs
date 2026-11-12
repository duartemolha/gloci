//! Python exceptions from a run's error (D90): the classes of `gloci.errors`, built from the
//! [`StageError`] the stage runtime judged while the run's in-memory names were registered, so
//! messages and attributes show display names, never a raw `\0gloci:` name (ruling R4).

use gloci::embed::{CommandSpec, ErrorSite, Outcome, StageError, rewrite_flags_as_kwargs};
use gloci_io::{IoError, Loc};
use pyo3::call::PyCallArgs;
use pyo3::exceptions::{PyKeyboardInterrupt, PyRuntimeError};
use pyo3::prelude::*;

/// The message of a run whose failure was logged rather than returned (`Outcome::Reported`),
/// when no error record was forwarded.
const LOGGED: &str = "the step's error was logged";

/// An instance of the class `name` of `gloci.errors`, made with `args`, as an error to raise.
fn raise<'py, A: PyCallArgs<'py>>(py: Python<'py>, name: &str, args: A) -> PyErr {
    let made = py
        .import("gloci.errors")
        .and_then(|m| m.getattr(name))
        .and_then(|cls| cls.call1(args));
    match made {
        Ok(e) => PyErr::from_value(e),
        // Only a broken installation gets here; say what was being raised.
        Err(e) => PyRuntimeError::new_err(format!("cannot raise gloci.errors.{name}: {e}")),
    }
}

/// `gloci.UsageError(msg)`.
pub(crate) fn usage(py: Python<'_>, msg: impl Into<String>) -> PyErr {
    raise(py, "UsageError", (msg.into(),))
}

/// `gloci.DataError(msg, path)`, with no line: an in-memory source's problem as a whole.
pub(crate) fn data_error(py: Python<'_>, msg: impl Into<String>, path: &str) -> PyErr {
    raise(py, "DataError", (msg.into(), path, None::<u64>))
}

/// `gloci.DataError(msg, path, line)`: a source's item or row `line` (1-based, D88).
pub(crate) fn data_at(py: Python<'_>, msg: impl Into<String>, path: &str, line: u64) -> PyErr {
    raise(py, "DataError", (msg.into(), path, line))
}

/// `gloci.GlociError(msg)`: an internal error.
pub(crate) fn internal(py: Python<'_>, msg: impl Into<String>) -> PyErr {
    raise(py, "GlociError", (msg.into(),))
}

/// The internal-error text of a panic in the host's own code (a stage's panic is turned into
/// the same text by the stage runtime).
pub(crate) fn panic_message(payload: &(dyn std::any::Any + Send)) -> String {
    let what = payload
        .downcast_ref::<&str>()
        .map(|s| s.to_string())
        .or_else(|| payload.downcast_ref::<String>().cloned())
        .unwrap_or_else(|| "a panic".into());
    IoError::Internal(format!("gloci-py panicked ({what})")).to_string()
}

fn number(loc: Loc) -> u64 {
    match loc {
        Loc::Line(n) | Loc::Record(n) => n,
    }
}

/// `gloci.DataError` or `gloci.UnsortedError` for a failure, with the site's fields.
fn data(py: Python<'_>, msg: &str, site: Option<&ErrorSite>) -> PyErr {
    match site {
        Some(ErrorSite {
            order: Some(o),
            path,
            ..
        }) => raise(
            py,
            "UnsortedError",
            (
                msg,
                path.as_str(),
                number(o.loc),
                (o.prev.as_str(), o.rec.as_str()),
                (number(o.prev_loc), number(o.loc)),
                o.hint.as_str(),
            ),
        ),
        Some(s) => raise(py, "DataError", (msg, s.path.as_str(), s.loc.map(number))),
        None => raise(py, "DataError", (msg, None::<&str>, None::<u64>)),
    }
}

/// Whether the error is an internal one: a panic, a lost thread, a bug.
fn is_internal(e: &anyhow::Error) -> bool {
    e.chain()
        .any(|c| matches!(c.downcast_ref::<IoError>(), Some(IoError::Internal(_))))
}

/// The exception for a failed run (D90):
/// - a usage problem is a `UsageError`, in keyword terms when the stage's command is known;
/// - an internal error (a panic) is a `GlociError`;
/// - any other failure is a `DataError`, or an `UnsortedError` for records out of order, with
///   `.path` and `.line` (and the order's fields) from the error's site;
/// - a failure the command logged is a `DataError` with the last error record forwarded
///   (`logged`);
/// - a cancelled run with no exception to raise again raises `KeyboardInterrupt`.
pub(crate) fn from_stage(
    py: Python<'_>,
    e: &StageError,
    spec: Option<&CommandSpec>,
    logged: Option<&str>,
) -> PyErr {
    match &e.outcome {
        Outcome::Usage(msg) => usage(
            py,
            spec.map_or_else(|| msg.clone(), |s| rewrite_flags_as_kwargs(msg, s)),
        ),
        Outcome::Failure(msg) if is_internal(&e.error) => internal(py, msg.as_str()),
        Outcome::Failure(msg) => data(py, msg, e.site.as_ref()),
        Outcome::Reported => data(py, logged.unwrap_or(LOGGED), None),
        Outcome::Cancelled => PyKeyboardInterrupt::new_err(()),
        // A broken pipe is a quiet end, never a set's error; report it rather than hide it.
        Outcome::BrokenPipe => internal(py, "a step ended with a broken pipe"),
    }
}
