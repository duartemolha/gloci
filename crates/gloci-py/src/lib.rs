//! `gloci._native`: the Python library's glue over [`gloci::embed`] (D83). It holds no command
//! logic: the commands, their arguments and how a step takes its input come from the CLI's own
//! definition, and every step runs through the CLI's parsing and dispatch.
//!
//! - [`command_specs`] and [`bindings`] describe the commands and the step methods.
//! - [`kwargs_to_argv`] checks a step's keyword arguments when the step is called (D85).
//! - [`step_argv`] and [`command_line`] build a stage's argv and its shell equivalent.
//! - [`run_plan`] runs a pipeline ([`run`]).

#[cfg(windows)]
compile_error!(
    "gloci does not build on Windows: use WSL (the Windows Subsystem for Linux) and install \
     gloci there with `pip install gloci`"
);

mod errors;
mod logs;
mod run;

use gloci::embed::{
    self, ArgSpec, Bind, Binding, CommandSpec, KwValue, Role, rewrite_flags_as_kwargs,
};
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyDict, PyList, PyTuple};
use std::ffi::OsString;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::sync::OnceLock;

/// The spec of the command at `path`, or a `UsageError`.
fn spec(py: Python<'_>, path: &[String]) -> PyResult<&'static CommandSpec> {
    embed::find_spec(path)
        .ok_or_else(|| errors::usage(py, format!("gloci {}: no such command", path.join(" "))))
}

fn arg_dict<'py>(py: Python<'py>, a: &ArgSpec) -> PyResult<Bound<'py, PyDict>> {
    let d = PyDict::new(py);
    d.set_item("id", &a.id)?;
    d.set_item("kwarg", a.kwarg())?;
    d.set_item("long", &a.long)?;
    d.set_item("short", a.short.map(String::from))?;
    d.set_item("aliases", &a.aliases)?;
    d.set_item("takes_value", a.takes_value)?;
    d.set_item("multiple", a.multiple)?;
    d.set_item("positional", a.positional)?;
    d.set_item("required", a.required)?;
    d.set_item("require_equals", a.require_equals)?;
    d.set_item("value_names", &a.value_names)?;
    d.set_item("possible_values", &a.possible_values)?;
    d.set_item("help", &a.help)?;
    Ok(d)
}

/// Every command that runs, with its arguments: `[{"path", "about", "args": [...]}]`.
#[pyfunction]
fn command_specs(py: Python<'_>) -> PyResult<Bound<'_, PyList>> {
    let out = PyList::empty(py);
    for s in embed::command_specs() {
        let d = PyDict::new(py);
        d.set_item("path", &s.path)?;
        d.set_item("about", &s.about)?;
        let args = PyList::empty(py);
        for a in &s.args {
            args.append(arg_dict(py, a)?)?;
        }
        d.set_item("args", args)?;
        out.append(d)?;
    }
    Ok(out)
}

/// The keyword argument that sets the format of a step's input (`a_format`, `format`), if the
/// command has one.
fn input_format(spec: &CommandSpec, b: &Binding) -> Option<String> {
    let name = match b.input {
        Bind::Flag(_) => format!("{}_format", embed::input_arg(spec, b)?.kwarg()),
        Bind::Positional { .. } | Bind::PositionalAfter(_) => "format".into(),
    };
    let a = spec.arg(&name)?;
    (!a.positional && a.value_names.iter().any(|v| v == "FMT")).then(|| a.kwarg())
}

/// The step methods (D87): `[{"method", "path", "input", "flag", "nary", "needs_file",
/// "input_format", "threads", "outputs"}]`.
/// - `input` is how the step takes the pipeline: `"flag"` (as the flag `flag`), `"positional"`
///   (its first positional; `nary` when it takes extra positional sources) or `"after"` (the
///   positional after `flag`).
/// - `input_format` is the keyword argument for the input's format, or None.
/// - `threads` says whether the command has `-t`; `outputs` lists the terminal options it has,
///   of `output`, `output_format`, `write_index` and `compression_level`.
#[pyfunction]
fn bindings(py: Python<'_>) -> PyResult<Bound<'_, PyList>> {
    let out = PyList::empty(py);
    for b in embed::BINDINGS {
        let spec = embed::find_spec(b.path)
            .ok_or_else(|| errors::internal(py, format!("no spec for {}", b.path.join(" "))))?;
        let has = |id: &str| spec.args.iter().any(|a| a.id == id);
        let (input, flag, nary) = match b.input {
            Bind::Flag(f) => ("flag", Some(f), false),
            Bind::Positional { nary } => ("positional", None, nary),
            Bind::PositionalAfter(p) => ("after", Some(p), false),
        };
        let d = PyDict::new(py);
        d.set_item("method", b.method)?;
        d.set_item("path", b.path)?;
        d.set_item("about", &spec.about)?;
        d.set_item("input", input)?;
        d.set_item("flag", flag)?;
        d.set_item("nary", nary)?;
        d.set_item("needs_file", b.needs_file)?;
        d.set_item("input_format", input_format(spec, b))?;
        d.set_item("threads", has("threads"))?;
        let outputs: Vec<&str> = [
            "output",
            "output_format",
            "write_index",
            "compression_level",
        ]
        .into_iter()
        .filter(|id| has(id))
        .collect();
        d.set_item("outputs", outputs)?;
        out.append(d)?;
    }
    Ok(out)
}

/// A keyword argument's value as text: a path-like object by its path, anything else by
/// `str()` (D86).
fn text(py: Python<'_>, name: &str, v: &Bound<'_, PyAny>) -> PyResult<String> {
    let s: String = if v.hasattr("__fspath__")? {
        let os = py.import("os")?;
        os.call_method1("fsdecode", (os.call_method1("fspath", (v,))?,))?
            .extract()?
    } else {
        v.str()?.extract()?
    };
    if s.contains('\0') {
        return Err(errors::usage(
            py,
            format!("{name}= holds a NUL character, which no argument can hold"),
        ));
    }
    Ok(s)
}

/// A keyword argument's value (D86): `True` gives the flag, `False` and `None` nothing, a list
/// or tuple the flag once per item, anything else its text.
fn kw_value(py: Python<'_>, name: &str, v: &Bound<'_, PyAny>) -> PyResult<KwValue> {
    if v.is_none() {
        return Ok(KwValue::Absent);
    }
    if let Ok(b) = v.cast::<PyBool>() {
        return Ok(KwValue::Flag(b.is_true()));
    }
    if v.is_instance_of::<PyList>() || v.is_instance_of::<PyTuple>() {
        let items = v
            .try_iter()?
            .map(|i| text(py, name, &i?))
            .collect::<PyResult<_>>()?;
        return Ok(KwValue::Values(items));
    }
    Ok(KwValue::Value(text(py, name, v)?))
}

fn kwargs(py: Python<'_>, kw: Option<&Bound<'_, PyDict>>) -> PyResult<Vec<(String, KwValue)>> {
    let Some(kw) = kw else {
        return Ok(Vec::new());
    };
    kw.iter()
        .map(|(k, v)| {
            let name: String = k.extract()?;
            let value = kw_value(py, &name, &v)?;
            Ok((name, value))
        })
        .collect()
}

fn strings(argv: Vec<OsString>) -> Vec<String> {
    argv.into_iter()
        .map(|a| a.to_string_lossy().into_owned())
        .collect()
}

/// Parse `argv` with the CLI's definition; its error is a `UsageError` in keyword terms.
fn check_parses(py: Python<'_>, spec: &CommandSpec, argv: &[OsString]) -> PyResult<()> {
    embed::parse(argv)
        .map(drop)
        .map_err(|e| errors::usage(py, rewrite_flags_as_kwargs(&e.0, spec)))
}

/// The argv of a step's keyword arguments (D86), checked as the step is called (D85): each
/// must be one of the command's, and clap must accept them with an input in place, so a bad
/// name or value is a `UsageError` at once, in keyword terms. For a module function's command
/// (`["vcf", "stats"]`), only the names and kinds are checked.
#[pyfunction]
#[pyo3(signature = (path, kwargs = None))]
fn kwargs_to_argv(
    py: Python<'_>,
    path: Vec<String>,
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Vec<String>> {
    let spec = spec(py, &path)?;
    let kw = self::kwargs(py, kwargs)?;
    let usage = |e: gloci_io::UsageError| errors::usage(py, e.0);
    let Some(b) = embed::binding(&path) else {
        if embed::module_function(&path).is_none() {
            return Err(errors::usage(
                py,
                format!("gloci {} is not offered in Python", path.join(" ")),
            ));
        }
        return embed::kwargs_to_argv_as(spec, &kw, Role::Function)
            .map(strings)
            .map_err(usage);
    };
    let argv = embed::kwargs_to_argv(spec, &kw).map_err(usage)?;
    let whole = embed::step_argv(b, &[OsString::from("-")], &kw, &[]).map_err(usage)?;
    check_parses(py, spec, &whole)?;
    Ok(strings(argv))
}

/// A stage's whole argv (`gloci` first, [`embed::step_argv`]): `inputs` where the binding puts
/// them, the keyword arguments, and the terminal's options (`-o`, `-t`, ...). With `check`,
/// clap must accept it, or a `UsageError` in keyword terms is raised.
#[pyfunction]
#[pyo3(signature = (path, inputs, kwargs = None, terminal = Vec::new(), check = true))]
fn step_argv(
    py: Python<'_>,
    path: Vec<String>,
    inputs: Vec<String>,
    kwargs: Option<&Bound<'_, PyDict>>,
    terminal: Vec<String>,
    check: bool,
) -> PyResult<Vec<String>> {
    let spec = spec(py, &path)?;
    let b = embed::binding(&path)
        .ok_or_else(|| errors::usage(py, format!("gloci {} is not a step", path.join(" "))))?;
    let kw = self::kwargs(py, kwargs)?;
    let os = |v: Vec<String>| v.into_iter().map(OsString::from).collect::<Vec<_>>();
    let argv =
        embed::step_argv(b, &os(inputs), &kw, &os(terminal)).map_err(|e| errors::usage(py, e.0))?;
    if check {
        check_parses(py, spec, &argv)?;
    }
    Ok(strings(argv))
}

/// The clap definition of the command at `path`, built.
fn clap_command(path: &[String]) -> Option<&'static clap::Command> {
    static ROOT: OnceLock<clap::Command> = OnceLock::new();
    let mut cmd = ROOT.get_or_init(|| {
        let mut root = gloci::command();
        root.build();
        root
    });
    for p in path {
        cmd = cmd.find_subcommand(p)?;
    }
    Some(cmd)
}

/// Whether the `--` of a step's argv can be left out of its shell form: no positional after it
/// looks like a flag, and the last option before it takes no more values than it was given, so
/// clap reads the positionals as positionals without it.
fn separator_needed(path: &[String], before: &[String], after: &[String]) -> bool {
    if after.iter().any(|a| a.starts_with('-') && a != "-") {
        return true;
    }
    let Some(cmd) = clap_command(path) else {
        return true;
    };
    let option = |a: &str| -> Option<(&clap::Arg, bool)> {
        if a == "-" || !a.starts_with('-') {
            return None;
        }
        let (name, attached) = match a.split_once('=') {
            Some((n, _)) => (n, true),
            None => (a, false),
        };
        let arg = if let Some(long) = name.strip_prefix("--") {
            cmd.get_arguments().find(|x| {
                x.get_long() == Some(long)
                    || x.get_all_aliases().is_some_and(|al| al.contains(&long))
            })
        } else {
            let mut c = name[1..].chars();
            match (c.next(), c.next()) {
                (Some(s), None) => cmd.get_arguments().find(|x| x.get_short() == Some(s)),
                _ => None,
            }
        };
        arg.map(|x| (x, attached))
    };
    let Some((at, (arg, attached))) = before
        .iter()
        .enumerate()
        .rev()
        .find_map(|(i, a)| option(a).map(|o| (i, o)))
    else {
        return false;
    };
    let given = before.len() - at - 1 + usize::from(attached);
    let max = arg
        .get_num_args()
        .map_or(usize::from(arg.get_action().takes_values()), |r| {
            r.max_values()
        });
    given < max
}

fn is_placeholder(a: &str) -> bool {
    a.strip_prefix("<python object ")
        .and_then(|r| r.strip_suffix('>'))
        .is_some_and(|n| !n.is_empty() && n.bytes().all(|c| c.is_ascii_digit()))
}

/// The shell form of a step's argv (D84): each word quoted as a shell needs it, except the
/// `<python object N>` that stands for an in-memory input, and without the `--` that
/// [`step_argv`] puts before positionals, where the shell form reads the same without it.
#[pyfunction]
fn command_line(path: Vec<String>, argv: Vec<String>) -> String {
    let mut words = argv;
    if let Some(sep) = words.iter().position(|a| a == "--")
        && !separator_needed(&path, &words[..sep], &words[sep + 1..])
    {
        words.remove(sep);
    }
    words
        .iter()
        .map(|w| {
            if is_placeholder(w) {
                w.clone()
            } else {
                gloci_io::shell_quote_args(&[OsString::from(w)])
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
}

/// Run a plan ([`run`]); a panic in this code raises `GlociError` (D90).
#[pyfunction]
fn run_plan(py: Python<'_>, plan: &Bound<'_, PyDict>) -> PyResult<Py<PyAny>> {
    catch_unwind(AssertUnwindSafe(|| run::run_plan(py, plan)))
        .unwrap_or_else(|p| Err(errors::internal(py, errors::panic_message(&*p))))
}

/// Panic in a stage (or, with `on_driver`, in the driver's own code), for the tests.
#[cfg(feature = "test-hooks")]
#[pyfunction]
#[pyo3(signature = (on_driver = false))]
fn _test_panic(py: Python<'_>, on_driver: bool) -> PyResult<()> {
    catch_unwind(AssertUnwindSafe(|| {
        if on_driver {
            panic!("test panic on the driver");
        }
        let set = gloci::embed::StageSet::spawn_fns(vec![Box::new(|| panic!("test panic"))]);
        run::drive(py, set).settle(py, |_| None).map(drop)
    }))
    .unwrap_or_else(|p| Err(errors::internal(py, errors::panic_message(&*p))))
}

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    logs::install();
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(command_specs, m)?)?;
    m.add_function(wrap_pyfunction!(bindings, m)?)?;
    m.add_function(wrap_pyfunction!(kwargs_to_argv, m)?)?;
    m.add_function(wrap_pyfunction!(step_argv, m)?)?;
    m.add_function(wrap_pyfunction!(command_line, m)?)?;
    m.add_function(wrap_pyfunction!(run_plan, m)?)?;
    #[cfg(feature = "test-hooks")]
    m.add_function(wrap_pyfunction!(_test_panic, m)?)?;
    Ok(())
}
