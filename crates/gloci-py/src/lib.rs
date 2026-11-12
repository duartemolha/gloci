//! `gloci._native`: the Python library's glue over [`gloci::embed`] (D83). It holds no command
//! logic: the commands, their arguments and how a step takes its input come from the CLI's own
//! definition, and every step runs through the CLI's parsing and dispatch.
//!
//! - [`command_specs`] and [`bindings`] describe the commands and the step methods.
//! - [`kwargs_to_argv`] checks a step's keyword arguments when the step is called (D85).
//! - [`step_argv`] and [`command_line`] build a stage's argv and its shell equivalent;
//!   [`source_argv`] and [`function_argv`] do the same for the module sources (`windows`,
//!   `random`) and the module functions (`vcf_stats`, `index`, ...), which
//!   [`module_functions`] lists.
//! - [`run_plan`] runs a pipeline ([`run`]), reading Python objects as sources ([`sources`],
//!   [`frames`]), and reading its output back as records ([`records`]: iteration, `count()`)
//!   or an Arrow table ([`arrow`]).

#[cfg(windows)]
compile_error!(
    "gloci does not build on Windows: use WSL (the Windows Subsystem for Linux) and install \
     gloci there with `pip install gloci`"
);

mod arrow;
mod errors;
mod frames;
mod logs;
mod records;
mod run;
mod sources;

use gloci::embed::{
    self, ArgSpec, Bind, Binding, CommandSpec, KwValue, Role, rewrite_flags_as_kwargs,
};
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyDict, PyList, PyTuple};
use std::ffi::OsString;
use std::panic::{AssertUnwindSafe, catch_unwind};

/// Run an entry point's body: a panic in it raises `GlociError` with the internal-error text
/// (D90), not pyo3's `PanicException`.
pub(crate) fn guarded<T>(py: Python<'_>, body: impl FnOnce() -> PyResult<T>) -> PyResult<T> {
    catch_unwind(AssertUnwindSafe(|| {
        test_point();
        body()
    }))
    .unwrap_or_else(|p| Err(errors::internal(py, errors::panic_message(&*p))))
}

/// With `test-hooks`, after `_test_panic("next")`: panic in the next entry point called.
#[cfg(feature = "test-hooks")]
static PANIC_NEXT: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

fn test_point() {
    #[cfg(feature = "test-hooks")]
    if PANIC_NEXT.swap(false, std::sync::atomic::Ordering::SeqCst) {
        panic!("test panic in an entry point");
    }
}

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
    d.set_item("takes_file", a.takes_file)?;
    d.set_item("help", &a.help)?;
    Ok(d)
}

/// Every command that runs, with its arguments: `[{"path", "about", "args": [...]}]`.
#[pyfunction]
fn command_specs(py: Python<'_>) -> PyResult<Bound<'_, PyList>> {
    guarded(py, || {
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
    })
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
    guarded(py, || {
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
    })
}

/// The keyword-only parameters of each step method and module function, from
/// [`embed::step_signature`] and [`embed::function_signature`] (the stub's source too):
/// `{name: {"lead", "nary", "out", "args", "skipped"}}`: `lead` the `(name, help)` of a function's positional
/// parameters, `out` the help of its `out=` (or None), `args` the keyword parameters in order
/// (`name`, the arg's clap `id`, `alias_of` its primary name for a short flag or alias, and
/// `required`), `skipped` the arguments (`name`, `switch`) whose own name is not a Python identifier.
#[pyfunction]
fn signatures(py: Python<'_>) -> PyResult<Bound<'_, PyDict>> {
    guarded(py, || {
        let put = |out: &Bound<'_, PyDict>, name: &str, sig: embed::Signature| -> PyResult<()> {
            let d = PyDict::new(py);
            d.set_item("lead", sig.lead.to_vec())?;
            d.set_item("nary", sig.nary)?;
            d.set_item("out", sig.out)?;
            let args = PyList::empty(py);
            for k in &sig.args {
                let a = PyDict::new(py);
                a.set_item("name", &k.name)?;
                a.set_item("id", &k.arg.id)?;
                a.set_item("alias_of", &k.alias_of)?;
                a.set_item("required", k.required)?;
                args.append(a)?;
            }
            d.set_item("args", args)?;
            let skipped = PyList::empty(py);
            for k in &sig.skipped {
                let a = PyDict::new(py);
                a.set_item("name", &k.name)?;
                a.set_item("switch", k.switch)?;
                skipped.append(a)?;
            }
            d.set_item("skipped", skipped)?;
            out.set_item(name, d)
        };
        let out = PyDict::new(py);
        for b in embed::BINDINGS {
            let sig = embed::step_signature(b)
                .ok_or_else(|| errors::internal(py, format!("no spec for {}", b.path.join(" "))))?;
            put(&out, b.method, sig)?;
        }
        for (_, f) in embed::MODULE_FUNCTIONS {
            let sig = embed::function_signature(f)
                .ok_or_else(|| errors::internal(py, format!("no spec for function {f}")))?;
            put(&out, f, sig)?;
        }
        Ok(out)
    })
}

/// Whether `s` is a word only a plan makes, which a keyword argument holds where a pipeline's
/// source goes (`b=`): a token for an input (`"\0src:<k>"`, `"\0pipe:<j>"`), or a shown
/// placeholder ([`PLACEHOLDER_OPEN`]). No value a user gives holds a NUL.
fn is_plan_token(s: &str) -> bool {
    let token = s
        .strip_prefix("\0src:")
        .or_else(|| s.strip_prefix("\0pipe:"))
        .is_some_and(|n| !n.is_empty() && n.bytes().all(|c| c.is_ascii_digit()));
    let placeholder = s
        .strip_prefix(PLACEHOLDER_OPEN)
        .and_then(|r| r.strip_suffix(PLACEHOLDER_CLOSE))
        .is_some_and(|inner| !inner.contains('\0'));
    token || placeholder
}

/// What a keyword argument's text may hold.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Tokens {
    /// A value the caller gave: no NUL.
    Refused,
    /// A plan's argv: a value may be a plan token.
    Allowed,
}

impl From<bool> for Tokens {
    fn from(allowed: bool) -> Tokens {
        if allowed {
            Tokens::Allowed
        } else {
            Tokens::Refused
        }
    }
}

/// A keyword argument's value as text: a path-like object by its path, anything else by
/// `str()` (D86).
fn text(py: Python<'_>, name: &str, v: &Bound<'_, PyAny>, tokens: Tokens) -> PyResult<String> {
    let s: String = if v.hasattr("__fspath__")? {
        let os = py.import("os")?;
        os.call_method1("fsdecode", (os.call_method1("fspath", (v,))?,))?
            .extract()?
    } else {
        v.str()?.extract()?
    };
    if s.contains('\0') && !(tokens == Tokens::Allowed && is_plan_token(&s)) {
        return Err(errors::usage(
            py,
            format!("{name}= holds a NUL character, which no argument can hold"),
        ));
    }
    Ok(s)
}

/// A keyword argument's value (D86): `True` gives the flag, `False` and `None` nothing, a list
/// or tuple the flag once per item, anything else its text.
fn kw_value(py: Python<'_>, name: &str, v: &Bound<'_, PyAny>, tokens: Tokens) -> PyResult<KwValue> {
    if v.is_none() {
        return Ok(KwValue::Absent);
    }
    if let Ok(b) = v.cast::<PyBool>() {
        return Ok(KwValue::Flag(b.is_true()));
    }
    if v.is_instance_of::<PyList>() || v.is_instance_of::<PyTuple>() {
        let items = v
            .try_iter()?
            .map(|i| text(py, name, &i?, tokens))
            .collect::<PyResult<_>>()?;
        return Ok(KwValue::Values(items));
    }
    Ok(KwValue::Value(text(py, name, v, tokens)?))
}

fn kwargs(
    py: Python<'_>,
    kw: Option<&Bound<'_, PyDict>>,
    tokens: Tokens,
) -> PyResult<Vec<(String, KwValue)>> {
    let Some(kw) = kw else {
        return Ok(Vec::new());
    };
    kw.iter()
        .map(|(k, v)| {
            let name: String = k.extract()?;
            let value = kw_value(py, &name, &v, tokens)?;
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
    guarded(py, || kwargs_to_argv_in(py, &path, kwargs))
}

fn kwargs_to_argv_in(
    py: Python<'_>,
    path: &[String],
    kwargs: Option<&Bound<'_, PyDict>>,
) -> PyResult<Vec<String>> {
    let spec = spec(py, path)?;
    let kw = self::kwargs(py, kwargs, Tokens::Refused)?;
    let usage = |e: gloci_io::UsageError| errors::usage(py, e.0);
    let Some(b) = embed::binding(path) else {
        if embed::module_function(path).is_none() {
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
    guarded(py, || {
        step_argv_in(py, &path, inputs, kwargs, terminal, check)
    })
}

fn step_argv_in(
    py: Python<'_>,
    path: &[String],
    inputs: Vec<String>,
    kwargs: Option<&Bound<'_, PyDict>>,
    terminal: Vec<String>,
    check: bool,
) -> PyResult<Vec<String>> {
    let spec = spec(py, path)?;
    let b = embed::binding(path)
        .ok_or_else(|| errors::usage(py, format!("gloci {} is not a step", path.join(" "))))?;
    let kw = self::kwargs(py, kwargs, Tokens::Allowed)?;
    let os = |v: Vec<String>| v.into_iter().map(OsString::from).collect::<Vec<_>>();
    let argv =
        embed::step_argv(b, &os(inputs), &kw, &os(terminal)).map_err(|e| errors::usage(py, e.0))?;
    if check {
        check_parses(py, spec, &argv)?;
    }
    Ok(strings(argv))
}

/// The argv of a module source (`windows`, `random`; [`embed::SOURCES`]), which starts a
/// pipeline: `gloci <path> <keyword arguments> <terminal>`. Its keyword arguments take the step
/// rules ([`Role::Step`]: the output and threads belong to the terminal). With `check`, clap
/// must accept it, or a `UsageError` in keyword terms is raised. With `tokens` false (a check
/// of the values a user gave), no value may be a plan's token.
#[pyfunction]
#[pyo3(signature = (path, kwargs = None, terminal = Vec::new(), check = true, tokens = true))]
fn source_argv(
    py: Python<'_>,
    path: Vec<String>,
    kwargs: Option<&Bound<'_, PyDict>>,
    terminal: Vec<String>,
    check: bool,
    tokens: bool,
) -> PyResult<Vec<String>> {
    guarded(py, || {
        let spec = spec(py, &path)?;
        let is_source = embed::module_function(&path).is_some_and(|f| embed::SOURCES.contains(&f));
        if !is_source {
            return Err(errors::usage(
                py,
                format!("gloci {} is not a module source", path.join(" ")),
            ));
        }
        let kw = self::kwargs(py, kwargs, Tokens::from(tokens))?;
        let words =
            embed::kwargs_to_argv_as(spec, &kw, Role::Step).map_err(|e| errors::usage(py, e.0))?;
        // The terminal's options go with the flags, before any `--`.
        let sep = words.iter().position(|w| w == "--").unwrap_or(words.len());
        let argv: Vec<OsString> = std::iter::once(OsString::from("gloci"))
            .chain(path.iter().map(OsString::from))
            .chain(words[..sep].iter().cloned())
            .chain(terminal.into_iter().map(OsString::from))
            .chain(words[sep..].iter().cloned())
            .collect();
        if check {
            check_parses(py, spec, &argv)?;
        }
        Ok(strings(argv))
    })
}

/// The argv of a module function (`vcf_stats`, `vcf_extract_matrix`, `vcf_compare`, `index`;
/// [`embed::MODULE_FUNCTIONS`]), which runs at once: `gloci <path> <keyword arguments> --
/// <inputs>`, its inputs being all of its positionals ([`Role::Function`]). With `check`, clap
/// must accept it, or a `UsageError` in keyword terms is raised. With `tokens` false (a check of
/// the values a user gave), no value may be a plan's token.
#[pyfunction]
#[pyo3(signature = (path, inputs, kwargs = None, check = true, tokens = true))]
fn function_argv(
    py: Python<'_>,
    path: Vec<String>,
    inputs: Vec<String>,
    kwargs: Option<&Bound<'_, PyDict>>,
    check: bool,
    tokens: bool,
) -> PyResult<Vec<String>> {
    guarded(py, || {
        let spec = spec(py, &path)?;
        let eager = embed::module_function(&path).is_some_and(|f| !embed::SOURCES.contains(&f));
        if !eager {
            return Err(errors::usage(
                py,
                format!("gloci {} is not a module function", path.join(" ")),
            ));
        }
        let kw = self::kwargs(py, kwargs, Tokens::from(tokens))?;
        let words = embed::kwargs_to_argv_as(spec, &kw, Role::Function)
            .map_err(|e| errors::usage(py, e.0))?;
        let argv: Vec<OsString> = std::iter::once(OsString::from("gloci"))
            .chain(path.iter().map(OsString::from))
            .chain(words)
            .chain(std::iter::once(OsString::from("--")))
            .chain(inputs.into_iter().map(OsString::from))
            .collect();
        if check {
            check_parses(py, spec, &argv)?;
        }
        Ok(strings(argv))
    })
}

/// The module functions and sources (D87): `[{"function", "path", "about", "source",
/// "inputs", "threads", "outputs"}]`. `source` says it starts a pipeline (`windows`, `random`)
/// rather than run at once; `inputs` names its positionals (the sources an eager function
/// takes); `threads` and `outputs` are as in [`bindings`].
#[pyfunction]
fn module_functions(py: Python<'_>) -> PyResult<Bound<'_, PyList>> {
    guarded(py, || {
        let out = PyList::empty(py);
        for (path, function) in embed::MODULE_FUNCTIONS {
            let spec = embed::find_spec(path)
                .ok_or_else(|| errors::internal(py, format!("no spec for {}", path.join(" "))))?;
            let has = |id: &str| spec.args.iter().any(|a| a.id == id);
            let d = PyDict::new(py);
            d.set_item("function", *function)?;
            d.set_item("path", *path)?;
            d.set_item("about", &spec.about)?;
            d.set_item("source", embed::SOURCES.contains(function))?;
            let inputs: Vec<&str> = spec
                .args
                .iter()
                .filter(|a| a.positional)
                .map(|a| a.id.as_str())
                .collect();
            d.set_item("inputs", inputs)?;
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
    })
}

/// The start and end of a placeholder in a shown argv: the text between them stands for an
/// in-memory input (`<python object N>`, or `<(gloci ...)` for a nested pipeline) and is shown
/// as it is. Only a plan makes them (gloci._pipeline); no value a user gives holds a NUL, so a
/// path or value that merely looks like a placeholder is quoted like any other word (I2).
const PLACEHOLDER_OPEN: &str = "\0{";
const PLACEHOLDER_CLOSE: &str = "\0}";

/// The shell form of a step's argv (D84): each word quoted as a shell needs it, except the
/// placeholders a plan marked ([`PLACEHOLDER_OPEN`]), shown as they are without their marks,
/// and without the `--` that [`step_argv`] puts before positionals where the line reads the
/// same without it ([`embed::drop_needless_separator`]). The result holds no NUL.
#[pyfunction]
fn command_line(py: Python<'_>, argv: Vec<String>) -> PyResult<String> {
    guarded(py, || Ok(shell_form(argv)))
}

fn shell_form(argv: Vec<String>) -> String {
    let argv: Vec<OsString> = argv.into_iter().map(OsString::from).collect();
    embed::drop_needless_separator(&argv)
        .iter()
        .map(|w| match w.to_str() {
            Some(word) if word.contains(PLACEHOLDER_OPEN) => shell_word(word),
            _ => gloci_io::shell_quote_args(std::slice::from_ref(w)),
        })
        .collect::<Vec<_>>()
        .join(" ")
}

/// One word holding marked placeholders (`--x=` and a placeholder, say): the placeholders as
/// they are, the rest quoted.
fn shell_word(word: &str) -> String {
    let quote = |s: &str| gloci_io::shell_quote_args(&[OsString::from(s)]);
    let mut out = String::new();
    let mut rest = word;
    while let Some(i) = rest.find(PLACEHOLDER_OPEN) {
        if i > 0 {
            out.push_str(&quote(&rest[..i]));
        }
        let after = &rest[i + PLACEHOLDER_OPEN.len()..];
        let Some(end) = after.find(PLACEHOLDER_CLOSE) else {
            // An open mark with no close: not a placeholder, so shown as text, without the NUL.
            rest = after;
            continue;
        };
        out.push_str(&after[..end].replace('\0', ""));
        rest = &after[end + PLACEHOLDER_CLOSE.len()..];
    }
    if !rest.is_empty() {
        out.push_str(&quote(&rest.replace('\0', "")));
    }
    out
}

/// Run a plan ([`run`]); a panic in this code raises `GlociError` (D90).
#[pyfunction]
fn run_plan(py: Python<'_>, plan: &Bound<'_, PyDict>) -> PyResult<Py<PyAny>> {
    guarded(py, || run::run_plan(py, plan))
}

/// Panic, for the tests: `"stage"` in a stage; `"driver"` in the driver's own code; `"unwind"`
/// in the driver while a stage runs that ends only once cancelled and then takes the GIL, so
/// the set must be cancelled and dropped with the GIL released while unwinding (R2); `"next"`
/// in the next entry point called.
#[cfg(feature = "test-hooks")]
#[pyfunction]
#[pyo3(signature = (place = "stage"))]
fn _test_panic(py: Python<'_>, place: &str) -> PyResult<()> {
    use gloci::embed::StageSet;
    catch_unwind(AssertUnwindSafe(|| match place {
        "stage" => {
            let set = StageSet::spawn_fns(vec![Box::new(|| panic!("test panic"))]);
            let mut driver = run::Driver::new(set, Vec::new());
            driver.finish(py);
            driver.settle(py, |_| None).map(drop)
        }
        "driver" => panic!("test panic on the driver"),
        "next" => {
            PANIC_NEXT.store(true, std::sync::atomic::Ordering::SeqCst);
            Ok(())
        }
        "unwind" => {
            let set = StageSet::spawn_fns(vec![Box::new(|| {
                while !gloci_io::cancel::is_cancelled() {
                    std::thread::sleep(std::time::Duration::from_millis(5));
                }
                Python::attach(|_| ());
                Ok(true)
            })]);
            let _running = run::Driver::new(set, Vec::new());
            panic!("test panic while a step runs");
        }
        other => Err(errors::usage(py, format!("no panic place {other:?}"))),
    }))
    .unwrap_or_else(|p| Err(errors::internal(py, errors::panic_message(&*p))))
}

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    logs::install();
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(command_specs, m)?)?;
    m.add_function(wrap_pyfunction!(bindings, m)?)?;
    m.add_function(wrap_pyfunction!(signatures, m)?)?;
    m.add_function(wrap_pyfunction!(kwargs_to_argv, m)?)?;
    m.add_function(wrap_pyfunction!(step_argv, m)?)?;
    m.add_function(wrap_pyfunction!(source_argv, m)?)?;
    m.add_function(wrap_pyfunction!(function_argv, m)?)?;
    m.add_function(wrap_pyfunction!(module_functions, m)?)?;
    m.add_function(wrap_pyfunction!(command_line, m)?)?;
    m.add_function(wrap_pyfunction!(run_plan, m)?)?;
    m.add_class::<records::Record>()?;
    m.add_class::<records::RecordIterator>()?;
    m.add_class::<arrow::ArrowStream>()?;
    #[cfg(feature = "test-hooks")]
    m.add_function(wrap_pyfunction!(_test_panic, m)?)?;
    Ok(())
}
