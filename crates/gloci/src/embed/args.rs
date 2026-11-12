//! Keyword arguments to an argv (D86): a host's `f=0.5, s=True, c=["a", "b"]` become
//! `-f 0.5 -s -c a -c b`, which clap then parses exactly as on the command line, so the
//! defaults, validation and messages are the CLI's.

use super::bindings::{Bind, Binding, binding, method_name};
use super::spec::find_spec;
use super::spec::{ArgSpec, CommandSpec};
use gloci_io::UsageError;
use std::ffi::OsString;

/// The value of one keyword argument.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KwValue {
    /// `True`: the switch is given; `False`: it is not.
    Flag(bool),
    /// One value (any other host value, as text).
    Value(String),
    /// A list: the flag once per value.
    Values(Vec<String>),
    /// `None`: as if not given.
    Absent,
}

/// What the keyword arguments are for, which decides the ones a host sets elsewhere.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Role {
    /// A step of a pipeline: the output, its format, index and compression, and the threads
    /// belong to the terminal method; the input is the pipeline the step is called on.
    Step,
    /// A module function (`vcf_stats`, `index`, ...), which runs at once: only its inputs, the
    /// command's positionals (both files of `vcf compare`), are set elsewhere. The sources
    /// (`windows`, `random`) start a pipeline, so they take the step rules.
    Function,
}

/// The arguments a step never takes as keyword arguments (by clap id), and why (D86).
pub const RESERVED_ON_STEPS: &[(&str, &str)] = &[
    ("output", "output is set by save(path)"),
    (
        "output_format",
        "output_format is set by save(path, output_format=...)",
    ),
    (
        "write_index",
        "write_index is set by save(path, write_index=True)",
    ),
    (
        "compression_level",
        "compression_level is set by save(path, compression_level=...)",
    ),
    (
        "threads",
        "threads belongs to the terminal method: pass threads= to save(), to_bytes(), count() \
or to_arrow()",
    ),
];

/// The arguments no keyword argument sets, in any role: the log level is the host's logger's.
pub const RESERVED_ALWAYS: &[(&str, &str)] = &[
    (
        "quiet",
        "quiet is the level of the logger 'gloci': logging.getLogger('gloci').setLevel(logging.ERROR)",
    ),
    (
        "verbose",
        "verbose is the level of the logger 'gloci': logging.getLogger('gloci').setLevel(logging.INFO)",
    ),
];

/// Why `arg` of the command `spec` is not a keyword argument in `role`, if it is not.
pub fn reserved(spec: &CommandSpec, arg: &ArgSpec, role: Role) -> Option<String> {
    reserved_as(spec, arg, &arg.kwarg(), role)
}

/// [`reserved`], saying `kw`, the name the host used.
fn reserved_as(spec: &CommandSpec, arg: &ArgSpec, kw: &str, role: Role) -> Option<String> {
    let listed = |list: &[(&str, &str)]| {
        list.iter()
            .find(|(id, _)| *id == arg.id)
            .map(|(_, why)| why.to_string())
    };
    if let Some(why) = listed(RESERVED_ALWAYS) {
        return Some(why);
    }
    match role {
        Role::Step => {
            if let Some(why) = listed(RESERVED_ON_STEPS) {
                return Some(why);
            }
            let b = binding(&spec.path)?;
            (input_arg(spec, b) == Some(arg)).then(|| match b.input {
                Bind::Positional { nary: true } => format!(
                    "{kw} is the input of {m}(): {m}() reads the pipeline it is called on, plus \
any sources passed as positional arguments",
                    m = b.method
                ),
                _ => format!(
                    "{kw} is the input of {m}(): {m}() reads the pipeline it is called on",
                    m = b.method
                ),
            })
        }
        Role::Function => arg
            .positional
            .then(|| format!("{kw} is an argument of the function, not a keyword")),
    }
}

/// The argument a step's input fills.
pub fn input_arg<'s>(spec: &'s CommandSpec, b: &Binding) -> Option<&'s ArgSpec> {
    let mut positionals = spec.args.iter().filter(|a| a.positional);
    match b.input {
        Bind::Flag(name) => spec.args.iter().find(|a| !a.positional && a.named(name)),
        Bind::Positional { .. } => positionals.next(),
        Bind::PositionalAfter(before) => {
            positionals.by_ref().find(|a| a.named(before))?;
            positionals.next()
        }
    }
}

/// The argv for keyword arguments of a pipeline step ([`Role::Step`]); see
/// [`kwargs_to_argv_as`]. This is the keywords' part only, for checking them: a step's whole
/// argv comes from [`step_argv`], which puts the inputs where the binding says.
pub fn kwargs_to_argv(
    spec: &CommandSpec,
    kwargs: &[(String, KwValue)],
) -> Result<Vec<OsString>, UsageError> {
    kwargs_to_argv_as(spec, kwargs, Role::Step)
}

/// The argv for keyword arguments (D86), in their order: the flags, then `--` and the
/// positionals given as keyword arguments (named by their clap id).
///
/// Where the inputs go depends on the command, so this is not a command line to prepend
/// `gloci <path> <input>` to: a step's input may be a flag (`-a`), or follow a positional
/// keyword (`reheader`'s input comes after `header=`), and nothing may follow the `--`. Build a
/// step's argv with [`step_argv`]. A module function's inputs are all of its positionals, so
/// `gloci <path> <inputs> <these>` is right for [`Role::Function`].
///
/// A kwarg is a long flag with `-` as `_`, or a one-letter short flag; the spelling follows:
/// `s=True` gives `-s`, `min_frac_a="50%"` gives `--min-frac-a 50%`. `True` gives the bare
/// flag, `False` and `Absent` nothing, a list the flag once per value. A value that starts with
/// `-` is attached (`--shift=-5`), so it is never read as a flag, and so is the value of an
/// argument that needs it (`--split=snps`). Values are not checked here:
/// clap checks them when the argv is parsed, and [`super::rewrite_flags_as_kwargs`] turns its
/// message into keyword terms.
pub fn kwargs_to_argv_as(
    spec: &CommandSpec,
    kwargs: &[(String, KwValue)],
    role: Role,
) -> Result<Vec<OsString>, UsageError> {
    let (mut argv, positionals) = split_kwargs(spec, kwargs, role)?;
    if !positionals.is_empty() {
        argv.push("--".into());
        argv.extend(positionals);
    }
    Ok(argv)
}

/// The whole argv of a step (`gloci` first), with its inputs where its binding puts them:
/// `inputs` are the paths or in-memory names of the pipeline the step is called on and, for an
/// n-ary command, of the extra sources, in order; `terminal` are the options the host adds for
/// the terminal method (`-o`, `-O`, `-t`, ...), which go with the flags, before any `--`:
/// - [`Bind::Flag`]: `gloci intersect -a <input> <flags> <terminal> [-- <positionals>]`;
/// - [`Bind::Positional`]: `gloci merge <flags> <terminal> -- <inputs> <positionals>`;
/// - [`Bind::PositionalAfter`]: `gloci reheader <flags> <terminal> -- <header> <input>`.
///
/// The keyword arguments follow [`kwargs_to_argv`]; a [`Bind::PositionalAfter`] step needs the
/// positional its input follows. Nothing may be appended to the result: after `--` it would
/// be read as an input.
pub fn step_argv(
    b: &Binding,
    inputs: &[OsString],
    kwargs: &[(String, KwValue)],
    terminal: &[OsString],
) -> Result<Vec<OsString>, UsageError> {
    let spec = find_spec(b.path)
        .ok_or_else(|| UsageError(format!("{}: no such command", b.path.join(" "))))?;
    let one = !matches!(b.input, Bind::Positional { nary: true });
    if inputs.is_empty() || (one && inputs.len() > 1) {
        return Err(UsageError(format!(
            "{}() takes {} input{}, not {}",
            b.method,
            if one { "one" } else { "one or more" },
            if one { "" } else { "s" },
            inputs.len()
        )));
    }
    let (flags, positionals) = split_kwargs(spec, kwargs, Role::Step)?;
    let mut argv: Vec<OsString> = std::iter::once("gloci")
        .chain(b.path.iter().copied())
        .map(OsString::from)
        .collect();
    match b.input {
        Bind::Flag(name) => {
            let f = if name.chars().count() == 1 {
                format!("-{name}")
            } else {
                input_arg(spec, b)
                    .and_then(ArgSpec::flag)
                    .unwrap_or_default()
            };
            let input = &inputs[0];
            if input
                .to_str()
                .is_some_and(|i| i.starts_with('-') && i != "-")
            {
                let mut joined = OsString::from(format!("{f}="));
                joined.push(input);
                argv.push(joined);
            } else {
                argv.push(f.into());
                argv.push(input.clone());
            }
            argv.extend(flags);
            argv.extend(terminal.iter().cloned());
            if !positionals.is_empty() {
                argv.push("--".into());
                argv.extend(positionals);
            }
        }
        Bind::Positional { .. } => {
            argv.extend(flags);
            argv.extend(terminal.iter().cloned());
            argv.push("--".into());
            argv.extend(inputs.iter().cloned());
            argv.extend(positionals);
        }
        Bind::PositionalAfter(_) => {
            argv.extend(flags);
            argv.extend(terminal.iter().cloned());
            argv.push("--".into());
            argv.extend(positionals);
            argv.extend(inputs.iter().cloned());
        }
    }
    Ok(argv)
}

/// The flags and the positionals of keyword arguments, checked.
fn split_kwargs(
    spec: &CommandSpec,
    kwargs: &[(String, KwValue)],
    role: Role,
) -> Result<(Vec<OsString>, Vec<OsString>), UsageError> {
    let mut argv: Vec<OsString> = Vec::new();
    let mut positionals: Vec<OsString> = Vec::new();
    let mut seen: Vec<(&str, &str)> = Vec::new();
    for (name, value) in kwargs {
        let Some(arg) = spec.arg(name) else {
            return Err(unknown(spec, name));
        };
        if let Some(&(_, first)) = seen.iter().find(|(id, _)| *id == arg.id) {
            return Err(UsageError(format!(
                "{first}= and {name}= set the same option; pass one of them"
            )));
        }
        seen.push((&arg.id, name));
        if let Some(why) = reserved_as(spec, arg, name, role) {
            return Err(UsageError(why));
        }
        let flag = if arg.positional {
            None
        } else if name.chars().count() == 1 && arg.short.is_some() {
            Some(format!("-{name}"))
        } else if let Some(alias) = arg.alias_named(name) {
            Some(format!("--{alias}"))
        } else {
            arg.flag()
        };
        let values: Vec<&String> = match value {
            KwValue::Absent | KwValue::Flag(false) => continue,
            KwValue::Flag(true) => match &flag {
                Some(f) => {
                    argv.push(f.into());
                    continue;
                }
                None => {
                    return Err(UsageError(format!(
                        "{name}= takes a value, not True; pass {name}=\"...\""
                    )));
                }
            },
            KwValue::Value(v) => vec![v],
            KwValue::Values(vs) => vs.iter().collect(),
        };
        if !arg.takes_value {
            return Err(UsageError(format!(
                "{name}= is a switch: pass {name}=True, or leave it out"
            )));
        }
        match &flag {
            None => positionals.extend(values.into_iter().map(OsString::from)),
            Some(f) => {
                for v in values {
                    if v.starts_with('-') || arg.require_equals {
                        argv.push(format!("{f}={v}").into());
                    } else {
                        argv.push(f.into());
                        argv.push(v.into());
                    }
                }
            }
        }
    }
    if role == Role::Step
        && let Some(b) = binding(&spec.path)
        && let Bind::PositionalAfter(before) = b.input
        && !seen
            .iter()
            .any(|(_, name)| spec.arg(name).is_some_and(|a| a.named(before)))
    {
        return Err(UsageError(format!(
            "{m}() needs {before}=: the command reads {before} before its input",
            m = b.method
        )));
    }
    Ok((argv, positionals))
}

/// `unexpected keyword argument 'fooo'`, with the closest keyword argument when one is within
/// two edits (and closer than its own length, so `D` does not suggest `a`), or a short flag
/// that differs only in case.
fn unknown(spec: &CommandSpec, name: &str) -> UsageError {
    let near = spec
        .args
        .iter()
        .filter(|a| a.long.is_some() || a.positional)
        .flat_map(|a| {
            std::iter::once(a.kwarg()).chain(a.aliases.iter().map(|l| l.replace('-', "_")))
        })
        .map(|k| (gloci_io::edit_distance(name, &k), k))
        .filter(|(d, k)| *d <= 2 && *d < k.chars().count())
        .min_by_key(|(d, _)| *d)
        .map(|(_, k)| k)
        .or_else(|| {
            let mut chars = name.chars();
            let (Some(c), None) = (chars.next(), chars.next()) else {
                return None;
            };
            let other = if c.is_ascii_lowercase() {
                c.to_ascii_uppercase()
            } else {
                c.to_ascii_lowercase()
            };
            spec.args
                .iter()
                .any(|a| a.short == Some(other))
                .then(|| other.to_string())
        });
    let mut msg = format!(
        "{}() got an unexpected keyword argument '{name}'",
        method_name(spec)
    );
    if let Some(k) = near {
        msg.push_str(&format!("; did you mean '{k}'?"));
    }
    UsageError(msg)
}

/// `msg` (a clap or gloci message about the command `spec`) with every flag of the command in
/// keyword terms (D86): `'--max-gap <N>'` becomes `'max_gap='`, `-s` becomes `s=`. clap's
/// `Usage:` line and its `-- --flag` tip, which describe a command line, are dropped.
pub fn rewrite_flags_as_kwargs(msg: &str, spec: &CommandSpec) -> String {
    let kept: Vec<&str> = msg
        .lines()
        .filter(|l| !l.starts_with("Usage: "))
        .filter(|l| !l.trim_start().starts_with("tip: to pass '"))
        .collect();
    let text = kept.join("\n");
    let mut text = text.trim_end().to_string();
    // A blank line left at the end, or doubled, by the dropped lines.
    while text.contains("\n\n\n") {
        text = text.replace("\n\n\n", "\n\n");
    }
    rewrite(&text, spec)
}

fn rewrite(text: &str, spec: &CommandSpec) -> String {
    let b = text.as_bytes();
    let mut out = String::with_capacity(text.len());
    let mut i = 0;
    while i < b.len() {
        let boundary = i == 0 || !(b[i - 1].is_ascii_alphanumeric() || b[i - 1] == b'-');
        let found = match b[i] {
            b'-' if boundary => flag_at(&text[i..], spec)
                .map(|(kw, len)| (kw, len + placeholder_len(&text[i + len..]))),
            b'<' | b'[' if boundary => positional_at(&text[i..], spec),
            _ => None,
        };
        if let Some((kw, len)) = found {
            out.push_str(&kw);
            out.push('=');
            i += len;
            continue;
        }
        let ch = text[i..].chars().next().unwrap_or('\u{fffd}');
        out.push(ch);
        i += ch.len_utf8();
    }
    out
}

/// The kwarg for a `--long` or `-s` flag of `spec` at the start of `s`, with its length.
fn flag_at(s: &str, spec: &CommandSpec) -> Option<(String, usize)> {
    let end_ok = |rest: &str| {
        !rest
            .chars()
            .next()
            .is_some_and(|c| c.is_ascii_alphanumeric() || c == '-' || c == '_')
    };
    if let Some(long) = s.strip_prefix("--") {
        let n = long
            .find(|c: char| !(c.is_ascii_alphanumeric() || c == '-'))
            .unwrap_or(long.len());
        let name = &long[..n];
        let arg = spec
            .args
            .iter()
            .find(|a| a.long.as_deref() == Some(name) || a.aliases.iter().any(|l| l == name))?;
        let kw = if arg.long.as_deref() == Some(name) {
            arg.kwarg()
        } else {
            name.replace('-', "_")
        };
        return Some((kw, 2 + n));
    }
    let c = s[1..].chars().next()?;
    if !c.is_ascii_alphabetic() || !end_ok(&s[2..]) {
        return None;
    }
    spec.args.iter().find(|a| a.short == Some(c))?;
    Some((c.to_string(), 2))
}

/// The length of a value placeholder after a flag: ` <N>`, `=<N>`, ` <N>...`; 0 when none.
fn placeholder_len(s: &str) -> usize {
    let body = s
        .strip_prefix(' ')
        .or_else(|| s.strip_prefix('='))
        .unwrap_or("");
    if !body.starts_with('<') {
        return 0;
    }
    let Some(close) = body.find('>') else {
        return 0;
    };
    let inner = &body[1..close];
    if inner.is_empty()
        || !inner
            .chars()
            .all(|c| c.is_ascii_uppercase() || c.is_ascii_digit() || "_-|,.+".contains(c))
    {
        return 0;
    }
    let mut n = 1 + close + 1;
    if body[close + 1..].starts_with("...") {
        n += 3;
    }
    n
}

/// The kwarg for a positional of `spec` written as clap does (`<IN>`, `[REGION]...`): by its
/// value name, or its id in capitals when it has none; with its length.
fn positional_at(s: &str, spec: &CommandSpec) -> Option<(String, usize)> {
    let close = if s.starts_with('<') { '>' } else { ']' };
    let end = s.find(close)?;
    let name = &s[1..end];
    let arg = spec.args.iter().find(|a| {
        a.positional
            && if a.value_names.is_empty() {
                a.id.to_uppercase() == name
            } else {
                a.value_names.iter().any(|v| v == name)
            }
    })?;
    let mut n = end + 1;
    if s[n..].starts_with("...") {
        n += 3;
    }
    Some((arg.kwarg(), n))
}
