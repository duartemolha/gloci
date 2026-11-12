//! Keyword arguments to an argv (D86): a host's `f=0.5, s=True, c=["a", "b"]` become
//! `-f 0.5 -s -c a -c b`, which clap then parses exactly as on the command line, so the
//! defaults, validation and messages are the CLI's.

use super::bindings::{Bind, Binding, binding, method_name};
use super::spec::{ArgSpec, CommandSpec};
use super::spec::{command_specs, find_spec};
use gloci_io::UsageError;
use std::ffi::OsString;
use std::sync::OnceLock;

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

/// Arguments one command's step never takes (by command, clap id and why): what a module
/// function of the same command covers, or what writes no output for a terminal to take.
pub const RESERVED_ON_STEP_OF: &[(&[&str], &str, &str)] = &[
    (
        &["vcf", "extract"],
        "genotype_matrix",
        "genotype_matrix is set by vcf_extract_matrix(path, genotype_matrix=...): the matrix is \
written to files, not to a pipeline",
    ),
    // Every terminal passes `-o`, which `sort --check` refuses: it writes nothing.
    (
        &["sort"],
        "check",
        "sort(check=True) writes no output; check sortedness with the CLI: gloci sort --check \
FILE",
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
            if let Some((_, _, why)) = RESERVED_ON_STEP_OF.iter().find(|(path, id, _)| {
                *id == arg.id
                    && path
                        .iter()
                        .copied()
                        .eq(spec.path.iter().map(String::as_str))
            }) {
                return Some(why.to_string());
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

/// A step's argv ([`step_argv`]) as its shell form shows it: without the `--` before the
/// positionals where clap reads the line the same without it, so `gloci merge --max-gap 50 -`
/// rather than `gloci merge --max-gap 50 -- -`. The `--` stays when a positional after it
/// starts with `-` (other than `-` itself), or when the last option before it was given fewer
/// values than it can take (`-b x y`), since it would take the positionals as its values.
pub fn drop_needless_separator(argv: &[OsString]) -> Vec<OsString> {
    let mut out = argv.to_vec();
    if let Some(sep) = argv.iter().position(|a| a == "--")
        && !separator_needed(&argv[..sep], &argv[sep + 1..])
    {
        out.remove(sep);
    }
    out
}

fn separator_needed(before: &[OsString], after: &[OsString]) -> bool {
    let flag_like = |a: &OsString| a.as_encoded_bytes().first() == Some(&b'-') && a != "-";
    if after.iter().any(flag_like) {
        return true;
    }
    // The command, by the subcommand words after `gloci`.
    static ROOT: OnceLock<clap::Command> = OnceLock::new();
    let mut cmd = ROOT.get_or_init(|| {
        let mut root = crate::command();
        root.build();
        root
    });
    for w in before.iter().skip(1) {
        match w.to_str().and_then(|w| cmd.find_subcommand(w)) {
            Some(sub) => cmd = sub,
            None => break,
        }
    }
    let option = |a: &OsString| -> Option<(&clap::Arg, bool)> {
        let a = a.to_str()?;
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
            return Err(unknown(spec, name, role));
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
/// that differs only in case. A keyword reserved in `role` is never suggested: it would only
/// be refused next.
fn unknown(spec: &CommandSpec, name: &str, role: Role) -> UsageError {
    let open = |a: &&ArgSpec| reserved(spec, a, role).is_none();
    let near = spec
        .args
        .iter()
        .filter(open)
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
                .filter(open)
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
        if boundary && let Some(n) = shell_command_at(&text[i..]) {
            // A gloci command line in a hint stays a command line: it is run in a shell.
            out.push_str(&text[i..i + n]);
            i += n;
            continue;
        }
        if b[i] == b'-'
            && boundary
            && let Some((arg, kw, len)) = flag_at(&text[i..], spec)
        {
            let rest = &text[i + len..];
            let (value, used) = match placeholder_len(rest) {
                0 if arg.takes_value => value_at(rest, arg).unwrap_or((None, 0)),
                n => (None, n),
            };
            match (arg.id.as_str(), value) {
                // The terminal's options (D86): save()'s path, and threads=.
                ("output", Some(v)) => out.push_str(&format!("save(\"{v}\")")),
                ("output", None) => out.push_str("save(path)"),
                (_, v) => {
                    out.push_str(&kw);
                    out.push('=');
                    out.push_str(v.unwrap_or(""));
                }
            }
            i += len + used;
            continue;
        }
        if matches!(b[i], b'<' | b'[')
            && boundary
            && let Some((kw, len)) = positional_at(&text[i..], spec)
        {
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

/// The length of a gloci command line at the start of `s` (`gloci sort -n x.bam -o y.bam`, a
/// hint's suggestion): `gloci`, a command, and the rest up to the end of the line, a `, `, a
/// `;`, a `)` or a backtick. `None` when `s` does not start with one.
fn shell_command_at(s: &str) -> Option<usize> {
    let rest = s.strip_prefix("gloci ")?;
    let word = rest
        .split(|c: char| !(c.is_ascii_alphanumeric() || c == '-'))
        .next()?;
    if !command_specs().iter().any(|c| c.path[0] == word) {
        return None;
    }
    let end = ["\n", ", ", ";", ")", "`"]
        .iter()
        .filter_map(|e| s.find(e))
        .min()
        .unwrap_or(s.len());
    Some(end)
}

/// The value written after a flag that takes one, at the start of `s`: ` VALUE` or `=VALUE`, up
/// to a space, a quote, a `,`, `;` or `)`, without a full stop or colon that ends a sentence.
/// Only what reads as a value is one: one of the flag's possible values, or a word with a
/// digit, a `.` or a `/` in it (`4194304`, `16M`, `x.byname.bam`); any other word (`--max-mem
/// needs a size`) is prose. The flag's own value name written bare (`-o FILE`) is a
/// placeholder: `None`. Its value and the length used.
fn value_at<'s>(s: &'s str, arg: &ArgSpec) -> Option<(Option<&'s str>, usize)> {
    let body = s.strip_prefix(' ').or_else(|| s.strip_prefix('='))?;
    let n = body
        .find(|c: char| c.is_whitespace() || "'\"`,;)".contains(c))
        .unwrap_or(body.len());
    let mut v = &body[..n];
    while let Some(t) = v.strip_suffix(['.', ':']) {
        v = t;
    }
    if arg.value_names.iter().any(|p| p == v) {
        return Some((None, 1 + v.len()));
    }
    let reads_as_value = arg.possible_values.iter().any(|p| p == v)
        || v.chars()
            .any(|c| c.is_ascii_digit() || c == '.' || c == '/');
    (!v.is_empty() && reads_as_value).then(|| (Some(v), 1 + v.len()))
}

/// The flag of `spec` at the start of `s`, `--long` or `-s`: its argument, its keyword, and
/// its length.
fn flag_at<'a>(s: &str, spec: &'a CommandSpec) -> Option<(&'a ArgSpec, String, usize)> {
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
        return Some((arg, kw, 2 + n));
    }
    let c = s[1..].chars().next()?;
    if !c.is_ascii_alphabetic() || !end_ok(&s[2..]) {
        return None;
    }
    let arg = spec.args.iter().find(|a| a.short == Some(c))?;
    // A one-letter flag is its keyword (D86), but the threads are the terminal's.
    let kw = if arg.id == "threads" {
        arg.kwarg()
    } else {
        c.to_string()
    };
    Some((arg, kw, 2))
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
