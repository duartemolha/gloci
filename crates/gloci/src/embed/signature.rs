//! The keyword arguments of a step method or module function, as the Python stub declares them
//! and the runtime's `__signature__` lists them (P2b Task 9): both come from here, so they agree.

use super::args::{Role, reserved};
use super::bindings::{Bind, Binding};
use super::spec::{ArgSpec, CommandSpec, find_spec};

/// The parameters of one method or function.
#[derive(Debug, Clone)]
pub struct Signature {
    /// The command's spec.
    pub spec: &'static CommandSpec,
    /// The leading positional parameters of a module function (`src`, or `a` and `b`), each
    /// with its help.
    pub lead: &'static [(&'static str, &'static str)],
    /// A step that takes extra sources as positional arguments (`merge(*sources)`).
    pub nary: bool,
    /// A function that writes files named by a prefix takes it as the required keyword `out`
    /// (the command's `--out-prefix`, which is not a keyword of its own); its help.
    pub out: Option<&'static str>,
    /// The keyword-only parameters, in the command's order: each argument by its [`ArgSpec::kwarg`]
    /// name, followed by its one-letter short flag and visible aliases (D86, D87).
    pub args: Vec<Kw>,
    /// Keyword names that are not Python identifiers (`genomecov`'s `5p`): they are accepted
    /// at runtime (`**{"5p": True}`) but have no parameter.
    pub skipped: Vec<String>,
}

/// One keyword-only parameter.
#[derive(Debug, Clone)]
pub struct Kw {
    pub name: String,
    pub arg: &'static ArgSpec,
    /// The primary name, when this is a short flag or an alias of it.
    pub alias_of: Option<String>,
    /// It has no default. An argument that is required but has other names is not: the
    /// signature cannot say "one of these", so the command line's check reports a missing one.
    pub required: bool,
}

const SOURCE_HELP: &str = "Any source: a path, a file object, an iterable of lines or tuples, \
an Arrow table, or a Pipeline.";

/// The help of the keyword `out`.
pub const OUT_HELP: &str =
    "The prefix of the files written (`--out-prefix`); the files are returned, sorted.";

/// A module function's leading positional parameters, with their help.
pub fn function_inputs(function: &str) -> &'static [(&'static str, &'static str)] {
    match function {
        "vcf_stats" | "vcf_extract_matrix" => &[("src", SOURCE_HELP)],
        "vcf_compare" => &[("a", SOURCE_HELP), ("b", SOURCE_HELP)],
        "index" => &[(
            "path",
            "The file to index: a path (an index is of a file, so no other source will do).",
        )],
        _ => &[],
    }
}

/// Whether `name` can be a Python parameter. `genomecov`'s `--5p` and `--3p` cannot: they stay
/// keyword arguments of the runtime (`**{"5p": True}`), but have no place in a signature.
fn is_identifier(name: &str) -> bool {
    const KEYWORDS: &[&str] = &[
        "False", "None", "True", "and", "as", "assert", "async", "await", "break", "class",
        "continue", "def", "del", "elif", "else", "except", "finally", "for", "from", "global",
        "if", "import", "in", "is", "lambda", "nonlocal", "not", "or", "pass", "raise", "return",
        "try", "while", "with", "yield",
    ];
    let mut chars = name.chars();
    chars
        .next()
        .is_some_and(|c| c.is_ascii_alphabetic() || c == '_')
        && chars.all(|c| c.is_ascii_alphanumeric() || c == '_')
        && !KEYWORDS.contains(&name)
}

fn keywords(spec: &'static CommandSpec, role: Role, taken: &[&str]) -> (Vec<Kw>, Vec<String>) {
    let args: Vec<&'static ArgSpec> = spec
        .args
        .iter()
        .filter(|a| reserved(spec, a, role).is_none())
        .filter(|a| !(role == Role::Function && a.id == "out_prefix"))
        .collect();
    let mut skipped = Vec::new();
    let mut used: Vec<String> = taken.iter().map(|s| s.to_string()).collect();
    let primary = |a: &ArgSpec| a.kwarg();
    for a in &args {
        used.push(primary(a));
    }
    let mut out = Vec::new();
    for a in args {
        let name = a.kwarg();
        let mut extras: Vec<String> = Vec::new();
        if a.long.is_some() {
            if let Some(c) = a.short {
                extras.push(c.to_string());
            }
        }
        extras.extend(a.aliases.iter().map(|l| l.replace('-', "_")));
        let extras: Vec<String> = extras
            .into_iter()
            .filter(|e| {
                if !is_identifier(e) {
                    skipped.push(e.clone());
                    false
                } else if used.contains(e) {
                    false
                } else {
                    used.push(e.clone());
                    true
                }
            })
            .collect();
        if !is_identifier(&name) {
            skipped.push(name);
            continue;
        }
        out.push(Kw {
            name: name.clone(),
            arg: a,
            alias_of: None,
            required: a.required && extras.is_empty(),
        });
        for e in extras {
            out.push(Kw {
                name: e,
                arg: a,
                alias_of: Some(name.clone()),
                required: false,
            });
        }
    }
    (out, skipped)
}

/// The signature of a step method.
pub fn step_signature(b: &Binding) -> Option<Signature> {
    let spec = find_spec(b.path)?;
    let (args, skipped) = keywords(spec, Role::Step, &[]);
    Some(Signature {
        spec,
        lead: &[],
        nary: b.input == Bind::Positional { nary: true },
        out: None,
        args,
        skipped,
    })
}

/// The signature of a module function (`vcf_stats`) or source (`windows`, `random`), by name.
pub fn function_signature(function: &str) -> Option<Signature> {
    let (path, _) = super::bindings::MODULE_FUNCTIONS
        .iter()
        .find(|(_, f)| *f == function)?;
    let spec = find_spec(path)?;
    let source = super::bindings::SOURCES.contains(&function);
    let role = if source { Role::Step } else { Role::Function };
    let lead = function_inputs(function);
    let mut taken: Vec<&str> = lead.iter().map(|(n, _)| *n).collect();
    taken.push("out");
    let (args, skipped) = keywords(spec, role, &taken);
    Some(Signature {
        spec,
        lead,
        nary: false,
        out: (!source && spec.args.iter().any(|a| a.id == "out_prefix")).then_some(OUT_HELP),
        args,
        skipped,
    })
}
