//! The commands and their arguments, read from the clap tree ([`crate::command`]), so a host
//! (the Python library) offers exactly what the CLI accepts (D83, D86).

use std::sync::OnceLock;

/// One argument of a command.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArgSpec {
    /// clap's id: the field name, which names a positional's keyword argument.
    pub id: String,
    pub long: Option<String>,
    pub short: Option<char>,
    /// Other long names it is shown with (`--diff-indv-map` for `--sample-map`); each is a
    /// keyword argument too.
    pub aliases: Vec<String>,
    /// It takes a value (`-f 0.5`), rather than being a switch (`-s`).
    pub takes_value: bool,
    /// It may be given more than once, or takes several values.
    pub multiple: bool,
    /// The most values one occurrence takes (clap's `num_args`): 0 for a switch, 1 for most
    /// options, None when there is no limit (`-b x y`). A host that reads an argv back into
    /// keyword arguments (the Python parity tests) needs it to tell `-b x y` from `-b x`
    /// followed by a positional.
    pub max_values: Option<usize>,
    pub positional: bool,
    pub required: bool,
    /// Its value must be attached (`--split=snps`).
    pub require_equals: bool,
    /// The placeholders clap prints for its values (`FILE`, `REF`, `IN`); clap names a
    /// positional by them in messages (`<IN>`, `[REGION]...`).
    pub value_names: Vec<String>,
    /// The values clap accepts, when it has a fixed list.
    pub possible_values: Vec<String>,
    /// Its values are files to read ([`FILE_VALUE_NAMES`]), which a host may give an in-memory
    /// source for (D88). Any other value is text: a label, a number, a directory or a prefix.
    pub takes_file: bool,
    /// The one-line help.
    pub help: String,
}

impl ArgSpec {
    /// The keyword argument that names it (D86): the long flag with `-` as `_`, else the
    /// one-letter short flag, else (a positional) the clap id.
    pub fn kwarg(&self) -> String {
        match (&self.long, self.short) {
            (Some(l), _) => l.replace('-', "_"),
            (None, Some(c)) => c.to_string(),
            (None, None) => self.id.clone(),
        }
    }

    /// Whether `name` names it: its long flag or an alias with `-` as `_`, its short letter,
    /// or (a positional) its clap id.
    pub fn named(&self, name: &str) -> bool {
        if self
            .long
            .as_deref()
            .is_some_and(|l| l.replace('-', "_") == name)
            || self.alias_named(name).is_some()
        {
            return true;
        }
        let mut chars = name.chars();
        if let (Some(c), None) = (chars.next(), chars.next())
            && self.short == Some(c)
        {
            return true;
        }
        self.positional && self.id == name
    }

    /// The alias that `name` (with `_` for `-`) names, if any.
    pub fn alias_named(&self, name: &str) -> Option<&str> {
        self.aliases
            .iter()
            .find(|l| l.replace('-', "_") == name)
            .map(String::as_str)
    }

    /// The command-line spelling: `--long`, else `-s`.
    pub fn flag(&self) -> Option<String> {
        match (&self.long, self.short) {
            (Some(l), _) => Some(format!("--{l}")),
            (None, Some(c)) => Some(format!("-{c}")),
            (None, None) => None,
        }
    }
}

/// One command that runs (a leaf of the command tree).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CommandSpec {
    /// The subcommand names under `gloci`: `["intersect"]`, `["vcf", "filter"]`.
    pub path: Vec<String>,
    pub about: String,
    /// Every argument but `--help`, `--version` and hidden ones, in clap's order (the global
    /// `-q`/`-v` included).
    pub args: Vec<ArgSpec>,
}

impl CommandSpec {
    /// The argument a keyword argument names, by [`ArgSpec::named`]; a one-letter name is a
    /// short flag before it is a positional's id.
    pub fn arg(&self, name: &str) -> Option<&ArgSpec> {
        let one_letter = name.chars().count() == 1;
        if one_letter && let Some(a) = self.args.iter().find(|a| !a.positional && a.named(name)) {
            return Some(a);
        }
        self.args.iter().find(|a| a.named(name))
    }

    /// The command as typed: `vcf filter`.
    pub fn name(&self) -> String {
        self.path.join(" ")
    }
}

/// One spec per leaf command: a command without subcommands, or one whose subcommand is
/// optional (`genomes` runs alone and has `genomes show`).
pub fn command_specs() -> Vec<CommandSpec> {
    specs().to_vec()
}

/// The spec of the command at `path`, if it is a leaf.
pub fn find_spec<S: AsRef<str>>(path: &[S]) -> Option<&'static CommandSpec> {
    specs().iter().find(|s| {
        s.path
            .iter()
            .map(String::as_str)
            .eq(path.iter().map(AsRef::as_ref))
    })
}

fn specs() -> &'static [CommandSpec] {
    static SPECS: OnceLock<Vec<CommandSpec>> = OnceLock::new();
    SPECS.get_or_init(|| {
        let mut root = crate::command();
        root.build();
        let mut out = Vec::new();
        for sub in root.get_subcommands() {
            walk(sub, &mut Vec::new(), &mut out);
        }
        out
    })
}

/// The first paragraph of the long about (the sentence that names every tool the command
/// replaces), without its final full stop; the short about when there is no long one.
fn long_summary(cmd: &clap::Command) -> String {
    let long = cmd
        .get_long_about()
        .or_else(|| cmd.get_about())
        .map(|s| s.to_string())
        .unwrap_or_default();
    let first = long.trim().split("\n\n").next().unwrap_or("");
    let one = first.split_whitespace().collect::<Vec<_>>().join(" ");
    one.strip_suffix('.').map(String::from).unwrap_or(one)
}

fn walk(cmd: &clap::Command, path: &mut Vec<String>, out: &mut Vec<CommandSpec>) {
    if cmd.get_name() == "help" {
        return;
    }
    path.push(cmd.get_name().to_string());
    let has_subs = cmd.get_subcommands().any(|s| s.get_name() != "help");
    if !has_subs || !cmd.is_subcommand_required_set() {
        out.push(CommandSpec {
            path: path.clone(),
            about: long_summary(cmd),
            args: cmd
                .get_arguments()
                .filter(|a| !a.is_hide_set())
                .filter(|a| !matches!(a.get_id().as_str(), "help" | "version"))
                .map(arg_spec)
                .collect(),
        });
    }
    for sub in cmd.get_subcommands() {
        walk(sub, path, out);
    }
    path.pop();
}

/// The value names (clap's placeholders) of the arguments whose values are files to read, besides
/// any name holding `FILE` (`FILE`, `FILES`, `NAME|FILE`): the inputs named by what they hold
/// (`IN`, `ALN`, `BED`, `FASTA`, `REF`, `SRC`, `HEADER`) and `vcf compare`'s `A` and `B`. Clap's
/// value types cannot tell: many file arguments are text (`-b`), and `--tmp-dir` and
/// `--out-prefix` are paths that are not read.
pub const FILE_VALUE_NAMES: &[&str] = &[
    "IN", "ALN", "A", "B", "BED", "FASTA", "REF", "SRC", "HEADER",
];

fn names_a_file(name: &str) -> bool {
    name.contains("FILE") || FILE_VALUE_NAMES.contains(&name)
}

/// Whether the values of `a` are files to read ([`ArgSpec::takes_file`]).
fn takes_file(a: &clap::Arg) -> bool {
    a.get_action().takes_values()
        && a.get_id() != "output"
        && a.get_value_names()
            .unwrap_or_default()
            .iter()
            .any(|n| names_a_file(n.as_str()))
}

/// The inputs a parsed command line gives `-` (standard input) to, in the order given: each
/// argument whose values are files to read ([`ArgSpec::takes_file`]) for each `-` the user wrote
/// for it, named as the user can find it: its short flag (`-a`), else its long one, else (a
/// positional) its value name (`A`), and the value's place when the argument was given several
/// (`-b (value 2)`). `matches` are those of `cmd`, the root, as clap built it.
pub(crate) fn stdin_inputs(cmd: &clap::Command, matches: &clap::ArgMatches) -> Vec<String> {
    let (mut cmd, mut m) = (cmd, matches);
    while let Some((name, sub)) = m.subcommand() {
        let Some(c) = cmd.find_subcommand(name) else {
            break;
        };
        (cmd, m) = (c, sub);
    }
    let mut found: Vec<(usize, String)> = Vec::new();
    for a in cmd.get_arguments().filter(|a| takes_file(a)) {
        let id = a.get_id().as_str();
        if m.value_source(id) != Some(clap::parser::ValueSource::CommandLine) {
            continue;
        }
        let (Ok(Some(values)), Some(at)) = (m.try_get_raw(id), m.indices_of(id)) else {
            continue;
        };
        let values: Vec<_> = values.zip(at).collect();
        let name = match (a.get_short(), a.get_long()) {
            (Some(c), _) => format!("-{c}"),
            (None, Some(l)) => format!("--{l}"),
            (None, None) => a
                .get_value_names()
                .and_then(|n| n.first())
                .map_or_else(|| id.to_uppercase(), |n| n.to_string()),
        };
        for (k, (v, i)) in values.iter().enumerate() {
            if *v == "-" {
                let label = if values.len() > 1 {
                    format!("{name} (value {})", k + 1)
                } else {
                    name.clone()
                };
                found.push((*i, label));
            }
        }
    }
    found.sort_by_key(|(i, _)| *i);
    found.into_iter().map(|(_, label)| label).collect()
}

fn arg_spec(a: &clap::Arg) -> ArgSpec {
    let action = a.get_action();
    let takes_value = action.takes_values();
    let many_values = a
        .get_num_args()
        .is_some_and(|r| r.max_values() > 1 && takes_value);
    ArgSpec {
        id: a.get_id().to_string(),
        long: a.get_long().map(String::from),
        short: a.get_short(),
        aliases: a
            .get_visible_aliases()
            .unwrap_or_default()
            .into_iter()
            .map(String::from)
            .collect(),
        takes_value,
        multiple: matches!(action, clap::ArgAction::Append | clap::ArgAction::Count) || many_values,
        max_values: if takes_value {
            match a.get_num_args().map(|r| r.max_values()) {
                Some(usize::MAX) => None,
                Some(n) => Some(n),
                None => Some(1),
            }
        } else {
            Some(0)
        },
        positional: a.is_positional(),
        required: a.is_required_set(),
        require_equals: a.is_require_equals_set(),
        value_names: a
            .get_value_names()
            .unwrap_or_default()
            .iter()
            .map(|n| n.to_string())
            .collect(),
        possible_values: a
            .get_possible_values()
            .iter()
            .filter(|v| !v.is_hide_set())
            .map(|v| v.get_name().to_string())
            .collect(),
        takes_file: takes_file(a),
        help: a.get_help().map(|s| s.to_string()).unwrap_or_default(),
    }
}
