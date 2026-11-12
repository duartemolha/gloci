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
    /// It takes a value (`-f 0.5`), rather than being a switch (`-s`).
    pub takes_value: bool,
    /// It may be given more than once, or takes several values.
    pub multiple: bool,
    pub positional: bool,
    pub required: bool,
    /// Its value must be attached (`--split=snps`).
    pub require_equals: bool,
    /// The values clap accepts, when it has a fixed list.
    pub possible_values: Vec<String>,
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

    /// Whether `name` names it: its long flag with `-` as `_`, its short letter, or (a
    /// positional) its clap id.
    pub fn named(&self, name: &str) -> bool {
        if self
            .long
            .as_deref()
            .is_some_and(|l| l.replace('-', "_") == name)
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

fn walk(cmd: &clap::Command, path: &mut Vec<String>, out: &mut Vec<CommandSpec>) {
    if cmd.get_name() == "help" {
        return;
    }
    path.push(cmd.get_name().to_string());
    let has_subs = cmd.get_subcommands().any(|s| s.get_name() != "help");
    if !has_subs || !cmd.is_subcommand_required_set() {
        out.push(CommandSpec {
            path: path.clone(),
            about: cmd.get_about().map(|s| s.to_string()).unwrap_or_default(),
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
        takes_value,
        multiple: matches!(action, clap::ArgAction::Append | clap::ArgAction::Count) || many_values,
        positional: a.is_positional(),
        required: a.is_required_set(),
        require_equals: a.is_require_equals_set(),
        possible_values: a
            .get_possible_values()
            .iter()
            .filter(|v| !v.is_hide_set())
            .map(|v| v.get_name().to_string())
            .collect(),
        help: a.get_help().map(|s| s.to_string()).unwrap_or_default(),
    }
}
