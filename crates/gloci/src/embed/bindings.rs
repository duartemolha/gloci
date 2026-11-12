//! Which commands a host offers, and how (D87, as amended in P2b Task 4). Every leaf command
//! is in exactly one of these groups, which a test checks against the clap tree, so a new
//! command cannot go unnoticed:
//! - [`BINDINGS`]: a step method of a pipeline, reading the pipeline as its input;
//! - [`MODULE_FUNCTIONS`]: a function that runs at once, or a source that starts a pipeline
//!   ([`SOURCES`]); `windows` and `vcf extract` are also steps (`gloci.windows(g=...)` is a
//!   source, `vcf_extract_matrix` writes the 012 files);
//! - [`CLI_ONLY`]: never offered;
//! - [`UNBOUND`]: commands not placed yet (empty; a new command goes here until it is).

use super::spec::CommandSpec;

/// How a step receives the pipeline it is called on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Bind {
    /// As this flag, named as a keyword argument (`a` for `-a`, `i` for `-i`, `fasta`).
    Flag(&'static str),
    /// As the command's first positional; `nary` when it takes several (`merge a b c`), which
    /// the method takes as extra positional sources.
    Positional { nary: bool },
    /// As the positional after the named one, which the step needs as a keyword argument
    /// (`reheader(header=...)`: `gloci reheader <header> <input>`).
    PositionalAfter(&'static str),
}

/// One step method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Binding {
    /// The command: `["vcf", "filter"]`.
    pub path: &'static [&'static str],
    /// The method: the path joined with `_`, with `-` as `_` (`pair_intersect`).
    pub method: &'static str,
    pub input: Bind,
    /// The input must be a file: the command reads it at random or twice (`faidx`,
    /// `reheader`), so it cannot follow another step or read an in-memory source.
    pub needs_file: bool,
}

const fn flag(path: &'static [&'static str], method: &'static str, f: &'static str) -> Binding {
    Binding {
        path,
        method,
        input: Bind::Flag(f),
        needs_file: false,
    }
}

const fn one(path: &'static [&'static str], method: &'static str) -> Binding {
    Binding {
        path,
        method,
        input: Bind::Positional { nary: false },
        needs_file: false,
    }
}

const fn nary(path: &'static [&'static str], method: &'static str) -> Binding {
    Binding {
        path,
        method,
        input: Bind::Positional { nary: true },
        needs_file: false,
    }
}

/// The step methods (D87).
pub const BINDINGS: &[Binding] = &[
    // relations: the pipeline is A
    flag(&["intersect"], "intersect", "a"),
    flag(&["subtract"], "subtract", "a"),
    flag(&["closest"], "closest", "a"),
    flag(&["map"], "map", "a"),
    flag(&["coverage"], "coverage", "a"),
    flag(&["pair-intersect"], "pair_intersect", "a"),
    flag(&["pair-pair"], "pair_pair", "a"),
    flag(&["stats", "jaccard"], "stats_jaccard", "a"),
    flag(&["stats", "fisher"], "stats_fisher", "a"),
    flag(&["stats", "reldist"], "stats_reldist", "a"),
    // several inputs: the pipeline, then extra positional sources
    nary(&["merge"], "merge"),
    nary(&["partition"], "partition"),
    nary(&["complement"], "complement"),
    nary(&["sort"], "sort"),
    nary(&["vcf", "concat"], "vcf_concat"),
    nary(&["vcf", "merge"], "vcf_merge"),
    nary(&["cat"], "cat"),
    nary(&["depth"], "depth"),
    nary(&["samples"], "samples"),
    // one input
    one(&["view"], "view"),
    one(&["convert"], "convert"),
    one(&["slop"], "slop"),
    one(&["flank"], "flank"),
    one(&["shift"], "shift"),
    one(&["cluster"], "cluster"),
    one(&["spacing"], "spacing"),
    one(&["groupby"], "groupby"),
    one(&["genomecov"], "genomecov"),
    one(&["vcf", "filter"], "vcf_filter"),
    one(&["vcf", "extract"], "vcf_extract"),
    one(&["vcf", "split"], "vcf_split"),
    one(&["vcf", "norm"], "vcf_norm"),
    one(&["vcf", "annotate"], "vcf_annotate"),
    one(&["getfasta"], "getfasta"),
    one(&["nuc"], "nuc"),
    one(&["sample"], "sample"),
    one(&["shuffle"], "shuffle"),
    Binding {
        path: &["faidx"],
        method: "faidx",
        input: Bind::Positional { nary: false },
        needs_file: true,
    },
    one(&["dict"], "dict"),
    one(&["flagstat"], "flagstat"),
    one(&["idxstats"], "idxstats"),
    one(&["head"], "head"),
    one(&["bedcov"], "bedcov"),
    // the input as a flag or after another positional
    flag(&["maskfasta"], "maskfasta", "fasta"),
    flag(&["windows"], "windows", "i"),
    Binding {
        path: &["reheader"],
        method: "reheader",
        input: Bind::PositionalAfter("header"),
        needs_file: true,
    },
];

/// The module functions (D87), by command and function name: `vcf_extract_matrix` is
/// `vcf extract --genotype-matrix`; `windows` (from `-g`) and `random` are sources.
pub const MODULE_FUNCTIONS: &[(&[&str], &str)] = &[
    (&["vcf", "stats"], "vcf_stats"),
    (&["vcf", "extract"], "vcf_extract_matrix"),
    (&["vcf", "compare"], "vcf_compare"),
    (&["index"], "index"),
    (&["windows"], "windows"),
    (&["random"], "random"),
];

/// The module functions that start a pipeline (`gloci.windows(g=...)`), rather than run at
/// once; their keyword arguments take the step rules ([`super::Role::Step`]).
pub const SOURCES: &[&str] = &["windows", "random"];

/// Commands that stay on the command line (D87), each with its subcommands: `flags` is a lookup
/// with no input.
pub const CLI_ONLY: &[&[&str]] = &[&["genomes"], &["check"], &["flags"]];

/// Commands not placed yet: none. A new command goes here until it is bound, made a module
/// function or added to [`CLI_ONLY`].
pub const UNBOUND: &[&[&str]] = &[];

fn same<S: AsRef<str>>(a: &[&str], b: &[S]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| *x == y.as_ref())
}

/// The step method of the command at `path`, if it has one.
pub fn binding<S: AsRef<str>>(path: &[S]) -> Option<&'static Binding> {
    BINDINGS.iter().find(|b| same(b.path, path))
}

/// The module function of the command at `path`, if it has one.
pub fn module_function<S: AsRef<str>>(path: &[S]) -> Option<&'static str> {
    MODULE_FUNCTIONS
        .iter()
        .find(|(p, _)| same(p, path))
        .map(|(_, f)| *f)
}

/// Whether the command at `path` (or a command it is under) stays on the command line.
pub fn is_cli_only<S: AsRef<str>>(path: &[S]) -> bool {
    CLI_ONLY
        .iter()
        .any(|p| p.len() <= path.len() && same(p, &path[..p.len()]))
}

/// Whether the command at `path` awaits a decision ([`UNBOUND`]).
pub fn is_unbound<S: AsRef<str>>(path: &[S]) -> bool {
    UNBOUND.iter().any(|p| same(p, path))
}

/// The name a host calls the command by: its method, else the path joined with `_`.
pub fn method_name(spec: &CommandSpec) -> String {
    binding(&spec.path)
        .map(|b| b.method.to_string())
        .unwrap_or_else(|| spec.path.join("_").replace('-', "_"))
}
