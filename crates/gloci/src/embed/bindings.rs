//! Which commands a host offers, and how (D87). Every leaf command is in exactly one of these
//! groups, which a test checks against the clap tree, so a new command cannot go unnoticed:
//! - [`BINDINGS`]: a step method of a pipeline, reading the pipeline as its input;
//! - [`MODULE_FUNCTIONS`]: a function that runs at once (`windows` and `vcf extract` are also
//!   steps: `gloci.windows(g=...)` is a source, `vcf_extract_matrix` writes the 012 files);
//! - [`CLI_ONLY`]: never offered;
//! - [`UNBOUND`]: commands D87 does not place yet (they came after it), awaiting a decision.

use super::spec::CommandSpec;

/// How a step receives the pipeline it is called on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Bind {
    /// As this flag, named as a keyword argument (`a` for `-a`, `i` for `-i`).
    Flag(&'static str),
    /// As the command's first positional; `nary` when it takes several (`merge a b c`), which
    /// the method takes as extra positional sources.
    Positional { nary: bool },
}

/// One step method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Binding {
    /// The command: `["vcf", "filter"]`.
    pub path: &'static [&'static str],
    /// The method: the path joined with `_`.
    pub method: &'static str,
    pub input: Bind,
}

const fn flag(path: &'static [&'static str], method: &'static str, f: &'static str) -> Binding {
    Binding {
        path,
        method,
        input: Bind::Flag(f),
    }
}

const fn one(path: &'static [&'static str], method: &'static str) -> Binding {
    Binding {
        path,
        method,
        input: Bind::Positional { nary: false },
    }
}

const fn nary(path: &'static [&'static str], method: &'static str) -> Binding {
    Binding {
        path,
        method,
        input: Bind::Positional { nary: true },
    }
}

/// The step methods (D87).
pub const BINDINGS: &[Binding] = &[
    flag(&["intersect"], "intersect", "a"),
    flag(&["subtract"], "subtract", "a"),
    flag(&["closest"], "closest", "a"),
    flag(&["map"], "map", "a"),
    flag(&["coverage"], "coverage", "a"),
    nary(&["merge"], "merge"),
    nary(&["partition"], "partition"),
    nary(&["complement"], "complement"),
    nary(&["sort"], "sort"),
    nary(&["vcf", "concat"], "vcf_concat"),
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
    flag(&["windows"], "windows", "i"),
];

/// The module functions (D87), by command and function name: `vcf_extract_matrix` is
/// `vcf extract --genotype-matrix`, and `windows` (from `-g`) is a source.
pub const MODULE_FUNCTIONS: &[(&[&str], &str)] = &[
    (&["vcf", "stats"], "vcf_stats"),
    (&["vcf", "extract"], "vcf_extract_matrix"),
    (&["index"], "index"),
    (&["windows"], "windows"),
];

/// Commands that stay on the command line (D87), each with its subcommands.
pub const CLI_ONLY: &[&[&str]] = &[&["genomes"], &["check"]];

/// Commands D87 does not place: P3 and the later VCF commands, which came after the plan. They
/// are neither offered nor refused until the controller decides (bind them, make them module
/// functions, or add them to [`CLI_ONLY`]).
pub const UNBOUND: &[&[&str]] = &[
    // P3a: sequences and BAM
    &["getfasta"],
    &["nuc"],
    &["maskfasta"],
    &["cat"],
    // P3b: statistics, randomness, pairs, vcf compare
    &["stats", "jaccard"],
    &["stats", "fisher"],
    &["stats", "reldist"],
    &["random"],
    &["sample"],
    &["shuffle"],
    &["pair-intersect"],
    &["pair-pair"],
    &["vcf", "compare"],
    // the VCF commands added with P3b (bcftools norm -m-, norm -f, annotate, merge)
    &["vcf", "split"],
    &["vcf", "norm"],
    &["vcf", "annotate"],
    &["vcf", "merge"],
    // P3c: samtools essentials
    &["faidx"],
    &["dict"],
    &["flags"],
    &["idxstats"],
    &["flagstat"],
    &["head"],
    &["samples"],
    &["depth"],
    &["bedcov"],
    &["reheader"],
];

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
