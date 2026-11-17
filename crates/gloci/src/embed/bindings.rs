//! Which commands a host offers, and how (D87, as amended). Every leaf command
//! is in exactly one of these groups, which a test checks against the clap tree, so a new
//! command cannot go unnoticed:
//! - [`BINDINGS`]: a step method of a pipeline, reading the pipeline as its input;
//! - [`MODULE_FUNCTIONS`]: a function that runs at once, or a source that starts a pipeline
//!   ([`SOURCES`]); `windows` and `vcf extract` are also steps (`gloci.windows(g=...)` is a
//!   source, `vcf_extract_matrix` writes the 012 files);
//! - [`CLI_ONLY`]: never offered;
//! - [`UNBOUND`]: commands not placed yet (a new command goes here until it is).

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

/// What a step's output is, for a host that reads it back as records (D89): every
/// terminal but `to_bytes()` and `save()` (iteration, `count()`, the Arrow tables) reads it the
/// same way, so they agree on what a record is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Shape {
    /// Records in a format gloci reads (BED, VCF, GFF, SAM, BAM, BCF, FASTA), told by the
    /// output's content; text gloci cannot read as records is a table with no header row.
    Records,
    /// A table of values, not records: its rows have no coordinates, even when they look like
    /// BED (`idxstats` writes `chr1  248956422  11  0`).
    Table,
    /// A table whose header row is its first line that is not blank and does not start with
    /// `#`: the header row names the columns and is not a row (`stats jaccard` writes
    /// `intersection  union  jaccard  n_intersections`).
    HeaderRow,
    /// A table whose header row is its first line that is not blank, whether or not it starts
    /// with `#` (`groupby --header` writes the input header's key names, with its `#` if it
    /// had one, and the `-c` labels).
    HeaderLine,
}

/// A step's [`Shape`], and the flags that change it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Output {
    pub shape: Shape,
    /// Flags (clap ids) that change the shape when given, in order: the first one given wins
    /// (`coverage --hist` writes a table, `stats reldist --detail` records).
    pub with: &'static [(&'static str, Shape)],
}

/// The output of most steps: records.
pub const RECORDS: Output = Output {
    shape: Shape::Records,
    with: &[],
};

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
    /// What its output is, read back ([`read_back`](super::read_back)).
    pub output: Output,
}

impl Binding {
    /// The binding whose output is `shape`, but as `with` says when one of its flags is given.
    const fn shaped(self, shape: Shape, with: &'static [(&'static str, Shape)]) -> Binding {
        Binding {
            output: Output { shape, with },
            ..self
        }
    }
}

const fn flag(path: &'static [&'static str], method: &'static str, f: &'static str) -> Binding {
    Binding {
        path,
        method,
        input: Bind::Flag(f),
        needs_file: false,
        output: RECORDS,
    }
}

const fn one(path: &'static [&'static str], method: &'static str) -> Binding {
    Binding {
        path,
        method,
        input: Bind::Positional { nary: false },
        needs_file: false,
        output: RECORDS,
    }
}

const fn nary(path: &'static [&'static str], method: &'static str) -> Binding {
    Binding {
        path,
        method,
        input: Bind::Positional { nary: true },
        needs_file: false,
        output: RECORDS,
    }
}

/// The step methods (D87). Each step whose output is a table rather than records is marked
/// with its [`Shape`], and the flags that change it: those that print an unprefixed header row
/// (found by running every binding's command on the golden fixtures; `tests/embed_outputs.rs`
/// keeps the list), and those whose rows would read back as wrong intervals. `groupby` writes
/// records only when its keys start with BED's chrom, start and end ([`super::read_back`]).
pub const BINDINGS: &[Binding] = &[
    // relations: the pipeline is A
    flag(&["intersect"], "intersect", "a"),
    flag(&["subtract"], "subtract", "a"),
    flag(&["closest"], "closest", "a"),
    flag(&["map"], "map", "a"),
    // --hist: A's lines with histogram columns, then `all` summary lines.
    flag(&["coverage"], "coverage", "a").shaped(Shape::Records, &[("hist", Shape::Table)]),
    flag(&["pair-intersect"], "pair_intersect", "a"),
    flag(&["pair-pair"], "pair_pair", "a"),
    flag(&["stats", "jaccard"], "stats_jaccard", "a").shaped(Shape::HeaderRow, &[]),
    // `#` lines (the contingency table), then the header row of the p-values.
    flag(&["stats", "fisher"], "stats_fisher", "a").shaped(Shape::HeaderRow, &[]),
    flag(&["stats", "reldist"], "stats_reldist", "a")
        .shaped(Shape::HeaderRow, &[("detail", Shape::Records)]),
    // several inputs: the pipeline, then extra positional sources
    nary(&["merge"], "merge"),
    nary(&["partition"], "partition"),
    nary(&["complement"], "complement"),
    nary(&["sort"], "sort"),
    nary(&["vcf", "concat"], "vcf_concat"),
    nary(&["vcf", "merge"], "vcf_merge"),
    nary(&["cat"], "cat"),
    // chrom, position, depth: a position is not a BED start, nor a depth an end.
    nary(&["depth"], "depth").shaped(Shape::Table, &[("summary", Shape::HeaderLine)]),
    nary(&["samples"], "samples").shaped(Shape::Table, &[]),
    // one input
    // --count: the number of records.
    one(&["view"], "view").shaped(Shape::Records, &[("count", Shape::Table)]),
    one(&["convert"], "convert"),
    one(&["slop"], "slop"),
    one(&["flank"], "flank"),
    one(&["shift"], "shift"),
    one(&["cluster"], "cluster"),
    one(&["spacing"], "spacing"),
    // --header: the first line names the columns. Without it, records only when keyed on
    // columns 1,2,3 (read_back).
    one(&["groupby"], "groupby").shaped(Shape::Records, &[("header", Shape::HeaderLine)]),
    // A histogram, or one depth per base, unless it writes bedGraph.
    one(&["genomecov"], "genomecov").shaped(
        Shape::Table,
        &[("bedgraph", Shape::Records), ("bga", Shape::Records)],
    ),
    one(&["vcf", "filter"], "vcf_filter"),
    one(&["vcf", "extract"], "vcf_extract").shaped(Shape::HeaderRow, &[]),
    // A format of the user's own, not records: a table, headed only under `print_header`.
    one(&["vcf", "query"], "vcf_query")
        .shaped(Shape::Table, &[("print_header", Shape::HeaderLine)]),
    one(&["vcf", "split"], "vcf_split"),
    one(&["vcf", "norm"], "vcf_norm"),
    one(&["vcf", "annotate"], "vcf_annotate"),
    // --tab: name and sequence.
    one(&["getfasta"], "getfasta").shaped(Shape::Records, &[("tab", Shape::Table)]),
    one(&["nuc"], "nuc"),
    one(&["sample"], "sample"),
    one(&["shuffle"], "shuffle"),
    Binding {
        path: &["faidx"],
        method: "faidx",
        input: Bind::Positional { nary: false },
        needs_file: true,
        output: RECORDS,
    },
    one(&["dict"], "dict"),
    one(&["flagstat"], "flagstat").shaped(Shape::Table, &[]),
    // chrom, length, mapped, unmapped.
    one(&["idxstats"], "idxstats").shaped(Shape::Table, &[]),
    one(&["head"], "head"),
    one(&["bedcov"], "bedcov"),
    one(&["tag"], "tag"),
    // An IGV batch script and an HTML page: text, not records.
    one(&["igv"], "igv").shaped(Shape::Table, &[]),
    one(&["links"], "links").shaped(Shape::Table, &[]),
    // the input as a flag or after another positional
    flag(&["maskfasta"], "maskfasta", "fasta"),
    flag(&["windows"], "windows", "i"),
    Binding {
        path: &["reheader"],
        method: "reheader",
        input: Bind::PositionalAfter("header"),
        needs_file: true,
        output: RECORDS,
    },
];

/// The module functions (D87), by command and function name: `vcf_extract_matrix` is
/// `vcf extract --genotype-matrix`; `windows` (from `-g`) and `random` are sources.
pub const MODULE_FUNCTIONS: &[(&[&str], &str)] = &[
    (&["vcf", "stats"], "vcf_stats"),
    (&["vcf", "extract"], "vcf_extract_matrix"),
    (&["vcf", "compare"], "vcf_compare"),
    (&["vcf", "popgen"], "vcf_popgen"),
    (&["depth"], "depth_files"),
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

/// Commands not placed yet. A new command goes here until it is bound, made a module function
/// or added to [`CLI_ONLY`]; every command is placed now.
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
