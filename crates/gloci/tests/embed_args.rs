//! `gloci::embed` (D83, D86, D87): parsing, in-process execution into in-memory outputs,
//! keyword arguments to argv, and the bindings checked against the clap tree.
use gloci::embed::{
    self, Bind, ExecOptions, KwValue, Outcome, Role, kwargs_to_argv, kwargs_to_argv_as,
    rewrite_flags_as_kwargs,
};
use gloci_io::cancel::{RunGuard, RunToken};
use gloci_io::memio::{self, Registered};
use std::ffi::OsString;
use std::io::{Cursor, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

fn os(v: &[&str]) -> Vec<OsString> {
    v.iter().map(OsString::from).collect()
}

fn strs(v: &[OsString]) -> Vec<String> {
    v.iter().map(|s| s.to_string_lossy().into_owned()).collect()
}

/// A sink whose bytes the test keeps.
#[derive(Clone, Default)]
struct Shared(Arc<Mutex<Vec<u8>>>);

impl Write for Shared {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

impl Shared {
    fn bytes(&self) -> Vec<u8> {
        self.0.lock().unwrap().clone()
    }
}

/// Parse and execute `argv` (without `gloci`), as a host does.
fn run(argv: &[OsString]) -> anyhow::Result<bool> {
    let mut full = os(&["gloci"]);
    full.extend_from_slice(argv);
    let cli = embed::parse(&full).map_err(anyhow::Error::new)?;
    let o = ExecOptions {
        command_line: gloci_io::shell_quote_args(&full),
        tmp_dir: None,
    };
    embed::execute(&cli, &o)
}

#[test]
fn parse_errors_are_usage_errors() {
    let e = embed::parse(&os(&["gloci", "merge", "--bogus", "a.bed"]))
        .err()
        .expect("--bogus is refused");
    assert!(e.0.contains("'--bogus'"), "{e}");
    assert!(!e.0.contains('\x1b'), "{e:?}");
    assert!(!e.0.contains("For more information"), "{e}");
    assert!(!e.0.starts_with("error"), "{e}");
    // An unknown subcommand too.
    let e = embed::parse(&os(&["gloci", "intersekt"])).err().unwrap();
    assert!(e.0.contains("intersekt") && !e.0.contains('\x1b'), "{e}");
    assert!(embed::parse(&os(&["gloci", "merge", "a.bed"])).is_ok());
}

fn golden_dir(case: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(case)
}

/// The case's argv, with each file of the case dir replaced by `input(path)`.
fn golden_argv(case: &str, mut input: impl FnMut(&Path) -> String) -> Vec<OsString> {
    let dir = golden_dir(case);
    let args = std::fs::read_to_string(dir.join("args")).unwrap();
    args.split_whitespace()
        .map(|a| {
            let p = dir.join(a);
            if p.is_file() {
                OsString::from(input(&p))
            } else {
                OsString::from(a)
            }
        })
        .collect()
}

const EXACT_GOLDENS: [&str; 4] = [
    "intersect-pairs-vs-bedtools",
    "merge-default-vs-bedtools",
    "closest-d-vs-bedtools",
    "partition-all-vs-bedops-intersect",
];

#[test]
fn execute_equals_golden() {
    for case in EXACT_GOLDENS {
        let dir = golden_dir(case);
        assert_eq!(
            std::fs::read_to_string(dir.join("compare")).unwrap().trim(),
            "exact"
        );
        let expected = std::fs::read(dir.join("expected")).unwrap();
        // Inputs as files, then as in-memory sources named like the files.
        for in_memory in [false, true] {
            let mut held = Vec::new();
            let mut argv = golden_argv(case, |p| {
                if in_memory {
                    let name = p.file_name().unwrap().to_str().unwrap();
                    let r = Registered::input(name, Box::new(std::fs::File::open(p).unwrap()));
                    let n = r.name().to_string();
                    held.push(r);
                    n
                } else {
                    p.to_str().unwrap().to_string()
                }
            });
            let sink = Shared::default();
            let out = Registered::output(&format!("<{case} output>"), Box::new(sink.clone()));
            argv.extend(os(&["-o", out.name()]));
            let ok = run(&argv).unwrap_or_else(|e| panic!("{case}: {e:#}"));
            assert!(ok, "{case}");
            assert_eq!(
                String::from_utf8_lossy(&sink.bytes()),
                String::from_utf8_lossy(&expected),
                "{case} (inputs in memory: {in_memory})"
            );
        }
    }
}

#[test]
fn execute_reports_errors_without_printing_and_under_display_names() {
    let src = Registered::input(
        "<list>",
        Box::new(Cursor::new(b"chr1\t1\t10\nchr1\tx\t20\n".to_vec())),
    );
    let sink = Shared::default();
    let out = Registered::output("<merge output>", Box::new(sink.clone()));
    let err = run(&os(&["merge", src.name(), "-o", out.name()])).unwrap_err();
    let Outcome::Failure(msg) = embed::classify(&err) else {
        panic!("{err:#}")
    };
    assert!(msg.starts_with("<list> line 2:"), "{msg}");
    let site = embed::error_site(&err).unwrap();
    assert_eq!(site.kind, embed::SiteKind::Parse);
    assert_eq!(site.path, "<list>");
    assert_eq!(site.loc, Some(gloci_io::Loc::Line(2)));
    // A usage error found by the command, not by clap.
    let src = Registered::input("<list>", Box::new(Cursor::new(Vec::new())));
    let err = run(&os(&["index", src.name()])).unwrap_err();
    assert_eq!(
        embed::classify(&err),
        Outcome::Usage("<list> is in memory; index needs a file".into())
    );
}

/// A source that cancels its run when it reaches the end of its bytes.
struct CancelAtEof {
    inner: Cursor<Vec<u8>>,
    token: RunToken,
}

impl Read for CancelAtEof {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.inner.read(buf)?;
        if n == 0 {
            self.token.cancel();
        }
        Ok(n)
    }
}

#[test]
fn a_cancel_stops_a_command_using_the_run_pool() {
    // sort -t 4 sorts its chunks on the run's rayon pool; the token fires when the input ends,
    // so the cancellation lands in the sort and the spill, on the pool.
    let mut bed = String::new();
    for i in (0..600_000u64).rev() {
        bed.push_str(&format!("chr1\t{}\t{}\n", i * 10, i * 10 + 5));
    }
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out.bed");
    let token = RunToken::new();
    let src = Registered::input(
        "<generator>",
        Box::new(CancelAtEof {
            inner: Cursor::new(bed.into_bytes()),
            token: token.clone(),
        }),
    );
    let argv = os(&[
        "sort",
        src.name(),
        "-t",
        "4",
        "--max-mem",
        "16M",
        "--tmp-dir",
        dir.path().to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    // The run's thread enters the token; the error is judged on another thread, by the token.
    let err = std::thread::scope(|s| {
        s.spawn(|| {
            let _g = RunGuard::enter(&token);
            run(&argv).unwrap_err()
        })
        .join()
        .unwrap()
    });
    assert_eq!(
        embed::classify_under(&err, Some(&token)),
        Outcome::Cancelled
    );
    let left: Vec<_> = std::fs::read_dir(dir.path()).unwrap().collect();
    assert!(left.is_empty(), "{left:?}");
    // The next run, with a fresh token, works.
    let src = Registered::input(
        "<list>",
        Box::new(Cursor::new(b"chr1\t5\t6\nchr1\t1\t2\n".to_vec())),
    );
    let sink = Shared::default();
    let o = Registered::output("<out>", Box::new(sink.clone()));
    let _g = RunGuard::enter(&RunToken::new());
    assert!(run(&os(&["sort", src.name(), "-t", "4", "-o", o.name()])).unwrap());
    assert_eq!(sink.bytes(), b"chr1\t1\t2\nchr1\t5\t6\n");
}

fn spec(path: &[&str]) -> &'static embed::CommandSpec {
    embed::find_spec(path).unwrap_or_else(|| panic!("no spec for {path:?}"))
}

fn kw(pairs: &[(&str, KwValue)]) -> Vec<(String, KwValue)> {
    pairs
        .iter()
        .map(|(k, v)| (k.to_string(), v.clone()))
        .collect()
}

fn argv_of(path: &[&str], pairs: &[(&str, KwValue)]) -> Result<Vec<String>, String> {
    kwargs_to_argv(spec(path), &kw(pairs))
        .map(|v| strs(&v))
        .map_err(|e| e.0)
}

fn val(s: &str) -> KwValue {
    KwValue::Value(s.into())
}

#[test]
fn kwargs_rules() {
    use KwValue::{Absent, Flag, Values};
    let i = &["intersect"][..];
    assert_eq!(argv_of(i, &[("f", val("0.5"))]).unwrap(), ["-f", "0.5"]);
    assert_eq!(
        argv_of(i, &[("min_frac_a", val("50%"))]).unwrap(),
        ["--min-frac-a", "50%"]
    );
    assert_eq!(argv_of(i, &[("s", Flag(true))]).unwrap(), ["-s"]);
    assert_eq!(
        argv_of(i, &[("split", Flag(false)), ("sort", Absent)]).unwrap(),
        Vec::<String>::new()
    );
    assert_eq!(
        argv_of(&["merge"], &[("c", Values(vec!["a".into(), "b".into()]))]).unwrap(),
        ["-c", "a", "-c", "b"]
    );

    // A bad value is clap's to judge, with clap's list, in keyword terms.
    let mut full = os(&["gloci", "intersect", "-a", "x.bed", "-b", "y.bed"]);
    full.extend(kwargs_to_argv(spec(i), &kw(&[("mode", val("bogus"))])).unwrap());
    let e = embed::parse(&full).err().unwrap();
    assert!(
        e.0.contains("[possible values: any, none, pairs, clip, count]"),
        "{e}"
    );
    let msg = rewrite_flags_as_kwargs(&e.0, spec(i));
    assert!(
        msg.starts_with("invalid value 'bogus' for 'mode='"),
        "{msg}"
    );

    // Reserved: set by the terminal method.
    for (k, want) in [
        ("output", "output is set by save(path)"),
        ("o", "output is set by save(path)"),
        ("output_format", "output_format is set by save("),
        ("write_index", "write_index is set by save("),
        ("compression_level", "compression_level is set by save("),
        ("threads", "threads belongs to the terminal method"),
        ("t", "threads belongs to the terminal method"),
        ("quiet", "logger 'gloci'"),
        ("verbose", "logger 'gloci'"),
    ] {
        let e = argv_of(i, &[(k, val("x"))]).unwrap_err();
        assert!(e.contains(want), "{k}: {e}");
    }
    // The step's own input.
    let e = argv_of(i, &[("a", val("x.bed"))]).unwrap_err();
    assert!(e.starts_with("a is the input of intersect()"), "{e}");
    let e = argv_of(&["merge"], &[("files", val("x.bed"))]).unwrap_err();
    assert!(e.contains("positional arguments"), "{e}");
    let e = argv_of(&["windows"], &[("i", val("x.bed"))]).unwrap_err();
    assert!(e.starts_with("i is the input of windows()"), "{e}");

    // Unknown names, with a suggestion only when one is close.
    let e = argv_of(i, &[("fooo", val("1"))]).unwrap_err();
    assert_eq!(e, "intersect() got an unexpected keyword argument 'fooo'");
    let e = argv_of(i, &[("min_frac", val("1"))]).unwrap_err();
    assert!(
        e.starts_with(
            "intersect() got an unexpected keyword argument 'min_frac'; did you mean 'min_frac_"
        ),
        "{e}"
    );
    let e = argv_of(i, &[("same_strnad", Flag(true))]).unwrap_err();
    assert!(e.ends_with("did you mean 'same_strand'?"), "{e}");
    // A one-letter name suggests the short flag of the other case, never a far long name.
    let e = argv_of(&["merge"], &[("C", val("1"))]).unwrap_err();
    assert!(e.ends_with("did you mean 'c'?"), "{e}");
    // ... unless that flag is reserved (q is quiet, the logger's level).
    let e = argv_of(&["closest"], &[("Q", Flag(true))]).unwrap_err();
    assert!(!e.contains("did you mean"), "{e}");
    let e = argv_of(&["closest"], &[("D", Flag(true))]).unwrap_err();
    assert!(!e.contains("'a'"), "{e}");
    let e = argv_of(i, &[("F", Flag(true))]);
    assert!(e.is_ok(), "-F exists: {e:?}");
}

#[test]
fn kwargs_spelling_and_shape() {
    use KwValue::{Flag, Values};
    // A value that starts with '-' is attached, so it is never read as a flag.
    assert_eq!(
        argv_of(&["merge"], &[("max_gap", val("-5"))]).unwrap(),
        ["--max-gap=-5"]
    );
    assert_eq!(
        argv_of(&["slop"], &[("b", val("-5")), ("g", val("hg38"))]).unwrap(),
        ["-b=-5", "-g", "hg38"]
    );
    // So is the value of an argument that needs `=`; its bare flag takes the default.
    let norm = &["vcf", "norm"][..];
    assert!(spec(norm).arg("split").unwrap().require_equals);
    assert_eq!(
        argv_of(norm, &[("split", val("snps"))]).unwrap(),
        ["--split=snps"]
    );
    assert_eq!(
        argv_of(norm, &[("split", Flag(true))]).unwrap(),
        ["--split"]
    );
    for kws in [&[("split", val("snps"))][..], &[("split", Flag(true))]] {
        let mut full = os(&["gloci", "vcf", "norm", "x.vcf", "-f", "ref.fa"]);
        full.extend(kwargs_to_argv(spec(norm), &kw(kws)).unwrap());
        assert!(embed::parse(&full).is_ok(), "{full:?}");
    }
    // Positionals given as keywords go last, after `--`.
    assert_eq!(
        argv_of(
            &["view"],
            &[
                ("regions", Values(vec!["chr1:1-10".into(), "-odd".into()])),
                ("header_only", Flag(true)),
            ]
        )
        .unwrap(),
        ["--header-only", "--", "chr1:1-10", "-odd"]
    );
    let mut full = os(&["gloci", "view", "x.bed"]);
    full.extend(
        kwargs_to_argv(
            spec(&["view"]),
            &kw(&[("regions", Values(vec!["chr1:1-10".into(), "chr2".into()]))]),
        )
        .unwrap(),
    );
    assert!(embed::parse(&full).is_ok());
    // Shape errors are caught before clap.
    let e = argv_of(&["merge"], &[("s", val("1"))]).unwrap_err();
    assert_eq!(e, "s= is a switch: pass s=True, or leave it out");
    let e = argv_of(&["view"], &[("regions", Flag(true))]).unwrap_err();
    assert!(e.starts_with("regions= takes a value"), "{e}");
    let e = argv_of(
        &["intersect"],
        &[("f", val("0.5")), ("min_frac_a", val("0.5"))],
    )
    .unwrap_err();
    assert_eq!(
        e,
        "f= and min_frac_a= set the same option; pass one of them"
    );
    // A module function runs at once: it takes threads, not its input as a keyword.
    let stats = spec(&["vcf", "stats"]);
    let a = kwargs_to_argv_as(stats, &kw(&[("threads", val("2"))]), Role::Function).unwrap();
    assert_eq!(strs(&a), ["--threads", "2"]);
    let e = kwargs_to_argv_as(stats, &kw(&[("file", val("x.vcf"))]), Role::Function).unwrap_err();
    assert!(e.0.contains("an argument of the function"), "{e}");
    let e = kwargs_to_argv_as(stats, &kw(&[("quiet", Flag(true))]), Role::Function).unwrap_err();
    assert!(e.0.contains("logger 'gloci'"), "{e}");
}

#[test]
fn rewrite_flags() {
    let merge = spec(&["merge"]);
    assert_eq!(
        rewrite_flags_as_kwargs("invalid value 'x' for '--max-gap <N>': bad", merge),
        "invalid value 'x' for 'max_gap=': bad"
    );
    let i = spec(&["intersect"]);
    assert_eq!(
        rewrite_flags_as_kwargs(
            "the argument '--same-strand' cannot be used with '--opposite-strand'\n\nUsage: gloci intersect --a <FILE> --b <FILE>... --same-strand",
            i
        ),
        "the argument 'same_strand=' cannot be used with 'opposite_strand='"
    );
    // Short flags and placeholders; words with dashes and unknown flags are left alone. The
    // output is save()'s path (FILE is -o's placeholder).
    assert_eq!(
        rewrite_flags_as_kwargs(
            "-f <F> and -s: half-open, --no-such-flag, chr1:1-100, -o FILE",
            i
        ),
        "min_frac_a= and same_strand=: half-open, --no-such-flag, chr1:1-100, save(path)"
    );
    assert_eq!(
        rewrite_flags_as_kwargs(
            "the following required arguments were not provided:\n  --b <FILE>...",
            i
        ),
        "the following required arguments were not provided:\n  b="
    );
    // The message cannot tell which spelling the caller used (`c=` or `column=`), so a short flag
    // is named by its long keyword, which says what it is; a flag with no long form keeps its letter.
    let map = spec(&["map"]);
    assert_eq!(
        rewrite_flags_as_kwargs("-c 5: expected COL:OP", map),
        "column=5: expected COL:OP"
    );
    // A flag's value is written after `=`: the messages the review saw.
    let sort = spec(&["sort"]);
    assert_eq!(
        rewrite_flags_as_kwargs(
            "--max-mem 4194304 bytes is too small; use at least 16M",
            sort
        ),
        "max_mem=4194304 bytes is too small; use at least 16M"
    );
    assert_eq!(
        rewrite_flags_as_kwargs("pass --max-mem=64M, or --chrom-order lex.", sort),
        "pass max_mem=64M, or chrom_order=lex."
    );
    // A word after a flag that is no value of it stays a word.
    assert_eq!(
        rewrite_flags_as_kwargs("--max-mem needs a size, and --unique cannot be used", sort),
        "max_mem= needs a size, and unique= cannot be used"
    );
    // The output is save()'s path, and the threads are the terminal's.
    assert_eq!(
        rewrite_flags_as_kwargs(
            "write it with -o x.sorted.bed, or --output y.bed; -t 4",
            sort
        ),
        "write it with save(\"x.sorted.bed\"), or save(\"y.bed\"); threads=4"
    );
    // A gloci command line in a hint stays a command line: it is run in a shell.
    let conv = spec(&["convert"]);
    assert_eq!(
        rewrite_flags_as_kwargs(
            "x.bam is sorted by coordinate (@HD SO:coordinate), and --bedpe pairs the reads that \
are next to each other by name; sort it by name first: gloci sort -n x.bam -o x.byname.bam",
            conv
        ),
        "x.bam is sorted by coordinate (@HD SO:coordinate), and bedpe= pairs the reads that are \
next to each other by name; sort it by name first: gloci sort -n x.bam -o x.byname.bam"
    );
    assert_eq!(
        rewrite_flags_as_kwargs(
            "hint: run: gloci sort in.bcf -o s.bcf, or pass --sort to sort the inputs first",
            i
        ),
        "hint: run: gloci sort in.bcf -o s.bcf, or pass sort= to sort the inputs first"
    );
    // Positionals as clap prints them: by their value names (<IN>, <REF>), from real errors.
    for (argv, path, want) in [
        (
            &["gloci", "reheader"][..],
            &["reheader"][..],
            &["header=", "input="][..],
        ),
        (&["gloci", "dict"], &["dict"], &["fasta="]),
        (&["gloci", "view"], &["view"], &["file="]),
    ] {
        let e = embed::parse(&os(argv)).err().unwrap();
        let msg = rewrite_flags_as_kwargs(&e.0, spec(path));
        for w in want {
            assert!(msg.contains(w), "{path:?}: {w} in {msg:?} (from {:?})", e.0);
        }
        assert!(!msg.contains('<'), "{path:?}: {msg:?}");
    }
}

/// Every runnable command under `cmd`, found by walking the clap tree.
fn leaves(cmd: &clap::Command, path: &mut Vec<String>, out: &mut Vec<Vec<String>>) {
    for sub in cmd.get_subcommands().filter(|s| s.get_name() != "help") {
        path.push(sub.get_name().to_string());
        let subs = sub.get_subcommands().any(|s| s.get_name() != "help");
        if !subs || !sub.is_subcommand_required_set() {
            out.push(path.clone());
        }
        leaves(sub, path, out);
        path.pop();
    }
}

fn clap_leaves() -> Vec<Vec<String>> {
    let mut root = gloci::command();
    root.build();
    let mut out = Vec::new();
    leaves(&root, &mut Vec::new(), &mut out);
    out
}

#[test]
fn specs_cover_every_subcommand_and_argument() {
    let mut root = gloci::command();
    root.build();
    let leaves = clap_leaves();
    let specs = embed::command_specs();
    let paths: Vec<Vec<String>> = specs.iter().map(|s| s.path.clone()).collect();
    assert_eq!(paths, leaves);
    // The alignment, statistics and samtools-style commands are among them.
    for p in [
        &["depth"][..],
        &["cat"],
        &["stats", "fisher"],
        &["vcf", "merge"],
    ] {
        assert!(embed::find_spec(p).is_some(), "{p:?}");
    }
    for s in &specs {
        let mut cmd = &root;
        for name in &s.path {
            cmd = cmd.find_subcommand(name).unwrap();
        }
        let want: Vec<&str> = cmd
            .get_arguments()
            .filter(|a| !a.is_hide_set() && !matches!(a.get_id().as_str(), "help" | "version"))
            .map(|a| a.get_id().as_str())
            .collect();
        let got: Vec<&str> = s.args.iter().map(|a| a.id.as_str()).collect();
        assert_eq!(got, want, "{}", s.name());
        assert!(!s.about.is_empty(), "{}", s.name());
        // Each keyword names its own argument, and no two arguments share one.
        for a in &s.args {
            assert_eq!(
                s.arg(&a.kwarg()).map(|x| &x.id),
                Some(&a.id),
                "{} {}",
                s.name(),
                a.id
            );
            if let Some(c) = a.short {
                assert_eq!(s.arg(&c.to_string()).map(|x| &x.id), Some(&a.id));
            }
            for alias in &a.aliases {
                let k = alias.replace('-', "_");
                assert_eq!(s.arg(&k).map(|x| &x.id), Some(&a.id), "{} {k}", s.name());
            }
        }
    }
    let i = spec(&["intersect"]);
    let mode = i.arg("mode").unwrap();
    assert_eq!(
        mode.possible_values,
        ["any", "none", "pairs", "clip", "count"]
    );
    assert!(mode.takes_value && !mode.multiple && !mode.positional);
    let b = i.arg("b").unwrap();
    assert!(b.required && b.multiple && b.takes_value);
    // `-b x y` takes any number of values at once; `-a` one; a switch none.
    assert_eq!(b.max_values, None);
    assert_eq!(i.arg("a").unwrap().max_values, Some(1));
    assert_eq!(i.arg("s").unwrap().max_values, Some(0));
    // `--split[=TYPES]`: at most one value, attached.
    let split = spec(&["vcf", "norm"]).arg("split").unwrap();
    assert!(split.require_equals && split.max_values == Some(1));
    assert!(!i.arg("s").unwrap().takes_value);
    assert!(
        i.arg("f").unwrap().help.contains("fraction"),
        "{:?}",
        i.arg("f")
    );
    let files = spec(&["merge"]).arg("files").unwrap();
    assert!(files.positional && files.multiple);
}

#[test]
fn every_cli_command_is_bound_or_excluded() {
    let mut unplaced = Vec::new();
    for leaf in clap_leaves() {
        let bound = embed::binding(&leaf).is_some();
        let function = embed::module_function(&leaf).is_some();
        let cli_only = embed::is_cli_only(&leaf);
        let unbound = embed::is_unbound(&leaf);
        let groups = [bound || function, cli_only, unbound];
        match groups.iter().filter(|&&g| g).count() {
            1 => {}
            0 => unplaced.push(leaf.join(" ")),
            _ => panic!("{leaf:?} is in more than one group: {groups:?}"),
        }
    }
    assert!(
        unplaced.is_empty(),
        "place these commands (bindings.rs): {unplaced:?}"
    );
    // Exactly depth, windows and vcf extract are both a step and a module function.
    let both: Vec<String> = embed::BINDINGS
        .iter()
        .filter(|b| embed::module_function(b.path).is_some())
        .map(|b| b.path.join(" "))
        .collect();
    assert_eq!(both, ["depth", "vcf extract", "windows"]);
    // Every entry names a real command, once.
    let leaves = clap_leaves();
    let real = |p: &[&str]| {
        leaves
            .iter()
            .any(|l| l.iter().map(String::as_str).eq(p.iter().copied()))
    };
    let mut seen = std::collections::HashSet::new();
    for p in embed::BINDINGS
        .iter()
        .map(|b| b.path)
        .chain(embed::UNBOUND.iter().copied())
        .chain(embed::CLI_ONLY.iter().copied())
    {
        assert!(real(p), "{p:?} is not a command");
        assert!(seen.insert(p), "{p:?} listed twice");
    }
    for (p, _) in embed::MODULE_FUNCTIONS {
        assert!(real(p), "{p:?} is not a command");
    }
    // Methods are the path joined with `_` (and `-` as `_`), and every bound input exists with
    // its shape.
    for b in embed::BINDINGS {
        assert_eq!(b.method, b.path.join("_").replace('-', "_"));
        let s = spec(b.path);
        let input = embed::input_arg(s, b)
            .unwrap_or_else(|| panic!("{}: no argument for {:?}", b.method, b.input));
        match b.input {
            Bind::Flag(name) => {
                assert_eq!(s.arg(name), Some(input));
                assert!(
                    input.takes_value && !input.positional,
                    "{}: {name}",
                    b.method
                );
            }
            Bind::Positional { nary } => {
                assert!(input.positional);
                assert_eq!(input.multiple, nary, "{}: {}", b.method, input.id);
            }
            Bind::PositionalAfter(before) => {
                let pos: Vec<&str> = s
                    .args
                    .iter()
                    .filter(|a| a.positional)
                    .map(|a| a.id.as_str())
                    .collect();
                let i = pos.iter().position(|id| *id == before).unwrap();
                assert_eq!(pos.get(i + 1), Some(&input.id.as_str()), "{}", b.method);
                assert!(!input.multiple);
            }
        }
    }
    // Only faidx and reheader need their input as a file.
    let needs_file: Vec<&str> = embed::BINDINGS
        .iter()
        .filter(|b| b.needs_file)
        .map(|b| b.method)
        .collect();
    assert_eq!(needs_file, ["faidx", "reheader"]);
    // Everything is placed; a new command fails above until it is.
    assert!(embed::UNBOUND.is_empty(), "{:?}", embed::UNBOUND);
    for f in embed::SOURCES {
        assert!(
            embed::MODULE_FUNCTIONS.iter().any(|(_, m)| m == f),
            "source {f} is not a module function"
        );
    }
    assert!(embed::is_cli_only(&["genomes", "show"]));
    assert!(embed::is_cli_only(&["flags"]));
}

#[test]
fn step_argv_puts_inputs_where_the_binding_says() {
    use KwValue::Values;
    let b = |p: &[&str]| embed::binding(p).unwrap();
    let argv = |p: &[&str], inputs: &[&str], k: &[(&str, KwValue)]| {
        embed::step_argv(b(p), &os(inputs), &kw(k), &[]).map(|v| strs(&v))
    };
    // A flag.
    assert_eq!(
        argv(&["intersect"], &["a.bed"], &[("b", val("b.bed"))]).unwrap(),
        ["gloci", "intersect", "-a", "a.bed", "-b", "b.bed"]
    );
    assert_eq!(
        argv(&["maskfasta"], &["ref.fa"], &[("b", val("m.bed"))]).unwrap(),
        ["gloci", "maskfasta", "--fasta", "ref.fa", "-b", "m.bed"]
    );
    // A path that looks like a flag is attached; stdin is not.
    assert_eq!(
        argv(&["intersect"], &["-odd.bed"], &[("b", val("b.bed"))]).unwrap()[2],
        "-a=-odd.bed"
    );
    assert_eq!(
        argv(&["intersect"], &["-"], &[("b", val("b.bed"))]).unwrap()[2..4],
        ["-a", "-"]
    );
    // Positionals: the inputs first, after `--`, then positional keywords.
    assert_eq!(
        argv(&["merge"], &["x.bed", "y.bed"], &[("max_gap", val("5"))]).unwrap(),
        ["gloci", "merge", "--max-gap", "5", "--", "x.bed", "y.bed"]
    );
    assert_eq!(
        argv(
            &["bedcov"],
            &["t.bed"],
            &[("files", Values(vec!["a.bam".into(), "b.bam".into()]))]
        )
        .unwrap(),
        ["gloci", "bedcov", "--", "t.bed", "a.bam", "b.bam"]
    );
    // reheader: the header first, then the input.
    let r = argv(&["reheader"], &["in.bam"], &[("header", val("new.sam"))]).unwrap();
    assert_eq!(r, ["gloci", "reheader", "--", "new.sam", "in.bam"]);
    let cli = embed::parse(&os(&r.iter().map(String::as_str).collect::<Vec<_>>())).unwrap();
    let gloci::cli::Commands::Reheader(a) = &cli.command else {
        panic!("not reheader")
    };
    assert_eq!(a.header.as_str(), "new.sam");
    assert_eq!(a.input, std::path::Path::new("in.bam"));
    let e = argv(&["reheader"], &["in.bam"], &[]).unwrap_err();
    assert_eq!(
        e.0,
        "reheader() needs header=: the command reads header before its input"
    );
    let e = kwargs_to_argv(spec(&["reheader"]), &[]).unwrap_err();
    assert!(e.0.starts_with("reheader() needs header="), "{e}");
    let e = argv(&["reheader"], &["in.bam"], &[("input", val("x"))]).unwrap_err();
    assert!(e.0.starts_with("input is the input of reheader()"), "{e}");
    // The number of inputs.
    let e = argv(&["intersect"], &["a", "b"], &[]).unwrap_err();
    assert_eq!(e.0, "intersect() takes one input, not 2");
    let e = argv(&["merge"], &[], &[]).unwrap_err();
    assert_eq!(e.0, "merge() takes one or more inputs, not 0");
    // Every binding's argv parses with a placeholder input and its required keywords, or
    // misses only an argument group clap requires (slop's --both or -l/-r), never the input.
    for bd in embed::BINDINGS {
        let s = spec(bd.path);
        let input = embed::input_arg(s, bd).unwrap();
        let required: Vec<(String, KwValue)> = s
            .args
            .iter()
            .filter(|a| a.required && a.id != input.id)
            .map(|a| {
                // A value every such argument parses: a number is also a path.
                let v = a.possible_values.first().map_or("1", String::as_str);
                (a.kwarg(), val(v))
            })
            .collect();
        let v = embed::step_argv(bd, &os(&["in"]), &required, &[])
            .unwrap_or_else(|e| panic!("{}: {e}", bd.method));
        match embed::parse(&v) {
            Ok(_) => {}
            Err(e) => {
                let msg = rewrite_flags_as_kwargs(&e.0, s);
                assert!(
                    msg.starts_with("the following required arguments were not provided"),
                    "{}: {v:?}: {msg}",
                    bd.method
                );
                assert!(
                    !msg.contains(&format!("{}=", input.kwarg())),
                    "{}: {msg}",
                    bd.method
                );
            }
        }
    }
}

#[test]
fn visible_aliases_are_keywords() {
    use KwValue::Flag;
    let cmp = spec(&["vcf", "compare"]);
    let map = cmp.arg("diff_indv_map").unwrap();
    assert_eq!(map.long.as_deref(), Some("sample-map"));
    assert!(map.aliases.iter().any(|a| a == "diff-indv-map"));
    let a = kwargs_to_argv_as(
        cmp,
        &kw(&[
            ("diff_indv_map", val("m.txt")),
            ("diff_site_discordance", Flag(true)),
        ]),
        Role::Function,
    )
    .unwrap();
    assert_eq!(
        strs(&a),
        ["--diff-indv-map", "m.txt", "--diff-site-discordance"]
    );
    let e = kwargs_to_argv_as(
        cmp,
        &kw(&[("sample_map", val("m")), ("diff_indv_map", val("m"))]),
        Role::Function,
    )
    .unwrap_err();
    assert!(e.0.contains("set the same option"), "{e}");
    assert_eq!(
        rewrite_flags_as_kwargs("'--diff-indv-map <FILE>' is bad", cmp),
        "'diff_indv_map=' is bad"
    );
    // Both inputs of a module function are its arguments.
    for k in ["a", "b"] {
        let e = kwargs_to_argv_as(cmp, &kw(&[(k, val("x.vcf"))]), Role::Function).unwrap_err();
        assert!(e.0.contains("an argument of the function"), "{k}: {e}");
    }
}

#[test]
fn execute_equals_p3_goldens() {
    // One exact golden for some of the commands placed in the follow-up, run in-process into
    // an in-memory output (inputs as files: CRAM, FASTA and indexes are read at random).
    for case in [
        "depth-all-positions-vs-samtools",
        "getfasta-bed-out-strand-vs-bedtools",
        "head-bcf-vs-bcftools",
        "idxstats-cram-vs-samtools",
        "maskfasta-absent-chromosome-vs-bedtools",
        "pairtobed-11-columns-vs-bedtools",
        "stats-jaccard-both-empty-vs-bedtools",
    ] {
        let dir = golden_dir(case);
        let expected = std::fs::read(dir.join("expected")).unwrap();
        let mut argv = golden_argv(case, |p| p.to_str().unwrap().to_string());
        let sink = Shared::default();
        let out = Registered::output("<output>", Box::new(sink.clone()));
        argv.extend(os(&["-o", out.name()]));
        assert!(
            run(&argv).unwrap_or_else(|e| panic!("{case}: {e:#}")),
            "{case}"
        );
        assert_eq!(
            String::from_utf8_lossy(&sink.bytes()),
            String::from_utf8_lossy(&expected),
            "{case}"
        );
    }
}

#[test]
fn registered_names_are_discarded_on_drop() {
    let r = Registered::input("<list>", Box::new(Cursor::new(b"x".to_vec())));
    let name = r.name().to_string();
    assert_eq!(memio::display_name(&name), "<list>");
    drop(r);
    assert_eq!(memio::display_name(&name), "<memory>");
    let e = memio::take_input(&name).err().unwrap();
    assert!(
        e.to_string().contains("never registered or was discarded"),
        "{e}"
    );
}

#[test]
fn log_records_render_display_names() {
    let r = Registered::input("<generator>", Box::new(std::io::empty()));
    let text = format!("{} has no chr prefix", r.name());
    let rendered = gloci::logging::render(
        &log::Record::builder()
            .args(format_args!("{text}"))
            .level(log::Level::Warn)
            .build(),
    );
    assert_eq!(rendered, "<generator> has no chr prefix");
}

#[test]
fn faidx_refuses_an_in_memory_fasta() {
    let src = Registered::input(
        "<view output>",
        Box::new(Cursor::new(b">chr1\nACGT\n".to_vec())),
    );
    let out = Registered::output("<out>", Box::new(std::io::sink()));
    let err = run(&os(&["faidx", src.name(), "chr1:1-2", "-o", out.name()])).unwrap_err();
    assert_eq!(
        embed::classify(&err),
        Outcome::Usage(
            "<view output> is in memory; faidx needs an indexed FASTA file: save it to a file \
and index it first"
                .into()
        )
    );
}

#[test]
fn reheader_refuses_an_in_memory_input() {
    let src = Registered::input("<view output>", Box::new(Cursor::new(b"BAM\x01".to_vec())));
    let out = Registered::output("<out>", Box::new(std::io::sink()));
    let err = run(&os(&["reheader", "new.sam", src.name(), "-o", out.name()])).unwrap_err();
    assert_eq!(
        embed::classify(&err),
        Outcome::Usage(
            "<view output> is in memory; reheader reads its input twice and needs a file: save \
it to a file first"
                .into()
        )
    );
}

/// A step's argv with the terminal's options (`-o`, `-t`), as a host builds it.
fn step_with_terminal(
    path: &[&str],
    inputs: &[OsString],
    kwargs: &[(String, KwValue)],
    terminal: &[OsString],
) -> Vec<OsString> {
    embed::step_argv(embed::binding(path).unwrap(), inputs, kwargs, terminal).unwrap()
}

#[test]
fn step_argv_takes_the_terminal_options_for_every_binding_shape() {
    let file = |case: &str, name: &str| golden_dir(case).join(name).into_os_string();
    let s = |o: &OsString| o.to_str().unwrap().to_string();
    // (step, inputs, kwargs, the same command for the CLI)
    let merge_in = file("merge-default-vs-bedtools", "input.bed");
    let a = file("intersect-pairs-vs-bedtools", "input.bed");
    let b = file("intersect-pairs-vs-bedtools", "rel-b.bed");
    let bam = file("reheader-bam-vs-samtools", "input.bam");
    let hdr = file("reheader-bam-vs-samtools", "new-header.sam");
    let fa = file("faidx-line-width-7-vs-samtools", "ref.fa");
    type Case<'a> = (
        &'a [&'a str],
        Vec<OsString>,
        Vec<(String, KwValue)>,
        Vec<String>,
    );
    let cases: Vec<Case> = vec![
        // Positional
        (
            &["merge"],
            vec![merge_in.clone()],
            kw(&[("chrom_order", val("lex"))]),
            vec![
                "merge".into(),
                "--chrom-order".into(),
                "lex".into(),
                s(&merge_in),
            ],
        ),
        // Flag
        (
            &["intersect"],
            vec![a.clone()],
            kw(&[
                ("b", val(&s(&b))),
                ("mode", val("pairs")),
                ("chrom_order", val("lex")),
            ]),
            [
                "intersect",
                "--chrom-order",
                "lex",
                "-a",
                &s(&a),
                "-b",
                &s(&b),
                "--mode",
                "pairs",
            ]
            .map(String::from)
            .to_vec(),
        ),
        // PositionalAfter
        (
            &["reheader"],
            vec![bam.clone()],
            kw(&[("header", val(&s(&hdr)))]),
            vec!["reheader".into(), s(&hdr), s(&bam)],
        ),
        // Positional, with positional keywords after it
        (
            &["faidx"],
            vec![fa.clone()],
            kw(&[
                ("line_width", val("7")),
                (
                    "regions",
                    KwValue::Values(vec!["chr1:5-70".into(), "chrM:1-20".into()]),
                ),
            ]),
            [
                "faidx",
                &s(&fa),
                "--line-width",
                "7",
                "chr1:5-70",
                "chrM:1-20",
            ]
            .map(String::from)
            .to_vec(),
        ),
    ];
    for (path, inputs, kwargs, cli) in cases {
        let want = assert_cmd::Command::cargo_bin("gloci")
            .unwrap()
            .args(&cli)
            .output()
            .unwrap();
        assert!(want.status.success(), "{cli:?}");
        let sink = Shared::default();
        let out = Registered::output("<save>", Box::new(sink.clone()));
        let terminal = os(&["-o", out.name(), "-t", "2"]);
        let argv = step_with_terminal(path, &inputs, &kwargs, &terminal);
        let dash = argv.iter().position(|a| a == "--");
        let o = argv.iter().position(|a| a == "-o").unwrap();
        assert!(
            dash.is_none_or(|d| o < d),
            "{path:?}: -o after --: {argv:?}"
        );
        // The CLI records its argv[0] (the binary's path) in reheader's @PG line.
        let bin = assert_cmd::cargo::cargo_bin("gloci");
        let cli_line: Vec<OsString> = std::iter::once(bin.into_os_string())
            .chain(cli.iter().map(OsString::from))
            .collect();
        let parsed = embed::parse(&argv).unwrap_or_else(|e| panic!("{path:?}: {argv:?}: {e}"));
        let ok = embed::execute(
            &parsed,
            &ExecOptions {
                command_line: gloci_io::shell_quote_args(&cli_line),
                tmp_dir: None,
            },
        )
        .unwrap_or_else(|e| panic!("{path:?}: {e:#}"));
        assert!(ok);
        assert!(
            sink.bytes() == want.stdout,
            "{path:?}: the bytes differ from the CLI's"
        );
    }
}

#[test]
fn suggestions_never_name_a_reserved_keyword() {
    // On a step, output and threads belong to the terminal and a is the input: never suggested.
    for typo in ["outptu", "thread", "aa"] {
        let e = argv_of(&["intersect"], &[(typo, val("1"))]).unwrap_err();
        assert_eq!(
            e,
            format!("intersect() got an unexpected keyword argument '{typo}'"),
            "{typo}"
        );
    }
    // A module function takes threads, so it is suggested there.
    let e = kwargs_to_argv_as(
        spec(&["vcf", "stats"]),
        &kw(&[("thread", val("2"))]),
        Role::Function,
    )
    .unwrap_err();
    assert!(e.0.ends_with("did you mean 'threads'?"), "{e}");
}

#[test]
fn the_vcf_extract_step_leaves_the_matrix_to_its_function() {
    let ve = &["vcf", "extract"][..];
    let e = argv_of(ve, &[("genotype_matrix", val("012"))]).unwrap_err();
    assert_eq!(
        e,
        "genotype_matrix is set by vcf_extract_matrix(path, genotype_matrix=...): the matrix is \
written to files, not to a pipeline"
    );
    // The module function takes it.
    let a = kwargs_to_argv_as(
        spec(ve),
        &kw(&[("genotype_matrix", val("012"))]),
        Role::Function,
    )
    .unwrap();
    assert_eq!(strs(&a), ["--genotype-matrix", "012"]);
    assert!(argv_of(ve, &[("info", val("DP"))]).is_ok());
}

#[test]
fn the_depth_step_leaves_the_files_to_its_function() {
    let d = &["depth"][..];
    let e = argv_of(d, &[("out_prefix", val("x"))]).unwrap_err();
    assert_eq!(
        e,
        "out_prefix is set by depth_files(src, out=...): mosdepth's files are written to disk, \
not to a pipeline"
    );
    // The module function takes them.
    let a = kwargs_to_argv_as(spec(d), &kw(&[("by", val("100"))]), Role::Function).unwrap();
    assert_eq!(strs(&a), ["--by", "100"]);
    assert!(argv_of(d, &[("summary", KwValue::Flag(true))]).is_ok());
}

#[test]
fn the_sort_step_refuses_check() {
    // Every terminal passes -o, and `sort --check` writes nothing: the check stays on the CLI.
    for v in [KwValue::Flag(true), KwValue::Flag(false)] {
        let e = argv_of(&["sort"], &[("check", v)]).unwrap_err();
        assert_eq!(
            e,
            "sort(check=True) writes no output; check sortedness with the CLI: gloci sort \
--check FILE"
        );
    }
    // Never suggested for a near miss.
    let e = argv_of(&["sort"], &[("chek", KwValue::Flag(true))]).unwrap_err();
    assert_eq!(e, "sort() got an unexpected keyword argument 'chek'");
}

#[test]
fn the_vcf_merge_step_refuses_file_list() {
    // The step's inputs are sources: open(a).vcf_merge(b, c), so the list file is not a keyword.
    let e = argv_of(&["vcf", "merge"], &[("file_list", val("shards.txt"))]).unwrap_err();
    assert_eq!(e, "pass the files as sources: open(a).vcf_merge(b, c)");
    // Not suggested for a near miss.
    let e = argv_of(&["vcf", "merge"], &[("file_lst", val("x"))]).unwrap_err();
    assert_eq!(
        e,
        "vcf_merge() got an unexpected keyword argument 'file_lst'"
    );
}

#[test]
fn file_arguments_are_known() {
    let takes = |path: &[&str], name: &str| spec(path).arg(name).unwrap().takes_file;
    // Files to read, whatever clap's value type.
    for (path, name) in [
        (&["intersect"][..], "b"),
        (&["intersect"], "genome"),
        (&["intersect"], "chrom_map"),
        (&["vcf", "compare"], "a"),
        (&["vcf", "compare"], "b"),
        (&["vcf", "annotate"], "annotations"),
        (&["vcf", "norm"], "fasta"),
        (&["getfasta"], "fasta"),
        (&["bedcov"], "files"),
        (&["reheader"], "header"),
        (&["shuffle"], "exclude"),
        (&["vcf", "popgen"], "fst_pop"),
    ] {
        assert!(takes(path, name), "{path:?} {name}");
    }
    // One population per value, given once per population.
    assert!(spec(&["vcf", "popgen"]).arg("fst_pop").unwrap().multiple);
    // Text, directories, prefixes, switches and the output are not.
    for (path, name) in [
        (&["intersect"][..], "names"),
        (&["intersect"], "mode"),
        (&["intersect"], "tmp_dir"),
        (&["intersect"], "output"),
        (&["merge"], "max_gap"),
        (&["sort"], "header"),
        (&["vcf", "stats"], "out_prefix"),
        (&["vcf", "filter"], "samples"),
        (&["vcf", "popgen"], "fst_window_size"),
        (&["vcf", "popgen"], "ld_window_bp"),
        (&["vcf", "popgen"], "min_r2"),
        (&["vcf", "popgen"], "max_mem"),
    ] {
        assert!(!takes(path, name), "{path:?} {name}");
    }
    // The LD and relatedness statistics are switches; the LD window and
    // memory take a value.
    let popgen = spec(&["vcf", "popgen"]);
    for name in [
        "geno_r2",
        "hap_r2",
        "interchrom_geno_r2",
        "interchrom_hap_r2",
        "relatedness",
        "relatedness2",
    ] {
        assert!(!popgen.arg(name).unwrap().takes_value, "{name}");
    }
    for name in [
        "ld_window",
        "ld_window_bp",
        "ld_window_min",
        "ld_window_bp_min",
        "min_r2",
        "max_mem",
    ] {
        assert!(popgen.arg(name).unwrap().takes_value, "{name}");
    }
}
