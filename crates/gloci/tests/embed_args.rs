//! `gloci::embed` (D83, D86, D87): parsing, in-process execution into in-memory outputs,
//! keyword arguments to argv, and the bindings checked against the clap tree.
use gloci::embed::{
    self, Bind, ExecOptions, KwValue, Outcome, Role, kwargs_to_argv, kwargs_to_argv_as,
    rewrite_flags_as_kwargs,
};
use gloci_io::cancel::{CancelToken, RunGuard};
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
    token: CancelToken,
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
    let token = CancelToken::new();
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
    let _g = RunGuard::enter(&CancelToken::new());
    assert!(run(&os(&["sort", src.name(), "-t", "4", "-o", o.name()])).unwrap());
    assert_eq!(sink.bytes(), b"chr1\t1\t2\nchr1\t5\t6\n");
}

#[test]
fn execute_leaves_the_global_rayon_pool_alone() {
    let src = Registered::input("<list>", Box::new(Cursor::new(b"chr1\t5\t6\n".to_vec())));
    let o = Registered::output("<out>", Box::new(std::io::sink()));
    assert!(run(&os(&["sort", src.name(), "-t", "3", "-o", o.name()])).unwrap());
    // Nothing in this test binary builds or uses the global pool, so a host still can.
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build_global()
        .expect("execute must not build the global rayon pool");
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
    let e = argv_of(&["closest"], &[("Q", Flag(true))]).unwrap_err();
    assert!(e.ends_with("did you mean 'q'?"), "{e}");
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
    assert!(e.0.contains("first argument"), "{e}");
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
    // Short flags and placeholders; words with dashes and unknown flags are left alone.
    assert_eq!(
        rewrite_flags_as_kwargs(
            "-f <F> and -s: half-open, --no-such-flag, chr1:1-100, -o FILE",
            i
        ),
        "f= and s=: half-open, --no-such-flag, chr1:1-100, o= FILE"
    );
    assert_eq!(
        rewrite_flags_as_kwargs(
            "the following required arguments were not provided:\n  --b <FILE>...",
            i
        ),
        "the following required arguments were not provided:\n  b="
    );
    assert_eq!(
        rewrite_flags_as_kwargs("invalid value for '[REGIONS]...'", spec(&["view"])),
        "invalid value for 'regions='"
    );
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
    // P3 commands among them.
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
    // Exactly windows and vcf extract are both a step and a module function.
    let both: Vec<String> = embed::BINDINGS
        .iter()
        .filter(|b| embed::module_function(b.path).is_some())
        .map(|b| b.path.join(" "))
        .collect();
    assert_eq!(both, ["vcf extract", "windows"]);
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
    // Methods are the path joined with `_`, and every bound input exists with its shape.
    for b in embed::BINDINGS {
        assert_eq!(b.method, b.path.join("_"));
        let s = spec(b.path);
        match b.input {
            Bind::Flag(name) => {
                let a = s
                    .arg(name)
                    .unwrap_or_else(|| panic!("{}: no argument {name}", b.method));
                assert!(a.takes_value && !a.positional, "{}: {name}", b.method);
            }
            Bind::Positional { nary } => {
                let a = s
                    .args
                    .iter()
                    .find(|a| a.positional)
                    .unwrap_or_else(|| panic!("{}: no positional", b.method));
                assert_eq!(a.multiple, nary, "{}: {}", b.method, a.id);
            }
        }
    }
    assert!(embed::is_cli_only(&["genomes", "show"]));
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
