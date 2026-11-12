//! What each step's output reads back as (D89; P2b Task 8, ruling R1): iteration, `count()` and
//! the Arrow tables of the Python library share one notion of a record, so a step whose output
//! starts with an unprefixed header row (`stats jaccard`: `intersection  union  jaccard ...`)
//! is marked in `BINDINGS`, and that row names the columns rather than being a record.
//!
//! Each step's output is taken from the golden cases: the CLI run on small fixtures, which
//! `golden.rs` checks. A step that prints an unprefixed header row without being marked, a
//! marked one that prints none, and a new step with no case to check all fail here.

use gloci::embed::{self, BINDINGS, Shape};
use gloci_core::VcfSpan;
use gloci_testkit::differential::split_args;
use std::ffi::OsString;
use std::path::Path;

/// A golden case whose output is text: its argv (`gloci` first) and gloci's output.
struct Case {
    name: String,
    argv: Vec<OsString>,
    out: Vec<u8>,
}

fn read(dir: &Path, file: &str) -> Option<String> {
    std::fs::read_to_string(dir.join(file))
        .ok()
        .map(|s| s.trim().to_string())
}

fn text_cases() -> Vec<Case> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    let mut cases: Vec<Case> = std::fs::read_dir(&root)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.is_dir())
        .filter(|dir| {
            // Binary output, and cases that check an exit status rather than the output.
            let compare = read(dir, "compare").unwrap_or_else(|| "exact".into());
            !matches!(compare.as_str(), "bam" | "cram" | "exit-status")
        })
        .map(|dir| {
            let args = split_args(&read(&dir, "args").unwrap()).unwrap();
            // gloci's own output where it differs from the reference tool's.
            let out = std::fs::read(dir.join("expected_gloci"))
                .or_else(|_| std::fs::read(dir.join("expected")))
                .unwrap();
            Case {
                name: dir.file_name().unwrap().to_string_lossy().into_owned(),
                argv: std::iter::once("gloci".to_string())
                    .chain(args)
                    .map(OsString::from)
                    .collect(),
                out,
            }
        })
        .collect();
    // Steps that write BAM only: their output, as `gloci view -O sam` shows it (as SAM text is
    // how it reads back).
    for (dir, args) in [
        ("cat-two-bams-vs-samtools", "cat input.bam cat-b.bam"),
        (
            "reheader-bam-vs-samtools",
            "reheader new-header.sam input.bam",
        ),
    ] {
        let gloci = || assert_cmd::Command::cargo_bin("gloci").unwrap();
        let bam = gloci()
            .current_dir(root.join(dir))
            .args(args.split(' '))
            .output()
            .unwrap();
        assert!(bam.status.success(), "gloci {args}: {bam:?}");
        let sam = gloci()
            .args(["view", "-O", "sam", "-"])
            .write_stdin(bam.stdout)
            .output()
            .unwrap();
        assert!(sam.status.success(), "gloci {args} | gloci view: {sam:?}");
        cases.push(Case {
            name: format!("{dir} (as SAM)"),
            argv: std::iter::once("gloci")
                .chain(args.split(' '))
                .map(OsString::from)
                .collect(),
            out: sam.stdout,
        });
    }
    cases.sort_by(|a, b| a.name.cmp(&b.name));
    cases
}

fn is_number(field: &[u8]) -> bool {
    std::str::from_utf8(field).is_ok_and(|s| s.parse::<f64>().is_ok())
}

/// The lines that are neither blank nor prefixed: `#` and `@` header lines, FASTA's `>`
/// lines, `track` and `browser` lines.
fn unprefixed(out: &[u8]) -> Vec<&[u8]> {
    out.split(|&b| b == b'\n')
        .map(|l| l.strip_suffix(b"\r").unwrap_or(l))
        .filter(|l| !l.iter().all(u8::is_ascii_whitespace))
        .filter(|l| {
            !(l.starts_with(b"#")
                || l.starts_with(b"@")
                || l.starts_with(b">")
                || l.starts_with(b"track")
                || l.starts_with(b"browser"))
        })
        .collect()
}

/// Whether the output starts with an unprefixed header row: a first line with no number in
/// it, above lines that have a number in one of its columns every time. `None` when there
/// are fewer than two lines to compare.
fn has_header_row(out: &[u8]) -> Option<bool> {
    let lines = unprefixed(out);
    let (first, rest) = lines.split_first()?;
    if rest.is_empty() {
        return None;
    }
    let fields = |l: &&[u8]| {
        l.split(|&b| b == b'\t')
            .map(<[u8]>::to_vec)
            .collect::<Vec<_>>()
    };
    let head = fields(first);
    if head.iter().any(|f| is_number(f)) {
        return Some(false);
    }
    let below: Vec<Vec<Vec<u8>>> = rest.iter().map(fields).collect();
    Some((0..head.len()).any(|i| below.iter().all(|r| r.get(i).is_some_and(|f| is_number(f)))))
}

fn command_of(argv: &[OsString]) -> Option<&'static embed::Binding> {
    let words: Vec<String> = argv[1..]
        .iter()
        .take(2)
        .map(|w| w.to_string_lossy().into_owned())
        .collect();
    (1..=words.len())
        .rev()
        .find_map(|n| embed::binding(&words[..n]))
}

#[test]
fn every_step_that_prints_a_header_row_is_marked() {
    let mut checked: Vec<&str> = Vec::new();
    let mut wrong = Vec::new();
    for case in text_cases() {
        let Some(b) = command_of(&case.argv) else {
            continue;
        };
        let back = embed::read_back(&case.argv)
            .unwrap_or_else(|e| panic!("{}: read_back: {}", case.name, e.0));
        let Some(seen) = has_header_row(&case.out) else {
            continue;
        };
        checked.push(b.method);
        let marked = back.shape == Shape::HeaderRow;
        if seen != marked {
            wrong.push(format!(
                "{} ({}): {}",
                case.name,
                b.method,
                if seen {
                    "prints an unprefixed header row, but its binding is not marked \
Shape::HeaderRow"
                } else {
                    "is marked Shape::HeaderRow, but prints no header row"
                }
            ));
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));

    // Every step's output was looked at, but dict's, which is SAM header lines only.
    let unchecked: Vec<&str> = BINDINGS
        .iter()
        .map(|b| b.method)
        .filter(|m| !checked.contains(m) && *m != "dict")
        .collect();
    assert!(
        unchecked.is_empty(),
        "no golden case shows two output lines of {unchecked:?}: add one, so that a header row \
would be seen"
    );
}

#[test]
fn the_marked_steps() {
    // Each step whose output is not plain records, with the flags that change it: tables whose
    // rows would read back as wrong intervals (or none), and tables with a header row.
    let marked: Vec<String> = BINDINGS
        .iter()
        .filter(|b| b.output != embed::RECORDS)
        .map(|b| {
            let with: Vec<String> = b
                .output
                .with
                .iter()
                .map(|(id, s)| format!("{id}: {s:?}"))
                .collect();
            let with = if with.is_empty() {
                String::new()
            } else {
                format!(" ({})", with.join(", "))
            };
            format!("{}: {:?}{with}", b.path.join(" "), b.output.shape)
        })
        .collect();
    assert_eq!(
        marked,
        [
            "coverage: Records (hist: Table)",
            "stats jaccard: HeaderRow",
            "stats fisher: HeaderRow",
            "stats reldist: HeaderRow (detail: Records)",
            "depth: Table",
            "samples: Table",
            "view: Records (count: Table)",
            "genomecov: Table (bedgraph: Records, bga: Records)",
            "vcf extract: HeaderRow",
            "getfasta: Records (tab: Table)",
            "flagstat: Table",
            "idxstats: Table",
        ]
    );
    // The flags are the command's own switches.
    for b in BINDINGS {
        let spec = embed::find_spec(b.path).unwrap();
        for (id, _) in b.output.with {
            assert!(
                spec.args.iter().any(|a| a.id == *id && !a.takes_value),
                "{}: no flag {id}",
                spec.name()
            );
        }
    }
}

fn back(args: &str) -> embed::ReadBack {
    let argv: Vec<OsString> = std::iter::once("gloci")
        .chain(args.split(' '))
        .map(OsString::from)
        .collect();
    embed::read_back(&argv).unwrap()
}

#[test]
fn read_back_follows_the_flags_given() {
    assert_eq!(
        back("stats reldist -a a.bed -b b.bed").shape,
        Shape::HeaderRow
    );
    assert_eq!(
        back("stats reldist -a a.bed -b b.bed --detail").shape,
        Shape::Records
    );
    assert_eq!(back("genomecov -g g.txt a.bed").shape, Shape::Table);
    assert_eq!(
        back("genomecov -g g.txt a.bed --per-base").shape,
        Shape::Table
    );
    assert_eq!(
        back("genomecov -g g.txt a.bed --bedgraph").shape,
        Shape::Records
    );
    assert_eq!(back("genomecov -g g.txt a.bed --bga").shape, Shape::Records);
    assert_eq!(back("intersect -a a.bed -b b.bed").shape, Shape::Records);
    // coverage --hist: A's lines with histogram columns, then `all` lines.
    assert_eq!(back("coverage -a a.bed -b b.bed").shape, Shape::Records);
    assert_eq!(
        back("coverage -a a.bed -b b.bed --hist").shape,
        Shape::Table
    );
    assert_eq!(
        back("coverage -a a.bed -b b.bed --per-base").shape,
        Shape::Records
    );
    // getfasta --tab: name and sequence; view --count: a number.
    assert_eq!(
        back("getfasta --fasta r.fa a.bed --tab").shape,
        Shape::Table
    );
    assert_eq!(
        back("getfasta --fasta r.fa a.bed --bed-out").shape,
        Shape::Records
    );
    assert_eq!(back("view x.bam --count").shape, Shape::Table);
    // groupby writes records only when its keys start with BED's chrom, start and end.
    assert_eq!(
        back("groupby a.bed -g 1,2,3 -c 5:sum").shape,
        Shape::Records
    );
    assert_eq!(
        back("groupby a.bed -g 1,2,3,6 -c 5:sum").shape,
        Shape::Records
    );
    assert_eq!(
        back("groupby a.bed -g 1 -c 5:sum -c 5:count").shape,
        Shape::Table
    );
    assert_eq!(back("groupby a.bed -g 1,2 -c 5:sum").shape, Shape::Table);
    assert_eq!(back("groupby a.bed -g 2,1,3 -c 5:sum").shape, Shape::Table);
    assert_eq!(
        back("groupby a.bed -g 6,1,2,3 -c 5:sum").shape,
        Shape::Table
    );
    // A module source is no step: it writes records.
    assert_eq!(back("windows -g g.txt -w 10").shape, Shape::Records);
    // The footprint of VCF records is the step's --vcf-span, else the allele span.
    assert_eq!(back("view x.vcf").vcf_span, VcfSpan::Allele);
    assert_eq!(back("view x.vcf --vcf-span ref").vcf_span, VcfSpan::Ref);
    assert_eq!(back("cluster x.vcf --vcf-span=ref").vcf_span, VcfSpan::Ref);
    assert_eq!(back("head x.vcf").vcf_span, VcfSpan::Allele);
    assert!(embed::read_back(&[OsString::from("gloci"), OsString::from("nope")]).is_err());
}
