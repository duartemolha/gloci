//! `gloci tag` (D168): formats, option errors, threads, stdin, indexes. The output itself is
//! checked against bedtools by the `tag-*-vs-bedtools` golden cases, whose fixtures are here.
mod common;
use common::*;
use predicates::prelude::PredicateBooleanExt;
use predicates::str::contains;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs/tag")
}

fn f(name: &str) -> String {
    dir().join(name).display().to_string()
}

fn s(p: &Path) -> String {
    p.display().to_string()
}

/// `gloci tag ARGS` that must succeed: its standard output.
fn tag(args: &[&str]) -> Vec<u8> {
    let o = gloci().arg("tag").args(args).output().unwrap();
    assert!(
        o.status.success(),
        "gloci tag {args:?}: {}",
        String::from_utf8_lossy(&o.stderr)
    );
    o.stdout
}

/// SAM text without its `@PG` lines (gloci's names the command line); BAM is read back first.
fn sam_text(bytes: &[u8]) -> String {
    let text = if bytes.starts_with(&[0x1f, 0x8b]) {
        let o = gloci()
            .args(["view", "-O", "sam", "-"])
            .write_stdin(bytes.to_vec())
            .output()
            .unwrap();
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
        o.stdout
    } else {
        bytes.to_vec()
    };
    String::from_utf8(text)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with("@PG\t"))
        .map(|l| format!("{l}\n"))
        .collect()
}

fn records(text: &str) -> Vec<&str> {
    text.lines().filter(|l| !l.starts_with('@')).collect()
}

/// A SAM of `n` reads on three references, in no order, with CIGARs that have introns,
/// deletions and clips, both strands, some unmapped; and a BED6 of features over them, some
/// zero-length, some unknown-strand, some long.
fn many(dir: &Path, n: usize) -> (PathBuf, PathBuf, PathBuf) {
    let mut seed = 0x9E37_79B9_7F4A_7C15u64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed
    };
    let refs = [("chr1", 400_000u64), ("chr2", 200_000), ("chr3", 50_000)];
    let mut sam = String::from("@HD\tVN:1.6\n");
    for (c, len) in refs {
        writeln!(sam, "@SQ\tSN:{c}\tLN:{len}").unwrap();
    }
    let cigars = [
        "50M",
        "20M300N30M",
        "10S40M",
        "25M5D25M",
        "30M2I18M",
        "5S40M5S",
    ];
    for i in 0..n {
        let r = next();
        let (c, len) = refs[(r % 3) as usize];
        let pos = 1 + (r >> 8) % (len - 1000);
        let cigar = cigars[(r >> 40) as usize % cigars.len()];
        let flag = match (r >> 50) % 10 {
            0 => 4,
            1 | 2 => 16,
            3 => 256,
            _ => 0,
        };
        writeln!(
            sam,
            "q{i}\t{flag}\t{c}\t{pos}\t60\t{cigar}\t*\t0\t0\t{}\t{}",
            "A".repeat(50),
            "I".repeat(50)
        )
        .unwrap();
    }
    let mut bed = String::new();
    for i in 0..n / 4 {
        let r = next();
        let (c, len) = refs[(r % 3) as usize];
        let start = (r >> 8) % (len - 40_000);
        let width = match (r >> 40) % 8 {
            0 => 0,
            1 => 30_000,
            _ => 1 + (r >> 44) % 400,
        };
        let strand = ["+", "-", "."][(r >> 52) as usize % 3];
        writeln!(
            bed,
            "{c}\t{start}\t{}\tf{i}\t{}\t{strand}",
            start + width,
            r % 1000
        )
        .unwrap();
    }
    let sam_path = plain_file(dir, "many.sam", &sam);
    let bam_path = dir.join("many.bam");
    gloci()
        .args(["view", &s(&sam_path), "-o", &s(&bam_path)])
        .assert()
        .success();
    (sam_path, bam_path, plain_file(dir, "many.bed", &bed))
}

#[test]
fn output_is_the_same_at_every_thread_count() {
    let t = tempdir();
    let (sam, bam, bed) = many(t.path(), 8_000);
    for input in [&sam, &bam] {
        for mode in [
            vec!["--labels", "f"],
            vec!["--names", "-S"],
            vec!["--intervals", "--labels", "f"],
        ] {
            let mut first: Option<(String, String)> = None;
            for th in ["1", "2", "8"] {
                let out = t.path().join(format!("o-{th}.bam"));
                let mut args = vec![s(input), "-b".into(), s(&bed), "-t".into(), th.into()];
                args.extend(mode.iter().map(|m| m.to_string()));
                let argv: Vec<&str> = args.iter().map(String::as_str).collect();
                let sam_out = sam_text(&tag(&[&argv[..], &["-O", "sam"]].concat()));
                tag(&[&argv[..], &["-o", &s(&out)]].concat());
                let bam_out = sam_text(&std::fs::read(&out).unwrap());
                assert_eq!(records(&sam_out).len(), 8_000);
                assert!(sam_out.contains("YB:Z:"), "{mode:?}: nothing tagged");
                match &first {
                    None => first = Some((sam_out, bam_out)),
                    Some((s1, b1)) => {
                        assert!(*s1 == sam_out, "{input:?} {mode:?} -t {th}: SAM differs");
                        assert!(*b1 == bam_out, "{input:?} {mode:?} -t {th}: BAM differs");
                    }
                }
            }
            // The BAM output reads back as the SAM output.
            let (s1, b1) = first.unwrap();
            assert!(s1 == b1, "{input:?} {mode:?}: BAM and SAM outputs differ");
        }
    }
}

#[test]
fn sam_and_bam_inputs_give_the_same_records() {
    let t = tempdir();
    let (sam, bam, bed) = many(t.path(), 3000);
    let a = sam_text(&tag(&[&s(&sam), "-b", &s(&bed), "--names"]));
    let b = sam_text(&tag(&[&s(&bam), "-b", &s(&bed), "--names", "-O", "sam"]));
    assert_eq!(records(&a), records(&b));
}

#[test]
fn output_format_follows_the_input_unless_chosen() {
    let t = tempdir();
    let args = |input: &str| {
        vec![
            input.to_string(),
            "-b".into(),
            f("genes.bed"),
            "--labels".into(),
            "G".into(),
        ]
    };
    // BAM in: BAM out; SAM in: SAM out; CRAM in: BAM out.
    let run = |input: &str, extra: &[&str]| {
        let mut a = args(input);
        a.extend(extra.iter().map(|e| e.to_string()));
        let argv: Vec<&str> = a.iter().map(String::as_str).collect();
        tag(&argv)
    };
    assert!(run(&f("reads.bam"), &[]).starts_with(&[0x1f, 0x8b]));
    assert!(run(&f("reads.sam"), &[]).starts_with(b"@HD\t"));
    assert!(run(&f("reads.cram"), &[]).starts_with(&[0x1f, 0x8b]));
    assert!(run(&f("reads.bam"), &["-O", "sam"]).starts_with(b"@HD\t"));
    assert!(run(&f("reads.sam"), &["-O", "bam"]).starts_with(&[0x1f, 0x8b]));
    // -o picks the format from its extension.
    for (name, magic) in [("x.sam", &b"@HD\t"[..]), ("x.bam", &[0x1f, 0x8b][..])] {
        let out = t.path().join(name);
        run(&f("reads.bam"), &["-o", &s(&out)]);
        assert!(std::fs::read(&out).unwrap().starts_with(magic), "{name}");
    }
    // Neither CRAM nor a non-alignment format can be written.
    for extra in [vec!["-o", "x.cram"], vec!["-O", "bed"], vec!["-o", "x.bed"]] {
        let mut a = args(&f("reads.bam"));
        a.extend(extra.iter().map(|e| e.to_string()));
        gloci()
            .current_dir(t.path())
            .arg("tag")
            .args(&a)
            .assert()
            .code(2)
            .stderr(contains("BAM or SAM").or(contains("CRAM")));
    }
    assert!(!t.path().join("x.cram").exists() && !t.path().join("x.bed").exists());
}

#[test]
fn option_mistakes_are_usage_errors_and_write_nothing() {
    let t = tempdir();
    let out = s(&t.path().join("out.bam"));
    let bam = f("reads.bam");
    let genes = f("genes.bed");
    let vcf = f("calls.vcf");
    let cases: Vec<(Vec<&str>, &str)> = vec![
        (
            vec!["-b", &genes, "--labels", "G", "--tag", "X"],
            "not a tag name",
        ),
        (
            vec!["-b", &genes, "--labels", "G", "--tag", "XYZ"],
            "not a tag name",
        ),
        (
            vec!["-b", &genes, "--labels", "G", "--tag", "1X"],
            "not a tag name",
        ),
        (
            vec!["-b", &genes, "--labels", "G,H"],
            "--labels gives 2 label(s) but there are 1 -b file(s)",
        ),
        (
            vec!["-b", &genes, &genes, "--labels", "G"],
            "--labels gives 1 label(s) but there are 2",
        ),
        (vec!["-b", &genes, "--labels", "G", "-f", "0"], "0"),
        (vec!["-b", &genes, "--labels", "G", "-f", "1.5"], "1.5"),
        (
            vec!["-b", &genes, "--labels", "G", "--fraction", "abc"],
            "abc",
        ),
        (vec!["-b", &genes, "--names", "--labels", "G"], "--labels"),
        (vec!["-b", &genes, "--scores", "--labels", "G"], "--labels"),
        (vec!["-b", &genes, "--names", "--scores"], "--scores"),
        (vec!["-b", &genes, "--intervals"], "--labels"),
        (vec!["-b", &genes], "--labels"),
        (
            vec!["-b", &genes, "--labels", "G", "-s", "-S"],
            "--opposite-strand",
        ),
        (vec!["-b", &vcf, "--labels", "V"], "gloci convert"),
        (vec!["-b", &bam, "--labels", "B"], "alignments"),
    ];
    for (args, msg) in cases {
        gloci()
            .arg("tag")
            .arg(&bam)
            .args(&args)
            .args(["-o", &out])
            .assert()
            .code(2)
            .stderr(contains(msg));
        assert!(!Path::new(&out).exists(), "{args:?} wrote {out}");
    }
}

#[test]
fn an_error_in_a_feature_file_writes_nothing() {
    let t = tempdir();
    let bad = plain_file(
        t.path(),
        "bad.bed",
        "chr1\t100\t150\tok\nchr1\t300\t200\tbackwards\n",
    );
    let out = t.path().join("out.bam");
    gloci()
        .args([
            "tag",
            &f("reads.bam"),
            "-b",
            &s(&bad),
            "--names",
            "-o",
            &s(&out),
        ])
        .assert()
        .code(1)
        .stderr(contains("bad.bed").and(contains("line 2")));
    assert!(!out.exists());
    gloci()
        .args([
            "tag",
            &f("reads.bam"),
            "-b",
            "missing.bed",
            "--names",
            "-o",
            &s(&out),
        ])
        .current_dir(t.path())
        .assert()
        .code(1)
        .stderr(contains("missing.bed"));
    assert!(!out.exists());
}

#[test]
fn a_fraction_may_be_a_percentage() {
    let base = [
        f("reads.bam"),
        "-b".into(),
        f("genes.bed"),
        "--names".into(),
        "-O".into(),
        "sam".into(),
    ];
    let run = |frac: &str| {
        let mut a: Vec<&str> = base.iter().map(String::as_str).collect();
        a.extend(["-f", frac]);
        sam_text(&tag(&a))
    };
    assert_eq!(run("50%"), run("0.5"));
    assert_ne!(run("50%"), run("51%"));
}

#[test]
fn cram_input_gives_what_its_records_as_bam_give() {
    let a = sam_text(&tag(&[
        &f("reads.cram"),
        "-b",
        &f("genes.bed"),
        "--names",
        "-O",
        "sam",
    ]));
    let b = sam_text(&tag(&[
        &f("reads-cram.bam"),
        "-b",
        &f("genes.bed"),
        "--names",
        "-O",
        "sam",
    ]));
    assert_eq!(records(&a), records(&b));
}

#[test]
fn standard_input_is_read_as_a_file_is() {
    let file = tag(&[
        &f("reads.bam"),
        "-b",
        &f("genes.bed"),
        "--intervals",
        "--labels",
        "G",
        "-O",
        "sam",
    ]);
    let o = gloci()
        .args([
            "tag",
            "-",
            "-b",
            &f("genes.bed"),
            "--intervals",
            "--labels",
            "G",
            "-O",
            "sam",
        ])
        .write_stdin(std::fs::read(f("reads.bam")).unwrap())
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert_eq!(records(&sam_text(&o.stdout)), records(&sam_text(&file)));
}

#[test]
fn compressed_features_and_several_files() {
    let t = tempdir();
    let text = std::fs::read_to_string(f("genes.bed")).unwrap();
    let gz = bgzf_file(t.path(), "genes.bed.gz", &text);
    let plain = tag(&[
        &f("reads.bam"),
        "-b",
        &f("genes.bed"),
        &f("genes2.bed"),
        "--names",
        "-O",
        "sam",
    ]);
    let packed = tag(&[
        &f("reads.bam"),
        "-b",
        &s(&gz),
        "-b",
        &f("genes2.bed"),
        "--names",
        "-O",
        "sam",
    ]);
    assert_eq!(records(&sam_text(&plain)), records(&sam_text(&packed)));
}

#[test]
fn write_index_needs_sorted_reads() {
    let t = tempdir();
    let sorted = t.path().join("sorted.bam");
    gloci()
        .args(["sort", &f("reads.bam"), "-o", &s(&sorted)])
        .assert()
        .success();
    let out = t.path().join("tagged.bam");
    tag(&[
        &s(&sorted),
        "-b",
        &f("genes.bed"),
        "--names",
        "-o",
        &s(&out),
        "--write-index",
    ]);
    assert!(t.path().join("tagged.bam.bai").exists());
    // The tagged file is the sorted one with the tags: the index finds its reads.
    let o = gloci()
        .args([
            "view",
            "-O",
            "sam",
            "--no-header",
            &s(&out),
            "chr2:5500-5600",
        ])
        .output()
        .unwrap();
    assert!(String::from_utf8_lossy(&o.stdout).contains("YB:Z:P,Q,L0"));
    // reads.bam is in no order: the index cannot be written, and neither is the BAM.
    let bad = t.path().join("unsorted.bam");
    gloci()
        .args([
            "tag",
            &f("reads.bam"),
            "-b",
            &f("genes.bed"),
            "--names",
            "-o",
            &s(&bad),
            "--write-index",
        ])
        .assert()
        .code(1)
        .stderr(contains("sort"));
    assert!(!bad.exists());
}

/// The tag of each read named in `names`, from `gloci tag reads.bam ARGS -O sam` (empty when the
/// read has none), and the standard error.
fn tags_of(args: &[&str], names: &[&str]) -> (Vec<String>, String) {
    let o = gloci()
        .arg("tag")
        .arg(f("reads.bam"))
        .args(args)
        .args(["-O", "sam"])
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let text = String::from_utf8(o.stdout).unwrap();
    let tags = names
        .iter()
        .map(|n| {
            let line = records(&text)
                .into_iter()
                .find(|l| l.split('\t').next() == Some(n))
                .unwrap_or_else(|| panic!("no read {n}"));
            line.split('\t')
                .skip(11)
                .filter(|t| t.starts_with("YB:Z:"))
                .collect::<Vec<_>>()
                .join(" ")
        })
        .collect();
    (tags, String::from_utf8(o.stderr).unwrap())
}

/// gloci's strand rule (spec section 4.3), not bedtools tag's: a feature whose strand is not +
/// or - (`.`, `?`, or a BED3 with no strand column) matches neither -s nor -S, and each file that
/// has such features gets one warning. `reverse` (16) and `forward` (0) sit on four copies of one
/// feature in strand.bed (+, -, ., ?) and on one in plain.bed (BED3).
#[test]
fn unknown_strands_match_neither_strand() {
    let reads = ["reverse", "forward", "r1"];
    let files = [f("strand.bed"), f("plain.bed")];
    let base = ["-b", files[0].as_str(), files[1].as_str(), "--names"];
    let (same, err) = tags_of(&[&base[..], &["-s"]].concat(), &reads);
    assert_eq!(same, ["YB:Z:minus", "YB:Z:plus", ""]);
    assert!(
        err.contains("strand.bed has records with unknown strand; they never match under -s or -S"),
        "{err}"
    );
    assert!(
        err.contains("plain.bed has records with unknown strand"),
        "{err}"
    );
    let (opposite, _) = tags_of(&[&base[..], &["-S"]].concat(), &reads);
    assert_eq!(opposite, ["YB:Z:plus", "YB:Z:minus", ""]);
    // Without a strand rule every copy counts, and nothing is said about strands.
    let (any, err) = tags_of(&base, &reads);
    assert_eq!(
        any,
        [
            "YB:Z:plus,minus,dot,query;",
            "YB:Z:plus,minus,dot,query;",
            "YB:Z:"
        ]
    );
    assert!(!err.contains("unknown strand"), "{err}");
    // Labels: plain.bed (strand unknown) never counts under -S, so no read gets P.
    let (labels, _) = tags_of(
        &[
            "-b",
            files[0].as_str(),
            files[1].as_str(),
            "--labels",
            "S,P",
            "-S",
        ],
        &reads,
    );
    assert_eq!(labels, ["YB:Z:S", "YB:Z:S", ""]);
}

#[test]
fn standard_input_is_read_once() {
    let bam = f("reads.bam");
    for args in [
        vec!["-", "-b", "-", "--labels", "x"],
        vec![bam.as_str(), "-b", "-", "-", "--names"],
    ] {
        gloci()
            .arg("tag")
            .args(&args)
            .write_stdin("chr1\t1\t2\n")
            .assert()
            .code(2)
            .stdout("")
            .stderr(contains("standard input for at most one of its files"));
    }
}
