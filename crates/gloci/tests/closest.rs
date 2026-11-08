//! `gloci closest`: distances, ties, `-k`, signed frames, strand rules,
//! `--each-side`, `--per-file`, `--max-distance`, VCF A, usage errors and
//! determinism.
mod common;
use common::*;
use predicates::str::contains;
use std::path::Path;

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

struct Fx {
    dir: tempfile::TempDir,
}

impl Fx {
    fn new() -> Fx {
        Fx { dir: tempdir() }
    }
    fn file(&self, name: &str, text: &str) -> String {
        s(&plain_file(self.dir.path(), name, text))
    }
    /// `closest -a a -b b extra...` over two BED texts, as output lines.
    fn run(&self, a: &str, b: &str, extra: &[&str]) -> Vec<String> {
        let (af, bf) = (self.file("a.bed", a), self.file("b.bed", b));
        let mut v = vec!["closest", "-a", &af, "-b", &bf];
        v.extend_from_slice(extra);
        let o = gloci().args(&v).assert().success();
        String::from_utf8(o.get_output().stdout.clone())
            .unwrap()
            .lines()
            .map(String::from)
            .collect()
    }
}

fn code(args: &[&str]) -> i32 {
    gloci().args(args).output().unwrap().status.code().unwrap()
}

const A: &str = "chr1\t100\t200\ta1\t0\t+\n";

/// The B names (4th column of the B part) of rows that have one.
fn b_names(rows: &[String], a_cols: usize) -> Vec<String> {
    rows.iter()
        .map(|r| r.split('\t').nth(a_cols + 3).unwrap().to_string())
        .collect()
}

#[test]
fn distances() {
    let fx = Fx::new();
    let one = |b: &str| fx.run(A, b, &["--distance"]);
    // overlap is 0
    assert_eq!(
        one("chr1\t150\t160\tb\t0\t+\n"),
        ["chr1\t100\t200\ta1\t0\t+\tchr1\t150\t160\tb\t0\t+\t0"]
    );
    // bookended is 1
    assert_eq!(
        one("chr1\t200\t210\tb\t0\t+\n"),
        ["chr1\t100\t200\ta1\t0\t+\tchr1\t200\t210\tb\t0\t+\t1"]
    );
    assert_eq!(
        one("chr1\t90\t100\tb\t0\t+\n"),
        ["chr1\t100\t200\ta1\t0\t+\tchr1\t90\t100\tb\t0\t+\t1"]
    );
    // a gap of 10 bases is 11
    assert_eq!(
        one("chr1\t210\t220\tb\t0\t+\n"),
        ["chr1\t100\t200\ta1\t0\t+\tchr1\t210\t220\tb\t0\t+\t11"]
    );
    // a point on A's edge does not overlap it: distance 1
    assert_eq!(
        one("chr1\t200\t200\tb\t0\t+\n"),
        ["chr1\t100\t200\ta1\t0\t+\tchr1\t200\t200\tb\t0\t+\t1"]
    );
    // with --inclusive-points the same point overlaps
    assert_eq!(
        fx.run(
            A,
            "chr1\t200\t200\tb\t0\t+\n",
            &["--distance", "--inclusive-points"]
        ),
        ["chr1\t100\t200\ta1\t0\t+\tchr1\t200\t200\tb\t0\t+\t0"]
    );
    // without --distance there is no distance column
    assert_eq!(
        fx.run(A, "chr1\t210\t220\tb\t0\t+\n", &[]),
        ["chr1\t100\t200\ta1\t0\t+\tchr1\t210\t220\tb\t0\t+"]
    );
}

#[test]
fn ties_and_k() {
    let fx = Fx::new();
    let b = "chr1\t50\t90\tL\t0\t+\nchr1\t210\t250\tR\t0\t+\nchr1\t900\t950\tfar\t0\t+\n";
    assert_eq!(b_names(&fx.run(A, b, &[]), 6), ["L", "R"]);
    assert_eq!(b_names(&fx.run(A, b, &["--ties", "all"]), 6), ["L", "R"]);
    assert_eq!(b_names(&fx.run(A, b, &["--ties", "first"]), 6), ["L"]);
    assert_eq!(b_names(&fx.run(A, b, &["--ties", "last"]), 6), ["R"]);

    // distances 1, 11, 11, far
    let b = "chr1\t200\t205\tx1\t0\t+\nchr1\t50\t90\tx2\t0\t+\nchr1\t210\t250\tx3\t0\t+\nchr1\t900\t950\tx4\t0\t+\n";
    // sorted by position: x2 first
    let b = {
        let mut l: Vec<&str> = b.lines().collect();
        l.sort_by_key(|x| x.split('\t').nth(1).unwrap().parse::<u64>().unwrap());
        l.join("\n") + "\n"
    };
    assert_eq!(b_names(&fx.run(A, &b, &["-k", "2"]), 6), ["x1", "x2", "x3"]);
    assert_eq!(
        b_names(&fx.run(A, &b, &["-k", "2", "--ties", "first"]), 6),
        ["x1", "x2"]
    );
    assert_eq!(
        b_names(&fx.run(A, &b, &["-k", "2", "--ties", "last"]), 6),
        ["x1", "x3"]
    );
    assert_eq!(
        b_names(&fx.run(A, &b, &["-k", "3", "--distance"]), 6),
        ["x1", "x2", "x3"]
    );
    assert_eq!(b_names(&fx.run(A, &b, &["-k", "10"]), 6).len(), 4);
    // an overlapping B is the nearest, then the rest
    let o = format!("{b}chr1\t150\t160\tov\t0\t+\n");
    let mut l: Vec<&str> = o.lines().collect();
    l.sort_by_key(|x| x.split('\t').nth(1).unwrap().parse::<u64>().unwrap());
    let o = l.join("\n") + "\n";
    assert_eq!(b_names(&fx.run(A, &o, &["-k", "2"]), 6), ["ov", "x1"]);
}

#[test]
fn no_candidate() {
    let fx = Fx::new();
    let a = "chr3\t5\t10\ta5\t0\t+\n";
    assert_eq!(
        fx.run(a, "chr1\t10\t50\tb\t0\t+\n", &["--distance"]),
        ["chr3\t5\t10\ta5\t0\t+\t.\t-1\t-1\t.\t-1\t.\t-1"]
    );
    assert_eq!(
        fx.run(a, "chr1\t10\t50\tb\t0\t+\n", &[]),
        ["chr3\t5\t10\ta5\t0\t+\t.\t-1\t-1\t.\t-1\t."]
    );
    // a BED3 B
    assert_eq!(
        fx.run(a, "chr1\t10\t50\n", &["--distance"]),
        ["chr3\t5\t10\ta5\t0\t+\t.\t-1\t-1\t-1"]
    );
}

#[test]
fn ignore_overlaps() {
    let fx = Fx::new();
    let b = "chr1\t50\t90\tL\t0\t+\nchr1\t150\t160\tov\t0\t+\nchr1\t210\t250\tR\t0\t+\n";
    assert_eq!(b_names(&fx.run(A, b, &[]), 6), ["ov"]);
    assert_eq!(
        b_names(&fx.run(A, b, &["--ignore-overlaps"]), 6),
        ["L", "R"]
    );
}

#[test]
fn signed_frames() {
    let fx = Fx::new();
    let b = "chr1\t50\t90\tL\t0\t+\nchr1\t210\t250\tR\t0\t+\n";
    let plus = "chr1\t100\t200\ta1\t0\t+\n";
    let minus = "chr1\t100\t200\ta1\t0\t-\n";
    let dist = |rows: Vec<String>| -> Vec<(String, String)> {
        rows.iter()
            .map(|r| {
                let f: Vec<&str> = r.split('\t').collect();
                (f[9].to_string(), f[12].to_string())
            })
            .collect()
    };
    let pair = |n: &str, d: &str| (n.to_string(), d.to_string());
    // ref: lower coordinates are upstream
    assert_eq!(
        dist(fx.run(plus, b, &["--signed", "ref"])),
        [pair("L", "-11"), pair("R", "11")]
    );
    assert_eq!(
        dist(fx.run(minus, b, &["--signed", "ref"])),
        [pair("L", "-11"), pair("R", "11")]
    );
    // a: upstream follows A's strand
    assert_eq!(
        dist(fx.run(plus, b, &["--signed", "a"])),
        [pair("L", "-11"), pair("R", "11")]
    );
    assert_eq!(
        dist(fx.run(minus, b, &["--signed", "a"])),
        [pair("L", "11"), pair("R", "-11")]
    );
    // --ignore-upstream --signed a on a minus-strand A drops R, not L
    assert_eq!(
        b_names(
            &fx.run(minus, b, &["--ignore-upstream", "--signed", "a"]),
            6
        ),
        ["L"]
    );
    assert_eq!(
        b_names(&fx.run(plus, b, &["--ignore-upstream", "--signed", "a"]), 6),
        ["R"]
    );
    assert_eq!(
        b_names(
            &fx.run(minus, b, &["--ignore-downstream", "--signed", "a"]),
            6
        ),
        ["R"]
    );
    // no --signed means ref
    assert_eq!(b_names(&fx.run(minus, b, &["--ignore-upstream"]), 6), ["R"]);
    // unknown strand counts as +
    assert_eq!(
        b_names(
            &fx.run(
                "chr1\t100\t200\ta1\t0\t.\n",
                b,
                &["--ignore-upstream", "--signed", "a"]
            ),
            6
        ),
        ["R"]
    );
    // --signed alone prints the signed distance
    assert_eq!(
        fx.run(plus, "chr1\t50\t90\tL\t0\t+\n", &["--signed", "ref"]),
        ["chr1\t100\t200\ta1\t0\t+\tchr1\t50\t90\tL\t0\t+\t-11"]
    );
    // an overlap is neither upstream nor downstream
    let ov = "chr1\t150\t160\tov\t0\t+\n";
    assert_eq!(
        b_names(
            &fx.run(plus, ov, &["--ignore-upstream", "--ignore-downstream"]),
            6
        ),
        ["ov"]
    );
}

/// `--signed b`: whether a B is upstream depends on B's own strand, so the
/// nearest upstream-filtered B may be far behind a nearer one.
#[test]
fn signed_b_frame_filters_per_b_strand() {
    let fx = Fx::new();
    // left: nearer `-` (B2, d=11), farther `+` (B1, d=31).
    // right: nearer `+` (B3, d=11), farther `-` (B4, d=41).
    let b = "chr1\t10\t70\tB1\t0\t+\n\
             chr1\t80\t90\tB2\t0\t-\n\
             chr1\t210\t220\tB3\t0\t+\n\
             chr1\t240\t250\tB4\t0\t-\n";
    let rows = |extra: &[&str]| -> Vec<(String, String)> {
        fx.run(A, b, extra)
            .iter()
            .map(|r| {
                let f: Vec<&str> = r.split('\t').collect();
                (f[9].to_string(), f[12].to_string())
            })
            .collect()
    };
    let pair = |n: &str, d: &str| (n.to_string(), d.to_string());
    // `+` lower B is downstream (+), `-` lower B is upstream (-);
    // `+` higher B is upstream (-), `-` higher B is downstream (+)
    assert_eq!(
        rows(&["--signed", "b", "-k", "4"]),
        [
            pair("B2", "-11"),
            pair("B3", "-11"),
            pair("B1", "31"),
            pair("B4", "41"),
        ]
    );
    assert_eq!(
        rows(&["--signed", "b", "--ignore-upstream"]),
        [pair("B1", "31")]
    );
    assert_eq!(
        rows(&["--signed", "b", "--ignore-downstream"]),
        [pair("B2", "-11"), pair("B3", "-11")]
    );
}

#[test]
fn each_side() {
    let fx = Fx::new();
    let b = "chr1\t10\t20\tu2\t0\t+\nchr1\t50\t90\tu1\t0\t+\nchr1\t210\t250\td1\t0\t+\nchr1\t400\t500\td2\t0\t+\n";
    assert_eq!(
        fx.run(A, b, &["--each-side", "--distance"]),
        ["chr1\t100\t200\ta1\t0\t+\tchr1\t50\t90\tu1\t0\t+\t-11\tchr1\t210\t250\td1\t0\t+\t11"]
    );
    assert_eq!(
        fx.run(A, b, &["--each-side"]),
        ["chr1\t100\t200\ta1\t0\t+\tchr1\t50\t90\tu1\t0\t+\tchr1\t210\t250\td1\t0\t+"]
    );
    // -k applies per side: one row per rank
    let rows = fx.run(A, b, &["--each-side", "-k", "2", "--distance"]);
    assert_eq!(
        rows,
        [
            "chr1\t100\t200\ta1\t0\t+\tchr1\t50\t90\tu1\t0\t+\t-11\tchr1\t210\t250\td1\t0\t+\t11",
            "chr1\t100\t200\ta1\t0\t+\tchr1\t10\t20\tu2\t0\t+\t-81\tchr1\t400\t500\td2\t0\t+\t201",
        ]
    );
    // an overlapping B counts on both sides
    let ov = format!("{b}chr1\t150\t160\tov\t0\t+\n");
    let mut l: Vec<&str> = ov.lines().collect();
    l.sort_by_key(|x| x.split('\t').nth(1).unwrap().parse::<u64>().unwrap());
    let ov = l.join("\n") + "\n";
    assert_eq!(
        fx.run(A, &ov, &["--each-side", "--distance"]),
        ["chr1\t100\t200\ta1\t0\t+\tchr1\t150\t160\tov\t0\t+\t0\tchr1\t150\t160\tov\t0\t+\t0"]
    );
    assert_eq!(
        fx.run(A, &ov, &["--each-side", "--distance", "--ignore-overlaps"]),
        fx.run(A, b, &["--each-side", "--distance"])
    );
    // a missing side is a placeholder with distance -1
    assert_eq!(
        fx.run(
            A,
            "chr1\t50\t90\tu1\t0\t+\n",
            &["--each-side", "--distance"]
        ),
        ["chr1\t100\t200\ta1\t0\t+\tchr1\t50\t90\tu1\t0\t+\t-11\t.\t-1\t-1\t.\t-1\t.\t-1"]
    );
    // the upstream side follows the --signed frame
    assert_eq!(
        fx.run(
            "chr1\t100\t200\ta1\t0\t-\n",
            b,
            &["--each-side", "--distance", "--signed", "a"]
        ),
        ["chr1\t100\t200\ta1\t0\t-\tchr1\t210\t250\td1\t0\t+\t-11\tchr1\t50\t90\tu1\t0\t+\t11"]
    );
}

#[test]
fn per_file() {
    let fx = Fx::new();
    let b1 = fx.file("b1.bed", "chr1\t50\t90\tx\t0\t+\nchr1\t210\t250\ty\t0\t+\n");
    let b2 = fx.file("b2.bed", "chr1\t120\t130\tz\t0\t+\nchr2\t5\t6\tw\t0\t+\n");
    let a = fx.file("a.bed", "chr1\t100\t200\ta1\t0\t+\nchr3\t0\t5\ta2\t0\t+\n");
    let run = |extra: &[&str]| -> Vec<String> {
        let mut v = vec!["closest", "-a", &a, "-b", &b1, &b2, "--names", "F,G"];
        v.extend_from_slice(extra);
        let o = gloci().args(&v).assert().success();
        String::from_utf8(o.get_output().stdout.clone())
            .unwrap()
            .lines()
            .map(String::from)
            .collect()
    };
    // across files: the overlap in the second file wins; label column
    assert_eq!(
        run(&["--distance"]),
        [
            "chr1\t100\t200\ta1\t0\t+\tG\tchr1\t120\t130\tz\t0\t+\t0",
            "chr3\t0\t5\ta2\t0\t+\t.\t.\t-1\t-1\t.\t-1\t.\t-1",
        ]
    );
    // per file: each file's rows in file order
    assert_eq!(
        run(&["--distance", "--per-file"]),
        [
            "chr1\t100\t200\ta1\t0\t+\tF\tchr1\t50\t90\tx\t0\t+\t11",
            "chr1\t100\t200\ta1\t0\t+\tF\tchr1\t210\t250\ty\t0\t+\t11",
            "chr1\t100\t200\ta1\t0\t+\tG\tchr1\t120\t130\tz\t0\t+\t0",
            "chr3\t0\t5\ta2\t0\t+\t.\t.\t-1\t-1\t.\t-1\t.\t-1",
        ]
    );
    // default labels are the file basenames
    let o = gloci()
        .args(["closest", "-a", &a, "-b", &b1, &b2, "--per-file"])
        .assert()
        .success();
    let out = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    assert!(out.contains("\tb1.bed\tchr1\t50\t90"), "{out}");
    assert!(out.contains("\tb2.bed\tchr1\t120\t130"), "{out}");
}

#[test]
fn max_distance() {
    let fx = Fx::new();
    let b = "chr1\t210\t220\tb\t0\t+\n";
    let none = "chr1\t100\t200\ta1\t0\t+\t.\t-1\t-1\t.\t-1\t.\t-1";
    assert_eq!(
        fx.run(A, b, &["--max-distance", "10", "--distance"]),
        [none]
    );
    assert_eq!(
        fx.run(A, b, &["--max-distance", "11", "--distance"]),
        ["chr1\t100\t200\ta1\t0\t+\tchr1\t210\t220\tb\t0\t+\t11"]
    );
    // an overlap is always within reach
    assert_eq!(
        b_names(
            &fx.run(A, "chr1\t150\t160\tb\t0\t+\n", &["--max-distance", "0"]),
            6
        ),
        ["b"]
    );
}

#[test]
fn strand_rules() {
    let fx = Fx::new();
    // the nearer B is on the opposite strand
    let b = "chr1\t50\t90\tnear\t0\t-\nchr1\t210\t250\tfar\t0\t+\nchr1\t260\t270\tfarther\t0\t-\n";
    assert_eq!(b_names(&fx.run(A, b, &[]), 6), ["near", "far"]);
    assert_eq!(b_names(&fx.run(A, b, &["-s"]), 6), ["far"]);
    assert_eq!(b_names(&fx.run(A, b, &["-S"]), 6), ["near"]);
    // unknown strand never matches under a rule
    let unk = "chr1\t100\t200\ta1\t0\t.\n";
    assert_eq!(
        fx.run(unk, b, &["-s", "--distance"]),
        ["chr1\t100\t200\ta1\t0\t.\t.\t-1\t-1\t.\t-1\t.\t-1"]
    );
}

#[test]
fn unsorted_matches_sorted() {
    let fx = Fx::new();
    let a = "chr1\t100\t200\ta1\t0\t+\nchr1\t5000\t5100\ta2\t0\t-\n";
    let b = "chr1\t50\t90\tb1\t0\t+\nchr1\t210\t250\tb2\t0\t-\nchr1\t4000\t4100\tb3\t0\t+\n";
    let shuffled = "chr1\t4000\t4100\tb3\t0\t+\nchr1\t50\t90\tb1\t0\t+\nchr1\t210\t250\tb2\t0\t-\n";
    for extra in [
        &["--distance"][..],
        &["-k", "2", "--distance"],
        &["--each-side", "--distance"],
    ] {
        let want = fx.run(a, b, extra);
        let mut e = extra.to_vec();
        e.push("--unsorted");
        // B in input order differs, so compare as sets of rows per A
        let mut got = fx.run(a, shuffled, &e);
        let mut want = want;
        got.sort();
        want.sort();
        assert_eq!(got, want, "{extra:?}");
    }
    // sorted input is rejected without --unsorted
    assert_eq!(
        code(&[
            "closest",
            "-a",
            &fx.file("a2.bed", a),
            "-b",
            &fx.file("b2.bed", shuffled)
        ]),
        1
    );
}

#[test]
fn vcf_a_gets_a_column_line() {
    let fx = Fx::new();
    let a = fx.file(
        "a.vcf",
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t150\t.\tA\tC\t.\t.\t.\n",
    );
    let b = fx.file("b.bed", "chr1\t10\t50\tb\t0\t+\n");
    let o = gloci()
        .args(["closest", "-a", &a, "-b", &b, "--distance"])
        .assert()
        .success();
    let out = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    assert_eq!(
        out,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tb_1\tb_2\tb_3\tb_4\tb_5\tb_6\tdistance\n\
         chr1\t150\t.\tA\tC\t.\t.\t.\tchr1\t10\t50\tb\t0\t+\t100\n"
    );
}

#[test]
fn usage_errors() {
    let fx = Fx::new();
    let (a, b) = (fx.file("a.bed", A), fx.file("b.bed", "chr1\t1\t2\n"));
    assert_eq!(code(&["closest", "-a", &a, "-b", &b, "-k", "0"]), 2);
    assert_eq!(code(&["closest", "-a", &a, "-b", &b, "--ties", "some"]), 2);
    assert_eq!(code(&["closest", "-a", &a, "-b", &b, "--signed", "x"]), 2);
    assert_eq!(code(&["closest", "-a", &a, "-b", &b, "-w", "5"]), 2);
    assert_eq!(code(&["closest", "-a", &a, "-b", &b, "-f", "0.5"]), 2);
    gloci()
        .args(["closest", "-a", &a, "-b", &b, "-k", "0"])
        .assert()
        .code(2)
        .stderr(contains("-k"));
}

#[test]
fn help_documents_the_strand_look_ahead() {
    let o = gloci().args(["closest", "--help"]).assert().success();
    let out = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    assert!(out.contains("--max-distance"), "{out}");
    assert!(out.contains("Examples:"), "{out}");
    assert!(
        out.contains("remaining") && out.contains("buffer"),
        "the look-ahead pathology is documented: {out}"
    );
    assert!(out.contains("--unsorted"), "{out}");
}

#[test]
fn determinism() {
    let fx = Fx::new();
    let a = fx.file(
        "a.bed",
        "chr1\t100\t200\ta1\t0\t+\nchr1\t300\t400\ta2\t0\t-\nchr2\t5\t15\ta3\t0\t+\nchr3\t0\t9\ta4\t0\t+\n",
    );
    let b = fx.file(
        "b.bed",
        "chr1\t10\t20\tb1\t0\t+\nchr1\t150\t160\tb2\t0\t-\nchr1\t450\t500\tb3\t0\t+\nchr2\t0\t8\tb4\t0\t-\n",
    );
    let b2 = fx.file("b2.bed", "chr1\t500\t600\tc1\t0\t+\n");
    for extra in [
        &[][..],
        &["--distance", "-k", "2"],
        &["--each-side", "--distance"],
        &["--per-file", "--signed", "a"],
        &["-s", "--distance"],
        &["--unsorted", "--distance"],
    ] {
        let run = |t: &str| {
            let mut v = vec!["closest", "-a", &a, "-b", &b, &b2, "-t", t];
            v.extend_from_slice(extra);
            gloci()
                .args(&v)
                .assert()
                .success()
                .get_output()
                .stdout
                .clone()
        };
        let one = run("1");
        assert!(!one.is_empty());
        assert_eq!(one, run("2"), "{extra:?}");
        assert_eq!(one, run("8"), "{extra:?}");
    }
}
