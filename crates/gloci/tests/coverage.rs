//! `gloci coverage`: summary, `--hist` with its `all` rows, `--per-base`,
//! `--mean`, several `-b` files (bedtools annotate), points (D56),
//! thresholds, `--split`, VCF A, usage errors and determinism.
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
    /// `coverage -a a -b b... extra...` over BED texts, as output lines.
    fn run(&self, a: &str, bs: &[&str], extra: &[&str]) -> Vec<String> {
        let af = self.file("a.bed", a);
        let bf: Vec<String> = bs
            .iter()
            .enumerate()
            .map(|(i, t)| self.file(&format!("b{i}.bed"), t))
            .collect();
        let mut v = vec!["coverage", "-a", &af, "-b"];
        v.extend(bf.iter().map(String::as_str));
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

fn f(n: u64, d: u64) -> String {
    (n as f64 / d as f64).to_string()
}

const A: &str = "chr1\t100\t200\ta1\t0\t+\nchr1\t300\t400\ta2\t0\t-\nchr2\t5\t15\ta3\t0\t+\nchr3\t10\t10\tpt\t0\t+\n";
const B: &str = "chr1\t150\t160\tb1\t2\t+\nchr1\t180\t260\tb2\t4\t-\nchr1\t395\t500\tb3\t1\t+\nchr3\t5\t15\tb4\t1\t+\n";

#[test]
fn summary() {
    let fx = Fx::new();
    assert_eq!(
        fx.run(A, &[B], &[]),
        [
            "chr1\t100\t200\ta1\t0\t+\t2\t30\t100\t0.3".to_string(),
            "chr1\t300\t400\ta2\t0\t-\t1\t5\t100\t0.05".to_string(),
            "chr2\t5\t15\ta3\t0\t+\t0\t0\t10\t0".to_string(),
            "chr3\t10\t10\tpt\t0\t+\t1\t0\t0\t.".to_string(),
        ]
    );
}

#[test]
fn hist_with_all() {
    let fx = Fx::new();
    let rows = fx.run(A, &[B], &["--hist"]);
    let want = vec![
        format!("chr1\t100\t200\ta1\t0\t+\t0\t70\t100\t{}", f(70, 100)),
        format!("chr1\t100\t200\ta1\t0\t+\t1\t30\t100\t{}", f(30, 100)),
        format!("chr1\t300\t400\ta2\t0\t-\t0\t95\t100\t{}", f(95, 100)),
        format!("chr1\t300\t400\ta2\t0\t-\t1\t5\t100\t{}", f(5, 100)),
        "chr2\t5\t15\ta3\t0\t+\t0\t10\t10\t1".to_string(),
        // the point has no rows and adds nothing to the totals
        format!("all\t0\t175\t210\t{}", f(175, 210)),
        format!("all\t1\t35\t210\t{}", f(35, 210)),
    ];
    assert_eq!(rows, want);
}

#[test]
fn hist_full_cover_has_no_depth_zero_row() {
    let fx = Fx::new();
    let rows = fx.run("chr1\t0\t10\n", &["chr1\t0\t6\nchr1\t4\t12\n"], &["--hist"]);
    assert_eq!(
        rows,
        [
            format!("chr1\t0\t10\t1\t8\t10\t{}", f(8, 10)),
            format!("chr1\t0\t10\t2\t2\t10\t{}", f(2, 10)),
            format!("all\t1\t8\t10\t{}", f(8, 10)),
            format!("all\t2\t2\t10\t{}", f(2, 10)),
        ]
    );
}

#[test]
fn per_base() {
    let fx = Fx::new();
    let rows = fx.run(
        "chr1\t10\t14\tx\n",
        &["chr1\t11\t13\nchr1\t12\t20\n"],
        &["--per-base"],
    );
    assert_eq!(
        rows,
        [
            "chr1\t10\t14\tx\t1\t0",
            "chr1\t10\t14\tx\t2\t1",
            "chr1\t10\t14\tx\t3\t2",
            "chr1\t10\t14\tx\t4\t1",
        ]
    );
    // the point has no rows
    let rows = fx.run(A, &[B], &["--per-base"]);
    assert_eq!(rows.len(), 100 + 100 + 10);
    assert!(rows.iter().all(|r| !r.contains("\tpt\t")));
    assert_eq!(rows[50], "chr1\t100\t200\ta1\t0\t+\t51\t1");
}

#[test]
fn mean() {
    let fx = Fx::new();
    assert_eq!(
        fx.run(A, &[B], &["--mean"]),
        [
            format!("chr1\t100\t200\ta1\t0\t+\t{}", f(30, 100)),
            format!("chr1\t300\t400\ta2\t0\t-\t{}", f(5, 100)),
            "chr2\t5\t15\ta3\t0\t+\t0".to_string(),
            "chr3\t10\t10\tpt\t0\t+\t.".to_string(),
        ]
    );
    // overlapping B add up: depth 2 over 2 bases of a 4-base A
    assert_eq!(
        fx.run(
            "chr1\t10\t14\n",
            &["chr1\t11\t13\nchr1\t12\t20\n"],
            &["--mean"]
        ),
        ["chr1\t10\t14\t1"]
    );
}

#[test]
fn annotate_multi_b() {
    let fx = Fx::new();
    let b2 = "chr1\t0\t120\nchr2\t0\t8\n";
    assert_eq!(
        fx.run(A, &[B, b2], &[]),
        [
            format!("chr1\t100\t200\ta1\t0\t+\t{}\t{}", f(30, 100), f(20, 100)),
            format!("chr1\t300\t400\ta2\t0\t-\t{}\t0", f(5, 100)),
            format!("chr2\t5\t15\ta3\t0\t+\t0\t{}", f(3, 10)),
            "chr3\t10\t10\tpt\t0\t+\t.\t.".to_string(),
        ]
    );
}

#[test]
fn points() {
    let fx = Fx::new();
    // B on the edge of the point A: no overlap unless --inclusive-points
    let a = "chr1\t10\t10\tp\n";
    let b = "chr1\t5\t15\nchr1\t10\t10\nchr1\t10\t20\n";
    assert_eq!(fx.run(a, &[b], &[]), ["chr1\t10\t10\tp\t2\t0\t0\t."]);
    assert_eq!(
        fx.run(a, &[b], &["--inclusive-points"]),
        ["chr1\t10\t10\tp\t3\t0\t0\t."]
    );
    assert_eq!(fx.run(a, &[b], &["--mean"]), ["chr1\t10\t10\tp\t."]);
    assert!(fx.run(a, &[b], &["--hist"]).is_empty());
    assert!(fx.run(a, &[b], &["--per-base"]).is_empty());
    // a point B counts as a hit and covers nothing
    assert_eq!(
        fx.run("chr1\t0\t10\n", &["chr1\t5\t5\n"], &[]),
        ["chr1\t0\t10\t1\t0\t10\t0"]
    );
}

#[test]
fn thresholds_and_strand() {
    let fx = Fx::new();
    // b2 covers 20% of a1, b1 10%
    assert_eq!(
        fx.run(A, &[B], &["-f", "0.2"])[0],
        format!("chr1\t100\t200\ta1\t0\t+\t1\t20\t100\t{}", f(20, 100))
    );
    // b1 is + like a1, b2 is -
    assert_eq!(
        fx.run(A, &[B], &["-s"])[0],
        format!("chr1\t100\t200\ta1\t0\t+\t1\t10\t100\t{}", f(10, 100))
    );
    assert_eq!(
        fx.run(A, &[B], &["-S"])[0],
        format!("chr1\t100\t200\ta1\t0\t+\t1\t20\t100\t{}", f(20, 100))
    );
    assert_eq!(
        fx.run(A, &[B], &["--min-bases", "15"])[0],
        format!("chr1\t100\t200\ta1\t0\t+\t1\t20\t100\t{}", f(20, 100))
    );
}

#[test]
fn split_blocks() {
    let fx = Fx::new();
    // A blocks 100-120 and 180-200; B blocks 100-120 and 160-180: one
    // shared block of 20 bases out of A's 40 block bases
    let a = "chr1\t100\t200\tsp\t0\t+\t100\t200\t0\t2\t20,20\t0,80\n";
    let b = "chr1\t100\t180\tsb\t0\t+\t100\t180\t0\t2\t20,20\t0,60\n";
    let row = |t: &str| format!("{}{t}", a.trim_end());
    assert_eq!(fx.run(a, &[b], &["--split"]), [row("\t1\t20\t40\t0.5")]);
    // without --split the spans cover 80 of A's 100 bases
    assert_eq!(fx.run(a, &[b], &[]), [row("\t1\t80\t100\t0.8")]);
    // a plain B has no blocks: it covers the part of A's blocks it meets
    assert_eq!(
        fx.run(a, &["chr1\t150\t190\n"], &["--split"]),
        [row("\t1\t10\t40\t0.25")]
    );
    // per-base offsets run along the blocks
    let rows = fx.run(a, &["chr1\t110\t190\n"], &["--split", "--per-base"]);
    assert_eq!(rows.len(), 40);
    assert_eq!(rows[9], row("\t10\t0"));
    assert_eq!(rows[10], row("\t11\t1"));
}

#[test]
fn split_hist_and_mean() {
    let fx = Fx::new();
    let a = "chr1\t100\t200\tsp\t0\t+\t100\t200\t0\t2\t20,20\t0,80\n";
    let row = |t: &str| format!("{}{t}", a.trim_end());
    // B1 blocks 100-120 and 160-180, B2 (plain) 110-190: over A's blocks
    // depth is 1 on 100-110, 2 on 110-120, 1 on 180-190, 0 on 190-200
    let b = "chr1\t100\t180\tsb\t0\t+\t100\t180\t0\t2\t20,20\t0,60\nchr1\t110\t190\tp\n";
    assert_eq!(
        fx.run(a, &[b], &["--split", "--hist"]),
        [
            row(&format!("\t0\t10\t40\t{}", f(10, 40))),
            row(&format!("\t1\t20\t40\t{}", f(20, 40))),
            row(&format!("\t2\t10\t40\t{}", f(10, 40))),
            format!("all\t0\t10\t40\t{}", f(10, 40)),
            format!("all\t1\t20\t40\t{}", f(20, 40)),
            format!("all\t2\t10\t40\t{}", f(10, 40)),
        ]
    );
    assert_eq!(fx.run(a, &[b], &["--split", "--mean"]), [row("\t1")]);
    // without --split B1 and B2 are plain spans over A's whole span
    assert_eq!(
        fx.run(a, &[b], &["--mean"]),
        [row(&format!("\t{}", f(80 + 80, 100)))]
    );
}

#[test]
fn vcf_a() {
    let fx = Fx::new();
    let vcf = fx.file(
        "a.vcf",
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t151\t.\tA\tC\t.\t.\t.\nchr1\t161\t.\tA\tC\t.\t.\t.\n",
    );
    let b = fx.file("b.bed", &B[..B.find("chr3").unwrap()]);
    let o = gloci()
        .args(["coverage", "-a", &vcf, "-b", &b])
        .assert()
        .success();
    let out = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    let lines: Vec<&str> = out.lines().collect();
    assert_eq!(
        lines,
        [
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\thits\tcovered_bases\ta_length\tcovered_fraction",
            "chr1\t151\t.\tA\tC\t.\t.\t.\t1\t1\t1\t1",
            "chr1\t161\t.\tA\tC\t.\t.\t.\t0\t0\t1\t0",
        ]
    );
}

#[test]
fn unsorted_matches_sorted() {
    let fx = Fx::new();
    let shuffled = "chr3\t5\t15\tb4\t1\t+\nchr1\t395\t500\tb3\t1\t+\nchr1\t180\t260\tb2\t4\t-\nchr1\t150\t160\tb1\t2\t+\n";
    for extra in [&[][..], &["--hist"], &["--mean"], &["--per-base"]] {
        let sorted = fx.run(A, &[B], extra);
        let mut v = extra.to_vec();
        v.push("--unsorted");
        assert_eq!(fx.run(A, &[shuffled], &v), sorted, "{extra:?}");
    }
}

#[test]
fn usage_and_data_errors() {
    let fx = Fx::new();
    let (a, b) = (fx.file("a.bed", A), fx.file("b.bed", B));
    // several -b files take only the default report
    for flag in ["--hist", "--per-base", "--mean"] {
        assert_eq!(
            code(&["coverage", "-a", &a, "-b", &b, &b, flag]),
            2,
            "{flag}"
        );
    }
    gloci()
        .args(["coverage", "-a", &a, "-b", &b, &b, "--hist"])
        .assert()
        .stderr(contains("one covered fraction per file"));
    // the reports exclude each other, and there are no window flags
    assert_eq!(
        code(&["coverage", "-a", &a, "-b", &b, "--hist", "--mean"]),
        2
    );
    assert_eq!(code(&["coverage", "-a", &a, "-b", &b, "-w", "10"]), 2);
    // an unsorted A is an order error, not silent
    let bad = fx.file("bad.bed", "chr1\t50\t60\nchr1\t10\t20\n");
    assert_eq!(code(&["coverage", "-a", &bad, "-b", &b]), 1);
}

/// A deterministic sorted BED over chr1-chr3: `n` records per chromosome,
/// spaced `step` apart with lengths from `len`..`len + spread`.
fn generated_bed(n: u64, step: u64, len: u64, spread: u64, seed: u64) -> String {
    let mut x = seed;
    let mut next = move || {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        x >> 33
    };
    let mut out = String::new();
    for chrom in ["chr1", "chr2", "chr3"] {
        for i in 0..n {
            let start = i * step + next() % step;
            let end = start + len + next() % spread;
            let strand = ["+", "-", "."][(next() % 3) as usize];
            out += &format!("{chrom}\t{start}\t{end}\tr{i}\t0\t{strand}\n");
        }
    }
    out
}

#[test]
fn determinism() {
    let dir = tempdir();
    // thousands of BGZF records, so -t 8 runs the threaded read and write
    // paths; A is short enough for --per-base
    let a = s(&bgzf_file(
        dir.path(),
        "a.bed.gz",
        &generated_bed(1500, 400, 100, 100, 1),
    ));
    let b = s(&bgzf_file(
        dir.path(),
        "b.bed.gz",
        &generated_bed(3000, 200, 30, 300, 2),
    ));
    for extra in [
        &[][..],
        &["--hist"],
        &["--per-base"],
        &["--mean", "--unsorted"],
        &["-f", "0.1", "-s"],
    ] {
        let run = |t: &str| {
            let mut v = vec!["coverage", "-a", &a, "-b", &b, "-t", t];
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
        assert!(one.len() > 100_000, "{extra:?}: {} bytes", one.len());
        assert_eq!(one, run("2"), "{extra:?}");
        assert_eq!(one, run("8"), "{extra:?}");
    }
    // BGZF output, compressed on several threads
    let out = |t: &str| {
        let path = dir.path().join(format!("out{t}.tsv.gz"));
        gloci()
            .args(["coverage", "-a", &a, "-b", &b, "--hist", "-t", t, "-o"])
            .arg(&path)
            .assert()
            .success();
        std::fs::read(path).unwrap()
    };
    let one = out("1");
    assert_eq!(one, out("8"));
}

#[test]
fn write_error_names_the_destination() {
    let fx = Fx::new();
    let (a, b) = (fx.file("a.bed", A), fx.file("b.bed", B));
    gloci()
        .args(["coverage", "-a", &a, "-b", &b, "-o", "/dev/full"])
        .assert()
        .failure()
        .stderr(contains("error writing /dev/full"));
}
