//! `gloci stats fisher`: the report text, merging, `--exclude`, strands,
//! empty inputs, `-g` requirements, order errors and determinism across
//! thread counts. The bedtools parity is in the golden cases
//! (`stats-fisher-*-vs-bedtools`).
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
    /// `stats fisher -a A -b B -g genome extra...` over BED texts, genome `chr1 100000`: stdout.
    fn run(&self, a: &str, b: &str, extra: &[&str]) -> String {
        let af = self.file("a.bed", a);
        let bf = self.file("b.bed", b);
        let g = self.file("g.txt", "chr1\t100000\n");
        let o = gloci()
            .args(["stats", "fisher", "-a", &af, "-b", &bf, "-g", &g])
            .args(extra)
            .assert()
            .success();
        String::from_utf8(o.get_output().stdout.clone()).unwrap()
    }
}

const A: &str = "chr1\t100\t200\nchr1\t300\t400\nchr1\t500\t600\n";
const B: &str = "chr1\t150\t250\nchr1\t550\t700\nchr1\t800\t900\n";

#[test]
fn fisher_report_text() {
    let got = Fx::new().run(A, B, &[]);
    let want = "\
# Number of query intervals: 3
# Number of db intervals: 3
# Number of overlaps: 2
# Number of possible intervals (estimated): 457
# phyper(2 - 1, 3, 457 - 3, 3, lower.tail=F)
# Contingency Table Of Counts
#_________________________________________
#           |  in -b       | not in -b    |
#     in -a | 2            | 1            |
# not in -a | 1            | 453          |
#_________________________________________
# p-values for fisher's exact test
left\tright\ttwo-tail\tratio
1\t8.6249e-05\t8.6249e-05\t906.000
";
    assert_eq!(got, want);
}

fn last_line(r: &str) -> &str {
    r.lines().last().unwrap()
}

fn estimate(r: &str) -> &str {
    r.lines().nth(3).unwrap()
}

#[test]
fn fisher_merge_overlaps_counts_merged_intervals() {
    let fx = Fx::new();
    // 10-20 and 15-30 merge: unmerged A has two records, merged A one.
    let a = "chr1\t10\t20\nchr1\t15\t30\n";
    let b = "chr1\t12\t14\n";
    assert!(
        fx.run(a, b, &[])
            .contains("# Number of query intervals: 2\n")
    );
    let m = fx.run(a, b, &["-m"]);
    assert!(m.contains("# Number of query intervals: 1\n"), "{m}");
    assert_eq!(m, fx.run(a, b, &["--merge-overlaps"]));
}

#[test]
fn fisher_exact_edge_text() {
    let fx = Fx::new();
    // All pairs overlap: n12 = n21 = 0, so the ratio is inf, and n11 > nA clamps n12 to 0.
    let r = fx.run(A, A, &[]);
    assert!(last_line(&r).ends_with("\tinf"), "{r}");
    // An empty A: ratio -nan, p-values 1, and the n22 estimate stays at n11 + n12 + n21.
    let r = fx.run("", B, &[]);
    assert_eq!(last_line(&r), "1\t1\t1\t-nan");
    assert_eq!(
        estimate(&r),
        "# Number of possible intervals (estimated): 3"
    );
    // p-values below 1e-5 print in exponent form.
    let r = fx.run(A, A, &[]);
    assert!(last_line(&r).contains("e-"), "{r}");
}

#[test]
fn fisher_exclude_shrinks_the_genome_once_per_base() {
    let fx = Fx::new();
    let base = fx.run(A, B, &[]);
    assert_eq!(
        estimate(&base),
        "# Number of possible intervals (estimated): 457"
    );
    // Merged exclude 0-50000 (three overlapping records) leaves 50000 bases: 228 intervals.
    let ex = fx.file(
        "ex.bed",
        "chr1\t0\t30000\nchr1\t10000\t40000\nchr1\t20000\t50000\n",
    );
    let r = fx.run(A, B, &["--exclude", &ex]);
    assert_eq!(
        estimate(&r),
        "# Number of possible intervals (estimated): 228"
    );
    // An unsorted exclude file is refused, and --sort sorts it.
    let bad = fx.file("bad.bed", "chr1\t1000\t2000\nchr1\t0\t1500\n");
    let af = fx.file("a2.bed", A);
    let bf = fx.file("b2.bed", B);
    let g = fx.file("g2.txt", "chr1\t100000\n");
    gloci()
        .args([
            "stats",
            "fisher",
            "-a",
            &af,
            "-b",
            &bf,
            "-g",
            &g,
            "--exclude",
            &bad,
        ])
        .assert()
        .code(1)
        .stderr(contains("not sorted"));
    let r = fx.run(A, B, &["--exclude", &bad, "--sort"]);
    assert_eq!(
        estimate(&r),
        "# Number of possible intervals (estimated): 448"
    );
}

#[test]
fn fisher_strand_rules() {
    let fx = Fx::new();
    let a = "chr1\t100\t200\t.\t0\t+\nchr1\t300\t400\t.\t0\t.\nchr1\t500\t600\t.\t0\t-\n";
    let b = "chr1\t150\t250\t.\t0\t+\nchr1\t350\t450\t.\t0\t.\nchr1\t550\t700\t.\t0\t+\n";
    // Unmerged: unknown-strand records still count in nA, nB, but never match.
    let r = fx.run(a, b, &["-s"]);
    assert!(r.contains("# Number of query intervals: 3\n"), "{r}");
    assert!(r.contains("# Number of overlaps: 1\n"), "{r}");
    // Merged per strand: unknown-strand clusters take no part.
    let r = fx.run(a, b, &["-s", "-m"]);
    assert!(r.contains("# Number of query intervals: 2\n"), "{r}");
    assert!(r.contains("# Number of db intervals: 2\n"), "{r}");
    // -S: opposite strands only (A - 500-600 and B + 550-700).
    let r = fx.run(a, b, &["-S"]);
    assert!(r.contains("# Number of overlaps: 1\n"), "{r}");
}

#[test]
fn fisher_needs_genome_lengths() {
    let fx = Fx::new();
    let af = fx.file("a.bed", A);
    let bf = fx.file("b.bed", B);
    gloci()
        .args(["stats", "fisher", "-a", &af, "-b", &bf])
        .assert()
        .code(2);
    // A genome without lengths (a VCF header with bare contigs) cannot size the genome.
    let vcf = fx.file(
        "h.vcf",
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    gloci()
        .args(["stats", "fisher", "-a", &af, "-b", &bf, "-g", &vcf])
        .assert()
        .code(2)
        .stderr(contains("no length for chr1"))
        .stderr(contains("-g"));
}

#[test]
fn fisher_unsorted_input_is_refused_and_sort_fixes_it() {
    let fx = Fx::new();
    let a = "chr1\t500\t600\nchr1\t100\t200\n";
    let af = fx.file("a.bed", a);
    let bf = fx.file("b.bed", B);
    let g = fx.file("g.txt", "chr1\t100000\n");
    gloci()
        .args(["stats", "fisher", "-a", &af, "-b", &bf, "-g", &g])
        .assert()
        .code(1)
        .stderr(contains("not sorted"))
        .stderr(contains("--sort"));
    let sorted = fx.run("chr1\t100\t200\nchr1\t500\t600\n", B, &[]);
    assert_eq!(fx.run(a, B, &["--sort"]), sorted);
}

#[test]
fn fisher_is_identical_at_every_thread_count() {
    let fx = Fx::new();
    let mut a = String::new();
    let mut b = String::new();
    for i in 0..600u64 {
        a += &format!(
            "chr1\t{}\t{}\t.\t0\t{}\n",
            i * 37,
            i * 37 + 50,
            ['+', '-'][(i % 2) as usize]
        );
        b += &format!(
            "chr1\t{}\t{}\t.\t0\t{}\n",
            i * 41,
            i * 41 + 30,
            ['+', '-'][((i / 3) % 2) as usize]
        );
    }
    let af = fx.file("a.bed", &a);
    let bf = fx.file("b.bed", &b);
    let g = fx.file("g.txt", "chr1\t1000000\n");
    for extra in [&[][..], &["-s"][..], &["-m"][..], &["-f", "0.3", "-r"][..]] {
        let outs: Vec<String> = ["1", "2", "8"]
            .iter()
            .map(|t| {
                let o = gloci()
                    .args(["stats", "fisher", "-a", &af, "-b", &bf, "-g", &g, "-t", t])
                    .args(extra)
                    .assert()
                    .success();
                String::from_utf8(o.get_output().stdout.clone()).unwrap()
            })
            .collect();
        assert_eq!(outs[0], outs[1]);
        assert_eq!(outs[0], outs[2]);
        assert!(outs[0].starts_with("# Number of query intervals: "));
    }
}

#[test]
fn fisher_empty_a_and_b_estimate() {
    let fx = Fx::new();
    // bedtools' value is undefined behaviour (differs by CPU); gloci gives n11 + n12 + n21.
    let r = fx.run("", B, &[]);
    assert_eq!(
        estimate(&r),
        "# Number of possible intervals (estimated): 3"
    );
    assert!(
        r.contains("# not in -a | 3            | 0            |\n"),
        "{r}"
    );
    assert_eq!(last_line(&r), "1\t1\t1\t-nan");
    let r = fx.run(A, "", &[]);
    assert_eq!(
        estimate(&r),
        "# Number of possible intervals (estimated): 3"
    );
    assert!(
        r.contains("#     in -a | 0            | 3            |\n"),
        "{r}"
    );
    assert_eq!(last_line(&r), "1\t1\t1\t-nan");
}
