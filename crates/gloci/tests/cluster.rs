//! `gloci cluster` and `gloci spacing`: gap table, strands, global ids,
//! points, VCF footprints, a property test against a naive version and
//! against `merge`, order errors and determinism.
mod common;
use common::*;
use predicates::str::contains;
use std::path::Path;

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn out(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

fn bed(rows: &[(u64, u64)]) -> String {
    rows.iter()
        .map(|(a, b)| format!("chr1\t{a}\t{b}\n"))
        .collect()
}

fn ids(text: &str) -> Vec<String> {
    text.lines()
        .map(|l| l.rsplit('\t').next().unwrap().to_string())
        .collect()
}

#[test]
fn gap_table() {
    let dir = tempdir();
    type Case = (&'static [(u64, u64)], &'static str, &'static [&'static str]);
    let cases: [Case; 8] = [
        (&[(1, 5), (5, 8)], "0", &["1", "1"]),
        (&[(1, 5), (6, 8)], "0", &["1", "2"]),
        (&[(1, 5), (6, 8)], "1", &["1", "1"]),
        (&[(1, 5), (6, 18), (16, 20)], "10", &["1", "1", "1"]),
        (&[(1, 5), (4, 8)], "-2", &["1", "2"]),
        (&[(1, 5), (3, 8)], "-2", &["1", "1"]),
        // the cluster end is the furthest end: (2, 3) is inside (1, 9)
        (&[(1, 9), (2, 3), (8, 12)], "-2", &["1", "1", "2"]),
        // points: a point joins a cluster it lies inside or touches
        (
            &[(1, 5), (3, 3), (5, 5), (9, 9)],
            "0",
            &["1", "1", "1", "2"],
        ),
    ];
    for (i, (rows, d, want)) in cases.iter().enumerate() {
        let f = s(&plain_file(dir.path(), &format!("g{i}.bed"), &bed(rows)));
        let got = out(&["cluster", &f, "-d", d]);
        assert_eq!(ids(&got), *want, "case {i}: {got}");
        assert_eq!(out(&["cluster", &f, &format!("--max-gap={d}")]), got);
    }
}

#[test]
fn input_lines_are_kept_byte_for_byte() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t1\t5\tname with space\t0.50\t+\textra\nchr1\t20\t30\n",
    ));
    assert_eq!(
        out(&["cluster", &f]),
        "chr1\t1\t5\tname with space\t0.50\t+\textra\t1\nchr1\t20\t30\t2\n"
    );
}

#[test]
fn ids_are_global_across_chromosomes() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t1\t5\nchr1\t4\t9\nchr2\t1\t5\nchr2\t100\t110\nchr3\t5\t9\n",
    ));
    assert_eq!(ids(&out(&["cluster", &f])), ["1", "1", "2", "3", "4"]);
}

#[test]
fn strands_get_ids_in_order_of_first_record() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t0\t100\ta\t1\t+\nchr1\t50\t150\tb\t1\t-\nchr1\t60\t70\tc\t1\t-\nchr1\t90\t200\td\t1\t+\n\
         chr1\t300\t310\te\t1\t.\nchr2\t1\t5\tf\t1\t-\n",
    ));
    // input order is kept, a `+` and a `-` cluster interleave
    let o = gloci().args(["cluster", &f, "-s"]).assert().success();
    let got = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    assert_eq!(ids(&got), ["1", "2", "2", "1", "3", "4"]);
    o.stderr(contains("no strand"));
    // without -s the strands share clusters
    assert_eq!(ids(&out(&["cluster", &f])), ["1", "1", "1", "1", "2", "3"]);
}

#[test]
fn vcf_input_uses_footprints_and_names_the_column() {
    let dir = tempdir();
    let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
        chr1\t5\t.\tA\tT\t.\t.\t.\nchr1\t6\t.\tA\tT\t.\t.\t.\n\
        chr1\t99\t.\tA\tT\t.\t.\t.\nchr1\t100\t.\tACG\tA\t.\t.\t.\n";
    let f = s(&plain_file(dir.path(), "a.vcf", vcf));
    // reference span of the deletion is [99, 102): bookended to the SNV at 99
    let o = out(&["cluster", &f, "--vcf-span", "ref"]);
    assert_eq!(
        o,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tcluster\n\
         chr1\t5\t.\tA\tT\t.\t.\t.\t1\nchr1\t6\t.\tA\tT\t.\t.\t.\t1\n\
         chr1\t99\t.\tA\tT\t.\t.\t.\t2\nchr1\t100\t.\tACG\tA\t.\t.\t.\t2\n"
    );
    // the allele footprint of the deletion is [100, 102): one base from the SNV
    let o = out(&["cluster", &f]);
    assert_eq!(ids(&o)[1..], ["1", "1", "2", "3"]);
}

#[test]
fn bcf_is_refused_with_the_fix() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "a.vcf",
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t5\t.\tA\tT\t.\t.\t.\n",
    ));
    let b = s(&dir.path().join("a.bcf"));
    gloci()
        .args(["view", &f, "-O", "bcf", "-o", &b])
        .assert()
        .success();
    for cmd in ["cluster", "spacing"] {
        gloci()
            .args([cmd, &b])
            .assert()
            .code(2)
            .stderr(contains("gloci view FILE -O vcf"));
    }
}

#[test]
fn header_flag_keeps_bed_header_lines() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "a.bed",
        "track name=x\nchr1\t1\t5\nchr1\t7\t9\n",
    ));
    assert_eq!(
        out(&["cluster", &f, "--header"]),
        "track name=x\nchr1\t1\t5\t1\nchr1\t7\t9\t2\n"
    );
    assert_eq!(out(&["cluster", &f]), "chr1\t1\t5\t1\nchr1\t7\t9\t2\n");
}

#[test]
fn unsorted_input_names_the_fix_and_sort_flag_works() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "u.bed",
        "chr1\t50\t60\nchr1\t0\t10\nchr1\t5\t12\n",
    ));
    for cmd in ["cluster", "spacing"] {
        gloci()
            .args([cmd, &f])
            .assert()
            .code(1)
            .stderr(contains("not sorted").and(contains("--sort")));
    }
    assert_eq!(
        out(&["cluster", &f, "--sort"]),
        "chr1\t0\t10\t1\nchr1\t5\t12\t1\nchr1\t50\t60\t2\n"
    );
    assert_eq!(
        out(&["spacing", &f, "--sort"]),
        "chr1\t0\t10\t.\nchr1\t5\t12\t-1\nchr1\t50\t60\t38\n"
    );
}

use predicates::prelude::*;

#[test]
fn spacing_table() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t0\t10\nchr1\t10\t20\nchr1\t25\t30\nchr1\t26\t28\nchr1\t40\t50\n\
         chr1\t0\t0\nchr2\t5\t9\nchr3\t0\t5\nchr3\t0\t5\nchr3\t0\t3\n",
    ));
    // sorting is not requested here: chr1 0 0 is out of order, so use --sort
    let o = out(&["spacing", &f, "--sort"]);
    let got: Vec<&str> = o.lines().map(|l| l.rsplit('\t').next().unwrap()).collect();
    // chr1 sorted: (0,0) (0,10) (10,20) (25,30) (26,28) (40,50)
    assert_eq!(
        got,
        [".", "0", "0", "5", "-1", "12", ".", ".", "-1", "-1"],
        "{o}"
    );
}

#[test]
fn spacing_first_bookended_overlap_and_chromosome_change() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t0\t10\ta\nchr1\t10\t20\tb\nchr1\t15\t18\tc\nchr1\t30\t40\td\nchr2\t0\t5\te\n",
    ));
    assert_eq!(
        out(&["spacing", &f]),
        "chr1\t0\t10\ta\t.\nchr1\t10\t20\tb\t0\nchr1\t15\t18\tc\t-1\nchr1\t30\t40\td\t12\nchr2\t0\t5\te\t.\n"
    );
}

/// A small deterministic generator (no external crates).
struct Lcg(u64);
impl Lcg {
    fn next(&mut self, n: u64) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) % n
    }
}

struct Rec {
    chrom: usize,
    start: u64,
    end: u64,
    strand: usize,
}

fn random_bed(seed: u64) -> (Vec<Rec>, String) {
    let mut g = Lcg(seed);
    let mut recs: Vec<Rec> = (0..60)
        .map(|_| {
            let start = g.next(300);
            Rec {
                chrom: g.next(3) as usize,
                start,
                end: start + g.next(25),
                strand: g.next(2) as usize,
            }
        })
        .collect();
    recs.sort_by_key(|r| (r.chrom, r.start, r.end));
    let text = recs
        .iter()
        .map(|r| {
            format!(
                "chr{}\t{}\t{}\tn\t0\t{}\n",
                r.chrom + 1,
                r.start,
                r.end,
                ['+', '-'][r.strand]
            )
        })
        .collect();
    (recs, text)
}

/// The definition, written out: per chromosome (and strand), a record joins
/// the cluster of the previous one in its group when it starts at most `d`
/// after the cluster's furthest end; ids count clusters by first record.
fn naive(recs: &[Rec], d: i64, per_strand: bool) -> Vec<u64> {
    let mut next = 0;
    let mut ids = Vec::new();
    // (chrom, strand slot) -> (furthest end, id)
    let mut open: std::collections::HashMap<(usize, usize), (u64, u64)> = Default::default();
    let mut chrom = usize::MAX;
    for r in recs {
        if r.chrom != chrom {
            open.clear();
            chrom = r.chrom;
        }
        let key = (r.chrom, if per_strand { r.strand } else { 0 });
        let joined = open
            .get(&key)
            .filter(|(end, _)| i128::from(r.start) <= i128::from(*end) + i128::from(d))
            .map(|&(end, id)| (end.max(r.end), id));
        let (end, id) = joined.unwrap_or_else(|| {
            next += 1;
            (r.end, next)
        });
        open.insert(key, (end, id));
        ids.push(id);
    }
    ids
}

#[test]
fn matches_the_naive_definition_and_merge() {
    let dir = tempdir();
    for seed in 0..25 {
        let (recs, text) = random_bed(seed);
        let f = s(&plain_file(dir.path(), "r.bed", &text));
        for d in [-5i64, -1, 0, 3, 40] {
            for per_strand in [false, true] {
                let mut args = vec!["cluster".to_string(), f.clone(), format!("--max-gap={d}")];
                if per_strand {
                    args.push("-s".into());
                }
                let a: Vec<&str> = args.iter().map(String::as_str).collect();
                let got = out(&a);
                let want = naive(&recs, d, per_strand);
                let got_ids: Vec<u64> = ids(&got).iter().map(|i| i.parse().unwrap()).collect();
                assert_eq!(got_ids, want, "seed {seed} d {d} strand {per_strand}");
                // the lines are the input's, in order
                let lines: Vec<&str> = got
                    .lines()
                    .map(|l| l.rsplit_once('\t').unwrap().0)
                    .collect();
                assert_eq!(lines, text.lines().collect::<Vec<_>>());
                // as many clusters as merge writes rows, with the same sizes
                let mut m = vec!["merge", &f, "--count", "-d"];
                let ds = d.to_string();
                m.push(&ds);
                if per_strand {
                    m.push("-s");
                }
                let merged = out(&m);
                let mut sizes: Vec<u64> = merged
                    .lines()
                    .map(|l| l.rsplit('\t').next().unwrap().parse().unwrap())
                    .collect();
                let mut counts = std::collections::BTreeMap::new();
                for i in &want {
                    *counts.entry(i).or_insert(0u64) += 1;
                }
                let mut from_ids: Vec<u64> = counts.into_values().collect();
                sizes.sort_unstable();
                from_ids.sort_unstable();
                assert_eq!(sizes, from_ids, "seed {seed} d {d} strand {per_strand}");
            }
        }
    }
}

#[test]
fn spacing_matches_the_definition() {
    let dir = tempdir();
    for seed in 0..25 {
        let (recs, text) = random_bed(seed + 100);
        let f = s(&plain_file(dir.path(), "r.bed", &text));
        let got = out(&["spacing", &f]);
        let want: Vec<String> = recs
            .iter()
            .enumerate()
            .map(|(i, r)| match i.checked_sub(1).map(|p| &recs[p]) {
                Some(p) if p.chrom == r.chrom => {
                    if r.start < p.end {
                        "-1".to_string()
                    } else {
                        (r.start - p.end).to_string()
                    }
                }
                _ => ".".to_string(),
            })
            .collect();
        assert_eq!(ids(&got), want, "seed {seed}");
    }
}

#[test]
fn identical_for_every_thread_count_and_through_stdin() {
    let dir = tempdir();
    let (_, text) = random_bed(7);
    let f = s(&plain_file(dir.path(), "r.bed", &text));
    for cmd in ["cluster", "spacing"] {
        let base = out(&[cmd, &f, "-t", "1"]);
        for t in ["2", "4", "8"] {
            assert_eq!(out(&[cmd, &f, "-t", t]), base, "{cmd} -t {t}");
        }
        let o = gloci()
            .args([cmd, "-"])
            .write_stdin(text.clone())
            .assert()
            .success();
        assert_eq!(
            String::from_utf8(o.get_output().stdout.clone()).unwrap(),
            base
        );
    }
}

#[test]
fn help_has_examples() {
    for cmd in ["cluster", "spacing"] {
        gloci()
            .args([cmd, "--help"])
            .assert()
            .success()
            .stdout(contains("Examples:"));
    }
}

/// Groups of record names: `cluster`'s (records sharing an id) and `merge`'s
/// (`-c COL:collapse` over the names), each as sorted lists of sorted names.
fn groups_from_cluster(text: &str, name_col: usize) -> Vec<Vec<String>> {
    let mut by_id: std::collections::BTreeMap<String, Vec<String>> = Default::default();
    for l in text.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = l.split('\t').collect();
        by_id
            .entry(f[f.len() - 1].to_string())
            .or_default()
            .push(f[name_col].to_string());
    }
    let mut g: Vec<Vec<String>> = by_id.into_values().collect();
    g.iter_mut().for_each(|v| v.sort());
    g.sort();
    g
}

fn groups_from_merge(text: &str, collapse_col: usize) -> Vec<Vec<String>> {
    let mut g: Vec<Vec<String>> = text
        .lines()
        .map(|l| {
            let mut v: Vec<String> = l
                .split('\t')
                .nth(collapse_col)
                .unwrap()
                .split(',')
                .map(String::from)
                .collect();
            v.sort();
            v
        })
        .collect();
    g.sort();
    g
}

fn random_vcf(seed: u64) -> String {
    let mut g = Lcg(seed);
    let mut recs: Vec<(u64, usize)> = (0..40)
        .map(|_| (90 + g.next(40), 1 + g.next(6) as usize))
        .collect();
    recs.sort();
    let mut t =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for (i, (pos, len)) in recs.iter().enumerate() {
        let refa = &"ACGTAC"[..*len];
        let alt = if *len > 1 || g.next(2) == 0 { "A" } else { "T" };
        t.push_str(&format!("chr1\t{pos}\tR{i}\t{refa}\t{alt}\t.\t.\t.\n"));
    }
    t
}

#[test]
fn groups_equal_merge_groups_on_random_bed_and_vcf() {
    let dir = tempdir();
    let mut inputs: Vec<(String, usize, Vec<&str>)> = Vec::new();
    for seed in 0..20 {
        let (_, text) = random_bed(seed + 500);
        // unique names in column 4
        let named: String = text
            .lines()
            .enumerate()
            .map(|(i, l)| {
                let mut f: Vec<&str> = l.split('\t').collect();
                let n = format!("n{i}");
                f[3] = &n;
                f.join("\t") + "\n"
            })
            .collect();
        inputs.push((format!("b{seed}.bed"), 3, vec![]));
        plain_file(dir.path(), &format!("b{seed}.bed"), &named);
        plain_file(dir.path(), &format!("v{seed}.vcf"), &random_vcf(seed + 900));
        inputs.push((format!("v{seed}.vcf"), 2, vec![]));
        inputs.push((format!("v{seed}.vcf"), 2, vec!["--vcf-span", "ref"]));
    }
    // the reviewer's repro
    plain_file(
        dir.path(),
        "repro.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         chr1\t95\tR0\tACGTAC\tA\t.\t.\t.\nchr1\t100\tR1\tAC\tA\t.\t.\t.\nchr1\t100\tR2\tA\tT\t.\t.\t.\n",
    );
    inputs.push(("repro.vcf".into(), 2, vec![]));
    for (name, col, extra) in &inputs {
        let f = s(&dir.path().join(name));
        let is_vcf = name.ends_with(".vcf");
        for d in ["-7", "-2", "-1", "0", "1", "5"] {
            for strand in [false, true] {
                if strand && is_vcf {
                    continue;
                }
                let mut c = vec!["cluster", &f, "-d", d];
                let mut m = vec!["merge".to_string(), f.clone(), "-d".into(), d.into()];
                let cc = format!("{}:collapse", col + 1);
                m.extend(["-c".into(), cc]);
                if strand {
                    c.push("-s");
                    m.push("-s".into());
                }
                c.extend(extra);
                m.extend(extra.iter().map(|s| s.to_string()));
                let m: Vec<&str> = m.iter().map(String::as_str).collect();
                let got = groups_from_cluster(&out(&c), *col);
                let want = groups_from_merge(&out(&m), 3);
                assert_eq!(got, want, "{name} d {d} strand {strand} {extra:?}");
            }
        }
    }
    // the repro explicitly: R0 and R2 together, R1 apart (as merge says)
    let f = s(&dir.path().join("repro.vcf"));
    let o = out(&["cluster", &f, "-d", "-1"]);
    assert_eq!(ids(&o)[1..], ["1", "2", "1"]);
}
