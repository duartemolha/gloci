//! `--sort` on the relation commands (D39): the result equals `gloci sort`
//! on A and on each B, followed by the command, byte for byte.
mod common;
use common::*;
use predicates::str::contains;
use std::path::Path;

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

/// Deterministic pseudo-random BED6 lines over a few chromosomes.
fn shuffled_bed(seed: u64, n: usize) -> String {
    let mut x = seed | 1;
    let mut next = move || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    };
    let mut out = String::new();
    for i in 0..n {
        let c = ["chr1", "chr2", "chr10", "chrX", "chrM"][(next() % 5) as usize];
        let start = next() % 20_000;
        let len = next() % 300 + 1;
        out += &format!(
            "{c}\t{start}\t{}\tr{i}\t{}\t{}\n",
            start + len,
            next() % 100,
            ["+", "-", "."][(next() % 3) as usize]
        );
    }
    out
}

fn run(args: &[&str]) -> Vec<u8> {
    gloci()
        .args(args)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone()
}

/// Shuffled A and two shuffled B files.
fn inputs(dir: &Path) -> (String, Vec<String>) {
    (
        s(&plain_file(dir, "a.bed", &shuffled_bed(1, 500))),
        vec![
            s(&plain_file(dir, "b1.bed", &shuffled_bed(2, 400))),
            s(&plain_file(dir, "b2.bed", &shuffled_bed(3, 300))),
        ],
    )
}

fn sorted_copy(dir: &Path, name: &str, input: &str) -> String {
    let p = plain_file(dir, name, "");
    std::fs::write(&p, run(&["sort", input])).unwrap();
    s(&p)
}

/// `gloci CMD --sort -a A -b B... EXTRA` against `gloci sort` on each input
/// and the plain command on the sorted files.
fn check(cmd: &str, extra: &[&str]) {
    check_n(cmd, extra, 2);
}

fn check_n(cmd: &str, extra: &[&str], n_b: usize) {
    let dir = tempdir();
    let (a, mut bs) = inputs(dir.path());
    bs.truncate(n_b);
    let sa = sorted_copy(dir.path(), "a.s.bed", &a);
    let sbs: Vec<String> = bs
        .iter()
        .enumerate()
        .map(|(i, b)| sorted_copy(dir.path(), &format!("b{i}.s.bed"), b))
        .collect();
    let build = |a: &str, bs: &[String], flag: Option<&str>| {
        let mut args: Vec<String> = vec![cmd.into()];
        args.extend(flag.map(String::from));
        args.extend(["-a".into(), a.into(), "-b".into()]);
        args.extend(bs.iter().cloned());
        args.extend(extra.iter().map(|x| x.to_string()));
        args
    };
    let run_args = |args: Vec<String>| run(&args.iter().map(String::as_str).collect::<Vec<_>>());
    let want = run_args(build(&sa, &sbs, None));
    assert!(!want.is_empty(), "{cmd} {extra:?} produced nothing");
    let got = run_args(build(&a, &bs, Some("--sort")));
    assert_eq!(
        String::from_utf8_lossy(&got),
        String::from_utf8_lossy(&want),
        "{cmd} {extra:?}"
    );
}

#[test]
fn sort_equals_presorted_pipeline() {
    check("intersect", &[]);
    check("intersect", &["--mode", "pairs", "--names", "X,Y"]);
    check(
        "intersect",
        &[
            "--mode",
            "pairs",
            "--keep-unmatched",
            "-s",
            "-f",
            "50%",
            "--names",
            "X,Y",
        ],
    );
    check("intersect", &["--mode", "count"]);
    check("subtract", &[]);
    check("closest", &["--distance", "-k", "2", "--names", "X,Y"]);
    check(
        "closest",
        &["--signed", "ref", "--per-file", "--names", "X,Y"],
    );
    check("map", &["-c", "5:sum", "-c", "0:count"]);
    check("coverage", &[]);
    check_n("coverage", &["--hist"], 1);
    check_n("coverage", &["--mean", "-s"], 1);
}

#[test]
fn sort_gives_sorted_output_unlike_unsorted() {
    let dir = tempdir();
    let (a, bs) = inputs(dir.path());
    let sorted = run(&["intersect", "--sort", "-a", &a, "-b", &bs[0]]);
    let plain = run(&["intersect", "--unsorted", "-a", &a, "-b", &bs[0]]);
    assert_ne!(sorted, plain, "--unsorted keeps A's order");
    let sa = sorted_copy(dir.path(), "a.s.bed", &a);
    let sb = sorted_copy(dir.path(), "b.s.bed", &bs[0]);
    assert_eq!(sorted, run(&["intersect", "-a", &sa, "-b", &sb]));
}

#[test]
fn sort_on_sorted_input_changes_nothing() {
    let dir = tempdir();
    let (a, bs) = inputs(dir.path());
    let sa = sorted_copy(dir.path(), "a.s.bed", &a);
    let sb = sorted_copy(dir.path(), "b.s.bed", &bs[0]);
    assert_eq!(
        run(&["map", "--sort", "-a", &sa, "-b", &sb, "-c", "5:max"]),
        run(&["map", "-a", &sa, "-b", &sb, "-c", "5:max"])
    );
}

const VCF_HEAD: &str = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n\
##contig=<ID=chr2,length=100000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

#[test]
fn sort_with_vcf_a_uses_the_users_span() {
    let dir = tempdir();
    let mut vcf = VCF_HEAD.to_string();
    for (c, pos, r, a) in [
        ("chr2", 120, "ATTTT", "A"),
        ("chr1", 300, "A", "G"),
        ("chr1", 98, "ATTTT", "A"),
        ("chr2", 40, "A", "ATT"),
        ("chr1", 101, "A", "G"),
    ] {
        vcf += &format!("{c}\t{pos}\t.\t{r}\t{a}\t.\t.\t.\n");
    }
    let a = s(&plain_file(dir.path(), "a.vcf", &vcf));
    let b = s(&plain_file(
        dir.path(),
        "b.bed",
        "chr2\t30\t130\nchr1\t99\t110\nchr1\t290\t310\n",
    ));
    let sa = s(&dir.path().join("a.s.vcf"));
    gloci().args(["sort", &a, "-o", &sa]).assert().success();
    let sb = sorted_copy(dir.path(), "b.s.bed", &b);
    for span in ["allele", "ref"] {
        for mode in ["any", "none"] {
            let want = run(&[
                "intersect",
                "-a",
                &sa,
                "-b",
                &sb,
                "--mode",
                mode,
                "--vcf-span",
                span,
            ]);
            let got = run(&[
                "intersect",
                "--sort",
                "-a",
                &a,
                "-b",
                &b,
                "--mode",
                mode,
                "--vcf-span",
                span,
            ]);
            // The ##glociCommand header line names the input files.
            let body = |v: Vec<u8>| {
                String::from_utf8(v)
                    .unwrap()
                    .lines()
                    .filter(|l| !l.starts_with("##glociCommand"))
                    .collect::<Vec<_>>()
                    .join("\n")
            };
            let (got, want) = (body(got), body(want));
            assert!(
                mode == "none" || want.lines().any(|l| !l.starts_with('#')),
                "{mode} {span}: {want}"
            );
            assert_eq!(got, want, "{mode} {span}");
        }
    }
}

#[test]
fn sort_and_unsorted_conflict() {
    let dir = tempdir();
    let (a, bs) = inputs(dir.path());
    for cmd in ["intersect", "subtract", "closest", "map", "coverage"] {
        let mut args = vec![cmd, "--sort", "--unsorted", "-a", &a, "-b", &bs[0]];
        if cmd == "map" {
            args.extend(["-c", "5:sum"]);
        }
        gloci()
            .args(&args)
            .assert()
            .code(2)
            .stderr(contains("cannot be used with"));
    }
}

#[test]
fn max_mem_and_tmp_dir_need_sort() {
    let dir = tempdir();
    let (a, bs) = inputs(dir.path());
    gloci()
        .args(["intersect", "-a", &a, "-b", &bs[0], "--max-mem", "32M"])
        .assert()
        .code(2);
}

/// 330k shuffled BED lines (about 9 MB) as A: with `--max-mem 16M` the
/// sorter spills at least once. B is small, so only A's sort needs to spill.
fn spill_bed() -> String {
    let mut s = String::with_capacity(10 << 20);
    let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut next = move || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    };
    for _ in 0..330_000 {
        let c = next() % 24 + 1;
        let start = next() % 200_000_000;
        s.push_str(&format!("chr{c}\t{start}\t{}\n", start + next() % 5000 + 1));
    }
    s
}

#[test]
fn sort_spills_within_max_mem() {
    let dir = tempdir();
    let big = s(&plain_file(dir.path(), "big.bed", &spill_bed()));
    let b = s(&plain_file(dir.path(), "b.bed", &shuffled_bed(9, 200)));
    let tmp = dir.path().join("spill");
    std::fs::create_dir(&tmp).unwrap();
    let sa = sorted_copy(dir.path(), "a.s.bed", &big);
    let sb = sorted_copy(dir.path(), "b.s.bed", &b);
    let small = gloci()
        .args([
            "intersect",
            "--sort",
            "--max-mem",
            "16M",
            "-v",
            "-a",
            &big,
            "-b",
            &b,
            "--tmp-dir",
        ])
        .arg(&tmp)
        .assert()
        .success()
        .get_output()
        .clone();
    let log = String::from_utf8_lossy(&small.stderr);
    assert!(
        log.lines()
            .any(|l| l.contains("sorted 330000 records") && !l.contains(" 0 chunks spilled")),
        "A's sort should have spilled: {log}"
    );
    let want = run(&["intersect", "-a", &sa, "-b", &sb]);
    assert!(!want.is_empty());
    assert_eq!(small.stdout, want);
    assert_eq!(std::fs::read_dir(&tmp).unwrap().count(), 0);
}

#[test]
fn order_error_hint_names_both_flags() {
    let dir = tempdir();
    let a = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t50\t60\nchr1\t10\t20\n",
    ));
    let b = s(&plain_file(dir.path(), "b.bed", "chr1\t0\t100\n"));
    for cmd in ["intersect", "subtract", "closest", "map", "coverage"] {
        let mut args = vec![cmd, "-a", &a, "-b", &b];
        if cmd == "map" {
            args.extend(["-c", "0:count"]);
        }
        gloci()
            .args(&args)
            .assert()
            .code(1)
            .stderr(contains("is not sorted"))
            .stderr(contains("--sort"))
            .stderr(contains("--unsorted"));
    }
}

#[test]
fn sort_determinism() {
    let dir = tempdir();
    let (a, bs) = inputs(dir.path());
    let run_t = |t: &str| {
        run(&[
            "intersect",
            "--sort",
            "--mode",
            "pairs",
            "-s",
            "-a",
            &a,
            "-b",
            &bs[0],
            &bs[1],
            "--names",
            "X,Y",
            "-t",
            t,
        ])
    };
    let base = run_t("1");
    assert!(!base.is_empty());
    assert_eq!(run_t("2"), base);
    assert_eq!(run_t("8"), base);
}

#[test]
fn default_labels_under_sort_are_the_users_files() {
    let dir = tempdir();
    let (a, bs) = inputs(dir.path());
    let out = run(&[
        "intersect",
        "--sort",
        "--mode",
        "pairs",
        "-a",
        &a,
        "-b",
        &bs[0],
        &bs[1],
    ]);
    let text = String::from_utf8(out).unwrap();
    assert!(!text.contains("sorted inputs"), "{text}");
    let labels: std::collections::BTreeSet<&str> = text
        .lines()
        .map(|l| l.split('\t').nth(6).unwrap())
        .collect();
    // The same labels as the unsorted run gives, and not temp-file names.
    let plain = run(&[
        "intersect",
        "--unsorted",
        "--mode",
        "pairs",
        "-a",
        &a,
        "-b",
        &bs[0],
        &bs[1],
    ]);
    let plain = String::from_utf8(plain).unwrap();
    let want: std::collections::BTreeSet<&str> = plain
        .lines()
        .map(|l| l.split('\t').nth(6).unwrap())
        .collect();
    assert_eq!(want.len(), 2);
    assert_eq!(labels, want);
}

#[test]
fn failure_after_the_sorts_exits_1_and_cleans_up() {
    let dir = tempdir();
    // Overlapping BED12 blocks: rejected mid-relation, with the user's file name.
    let bad = s(&plain_file(
        dir.path(),
        "bad.bed",
        "chr1\t0\t100\tn\t0\t+\t0\t100\t0\t2\t50,50\t0,40\n",
    ));
    let b = s(&plain_file(dir.path(), "b.bed", "chr1\t0\t200\n"));
    let tmp = dir.path().join("spill");
    std::fs::create_dir(&tmp).unwrap();
    gloci()
        .args([
            "intersect",
            "--sort",
            "--split",
            "-a",
            &bad,
            "-b",
            &b,
            "--tmp-dir",
        ])
        .arg(&tmp)
        .assert()
        .code(1)
        .stderr(contains("bad.bed"))
        .stderr(contains("invalid BED12 blocks"));
    assert_eq!(std::fs::read_dir(&tmp).unwrap().count(), 0);
}

#[test]
fn tmp_dir_needs_sort() {
    let dir = tempdir();
    let (a, bs) = inputs(dir.path());
    gloci()
        .args(["intersect", "-a", &a, "-b", &bs[0], "--tmp-dir"])
        .arg(dir.path())
        .assert()
        .code(2);
}

#[test]
fn usage_errors_come_before_any_sorting() {
    let dir = tempdir();
    let (a, bs) = inputs(dir.path());
    let no_sort = |args: &[&str]| {
        let out = gloci()
            .args(args)
            .args(["-v", "--sort", "-a", &a, "-b", &bs[0]])
            .assert()
            .code(2)
            .get_output()
            .clone();
        let err = String::from_utf8_lossy(&out.stderr).into_owned();
        assert!(
            !err.contains("--sort: sorted"),
            "{args:?} sorted first: {err}"
        );
    };
    no_sort(&["intersect", "--overlap-bases"]);
    no_sort(&["map", "-c", "5:nonsense"]);
    no_sort(&["coverage", "--hist", "--per-base"]);
    no_sort(&["subtract", "--sum-overlaps", "--whole"]);
    no_sort(&["closest", "-k", "0"]);
}

#[test]
fn warnings_under_sort_name_the_users_files() {
    let dir = tempdir();
    let (a, bs) = inputs(dir.path());
    let out = gloci()
        .args(["intersect", "--sort", "-s", "-a", &a, "-b", &bs[0]])
        .assert()
        .success()
        .get_output()
        .clone();
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        err.contains(&format!("{a} has records with unknown strand")),
        "{err}"
    );
    assert!(!err.contains("sorted inputs"), "{err}");
}
