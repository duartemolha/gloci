//! `--sort` on `merge` and `complement` (D39): the result equals
//! `gloci sort INPUTS | gloci <cmd> -`, byte for byte, within `--max-mem`.
mod common;
use common::*;
use predicates::str::contains;
use std::path::Path;

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

/// Deterministic pseudo-random BED lines, `n` of them, over a few chromosomes.
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
        let start = next() % 100_000;
        let len = next() % 300;
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

fn files(dir: &Path) -> Vec<String> {
    vec![
        s(&plain_file(dir, "a.bed", &shuffled_bed(1, 400))),
        s(&plain_file(dir, "b.bed", &shuffled_bed(2, 300))),
        s(&plain_file(dir, "c.bed", &shuffled_bed(3, 200))),
    ]
}

/// `gloci sort INPUTS > sorted.bed; gloci CMD ARGS sorted.bed`.
fn via_sort_pipe(dir: &Path, inputs: &[String], cmd: &[&str]) -> Vec<u8> {
    let mut sort_args = vec!["sort"];
    sort_args.extend(inputs.iter().map(String::as_str));
    // The same chromosome order as the command: -g goes to both.
    if let Some(i) = cmd.iter().position(|a| *a == "-g") {
        sort_args.extend(&cmd[i..=i + 1]);
    }
    let sorted = plain_file(dir, "sorted.bed", "");
    std::fs::write(&sorted, run(&sort_args)).unwrap();
    let mut args = cmd.to_vec();
    let sorted = s(&sorted);
    args.insert(1, &sorted);
    run(&args)
}

#[test]
fn merge_sort_equals_sort_pipe() {
    let dir = tempdir();
    let inputs = files(dir.path());
    for extra in [
        vec![],
        vec!["-d", "50"],
        vec!["-s", "--count", "-c", "5:max"],
    ] {
        let mut cmd = vec!["merge"];
        cmd.extend(&extra);
        let want = via_sort_pipe(dir.path(), &inputs, &cmd);
        assert!(!want.is_empty());
        let mut args = vec!["merge", "--sort"];
        args.extend(inputs.iter().map(String::as_str));
        args.extend(&extra);
        assert_eq!(run(&args), want, "merge {extra:?}");
    }
}

#[test]
fn complement_sort_equals_sort_pipe() {
    let dir = tempdir();
    let inputs = files(dir.path());
    let g = s(&plain_file(
        dir.path(),
        "g.genome",
        "chrX\t200000\nchr10\t200000\nchr2\t200000\nchr1\t200000\nchrM\t200000\n",
    ));
    for extra in [vec![], vec!["-g", g.as_str()]] {
        let mut cmd = vec!["complement"];
        cmd.extend(&extra);
        let want = via_sort_pipe(dir.path(), &inputs, &cmd);
        assert!(!want.is_empty());
        let mut args = vec!["complement", "--sort"];
        args.extend(inputs.iter().map(String::as_str));
        args.extend(&extra);
        assert_eq!(run(&args), want, "complement {extra:?}");
    }
}

#[test]
fn sort_flag_on_sorted_input_changes_nothing() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t1\t5\nchr1\t9\t12\n",
    ));
    assert_eq!(run(&["merge", &f, "--sort"]), run(&["merge", &f]));
    assert_eq!(run(&["complement", &f, "--sort"]), run(&["complement", &f]));
}

/// About 50 MB of shuffled BED lines.
fn big_bed() -> String {
    let mut s = String::with_capacity(52 << 20);
    let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut next = move || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    };
    while s.len() < 50 << 20 {
        let c = next() % 24 + 1;
        let start = next() % 200_000_000;
        s.push_str(&format!("chr{c}\t{start}\t{}\n", start + next() % 5000 + 1));
    }
    s
}

#[test]
fn sort_spills_within_max_mem() {
    let dir = tempdir();
    let big = s(&plain_file(dir.path(), "big.bed", &big_bed()));
    let tmp = dir.path().join("spill");
    std::fs::create_dir(&tmp).unwrap();
    for cmd in ["merge", "complement"] {
        let small = gloci()
            .args([cmd, "--sort", "--max-mem", "16M", "-v", &big, "--tmp-dir"])
            .arg(&tmp)
            .assert()
            .success();
        let small = small.get_output().clone();
        let log = String::from_utf8_lossy(&small.stderr);
        assert!(
            log.contains("chunks spilled") && !log.contains(" 0 chunks spilled"),
            "{cmd} should have spilled: {log}"
        );
        let roomy = run(&[cmd, "--sort", &big]);
        assert_eq!(small.stdout, roomy, "{cmd}");
        assert!(!roomy.is_empty());
        // Nothing is left in the temp directory.
        assert_eq!(std::fs::read_dir(&tmp).unwrap().count(), 0, "{cmd}");
    }
}

#[test]
fn max_mem_too_small_is_usage_error() {
    let dir = tempdir();
    let f = s(&plain_file(dir.path(), "a.bed", "chr1\t1\t5\n"));
    gloci()
        .args(["merge", &f, "--sort", "--max-mem", "1M"])
        .assert()
        .code(2)
        .stderr(contains("at least 16M"));
}

#[test]
fn mixed_formats_with_sort_is_usage_error() {
    let dir = tempdir();
    let bed = s(&plain_file(dir.path(), "a.bed", "chr1\t1\t5\n"));
    let vcf = s(&plain_file(
        dir.path(),
        "a.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t5\t.\tA\tG\t.\t.\t.\n",
    ));
    gloci()
        .args(["merge", &bed, &vcf, "--sort"])
        .assert()
        .code(2)
        .stderr(contains("one format"));
}
