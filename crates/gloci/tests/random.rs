//! `gloci random` (D150): the pinned stream, uniformity, constraints, errors and thread
//! independence. bedtools' random output is not comparable (docs/differences.md#rng-stream).
mod common;
use common::*;
use predicates::str::contains;
use std::time::Instant;

const G3: &str = "chr1\t100\nchr2\t200\nchr3\t700\n";

fn genome(dir: &std::path::Path) -> String {
    plain_file(dir, "g.genome", G3)
        .to_string_lossy()
        .into_owned()
}

fn run(args: &[&str]) -> String {
    let o = gloci().arg("random").args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

#[test]
fn random_pinned() {
    let d = tempdir();
    let g = genome(d.path());
    let out = run(&["-g", &g, "-l", "10", "-n", "5", "--seed", "42"]);
    assert_eq!(
        out,
        "chr3\t56\t66\t1\t10\t-\nchr3\t473\t483\t2\t10\t-\nchr2\t130\t140\t3\t10\t+\n\
chr3\t162\t172\t4\t10\t+\nchr3\t465\t475\t5\t10\t+\n"
    );
}

#[test]
fn random_uniform() {
    let d = tempdir();
    let g = genome(d.path());
    let n = 100_000usize;
    let out = run(&["-g", &g, "-l", "1", "-n", &n.to_string(), "--seed", "7"]);
    let mut per_chrom = [0f64; 3];
    let mut bins = [0f64; 10];
    let mut ids = 0u64;
    for line in out.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let (start, end): (u64, u64) = (f[1].parse().unwrap(), f[2].parse().unwrap());
        ids += 1;
        assert_eq!(f[3].parse::<u64>().unwrap(), ids);
        assert_eq!(f[4], "1");
        assert!(f[5] == "+" || f[5] == "-");
        assert_eq!(end, start + 1);
        let (c, len) = match f[0] {
            "chr1" => (0, 100),
            "chr2" => (1, 200),
            "chr3" => (2, 700),
            other => panic!("unexpected chromosome {other}"),
        };
        assert!(end <= len);
        per_chrom[c] += 1.0;
        if c == 2 {
            bins[(start / 70) as usize] += 1.0;
        }
    }
    assert_eq!(ids as usize, n);
    let chi2 = |obs: &[f64], exp: &[f64]| -> f64 {
        obs.iter()
            .zip(exp)
            .map(|(o, e)| (o - e) * (o - e) / e)
            .sum()
    };
    let nf = n as f64;
    let c = chi2(&per_chrom, &[nf * 0.1, nf * 0.2, nf * 0.7]);
    assert!(c < 18.42, "chromosome counts chi2 {c}");
    let total3: f64 = per_chrom[2];
    let b = chi2(&bins, &[total3 / 10.0; 10]);
    assert!(b < 33.72, "chr3 start bins chi2 {b}");
}

#[test]
fn random_stays_inside_when_some_chromosomes_are_short() {
    let d = tempdir();
    let g = genome(d.path());
    let out = run(&["-g", &g, "-l", "150", "-n", "2000", "--seed", "3"]);
    for line in out.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        assert_ne!(f[0], "chr1", "{line}");
        let len = if f[0] == "chr2" { 200 } else { 700 };
        assert!(f[2].parse::<u64>().unwrap() <= len, "{line}");
        assert_eq!(
            f[2].parse::<u64>().unwrap() - f[1].parse::<u64>().unwrap(),
            150
        );
    }
}

#[test]
fn random_too_long() {
    let d = tempdir();
    let g = genome(d.path());
    gloci()
        .args(["random", "-g", &g, "-l", "800", "-n", "5", "--seed", "1"])
        .assert()
        .code(2)
        .stderr(contains("chr3 (700 bp)"));
}

#[test]
fn random_needs_lengths() {
    let d = tempdir();
    // A VCF header whose contigs carry no length.
    let g = plain_file(
        d.path(),
        "nolen.vcf",
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    gloci()
        .args(["random", "-g"])
        .arg(&g)
        .args(["-n", "5", "--seed", "1"])
        .assert()
        .code(2)
        .stderr(contains("no length for chr1"));
}

#[test]
fn random_huge_contig_count() {
    let d = tempdir();
    let mut text = String::new();
    for i in 0..200_000 {
        text.push_str(&format!("ctg{i}\t{}\n", 1000 + i % 7));
    }
    let g = plain_file(d.path(), "many.genome", &text);
    let t = Instant::now();
    let o = gloci()
        .args(["random", "-g"])
        .arg(&g)
        .args(["-l", "100", "-n", "10000", "--seed", "5"])
        .assert()
        .success();
    assert_eq!(o.get_output().stdout.split(|&b| b == b'\n').count(), 10_001);
    assert!(t.elapsed().as_secs_f64() < 2.0, "took {:?}", t.elapsed());
}

#[test]
fn random_same_bytes_at_any_thread_count() {
    let d = tempdir();
    let g = genome(d.path());
    let base = run(&[
        "-g", &g, "-l", "20", "-n", "300000", "--seed", "9", "-t", "1",
    ]);
    for t in ["2", "8"] {
        let other = run(&["-g", &g, "-l", "20", "-n", "300000", "--seed", "9", "-t", t]);
        assert!(base == other, "-t {t} differs from -t 1");
    }
}

#[test]
fn seed_warning_reproduces() {
    let d = tempdir();
    let g = genome(d.path());
    let o = gloci()
        .args(["random", "-g", &g, "-l", "10", "-n", "50"])
        .assert()
        .success();
    let first = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    let err = String::from_utf8(o.get_output().stderr.clone()).unwrap();
    let seed = err
        .split("using --seed ")
        .nth(1)
        .and_then(|s| s.split(' ').next())
        .unwrap_or_else(|| panic!("no seed in {err:?}"))
        .to_string();
    let again = run(&["-g", &g, "-l", "10", "-n", "50", "--seed", &seed]);
    assert_eq!(first, again);
}

#[test]
fn random_zero_count_and_zero_length() {
    let d = tempdir();
    let g = genome(d.path());
    assert_eq!(run(&["-g", &g, "-n", "0", "--seed", "1"]), "");
    gloci()
        .args(["random", "-g", &g, "-l", "0", "--seed", "1"])
        .assert()
        .code(2)
        .stderr(contains("-l"));
}

#[test]
fn random_refuses_write_index_up_front() {
    let d = tempdir();
    let g = genome(d.path());
    let out = d.path().join("x.bed.gz");
    gloci()
        .args([
            "random",
            "-g",
            &g,
            "-n",
            "5",
            "--seed",
            "1",
            "--write-index",
            "-o",
        ])
        .arg(&out)
        .assert()
        .code(2)
        .stderr(contains("random output is unsorted"))
        .stderr(contains("gloci sort"));
    assert!(!out.exists(), "nothing is created");
}
