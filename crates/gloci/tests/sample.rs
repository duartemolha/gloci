//! `gloci sample` (D151): uniform subsets, input order, all records, too few, strand and VCF.
//! bedtools parity of the deterministic cases is in the golden cases (`sample-*-vs-bedtools`).
mod common;
use common::*;
use gloci_core::ChromTable;
use gloci_io::{IoError, ReaderOptions, RecordWriter, open_reader};
use gloci_ops::sample::{SampleOptions, sample};
use predicates::str::contains;

fn bed(n: usize) -> String {
    (0..n)
        .map(|i| format!("chr1\t{}\t{}\tr{i}\t{i}\t+\n", i * 10, i * 10 + 5))
        .collect()
}

#[derive(Default)]
struct Collect(Vec<Vec<u8>>);

impl RecordWriter for Collect {
    fn write_record(&mut self, raw: &[u8]) -> std::io::Result<()> {
        self.0.push(raw.to_vec());
        Ok(())
    }
    fn finish(self: Box<Self>) -> Result<(), IoError> {
        Ok(())
    }
}

fn lib_sample(path: &str, count: u64, seed: u64) -> Vec<Vec<u8>> {
    let ro = ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: None,
        threads: std::num::NonZeroUsize::MIN,
        format: None,
        vcf_span: gloci_core::VcfSpan::Allele,
        invalid: gloci_io::InvalidPolicy::Fail,
    };
    let r = open_reader(path, &ro).unwrap();
    let mut chroms = ChromTable::first_appearance();
    let mut w = Collect::default();
    let o = SampleOptions {
        count,
        seed,
        strand: None,
    };
    sample(r, &mut chroms, &o, &mut w).unwrap();
    w.0
}

fn run(args: &[&str]) -> String {
    let o = gloci().arg("sample").args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

#[test]
fn sample_is_uniform() {
    let d = tempdir();
    let f = plain_file(d.path(), "a.bed", &bed(1000));
    let f = f.to_string_lossy().into_owned();
    let mut counts = vec![0f64; 1000];
    for seed in 1..=2000u64 {
        let kept = lib_sample(&f, 100, seed);
        assert_eq!(kept.len(), 100);
        for r in kept {
            let s = String::from_utf8(r).unwrap();
            let i: usize = s.split('\t').nth(4).unwrap().parse().unwrap();
            counts[i] += 1.0;
        }
    }
    let chi2: f64 = counts
        .iter()
        .map(|c| (c - 200.0) * (c - 200.0) / 200.0)
        .sum();
    assert!(chi2 < 1174.0, "inclusion chi2 {chi2}");
}

#[test]
fn sample_input_order_and_determinism() {
    let d = tempdir();
    let f = plain_file(d.path(), "a.bed", &bed(300));
    let f = f.to_string_lossy().into_owned();
    let a = run(&[&f, "-n", "40", "--seed", "11", "-t", "1"]);
    let ids: Vec<usize> = a
        .lines()
        .map(|l| l.split('\t').nth(4).unwrap().parse().unwrap())
        .collect();
    assert_eq!(ids.len(), 40);
    assert!(ids.windows(2).all(|w| w[0] < w[1]), "{ids:?}");
    for t in ["2", "8"] {
        assert_eq!(a, run(&[&f, "-n", "40", "--seed", "11", "-t", t]));
    }
    assert_ne!(a, run(&[&f, "-n", "40", "--seed", "12"]));
}

#[test]
fn sample_all_records() {
    let d = tempdir();
    let text = format!("track name=x\n#c\n{}", bed(25));
    let f = plain_file(d.path(), "a.bed", &text);
    let out = run(&[f.to_str().unwrap(), "-n", "25", "--seed", "1"]);
    assert_eq!(out, bed(25));
    let out = run(&[f.to_str().unwrap(), "-n", "25", "--seed", "1", "--header"]);
    assert_eq!(out, text);
}

#[test]
fn sample_too_few() {
    let d = tempdir();
    let f = plain_file(d.path(), "a.bed", &bed(37));
    let out = d.path().join("out.bed");
    gloci()
        .args(["sample"])
        .arg(&f)
        .args(["-n", "100", "--seed", "1", "-o"])
        .arg(&out)
        .assert()
        .code(1)
        .stderr(contains("37 records"));
    assert!(!out.exists(), "no output file on error");
}

#[test]
fn sample_of_empty_input_writes_nothing_and_warns_once() {
    let d = tempdir();
    let f = plain_file(d.path(), "e.bed", "");
    let o = gloci()
        .args(["sample"])
        .arg(&f)
        .args(["-n", "5", "--seed", "1"])
        .assert()
        .success();
    assert!(o.get_output().stdout.is_empty());
    let e = String::from_utf8(o.get_output().stderr.clone()).unwrap();
    assert_eq!(e.lines().count(), 1, "{e}");
    assert!(e.contains("e.bed is empty: it has no records"), "{e}");
    assert!(!e.contains("only 0"), "{e}");
}

#[test]
fn sample_strand_filter_emptying_a_non_empty_input_is_still_too_few() {
    let d = tempdir();
    let f = plain_file(
        d.path(),
        "minus.bed",
        "chr1\t0\t5\ta\t0\t-\nchr1\t5\t9\tb\t0\t-\n",
    );
    gloci()
        .args(["sample"])
        .arg(&f)
        .args(["-n", "1", "--strand", "+", "--seed", "1"])
        .assert()
        .code(1)
        .stderr(contains("only 0 records on the + strand"));
}

#[test]
fn sample_strand() {
    let d = tempdir();
    let text = "chr1\t0\t5\ta\t0\t+\nchr1\t5\t9\tb\t0\t-\nchr1\t9\t12\tc\t0\t.\n\
chr1\t12\t15\td\t0\t+\nchr1\t15\t19\te\t0\t-\n";
    let f = plain_file(d.path(), "a.bed", text);
    let f = f.to_str().unwrap();
    assert_eq!(
        run(&[f, "-n", "2", "--strand", "+", "--seed", "1"]),
        "chr1\t0\t5\ta\t0\t+\nchr1\t12\t15\td\t0\t+\n"
    );
    assert_eq!(
        run(&[f, "-n", "2", "--strand", "-", "--seed", "1"]),
        "chr1\t5\t9\tb\t0\t-\nchr1\t15\t19\te\t0\t-\n"
    );
    // Unknown never complies: 2 of the 5 records are on each strand.
    gloci()
        .args(["sample", f, "-n", "3", "--strand", "+", "--seed", "1"])
        .assert()
        .code(1)
        .stderr(contains("2 records"));
}

#[test]
fn sample_vcf_keeps_header() {
    let d = tempdir();
    let mut text = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    for i in 1..=6 {
        text.push_str(&format!("chr1\t{}\t.\tA\tC\t.\t.\t.\n", i * 10));
    }
    let f = plain_file(d.path(), "a.vcf", &text);
    let out = run(&[f.to_str().unwrap(), "-n", "2", "--seed", "4"]);
    assert!(out.starts_with("##fileformat=VCFv4.2\n"), "{out}");
    assert!(out.contains("#CHROM\tPOS"), "{out}");
    assert_eq!(out.lines().filter(|l| !l.starts_with('#')).count(), 2);
}
