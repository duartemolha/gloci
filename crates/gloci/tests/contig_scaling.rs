//! Streaming relation commands must scale linearly with the number of contigs
//! (an assembly with hundreds of thousands of scaffolds). Ignored in normal
//! runs because it times commands; run it in release:
//! `cargo test --release -p gloci --test contig_scaling -- --ignored`.
//! Before the fix, 80,000 contigs took 8 to 13 s and 400,000 did not finish
//! in 100 s (a linear scan of the names seen so far, once per contig).
mod common;
use common::*;
use std::fmt::Write as _;
use std::time::Instant;

const CONTIGS: usize = 200_000;

#[test]
#[ignore = "timing test: run with --release -- --ignored"]
fn two_hundred_thousand_contigs_complete_in_seconds() {
    let dir = tempdir();
    let mut names: Vec<String> = (0..CONTIGS).map(|i| format!("ctg{i}")).collect();
    names.sort();
    let (mut a, mut b) = (String::new(), String::new());
    for n in &names {
        writeln!(a, "{n}\t100\t200").unwrap();
        writeln!(b, "{n}\t150\t250").unwrap();
    }
    let a = plain_file(dir.path(), "a.bed", &a);
    let b = plain_file(dir.path(), "b.bed", &b);
    let (a, b) = (a.to_str().unwrap(), b.to_str().unwrap());
    let cases: &[&[&str]] = &[
        &["intersect", "-t", "1", "-a", a, "-b", b],
        &["merge", a, b],
        &["partition", a, b],
        &["subtract", "-t", "1", "-a", a, "-b", b],
        &["closest", "-t", "1", "-a", a, "-b", b],
        &["map", "-t", "1", "-c", "2:sum", "-a", a, "-b", b],
        &["coverage", "-t", "1", "-a", a, "-b", b],
    ];
    for args in cases {
        let t0 = Instant::now();
        let out = gloci()
            .args(args[..1].iter().copied())
            .args(["--chrom-order", "lex"])
            .args(&args[1..])
            .assert()
            .success()
            .get_output()
            .stdout
            .clone();
        let secs = t0.elapsed().as_secs_f64();
        eprintln!("{}: {secs:.2} s", args[0]);
        assert!(!out.is_empty(), "{} wrote nothing", args[0]);
        // A generous bound: the quadratic code needed minutes at this size.
        assert!(
            secs < 30.0,
            "{} took {secs:.1} s for {CONTIGS} contigs",
            args[0]
        );
    }
}
