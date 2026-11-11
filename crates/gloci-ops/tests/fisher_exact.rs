//! `stats fisher`: klib's exact test against a direct hypergeometric sum,
//! bedtools' n22 estimate, and the genome size.

use gloci_core::{Genome, OverlapCriteria, VcfSpan};
use gloci_io::{Format, InvalidPolicy, ReaderOptions};
use gloci_ops::RelateInputs;
use gloci_ops::stats::{
    FisherCounts, FisherOptions, fisher_counts, fisher_table, kt_fisher_exact, write_fisher,
};
use gloci_testkit::mem::mem_reader;
use std::num::NonZeroUsize;

fn reader_opts() -> ReaderOptions {
    ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: None,
        threads: NonZeroUsize::MIN,
        format: Some(Format::Bed),
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    }
}

/// ln(k!) from a running sum of logs, independent of the port's `lgamma`.
fn ln_factorials(n: usize) -> Vec<f64> {
    let mut v = vec![0.0; n + 1];
    for i in 1..=n {
        v[i] = v[i - 1] + (i as f64).ln();
    }
    v
}

/// The hypergeometric probabilities of n11 = lo..=hi for the table's margins, and the observed one.
fn pmf(t: [i64; 4], lf: &[f64]) -> (i64, Vec<f64>, f64) {
    let [n11, n12, n21, n22] = t;
    let (r1, c1, n) = (n11 + n12, n11 + n21, n11 + n12 + n21 + n22);
    let lo = (r1 + c1 - n).max(0);
    let hi = r1.min(c1);
    let lc = |n: i64, k: i64| lf[n as usize] - lf[k as usize] - lf[(n - k) as usize];
    let p = |k: i64| (lc(r1, k) + lc(n - r1, c1 - k) - lc(n, c1)).exp();
    ((lo), (lo..=hi).map(p).collect(), p(n11))
}

fn close(got: f64, want: f64) -> bool {
    (got - want).abs() <= 1e-9 * want.abs().max(1e-12)
}

#[test]
fn kt_fisher_exact_matches_enumeration() {
    let lf = ln_factorials(200);
    let cells: Vec<i64> = (0..=40).step_by(3).collect();
    let mut checked = 0;
    for &a in &cells {
        for &b in &cells {
            for &c in &cells {
                for &d in &cells {
                    let (lo, ps, q) = pmf([a, b, c, d], &lf);
                    let hi = lo + ps.len() as i64 - 1;
                    let (left, right, two) = kt_fisher_exact(a, b, c, d);
                    if lo == hi {
                        // klib's `min == max` shortcut: nothing to test.
                        assert_eq!((left, right, two), (1.0, 1.0, 1.0), "{a} {b} {c} {d}");
                        continue;
                    }
                    let k = (a - lo) as usize;
                    let want_left: f64 = ps[..=k].iter().sum();
                    let want_right: f64 = ps[k..].iter().sum();
                    let want_two: f64 = ps.iter().filter(|&&p| p <= q * (1.0 + 1e-7)).sum();
                    let ctx = format!("table {a} {b} {c} {d}");
                    assert!(
                        close(left, want_left.min(1.0)),
                        "{ctx}: left {left} vs {want_left}"
                    );
                    assert!(
                        close(right, want_right.min(1.0)),
                        "{ctx}: right {right} vs {want_right}"
                    );
                    assert!(
                        close(two, want_two.min(1.0)),
                        "{ctx}: two {two} vs {want_two}"
                    );
                    checked += 1;
                }
            }
        }
    }
    assert!(checked > 30_000, "{checked}");
}

/// Fisher's tea-tasting table: P(X >= 3) = 17/70, two-sided 34/70.
#[test]
fn kt_fisher_exact_published_table() {
    let (left, right, two) = kt_fisher_exact(3, 1, 1, 3);
    assert!(close(right, 17.0 / 70.0), "{right}");
    assert!(close(left, 69.0 / 70.0), "{left}");
    assert!(close(two, 34.0 / 70.0), "{two}");
}

/// klib's `q == 0` branch (the table's probability underflows): the 1 goes to the side the mode
/// is on, as in bedtools 2.31.1's kfunc.cpp.
#[test]
fn kt_fisher_exact_underflow_branch() {
    // Enriched: all 3000 of A's and B's intervals shared, in a small universe.
    assert_eq!(kt_fisher_exact(3000, 0, 0, 100_000), (1.0, 0.0, 0.0));
    // Depleted: none shared although 3000 each would share about 1475.
    assert_eq!(kt_fisher_exact(0, 3000, 3000, 100), (0.0, 1.0, 0.0));
    // The reviewer's table: 5000 x 5000 with 300 shared on 10^8 bases.
    assert_eq!(
        kt_fisher_exact(300, 4700, 4700, 99_990_300),
        (1.0, 0.0, 0.0)
    );
}

fn counts(n_a: u64, len_a: u64, n_b: u64, len_b: u64, overlaps: u64, g: u64) -> FisherCounts {
    FisherCounts {
        n_a,
        n_b,
        overlaps,
        len_a,
        len_b,
        genome_size: g,
    }
}

#[test]
fn fisher_table_heuristic() {
    let t = fisher_table(&counts(10, 1000, 20, 1000, 4, 100_000));
    assert_eq!((t.n11, t.n12, t.n21, t.n22), (4, 6, 16, 631));
    assert_eq!(t.n22_full, 657);
    // 1000 / (2 + 4/3 + 5/3) is exactly 200 (a binary estimate may give 199).
    let t = fisher_table(&counts(3, 4, 3, 5, 0, 1000));
    assert_eq!(t.n22_full, 200);
    // An empty A (or B): the estimate is NaN, so the floor n11 + n12 + n21 stays.
    let t = fisher_table(&counts(0, 0, 7, 70, 0, 1000));
    assert_eq!((t.n22_full, t.n22), (7, 0));
    let t = fisher_table(&counts(5, 50, 0, 0, 0, 1000));
    assert_eq!((t.n22_full, t.n22), (5, 0));
    // A tiny genome never makes n22 negative.
    let t = fisher_table(&counts(5, 500, 5, 500, 2, 10));
    assert_eq!((t.n22_full, t.n22), (8, 0));
    // The u128 numerator overflows: the f64 fallback still gives the quotient.
    let big = 1u64 << 62;
    let t = fisher_table(&counts(big, big, big, big, 0, u64::MAX));
    assert!(t.n22_full > 0 && t.n22 >= 0);
}

fn run(a: &str, b: &str, genome: &[(&str, Option<u64>)], exclude: Option<&str>) -> FisherCounts {
    let ro = reader_opts();
    let g = Genome::new(
        genome
            .iter()
            .map(|(n, l)| (n.as_bytes().to_vec(), *l))
            .collect(),
    )
    .unwrap();
    let dir = tempfile::tempdir().unwrap();
    let exclude = exclude.map(|t| {
        let p = dir.path().join("x.bed");
        std::fs::write(&p, t).unwrap();
        p.to_string_lossy().into_owned()
    });
    let o = FisherOptions {
        criteria: OverlapCriteria::default(),
        merge: false,
        exclude,
        sort: None,
        genome: None,
    };
    fisher_counts(
        RelateInputs {
            a: mem_reader("a.bed", a, &ro),
            b: vec![mem_reader("b.bed", b, &ro)],
        },
        &g,
        &o,
    )
    .unwrap()
}

#[test]
fn fisher_counts_and_exclude() {
    let g = [("chr1", Some(1000)), ("chr2", Some(500))];
    let c = run("chr1\t0\t10\nchr1\t20\t40\n", "chr1\t5\t30\n", &g, None);
    assert_eq!(
        (c.n_a, c.n_b, c.overlaps, c.len_a, c.len_b, c.genome_size),
        (2, 1, 2, 30, 25, 1500)
    );
    // Overlapping and bookended exclude records count once: 100-300 and
    // 300-400 and 350-500 merge into 100-500.
    let ex = "chr1\t100\t300\nchr1\t300\t400\nchr1\t350\t500\nchr2\t0\t100\n";
    let c = run("chr1\t0\t10\n", "chr1\t5\t30\n", &g, Some(ex));
    assert_eq!(c.genome_size, 1500 - 400 - 100);
}

#[test]
fn fisher_genome_without_lengths_is_a_usage_error() {
    let ro = reader_opts();
    let g = Genome::new(vec![(b"chr1".to_vec(), None)]).unwrap();
    let o = FisherOptions {
        criteria: OverlapCriteria::default(),
        merge: false,
        exclude: None,
        sort: None,
        genome: None,
    };
    let e = fisher_counts(
        RelateInputs {
            a: mem_reader("a.bed", "chr1\t0\t10\n", &ro),
            b: vec![mem_reader("b.bed", "chr1\t5\t30\n", &ro)],
        },
        &g,
        &o,
    )
    .unwrap_err();
    assert!(e.to_string().contains("no length for chr1"), "{e}");
}

#[test]
fn fisher_report_text() {
    let c = counts(3, 300, 3, 350, 2, 100_000);
    let t = fisher_table(&c);
    let mut out = Vec::new();
    write_fisher(&c, &t, &mut out);
    let s = String::from_utf8(out).unwrap();
    assert!(s.starts_with("# Number of query intervals: 3\n"), "{s}");
    assert!(
        s.contains("# not in -a | 1            | 453          |\n"),
        "{s}"
    );
    assert!(
        s.ends_with("left\tright\ttwo-tail\tratio\n1\t8.6249e-05\t8.6249e-05\t906.000\n"),
        "{s}"
    );
}
