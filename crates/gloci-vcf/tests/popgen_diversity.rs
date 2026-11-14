//! Nucleotide diversity and Tajima's D (decision D159) on hand-made sites;
//! the vcftools 0.1.16 comparison is `crates/gloci/tests/vcf_popgen.rs`.

use gloci_core::Interval;
use gloci_vcf::popgen::{SitePi, TajimaD, WindowPi};
use gloci_vcf::{MaskRules, SiteView, Stat, StatSink, StatWrite, site_stats, text_genotypes};
use std::cell::RefCell;
use std::io::{self, Write};
use std::rc::Rc;

/// Captures the one file a statistic writes.
#[derive(Default, Clone)]
struct MemSink(Rc<RefCell<Vec<u8>>>);

struct MemFile(MemSink, Vec<u8>);

impl Write for MemFile {
    fn write(&mut self, b: &[u8]) -> io::Result<usize> {
        self.1.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

impl StatWrite for MemFile {
    fn commit(self: Box<Self>) -> io::Result<()> {
        *self.0.0.borrow_mut() = self.1;
        Ok(())
    }
}

impl StatSink for MemSink {
    fn create(&mut self, _ext: &str) -> io::Result<Box<dyn StatWrite>> {
        Ok(Box::new(MemFile(self.clone(), Vec::new())))
    }
}

const NO_MASK: MaskRules = MaskRules {
    min_gq: None,
    min_dp: None,
    max_dp: None,
    mask_filtered: false,
};

/// A record `CHROM POS . A ALTS . . . GT:DP SAMPLES...`; each sample is
/// `GT` or `GT:DP`.
fn record(chrom: &str, pos: u64, alts: &str, samples: &[&str]) -> String {
    format!(
        "{chrom}\t{pos}\t.\tA\t{alts}\t.\t.\t.\tGT:DP\t{}",
        samples.join("\t")
    )
}

/// Feed `lines` to `stat` under `mask` and return the file it writes.
fn run(mut stat: impl Stat, lines: &[String], mask: &MaskRules) -> String {
    let sink = MemSink::default();
    stat.start(&mut sink.clone()).unwrap();
    let mut masked = Vec::new();
    for line in lines {
        let pos: u64 = line.split('\t').nth(1).unwrap().parse().unwrap();
        let site =
            SiteView::from_text(line.as_bytes(), Interval::new(pos - 1, pos).unwrap()).unwrap();
        let n_alleles = if site.alts == b"." {
            1
        } else {
            1 + site.alts.split(|&b| b == b',').count()
        };
        let g = text_genotypes(line.as_bytes()).unwrap();
        let st = site_stats(&g, None, mask, n_alleles, &mut masked).unwrap();
        stat.site(&site, None, Some(&st)).unwrap();
    }
    stat.finish().unwrap();
    String::from_utf8(sink.0.take()).unwrap()
}

#[test]
fn site_pi_table() {
    let lines = [
        record("chr1", 1, "G", &["0/0", "0/1"]), // counts (3, 1)
        record("chr1", 2, "G", &["0/0", "0/0"]), // monomorphic
        record("chr1", 3, "G", &["./1", "./."]), // one called allele
        record("chr1", 4, "G", &["1", "0/1"]),   // a haploid sample: no row
        // the masked haploid sample is left out: its site counts, its allele not
        record("chr1", 5, "G", &["1:3", "0/1:9", "0/0:9"]),
    ];
    let mask = MaskRules {
        min_dp: Some(5),
        ..NO_MASK
    };
    assert_eq!(
        run(SitePi::new(), &lines, &mask),
        "CHROM\tPOS\tPI\nchr1\t1\t0.5\nchr1\t2\t0\nchr1\t3\tnan\nchr1\t5\t0.5\n"
    );
    // Unmasked, the haploid sample keeps the last site out.
    assert_eq!(
        run(SitePi::new(), &lines, &NO_MASK),
        "CHROM\tPOS\tPI\nchr1\t1\t0.5\nchr1\t2\t0\nchr1\t3\tnan\n"
    );
}

#[test]
fn window_pi_table() {
    let head = "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI\n";
    // 2 samples: 6 differing pairs of 12, plus 9 monomorphic positions of
    // 4 × 3 pairs each: 6 / 120.
    let at5 = [record("chr1", 5, "G", &["0/0", "0/1"])];
    assert_eq!(
        run(WindowPi::new(10, 10), &at5, &NO_MASK),
        format!("{head}chr1\t1\t10\t1\t0.05\n")
    );
    // Step 5: POS 7 is in 1-10 and 6-15 only (not 11-20, nor a window
    // starting below 1).
    let at7 = [record("chr1", 7, "G", &["0/0", "0/1"])];
    assert_eq!(
        run(WindowPi::new(10, 5), &at7, &NO_MASK),
        format!("{head}chr1\t1\t10\t1\t0.05\nchr1\t6\t15\t1\t0.05\n")
    );
    // A fixed or monomorphic site and a haploid one add nothing: no window
    // holds a polymorphic site, so there is no row.
    let none = [
        record("chr1", 3, "G", &["1/1", "1/1"]),
        record("chr1", 13, "G", &["0/0", "0/0"]),
        record("chr1", 23, "G", &["0", "0/1"]),
    ];
    assert_eq!(run(WindowPi::new(10, 10), &none, &NO_MASK), head);
}

#[test]
fn tajima_bin_reciprocal() {
    // 49 × (1/49) is just below 1 in floating point: bin 0, and 98 is in
    // bin 1, where integer division would give bins 1 and 2.
    let two = ["0/1", "0/0"];
    let at49 = [record("chr1", 49, "G", &two)];
    let got = run(TajimaD::new(49, 2), &at49, &NO_MASK);
    assert!(got.contains("\nchr1\t0\t1\t"), "{got}");
    let at98 = [record("chr1", 98, "G", &two)];
    let got = run(TajimaD::new(49, 2), &at98, &NO_MASK);
    assert!(got.contains("\nchr1\t49\t1\t"), "{got}");
    assert_eq!(got.lines().count(), 2, "{got}");
}

#[test]
fn tajima_trailing_bins() {
    let lines = [
        record("chr1", 5, "G", &["0/0", "0/0"]), // bin 0, monomorphic: no row
        record("chr1", 12, "G", &["0/0", "0/0"]), // bin 1, monomorphic: no row
        record("chr1", 25, "G", &["0/1", "0/0"]), // bin 2, S = 1
        record("chr1", 52, "G", &["1/1", "1/1"]), // bin 5, fixed: ends the rows
        record("chr1", 75, "G,C", &["0/1", "0/2"]), // multi-allelic: not a bin
        record("chr1", 85, "G", &["1", "0/1"]),  // a haploid sample: not a bin
    ];
    let got = run(TajimaD::new(10, 2), &lines, &NO_MASK);
    let rows: Vec<&str> = got.lines().skip(1).collect();
    assert_eq!(rows.len(), 4, "{got}");
    assert!(rows[0].starts_with("chr1\t20\t1\t"), "{got}");
    assert_eq!(
        &rows[1..],
        ["chr1\t30\t0\tnan", "chr1\t40\t0\tnan", "chr1\t50\t0\tnan"]
    );
    // A second chromosome starts its own rows at its first S > 0 bin.
    let mut two = lines.to_vec();
    two.push(record("chr2", 3, "G", &["0/0", "0/0"]));
    two.push(record("chr2", 31, "G", &["0/1", "1/1"]));
    let got = run(TajimaD::new(10, 2), &two, &NO_MASK);
    let last: Vec<&str> = got.lines().skip(5).collect();
    assert_eq!(last.len(), 1, "{got}");
    assert!(last[0].starts_with("chr2\t30\t1\t"), "{got}");
}

/// Tajima (1989), equations 3, 5 and 38 to 41, with π the mean number of
/// differences over all pairs of sequences.
fn tajima_1989(seqs: &[&[u8]]) -> f64 {
    let n = seqs.len();
    let sites = seqs[0].len();
    let s = (0..sites)
        .filter(|&j| seqs.iter().any(|q| q[j] != seqs[0][j]))
        .count() as f64;
    let mut diffs = 0usize;
    for i in 0..n {
        for k in i + 1..n {
            diffs += (0..sites).filter(|&j| seqs[i][j] != seqs[k][j]).count();
        }
    }
    let pi = diffs as f64 / (n * (n - 1) / 2) as f64;
    let nf = n as f64;
    let a1: f64 = (1..n).map(|i| 1.0 / i as f64).sum();
    let a2: f64 = (1..n).map(|i| 1.0 / (i * i) as f64).sum();
    let b1 = (nf + 1.0) / (3.0 * (nf - 1.0));
    let b2 = 2.0 * (nf * nf + nf + 3.0) / (9.0 * nf * (nf - 1.0));
    let c1 = b1 - 1.0 / a1;
    let c2 = b2 - (nf + 2.0) / (a1 * nf) + a2 / (a1 * a1);
    let e1 = c1 / a1;
    let e2 = c2 / (a1 * a1 + a2);
    (pi - s / a1) / (e1 * s + e2 * s * (s - 1.0)).sqrt()
}

#[test]
fn tajima_oracle() {
    // Ten chromosomes (five diploid samples, nothing missing), six sites;
    // site 4 is monomorphic.
    let seqs: [&[u8]; 10] = [
        b"ACGTAC", b"ACGTAC", b"TCGTAC", b"TCGTAA", b"TAGTAA", b"TAGTAA", b"TACTAA", b"AACTAC",
        b"ACCTAC", b"ACGTGC",
    ];
    let want = tajima_1989(&seqs);
    // The statistic's inputs: S and Σ p (1 − p), p the first sequence's
    // base frequency at each site.
    let (mut s, mut sum_pq) = (0u64, 0.0f64);
    for j in 0..seqs[0].len() {
        let same = seqs.iter().filter(|q| q[j] == seqs[0][j]).count();
        let p = same as f64 / seqs.len() as f64;
        if p < 1.0 {
            s += 1;
            sum_pq += p * (1.0 - p);
        }
    }
    let got = TajimaD::new(100, 5).d(s, sum_pq);
    assert!((got - want).abs() < 1e-12, "{got} vs {want}");
    assert!(want.is_finite() && want != 0.0);
}
