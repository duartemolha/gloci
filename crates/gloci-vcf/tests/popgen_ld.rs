//! Linkage disequilibrium (decision D161) on hand-made sites: genotype and
//! haplotype r² against oracles written from the published formulas, the
//! pair windows, the genotype rules, the memory guard and the order of the
//! interchromosomal pairs. The vcftools 0.1.16 comparison is
//! `crates/gloci/tests/vcf_popgen.rs`.

use gloci_core::Interval;
use gloci_vcf::popgen::{GenoR2, HapR2, InterchromR2, LdBudget, LdWindow};
use gloci_vcf::{
    MaskRules, MaskedGenotypes, SiteView, Stat, StatSink, StatWrite, site_stats, text_genotypes,
};
use std::cell::RefCell;
use std::collections::HashMap;
use std::io::{self, Write};
use std::rc::Rc;

/// Captures every file a statistic writes, by extension.
#[derive(Default, Clone)]
struct MemSink(Rc<RefCell<HashMap<String, Vec<u8>>>>);

struct MemFile(MemSink, String, Vec<u8>);

impl Write for MemFile {
    fn write(&mut self, b: &[u8]) -> io::Result<usize> {
        self.2.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

impl StatWrite for MemFile {
    fn commit(self: Box<Self>) -> io::Result<()> {
        let MemFile(sink, ext, data) = *self;
        sink.0.borrow_mut().insert(ext, data);
        Ok(())
    }
}

impl StatSink for MemSink {
    fn create(&mut self, ext: &str) -> io::Result<Box<dyn StatWrite>> {
        Ok(Box::new(MemFile(self.clone(), ext.to_string(), Vec::new())))
    }
}

const NO_MASK: MaskRules = MaskRules {
    min_gq: None,
    min_dp: None,
    max_dp: None,
    mask_filtered: false,
};

const DP5: MaskRules = MaskRules {
    min_dp: Some(5),
    ..NO_MASK
};

/// A record `CHROM POS . A ALTS . . . GT:DP SAMPLES...`; each sample is
/// `GT` or `GT:DP`.
fn record(chrom: &str, pos: u64, alts: &str, samples: &[&str]) -> String {
    format!(
        "{chrom}\t{pos}\t.\tA\t{alts}\t.\t.\t.\tGT:DP\t{}",
        samples.join("\t")
    )
}

/// Feed `lines` to `stat` under `mask`; the file it writes, or the first
/// error.
fn try_run(stat: &mut dyn Stat, lines: &[String], mask: &MaskRules) -> Result<String, String> {
    let sink = MemSink::default();
    stat.start(&mut sink.clone()).unwrap();
    let mut masked = Vec::new();
    for line in lines {
        let pos: u64 = line.split('\t').nth(1).unwrap().parse().unwrap();
        let site =
            SiteView::from_text(line.as_bytes(), Interval::new(pos, pos + 1).unwrap()).unwrap();
        let n_alleles = if site.alts == b"." {
            1
        } else {
            1 + site.alts.split(|&b| b == b',').count()
        };
        let g = text_genotypes(line.as_bytes()).unwrap();
        let st = site_stats(&g, None, mask, n_alleles, &mut masked).unwrap();
        let mg = MaskedGenotypes {
            inner: &g,
            masked: &masked,
        };
        stat.site(&site, Some(&mg), Some(&st))?;
    }
    stat.finish().map_err(|e| e.to_string())?;
    let files = sink.0.take();
    assert_eq!(files.len(), 1, "one file: {:?}", files.keys());
    let (_, text) = files.into_iter().next().unwrap();
    Ok(String::from_utf8(text).unwrap())
}

fn run(stat: &mut dyn Stat, lines: &[String], mask: &MaskRules) -> String {
    try_run(stat, lines, mask).unwrap()
}

fn samples(n: usize) -> Vec<usize> {
    (0..n).collect()
}

fn budget() -> LdBudget {
    LdBudget::new(1 << 30)
}

fn geno(n: usize, window: LdWindow) -> GenoR2 {
    GenoR2::new(samples(n), window, budget())
}

fn hap(n: usize, window: LdWindow) -> HapR2 {
    HapR2::new(samples(n), window, budget())
}

/// The data rows of a file, without its header.
fn rows(text: &str) -> Vec<Vec<String>> {
    text.lines()
        .skip(1)
        .map(|l| l.split('\t').map(str::to_string).collect())
        .collect()
}

/// The (POS1, POS2) of every row of a windowed file.
fn pairs(text: &str) -> Vec<(u64, u64)> {
    rows(text)
        .iter()
        .map(|r| (r[1].parse().unwrap(), r[2].parse().unwrap()))
        .collect()
}

/// A deterministic genotype for the generated cases (no platform RNG).
struct Lcg(u64);

impl Lcg {
    fn next(&mut self, below: u64) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 33) % below
    }
}

/// The squared Pearson correlation of the ALT allele counts (0, 1, 2) of
/// two sites over the samples with both alleles called at both, and that
/// count: r = Σ(x − x̄)(y − ȳ) / √(Σ(x − x̄)² Σ(y − ȳ)²).
fn pearson_r2(a: &[&str], b: &[&str]) -> (usize, f64) {
    let dosage = |g: &str| -> Option<f64> {
        let (x, y) = g.split_once(['/', '|'])?;
        Some(x.parse::<f64>().ok()? + y.parse::<f64>().ok()?)
    };
    let both: Vec<(f64, f64)> = a
        .iter()
        .zip(b)
        .filter_map(|(x, y)| Some((dosage(x)?, dosage(y)?)))
        .collect();
    let n = both.len() as f64;
    let mx = both.iter().map(|p| p.0).sum::<f64>() / n;
    let my = both.iter().map(|p| p.1).sum::<f64>() / n;
    let sxy: f64 = both.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
    let sxx: f64 = both.iter().map(|p| (p.0 - mx).powi(2)).sum();
    let syy: f64 = both.iter().map(|p| (p.1 - my).powi(2)).sum();
    (both.len(), sxy * sxy / (sxx * syy))
}

/// Lewontin's D, D′ and r² of two biallelic sites from the haplotypes
/// called at both: with p₁ the REF frequency at the first site, q₁ at the
/// second and p₁₁ that of REF-REF haplotypes, D = p₁₁ − p₁q₁; D′ = D / Dmax,
/// Dmax = min(p₁q₂, p₂q₁) for D > 0 and min(p₁q₁, p₂q₂) for D < 0 (Lewontin
/// 1964); r² = D² / (p₁p₂q₁q₂). Haploid genotypes are one haplotype.
fn lewontin(a: &[&str], b: &[&str]) -> (usize, f64, f64, f64) {
    let haps = |g: &str| -> [Option<u32>; 2] {
        let mut it = g.split('|').map(|x| x.parse::<u32>().ok());
        [it.next().flatten(), it.next().flatten()]
    };
    let (mut n, mut n11, mut n1, mut m1) = (0usize, 0usize, 0usize, 0usize);
    for (x, y) in a.iter().zip(b) {
        for (u, v) in haps(x).into_iter().zip(haps(y)) {
            let (Some(u), Some(v)) = (u, v) else { continue };
            n += 1;
            n1 += usize::from(u == 0);
            m1 += usize::from(v == 0);
            n11 += usize::from(u == 0 && v == 0);
        }
    }
    let nf = n as f64;
    let (p1, q1, p11) = (n1 as f64 / nf, m1 as f64 / nf, n11 as f64 / nf);
    let (p2, q2) = (1.0 - p1, 1.0 - q1);
    let d = p11 - p1 * q1;
    let dmax = if d > 0.0 {
        (p1 * q2).min(p2 * q1)
    } else {
        (p1 * q1).min(p2 * q2)
    };
    (n, d * d / (p1 * p2 * q1 * q2), d, d / dmax)
}

fn close(got: &str, want: f64, what: &str) {
    if want.is_nan() {
        assert_eq!(got, "nan", "{what}");
        return;
    }
    let g: f64 = got.parse().unwrap_or_else(|_| panic!("{what}: {got}"));
    // Six printed digits: a relative error of 5e-6 at most, and noise of
    // 1e-12 around 0.
    assert!(
        (g - want).abs() <= 5e-6 * want.abs() + 1e-12,
        "{what}: printed {got}, formula {want}"
    );
}

#[test]
fn geno_r2_oracle() {
    let mut rng = Lcg(7);
    let gts = ["0/0", "0/1", "1/0", "1/1", "./.", "./1"];
    let n = 23;
    let mut sites: Vec<Vec<&str>> = Vec::new();
    for _ in 0..40 {
        // skew some sites so that alleles are rare or fixed
        let skew = rng.next(4);
        sites.push(
            (0..n)
                .map(|_| {
                    let g = rng.next(gts.len() as u64) as usize;
                    if skew == 0 && g > 0 && g < 4 {
                        "0/0"
                    } else {
                        gts[g]
                    }
                })
                .collect(),
        );
    }
    let lines: Vec<String> = sites
        .iter()
        .enumerate()
        .map(|(i, s)| record("1", 10 * (i as u64 + 1), "G", s))
        .collect();
    let out = run(&mut geno(n, LdWindow::default()), &lines, &NO_MASK);
    assert!(out.starts_with("CHR\tPOS1\tPOS2\tN_INDV\tR^2\n"), "{out}");
    let rows = rows(&out);
    assert_eq!(rows.len(), 40 * 39 / 2, "every pair once");
    let mut finite = 0;
    for r in &rows {
        let (i, j) = (
            r[1].parse::<usize>().unwrap() / 10 - 1,
            r[2].parse::<usize>().unwrap() / 10 - 1,
        );
        let (count, r2) = pearson_r2(&sites[i], &sites[j]);
        assert_eq!(r[3], count.to_string(), "N_INDV of {i},{j}");
        close(&r[4], r2, &format!("r2 of {i},{j}"));
        finite += usize::from(r2.is_finite());
    }
    assert!(finite > 400, "{finite} pairs with a value");
}

#[test]
fn hap_r2_oracle() {
    let mut rng = Lcg(11);
    let gts = ["0|0", "0|1", "1|0", "1|1", ".|1", "0|.", ".|.", "0", "1"];
    let n = 17;
    let mut sites: Vec<Vec<&str>> = Vec::new();
    for _ in 0..30 {
        sites.push(
            (0..n)
                .map(|_| gts[rng.next(gts.len() as u64) as usize])
                .collect(),
        );
    }
    // a site in perfect LD with the first, and one with REF only
    sites.push(sites[0].clone());
    sites.push(vec!["0|0"; n]);
    let lines: Vec<String> = sites
        .iter()
        .enumerate()
        .map(|(i, s)| record("1", i as u64 + 1, "G", s))
        .collect();
    let out = run(&mut hap(n, LdWindow::default()), &lines, &NO_MASK);
    assert!(
        out.starts_with("CHR\tPOS1\tPOS2\tN_CHR\tR^2\tD\tDprime\n"),
        "{out}"
    );
    let rows = rows(&out);
    assert_eq!(rows.len(), 32 * 31 / 2);
    let mut negative = 0;
    for r in &rows {
        let (i, j) = (
            r[1].parse::<usize>().unwrap() - 1,
            r[2].parse::<usize>().unwrap() - 1,
        );
        let (count, r2, d, dprime) = lewontin(&sites[i], &sites[j]);
        assert_eq!(r[3], count.to_string(), "N_CHR of {i},{j}");
        close(&r[4], r2, &format!("r2 of {i},{j}"));
        close(&r[5], d, &format!("D of {i},{j}"));
        if d.abs() > 1e-12 {
            close(&r[6], dprime, &format!("D' of {i},{j}"));
        }
        negative += usize::from(d < -1e-9);
    }
    assert!(negative > 50, "{negative} pairs with D < 0");
    let perfect = rows.iter().find(|r| r[1] == "1" && r[2] == "31").unwrap();
    assert_eq!(
        (&perfect[4][..], &perfect[6][..]),
        ("1", "1"),
        "{perfect:?}"
    );
    // REF only: no variance, D = 0 and Dmax = 0
    let fixed = rows.iter().find(|r| r[1] == "1" && r[2] == "32").unwrap();
    assert_eq!(&fixed[4..], ["nan", "0", "nan"], "{fixed:?}");
}

/// Sites at 1:100, 1:105, 1:120, 1:200, 1:300, 2:10, 2:15; one multi-allelic
/// and one ALT `.` site between, which take no index.
fn window_sites() -> Vec<String> {
    // r² of s and t is 0, of s and s or t and t 1
    let s = ["0/0", "0/1", "1/1", "0/1"];
    let t = ["0/1", "1/1", "0/1", "0/0"];
    let mono = ["0/0", "0/0", "0/0", "0/0"];
    vec![
        record("1", 100, "G", &s),
        record("1", 105, "G", &t),
        record("1", 110, "G,T", &["0/2", "1/1", "0/0", "2/2"]),
        record("1", 115, ".", &mono),
        record("1", 120, "G", &s),
        record("1", 200, "G", &t),
        record("1", 300, "G", &mono),
        record("2", 10, "G", &s),
        record("2", 15, "G", &t),
    ]
}

#[test]
fn ld_window_table() {
    let lines = window_sites();
    let w = |f: fn(&mut LdWindow)| {
        let mut w = LdWindow::default();
        f(&mut w);
        pairs(&run(&mut geno(4, w), &lines, &NO_MASK))
    };
    // Every pair on one chromosome, (i, j) order, none across chromosomes.
    let all = vec![
        (100, 105),
        (100, 120),
        (100, 200),
        (100, 300),
        (105, 120),
        (105, 200),
        (105, 300),
        (120, 200),
        (120, 300),
        (200, 300),
        (10, 15),
    ];
    assert_eq!(w(|_| {}), all);
    // j − i ≤ 2 sites: the two skipped sites take no index, so 105 and 200
    // are 2 apart.
    assert_eq!(
        w(|w| w.snps = 2),
        vec![
            (100, 105),
            (100, 120),
            (105, 120),
            (105, 200),
            (120, 200),
            (120, 300),
            (200, 300),
            (10, 15)
        ]
    );
    assert_eq!(
        w(|w| w.snps = 1),
        vec![(100, 105), (105, 120), (120, 200), (200, 300), (10, 15)]
    );
    // Δpos ≤ 100 (100 to 200 is in, 105 to 300 is not).
    assert_eq!(
        w(|w| w.bp = 100),
        vec![
            (100, 105),
            (100, 120),
            (100, 200),
            (105, 120),
            (105, 200),
            (120, 200),
            (200, 300),
            (10, 15)
        ]
    );
    // j − i ≥ 2; and a minimum of 0 is a minimum of 1.
    assert_eq!(
        w(|w| w.min_snps = 2),
        vec![
            (100, 120),
            (100, 200),
            (100, 300),
            (105, 200),
            (105, 300),
            (120, 300)
        ]
    );
    assert_eq!(w(|w| w.min_snps = 0), all);
    // Δpos ≥ 20.
    assert_eq!(
        w(|w| w.min_bp = 20),
        vec![
            (100, 120),
            (100, 200),
            (100, 300),
            (105, 200),
            (105, 300),
            (120, 200),
            (120, 300),
            (200, 300)
        ]
    );
    // All four at once.
    assert_eq!(
        w(|w| {
            w.snps = 3;
            w.bp = 150;
            w.min_snps = 2;
            w.min_bp = 50;
        }),
        vec![(100, 200), (105, 200)]
    );
    // --min-r2 drops low values and nan (the REF-only site at 300).
    let all_rows = rows(&run(&mut geno(4, LdWindow::default()), &lines, &NO_MASK));
    let kept: Vec<(u64, u64)> = all_rows
        .iter()
        .filter(|r| r[4] != "nan" && r[4].parse::<f64>().unwrap() >= 0.2)
        .map(|r| (r[1].parse().unwrap(), r[2].parse().unwrap()))
        .collect();
    assert_eq!(kept, [(100, 120), (105, 200)]);
    assert_eq!(w(|w| w.min_r2 = Some(0.2)), kept);
    // A minimum of 0 keeps everything, nan included, as vcftools does.
    assert_eq!(w(|w| w.min_r2 = Some(0.0)), all);
}

#[test]
fn ld_genotype_rules() {
    // s1 haploid at the second site, s2 half-missing, s3 masked (DP 3), s4
    // missing; s5 and s6 always called.
    let lines = vec![
        record("1", 1, "G", &["0/1", "0/1", "1/1:3", "0/0", "0/1", "1/1"]),
        record("1", 2, "G", &["1", "./1", "0/1", "./.", "0/0", "1/1"]),
        record("1", 3, "G", &["0/1", "1/1", "0/0:3", "1/1", "0/1", "0/0"]),
    ];
    let out = run(&mut geno(6, LdWindow::default()), &lines, &DP5);
    let n_indv: Vec<String> = rows(&out).iter().map(|r| r[3].clone()).collect();
    // 1-2: s5, s6; 1-3: s1, s2, s4, s5, s6; 2-3: s5, s6.
    assert_eq!(n_indv, ["2", "5", "2"], "{out}");
    // Without the mask s3 counts.
    let out = run(&mut geno(6, LdWindow::default()), &lines, &NO_MASK);
    let n_indv: Vec<String> = rows(&out).iter().map(|r| r[3].clone()).collect();
    assert_eq!(n_indv, ["3", "6", "3"], "{out}");

    // Haplotypes: a haploid genotype is one haplotype, `.|1` one; a masked
    // genotype none.
    let lines = vec![
        record("1", 1, "G", &["0|1", "0|1", "1|1:3", "0|0"]),
        record("1", 2, "G", &["1", ".|1", "0|1", "1|0"]),
        record("1", 3, "G", &["0|1", "1|1", "0|0:3", "1|1"]),
    ];
    let out = run(&mut hap(4, LdWindow::default()), &lines, &DP5);
    let n_chr: Vec<String> = rows(&out).iter().map(|r| r[3].clone()).collect();
    assert_eq!(n_chr, ["4", "6", "4"], "{out}");

    // hap-r2 uses only sites where every genotype is phased, `.|.`
    // included and `./.` not, and the others take no index; geno-r2 uses
    // them all.
    let lines = vec![
        record("1", 1, "G", &["0|1", "1|1", "0|0"]),
        record("1", 2, "G", &["0|1", "0/1", "1|1"]),
        record("1", 3, "G", &["0|0", "./.", "1|1"]),
        record("1", 4, "G", &["1|1", ".|.", "0|1"]),
        record("1", 5, "G", &["1|1", "0|1", "0|1"]),
    ];
    let one = LdWindow {
        snps: 1,
        ..LdWindow::default()
    };
    assert_eq!(
        pairs(&run(&mut hap(3, one), &lines, &NO_MASK)),
        [(1, 4), (4, 5)]
    );
    assert_eq!(
        pairs(&run(&mut geno(3, one), &lines, &NO_MASK)),
        [(1, 2), (2, 3), (3, 4), (4, 5)]
    );
    // A masked unphased genotype still makes the site unphased (vcftools
    // reads the phase before it masks).
    let lines = vec![
        record("1", 1, "G", &["0|1", "1|1", "0|0"]),
        record("1", 2, "G", &["0|1", "0/1:2", "1|1"]),
        record("1", 3, "G", &["1|1", "0|1", "0|1"]),
    ];
    assert_eq!(
        pairs(&run(&mut hap(3, LdWindow::default()), &lines, &DP5)),
        [(1, 3)]
    );
    // A sample without GT is unphased too (vcftools reads it as `./.`): a
    // FORMAT without GT, or a column that stops before it.
    let lines = vec![
        record("1", 1, "G", &["0|1", "1|1", "0|0"]),
        "1\t2\t.\tA\tG\t.\t.\t.\tDP\t9\t9\t9".to_string(),
        "1\t3\t.\tA\tG\t.\t.\t.\tDP:GT\t9:0|1\t9\t9:1|1".to_string(),
        record("1", 4, "G", &["1|1", "0|1", "0|1"]),
    ];
    assert_eq!(pairs(&run(&mut hap(3, one), &lines, &NO_MASK)), [(1, 4)]);

    // A genotype of ploidy above two at a biallelic site stops both, as in
    // vcftools; at a multi-allelic site LD does not read it.
    let poly = vec![
        record("1", 1, "G", &["0/1", "1/1"]),
        record("1", 2, "G", &["0/1/1", "1/1"]),
    ];
    for stat in [
        &mut geno(2, LdWindow::default()) as &mut dyn Stat,
        &mut hap(2, LdWindow::default()),
    ] {
        let e = try_run(stat, &poly, &NO_MASK).unwrap_err();
        assert!(e.contains("polyploid genotypes are not supported"), "{e}");
    }
    let multi = vec![
        record("1", 1, "G", &["0/1", "1/1"]),
        record("1", 2, "G,T", &["0/1/2", "1/1"]),
        record("1", 3, "G", &["0/0", "0/1"]),
    ];
    assert_eq!(
        pairs(&run(&mut geno(2, LdWindow::default()), &multi, &NO_MASK)),
        [(1, 3)]
    );
}

/// 300 phased sites of 8 samples on chromosome 1.
fn many_sites() -> Vec<String> {
    let mut rng = Lcg(3);
    let gts = ["0|0", "0|1", "1|1"];
    (1..=300)
        .map(|i| {
            let s: Vec<&str> = (0..8).map(|_| gts[rng.next(3) as usize]).collect();
            record("1", i * 10, "G", &s)
        })
        .collect()
}

#[test]
fn ld_memory_guard() {
    let lines = many_sites();
    // 8 samples take one 64-sample word per plane: with its index and
    // position a site takes 40 bytes for genotypes, 48 for haplotypes, so
    // 1K holds at most 25 or 21 sites.
    let e = try_run(
        &mut GenoR2::new(samples(8), LdWindow::default(), LdBudget::new(1024)),
        &lines,
        &NO_MASK,
    )
    .unwrap_err();
    for fix in ["--ld-window-bp", "--ld-window", "--max-mem"] {
        assert!(e.contains(fix), "{fix}: {e}");
    }
    // A window keeps the held sites below it.
    let w = LdWindow {
        bp: 200,
        ..LdWindow::default()
    };
    let mut g = GenoR2::new(samples(8), w, LdBudget::new(1024));
    let out = run(&mut g, &lines, &NO_MASK);
    assert_eq!(pairs(&out).len(), 20 * 280 + (20 * 19) / 2);
    assert_eq!(g.peak_held(), 21, "a site and the 20 within 200 bases");
    // The budget is shared: two statistics each holding 21 sites take
    // slots of 40 and 48 bytes, more than 1K together.
    let budget = LdBudget::new(1024);
    let mut g = GenoR2::new(samples(8), w, budget.clone());
    let mut h = HapR2::new(samples(8), w, budget.clone());
    let sink = MemSink::default();
    g.start(&mut sink.clone()).unwrap();
    h.start(&mut sink.clone()).unwrap();
    let mut masked = Vec::new();
    let mut err = None;
    for line in &lines {
        let pos: u64 = line.split('\t').nth(1).unwrap().parse().unwrap();
        let site =
            SiteView::from_text(line.as_bytes(), Interval::new(pos, pos + 1).unwrap()).unwrap();
        let gt = text_genotypes(line.as_bytes()).unwrap();
        let st = site_stats(&gt, None, &NO_MASK, 2, &mut masked).unwrap();
        let mg = MaskedGenotypes {
            inner: &gt,
            masked: &masked,
        };
        if let Err(e) = g
            .site(&site, Some(&mg), Some(&st))
            .and_then(|()| h.site(&site, Some(&mg), Some(&st)))
        {
            err = Some(e);
            break;
        }
    }
    assert!(err.is_some_and(|e| e.contains("--max-mem")));
    // Interchromosomal pairs hold every site: the fix is --max-mem.
    let e = try_run(
        &mut InterchromR2::new(samples(8), false, None, LdBudget::new(1024)),
        &lines,
        &NO_MASK,
    )
    .unwrap_err();
    assert!(e.contains("--max-mem") && !e.contains("--ld-window"), "{e}");
}

#[test]
fn interchrom_order() {
    // Chromosome 1 comes back after 2: its sites pair with 2's, never with
    // each other; the pairs are written in (i, j) input order.
    let a = ["0/0", "0/1", "1/1", "0/1"];
    let b = ["0/1", "0/1", "1/1", "0/0"];
    let lines = vec![
        record("1", 10, "G", &a),
        record("2", 5, "G", &b),
        record("1", 20, "G,T", &["0/2", "0/1", "1/1", "0/0"]),
        record("3", 7, "G", &a),
        record("1", 30, "G", &b),
        record("2", 1, "G", &a),
    ];
    let out = run(
        &mut InterchromR2::new(samples(4), false, None, budget()),
        &lines,
        &NO_MASK,
    );
    assert!(
        out.starts_with("CHR1\tPOS1\tCHR2\tPOS2\tN_INDV\tR^2\n"),
        "{out}"
    );
    let got: Vec<String> = rows(&out)
        .iter()
        .map(|r| format!("{}:{}-{}:{}", r[0], r[1], r[2], r[3]))
        .collect();
    assert_eq!(
        got,
        [
            "1:10-2:5", "1:10-3:7", "1:10-2:1", "2:5-3:7", "2:5-1:30", "3:7-1:30", "3:7-2:1",
            "1:30-2:1"
        ]
    );
    // The haplotype form: N_CHR and R^2 only, phased sites only.
    let lines = vec![
        record("1", 10, "G", &["0|0", "0|1"]),
        record("1", 20, "G", &["0/0", "0|1"]),
        record("2", 5, "G", &["0|1", "1|1"]),
    ];
    let out = run(
        &mut InterchromR2::new(samples(2), true, None, budget()),
        &lines,
        &NO_MASK,
    );
    assert_eq!(rows(&out).len(), 1, "{out}");
    assert!(
        out.starts_with("CHR1\tPOS1\tCHR2\tPOS2\tN_CHR\tR^2\n1\t10\t2\t5\t4\t"),
        "{out}"
    );
    // --min-r2 applies.
    let out = run(
        &mut InterchromR2::new(samples(2), true, Some(2.0), budget()),
        &lines,
        &NO_MASK,
    );
    assert_eq!(rows(&out).len(), 0, "{out}");
}

#[test]
fn ld_windows_need_sorted_input() {
    let s = ["0/0", "0/1"];
    for (lines, why) in [
        (
            vec![
                record("1", 50, "G", &s),
                record("2", 5, "G", &s),
                record("1", 60, "G", &s),
            ],
            "the chromosome reappears at 1:60",
        ),
        (
            vec![record("1", 50, "G", &s), record("1", 40, "G", &s)],
            "a position goes backwards at 1:40",
        ),
    ] {
        for (stat, flag) in [
            (
                &mut geno(2, LdWindow::default()) as &mut dyn Stat,
                "--geno-r2",
            ),
            (&mut hap(2, LdWindow::default()), "--hap-r2"),
        ] {
            let e = try_run(stat, &lines, &NO_MASK).unwrap_err();
            assert!(
                e.contains(&format!("{flag} needs input sorted")) && e.contains(why),
                "{e}"
            );
            assert!(e.contains("gloci sort"), "{e}");
        }
        // The interchromosomal pairs take any order.
        run(
            &mut InterchromR2::new(samples(2), false, None, budget()),
            &lines,
            &NO_MASK,
        );
    }
}

/// A statistic and the bytes it has allocated for the sites it holds.
trait Holds: Stat {
    fn held(&self) -> u64;
}

impl Holds for GenoR2 {
    fn held(&self) -> u64 {
        self.held_bytes()
    }
}

impl Holds for HapR2 {
    fn held(&self) -> u64 {
        self.held_bytes()
    }
}

impl Holds for InterchromR2 {
    fn held(&self) -> u64 {
        self.held_bytes()
    }
}

/// Feed `lines` to every statistic in turn, site by site, until one fails;
/// after each site the bytes they hold together, by what is allocated, must
/// be within `max`. The number of sites fed.
fn within_budget(stats: &mut [&mut dyn Holds], lines: &[String], max: u64) -> usize {
    let sink = MemSink::default();
    for s in stats.iter_mut() {
        s.start(&mut sink.clone()).unwrap();
    }
    let mut masked = Vec::new();
    for (n, line) in lines.iter().enumerate() {
        let pos: u64 = line.split('\t').nth(1).unwrap().parse().unwrap();
        let site =
            SiteView::from_text(line.as_bytes(), Interval::new(pos, pos + 1).unwrap()).unwrap();
        let gt = text_genotypes(line.as_bytes()).unwrap();
        let st = site_stats(&gt, None, &NO_MASK, 2, &mut masked).unwrap();
        let mg = MaskedGenotypes {
            inner: &gt,
            masked: &masked,
        };
        for s in stats.iter_mut() {
            if s.site(&site, Some(&mg), Some(&st)).is_err() {
                return n;
            }
        }
        let now: u64 = stats.iter().map(|s| s.held()).sum();
        assert!(
            now <= max,
            "site {n}: {now} bytes held, past the budget of {max}"
        );
    }
    lines.len()
}

/// `n` sites of 8 samples on `chrom`, phased or not.
fn sites_on(chrom: &str, n: u64, phased: bool, seed: u64) -> Vec<String> {
    let mut rng = Lcg(seed);
    let gts: [&str; 3] = if phased {
        ["0|0", "0|1", "1|1"]
    } else {
        ["0/0", "0/1", "1/1"]
    };
    (1..=n)
        .map(|i| {
            let s: Vec<&str> = (0..8).map(|_| gts[rng.next(3) as usize]).collect();
            record(chrom, i * 10, "G", &s)
        })
        .collect()
}

/// The budget bounds what the LD statistics hold together, by the bytes
/// allocated: sites one statistic has let go of are not free for another
/// while it keeps their memory, and the interchromosomal store does not
/// grow past it.
#[test]
fn ld_budget_bounds_the_bytes_held() {
    // A phased chromosome both hold whole, then an unphased one only
    // --geno-r2 holds.
    let mut lines = sites_on("1", 30, true, 5);
    lines.extend(sites_on("2", 200, false, 6));
    let max = 4096;
    let budget = LdBudget::new(max);
    let mut g = GenoR2::new(samples(8), LdWindow::default(), budget.clone());
    let mut h = HapR2::new(samples(8), LdWindow::default(), budget.clone());
    let fed = within_budget(&mut [&mut g, &mut h], &lines, max);
    assert!(
        fed > 30 && fed < lines.len(),
        "{fed}: the budget stops chromosome 2"
    );
    // --geno-r2 beside the interchromosomal pairs.
    let budget = LdBudget::new(max);
    let w = LdWindow {
        snps: 20,
        ..LdWindow::default()
    };
    let mut g = GenoR2::new(samples(8), w, budget.clone());
    let mut x = InterchromR2::new(samples(8), false, None, budget.clone());
    let fed = within_budget(&mut [&mut g, &mut x], &lines, max);
    assert!(
        fed > 30 && fed < lines.len(),
        "{fed}: the budget stops the pairs"
    );
}
