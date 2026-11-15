//! Relatedness (decision D162) on hand-made sites: Manichaikul's φ and
//! Yang's unadjusted Ajk against oracles written from the published
//! formulas, the genotype rules, the 64-site blocks against a site-by-site
//! sum, the rows and the memory the matrices hold. The vcftools 0.1.16
//! comparison is `crates/gloci/tests/vcf_popgen.rs`.
// The oracles are matrices indexed by sample, beside rows named by index.
#![allow(clippy::needless_range_loop)]

use gloci_core::Interval;
use gloci_core::cfmt::write_g6;
use gloci_vcf::popgen::{RelatednessManichaikul, RelatednessYang};
use gloci_vcf::{
    MaskRules, MaskedGenotypes, SiteView, Stat, StatSink, StatWrite, site_stats, text_genotypes,
};
use std::cell::RefCell;
use std::collections::HashMap;
use std::io::{self, Write};
use std::num::NonZeroUsize;
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
fn record(pos: u64, alts: &str, samples: &[&str]) -> String {
    format!(
        "1\t{pos}\t.\tA\t{alts}\t.\t.\t.\tGT:DP\t{}",
        samples.join("\t")
    )
}

/// Feed `lines` to `stat` under `mask`; the file it writes, or the first
/// error.
fn try_run(stat: &mut dyn Stat, lines: &[String], mask: &MaskRules) -> Result<String, String> {
    try_run_kept(stat, lines, mask, None)
}

/// [`try_run`] with the site statistics over the samples `keep` selects,
/// as the selection stage computes them.
fn try_run_kept(
    stat: &mut dyn Stat,
    lines: &[String],
    mask: &MaskRules,
    keep: Option<&[bool]>,
) -> Result<String, String> {
    let sink = MemSink::default();
    stat.start(&mut sink.clone()).map_err(|e| e.to_string())?;
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
        let st = site_stats(&g, keep, mask, n_alleles, &mut masked).unwrap();
        let mg = MaskedGenotypes {
            inner: &g,
            masked: &masked,
        };
        stat.site(&site, Some(&mg), Some(&st))?;
    }
    stat.finish().map_err(|e| e.to_string())?;
    let files = sink.0.take();
    assert_eq!(files.len(), 1, "one file: {:?}", files.keys());
    let (ext, text) = files.into_iter().next().unwrap();
    assert_eq!(ext, stat.ext());
    Ok(String::from_utf8(text).unwrap())
}

fn run(stat: &mut dyn Stat, lines: &[String], mask: &MaskRules) -> String {
    try_run(stat, lines, mask).unwrap()
}

fn threads(n: usize) -> NonZeroUsize {
    NonZeroUsize::new(n).unwrap()
}

/// Samples 0..n named s0, s1, ...
fn named(n: usize) -> Vec<(usize, String)> {
    (0..n).map(|i| (i, format!("s{i}"))).collect()
}

fn yang(n: usize) -> RelatednessYang {
    RelatednessYang::new(named(n), threads(1)).unwrap()
}

fn mani(n: usize) -> RelatednessManichaikul {
    RelatednessManichaikul::new(named(n), threads(1)).unwrap()
}

/// The data rows of a file, without its header.
fn rows(text: &str) -> Vec<Vec<String>> {
    text.lines()
        .skip(1)
        .map(|l| l.split('\t').map(str::to_string).collect())
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

/// `n` samples at `count` biallelic sites: called genotypes at a frequency
/// that changes from site to site, missing, half-called (`./1`, `0/.`) and
/// low-depth (DP 3) ones, and now and then a site where every call is REF
/// or ALT.
fn generated(seed: u64, n: usize, count: usize) -> Vec<Vec<String>> {
    let mut rng = Lcg(seed);
    let mut sites = Vec::new();
    for _ in 0..count {
        let p = rng.next(1000);
        let fixed = rng.next(12);
        let site = (0..n)
            .map(|_| {
                let r = rng.next(100);
                let a = |rng: &mut Lcg| u32::from(rng.next(1000) < p);
                let g = match r {
                    0..=5 => "./.".to_string(),
                    6..=8 => "./1".to_string(),
                    9..=10 => "0/.".to_string(),
                    _ => {
                        let (x, y) = (a(&mut rng), a(&mut rng));
                        match fixed {
                            0 => "0/0".to_string(),
                            1 => "1/1".to_string(),
                            _ => format!("{x}/{y}"),
                        }
                    }
                };
                let dp = if rng.next(10) == 0 { 3 } else { 9 };
                format!("{g}:{dp}")
            })
            .collect();
        sites.push(site);
    }
    sites
}

fn lines_of(sites: &[Vec<String>]) -> Vec<String> {
    sites
        .iter()
        .enumerate()
        .map(|(i, s)| {
            let s: Vec<&str> = s.iter().map(String::as_str).collect();
            record(10 * (i as u64 + 1), "G", &s)
        })
        .collect()
}

/// A sample `GT:DP` under `mask`: the ALT allele count of a fully called
/// diploid genotype, `None` for a missing, half-called or masked one.
fn dosage(g: &str, mask: &MaskRules) -> Option<u32> {
    let (gt, dp) = g.split_once(':').unwrap_or((g, "9"));
    if mask.min_dp.is_some_and(|m| dp.parse::<u32>().unwrap() < m) {
        return None;
    }
    let (x, y) = gt.split_once(['/', '|'])?;
    Some(x.parse::<u32>().ok()? + y.parse::<u32>().ok()?)
}

/// The ALT allele frequency among the called alleles of the unmasked
/// samples, half-called genotypes included.
fn alt_frequency(site: &[String], mask: &MaskRules) -> f64 {
    let (mut alt, mut called) = (0u32, 0u32);
    for g in site {
        let (gt, dp) = g.split_once(':').unwrap_or((g, "9"));
        if mask.min_dp.is_some_and(|m| dp.parse::<u32>().unwrap() < m) {
            continue;
        }
        for a in gt.split(['/', '|']) {
            if let Ok(a) = a.parse::<u32>() {
                called += 1;
                alt += a;
            }
        }
    }
    f64::from(alt) / f64::from(called)
}

/// Manichaikul et al. (2010) from the genotypes: per pair the sites where
/// both are heterozygous, where they are opposite homozygotes, each one's
/// heterozygous sites, and φ = (N_AaAa − 2 N_AAaa) / (N_Aa,i + N_Aa,j).
fn kinship_oracle(sites: &[Vec<String>], mask: &MaskRules) -> Vec<Vec<(u32, u32, u32, f64)>> {
    let n = sites[0].len();
    let calls: Vec<Vec<Option<u32>>> = sites
        .iter()
        .map(|s| s.iter().map(|g| dosage(g, mask)).collect())
        .collect();
    let het = |i: usize| calls.iter().filter(|s| s[i] == Some(1)).count() as u32;
    (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    let both_het = calls
                        .iter()
                        .filter(|s| s[i] == Some(1) && s[j] == Some(1))
                        .count() as u32;
                    let opposite = calls
                        .iter()
                        .filter(|s| matches!((s[i], s[j]), (Some(0), Some(2)) | (Some(2), Some(0))))
                        .count() as u32;
                    let phi = (f64::from(both_het) - 2.0 * f64::from(opposite))
                        / f64::from(het(i) + het(j));
                    (both_het, opposite, het(i), phi)
                })
                .collect()
        })
        .collect()
}

/// Yang et al. (2010), equation 6, unadjusted: over the sites whose ALT
/// frequency p is strictly between 0 and 1 and the samples called at both,
/// the mean of (x_j − 2p)(x_k − 2p) / (2p(1 − p)), and for a sample with
/// itself 1 + the mean of (x² − (1 + 2p)x + 2p²) / (2p(1 − p)).
fn yang_oracle(sites: &[Vec<String>], mask: &MaskRules) -> Vec<Vec<f64>> {
    let n = sites[0].len();
    let mut sum = vec![vec![0.0f64; n]; n];
    let mut count = vec![vec![0u32; n]; n];
    for site in sites {
        let p = alt_frequency(site, mask);
        if !(p > 0.0 && p < 1.0) {
            continue;
        }
        let x: Vec<Option<f64>> = site
            .iter()
            .map(|g| dosage(g, mask).map(f64::from))
            .collect();
        let v = 2.0 * p * (1.0 - p);
        for j in 0..n {
            let Some(xj) = x[j] else { continue };
            sum[j][j] += (xj * xj - (1.0 + 2.0 * p) * xj + 2.0 * p * p) / v;
            count[j][j] += 1;
            for k in j + 1..n {
                let Some(xk) = x[k] else { continue };
                sum[j][k] += (xj - 2.0 * p) * (xk - 2.0 * p) / v;
                count[j][k] += 1;
            }
        }
    }
    (0..n)
        .map(|j| {
            (0..n)
                .map(|k| {
                    if k < j {
                        return f64::NAN;
                    }
                    let mean = sum[j][k] / f64::from(count[j][k]);
                    if j == k { 1.0 + mean } else { mean }
                })
                .collect()
        })
        .collect()
}

fn g6(x: f64) -> String {
    let mut b = Vec::new();
    write_g6(&mut b, x);
    String::from_utf8(b).unwrap()
}

fn assert_close(got: f64, want: f64, what: &str) {
    if want.is_nan() {
        assert!(got.is_nan(), "{what}: {got}, formula nan");
        return;
    }
    assert!(
        (got - want).abs() <= 1e-12 * want.abs().max(1.0),
        "{what}: {got}, formula {want}"
    );
}

#[test]
fn relatedness2_oracle() {
    let n = 23;
    let sites = generated(5, n, 150);
    for (mask, tag) in [(NO_MASK, "unmasked"), (DP5, "--min-dp 5")] {
        let mut stat = mani(n);
        let out = run(&mut stat, &lines_of(&sites), &mask);
        assert!(
            out.starts_with("INDV1\tINDV2\tN_AaAa\tN_AAaa\tN1_Aa\tN2_Aa\tRELATEDNESS_PHI\n"),
            "{out}"
        );
        let want = kinship_oracle(&sites, &mask);
        let rows = rows(&out);
        assert_eq!(rows.len(), n * n, "{tag}: every ordered pair");
        let mut k = 0;
        for i in 0..n {
            for j in 0..n {
                let r = &rows[k];
                k += 1;
                let (aaaa, aa_aa, het, phi) = want[i][j];
                assert_eq!(r[..2], [format!("s{i}"), format!("s{j}")], "{tag}");
                assert_eq!(
                    r[2..6],
                    [
                        aaaa.to_string(),
                        aa_aa.to_string(),
                        het.to_string(),
                        want[j][j].2.to_string()
                    ],
                    "{tag}: counts of s{i}, s{j}"
                );
                assert_close(stat.phi(i, j), phi, &format!("{tag}: phi of s{i}, s{j}"));
                assert_eq!(r[6], g6(stat.phi(i, j)), "{tag}: printed phi");
            }
        }
        // the sites hold rare and common alleles alike
        assert!(rows.iter().any(|r| r[6].starts_with('-')), "{tag}: phi < 0");
    }
}

#[test]
fn relatedness_oracle() {
    let n = 23;
    let sites = generated(9, n, 150);
    for (mask, tag) in [(NO_MASK, "unmasked"), (DP5, "--min-dp 5")] {
        let mut stat = yang(n);
        let out = run(&mut stat, &lines_of(&sites), &mask);
        assert!(out.starts_with("INDV1\tINDV2\tRELATEDNESS_AJK\n"), "{out}");
        let want = yang_oracle(&sites, &mask);
        let rows = rows(&out);
        assert_eq!(rows.len(), n * (n + 1) / 2, "{tag}: i <= j");
        let mut k = 0;
        for i in 0..n {
            for j in i..n {
                let r = &rows[k];
                k += 1;
                assert_eq!(r[..2], [format!("s{i}"), format!("s{j}")], "{tag}");
                assert_close(stat.ajk(i, j), want[i][j], &format!("{tag}: s{i}, s{j}"));
                assert_eq!(r[2], g6(stat.ajk(i, j)), "{tag}: printed");
            }
        }
    }
}

/// Four samples; the first site is every genotype rule at once.
#[test]
fn relatedness_genotype_rules() {
    let base = [
        record(10, "G", &["0/1:9", "1/1:9", "0/0:9", "0/1:9"]),
        record(20, "G", &["1/1:9", "0/1:9", "0/0:9", "0/0:9"]),
        record(30, "G", &["0/0:9", "0/1:9", "1/1:9", "0/1:9"]),
    ];
    let with = |first: &[&str]| -> Vec<String> {
        let mut v = base.to_vec();
        v.push(record(40, "G", first));
        v
    };
    let yang_of = |lines: &[String], mask: &MaskRules| -> RelatednessYang {
        let mut s = yang(4);
        run(&mut s, lines, mask);
        s
    };
    let mani_of = |lines: &[String], mask: &MaskRules| -> RelatednessManichaikul {
        let mut s = mani(4);
        run(&mut s, lines, mask);
        s
    };
    let same_yang = |a: &RelatednessYang, b: &RelatednessYang, what: &str| {
        for i in 0..4 {
            for j in i..4 {
                let (x, y) = (a.ajk(i, j), b.ajk(i, j));
                assert!(
                    x.to_bits() == y.to_bits() || (x.is_nan() && y.is_nan()),
                    "{what}: s{i}, s{j}: {x} and {y}"
                );
            }
        }
    };
    let without = yang_of(&base, &NO_MASK);

    // A half-called genotype is missing: ./1 is not 0 (vcftools' reading)
    // and its called allele still counts in the frequency, 4/7.
    let half = yang_of(&with(&["./1:9", "1/1:9", "0/0:9", "0/1:9"]), &NO_MASK);
    let sites: Vec<Vec<String>> = with(&["./1:9", "1/1:9", "0/0:9", "0/1:9"])
        .iter()
        .map(|l| l.split('\t').skip(9).map(str::to_string).collect())
        .collect();
    assert_eq!(alt_frequency(&sites[3], &NO_MASK), 4.0 / 7.0);
    let want = yang_oracle(&sites, &NO_MASK);
    for j in 0..4 {
        assert_close(half.ajk(0, j), want[0][j], &format!("./1 with s{j}"));
    }
    // ... and in --relatedness2
    let m_half = mani_of(&with(&["./1:9", "0/1:9", "0/0:9", "1/1:9"]), &NO_MASK);
    let m_none = mani_of(&with(&["./.:9", "0/1:9", "0/0:9", "1/1:9"]), &NO_MASK);
    assert_eq!(m_half.counts(0, 1), m_none.counts(0, 1));

    // A masked genotype is missing in both statistics, the frequency
    // included: the same as ./. there.
    let masked_first = with(&["0/1:3", "1/1:9", "0/0:9", "1/1:9"]);
    let missing_first = with(&["./.:9", "1/1:9", "0/0:9", "1/1:9"]);
    same_yang(
        &yang_of(&masked_first, &DP5),
        &yang_of(&missing_first, &NO_MASK),
        "masked",
    );
    let (mm, mn) = (
        mani_of(&masked_first, &DP5),
        mani_of(&missing_first, &NO_MASK),
    );
    let masked_hom = with(&["0/0:3", "1/1:9", "0/0:9", "1/1:9"]);
    let mh = mani_of(&masked_hom, &DP5);
    for i in 0..4 {
        for j in 0..4 {
            assert_eq!(mm.counts(i, j), mn.counts(i, j), "masked het s{i}, s{j}");
            assert_eq!(mh.counts(i, j), mn.counts(i, j), "masked hom s{i}, s{j}");
        }
    }
    // unmasked, that opposite homozygote counts
    assert_eq!(
        mani_of(&masked_hom, &NO_MASK).counts(0, 1).1,
        mn.counts(0, 1).1 + 1
    );

    // A haploid genotype drops the whole site, unless masked; a lone '.'
    // too. Multi-allelic and ALT '.' sites are not used; neither are sites
    // where every call is REF or ALT (--relatedness).
    for (first, mask, what) in [
        (&["0/1:9", "1:9", "0/0:9", "0/1:9"][..], NO_MASK, "haploid"),
        (&["0/1:9", ".:9", "0/0:9", "0/1:9"], NO_MASK, "lone ."),
        (&["0/0:9", "0/0:9", "0/0:9", "0/0:9"], NO_MASK, "REF only"),
        (&["1/1:9", "1/1:9", "./.:9", "1/1:9"], NO_MASK, "ALT only"),
        (
            &["1/1:9", "1/1:9", "1/1:9", "0/0:3"],
            DP5,
            "ALT only after masking",
        ),
    ] {
        same_yang(&yang_of(&with(first), &mask), &without, what);
    }
    let masked_haploid = yang_of(&with(&["0/1:9", "1:3", "0/0:9", "0/1:9"]), &DP5);
    let missing_one = yang_of(&with(&["0/1:9", "./.:9", "0/0:9", "0/1:9"]), &NO_MASK);
    same_yang(&masked_haploid, &missing_one, "masked haploid");
    let mut v = base.to_vec();
    v.push(record(40, "G,T", &["0/2:9", "1/1:9", "0/0:9", "2/2:9"]));
    v.push(record(50, ".", &["0/0:9", "0/0:9", "0/0:9", "0/0:9"]));
    same_yang(&yang_of(&v, &NO_MASK), &without, "multi-allelic, ALT .");
    let m_without = mani_of(&base, &NO_MASK);
    for (lines, what) in [
        (v.clone(), "multi-allelic"),
        (with(&["0/1:9", "1:9", "0/0:9", "1/1:9"]), "haploid"),
    ] {
        let m = mani_of(&lines, &NO_MASK);
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(m.counts(i, j), m_without.counts(i, j), "{what}");
            }
        }
    }

    // A sample without GT is missing and keeps the site.
    let mut nogt = base.to_vec();
    nogt.push("1\t40\t.\tA\tG\t.\t.\t.\tDP:GT\t9:0/1\t9\t9:0/0\t9:0/1".to_string());
    same_yang(
        &yang_of(&nogt, &NO_MASK),
        &yang_of(&with(&["0/1:9", "./.:9", "0/0:9", "0/1:9"]), &NO_MASK),
        "no GT",
    );

    // A genotype of ploidy above two stops both at a biallelic site, with
    // the hint to leave the sample out; not at a multi-allelic one, nor
    // when masked.
    for stat in [&mut yang(4) as &mut dyn Stat, &mut mani(4)] {
        let e = try_run(
            stat,
            &with(&["0/1/1:9", "1/1:9", "0/0:9", "0/1:9"]),
            &NO_MASK,
        )
        .unwrap_err();
        assert!(e.contains("polyploid"), "{e}");
        assert!(e.contains("--exclude-samples"), "{e}");
    }
    let mut v = base.to_vec();
    v.push(record(40, "G,T", &["0/1/2:9", "1/1:9", "0/0:9", "0/1:9"]));
    run(&mut yang(4), &v, &NO_MASK);
    run(&mut mani(4), &v, &NO_MASK);
    run(
        &mut yang(4),
        &with(&["0/1/1:3", "1/1:9", "0/0:9", "0/1:9"]),
        &DP5,
    );
}

/// The site-by-site sums, in vcftools' order of operations: what the
/// 64-site blocks must give bit for bit.
fn site_by_site(sites: &[Vec<String>]) -> (Vec<Vec<f64>>, Vec<Vec<u32>>) {
    let n = sites[0].len();
    let mut sum = vec![vec![0.0f64; n]; n];
    let mut count = vec![vec![0u32; n]; n];
    for site in sites {
        let freq = alt_frequency(site, &NO_MASK);
        if freq <= f64::EPSILON || freq >= 1.0 - f64::EPSILON || freq.is_nan() {
            continue;
        }
        let x: Vec<Option<f64>> = site
            .iter()
            .map(|g| dosage(g, &NO_MASK).map(f64::from))
            .collect();
        let div = 1.0 / (2.0 * freq * (1.0 - freq));
        for i in 0..n {
            let Some(xi) = x[i] else { continue };
            sum[i][i] += (xi * xi - (1.0 + 2.0 * freq) * xi + 2.0 * freq * freq) * div;
            count[i][i] += 1;
            for j in i + 1..n {
                let Some(xj) = x[j] else { continue };
                sum[i][j] += (xi - 2.0 * freq) * (xj - 2.0 * freq) * div;
                count[i][j] += 1;
            }
        }
    }
    (sum, count)
}

#[test]
fn relatedness_blocks() {
    // 300 samples: enough work per block for every thread to take rows.
    let n = 300;
    let sites = generated(13, n, 130);
    let (sum, count) = site_by_site(&sites);
    let lines = lines_of(&sites);
    let mut first: Option<String> = None;
    for t in [1, 2, 3, 8] {
        let mut stat = RelatednessYang::new(named(n), threads(t)).unwrap();
        let out = run(&mut stat, &lines, &NO_MASK);
        for i in 0..n {
            for j in i..n {
                let mean = sum[i][j] / f64::from(count[i][j]);
                let want = if i == j { 1.0 + mean } else { mean };
                let got = stat.ajk(i, j);
                assert!(
                    got.to_bits() == want.to_bits() || (got.is_nan() && want.is_nan()),
                    "-t {t}: s{i}, s{j}: {got}, site by site {want}"
                );
            }
        }
        let mut m = RelatednessManichaikul::new(named(n), threads(t)).unwrap();
        let out2 = run(&mut m, &lines, &NO_MASK);
        let both = format!("{out}{out2}");
        match &first {
            None => first = Some(both),
            Some(f) => assert!(*f == both, "-t {t} wrote other files than -t 1"),
        }
    }
    // the 130 sites hold two full blocks and a partial one, and some sites
    // are not used, so used sites and blocks do not line up
    assert!(count[0][0] > 64 && count[0][0] < 130, "{}", count[0][0]);
}

/// `--samples` picks the rows: only the selected samples, in header order,
/// each read from its own column.
#[test]
fn relatedness_rows_follow_the_selection() {
    let n = 6;
    let sites = generated(21, n, 40);
    let lines = lines_of(&sites);
    let pick = vec![
        (1, "b".to_string()),
        (3, "d".to_string()),
        (4, "e".to_string()),
    ];
    let picked: Vec<Vec<String>> = sites
        .iter()
        .map(|s| pick.iter().map(|(i, _)| s[*i].clone()).collect())
        .collect();
    let keep = [false, true, false, true, true, false];
    let mut y = RelatednessYang::new(pick.clone(), threads(2)).unwrap();
    let out = try_run_kept(&mut y, &lines, &NO_MASK, Some(&keep)).unwrap();
    let names: Vec<(String, String)> = rows(&out)
        .into_iter()
        .map(|r| (r[0].clone(), r[1].clone()))
        .collect();
    let pairs = |v: &[(&str, &str)]| -> Vec<(String, String)> {
        v.iter()
            .map(|(a, b)| (a.to_string(), b.to_string()))
            .collect()
    };
    assert_eq!(
        names,
        pairs(&[
            ("b", "b"),
            ("b", "d"),
            ("b", "e"),
            ("d", "d"),
            ("d", "e"),
            ("e", "e")
        ])
    );
    let want = yang_oracle(&picked, &NO_MASK);
    for i in 0..3 {
        for j in i..3 {
            assert_close(y.ajk(i, j), want[i][j], &format!("{i}, {j}"));
        }
    }
    let mut m = RelatednessManichaikul::new(pick, threads(2)).unwrap();
    let out = try_run_kept(&mut m, &lines, &NO_MASK, Some(&keep)).unwrap();
    let names: Vec<(String, String)> = rows(&out)
        .into_iter()
        .map(|r| (r[0].clone(), r[1].clone()))
        .collect();
    assert_eq!(
        names,
        pairs(&[
            ("b", "b"),
            ("b", "d"),
            ("b", "e"),
            ("d", "b"),
            ("d", "d"),
            ("d", "e"),
            ("e", "b"),
            ("e", "d"),
            ("e", "e")
        ])
    );
    let want = kinship_oracle(&picked, &NO_MASK);
    for i in 0..3 {
        for j in 0..3 {
            let (aaaa, aa_aa, _, phi) = want[i][j];
            assert_eq!(m.counts(i, j), (aaaa, aa_aa), "{i}, {j}");
            assert_close(m.phi(i, j), phi, &format!("{i}, {j}"));
        }
    }
}

/// No site used: every estimate is nan, and every row is written.
#[test]
fn relatedness_without_sites() {
    let lines = vec![record(10, "G", &["0/0:9", "0/0:9", "0/0:9"])];
    let out = run(&mut yang(3), &lines, &NO_MASK);
    assert_eq!(rows(&out).len(), 6);
    assert!(rows(&out).iter().all(|r| r[2] == "nan"), "{out}");
    let out = run(&mut mani(3), &[], &NO_MASK);
    assert_eq!(rows(&out).len(), 9);
    assert!(
        rows(&out)
            .iter()
            .all(|r| r[2..] == ["0", "0", "0", "0", "nan"])
    );
    // opposite homozygotes and no heterozygote: -inf
    let lines = vec![record(10, "G", &["0/0:9", "1/1:9"])];
    let out = run(&mut mani(2), &lines, &NO_MASK);
    assert_eq!(rows(&out)[1][2..], ["0", "1", "0", "0", "-inf"]);
}

/// The matrices take 12 bytes (--relatedness) or 8 (--relatedness2) per
/// pair of samples with itself, and one block of 64 sites beside them
/// (decision D162).
#[test]
fn relatedness_memory_bound() {
    for n in [1usize, 7, 64, 65, 300] {
        let pairs = (n * (n + 1) / 2) as u64;
        let sites = generated(3, n, 70);
        let mut y = yang(n);
        run(&mut y, &lines_of(&sites), &NO_MASK);
        let (matrix, block) = y.held_bytes();
        assert_eq!(matrix, 12 * pairs, "{n} samples");
        // the centred genotypes of 64 sites, the called samples of each
        // and the site constants
        assert!(
            block <= 64 * 8 * n as u64 + 8 * n as u64 + 64 * 32,
            "{n}: {block}"
        );
        let mut m = mani(n);
        run(&mut m, &lines_of(&sites), &NO_MASK);
        let (matrix, block) = m.held_bytes();
        assert_eq!(matrix, 8 * pairs, "{n} samples");
        // three 64-site masks per sample
        assert!(block <= 3 * 8 * n as u64, "{n}: {block}");
    }
}
