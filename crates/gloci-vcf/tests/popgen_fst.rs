//! Weir and Cockerham's Fst (decision D160) on hand-made sites: the oracle
//! from the published formula, the population files, all-missing
//! populations, masking and the streamed windows. The vcftools 0.1.16
//! comparison is `crates/gloci/tests/vcf_popgen.rs`.

use gloci_core::Interval;
use gloci_core::cfmt::write_g6;
use gloci_vcf::popgen::{Populations, WeirFst};
use gloci_vcf::{
    MaskRules, MaskedGenotypes, SiteView, Stat, StatSink, StatWrite, site_stats, text_genotypes,
};
use std::cell::RefCell;
use std::collections::{BTreeMap, HashMap};
use std::io::{self, Write};
use std::path::PathBuf;
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

/// A record `CHROM POS . A ALTS . . . GT:DP SAMPLES...`; each sample is
/// `GT` or `GT:DP`.
fn record(chrom: &str, pos: u64, alts: &str, samples: &[&str]) -> String {
    format!(
        "{chrom}\t{pos}\t.\tA\t{alts}\t.\t.\t.\tGT:DP\t{}",
        samples.join("\t")
    )
}

/// Populations of consecutive samples: `sizes[i]` samples each, in order.
fn pops(sizes: &[usize]) -> Populations {
    let mut at = 0;
    let mut members = Vec::new();
    for &n in sizes {
        members.push((at..at + n).collect());
        at += n;
    }
    Populations {
        names: (1..=sizes.len()).map(|i| format!("pop{i}.txt")).collect(),
        members,
    }
}

/// Feed `lines` to `stat` under `mask`; the files it writes, or the first
/// error.
fn try_run(
    stat: &mut dyn Stat,
    lines: &[String],
    mask: &MaskRules,
) -> Result<HashMap<String, String>, String> {
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
    stat.finish().unwrap();
    Ok(sink
        .0
        .take()
        .into_iter()
        .map(|(k, v)| (k, String::from_utf8(v).unwrap()))
        .collect())
}

fn run(stat: &mut dyn Stat, lines: &[String], mask: &MaskRules) -> HashMap<String, String> {
    try_run(stat, lines, mask).unwrap()
}

fn g6(x: f64) -> String {
    let mut b = Vec::new();
    write_g6(&mut b, x);
    String::from_utf8(b).unwrap()
}

/// Weir and Cockerham (1984), "Estimating F-statistics for the analysis of
/// population structure", equations 2 to 4 and 10: for each allele A, over
/// the r populations with n_i complete genotypes,
/// n̄ = Σ n_i / r, n_c = (r n̄ − Σ n_i² / (r n̄)) / (r − 1),
/// p̄ = Σ n_i p_i / (r n̄), s² = Σ n_i (p_i − p̄)² / ((r − 1) n̄),
/// h̄ = Σ n_i h_i / (r n̄) (h_i the share of individuals heterozygous for A),
/// a = n̄/n_c · (s² − (p̄(1 − p̄) − (r − 1)/r · s² − h̄/4) / (n̄ − 1)),
/// b = n̄/(n̄ − 1) · (p̄(1 − p̄) − (r − 1)/r · s² − (2n̄ − 1)/(4n̄) · h̄),
/// c = h̄/2, and θ = Σ_A a / Σ_A (a + b + c). Genotypes with a missing
/// allele are not complete.
fn weir_cockerham_1984(pops: &[Vec<&str>], n_alleles: usize) -> f64 {
    let r = pops.len() as f64;
    // Per population: complete genotypes as allele pairs.
    let called: Vec<Vec<(usize, usize)>> = pops
        .iter()
        .map(|p| {
            p.iter()
                .filter_map(|g| {
                    let (x, y) = g.split_once('/')?;
                    Some((x.parse().ok()?, y.parse().ok()?))
                })
                .collect()
        })
        .collect();
    let n: Vec<f64> = called.iter().map(|c| c.len() as f64).collect();
    let n_tot: f64 = n.iter().sum();
    let n_bar = n_tot / r;
    let n_c = (r * n_bar - n.iter().map(|x| x * x).sum::<f64>() / (r * n_bar)) / (r - 1.0);
    let (mut sum_a, mut sum_abc) = (0.0, 0.0);
    for allele in 0..n_alleles {
        let p: Vec<f64> = called
            .iter()
            .map(|c| {
                let copies: usize = c
                    .iter()
                    .map(|&(x, y)| usize::from(x == allele) + usize::from(y == allele))
                    .sum();
                copies as f64 / (2.0 * c.len() as f64)
            })
            .collect();
        let h: Vec<f64> = called
            .iter()
            .map(|c| {
                let het = c
                    .iter()
                    .filter(|&&(x, y)| x != y && (x == allele || y == allele))
                    .count();
                het as f64 / c.len() as f64
            })
            .collect();
        let p_bar = (0..pops.len()).map(|i| n[i] * p[i]).sum::<f64>() / (r * n_bar);
        let s2 = (0..pops.len())
            .map(|i| n[i] * (p[i] - p_bar).powi(2))
            .sum::<f64>()
            / ((r - 1.0) * n_bar);
        let h_bar = (0..pops.len()).map(|i| n[i] * h[i]).sum::<f64>() / (r * n_bar);
        let pq = p_bar * (1.0 - p_bar);
        let a = n_bar / n_c * (s2 - (pq - (r - 1.0) / r * s2 - h_bar / 4.0) / (n_bar - 1.0));
        let b = n_bar / (n_bar - 1.0)
            * (pq - (r - 1.0) / r * s2 - (2.0 * n_bar - 1.0) / (4.0 * n_bar) * h_bar);
        let c = h_bar / 2.0;
        sum_a += a;
        sum_abc += a + b + c;
    }
    sum_a / sum_abc
}

/// The homozygote and heterozygote counts per population and allele that
/// `WeirFst::terms` takes, from the genotypes.
fn counts(pops: &[Vec<&str>], n_alleles: usize) -> (Vec<Vec<u64>>, Vec<Vec<u64>>) {
    let mut hom = vec![vec![0; n_alleles]; pops.len()];
    let mut het = vec![vec![0; n_alleles]; pops.len()];
    for (i, p) in pops.iter().enumerate() {
        for g in p {
            let Some((x, y)) = g.split_once('/') else {
                continue;
            };
            let (Ok(x), Ok(y)) = (x.parse::<usize>(), y.parse::<usize>()) else {
                continue;
            };
            if x == y {
                hom[i][x] += 1;
            } else {
                het[i][x] += 1;
                het[i][y] += 1;
            }
        }
    }
    (hom, het)
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

#[test]
fn fst_oracle() {
    let cases: Vec<(usize, Vec<Vec<&str>>)> = vec![
        // biallelic, 2 populations of different sizes, a missing and a
        // half-missing genotype (neither is complete)
        (
            2,
            vec![
                vec!["0/0", "0/1", "0/1", "1/1", "./."],
                vec!["1/1", "0/1", "1/1", "./1"],
            ],
        ),
        // biallelic, 3 populations
        (
            2,
            vec![
                vec!["0/0", "0/1", "0/1", "1/1", "./."],
                vec!["1/1", "0/1", "1/1"],
                vec!["0/0", "0/0", "0/1", "0/1"],
            ],
        ),
        // triallelic, 2 populations, heterozygotes without REF
        (
            3,
            vec![
                vec!["0/1", "1/2", "2/2", "0/0"],
                vec!["0/2", "0/0", "1/1", "./1", "2/1"],
            ],
        ),
        // triallelic, 3 populations
        (
            3,
            vec![
                vec!["0/1", "1/2", "2/2", "0/0"],
                vec!["0/2", "0/0", "1/1"],
                vec!["2/2", "1/2", "0/0", "0/0", "1/1"],
            ],
        ),
        // the same frequencies everywhere: a negative estimate
        (
            2,
            vec![vec!["0/0", "0/1", "1/1"], vec!["0/0", "0/1", "1/1"]],
        ),
    ];
    let mut generated = Vec::new();
    let mut rng = Lcg(7);
    for _ in 0..300 {
        let n_alleles = 2 + rng.next(3) as usize;
        let n_pops = 2 + rng.next(3) as usize;
        let pops: Vec<Vec<String>> = (0..n_pops)
            .map(|_| {
                (0..2 + rng.next(6))
                    .map(|_| match rng.next(12) {
                        0 => "./.".to_string(),
                        1 => format!("./{}", rng.next(n_alleles as u64)),
                        _ => format!(
                            "{}/{}",
                            rng.next(n_alleles as u64),
                            rng.next(n_alleles as u64)
                        ),
                    })
                    .collect()
            })
            .collect();
        generated.push((n_alleles, pops));
    }
    let generated: Vec<(usize, Vec<Vec<&str>>)> = generated
        .iter()
        .map(|(a, p)| {
            (
                *a,
                p.iter()
                    .map(|g| g.iter().map(String::as_str).collect())
                    .collect(),
            )
        })
        .collect();
    let mut checked = 0;
    for (n_alleles, case) in cases.iter().chain(&generated) {
        let want = weir_cockerham_1984(case, *n_alleles);
        let (hom, het) = counts(case, *n_alleles);
        let (sum_a, sum_abc) = WeirFst::terms(&hom, &het);
        let got = sum_a / sum_abc;
        if want.is_nan() || !want.is_finite() {
            // A population with no complete genotype, or no variation.
            assert!(!got.is_finite() || got == want, "{case:?}: {got} vs {want}");
            continue;
        }
        assert!(
            (got - want).abs() < 1e-12,
            "{case:?}: {got} vs the 1984 formula {want}"
        );
        // The statistic counts the genotypes of a record the same way. (Its
        // printed value is the port's: near 0, the two orders of operations
        // leave different rounding noise, such as -2.03541e-16 for
        // -1.85037e-16.)
        let samples: Vec<&str> = case.iter().flatten().copied().collect();
        let alts = ["G", "G,T", "G,T,C", "G,T,C,A"][n_alleles - 2];
        let sizes: Vec<usize> = case.iter().map(Vec::len).collect();
        let files = run(
            &mut WeirFst::new(pops(&sizes), None),
            &[record("chr1", 5, alts, &samples)],
            &NO_MASK,
        );
        assert_eq!(
            files["weir.fst"],
            format!("CHROM\tPOS\tWEIR_AND_COCKERHAM_FST\nchr1\t5\t{}\n", g6(got)),
            "{case:?}"
        );
        checked += 1;
    }
    assert!(checked > 200, "{checked} finite cases");
}

/// Population 2 has no called genotype at the second site: that row is
/// nan, and neither its window nor the summary counts it.
#[test]
fn fst_population_all_missing() {
    let lines = [
        record("chr1", 10, "G", &["0/0", "0/1", "1/1", "1/1"]),
        record("chr1", 20, "G", &["0/1", "0/0", "./.", "./."]),
        record("chr1", 30, "G", &["0/0", "0/0", "0/1", "1/1"]),
    ];
    let files = run(
        &mut WeirFst::new(pops(&[2, 2]), Some((100, 100))),
        &lines,
        &NO_MASK,
    );
    let rows: Vec<&str> = files["weir.fst"].lines().collect();
    assert_eq!(rows[0], "CHROM\tPOS\tWEIR_AND_COCKERHAM_FST");
    assert_eq!(rows[2], "chr1\t20\tnan");
    let fst = |gts: [&str; 4]| {
        let p = vec![vec![gts[0], gts[1]], vec![gts[2], gts[3]]];
        let (hom, het) = counts(&p, 2);
        WeirFst::terms(&hom, &het)
    };
    let (a1, t1) = fst(["0/0", "0/1", "1/1", "1/1"]);
    let (a3, t3) = fst(["0/0", "0/0", "0/1", "1/1"]);
    assert_eq!(rows[1], format!("chr1\t10\t{}", g6(a1 / t1)));
    assert_eq!(rows[3], format!("chr1\t30\t{}", g6(a3 / t3)));
    assert_eq!(rows.len(), 4);
    let mean = (a1 / t1 + a3 / t3) / 2.0;
    let weighted = (a1 + a3) / (t1 + t3);
    assert_eq!(
        files["windowed.weir.fst"],
        format!(
            "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST\nchr1\t1\t100\t2\t{}\t{}\n",
            g6(weighted),
            g6(mean)
        )
    );
    assert_eq!(
        files["weir.fst.summary"],
        format!(
            "N_SITES\tMEAN_FST\tWEIGHTED_FST\n2\t{}\t{}\n",
            g6(mean),
            g6(weighted)
        )
    );
    // No site with an estimate at all: the summary is 0 nan nan and there is
    // no window.
    let files = run(
        &mut WeirFst::new(pops(&[2, 2]), Some((100, 100))),
        &lines[1..2],
        &NO_MASK,
    );
    assert_eq!(
        files["weir.fst.summary"],
        "N_SITES\tMEAN_FST\tWEIGHTED_FST\n0\tnan\tnan\n"
    );
    assert_eq!(
        files["windowed.weir.fst"],
        "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST\n"
    );
}

/// Fst uses a site when every population member that is not masked is
/// diploid; samples in no population do not count, and a half-called
/// genotype is missing.
#[test]
fn fst_diploid_members_only() {
    // samples: two populations of two, then one sample in neither
    let p = || WeirFst::new(pops(&[2, 2]), None);
    let lines = [
        record("chr1", 1, "G", &["0/0", "0/1", "1/1", "1/1", "1"]), // non-member haploid
        record("chr1", 2, "G", &["0/0", "1", "1/1", "1/1", "0/0"]), // member haploid
        record("chr1", 3, "G", &["0/0", "1:3", "1/1", "1/1", "0/0"]), // masked below
        record("chr1", 4, "G", &["0/0", "0/0", "1/1", "1/1", "0/1/1"]), // non-member triploid
        record("chr1", 5, "G", &["0/0", "./1", "1/1", "0/1", "0/0"]), // half-called
    ];
    let got = run(&mut p(), &lines, &NO_MASK)["weir.fst"].clone();
    let positions: Vec<&str> = got
        .lines()
        .skip(1)
        .map(|l| l.split('\t').nth(1).unwrap())
        .collect();
    assert_eq!(positions, ["1", "4", "5"], "{got}");
    let mask = MaskRules {
        min_dp: Some(5),
        ..NO_MASK
    };
    let got = run(&mut p(), &lines, &mask)["weir.fst"].clone();
    let rows: Vec<&str> = got.lines().skip(1).collect();
    assert_eq!(rows.len(), 4, "{got}");
    // the masked genotype is missing: population 1 has one called sample
    let one = |gts: Vec<Vec<&str>>| {
        let (hom, het) = counts(&gts, 2);
        let (a, t) = WeirFst::terms(&hom, &het);
        g6(a / t)
    };
    assert_eq!(
        rows[1],
        format!("chr1\t3\t{}", one(vec![vec!["0/0"], vec!["1/1", "1/1"]]))
    );
    // ./1 is left out of the counts, as vcftools does (pi counts its allele)
    assert_eq!(
        rows[3],
        format!("chr1\t5\t{}", one(vec![vec!["0/0"], vec!["1/1", "0/1"]]))
    );
}

/// A member's genotype of ploidy above two stops Fst, as vcftools does; a
/// masked one does not, and neither does a sample in no population.
#[test]
fn fst_polyploid_member_stops() {
    let lines = [
        record("chr1", 1, "G", &["0/0", "0/1", "1/1", "1/1", "0/0"]),
        // a haploid member before the triploid one does not hide it
        record("chr1", 2, "G", &["1", "0/0/1:3", "1/1", "1/1", "0/0"]),
    ];
    let e = try_run(&mut WeirFst::new(pops(&[2, 2]), None), &lines, &NO_MASK).unwrap_err();
    assert!(
        e.starts_with("polyploid genotypes are not supported (vcftools stops here too)"),
        "{e}"
    );
    let mask = MaskRules {
        min_dp: Some(5),
        ..NO_MASK
    };
    let files = run(&mut WeirFst::new(pops(&[2, 2]), None), &lines, &mask);
    // the haploid member still drops the second site
    assert_eq!(files["weir.fst"].lines().count(), 2, "{files:?}");
}

/// `--min-dp` masks genotypes, which changes the number of called
/// genotypes of their population and so the estimate.
#[test]
fn fst_mask_changes_n() {
    let line = [record(
        "chr1",
        7,
        "G",
        &["0/0:9", "0/1:2", "0/0:9", "1/1:9", "0/1:4", "1/1:9"],
    )];
    let mask = MaskRules {
        min_dp: Some(5),
        ..NO_MASK
    };
    let plain = run(&mut WeirFst::new(pops(&[3, 3]), None), &line, &NO_MASK)["weir.fst"].clone();
    let masked = run(&mut WeirFst::new(pops(&[3, 3]), None), &line, &mask)["weir.fst"].clone();
    let est = |gts: Vec<Vec<&str>>| {
        let (hom, het) = counts(&gts, 2);
        let (a, t) = WeirFst::terms(&hom, &het);
        g6(a / t)
    };
    assert!(plain.ends_with(&format!(
        "\t{}\n",
        est(vec![vec!["0/0", "0/1", "0/0"], vec!["1/1", "0/1", "1/1"]])
    )));
    assert!(masked.ends_with(&format!(
        "\t{}\n",
        est(vec![vec!["0/0", "0/0"], vec!["1/1", "1/1"]])
    )));
    assert_ne!(plain, masked);
}

/// The windows vcftools writes at the end of the file, from every site's
/// terms: `[k·step + 1, k·step + size]` for k from ⌈(pos − size)/step⌉ (at
/// least 0) to ⌈pos/step⌉ − 1, summing a, a + b + c, fst and the sites
/// whose fst is not nan; a window is written when its Σ(a + b + c) is not
/// 0 and it has a site.
fn windows_at_end(sites: &[(&str, u64, f64, f64)], size: u64, step: u64) -> String {
    let mut chroms: Vec<&str> = Vec::new();
    let mut bins: BTreeMap<(usize, u64), (f64, f64, f64, u64)> = BTreeMap::new();
    for &(chrom, pos, a, t) in sites {
        let fst = a / t;
        if fst.is_nan() {
            continue;
        }
        if chroms.last() != Some(&chrom) {
            chroms.push(chrom);
        }
        let c = chroms.len() - 1;
        let first = (pos as f64 - size as f64) / step as f64;
        let first = first.ceil().max(0.0) as u64;
        let last = (pos as f64 / step as f64).ceil() as u64;
        for k in first..last {
            let w = bins.entry((c, k)).or_default();
            w.0 += a;
            w.1 += t;
            w.2 += fst;
            w.3 += 1;
        }
    }
    let mut out = "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST\n".to_string();
    for (&(c, k), &(a, t, f, n)) in &bins {
        if t != 0.0 && !a.is_nan() && !t.is_nan() && n > 0 {
            out += &format!(
                "{}\t{}\t{}\t{}\t{}\t{}\n",
                chroms[c],
                k * step + 1,
                k * step + size,
                n,
                g6(a / t),
                g6(f / n as f64)
            );
        }
    }
    out
}

#[test]
fn fst_windows_progressive() {
    // Sites on three chromosomes, some at one position, some with no
    // estimate (monomorphic, or a population with no call), gaps longer
    // than a window.
    let mut rng = Lcg(11);
    let mut lines = Vec::new();
    let mut sites = Vec::new();
    for (chrom, n) in [("chr1", 60u64), ("chr2", 5), ("chr3", 40)] {
        let mut pos = 1 + rng.next(30);
        for _ in 0..n {
            let gts: Vec<String> = (0..6)
                .map(|_| match rng.next(9) {
                    0 => "./.".to_string(),
                    k => format!("{}/{}", u64::from(k > 4), u64::from(k > 6)),
                })
                .collect();
            let refs: Vec<&str> = gts.iter().map(String::as_str).collect();
            let p = vec![refs[..3].to_vec(), refs[3..].to_vec()];
            let (hom, het) = counts(&p, 2);
            let (a, t) = WeirFst::terms(&hom, &het);
            lines.push(record(chrom, pos, "G", &refs));
            sites.push((chrom, pos, a, t));
            pos += match rng.next(10) {
                0 => 0,
                1 => 250,
                _ => 1 + rng.next(40),
            };
        }
    }
    for (size, step) in [
        (100, 100),
        (100, 50),
        (100, 30),
        (37, 1),
        (1, 1),
        (500, 125),
    ] {
        let mut stat = WeirFst::new(pops(&[3, 3]), Some((size, step)));
        let files = run(&mut stat, &lines, &NO_MASK);
        assert_eq!(
            files["windowed.weir.fst"],
            windows_at_end(&sites, size, step),
            "size {size} step {step}"
        );
        assert!(
            stat.peak_open() as u64 <= size.div_ceil(step) + 1,
            "size {size} step {step}: {} windows open",
            stat.peak_open()
        );
    }
    // Sites every 10 bases keep open all the windows a position falls in,
    // ⌈100/30⌉ = 4 (within the bound of 5).
    let dense: Vec<String> = (1..=40)
        .map(|i| record("chr1", i * 10, "G", &["0/0", "0/1", "1/1", "1/1"]))
        .collect();
    let mut stat = WeirFst::new(pops(&[2, 2]), Some((100, 30)));
    run(&mut stat, &dense, &NO_MASK);
    assert_eq!(stat.peak_open(), 4);

    // Windows stream their rows, so they need sorted input; the per-site
    // file alone takes any order.
    let unsorted = [
        record("chr1", 50, "G", &["0/0", "0/1", "1/1", "1/1"]),
        record("chr2", 10, "G", &["0/0", "0/1", "1/1", "1/1"]),
        record("chr1", 60, "G", &["0/0", "0/1", "1/1", "1/1"]),
    ];
    let e = try_run(
        &mut WeirFst::new(pops(&[2, 2]), Some((100, 100))),
        &unsorted,
        &NO_MASK,
    )
    .unwrap_err();
    assert!(
        e.contains("--fst-window-size needs input sorted") && e.contains("chr1:60"),
        "{e}"
    );
    let backwards = [
        record("chr1", 50, "G", &["0/0", "0/1", "1/1", "1/1"]),
        record("chr1", 40, "G", &["0/0", "0/1", "1/1", "1/1"]),
    ];
    let e = try_run(
        &mut WeirFst::new(pops(&[2, 2]), Some((100, 100))),
        &backwards,
        &NO_MASK,
    )
    .unwrap_err();
    assert!(e.contains("a position goes backwards at chr1:40"), "{e}");
    let files = run(&mut WeirFst::new(pops(&[2, 2]), None), &unsorted, &NO_MASK);
    assert_eq!(files["weir.fst"].lines().count(), 4);
    assert!(!files.contains_key("windowed.weir.fst"));
}

fn write(dir: &std::path::Path, name: &str, text: &str) -> PathBuf {
    let p = dir.join(name);
    std::fs::write(&p, text).unwrap();
    p
}

fn selected(names: &[&str]) -> Vec<(usize, String)> {
    names
        .iter()
        .enumerate()
        .map(|(i, n)| (i * 2, n.to_string()))
        .collect()
}

/// One name per line: its first word, blank lines skipped, a CR before the
/// newline dropped; names that are not selected samples are left out.
#[test]
fn populations_load() {
    let d = tempfile::tempdir().unwrap();
    let p1 = write(d.path(), "p1.txt", "s1\ns2 more words\n\n  s3\tx\n");
    let p2 = write(d.path(), "p2.txt", "s4\r\nnobody\ns5\ns5\nunselected");
    let sel = selected(&["s1", "s2", "s3", "s4", "s5"]);
    let pops = Populations::load(&[p1.clone(), p2.clone()], &sel).unwrap();
    assert_eq!(
        pops.names,
        [p1.display().to_string(), p2.display().to_string()]
    );
    // indices in the record, as `selected` gives them
    assert_eq!(pops.members, [vec![0, 2, 4], vec![6, 8]]);

    let usage = |files: &[PathBuf], sel: &[(usize, String)]| match Populations::load(files, sel) {
        Err(gloci_io::IoError::Usage(e)) => e.0,
        other => panic!("expected a usage error, got {other:?}"),
    };
    let e = usage(std::slice::from_ref(&p1), &sel);
    assert!(e.contains("at least two"), "{e}");
    // a population with no selected sample
    let p3 = write(d.path(), "p3.txt", "nobody\nelse\n");
    let e = usage(&[p1.clone(), p3.clone()], &sel);
    assert!(
        e.contains(&p3.display().to_string()) && e.contains("no selected sample"),
        "{e}"
    );
    // after selection: s4 and s5 are not selected
    let e = usage(&[p1.clone(), p2.clone()], &selected(&["s1", "s2", "s3"]));
    assert!(e.contains(&p2.display().to_string()), "{e}");
    // a sample in two files
    let p4 = write(d.path(), "p4.txt", "nobody\ns3\n");
    let e = usage(&[p1.clone(), p2.clone(), p4.clone()], &sel);
    assert!(
        e.contains("s3")
            && e.contains(&p1.display().to_string())
            && e.contains(&p4.display().to_string()),
        "{e}"
    );
    // a file that does not exist
    match Populations::load(&[p1, d.path().join("absent.txt")], &sel) {
        Err(gloci_io::IoError::Open { path, .. }) => assert!(path.ends_with("absent.txt")),
        other => panic!("expected an open error, got {other:?}"),
    }
}
