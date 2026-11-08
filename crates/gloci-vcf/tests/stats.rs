use gloci_vcf::{Genotype, GenotypeSource, MaskRules, SampleFields, SiteStats, site_stats};

struct Fake(Vec<SampleFields>);
impl GenotypeSource for Fake {
    fn n_samples(&self) -> usize {
        self.0.len()
    }
    fn sample(&self, i: usize, out: &mut SampleFields) -> Result<(), String> {
        let s = &self.0[i];
        out.gt = s.gt.as_ref().map(|g| Genotype {
            alleles: g.alleles.clone(),
            phased: g.phased,
        });
        out.dp = s.dp;
        out.gq = s.gq;
        out.ft_pass = s.ft_pass;
        Ok(())
    }
}

fn sf(a: &[Option<u32>], phased: bool, dp: Option<u32>, gq: Option<u32>, ft: bool) -> SampleFields {
    SampleFields {
        gt: Some(Genotype {
            alleles: a.to_vec(),
            phased,
        }),
        dp,
        gq,
        ft_pass: ft,
    }
}
fn g(a: &[Option<u32>]) -> SampleFields {
    sf(a, false, None, None, true)
}
const NONE: MaskRules = MaskRules {
    min_gq: None,
    min_dp: None,
    max_dp: None,
    mask_filtered: false,
};

fn stats_of(counts: Vec<u64>) -> SiteStats {
    SiteStats {
        called_alleles: counts.iter().sum(),
        allele_counts: counts,
        samples: 1,
        missing: 0,
        dp_sum: 0,
        dp_n: 0,
        hwe_counts: None,
        all_phased: true,
    }
}

#[test]
fn stats() {
    let s = stats_of(vec![7, 3]);
    assert!((s.maf().unwrap() - 0.3).abs() < 1e-12);
    assert_eq!(s.mac(), Some(3));
    let s = stats_of(vec![5, 3, 2]);
    assert!((s.maf().unwrap() - 0.3).abs() < 1e-12);
    assert_eq!(s.mac(), Some(3));
    assert_eq!(stats_of(vec![0, 0]).maf(), None);
    assert_eq!(stats_of(vec![10]).mac(), Some(0));
}

#[test]
fn stats_from_genotypes() {
    let c = Some(0);
    let a = Some(1);
    let src = Fake(vec![
        sf(&[c, a], true, Some(10), None, true),
        sf(&[a, a], true, Some(20), None, true),
        g(&[None, a]), // missing, but its allele 1 still counts
        g(&[None, None]),
        SampleFields {
            gt: None,
            dp: None,
            gq: None,
            ft_pass: true,
        },
        g(&[c]),       // haploid: counts alleles, not HWE, not unphased
        g(&[c, c, a]), // triploid
    ]);
    let mut masked = Vec::new();
    let s = site_stats(&src, None, &NONE, 2, &mut masked).unwrap();
    assert_eq!(s.samples, 7);
    assert_eq!(s.missing, 3);
    assert_eq!(s.allele_counts, vec![1 + 1 + 2, 1 + 2 + 1 + 1]);
    assert_eq!(s.called_alleles, 9);
    assert_eq!(s.hwe_counts, Some((0, 1, 1)));
    assert_eq!((s.dp_sum, s.dp_n), (30, 2));
    assert_eq!(s.mean_dp(), Some(15.0));
    assert!((s.missing_fraction() - 3.0 / 7.0).abs() < 1e-12);
    assert!(!s.all_phased); // the triploid and ./1 are unphased, ploidy >= 2
}

#[test]
fn stats_phased_and_hwe_none() {
    let src = Fake(vec![
        sf(&[Some(0), Some(1)], true, None, None, true),
        g(&[Some(0)]),
    ]);
    let mut m = Vec::new();
    let s = site_stats(&src, None, &NONE, 2, &mut m).unwrap();
    assert!(s.all_phased);
    let s = site_stats(&src, None, &NONE, 3, &mut m).unwrap();
    assert_eq!(s.hwe_counts, None);
    assert_eq!(s.allele_counts.len(), 3);
    assert!(site_stats(&Fake(vec![g(&[Some(5)])]), None, &NONE, 2, &mut m).is_err());
}

#[test]
fn masking() {
    let a = Some(1);
    let src = Fake(vec![
        sf(&[a, a], false, Some(5), Some(50), true),   // low DP
        sf(&[a, a], false, None, None, true),          // missing DP/GQ: not masked
        sf(&[a, a], false, Some(50), Some(5), true),   // low GQ
        sf(&[a, a], false, Some(500), Some(50), true), // high DP
        sf(&[a, a], false, Some(20), Some(50), false), // FT fail
        sf(&[a, a], false, Some(20), Some(50), true),
    ]);
    let mut masked = Vec::new();
    let r = MaskRules {
        min_gq: None,
        min_dp: Some(10),
        max_dp: None,
        mask_filtered: false,
    };
    let s = site_stats(&src, None, &r, 2, &mut masked).unwrap();
    assert_eq!(masked, vec![true, false, false, false, false, false]);
    assert_eq!(s.missing, 1);
    assert_eq!(s.called_alleles, 10);
    assert_eq!((s.dp_sum, s.dp_n), (500 + 20 + 50 + 20, 4));

    let r = MaskRules {
        min_gq: Some(10),
        min_dp: None,
        max_dp: Some(100),
        mask_filtered: true,
    };
    site_stats(&src, None, &r, 2, &mut masked).unwrap();
    assert_eq!(masked, vec![false, false, true, true, true, false]);

    // Sample selection: unselected samples are not counted at all.
    let keep = [true, false, true, false, false, true];
    let s = site_stats(&src, Some(&keep), &NONE, 2, &mut masked).unwrap();
    assert_eq!(s.samples, 3);
    assert_eq!(s.called_alleles, 6);
    assert_eq!(masked.len(), 6);
}
