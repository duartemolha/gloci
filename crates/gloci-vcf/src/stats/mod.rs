//! Per-site statistics over the selected, unmasked samples.

mod engine;
mod report;
mod sample;
mod site;
mod tally;
mod whole;

pub use engine::{Out, Stat, StatSink, StatWrite, StatsEngine};
pub use report::{SUMMARY_SCHEMA_VERSION, Summary};
pub use sample::{Heterozygosity, MissingSample, SampleDepth};
pub use site::{
    Freq, GenoDepth, Hardy, MissingSite, Singletons, SiteDepth, SiteMeanDepth, SiteQuality,
};
pub use tally::{SampleTally, TallyWants};
pub use whole::{FilterSummary, IndelHist, SnpDensity, TsTv};

use crate::gt::{GenotypeSource, SampleFields};

/// Genotype masking thresholds (D50). A masked sample counts as missing.
/// Missing DP or GQ never masks.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MaskRules {
    pub min_gq: Option<u32>,
    pub min_dp: Option<u32>,
    pub max_dp: Option<u32>,
    /// Mask samples whose FORMAT/FT is not `PASS`, `.` or absent.
    pub mask_filtered: bool,
}

impl MaskRules {
    /// True when no genotype masking rule is set.
    pub fn is_empty(&self) -> bool {
        self.min_gq.is_none()
            && self.min_dp.is_none()
            && self.max_dp.is_none()
            && !self.mask_filtered
    }

    fn masks(&self, f: &SampleFields) -> bool {
        self.min_gq.is_some_and(|m| f.gq.is_some_and(|g| g < m))
            || self.min_dp.is_some_and(|m| f.dp.is_some_and(|d| d < m))
            || self.max_dp.is_some_and(|m| f.dp.is_some_and(|d| d > m))
            || (self.mask_filtered && !f.ft_pass)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SiteStats {
    /// Called alleles per allele index (a `./1` genotype still counts its 1).
    pub allele_counts: Vec<u64>,
    pub called_alleles: u64,
    /// Selected samples.
    pub samples: u64,
    /// Selected samples that are masked or whose genotype is missing (D50).
    pub missing: u64,
    pub dp_sum: u64,
    pub dp_n: u64,
    /// Sum of squared depths over the same samples as `dp_sum`.
    pub dp_sumsq: u64,
    /// Selected samples the mask rules masked (`vcftools` N_GENOTYPE_FILTERED).
    pub filtered: u64,
    /// Alleles of the selected, unmasked samples: their ploidy summed, where
    /// `.` is one allele and a sample without GT counts as two.
    pub data_alleles: u64,
    /// Missing alleles among `data_alleles`.
    pub missing_alleles: u64,
    /// An unmasked selected sample has a one-allele genotype (`.` too).
    pub haploid: bool,
    /// An unmasked selected sample has more than two alleles.
    pub polyploid: bool,
    /// Fully called diploid genotypes as (0/0, 0/1, 1/1); `None` unless the
    /// site is biallelic, in which case `--hwe` passes it.
    pub hwe_counts: Option<(u64, u64, u64)>,
    /// False if any genotype of ploidy 2 or more, called or not, is written
    /// unphased (a `/` separator).
    pub all_phased: bool,
}

/// An empty site: no samples, so `all_phased` is true, as `site_stats`
/// leaves it for a site without unphased genotypes.
impl Default for SiteStats {
    fn default() -> SiteStats {
        SiteStats {
            allele_counts: Vec::new(),
            called_alleles: 0,
            samples: 0,
            missing: 0,
            dp_sum: 0,
            dp_n: 0,
            dp_sumsq: 0,
            filtered: 0,
            data_alleles: 0,
            missing_alleles: 0,
            haploid: false,
            polyploid: false,
            hwe_counts: None,
            all_phased: true,
        }
    }
}

impl SiteStats {
    /// Second-highest allele count over called alleles.
    pub fn maf(&self) -> Option<f64> {
        if self.called_alleles == 0 {
            return None;
        }
        Some(self.mac()? as f64 / self.called_alleles as f64)
    }

    /// Second-highest allele count; `None` with no called alleles.
    pub fn mac(&self) -> Option<u64> {
        if self.called_alleles == 0 {
            return None;
        }
        let (mut top, mut second) = (0u64, 0u64);
        for &c in &self.allele_counts {
            if c > top {
                second = top;
                top = c;
            } else if c > second {
                second = c;
            }
        }
        Some(second)
    }

    /// Exact HWE p-value of the called diploid genotypes; `None` for a site
    /// that is not biallelic (such a site passes `--hwe`). This is the only
    /// place the `hwe_counts` order is mapped onto `hwe_p`'s arguments.
    pub fn hwe_p(&self) -> Option<f64> {
        let (hom_ref, het, hom_alt) = self.hwe_counts?;
        Some(crate::hwe::hwe_p(het, hom_ref, hom_alt))
    }

    /// Called count per ALT allele (index 0 of the result is ALT 1); the
    /// `AC` of filter expressions. Empty for a site with no ALT.
    pub fn ac(&self) -> &[u64] {
        self.allele_counts.get(1..).unwrap_or(&[])
    }

    /// Called alleles; the `AN` of filter expressions.
    pub fn an(&self) -> u64 {
        self.called_alleles
    }

    /// bcftools' `MAF`: `min(f, 1 - f)` for the frequency `f` of the first
    /// ALT among called alleles. `None` with no called allele or no ALT.
    /// This is not [`SiteStats::maf`] (the `--min-maf` definition).
    pub fn expr_maf(&self) -> Option<f64> {
        if self.called_alleles == 0 || self.allele_counts.len() < 2 {
            return None;
        }
        let f = self.allele_counts[1] as f64 / self.called_alleles as f64;
        Some(f.min(1.0 - f))
    }

    pub fn missing_fraction(&self) -> f64 {
        if self.samples == 0 {
            0.0
        } else {
            self.missing as f64 / self.samples as f64
        }
    }

    pub fn mean_dp(&self) -> Option<f64> {
        (self.dp_n > 0).then(|| self.dp_sum as f64 / self.dp_n as f64)
    }

    /// Sample variance (n - 1) of the depths; `None` with fewer than two.
    /// The numerator is exact in integers, so only the final division rounds.
    pub fn var_dp(&self) -> Option<f64> {
        if self.dp_n < 2 {
            return None;
        }
        let n = u128::from(self.dp_n);
        let num = n * u128::from(self.dp_sumsq) - u128::from(self.dp_sum) * u128::from(self.dp_sum);
        Some(num as f64 / (n * (n - 1)) as f64)
    }
}

/// Compute the statistics of one site. `keep` selects samples (`None` = all);
/// `masked` is resized to the sample count and says which samples the mask
/// rules masked (selected or not, only selected ones are ever set).
pub fn site_stats(
    src: &dyn GenotypeSource,
    keep: Option<&[bool]>,
    mask: &MaskRules,
    n_alleles: usize,
    masked: &mut Vec<bool>,
) -> Result<SiteStats, String> {
    let mut s = SiteStats::default();
    let mut f = SampleFields::default();
    site_stats_into(src, keep, mask, n_alleles, masked, &mut f, &mut s)?;
    Ok(s)
}

/// [`site_stats`] into a reused `SiteStats` and sample buffer `f`:
/// allocation-free once `out`, `masked` and `f` have grown to the largest
/// site.
pub fn site_stats_into(
    src: &dyn GenotypeSource,
    keep: Option<&[bool]>,
    mask: &MaskRules,
    n_alleles: usize,
    masked: &mut Vec<bool>,
    f: &mut SampleFields,
    out: &mut SiteStats,
) -> Result<(), String> {
    let n = src.n_samples();
    if let Some(k) = keep
        && k.len() != n
    {
        return Err(format!(
            "the sample selection has {} entries; the record has {n} samples",
            k.len()
        ));
    }
    masked.clear();
    masked.resize(n, false);
    let mut counts = std::mem::take(&mut out.allele_counts);
    counts.clear();
    counts.resize(n_alleles, 0);
    *out = SiteStats {
        allele_counts: counts,
        hwe_counts: (n_alleles == 2).then_some((0, 0, 0)),
        ..SiteStats::default()
    };
    let s = out;
    for i in 0..n {
        if keep.is_some_and(|k| !k[i]) {
            continue;
        }
        s.samples += 1;
        src.sample(i, f)?;
        // Phasing is read from the separators, called or not (vcftools
        // --phased): `.|.` is phased, `./.` is not. Masking keeps them.
        if let Some(g) = &f.gt
            && g.alleles.len() >= 2
            && !g.phased
        {
            s.all_phased = false;
        }
        if mask.masks(f) {
            masked[i] = true;
            s.missing += 1;
            s.filtered += 1;
            continue;
        }
        if let Some(dp) = f.dp {
            s.dp_sum += u64::from(dp);
            s.dp_sumsq += u64::from(dp) * u64::from(dp);
            s.dp_n += 1;
        }
        let Some(g) = &f.gt else {
            s.missing += 1;
            s.data_alleles += 2;
            s.missing_alleles += 2;
            continue;
        };
        s.data_alleles += g.alleles.len() as u64;
        s.haploid |= g.alleles.len() == 1;
        s.polyploid |= g.alleles.len() > 2;
        let mut fully_called = true;
        for a in &g.alleles {
            match a {
                Some(a) => {
                    let slot = s.allele_counts.get_mut(*a as usize).ok_or_else(|| {
                        format!("genotype allele {a} but the site lists {n_alleles} alleles")
                    })?;
                    *slot += 1;
                    s.called_alleles += 1;
                }
                None => {
                    fully_called = false;
                    s.missing_alleles += 1;
                }
            }
        }
        if !fully_called {
            s.missing += 1;
            continue;
        }
        if let (Some(c), [Some(a), Some(b)]) = (&mut s.hwe_counts, &g.alleles[..]) {
            match a + b {
                0 => c.0 += 1,
                1 => c.1 += 1,
                _ => c.2 += 1,
            }
        }
    }
    Ok(())
}
