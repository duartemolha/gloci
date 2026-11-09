//! Per-site statistics over the selected, unmasked samples.

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
    /// Fully called diploid genotypes as (0/0, 0/1, 1/1); `None` unless the
    /// site is biallelic, in which case `--hwe` passes it.
    pub hwe_counts: Option<(u64, u64, u64)>,
    /// False if any genotype of ploidy 2 or more, called or not, is written
    /// unphased (a `/` separator).
    pub all_phased: bool,
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
    let mut s = SiteStats {
        allele_counts: vec![0; n_alleles],
        called_alleles: 0,
        samples: 0,
        missing: 0,
        dp_sum: 0,
        dp_n: 0,
        hwe_counts: (n_alleles == 2).then_some((0, 0, 0)),
        all_phased: true,
    };
    let mut f = SampleFields::default();
    for i in 0..n {
        if keep.is_some_and(|k| !k[i]) {
            continue;
        }
        s.samples += 1;
        src.sample(i, &mut f)?;
        // Phasing is read from the separators, called or not (vcftools
        // --phased): `.|.` is phased, `./.` is not. Masking keeps them.
        if let Some(g) = &f.gt
            && g.alleles.len() >= 2
            && !g.phased
        {
            s.all_phased = false;
        }
        if mask.masks(&f) {
            masked[i] = true;
            s.missing += 1;
            continue;
        }
        if let Some(dp) = f.dp {
            s.dp_sum += u64::from(dp);
            s.dp_n += 1;
        }
        let Some(g) = &f.gt else {
            s.missing += 1;
            continue;
        };
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
                None => fully_called = false,
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
    Ok(s)
}
