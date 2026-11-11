//! Per-sample counters that selection workers keep per batch.
//!
//! `--depth`, `--missing-indv` and `--het` keep integer counters per sample.
//! Integer sums come out the same in any grouping, so a worker adds the
//! records of its batch and the driver adds the batches. `--het`'s expected
//! homozygosity is a floating-point sum per sample, which must be added in
//! record order: the worker only records which samples were fully called at
//! each contributing site (a bitset), and the driver does the additions.

use crate::gt::{GenotypeSource, SampleFields};
use crate::stats::SiteStats;

/// Which per-sample counters to keep.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct TallyWants {
    pub depth: bool,
    pub missing: bool,
    pub het: bool,
}

impl TallyWants {
    pub fn is_empty(&self) -> bool {
        !(self.depth || self.missing || self.het)
    }

    pub fn or(self, o: TallyWants) -> TallyWants {
        TallyWants {
            depth: self.depth || o.depth,
            missing: self.missing || o.missing,
            het: self.het || o.het,
        }
    }
}

/// Counters per selected sample (by position in the selection).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SampleTally {
    pub dp_sum: Vec<u64>,
    pub dp_n: Vec<u64>,
    pub data: Vec<u64>,
    pub filtered: Vec<u64>,
    pub miss: Vec<u64>,
    pub o_hom: Vec<u64>,
    pub n_sites: Vec<u64>,
}

impl SampleTally {
    fn columns(&mut self) -> [&mut Vec<u64>; 7] {
        [
            &mut self.dp_sum,
            &mut self.dp_n,
            &mut self.data,
            &mut self.filtered,
            &mut self.miss,
            &mut self.o_hom,
            &mut self.n_sites,
        ]
    }

    /// Zero every counter, keeping the buffers.
    pub fn clear(&mut self) {
        for c in self.columns() {
            c.iter_mut().for_each(|v| *v = 0);
        }
    }

    /// Size the columns of the requested counters for `k` samples; the
    /// others stay empty.
    fn ensure(&mut self, k: usize, wants: TallyWants) {
        let used = [
            (wants.depth, 0),
            (wants.depth, 1),
            (wants.missing, 2),
            (wants.missing, 3),
            (wants.missing, 4),
            (wants.het, 5),
            (wants.het, 6),
        ];
        let cols = self.columns();
        for (want, i) in used {
            let c = &mut *cols[i];
            if want && c.len() != k {
                c.clear();
                c.resize(k, 0);
            }
        }
    }

    /// Add `o`'s counters to these.
    pub fn add(&mut self, o: &SampleTally) {
        let mut src = [
            &o.dp_sum,
            &o.dp_n,
            &o.data,
            &o.filtered,
            &o.miss,
            &o.o_hom,
            &o.n_sites,
        ]
        .into_iter();
        for c in self.columns() {
            let s = src.next().expect("seven columns");
            if c.len() < s.len() {
                c.resize(s.len(), 0);
            }
            for (a, b) in c.iter_mut().zip(s) {
                *a += *b;
            }
        }
    }

    /// Count one record. `samples` are the selected samples' indices in the
    /// record, `src` reads them with the masked ones as missing (so `f.masked`
    /// is set for them), `st` is the record's site statistics. When the
    /// site contributes to `--het`, the words of a bitset (bit k set when
    /// selected sample k has a fully called diploid genotype) are appended to
    /// `called`, and `true` is returned.
    ///
    /// A polyploid site adds nothing to `--missing-indv` and `--het` (the
    /// statistics refuse the record).
    pub fn add_record(
        &mut self,
        wants: TallyWants,
        samples: &[usize],
        src: &dyn GenotypeSource,
        f: &mut SampleFields,
        st: &SiteStats,
        called: &mut Vec<u64>,
    ) -> Result<bool, String> {
        let k = samples.len();
        self.ensure(k, wants);
        let diploid_ok = !st.polyploid;
        let het_site = wants.het
            && diploid_ok
            && st.allele_counts.len() == 2
            && !st.haploid
            && st.allele_counts[0] != 0
            && st.allele_counts[1] != 0;
        let bits_at = called.len();
        if het_site {
            called.resize(bits_at + k.div_ceil(64), 0);
        }
        for (j, &i) in samples.iter().enumerate() {
            src.sample(i, f)?;
            if wants.depth
                && let Some(d) = f.dp
            {
                self.dp_sum[j] += u64::from(d);
                self.dp_n[j] += 1;
            }
            if wants.missing && diploid_ok {
                if f.masked {
                    self.filtered[j] += 1;
                } else {
                    self.data[j] += 1;
                    if f.gt
                        .as_ref()
                        .is_none_or(|g| g.alleles.first().is_none_or(Option::is_none))
                    {
                        self.miss[j] += 1;
                    }
                }
            }
            if het_site
                && let Some(g) = &f.gt
                && let [Some(a), Some(b)] = g.alleles[..]
            {
                self.n_sites[j] += 1;
                if a == b {
                    self.o_hom[j] += 1;
                }
                called[bits_at + j / 64] |= 1 << (j % 64);
            }
        }
        Ok(het_site)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn only_the_requested_columns_are_sized() {
        let mut t = SampleTally::default();
        t.ensure(
            5,
            TallyWants {
                het: true,
                ..TallyWants::default()
            },
        );
        assert_eq!((t.o_hom.len(), t.n_sites.len()), (5, 5));
        assert!(t.dp_sum.is_empty() && t.dp_n.is_empty());
        assert!(t.data.is_empty() && t.filtered.is_empty() && t.miss.is_empty());
    }
}
