//! The per-sample statistics: one row per selected sample, written when the
//! input ends. Each keeps a few counters per sample, so memory is bounded by
//! the number of samples, never by the number of records.

use super::engine::{Out, Stat, StatSink};
use crate::filter::SiteView;
use crate::gt::{GenotypeSource, SampleFields};
use crate::stats::SiteStats;
use gloci_core::cfmt::{write_fixed, write_g6, write_u64};
use std::io;

fn stats_of(st: Option<&SiteStats>) -> Result<&SiteStats, String> {
    st.ok_or_else(|| "internal error: this statistic needs the site statistics".to_string())
}

fn genos(gt: Option<&dyn GenotypeSource>) -> Result<&dyn GenotypeSource, String> {
    gt.ok_or_else(|| "internal error: this statistic needs genotypes".to_string())
}

fn header_row(out: &mut Out, line: &str) {
    out.buf.extend_from_slice(line.as_bytes());
    out.buf.push(b'\n');
}

/// `--depth` (`.idepth`): per sample, the number of sites that have a DP
/// and the mean of those depths. A masked genotype has no depth.
pub struct SampleDepth {
    samples: Vec<(usize, String)>,
    sum: Vec<u64>,
    n: Vec<u64>,
    f: SampleFields,
    out: Out,
}

impl SampleDepth {
    /// `samples`: (index in the record, name) of every selected sample.
    pub fn new(samples: Vec<(usize, String)>) -> SampleDepth {
        let k = samples.len();
        SampleDepth {
            samples,
            sum: vec![0; k],
            n: vec![0; k],
            f: SampleFields::default(),
            out: Out::closed(),
        }
    }
}

impl Stat for SampleDepth {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, self.ext())
    }
    fn site(
        &mut self,
        _s: &SiteView,
        gt: Option<&dyn GenotypeSource>,
        _st: Option<&SiteStats>,
    ) -> Result<(), String> {
        let gt = genos(gt)?;
        for (k, (i, _)) in self.samples.iter().enumerate() {
            gt.sample(*i, &mut self.f)?;
            if let Some(d) = self.f.dp {
                self.sum[k] += u64::from(d);
                self.n[k] += 1;
            }
        }
        Ok(())
    }
    fn finish(&mut self) -> io::Result<()> {
        header_row(&mut self.out, "INDV\tN_SITES\tMEAN_DEPTH");
        for (k, (_, name)) in self.samples.iter().enumerate() {
            let b = &mut self.out.buf;
            b.extend_from_slice(name.as_bytes());
            b.push(b'\t');
            write_u64(b, self.n[k]);
            b.push(b'\t');
            // 0 sites is nan, as vcftools prints.
            write_g6(b, self.sum[k] as f64 / self.n[k] as f64);
            b.push(b'\n');
        }
        self.out.row().map_err(io::Error::other)?;
        self.out.commit()
    }
    fn needs_genotypes(&self) -> bool {
        true
    }
    fn ext(&self) -> &'static str {
        "idepth"
    }
}

/// `--missing-indv` (`.imiss`): per sample, the genotypes that were not
/// masked (`N_DATA`), those masked (`N_GENOTYPES_FILTERED`), and how many
/// of the unmasked ones are missing: the first allele is `.` (`./1` and a
/// lone `.` are, `1/.` is not, as in vcftools), or the sample has no GT.
pub struct MissingSample {
    samples: Vec<(usize, String)>,
    data: Vec<u64>,
    filtered: Vec<u64>,
    miss: Vec<u64>,
    f: SampleFields,
    out: Out,
}

impl MissingSample {
    /// `samples`: (index in the record, name) of every selected sample.
    pub fn new(samples: Vec<(usize, String)>) -> MissingSample {
        let k = samples.len();
        MissingSample {
            samples,
            data: vec![0; k],
            filtered: vec![0; k],
            miss: vec![0; k],
            f: SampleFields::default(),
            out: Out::closed(),
        }
    }
}

impl Stat for MissingSample {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, self.ext())
    }
    fn site(
        &mut self,
        _s: &SiteView,
        gt: Option<&dyn GenotypeSource>,
        st: Option<&SiteStats>,
    ) -> Result<(), String> {
        super::site::diploid_only(stats_of(st)?)?;
        let gt = genos(gt)?;
        for (k, (i, _)) in self.samples.iter().enumerate() {
            gt.sample(*i, &mut self.f)?;
            if self.f.masked {
                self.filtered[k] += 1;
                continue;
            }
            self.data[k] += 1;
            if self
                .f
                .gt
                .as_ref()
                .is_none_or(|g| g.alleles.first().is_none_or(Option::is_none))
            {
                self.miss[k] += 1;
            }
        }
        Ok(())
    }
    fn finish(&mut self) -> io::Result<()> {
        header_row(
            &mut self.out,
            "INDV\tN_DATA\tN_GENOTYPES_FILTERED\tN_MISS\tF_MISS",
        );
        for (k, (_, name)) in self.samples.iter().enumerate() {
            let b = &mut self.out.buf;
            b.extend_from_slice(name.as_bytes());
            for v in [self.data[k], self.filtered[k], self.miss[k]] {
                b.push(b'\t');
                write_u64(b, v);
            }
            b.push(b'\t');
            // No unmasked genotype at all is nan, as vcftools prints.
            write_g6(b, self.miss[k] as f64 / self.data[k] as f64);
            b.push(b'\n');
        }
        self.out.row().map_err(io::Error::other)?;
        self.out.commit()
    }
    fn needs_genotypes(&self) -> bool {
        true
    }
    fn ext(&self) -> &'static str {
        "imiss"
    }
}

/// `--het` (`.het`): per sample, observed and expected homozygous sites and
/// the inbreeding coefficient F = (O - E) / (N - E).
///
/// A site contributes only when it has exactly two alleles, both of which
/// are seen among the called alleles, and no unmasked sample is haploid. On
/// such a site, each sample with a fully called genotype adds one to
/// `N_SITES`, one to `O(HOM)` when homozygous, and the expected
/// homozygosity `1 - 2pq * n / (n - 1)` to `E(HOM)` (p the REF frequency, n
/// the called alleles). A sample with no contributing site has no row.
pub struct Heterozygosity {
    samples: Vec<(usize, String)>,
    o_hom: Vec<u64>,
    e_hom: Vec<f64>,
    n_sites: Vec<u64>,
    f: SampleFields,
    out: Out,
}

impl Heterozygosity {
    /// `samples`: (index in the record, name) of every selected sample.
    pub fn new(samples: Vec<(usize, String)>) -> Heterozygosity {
        let k = samples.len();
        Heterozygosity {
            samples,
            o_hom: vec![0; k],
            e_hom: vec![0.0; k],
            n_sites: vec![0; k],
            f: SampleFields::default(),
            out: Out::closed(),
        }
    }
}

impl Stat for Heterozygosity {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, self.ext())
    }
    fn site(
        &mut self,
        _s: &SiteView,
        gt: Option<&dyn GenotypeSource>,
        st: Option<&SiteStats>,
    ) -> Result<(), String> {
        let st = stats_of(st)?;
        super::site::diploid_only(st)?;
        let gt = genos(gt)?;
        if st.allele_counts.len() != 2
            || st.haploid
            || st.allele_counts[0] == 0
            || st.allele_counts[1] == 0
        {
            return Ok(());
        }
        let n = st.an() as f64;
        let p = st.allele_counts[0] as f64 / n;
        let e = 1.0 - 2.0 * p * (1.0 - p) * (n / (n - 1.0));
        for (k, (i, _)) in self.samples.iter().enumerate() {
            gt.sample(*i, &mut self.f)?;
            let Some(g) = &self.f.gt else { continue };
            let [Some(a), Some(b)] = g.alleles[..] else {
                continue;
            };
            self.n_sites[k] += 1;
            self.e_hom[k] += e;
            if a == b {
                self.o_hom[k] += 1;
            }
        }
        Ok(())
    }
    fn finish(&mut self) -> io::Result<()> {
        header_row(&mut self.out, "INDV\tO(HOM)\tE(HOM)\tN_SITES\tF");
        for (k, (_, name)) in self.samples.iter().enumerate() {
            if self.n_sites[k] == 0 {
                continue;
            }
            let (o, e, n) = (self.o_hom[k] as f64, self.e_hom[k], self.n_sites[k] as f64);
            let b = &mut self.out.buf;
            b.extend_from_slice(name.as_bytes());
            b.push(b'\t');
            write_u64(b, self.o_hom[k]);
            b.push(b'\t');
            write_fixed(b, e, 1);
            b.push(b'\t');
            write_u64(b, self.n_sites[k]);
            b.push(b'\t');
            write_fixed(b, (o - e) / (n - e), 5);
            b.push(b'\n');
        }
        self.out.row().map_err(io::Error::other)?;
        self.out.commit()
    }
    fn needs_genotypes(&self) -> bool {
        true
    }
    fn ext(&self) -> &'static str {
        "het"
    }
}
