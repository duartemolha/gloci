//! The per-site statistics: one row per record, in input order.

use super::engine::{Out, Stat, StatSink};
use super::fmt::{write_e6, write_f2, write_g6, write_u64};
use crate::filter::SiteView;
use crate::gt::{GenotypeSource, SampleFields};
use crate::hwe::hwe_tails;
use crate::stats::SiteStats;
use std::io;

fn stats_of(st: Option<&SiteStats>) -> Result<&SiteStats, String> {
    st.ok_or_else(|| "internal error: this statistic needs the site statistics".to_string())
}

fn genos(gt: Option<&dyn GenotypeSource>) -> Result<&dyn GenotypeSource, String> {
    gt.ok_or_else(|| "internal error: this statistic needs genotypes".to_string())
}

/// vcftools stops on genotypes of ploidy above two.
pub(super) fn diploid_only(st: &SiteStats) -> Result<(), String> {
    if st.polyploid {
        return Err(
            "polyploid genotypes are not supported (vcftools stops here too); drop or fix those \
             records, or request only depth, quality, --FILTER-summary, --TsTv-summary or \
             --indel-hist"
                .into(),
        );
    }
    Ok(())
}

fn prefix(buf: &mut Vec<u8>, s: &SiteView) {
    buf.extend_from_slice(s.chrom);
    buf.push(b'\t');
    write_u64(buf, s.pos);
}

fn tab(buf: &mut Vec<u8>) {
    buf.push(b'\t');
}

fn header(out: &mut Out, line: &str) -> io::Result<()> {
    out.buf.extend_from_slice(line.as_bytes());
    out.buf.push(b'\n');
    out.row().map_err(io::Error::other)
}

macro_rules! plumbing {
    ($ext:expr, $genos:expr) => {
        fn finish(&mut self) -> io::Result<()> {
            self.out.commit()
        }
        fn needs_genotypes(&self) -> bool {
            $genos
        }
        fn ext(&self) -> &'static str {
            $ext
        }
    };
}

/// `--freq` (`.frq`) and `--counts` (`.frq.count`): per-allele frequency or
/// count over the called alleles of the selected, unmasked samples.
pub struct Freq {
    counts: bool,
    out: Out,
}

impl Freq {
    pub fn frequencies() -> Freq {
        Freq {
            counts: false,
            out: Out::closed(),
        }
    }
    pub fn counts() -> Freq {
        Freq {
            counts: true,
            out: Out::closed(),
        }
    }
}

impl Stat for Freq {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, self.ext())?;
        header(
            &mut self.out,
            if self.counts {
                "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}"
            } else {
                "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}"
            },
        )
    }

    fn site(
        &mut self,
        s: &SiteView,
        _gt: Option<&dyn GenotypeSource>,
        st: Option<&SiteStats>,
    ) -> Result<(), String> {
        let st = stats_of(st)?;
        diploid_only(st)?;
        let b = &mut self.out.buf;
        prefix(b, s);
        tab(b);
        write_u64(b, st.allele_counts.len() as u64);
        tab(b);
        write_u64(b, st.an());
        let mut alts = s.alts.split(|&c| c == b',');
        for (k, &c) in st.allele_counts.iter().enumerate() {
            tab(b);
            b.extend_from_slice(if k == 0 {
                s.ref_
            } else {
                alts.next().unwrap_or(b"")
            });
            b.push(b':');
            if self.counts {
                write_u64(b, c);
            } else {
                write_g6(b, c as f64 / st.an() as f64);
            }
        }
        b.push(b'\n');
        self.out.row()
    }

    fn finish(&mut self) -> io::Result<()> {
        self.out.commit()
    }
    fn needs_genotypes(&self) -> bool {
        true
    }
    fn ext(&self) -> &'static str {
        if self.counts { "frq.count" } else { "frq" }
    }
}

/// `--site-depth` (`.ldepth`): the sum and sum of squares of the FORMAT/DP
/// values of the selected, unmasked samples.
pub struct SiteDepth {
    out: Out,
}

impl SiteDepth {
    pub fn new() -> SiteDepth {
        SiteDepth { out: Out::closed() }
    }
}

impl Default for SiteDepth {
    fn default() -> Self {
        Self::new()
    }
}

impl Stat for SiteDepth {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, "ldepth")?;
        header(&mut self.out, "CHROM\tPOS\tSUM_DEPTH\tSUMSQ_DEPTH")
    }
    fn site(
        &mut self,
        s: &SiteView,
        _gt: Option<&dyn GenotypeSource>,
        st: Option<&SiteStats>,
    ) -> Result<(), String> {
        let st = stats_of(st)?;
        let b = &mut self.out.buf;
        prefix(b, s);
        tab(b);
        write_u64(b, st.dp_sum);
        tab(b);
        write_u64(b, st.dp_sumsq);
        b.push(b'\n');
        self.out.row()
    }
    fn finish(&mut self) -> io::Result<()> {
        self.out.commit()
    }
    fn needs_genotypes(&self) -> bool {
        true
    }
    fn ext(&self) -> &'static str {
        "ldepth"
    }
}

/// `--site-mean-depth` (`.ldepth.mean`): the mean and sample variance of
/// the depths of the samples that have one.
pub struct SiteMeanDepth {
    out: Out,
}

impl SiteMeanDepth {
    pub fn new() -> SiteMeanDepth {
        SiteMeanDepth { out: Out::closed() }
    }
}

impl Default for SiteMeanDepth {
    fn default() -> Self {
        Self::new()
    }
}

impl Stat for SiteMeanDepth {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, self.ext())?;
        header(&mut self.out, "CHROM\tPOS\tMEAN_DEPTH\tVAR_DEPTH")
    }
    fn site(
        &mut self,
        s: &SiteView,
        _gt: Option<&dyn GenotypeSource>,
        st: Option<&SiteStats>,
    ) -> Result<(), String> {
        let st = stats_of(st)?;
        let b = &mut self.out.buf;
        prefix(b, s);
        tab(b);
        write_g6(b, st.mean_dp().unwrap_or(f64::NAN));
        tab(b);
        write_g6(b, st.var_dp().unwrap_or(f64::NAN));
        b.push(b'\n');
        self.out.row()
    }
    plumbing!("ldepth.mean", true);
}

/// `--geno-depth` (`.gdepth`): every selected sample's DP, -1 when it has
/// none or the sample is masked. Columns carry the sample names.
pub struct GenoDepth {
    samples: Vec<(usize, String)>,
    f: SampleFields,
    out: Out,
}

impl GenoDepth {
    /// `samples`: (index in the record, name) of every selected sample.
    pub fn new(samples: Vec<(usize, String)>) -> GenoDepth {
        GenoDepth {
            samples,
            f: SampleFields::default(),
            out: Out::closed(),
        }
    }
}

impl Stat for GenoDepth {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, self.ext())?;
        self.out.buf.extend_from_slice(b"CHROM\tPOS");
        for (_, n) in &self.samples {
            self.out.buf.push(b'\t');
            self.out.buf.extend_from_slice(n.as_bytes());
        }
        self.out.buf.push(b'\n');
        self.out.row().map_err(io::Error::other)
    }
    fn site(
        &mut self,
        s: &SiteView,
        gt: Option<&dyn GenotypeSource>,
        _st: Option<&SiteStats>,
    ) -> Result<(), String> {
        let gt = genos(gt)?;
        prefix(&mut self.out.buf, s);
        for (i, _) in &self.samples {
            gt.sample(*i, &mut self.f)?;
            tab(&mut self.out.buf);
            match self.f.dp {
                Some(d) => write_u64(&mut self.out.buf, u64::from(d)),
                None => self.out.buf.extend_from_slice(b"-1"),
            }
        }
        self.out.buf.push(b'\n');
        self.out.row()
    }
    plumbing!("gdepth", true);
}

/// `--site-quality` (`.lqual`): QUAL, -1 when missing.
pub struct SiteQuality {
    out: Out,
}

impl SiteQuality {
    pub fn new() -> SiteQuality {
        SiteQuality { out: Out::closed() }
    }
}

impl Default for SiteQuality {
    fn default() -> Self {
        Self::new()
    }
}

impl Stat for SiteQuality {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, self.ext())?;
        header(&mut self.out, "CHROM\tPOS\tQUAL")
    }
    fn site(
        &mut self,
        s: &SiteView,
        _gt: Option<&dyn GenotypeSource>,
        _st: Option<&SiteStats>,
    ) -> Result<(), String> {
        let b = &mut self.out.buf;
        prefix(b, s);
        tab(b);
        write_g6(b, s.qual.unwrap_or(-1.0));
        b.push(b'\n');
        self.out.row()
    }
    plumbing!("lqual", false);
}

/// `--missing-site` (`.lmiss`): alleles with data, genotypes the mask
/// removed (counted in samples), missing alleles and their fraction.
pub struct MissingSite {
    out: Out,
}

impl MissingSite {
    pub fn new() -> MissingSite {
        MissingSite { out: Out::closed() }
    }
}

impl Default for MissingSite {
    fn default() -> Self {
        Self::new()
    }
}

impl Stat for MissingSite {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, self.ext())?;
        header(
            &mut self.out,
            "CHR\tPOS\tN_DATA\tN_GENOTYPE_FILTERED\tN_MISS\tF_MISS",
        )
    }
    fn site(
        &mut self,
        s: &SiteView,
        _gt: Option<&dyn GenotypeSource>,
        st: Option<&SiteStats>,
    ) -> Result<(), String> {
        let st = stats_of(st)?;
        diploid_only(st)?;
        let b = &mut self.out.buf;
        prefix(b, s);
        tab(b);
        write_u64(b, st.data_alleles);
        tab(b);
        write_u64(b, st.filtered);
        tab(b);
        write_u64(b, st.missing_alleles);
        tab(b);
        write_g6(b, st.missing_alleles as f64 / st.data_alleles as f64);
        b.push(b'\n');
        self.out.row()
    }
    plumbing!("lmiss", true);
}

/// `--hardy` (`.hwe`): observed and expected genotype counts, the chi-square
/// statistic and the exact test's p-values. Only biallelic sites without
/// haploid genotypes get a row.
pub struct Hardy {
    out: Out,
}

impl Hardy {
    pub fn new() -> Hardy {
        Hardy { out: Out::closed() }
    }
}

impl Default for Hardy {
    fn default() -> Self {
        Self::new()
    }
}

impl Stat for Hardy {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, self.ext())?;
        header(
            &mut self.out,
            "CHR\tPOS\tOBS(HOM1/HET/HOM2)\tE(HOM1/HET/HOM2)\tChiSq_HWE\tP_HWE\tP_HET_DEFICIT\tP_HET_EXCESS",
        )
    }
    fn site(
        &mut self,
        s: &SiteView,
        _gt: Option<&dyn GenotypeSource>,
        st: Option<&SiteStats>,
    ) -> Result<(), String> {
        let st = stats_of(st)?;
        diploid_only(st)?;
        let Some((hom1, het, hom2)) = st.hwe_counts else {
            return Ok(());
        };
        if st.haploid {
            return Ok(());
        }
        // Expected counts use the allele frequencies over every called
        // allele (a half-called genotype adds its allele) but the number of
        // fully called genotypes, as vcftools does.
        let n = (hom1 + het + hom2) as f64;
        let an = st.an() as f64;
        let p = st.allele_counts[0] as f64 / an;
        let q = 1.0 - p;
        // This operand order reproduces vcftools' rounding of exact halves
        // (2.025 prints as 2.03).
        let e = [p * p * n, 2.0 * p * q * n, q * q * n];
        let o = [hom1 as f64, het as f64, hom2 as f64];
        let chisq: f64 = (0..3).map(|i| (o[i] - e[i]) * (o[i] - e[i]) / e[i]).sum();
        let t = hwe_tails(het, hom1, hom2);
        let b = &mut self.out.buf;
        prefix(b, s);
        tab(b);
        write_u64(b, hom1);
        b.push(b'/');
        write_u64(b, het);
        b.push(b'/');
        write_u64(b, hom2);
        tab(b);
        for (i, v) in e.iter().enumerate() {
            if i > 0 {
                b.push(b'/');
            }
            write_f2(b, *v);
        }
        for v in [chisq, t.two_sided, t.het_deficit, t.het_excess] {
            tab(b);
            write_e6(b, v);
        }
        b.push(b'\n');
        self.out.row()
    }
    plumbing!("hwe", true);
}

/// `--singletons` (`.singletons`): alleles seen once (`S`) or twice in one
/// homozygote (`D`), with the carrying sample.
pub struct Singletons {
    samples: Vec<(usize, String)>,
    f: SampleFields,
    out: Out,
}

impl Singletons {
    /// `samples`: (index in the record, name) of every selected sample.
    pub fn new(samples: Vec<(usize, String)>) -> Singletons {
        Singletons {
            samples,
            f: SampleFields::default(),
            out: Out::closed(),
        }
    }
}

impl Stat for Singletons {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, self.ext())?;
        header(
            &mut self.out,
            "CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV",
        )
    }
    fn site(
        &mut self,
        s: &SiteView,
        gt: Option<&dyn GenotypeSource>,
        st: Option<&SiteStats>,
    ) -> Result<(), String> {
        let st = stats_of(st)?;
        let gt = genos(gt)?;
        diploid_only(st)?;
        for (k, &count) in st.allele_counts.iter().enumerate() {
            if count == 0 || count > 2 {
                continue;
            }
            // Find the carrier: the first sample with the allele for a
            // singleton, a homozygote for a doubleton.
            let mut hit = None;
            for (idx, (i, _)) in self.samples.iter().enumerate() {
                gt.sample(*i, &mut self.f)?;
                let Some(g) = &self.f.gt else { continue };
                let copies = g.alleles.iter().filter(|&&a| a == Some(k as u32)).count();
                if copies == 0 {
                    continue;
                }
                if count == 1 || copies == 2 {
                    hit = Some(idx);
                    break;
                }
            }
            let Some(idx) = hit else { continue };
            let b = &mut self.out.buf;
            prefix(b, s);
            tab(b);
            b.push(if count == 1 { b'S' } else { b'D' });
            tab(b);
            b.extend_from_slice(if k == 0 {
                s.ref_
            } else {
                s.alts.split(|&c| c == b',').nth(k - 1).unwrap_or(b"")
            });
            tab(b);
            b.extend_from_slice(self.samples[idx].1.as_bytes());
            b.push(b'\n');
            self.out.row()?;
        }
        Ok(())
    }
    plumbing!("singletons", true);
}
