//! Relatedness of every pair of selected samples as vcftools 0.1.16
//! computes it: Manichaikul et al.'s kinship φ (`--relatedness2`) and Yang
//! et al.'s unadjusted Ajk (`--relatedness`), over the biallelic sites where
//! every selected, unmasked sample is diploid. A genotype counts only when
//! both its alleles are called and it is not masked.
//!
//! The pairs of a sample with itself and with each later sample live in a
//! triangle, row i holding (i, i), (i, i + 1), …, (i, N − 1): 8 bytes a pair
//! for φ's two counts, 12 for Ajk's sum and count. Sites are read into
//! blocks of 64. A block keeps, per sample, bit masks of its sites (called;
//! or heterozygous, homozygous REF, homozygous ALT), whose popcounts give
//! the counts; for Ajk it also keeps each site's centred genotypes x − 2f,
//! which a pair adds site by site, in input order, as vcftools does. A full
//! block is added to the triangle with its rows split between threads, each
//! pair taken by one thread, so every thread count gives the same bits.

use super::diploid;
use crate::compare::Gt;
use crate::filter::SiteView;
use crate::gt::{GenotypeSource, SampleFields};
use crate::stats::{Out, SiteStats, Stat, StatSink};
use gloci_core::cfmt::{write_g6, write_u64};
use std::io;
use std::num::NonZeroUsize;
use std::ops::Range;

/// Sites per block: one bit of a 64-bit mask each.
const BLOCK: usize = 64;

/// The pair-site steps below which a thread is not worth starting.
const MIN_WORK: usize = 1 << 18;

/// The pairs of `n` samples with themselves and each later sample, row by
/// row.
#[derive(Clone, Copy, Debug)]
struct Triangle {
    n: usize,
}

impl Triangle {
    /// `n(n + 1)/2`, or `None` when it does not fit in memory sizes.
    fn pairs(n: usize) -> Option<usize> {
        n.checked_mul(n.checked_add(1)?).map(|x| x / 2)
    }

    /// The index of (i, i): the rows before i hold n, n − 1, … pairs.
    fn row(self, i: usize) -> usize {
        i * self.n - i * i.saturating_sub(1) / 2
    }

    /// The index of the pair of samples i and j, in either order.
    fn at(self, i: usize, j: usize) -> usize {
        let (a, b) = if i <= j { (i, j) } else { (j, i) };
        self.row(a) + (b - a)
    }

    /// Up to `parts` ranges of rows, in order, holding about as many pairs
    /// each.
    fn split(self, parts: usize) -> Vec<Range<usize>> {
        let total = self.row(self.n);
        let mut out = Vec::with_capacity(parts);
        let mut start = 0;
        for p in 1..=parts {
            let target = total * p / parts;
            let mut end = start;
            while end < self.n && self.row(end) < target {
                end += 1;
            }
            if p == parts {
                end = self.n;
            }
            if end > start {
                out.push(start..end);
                start = end;
            }
        }
        out
    }
}

/// The matrix of `n(n + 1)/2` zeroes, or the message for a run with too
/// many samples to hold it.
fn matrix<T: Clone + Default>(
    flag: &str,
    n: usize,
    per_pair: usize,
    pairs: Option<usize>,
) -> Result<Vec<T>, String> {
    let too_many = || {
        let bytes = pairs.map_or(u128::MAX, |p| p as u128 * per_pair as u128);
        format!(
            "{flag} holds {per_pair} bytes for each pair of the {n} selected samples, {bytes} \
             bytes in all, more than this machine can allocate; select fewer samples with \
             --samples or --samples-file"
        )
    };
    let pairs = pairs.ok_or_else(too_many)?;
    let mut v = Vec::new();
    v.try_reserve_exact(pairs).map_err(|_| too_many())?;
    v.resize(pairs, T::default());
    Ok(v)
}

/// A vector of `n` zeroes with exactly that capacity.
fn zeroes<T: Clone + Default>(n: usize) -> Vec<T> {
    let mut v = Vec::with_capacity(n);
    v.resize(n, T::default());
    v
}

/// The threads worth starting for `work` pair-site steps, at least
/// `min_work` each.
fn threads_for(work: usize, threads: NonZeroUsize, min_work: usize) -> usize {
    (work / min_work.max(1)).clamp(1, threads.get())
}

/// How a statistic splits its blocks between threads.
#[derive(Clone, Copy, Debug)]
struct Split {
    threads: NonZeroUsize,
    /// The pair-site steps a thread must have to be started ([`MIN_WORK`];
    /// tests lower it to split small blocks).
    min_work: usize,
    /// The most row ranges a block was split into so far.
    most: usize,
}

impl Split {
    fn new(threads: NonZeroUsize) -> Split {
        Split {
            threads,
            min_work: MIN_WORK,
            most: 0,
        }
    }

    /// The threads for a block of `work` pair-site steps.
    fn threads(&self, work: usize) -> usize {
        threads_for(work, self.threads, self.min_work)
    }
}

/// Run `f` over contiguous ranges of rows of `tri`, each with the slices of
/// `a` and `b` that hold its pairs, on up to `threads` threads; the number
/// of ranges.
fn split_rows<A: Send, B: Send>(
    tri: Triangle,
    a: &mut [A],
    b: &mut [B],
    threads: usize,
    f: impl Fn(Range<usize>, &mut [A], &mut [B]) + Sync,
) -> usize {
    let parts = tri.split(threads);
    let used = parts.len();
    if used <= 1 {
        f(0..tri.n, a, b);
        return used;
    }
    std::thread::scope(|scope| {
        let (mut a, mut b) = (a, b);
        let f = &f;
        let last = parts.len() - 1;
        for (k, rows) in parts.into_iter().enumerate() {
            let len = tri.row(rows.end) - tri.row(rows.start);
            let (ha, ta) = std::mem::take(&mut a).split_at_mut(len);
            let (hb, tb) = std::mem::take(&mut b).split_at_mut(len);
            (a, b) = (ta, tb);
            if k == last {
                f(rows, ha, hb);
            } else {
                scope.spawn(move || f(rows, ha, hb));
            }
        }
    });
    used
}

/// A sample's ALT allele count (0, 1 or 2), `None` when its genotype is
/// missing, half-called or masked.
fn dosage(g: Gt) -> Result<Option<usize>, String> {
    match g {
        Gt::Dip(a, b) if a <= 1 && b <= 1 => Ok(Some((a + b) as usize)),
        Gt::Dip(a, b) => Err(format!(
            "genotype allele {} but the site lists 2 alleles",
            a.max(b)
        )),
        _ => Ok(None),
    }
}

/// The site statistics of a site the relatedness statistics use (biallelic,
/// every selected and unmasked sample diploid), `None` for any other. A
/// genotype of ploidy above two among them stops the run, as vcftools stops
/// when it reads one; vcftools drops a multi-allelic site before.
fn site_used(st: Option<&SiteStats>) -> Result<Option<&SiteStats>, String> {
    let st = st.ok_or("internal error: this statistic needs the site statistics")?;
    if st.allele_counts.len() != 2 || !diploid(st)? {
        return Ok(None);
    }
    Ok(Some(st))
}

fn genotypes(gt: Option<&dyn GenotypeSource>) -> Result<&dyn GenotypeSource, String> {
    gt.ok_or_else(|| "internal error: this statistic needs genotypes".to_string())
}

/// The counts are 32-bit: the sites used must fit.
fn count_site(used: &mut u64, flag: &str) -> Result<(), String> {
    if *used == u64::from(u32::MAX) {
        return Err(format!(
            "{flag} counts the sites of each pair in 32 bits, and this input has more than \
             {} sites it uses; use fewer sites, with --region or the site filters",
            u32::MAX
        ));
    }
    *used += 1;
    Ok(())
}

/// A count as vcftools prints a double (`%g`): a whole number below a
/// million as it is, a larger one with an exponent (`1e+06`).
fn write_count(o: &mut Vec<u8>, c: u32) {
    if c < 1_000_000 {
        write_u64(o, u64::from(c));
    } else {
        write_g6(o, f64::from(c));
    }
}

fn header(out: &mut Out, line: &str) -> io::Result<()> {
    out.buf.extend_from_slice(line.as_bytes());
    out.buf.push(b'\n');
    out.row().map_err(io::Error::other)
}

/// `--relatedness2` (`.relatedness2`), Manichaikul et al. (2010): for every
/// ordered pair of selected samples, itself included, in header order, the
/// sites where both are heterozygous (`N_AaAa`), where they are homozygous
/// for different alleles (`N_AAaa`), each one's heterozygous sites (`N1_Aa`,
/// `N2_Aa`) and φ = (N_AaAa − 2 N_AAaa) / (N1_Aa + N2_Aa). A pair neither of
/// which is heterozygous anywhere is `nan`, or `-inf` when they are opposite
/// homozygotes somewhere. Counts print as vcftools prints a double (`1e+06`
/// from a million).
pub struct RelatednessManichaikul {
    samples: Vec<(usize, String)>,
    split: Split,
    tri: Triangle,
    /// Per pair: both heterozygous, and opposite homozygotes.
    both_het: Vec<u32>,
    opposite: Vec<u32>,
    /// The block: per sample, a bit per site for each called genotype.
    het: Vec<u64>,
    hom_ref: Vec<u64>,
    hom_alt: Vec<u64>,
    filled: usize,
    used: u64,
    scratch: SampleFields,
    out: Out,
}

impl RelatednessManichaikul {
    /// Over the `samples` (index in the record and name, in header order),
    /// adding blocks on up to `threads` threads; an error when the matrix
    /// of their pairs cannot be allocated.
    pub fn new(
        samples: Vec<(usize, String)>,
        threads: NonZeroUsize,
    ) -> Result<RelatednessManichaikul, String> {
        let n = samples.len();
        let pairs = Triangle::pairs(n);
        Ok(RelatednessManichaikul {
            both_het: matrix("--relatedness2", n, 8, pairs)?,
            opposite: matrix("--relatedness2", n, 8, pairs)?,
            het: zeroes(n),
            hom_ref: zeroes(n),
            hom_alt: zeroes(n),
            tri: Triangle { n },
            samples,
            split: Split::new(threads),
            filled: 0,
            used: 0,
            scratch: SampleFields::default(),
            out: Out::closed(),
        })
    }

    /// Add the block to the matrix and empty it.
    fn flush(&mut self) {
        if self.filled == 0 {
            return;
        }
        let (tri, het, hom_ref, hom_alt) = (self.tri, &self.het, &self.hom_ref, &self.hom_alt);
        let work = tri.row(tri.n).saturating_mul(4);
        let used = split_rows(
            tri,
            &mut self.both_het,
            &mut self.opposite,
            self.split.threads(work),
            |rows, both, opp| {
                let mut at = 0;
                for i in rows {
                    let len = tri.n - i;
                    let (h, r, a) = (het[i], hom_ref[i], hom_alt[i]);
                    if h | r | a != 0 {
                        let both = &mut both[at..at + len];
                        let opp = &mut opp[at..at + len];
                        let j = i..tri.n;
                        for (k, j) in j.enumerate() {
                            both[k] += (h & het[j]).count_ones();
                            opp[k] += (r & hom_alt[j]).count_ones() + (a & hom_ref[j]).count_ones();
                        }
                    }
                    at += len;
                }
            },
        );
        self.split.most = self.split.most.max(used);
        self.het.fill(0);
        self.hom_ref.fill(0);
        self.hom_alt.fill(0);
        self.filled = 0;
    }

    /// `N_AaAa` and `N_AAaa` of the selected samples i and j (positions in
    /// the selection). For tests.
    #[doc(hidden)]
    pub fn counts(&self, i: usize, j: usize) -> (u32, u32) {
        let p = self.tri.at(i, j);
        (self.both_het[p], self.opposite[p])
    }

    /// φ of the selected samples i and j. For tests.
    #[doc(hidden)]
    pub fn phi(&self, i: usize, j: usize) -> f64 {
        let (both, opp) = self.counts(i, j);
        let het = |k: usize| f64::from(self.both_het[self.tri.row(k)]);
        (f64::from(both) - 2.0 * f64::from(opp)) / (het(i) + het(j))
    }

    /// Split a block between threads from `steps` pair-site steps a
    /// thread instead of the default. For tests.
    #[doc(hidden)]
    pub fn set_min_work(&mut self, steps: usize) {
        self.split.min_work = steps;
    }

    /// The most threads that took rows of one block so far. For tests.
    #[doc(hidden)]
    pub fn most_threads(&self) -> usize {
        self.split.most
    }

    /// The bytes of the matrix and of the block. For tests.
    #[doc(hidden)]
    pub fn held_bytes(&self) -> (u64, u64) {
        let m = (self.both_het.capacity() + self.opposite.capacity()) * size_of::<u32>();
        let b = (self.het.capacity() + self.hom_ref.capacity() + self.hom_alt.capacity())
            * size_of::<u64>();
        (m as u64, b as u64)
    }
}

impl Stat for RelatednessManichaikul {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, self.ext())?;
        header(
            &mut self.out,
            "INDV1\tINDV2\tN_AaAa\tN_AAaa\tN1_Aa\tN2_Aa\tRELATEDNESS_PHI",
        )
    }

    fn site(
        &mut self,
        _s: &SiteView,
        gt: Option<&dyn GenotypeSource>,
        st: Option<&SiteStats>,
    ) -> Result<(), String> {
        if site_used(st)?.is_none() {
            return Ok(());
        }
        let gt = genotypes(gt)?;
        count_site(&mut self.used, "--relatedness2")?;
        let bit = 1u64 << self.filled;
        for (k, (i, _)) in self.samples.iter().enumerate() {
            match dosage(gt.gt(*i, &mut self.scratch)?)? {
                Some(0) => self.hom_ref[k] |= bit,
                Some(1) => self.het[k] |= bit,
                Some(_) => self.hom_alt[k] |= bit,
                None => {}
            }
        }
        self.filled += 1;
        if self.filled == BLOCK {
            self.flush();
        }
        Ok(())
    }

    fn finish(&mut self) -> io::Result<()> {
        self.flush();
        for (i, (_, a)) in self.samples.iter().enumerate() {
            let het_i = self.both_het[self.tri.row(i)];
            for (j, (_, b)) in self.samples.iter().enumerate() {
                let (both, opp) = self.counts(i, j);
                let het_j = self.both_het[self.tri.row(j)];
                let o = &mut self.out.buf;
                o.extend_from_slice(a.as_bytes());
                o.push(b'\t');
                o.extend_from_slice(b.as_bytes());
                for c in [both, opp, het_i, het_j] {
                    o.push(b'\t');
                    write_count(o, c);
                }
                o.push(b'\t');
                let phi = (f64::from(both) - 2.0 * f64::from(opp))
                    / (f64::from(het_i) + f64::from(het_j));
                write_g6(o, phi);
                o.push(b'\n');
                self.out.row().map_err(io::Error::other)?;
            }
        }
        self.out.commit()
    }

    fn needs_genotypes(&self) -> bool {
        true
    }

    fn reads_genotypes(&self) -> bool {
        true
    }

    fn ext(&self) -> &'static str {
        "relatedness2"
    }
}

/// `--relatedness` (`.relatedness`), Yang et al. (2010) equation 6,
/// unadjusted: over the sites whose ALT frequency f among the called
/// alleles (a half-called genotype's included) is strictly between 0 and 1
/// (beyond a double's epsilon), for each pair of selected samples in header
/// order, i before j, the mean over the sites where both are called of
/// (x_i − 2f)(x_j − 2f) / (2f(1 − f)), x the ALT allele count; for a sample
/// with itself 1 + the mean of (x² − (1 + 2f)x + 2f²) / (2f(1 − f)). A pair
/// called together at no site is `nan`.
pub struct RelatednessYang {
    samples: Vec<(usize, String)>,
    split: Split,
    tri: Triangle,
    /// Per pair: the sum and the sites called at both.
    sums: Vec<f64>,
    counts: Vec<u32>,
    /// The block: per site, each sample's x − 2f (0 when not called) ...
    centred: Vec<f64>,
    /// ... and 1/(2f(1 − f)); per sample, a bit per site where it is
    /// called.
    div: Vec<f64>,
    called: Vec<u64>,
    filled: usize,
    used: u64,
    scratch: SampleFields,
    out: Out,
}

impl RelatednessYang {
    /// Over the `samples` (index in the record and name, in header order),
    /// adding blocks on up to `threads` threads; an error when the matrix
    /// of their pairs cannot be allocated.
    pub fn new(
        samples: Vec<(usize, String)>,
        threads: NonZeroUsize,
    ) -> Result<RelatednessYang, String> {
        let n = samples.len();
        let pairs = Triangle::pairs(n);
        let centred = n
            .checked_mul(BLOCK)
            .ok_or_else(|| format!("{n} samples are too many for --relatedness"))?;
        Ok(RelatednessYang {
            sums: matrix("--relatedness", n, 12, pairs)?,
            counts: matrix("--relatedness", n, 12, pairs)?,
            centred: zeroes(centred),
            div: zeroes(BLOCK),
            called: zeroes(n),
            tri: Triangle { n },
            samples,
            split: Split::new(threads),
            filled: 0,
            used: 0,
            scratch: SampleFields::default(),
            out: Out::closed(),
        })
    }

    /// Add the block to the matrix and empty it: per pair, the count of
    /// sites called at both, and the terms of those sites in their order.
    fn flush(&mut self) {
        if self.filled == 0 {
            return;
        }
        let (tri, called, centred, div) = (self.tri, &self.called, &self.centred, &self.div);
        let work = tri.row(tri.n).saturating_mul(self.filled);
        let used = split_rows(
            tri,
            &mut self.sums,
            &mut self.counts,
            self.split.threads(work),
            |rows, sums, counts| {
                let n = tri.n;
                let mut at = 0;
                for i in rows {
                    let len = n - i;
                    let ci = called[i];
                    if ci != 0 {
                        // (i, i) took its terms as the sites were read
                        let sums = &mut sums[at + 1..at + len];
                        let counts = &mut counts[at + 1..at + len];
                        for (c, &cj) in counts.iter_mut().zip(&called[i + 1..]) {
                            *c += (ci & cj).count_ones();
                        }
                        let mut bits = ci;
                        while bits != 0 {
                            let s = bits.trailing_zeros() as usize;
                            bits &= bits - 1;
                            let row = &centred[s * n..(s + 1) * n];
                            let (xi, d) = (row[i], div[s]);
                            // A sample not called is 0 here, and adding the
                            // ±0 its term gives leaves the sum as it is
                            // (it starts at +0 and is never -0).
                            for (sum, &xj) in sums.iter_mut().zip(&row[i + 1..]) {
                                *sum += xi * xj * d;
                            }
                        }
                    }
                    at += len;
                }
            },
        );
        self.split.most = self.split.most.max(used);
        self.called.fill(0);
        self.filled = 0;
    }

    /// Ajk of the selected samples i ≤ j (positions in the selection). For
    /// tests.
    #[doc(hidden)]
    pub fn ajk(&self, i: usize, j: usize) -> f64 {
        let p = self.tri.at(i, j);
        let mean = self.sums[p] / f64::from(self.counts[p]);
        if i == j { 1.0 + mean } else { mean }
    }

    /// Split a block between threads from `steps` pair-site steps a
    /// thread instead of the default. For tests.
    #[doc(hidden)]
    pub fn set_min_work(&mut self, steps: usize) {
        self.split.min_work = steps;
    }

    /// The most threads that took rows of one block so far. For tests.
    #[doc(hidden)]
    pub fn most_threads(&self) -> usize {
        self.split.most
    }

    /// The bytes of the matrix and of the block. For tests.
    #[doc(hidden)]
    pub fn held_bytes(&self) -> (u64, u64) {
        let m = self.sums.capacity() * size_of::<f64>() + self.counts.capacity() * size_of::<u32>();
        let b = (self.centred.capacity() + self.div.capacity()) * size_of::<f64>()
            + self.called.capacity() * size_of::<u64>();
        (m as u64, b as u64)
    }
}

impl Stat for RelatednessYang {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, self.ext())?;
        header(&mut self.out, "INDV1\tINDV2\tRELATEDNESS_AJK")
    }

    fn site(
        &mut self,
        _s: &SiteView,
        gt: Option<&dyn GenotypeSource>,
        st: Option<&SiteStats>,
    ) -> Result<(), String> {
        let Some(st) = site_used(st)? else {
            return Ok(());
        };
        if st.called_alleles == 0 {
            return Ok(());
        }
        // vcftools' operations, in its order.
        let freq = st.allele_counts[1] as f64 / st.called_alleles as f64;
        if freq <= f64::EPSILON || freq >= 1.0 - f64::EPSILON {
            return Ok(());
        }
        let gt = genotypes(gt)?;
        count_site(&mut self.used, "--relatedness")?;
        let div = 1.0 / (2.0 * freq * (1.0 - freq));
        let centred = [0.0, 1.0, 2.0].map(|x: f64| x - 2.0 * freq);
        let own = [0.0, 1.0, 2.0]
            .map(|x: f64| (x * x - (1.0 + 2.0 * freq) * x + 2.0 * freq * freq) * div);
        let (s, n) = (self.filled, self.tri.n);
        let bit = 1u64 << s;
        let row = &mut self.centred[s * n..(s + 1) * n];
        for (k, (i, _)) in self.samples.iter().enumerate() {
            match dosage(gt.gt(*i, &mut self.scratch)?)? {
                Some(x) => {
                    row[k] = centred[x];
                    self.called[k] |= bit;
                    let d = self.tri.row(k);
                    self.sums[d] += own[x];
                    self.counts[d] += 1;
                }
                None => row[k] = 0.0,
            }
        }
        self.div[s] = div;
        self.filled += 1;
        if self.filled == BLOCK {
            self.flush();
        }
        Ok(())
    }

    fn finish(&mut self) -> io::Result<()> {
        self.flush();
        for (i, (_, a)) in self.samples.iter().enumerate() {
            for (j, (_, b)) in self.samples.iter().enumerate().skip(i) {
                let v = self.ajk(i, j);
                let o = &mut self.out.buf;
                o.extend_from_slice(a.as_bytes());
                o.push(b'\t');
                o.extend_from_slice(b.as_bytes());
                o.push(b'\t');
                write_g6(o, v);
                o.push(b'\n');
                self.out.row().map_err(io::Error::other)?;
            }
        }
        self.out.commit()
    }

    fn needs_genotypes(&self) -> bool {
        true
    }

    fn reads_genotypes(&self) -> bool {
        true
    }

    fn ext(&self) -> &'static str {
        "relatedness"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn triangle_rows_and_splits() {
        let t = Triangle { n: 5 };
        assert_eq!((t.row(0), t.row(1), t.row(2), t.row(5)), (0, 5, 9, 15));
        assert_eq!((t.at(1, 3), t.at(3, 1), t.at(4, 4)), (7, 7, 14));
        assert_eq!(Triangle::pairs(5), Some(15));
        assert_eq!(Triangle::pairs(usize::MAX), None);
        for parts in 1..8 {
            let s = t.split(parts);
            assert_eq!(s.first().unwrap().start, 0);
            assert_eq!(s.last().unwrap().end, 5);
            assert!(
                s.windows(2)
                    .all(|w| w[0].end == w[1].start && w[0].start < w[0].end)
            );
            assert!(s.len() <= parts);
        }
        // pairs, not rows, are shared: the first rows are the longest
        let big = Triangle { n: 1000 };
        let s = big.split(4);
        assert_eq!(s.len(), 4);
        assert!(s[0].len() < s[3].len());
        assert_eq!(Triangle { n: 0 }.split(3), Vec::<Range<usize>>::new());
    }

    #[test]
    fn counts_print_as_doubles() {
        for c in [0, 1, 999_999, 1_000_000, 1_234_567, 99_999_999, u32::MAX] {
            let (mut a, mut b) = (Vec::new(), Vec::new());
            write_count(&mut a, c);
            write_g6(&mut b, f64::from(c));
            assert_eq!(a, b, "{c}");
        }
    }

    #[test]
    fn work_decides_the_threads() {
        let t = NonZeroUsize::new(8).unwrap();
        assert_eq!(threads_for(0, t, MIN_WORK), 1);
        assert_eq!(threads_for(MIN_WORK * 3, t, MIN_WORK), 3);
        assert_eq!(threads_for(usize::MAX, t, MIN_WORK), 8);
        assert_eq!(threads_for(5, t, 1), 5);
        assert_eq!(threads_for(5, t, 0), 5);
    }
}
