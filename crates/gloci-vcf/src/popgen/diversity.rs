//! Nucleotide diversity per site (`--site-pi`) and in windows
//! (`--window-pi`), and Tajima's D in bins (`--TajimaD`), as vcftools 0.1.16
//! computes them. Products vcftools takes in 32-bit integers are taken in
//! u64 here, so they never wrap.

use super::{SortedInput, diploid};
use crate::filter::SiteView;
use crate::gt::GenotypeSource;
use crate::stats::{Out, SiteStats, Stat, StatSink};
use gloci_core::cfmt::{write_g6, write_u64};
use std::collections::VecDeque;
use std::io;

fn stats_of(st: Option<&SiteStats>) -> Result<&SiteStats, String> {
    st.ok_or_else(|| "internal error: this statistic needs the site statistics".to_string())
}

/// `Σ c (n − c)` over the allele counts `c` of the `n` called alleles: the
/// ordered pairs of called alleles that differ.
fn mismatches(st: &SiteStats) -> u64 {
    let n = st.called_alleles;
    st.allele_counts.iter().map(|&c| c * (n - c)).sum()
}

fn header(out: &mut Out, line: &str) -> io::Result<()> {
    out.buf.extend_from_slice(line.as_bytes());
    out.buf.push(b'\n');
    out.row().map_err(io::Error::other)
}

/// `--site-pi` (`.sites.pi`): per diploid site, the chance that two called
/// alleles drawn without replacement differ, `Σ c (n − c) / (n (n − 1))`;
/// `nan` with fewer than two called alleles. A site where a selected,
/// unmasked sample is not diploid has no row.
pub struct SitePi {
    out: Out,
}

impl SitePi {
    pub fn new() -> SitePi {
        SitePi { out: Out::closed() }
    }
}

impl Default for SitePi {
    fn default() -> Self {
        Self::new()
    }
}

impl Stat for SitePi {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, self.ext())?;
        header(&mut self.out, "CHROM\tPOS\tPI")
    }
    fn site(
        &mut self,
        s: &SiteView,
        _gt: Option<&dyn GenotypeSource>,
        st: Option<&SiteStats>,
    ) -> Result<(), String> {
        let st = stats_of(st)?;
        if !diploid(st) {
            return Ok(());
        }
        let n = st.called_alleles;
        let pi = if n < 2 {
            f64::NAN
        } else {
            mismatches(st) as f64 / (n * (n - 1)) as f64
        };
        let b = &mut self.out.buf;
        b.extend_from_slice(s.chrom);
        b.push(b'\t');
        write_u64(b, s.pos);
        b.push(b'\t');
        write_g6(b, pi);
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
        "sites.pi"
    }
}

/// The counts of one window of `--window-pi`.
#[derive(Clone, Copy, Debug, Default)]
struct Window {
    /// Polymorphic sites placed in the window. Each has a called allele other
    /// than REF, so this is also vcftools' N_VARIANTS.
    sites: u64,
    /// Their pairs of called alleles, `Σ n (n − 1)`.
    pairs: u64,
    /// Their differing pairs, `Σ Σ c (n − c)`.
    mismatches: u64,
}

/// `--window-pi SIZE` with `--window-pi-step STEP` (`.windowed.pi`):
/// windows `[k·STEP + 1, k·STEP + SIZE]` of each chromosome. A diploid
/// polymorphic site adds its pairs of called alleles and the pairs that
/// differ to every window holding it; every other position of a window is
/// counted as monomorphic with all `N (N − 1)` pairs called, N = 2 × the
/// selected samples. π is the differing pairs over all pairs. Windows
/// without a polymorphic site have no row.
///
/// Input must be sorted by chromosome and position. Only windows that hold
/// a polymorphic site are open, and each is written once a polymorphic site
/// lies past it or its chromosome ends, so at most `⌈SIZE/STEP⌉ + 1` are
/// open.
pub struct WindowPi {
    size: u64,
    step: u64,
    sorted: SortedInput,
    chrom: Vec<u8>,
    /// The index `k` of the first open window.
    base: u64,
    open: VecDeque<Window>,
    /// The most windows open at once.
    peak: usize,
    /// 2 × the selected samples.
    n_chr: u64,
    out: Out,
}

impl WindowPi {
    /// `size` and `step` must be positive, `step` at most `size`.
    pub fn new(size: u64, step: u64) -> WindowPi {
        assert!(
            size > 0 && step > 0 && step <= size,
            "a pi window needs 0 < step <= size"
        );
        WindowPi {
            size,
            step,
            sorted: SortedInput::new("--window-pi"),
            chrom: Vec::new(),
            base: 0,
            open: VecDeque::new(),
            peak: 0,
            n_chr: 0,
            out: Out::closed(),
        }
    }

    /// The most windows held open at once so far: at most
    /// `⌈size/step⌉ + 1`. For tests.
    #[doc(hidden)]
    pub fn peak_open(&self) -> usize {
        self.peak
    }

    /// Write the open windows below index `end` (all of them for
    /// `u64::MAX`) and drop them.
    fn close_before(&mut self, end: u64) -> Result<(), String> {
        while self.base < end {
            let Some(w) = self.open.pop_front() else {
                break;
            };
            let k = self.base;
            self.base += 1;
            // vcftools' unsigned arithmetic: with more records than positions
            // in a window (several at one position) the monomorphic count
            // goes below zero and wraps, and the sum wraps back.
            let mono = self.size.wrapping_sub(w.sites);
            let all_pairs = self.n_chr * self.n_chr.saturating_sub(1);
            let pairs = w.pairs.wrapping_add(mono.wrapping_mul(all_pairs));
            let pi = w.mismatches as f64 / pairs as f64;
            let b = &mut self.out.buf;
            b.extend_from_slice(&self.chrom);
            b.push(b'\t');
            write_u64(b, k * self.step + 1);
            b.push(b'\t');
            write_u64(b, k * self.step + self.size);
            b.push(b'\t');
            write_u64(b, w.sites);
            b.push(b'\t');
            write_g6(b, pi);
            b.push(b'\n');
            self.out.row()?;
        }
        Ok(())
    }
}

impl Stat for WindowPi {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, self.ext())?;
        header(&mut self.out, "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI")
    }
    fn site(
        &mut self,
        s: &SiteView,
        _gt: Option<&dyn GenotypeSource>,
        st: Option<&SiteStats>,
    ) -> Result<(), String> {
        let st = stats_of(st)?;
        if self.sorted.check(s)? {
            self.close_before(u64::MAX)?;
            self.chrom.clear();
            self.chrom.extend_from_slice(s.chrom);
        }
        self.n_chr = 2 * st.samples;
        if !diploid(st) {
            return Ok(());
        }
        let diffs = mismatches(st);
        if diffs == 0 {
            return Ok(());
        }
        let n = st.called_alleles;
        // The windows holding the position: k·step + 1 <= pos <= k·step + size.
        let first = if s.pos > self.size {
            (s.pos - self.size).div_ceil(self.step)
        } else {
            0
        };
        let last = s.pos.div_ceil(self.step);
        self.close_before(first)?;
        if self.open.is_empty() {
            self.base = first;
        }
        let have = self.base + self.open.len() as u64;
        if last > have {
            self.open
                .resize(self.open.len() + (last - have) as usize, Window::default());
            self.peak = self.peak.max(self.open.len());
        }
        for w in self
            .open
            .range_mut((first - self.base) as usize..(last - self.base) as usize)
        {
            w.sites += 1;
            w.pairs += n * (n - 1);
            w.mismatches += diffs;
        }
        Ok(())
    }
    fn finish(&mut self) -> io::Result<()> {
        self.close_before(u64::MAX).map_err(io::Error::other)?;
        self.out.commit()
    }
    fn needs_genotypes(&self) -> bool {
        true
    }
    fn ext(&self) -> &'static str {
        "windowed.pi"
    }
}

/// `--TajimaD BIN` (`.Tajima.D`): Tajima's D per bin of each chromosome,
/// over the biallelic diploid sites. A position falls in bin
/// `⌊POS × (1/BIN)⌋`, the product taken in floating point as vcftools does
/// (POS 49 with BIN 49 is in bin 0). Every such site extends its
/// chromosome's bins; a polymorphic one adds 1 to S and `p (1 − p)` to its
/// bin, p being the REF frequency among the called alleles. Rows run from
/// the first bin with S > 0 to the chromosome's last bin, empty bins as
/// `0 nan`. The sample size is n = 2 × the selected samples throughout,
/// whatever is missing at a site.
///
/// Input must be sorted by chromosome and position; one bin is held.
pub struct TajimaD {
    bin: u64,
    /// 1 / BIN: vcftools multiplies by it.
    inv: f64,
    /// n, as f64.
    n: f64,
    a1: f64,
    e1: f64,
    e2: f64,
    sorted: SortedInput,
    chrom: Vec<u8>,
    /// The current bin: its index, S and Σ p (1 − p).
    cur: Option<(u64, u64, f64)>,
    /// A bin with S > 0 has been written on this chromosome.
    started: bool,
    out: Out,
}

impl TajimaD {
    /// `bin` must be positive. The constants of Tajima (1989) are taken for
    /// n = 2 × `n_samples` chromosomes.
    pub fn new(bin: u64, n_samples: u64) -> TajimaD {
        assert!(bin > 0, "the Tajima's D bin size must be positive");
        let n = 2 * n_samples;
        let (mut a1, mut a2) = (0.0f64, 0.0f64);
        for i in 1..n {
            a1 += 1.0 / i as f64;
            a2 += 1.0 / (i * i) as f64;
        }
        let nf = n as f64;
        let b1 = (n + 1) as f64 / 3.0 / n.saturating_sub(1) as f64;
        let b2 = 2.0 * (n * n + n + 3) as f64 / 9.0 / nf / n.saturating_sub(1) as f64;
        let c1 = b1 - 1.0 / a1;
        let c2 = b2 - (n + 2) as f64 / (a1 * nf) + a2 / a1 / a1;
        TajimaD {
            bin,
            inv: 1.0 / bin as f64,
            n: nf,
            a1,
            e1: c1 / a1,
            e2: c2 / (a1 * a1 + a2),
            sorted: SortedInput::new("--tajima-d"),
            chrom: Vec::new(),
            cur: None,
            started: false,
            out: Out::closed(),
        }
    }

    /// D of a bin with `s` segregating sites whose `p (1 − p)` sum to
    /// `sum_pq`, in vcftools' operation order.
    pub fn d(&self, s: u64, sum_pq: f64) -> f64 {
        let sf = s as f64;
        let pi = 2.0 * sum_pq * self.n / (self.n - 1.0);
        let tw = sf / self.a1;
        let var = self.e1 * sf + self.e2 * sf * (sf - 1.0);
        (pi - tw) / var.sqrt()
    }

    fn row(&mut self, bin: u64, s: u64, sum_pq: f64) -> Result<(), String> {
        let d = if s > 0 { self.d(s, sum_pq) } else { f64::NAN };
        let b = &mut self.out.buf;
        b.extend_from_slice(&self.chrom);
        b.push(b'\t');
        write_u64(b, bin * self.bin);
        b.push(b'\t');
        write_u64(b, s);
        b.push(b'\t');
        write_g6(b, d);
        b.push(b'\n');
        self.out.row()
    }

    /// Write the current bin and, while rows have started, the empty bins
    /// after it up to `next` (exclusive).
    fn close_bin(&mut self, next: Option<u64>) -> Result<(), String> {
        let Some((bin, s, sum_pq)) = self.cur.take() else {
            return Ok(());
        };
        self.started |= s > 0;
        if !self.started {
            return Ok(());
        }
        self.row(bin, s, sum_pq)?;
        for empty in bin + 1..next.unwrap_or(bin + 1) {
            self.row(empty, 0, 0.0)?;
        }
        Ok(())
    }
}

impl Stat for TajimaD {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, self.ext())?;
        header(&mut self.out, "CHROM\tBIN_START\tN_SNPS\tTajimaD")
    }
    fn site(
        &mut self,
        s: &SiteView,
        _gt: Option<&dyn GenotypeSource>,
        st: Option<&SiteStats>,
    ) -> Result<(), String> {
        let st = stats_of(st)?;
        if self.sorted.check(s)? {
            self.close_bin(None)?;
            self.started = false;
            self.chrom.clear();
            self.chrom.extend_from_slice(s.chrom);
        }
        let &[c_ref, _] = st.allele_counts.as_slice() else {
            return Ok(());
        };
        if !diploid(st) {
            return Ok(());
        }
        let bin = (s.pos as f64 * self.inv) as u64;
        if self.cur.is_some_and(|c| c.0 != bin) {
            self.close_bin(Some(bin))?;
        }
        let cur = self.cur.get_or_insert((bin, 0, 0.0));
        let p = c_ref as f64 / st.called_alleles as f64;
        if p > 0.0 && p < 1.0 {
            cur.1 += 1;
            cur.2 += p * (1.0 - p);
        }
        Ok(())
    }
    fn finish(&mut self) -> io::Result<()> {
        self.close_bin(None).map_err(io::Error::other)?;
        self.out.commit()
    }
    fn needs_genotypes(&self) -> bool {
        true
    }
    fn ext(&self) -> &'static str {
        "Tajima.D"
    }
}
