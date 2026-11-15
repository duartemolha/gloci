//! Linkage disequilibrium between pairs of biallelic sites as vcftools
//! 0.1.16 computes it: the squared correlation of genotype dosages
//! (`--geno-r2`), Lewontin's D, D′ and r² of phased haplotypes
//! (`--hap-r2`), each over a window of pairs on one chromosome, and both
//! over every pair on different chromosomes (`--interchrom-geno-r2`,
//! `--interchrom-hap-r2`).
//!
//! A held site keeps its selected samples' genotypes as bit planes, one
//! 64-bit word per plane for every 64 samples: three planes for genotypes
//! (both alleles called, REF/REF, REF/ALT) and four for haplotypes (each
//! haplotype called, and ALT), so 3 or 4 bits a sample in whole words of 64
//! samples, and 16 bytes for its index and position. The sums a pair needs
//! are popcounts of their intersections; they are whole numbers, so the
//! floating-point steps that follow are vcftools' own, in its order.

use super::POLYPLOID_HINT;
use crate::filter::SiteView;
use crate::gt::{GenotypeSource, SampleFields};
use crate::stats::{Out, SiteStats, SortedInput, Stat, StatSink, no_polyploid};
use gloci_core::cfmt::{write_g6, write_u64};
use std::cell::Cell;
use std::collections::HashMap;
use std::io;
use std::rc::Rc;

/// The pairs `--geno-r2` and `--hap-r2` write: sites i < j of one
/// chromosome, counted among the sites the statistic uses, with
/// `max(1, min_snps) ≤ j − i ≤ snps` and `min_bp ≤ POSj − POSi ≤ bp`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LdWindow {
    /// `--ld-window`: at most this many sites apart.
    pub snps: u64,
    /// `--ld-window-bp`: at most this many bases apart.
    pub bp: u64,
    /// `--ld-window-min`: at least this many sites apart (0 is 1).
    pub min_snps: u64,
    /// `--ld-window-bp-min`: at least this many bases apart.
    pub min_bp: u64,
    /// `--min-r2`: when above 0, a pair whose r² is smaller, or `nan`, is
    /// left out.
    pub min_r2: Option<f64>,
}

impl Default for LdWindow {
    /// Every pair on one chromosome, as vcftools' defaults give.
    fn default() -> LdWindow {
        LdWindow {
            snps: u64::MAX,
            bp: u64::MAX,
            min_snps: 1,
            min_bp: 0,
            min_r2: None,
        }
    }
}

impl LdWindow {
    /// True when a pair `d` sites and `dpos` bases apart, inside the
    /// window's maxima, meets its minima.
    fn admits(&self, d: u64, dpos: u64) -> bool {
        d >= self.min_snps.max(1) && dpos >= self.min_bp
    }

    /// True when a site `d` sites and `dpos` bases after another is inside
    /// the window's maxima.
    fn holds(&self, d: u64, dpos: u64) -> bool {
        d <= self.snps && dpos <= self.bp
    }
}

/// vcftools' `--min-r2` test: a value above 0 drops a smaller r² and `nan`.
fn keeps(min_r2: Option<f64>, r2: f64) -> bool {
    match min_r2 {
        Some(x) if x > 0.0 => !(r2 < x || r2.is_nan()),
        _ => true,
    }
}

/// The memory the LD statistics of a run share for the sites they hold
/// (`--max-mem`). Clones share it. A statistic takes memory as its store of
/// sites grows and keeps it to the end of the run, so the bytes taken are
/// the bytes allocated.
#[derive(Clone, Debug)]
pub struct LdBudget {
    max: u64,
    used: Rc<Cell<u64>>,
}

impl LdBudget {
    /// A budget of `max` bytes.
    pub fn new(max: u64) -> LdBudget {
        LdBudget {
            max,
            used: Rc::new(Cell::new(0)),
        }
    }

    /// Take `bytes` if they fit.
    fn take(&self, bytes: u64) -> bool {
        match self.used.get().checked_add(bytes) {
            Some(t) if t <= self.max => {
                self.used.set(t);
                true
            }
            _ => false,
        }
    }

    /// The bytes not taken yet.
    fn room(&self) -> u64 {
        self.max - self.used.get()
    }

    /// The budget as `--max-mem` takes it (`2G`, `512M`, or bytes).
    fn shown(&self) -> String {
        for (shift, unit) in [(40, 'T'), (30, 'G'), (20, 'M'), (10, 'K')] {
            if self.max >= 1 << shift && self.max.is_multiple_of(1 << shift) {
                return format!("{}{unit}", self.max >> shift);
            }
        }
        format!("{} bytes", self.max)
    }
}

/// Whether a statistic scores genotypes or haplotypes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    Geno,
    Hap,
}

impl Kind {
    /// Bit planes per word of 64 samples.
    fn planes(self) -> usize {
        match self {
            Kind::Geno => 3,
            Kind::Hap => 4,
        }
    }
}

/// Reads the selected samples' genotypes of a biallelic site into bit
/// planes.
struct Encoder {
    kind: Kind,
    /// The selected samples, as their indices in the record.
    samples: Vec<usize>,
    /// Words of 64 samples.
    words: usize,
    f: SampleFields,
}

impl Encoder {
    fn new(kind: Kind, samples: Vec<usize>) -> Encoder {
        Encoder {
            kind,
            words: samples.len().div_ceil(64),
            samples,
            f: SampleFields::default(),
        }
    }

    /// The words of one site.
    fn stride(&self) -> usize {
        self.words * self.kind.planes()
    }

    /// Encode the site's genotypes into `out`, word by word, each word's
    /// planes together. Genotypes: both alleles called (a missing, masked,
    /// half-called or haploid genotype is not), REF/REF, and REF/ALT; the
    /// rest of the called ones are ALT/ALT. Haplotypes: for each of the two,
    /// called and ALT; a haploid genotype is one haplotype, as vcftools
    /// reads it, and a masked one none.
    fn encode(&mut self, gt: &dyn GenotypeSource, out: &mut [u64]) -> Result<(), String> {
        let p = self.kind.planes();
        out.fill(0);
        for (k, &i) in self.samples.iter().enumerate() {
            gt.sample(i, &mut self.f)?;
            let Some(g) = &self.f.gt else { continue };
            let (base, bit) = ((k / 64) * p, 1u64 << (k % 64));
            match self.kind {
                Kind::Geno => {
                    let &[Some(a), Some(b)] = g.alleles.as_slice() else {
                        continue;
                    };
                    out[base] |= bit;
                    if a != b {
                        out[base + 2] |= bit;
                    } else if a == 0 {
                        out[base + 1] |= bit;
                    }
                }
                Kind::Hap => {
                    let (first, second) = match *g.alleles.as_slice() {
                        [x] => (x, None),
                        [x, y] => (x, y),
                        _ => continue,
                    };
                    for (plane, allele) in [(0, first), (2, second)] {
                        if let Some(a) = allele {
                            out[base + plane] |= bit;
                            if a != 0 {
                                out[base + plane + 1] |= bit;
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

fn pc(x: u64) -> u64 {
    u64::from(x.count_ones())
}

/// vcftools' `calc_geno_r2` for two sites' planes: the samples with both
/// genotypes called and the squared correlation of their REF allele counts
/// (2, 1, 0), `cov² / (var₁ var₂)`.
fn geno_r2(a: &[u64], b: &[u64]) -> (u64, f64) {
    let (mut n, mut sx, mut sx2, mut sy, mut sy2, mut sxy) = (0u64, 0u64, 0u64, 0u64, 0u64, 0u64);
    for (wa, wb) in a.as_chunks::<3>().0.iter().zip(b.as_chunks::<3>().0) {
        let v = wa[0] & wb[0];
        if v == 0 {
            continue;
        }
        let (ha, ea, hb, eb) = (wa[1] & v, wa[2] & v, wb[1] & v, wb[2] & v);
        let (nha, nea, nhb, neb) = (pc(ha), pc(ea), pc(hb), pc(eb));
        n += pc(v);
        sx += 2 * nha + nea;
        sx2 += 4 * nha + nea;
        sy += 2 * nhb + neb;
        sy2 += 4 * nhb + neb;
        sxy += 4 * pc(ha & hb) + 2 * pc(ha & eb) + 2 * pc(ea & hb) + pc(ea & eb);
    }
    // vcftools adds the counts one sample at a time into doubles: whole
    // numbers, so the same values.
    let nf = n as f64;
    let x = sx as f64 / nf;
    let x2 = sx2 as f64 / nf;
    let y = sy as f64 / nf;
    let y2 = sy2 as f64 / nf;
    let xy = sxy as f64 / nf;
    let var1 = x2 - x * x;
    let var2 = y2 - y * y;
    let cov12 = xy - x * y;
    (n, cov12 * cov12 / (var1 * var2))
}

/// C++'s `std::min`: `b` when `b < a`, else `a` (so a `nan` `a` stays).
fn cmin(a: f64, b: f64) -> f64 {
    if b < a { b } else { a }
}

/// vcftools' `calc_hap_r2` for two sites' planes: the haplotypes called at
/// both, r², D and D′.
fn hap_r2(a: &[u64], b: &[u64]) -> (u64, f64, f64, f64) {
    // x11: REF at both; x12: REF then ALT; x21: ALT then REF; x22: ALT at both.
    let (mut x11, mut x12, mut x21, mut x22) = (0u64, 0u64, 0u64, 0u64);
    for (wa, wb) in a.as_chunks::<4>().0.iter().zip(b.as_chunks::<4>().0) {
        for h in [0, 2] {
            let v = wa[h] & wb[h];
            if v == 0 {
                continue;
            }
            let (alt_a, alt_b) = (wa[h + 1] & v, wb[h + 1] & v);
            let (n, na, nb, nab) = (pc(v), pc(alt_a), pc(alt_b), pc(alt_a & alt_b));
            x22 += nab;
            x21 += na - nab;
            x12 += nb - nab;
            x11 += n + nab - na - nb;
        }
    }
    let chr_count = x11 + x12 + x21 + x22;
    let n = chr_count as f64;
    let (x11, x12, x21, x22) = (x11 as f64, x12 as f64, x21 as f64, x22 as f64);
    let rel_x11 = x11 / n;
    let p1 = (x11 + x12) / n;
    let p2 = (x21 + x22) / n;
    let q1 = (x11 + x21) / n;
    let q2 = (x12 + x22) / n;
    let d = rel_x11 - p1 * q1;
    let dmax = if d < 0.0 {
        cmin(p1 * q1, p2 * q2)
    } else {
        cmin(p1 * q2, p2 * q1)
    };
    let dprime = d / dmax;
    // vcftools' X, X2, Y, Y2 and XY: REF is 1, ALT 0, so X = X2 = x11 + x12,
    // Y = Y2 = x11 + x21 and XY = x11, each divided by the count.
    let (x, x2) = ((x11 + x12) / n, (x11 + x12) / n);
    let (y, y2) = ((x11 + x21) / n, (x11 + x21) / n);
    let xy = x11 / n;
    let var1 = x2 - x * x;
    let var2 = y2 - y * y;
    let cov12 = xy - x * y;
    (chr_count, cov12 * cov12 / (var1 * var2), d, dprime)
}

fn header(out: &mut Out, line: &str) -> io::Result<()> {
    out.buf.extend_from_slice(line.as_bytes());
    out.buf.push(b'\n');
    out.row().map_err(io::Error::other)
}

/// The site statistics and genotypes of a site an LD statistic uses: only
/// biallelic ones, and for haplotypes only those where every selected
/// genotype is phased (vcftools' implied `--phased`, which reads the phase
/// before masking). `None` for a site it does not use. A genotype of
/// ploidy above two stops it, as vcftools stops when it parses one.
fn used<'a>(
    kind: Kind,
    st: Option<&SiteStats>,
    gt: Option<&'a dyn GenotypeSource>,
) -> Result<Option<&'a dyn GenotypeSource>, String> {
    let st = st.ok_or("internal error: this statistic needs the site statistics")?;
    if st.allele_counts.len() != 2 {
        return Ok(None);
    }
    no_polyploid(st.polyploid, POLYPLOID_HINT)?;
    if kind == Kind::Hap && !st.all_phased {
        return Ok(None);
    }
    gt.map(Some)
        .ok_or_else(|| "internal error: this statistic needs genotypes".to_string())
}

/// The sites a statistic holds, each with a key (its index or chromosome
/// id, and its position) and its planes, in a ring of slots. The slots are
/// allocated whole and taken from the budget as the ring grows (doubling,
/// but never past what the budget has left), and kept to the end, so the
/// memory the statistic holds is what the budget gave it.
struct Ring {
    /// Words of planes per site.
    stride: usize,
    keys: Vec<(u64, u64)>,
    planes: Vec<u64>,
    /// The slot of the first site held, and how many are held.
    head: usize,
    len: usize,
}

impl Ring {
    fn new(stride: usize) -> Ring {
        Ring {
            stride,
            keys: Vec::new(),
            planes: Vec::new(),
            head: 0,
            len: 0,
        }
    }

    /// The bytes of one slot.
    fn slot_bytes(&self) -> u64 {
        (self.stride * size_of::<u64>() + size_of::<(u64, u64)>()) as u64
    }

    /// The bytes allocated.
    fn bytes(&self) -> u64 {
        (self.keys.capacity() * size_of::<(u64, u64)>() + self.planes.capacity() * size_of::<u64>())
            as u64
    }

    fn slot(&self, k: usize) -> usize {
        (self.head + k) % self.keys.len()
    }

    /// The key of the `k`-th site held.
    fn key(&self, k: usize) -> (u64, u64) {
        self.keys[self.slot(k)]
    }

    /// The planes of the `k`-th site held.
    fn planes(&self, k: usize) -> &[u64] {
        let s = self.slot(k);
        &self.planes[s * self.stride..(s + 1) * self.stride]
    }

    fn pop_front(&mut self) {
        if self.len > 0 {
            self.head = (self.head + 1) % self.keys.len();
            self.len -= 1;
        }
    }

    /// The slot for one more site with `key`, whose planes the caller
    /// writes; `None` when the ring is full and `budget` has no room for
    /// another slot.
    fn push(&mut self, key: (u64, u64), budget: &LdBudget) -> Option<&mut [u64]> {
        if self.len == self.keys.len() && !self.grow(budget) {
            return None;
        }
        let s = (self.head + self.len) % self.keys.len();
        self.len += 1;
        self.keys[s] = key;
        Some(&mut self.planes[s * self.stride..(s + 1) * self.stride])
    }

    /// Double the slots, or take what the budget has left; false when it
    /// has no room for one more.
    fn grow(&mut self, budget: &LdBudget) -> bool {
        let cap = self.keys.len();
        let room = usize::try_from(budget.room() / self.slot_bytes()).unwrap_or(usize::MAX);
        let new = (cap.saturating_mul(2)).max(4).min(cap.saturating_add(room));
        if new == cap || !budget.take((new - cap) as u64 * self.slot_bytes()) {
            return false;
        }
        let mut keys = vec![(0, 0); new];
        let mut planes = vec![0; new * self.stride];
        for k in 0..self.len {
            keys[k] = self.key(k);
            planes[k * self.stride..(k + 1) * self.stride].copy_from_slice(self.planes(k));
        }
        self.keys = keys;
        self.planes = planes;
        self.head = 0;
        true
    }
}

/// `--geno-r2` and `--hap-r2`: the pairs of a window on each chromosome.
/// Sites stream through: a site is held until one past its window arrives
/// (or its chromosome ends), and then its pairs are written, so the rows
/// come in (i, j) order and only the sites of one window are held.
struct WindowLd {
    flag: &'static str,
    window: LdWindow,
    enc: Encoder,
    budget: LdBudget,
    sorted: SortedInput,
    chrom: Vec<u8>,
    /// The index of the next site used.
    next: u64,
    /// The sites held, keyed by index and position.
    ring: Ring,
    peak: usize,
    out: Out,
}

impl WindowLd {
    fn new(kind: Kind, samples: Vec<usize>, window: LdWindow, budget: LdBudget) -> WindowLd {
        let flag = match kind {
            Kind::Geno => "--geno-r2",
            Kind::Hap => "--hap-r2",
        };
        let enc = Encoder::new(kind, samples);
        WindowLd {
            flag,
            window,
            ring: Ring::new(enc.stride()),
            enc,
            budget,
            sorted: SortedInput::new(flag, "statistic"),
            chrom: Vec::new(),
            next: 0,
            peak: 0,
            out: Out::closed(),
        }
    }

    /// Write the pairs of the first held site with the sites held after it,
    /// and drop it.
    fn close_front(&mut self) -> Result<(), String> {
        if self.ring.len == 0 {
            return Ok(());
        }
        let (i_idx, i_pos) = self.ring.key(0);
        let a = self.ring.planes(0);
        for k in 1..self.ring.len {
            let (j_idx, j_pos) = self.ring.key(k);
            if !self.window.admits(j_idx - i_idx, j_pos - i_pos) {
                continue;
            }
            let b = self.ring.planes(k);
            let o = &mut self.out.buf;
            match self.enc.kind {
                Kind::Geno => {
                    let (n, r2) = geno_r2(a, b);
                    if !keeps(self.window.min_r2, r2) {
                        continue;
                    }
                    pair_prefix(o, &self.chrom, i_pos, j_pos, n);
                    write_g6(o, r2);
                }
                Kind::Hap => {
                    let (n, r2, d, dprime) = hap_r2(a, b);
                    if !keeps(self.window.min_r2, r2) {
                        continue;
                    }
                    pair_prefix(o, &self.chrom, i_pos, j_pos, n);
                    write_g6(o, r2);
                    o.push(b'\t');
                    write_g6(o, d);
                    o.push(b'\t');
                    write_g6(o, dprime);
                }
            }
            o.push(b'\n');
            self.out.row()?;
        }
        self.ring.pop_front();
        Ok(())
    }

    fn close_all(&mut self) -> Result<(), String> {
        while self.ring.len > 0 {
            self.close_front()?;
        }
        Ok(())
    }

    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, self.ext())?;
        header(
            &mut self.out,
            match self.enc.kind {
                Kind::Geno => "CHR\tPOS1\tPOS2\tN_INDV\tR^2",
                Kind::Hap => "CHR\tPOS1\tPOS2\tN_CHR\tR^2\tD\tDprime",
            },
        )
    }

    fn site(
        &mut self,
        s: &SiteView,
        gt: Option<&dyn GenotypeSource>,
        st: Option<&SiteStats>,
    ) -> Result<(), String> {
        if self.sorted.check(s)? {
            self.close_all()?;
            self.chrom.clear();
            self.chrom.extend_from_slice(s.chrom);
        }
        let Some(gt) = used(self.enc.kind, st, gt)? else {
            return Ok(());
        };
        let idx = self.next;
        self.next += 1;
        while self.ring.len > 0 {
            let (f_idx, f_pos) = self.ring.key(0);
            if self.window.holds(idx - f_idx, s.pos - f_pos) {
                break;
            }
            self.close_front()?;
        }
        let Some(slot) = self.ring.push((idx, s.pos), &self.budget) else {
            return Err(format!(
                "{} holds the sites of a pair window in memory, and {} sites of {} samples pass \
                 --max-mem {}; narrow the window with --ld-window-bp or --ld-window, or raise \
                 --max-mem",
                self.flag,
                self.ring.len + 1,
                self.enc.samples.len(),
                self.budget.shown()
            ));
        };
        self.enc.encode(gt, slot)?;
        self.peak = self.peak.max(self.ring.len);
        Ok(())
    }

    fn finish(&mut self) -> io::Result<()> {
        self.close_all().map_err(io::Error::other)?;
        self.out.commit()
    }

    fn ext(&self) -> &'static str {
        match self.enc.kind {
            Kind::Geno => "geno.ld",
            Kind::Hap => "hap.ld",
        }
    }
}

/// `CHROM POS1 POS2 N ` of a windowed row.
fn pair_prefix(b: &mut Vec<u8>, chrom: &[u8], pos1: u64, pos2: u64, n: u64) {
    b.extend_from_slice(chrom);
    b.push(b'\t');
    write_u64(b, pos1);
    b.push(b'\t');
    write_u64(b, pos2);
    b.push(b'\t');
    write_u64(b, n);
    b.push(b'\t');
}

macro_rules! window_stat {
    ($name:ident, $kind:expr) => {
        impl $name {
            /// Over the `samples` (their indices in the record), the pairs
            /// of `window`, holding sites within `budget`.
            pub fn new(samples: Vec<usize>, window: LdWindow, budget: LdBudget) -> $name {
                $name(WindowLd::new($kind, samples, window, budget))
            }

            /// The most sites held at once so far. For tests.
            #[doc(hidden)]
            pub fn peak_held(&self) -> usize {
                self.0.peak
            }

            /// The bytes allocated for held sites now. For tests.
            #[doc(hidden)]
            pub fn held_bytes(&self) -> u64 {
                self.0.ring.bytes()
            }
        }

        impl Stat for $name {
            fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
                self.0.start(sink)
            }
            fn site(
                &mut self,
                s: &SiteView,
                gt: Option<&dyn GenotypeSource>,
                st: Option<&SiteStats>,
            ) -> Result<(), String> {
                self.0.site(s, gt, st)
            }
            fn finish(&mut self) -> io::Result<()> {
                self.0.finish()
            }
            fn needs_genotypes(&self) -> bool {
                true
            }
            fn reads_genotypes(&self) -> bool {
                true
            }
            fn ext(&self) -> &'static str {
                self.0.ext()
            }
            fn set_input(&mut self, name: &str) {
                self.0.sorted.set_input(name);
            }
        }
    };
}

/// `--geno-r2` (`.geno.ld`): for each pair of biallelic sites in the
/// window, the samples whose genotypes are called at both (`N_INDV`) and
/// the squared correlation of their allele counts. A missing, half-called
/// (`./1`), haploid or masked genotype is left out of the pair. Input must
/// be sorted by chromosome and position.
pub struct GenoR2(WindowLd);

/// `--hap-r2` (`.hap.ld`): for each pair of biallelic sites in the window
/// where every selected genotype is phased, the haplotypes called at both
/// (`N_CHR`), r², Lewontin's D and D′. A haploid genotype is one haplotype,
/// a half-called one (`.|1`) its called haplotype, and a masked one none.
/// Sites with an unphased genotype (`./.` included) are not used and take
/// no place in `--ld-window`. Input must be sorted by chromosome and
/// position.
pub struct HapR2(WindowLd);

window_stat!(GenoR2, Kind::Geno);
window_stat!(HapR2, Kind::Hap);

/// `--interchrom-geno-r2` and `--interchrom-hap-r2`
/// (`.interchrom.geno.ld`, `.interchrom.hap.ld`): every pair of sites on
/// different chromosomes, with `--geno-r2`'s or `--hap-r2`'s values (no D
/// or D′), in input order (i, j), written at the end. Every site used is
/// held; the input may be in any order, and sites of a chromosome name
/// never pair with each other.
pub struct InterchromR2 {
    flag: &'static str,
    enc: Encoder,
    min_r2: Option<f64>,
    budget: LdBudget,
    /// The chromosome names, and the id of each.
    names: Vec<Vec<u8>>,
    ids: HashMap<Vec<u8>, u32>,
    /// Every site used, keyed by chromosome id and position, in input order.
    ring: Ring,
    out: Out,
}

impl InterchromR2 {
    /// Over the `samples` (their indices in the record): haplotypes when
    /// `haplotypes`, else genotypes; `min_r2` as in [`LdWindow`]; every
    /// site held within `budget`.
    pub fn new(
        samples: Vec<usize>,
        haplotypes: bool,
        min_r2: Option<f64>,
        budget: LdBudget,
    ) -> InterchromR2 {
        let (kind, flag) = if haplotypes {
            (Kind::Hap, "--interchrom-hap-r2")
        } else {
            (Kind::Geno, "--interchrom-geno-r2")
        };
        let enc = Encoder::new(kind, samples);
        InterchromR2 {
            flag,
            ring: Ring::new(enc.stride()),
            enc,
            min_r2,
            budget,
            names: Vec::new(),
            ids: HashMap::new(),
            out: Out::closed(),
        }
    }

    /// The bytes allocated for held sites now. For tests.
    #[doc(hidden)]
    pub fn held_bytes(&self) -> u64 {
        self.ring.bytes()
    }

    /// The id of chromosome `name`, given one if it has none.
    fn chrom_id(&mut self, name: &[u8]) -> u32 {
        // Most sites are on the chromosome of the site before.
        if let Some(last) = self.ring.len.checked_sub(1) {
            let id = self.ring.key(last).0 as u32;
            if self.names[id as usize] == name {
                return id;
            }
        }
        if let Some(&id) = self.ids.get(name) {
            return id;
        }
        let id = self.names.len() as u32;
        self.names.push(name.to_vec());
        self.ids.insert(name.to_vec(), id);
        id
    }

    fn write_pairs(&mut self) -> Result<(), String> {
        for i in 0..self.ring.len {
            let (ci, pi) = self.ring.key(i);
            let a = self.ring.planes(i);
            for j in i + 1..self.ring.len {
                let (cj, pj) = self.ring.key(j);
                if ci == cj {
                    continue;
                }
                let b = self.ring.planes(j);
                let (n, r2) = match self.enc.kind {
                    Kind::Geno => geno_r2(a, b),
                    Kind::Hap => {
                        let (n, r2, _, _) = hap_r2(a, b);
                        (n, r2)
                    }
                };
                if !keeps(self.min_r2, r2) {
                    continue;
                }
                let o = &mut self.out.buf;
                o.extend_from_slice(&self.names[ci as usize]);
                o.push(b'\t');
                write_u64(o, pi);
                o.push(b'\t');
                o.extend_from_slice(&self.names[cj as usize]);
                o.push(b'\t');
                write_u64(o, pj);
                o.push(b'\t');
                write_u64(o, n);
                o.push(b'\t');
                write_g6(o, r2);
                o.push(b'\n');
                self.out.row()?;
            }
        }
        Ok(())
    }
}

impl Stat for InterchromR2 {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, self.ext())?;
        header(
            &mut self.out,
            match self.enc.kind {
                Kind::Geno => "CHR1\tPOS1\tCHR2\tPOS2\tN_INDV\tR^2",
                Kind::Hap => "CHR1\tPOS1\tCHR2\tPOS2\tN_CHR\tR^2",
            },
        )
    }
    fn site(
        &mut self,
        s: &SiteView,
        gt: Option<&dyn GenotypeSource>,
        st: Option<&SiteStats>,
    ) -> Result<(), String> {
        let Some(gt) = used(self.enc.kind, st, gt)? else {
            return Ok(());
        };
        let id = self.chrom_id(s.chrom);
        let Some(slot) = self.ring.push((u64::from(id), s.pos), &self.budget) else {
            return Err(format!(
                "{} holds every site it uses until the end of the input, and {} sites of {} \
                 samples pass --max-mem {}; raise --max-mem, or use fewer sites or samples \
                 (--region, --samples)",
                self.flag,
                self.ring.len + 1,
                self.enc.samples.len(),
                self.budget.shown()
            ));
        };
        self.enc.encode(gt, slot)?;
        Ok(())
    }
    fn finish(&mut self) -> io::Result<()> {
        self.write_pairs().map_err(io::Error::other)?;
        self.out.commit()
    }
    fn needs_genotypes(&self) -> bool {
        true
    }
    fn reads_genotypes(&self) -> bool {
        true
    }
    fn ext(&self) -> &'static str {
        match self.enc.kind {
            Kind::Geno => "interchrom.geno.ld",
            Kind::Hap => "interchrom.hap.ld",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Planes of one word from per-sample genotypes (`None` missing).
    fn geno_planes(g: &[Option<(u32, u32)>]) -> Vec<u64> {
        let mut p = vec![0u64; 3];
        for (k, x) in g.iter().enumerate() {
            let Some((a, b)) = *x else { continue };
            p[0] |= 1 << k;
            if a != b {
                p[2] |= 1 << k;
            } else if a == 0 {
                p[1] |= 1 << k;
            }
        }
        p
    }

    #[test]
    fn geno_sums_match_a_sample_loop() {
        // vcftools' loop over the samples called at both, by hand.
        let a = [Some((0, 0)), Some((0, 1)), Some((1, 1)), None, Some((1, 0))];
        let b = [Some((0, 1)), Some((1, 1)), Some((1, 1)), Some((0, 0)), None];
        let (n, r2) = geno_r2(&geno_planes(&a), &geno_planes(&b));
        // REF counts at both: (2, 1), (1, 0), (0, 0).
        let (x, y) = ([2.0, 1.0, 0.0], [1.0, 0.0, 0.0]);
        let m = |v: &[f64]| v.iter().sum::<f64>() / 3.0;
        let (mx, my) = (m(&x), m(&y));
        let xy = x.iter().zip(&y).map(|(p, q)| p * q).sum::<f64>() / 3.0;
        let x2 = x.iter().map(|p| p * p).sum::<f64>() / 3.0;
        let y2 = y.iter().map(|q| q * q).sum::<f64>() / 3.0;
        let want = (xy - mx * my).powi(2) / ((x2 - mx * mx) * (y2 - my * my));
        assert_eq!(n, 3);
        assert_eq!(r2.to_bits(), want.to_bits());
    }

    #[test]
    fn min_r2_follows_vcftools() {
        assert!(keeps(None, f64::NAN));
        assert!(keeps(Some(0.0), f64::NAN));
        assert!(keeps(Some(-1.0), 0.0));
        assert!(!keeps(Some(0.2), f64::NAN));
        assert!(!keeps(Some(0.2), 0.19));
        assert!(keeps(Some(0.2), 0.2));
    }

    #[test]
    fn budget_shows_as_max_mem_is_written() {
        assert_eq!(LdBudget::new(2 << 30).shown(), "2G");
        assert_eq!(LdBudget::new(1024).shown(), "1K");
        assert_eq!(LdBudget::new(1536 << 20).shown(), "1536M");
        assert_eq!(LdBudget::new(1000).shown(), "1000 bytes");
        let b = LdBudget::new(100);
        let c = b.clone();
        assert!(b.take(60) && !c.take(41) && c.take(40));
        assert_eq!(b.room(), 0);
    }
}
