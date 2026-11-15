//! Weir and Cockerham's Fst between populations (`--weir-fst-pop`, with
//! `--fst-window-size` and `--fst-window-step`) as vcftools 0.1.16 computes
//! it: the multi-allele estimator of Weir and Cockerham (1984) per site, in
//! windows, and its mean and weighted mean over the run.

use super::window::OpenWindows;
use crate::filter::SiteView;
use crate::gt::{GenotypeSource, SampleFields};
use crate::stats::{Out, SiteStats, SortedInput, Stat, StatSink};
use gloci_core::cfmt::{write_g6, write_u64};
use gloci_io::{IoError, UsageError};
use std::collections::{HashMap, HashSet};
use std::io;
use std::path::PathBuf;

fn usage(msg: String) -> IoError {
    UsageError(msg).into()
}

/// The populations of `--fst-pop`, among the selected samples.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Populations {
    /// The population files as given, for messages.
    pub names: Vec<String>,
    /// Each population's samples, as their indices in the record, ascending.
    pub members: Vec<Vec<usize>>,
}

/// C's `isspace`, which splits vcftools' `>>` words.
fn space(b: &u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\n' | b'\x0b' | b'\x0c' | b'\r')
}

impl Populations {
    /// Read one population per file: the first word of each line names a
    /// sample (blank lines are skipped). `selected` are the selected samples
    /// as (index in the record, name). A name that is not a selected sample
    /// is left out, and all of them are reported in one warning with their
    /// count. Fewer than two files, a file that names no selected sample, or
    /// a sample in two files is a usage error naming it.
    pub fn load(files: &[PathBuf], selected: &[(usize, String)]) -> Result<Populations, IoError> {
        if files.len() < 2 {
            return Err(usage(TWO_FILES.to_string()));
        }
        let index: HashMap<&[u8], usize> = selected
            .iter()
            .map(|(i, name)| (name.as_bytes(), *i))
            .collect();
        let names: Vec<String> = files.iter().map(|f| f.display().to_string()).collect();
        // The population of each sample placed so far.
        let mut owner: HashMap<usize, usize> = HashMap::new();
        let mut unknown: Vec<Vec<u8>> = Vec::new();
        let mut seen_unknown: HashSet<Vec<u8>> = HashSet::new();
        let mut members = Vec::with_capacity(files.len());
        for (p, file) in files.iter().enumerate() {
            let text = std::fs::read(file).map_err(|source| IoError::Open {
                path: names[p].clone(),
                source,
            })?;
            let mut m = Vec::new();
            for line in text.split(|&b| b == b'\n') {
                let Some(word) = line.split(space).find(|w| !w.is_empty()) else {
                    continue;
                };
                let Some(&i) = index.get(word) else {
                    if seen_unknown.insert(word.to_vec()) {
                        unknown.push(word.to_vec());
                    }
                    continue;
                };
                match owner.get(&i) {
                    Some(&q) if q != p => {
                        return Err(usage(format!(
                            "sample {} is in two population files, {} and {}; a sample belongs \
                             to one population, so remove it from one of them",
                            String::from_utf8_lossy(word),
                            names[q],
                            names[p]
                        )));
                    }
                    Some(_) => {}
                    None => {
                        owner.insert(i, p);
                        m.push(i);
                    }
                }
            }
            if m.is_empty() {
                return Err(usage(format!(
                    "population file {} names no selected sample, and Fst needs at least one in \
                     each population; check its names against the input's sample names, \
                     --samples and --exclude-samples",
                    names[p]
                )));
            }
            m.sort_unstable();
            members.push(m);
        }
        if let Some(first) = unknown.first() {
            let n = unknown.len();
            log::warn!(
                "vcf popgen: {n} name{} in the --fst-pop files {} not a selected sample and {} \
                 left out of Fst (the first is {})",
                if n == 1 { "" } else { "s" },
                if n == 1 { "is" } else { "are" },
                if n == 1 { "is" } else { "are" },
                String::from_utf8_lossy(first)
            );
        }
        Ok(Populations { names, members })
    }
}

/// The message for fewer than two population files.
pub const TWO_FILES: &str = "--fst-pop needs at least two population files, one per \
                             population; give --fst-pop once for each population";

/// Reusable buffers of [`terms`].
#[derive(Default)]
struct Scratch {
    n: Vec<f64>,
    p: Vec<f64>,
    pbar: Vec<f64>,
    hbar: Vec<f64>,
    ssqr: Vec<f64>,
}

fn reset(v: &mut Vec<f64>, len: usize) {
    v.clear();
    v.resize(len, 0.0);
}

/// Σ a and Σ (a + b + c) over the alleles for which none of a, b and c is
/// NaN, for `n_pops` populations with `hom[i·k + j]` homozygotes of allele
/// j and `het[i·k + j]` heterozygotes carrying it (k alleles), in vcftools'
/// order of operations.
fn terms(hom: &[u64], het: &[u64], n_pops: usize, k: usize, s: &mut Scratch) -> (f64, f64) {
    let r = n_pops as f64;
    reset(&mut s.n, n_pops);
    reset(&mut s.p, n_pops * k);
    reset(&mut s.pbar, k);
    reset(&mut s.hbar, k);
    reset(&mut s.ssqr, k);
    let mut sum_nsqr = 0.0;
    for i in 0..n_pops {
        for j in 0..k {
            let (ho, he) = (hom[i * k + j], het[i * k + j]);
            s.n[i] += ho as f64 + 0.5 * he as f64;
            s.p[i * k + j] = (he + 2 * ho) as f64;
            s.pbar[j] += s.p[i * k + j];
            s.hbar[j] += he as f64;
        }
        for j in 0..k {
            s.p[i * k + j] /= 2.0 * s.n[i];
        }
        sum_nsqr += s.n[i] * s.n[i];
    }
    // vcftools adds the n into an int: each is a whole number of genotypes,
    // so the sum is the same.
    let n_sum: f64 = s.n.iter().sum();
    let nbar = n_sum / r;
    for j in 0..k {
        s.pbar[j] /= n_sum * 2.0;
        s.hbar[j] /= n_sum;
    }
    for j in 0..k {
        for i in 0..n_pops {
            let d = s.p[i * k + j] - s.pbar[j];
            s.ssqr[j] += s.n[i] * d * d;
        }
        s.ssqr[j] /= (r - 1.0) * nbar;
    }
    let nc = (n_sum - (sum_nsqr / n_sum)) / (r - 1.0);
    let (mut sum_a, mut sum_all) = (0.0, 0.0);
    for j in 0..k {
        let (pbar, hbar, ssqr) = (s.pbar[j], s.hbar[j], s.ssqr[j]);
        let a = (ssqr
            - (pbar * (1.0 - pbar) - (((r - 1.0) * ssqr) / r) - (hbar / 4.0)) / (nbar - 1.0))
            * nbar
            / nc;
        let b = (pbar * (1.0 - pbar)
            - (ssqr * (r - 1.0) / r)
            - hbar * (((2.0 * nbar) - 1.0) / (4.0 * nbar)))
            * nbar
            / (nbar - 1.0);
        let c = hbar / 2.0;
        if !a.is_nan() && !b.is_nan() && !c.is_nan() {
            sum_a += a;
            sum_all += a + b + c;
        }
    }
    (sum_a, sum_all)
}

/// The sums of one window of `--fst-window-size`.
#[derive(Clone, Copy, Debug, Default)]
struct Window {
    sum_a: f64,
    sum_all: f64,
    sum_fst: f64,
    /// The sites with an estimate.
    sites: u64,
}

/// `.windowed.weir.fst`: the windows of D159 on each chromosome.
struct Windows {
    windows: OpenWindows<Window>,
    sorted: SortedInput,
    chrom: Vec<u8>,
    out: Out,
}

impl Windows {
    /// Write the open windows below index `end` (all of them for
    /// `u64::MAX`), each that has a site and a Σ(a + b + c) that is a
    /// number other than 0, as vcftools does.
    fn close_before(&mut self, end: u64) -> Result<(), String> {
        while let Some((k, w)) = self.windows.pop_before(end) {
            if w.sum_all == 0.0 || w.sum_a.is_nan() || w.sum_all.is_nan() || w.sites == 0 {
                continue;
            }
            let (start, last) = self.windows.bounds(k);
            let b = &mut self.out.buf;
            b.extend_from_slice(&self.chrom);
            b.push(b'\t');
            write_u64(b, start);
            b.push(b'\t');
            write_u64(b, last);
            b.push(b'\t');
            // vcftools counts in a double and prints it as one.
            write_g6(b, w.sites as f64);
            b.push(b'\t');
            write_g6(b, w.sum_a / w.sum_all);
            b.push(b'\t');
            write_g6(b, w.sum_fst / w.sites as f64);
            b.push(b'\n');
            self.out.row()?;
        }
        Ok(())
    }
}

/// `--fst-pop FILE` (two or more): Weir and Cockerham's Fst between the
/// populations, with vcftools' multi-allele form. Per site where every
/// population member that is not masked is diploid, each population's
/// fully called genotypes give its size n_i, allele frequencies and
/// heterozygote frequencies; per allele these give vcftools' a, b and c,
/// and the site's Fst is Σa / Σ(a + b + c) over the alleles where none is
/// NaN. Every such site has a row in `.weir.fst` (`nan` when a population
/// has no call, or nothing varies). `.weir.fst.summary` gives the sites
/// with an estimate, the mean of their Fst and the weighted Fst
/// Σa / Σ(a + b + c) over them, which vcftools only logs.
///
/// With a window (`--fst-window-size`, `--fst-window-step`),
/// `.windowed.weir.fst` sums a, a + b + c and Fst over the sites with an
/// estimate in each window `[k·STEP + 1, k·STEP + SIZE]`. Windows stream
/// as positions pass them, so the input must be sorted by chromosome and
/// position, and at most `⌈SIZE/STEP⌉ + 1` are open.
pub struct WeirFst {
    pops: Populations,
    /// Per population and allele (`i·k + j`): homozygotes and heterozygotes.
    hom: Vec<u64>,
    het: Vec<u64>,
    scratch: Scratch,
    f: SampleFields,
    /// Over the sites with an estimate: Σa, Σ(a + b + c), ΣFst and their
    /// count.
    total: (f64, f64, f64, u64),
    windows: Option<Windows>,
    out: Out,
    summary: Out,
}

impl WeirFst {
    /// `window` is (SIZE, STEP), both positive and STEP at most SIZE.
    pub fn new(pops: Populations, window: Option<(u64, u64)>) -> WeirFst {
        WeirFst {
            pops,
            hom: Vec::new(),
            het: Vec::new(),
            scratch: Scratch::default(),
            f: SampleFields::default(),
            total: (0.0, 0.0, 0.0, 0),
            windows: window.map(|(size, step)| Windows {
                windows: OpenWindows::new(size, step),
                sorted: SortedInput::new("--fst-window-size", "statistic"),
                chrom: Vec::new(),
                out: Out::closed(),
            }),
            out: Out::closed(),
            summary: Out::closed(),
        }
    }

    /// Σa and Σ(a + b + c) of a site whose population i has `hom[i][j]`
    /// homozygotes of allele j and `het[i][j]` heterozygotes carrying it, in
    /// vcftools' order of operations; Fst is their ratio.
    pub fn terms(hom: &[Vec<u64>], het: &[Vec<u64>]) -> (f64, f64) {
        let k = hom.first().map_or(0, Vec::len);
        let flat = |v: &[Vec<u64>]| v.iter().flatten().copied().collect::<Vec<u64>>();
        terms(
            &flat(hom),
            &flat(het),
            hom.len(),
            k,
            &mut Scratch::default(),
        )
    }

    /// The most windows held open at once so far: at most
    /// `⌈SIZE/STEP⌉ + 1`. For tests.
    #[doc(hidden)]
    pub fn peak_open(&self) -> usize {
        self.windows.as_ref().map_or(0, |w| w.windows.peak())
    }

    /// Count the members' fully called genotypes into `hom` and `het`, for
    /// `k` alleles. False when a member that is not masked is not diploid
    /// (a masked genotype reads as missing).
    fn count(&mut self, gt: &dyn GenotypeSource, k: usize) -> Result<bool, String> {
        let n = self.pops.members.len() * k;
        self.hom.clear();
        self.hom.resize(n, 0);
        self.het.clear();
        self.het.resize(n, 0);
        for (p, members) in self.pops.members.iter().enumerate() {
            for &i in members {
                gt.sample(i, &mut self.f)?;
                // No GT (or a masked one) is a missing diploid genotype.
                let Some(g) = &self.f.gt else { continue };
                let &[x, y] = g.alleles.as_slice() else {
                    return Ok(false);
                };
                let (Some(x), Some(y)) = (x, y) else { continue };
                let (x, y) = (x as usize, y as usize);
                if x.max(y) >= k {
                    return Err(format!(
                        "genotype allele {} but the site lists {k} alleles",
                        x.max(y)
                    ));
                }
                if x == y {
                    self.hom[p * k + x] += 1;
                } else {
                    self.het[p * k + x] += 1;
                    self.het[p * k + y] += 1;
                }
            }
        }
        Ok(true)
    }
}

fn header(out: &mut Out, line: &str) -> io::Result<()> {
    out.buf.extend_from_slice(line.as_bytes());
    out.buf.push(b'\n');
    out.row().map_err(io::Error::other)
}

impl Stat for WeirFst {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, "weir.fst")?;
        header(&mut self.out, "CHROM\tPOS\tWEIR_AND_COCKERHAM_FST")?;
        self.summary.open(sink, "weir.fst.summary")?;
        if let Some(w) = &mut self.windows {
            w.out.open(sink, "windowed.weir.fst")?;
            header(
                &mut w.out,
                "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST",
            )?;
        }
        Ok(())
    }
    fn site(
        &mut self,
        s: &SiteView,
        gt: Option<&dyn GenotypeSource>,
        st: Option<&SiteStats>,
    ) -> Result<(), String> {
        let st =
            st.ok_or("internal error: this statistic needs the site statistics".to_string())?;
        let gt = gt.ok_or("internal error: this statistic needs genotypes".to_string())?;
        if let Some(w) = &mut self.windows {
            if w.sorted.check(s)? {
                w.close_before(u64::MAX)?;
                w.chrom.clear();
                w.chrom.extend_from_slice(s.chrom);
            }
            let first = w.windows.first_holding(s.pos);
            w.close_before(first)?;
        }
        let k = st.allele_counts.len();
        if !self.count(gt, k)? {
            return Ok(());
        }
        let n_pops = self.pops.members.len();
        let (sum_a, sum_all) = terms(&self.hom, &self.het, n_pops, k, &mut self.scratch);
        let fst = sum_a / sum_all;
        if !fst.is_nan() {
            self.total.0 += sum_a;
            self.total.1 += sum_all;
            self.total.2 += fst;
            self.total.3 += 1;
            if let Some(w) = &mut self.windows {
                w.windows.add(s.pos, |x| {
                    x.sum_a += sum_a;
                    x.sum_all += sum_all;
                    x.sum_fst += fst;
                    x.sites += 1;
                });
            }
        }
        let b = &mut self.out.buf;
        b.extend_from_slice(s.chrom);
        b.push(b'\t');
        write_u64(b, s.pos);
        b.push(b'\t');
        write_g6(b, fst);
        b.push(b'\n');
        self.out.row()
    }
    fn finish(&mut self) -> io::Result<()> {
        self.out.commit()?;
        let (sum_a, sum_all, sum_fst, n) = self.total;
        let b = &mut self.summary.buf;
        b.extend_from_slice(b"N_SITES\tMEAN_FST\tWEIGHTED_FST\n");
        write_u64(b, n);
        b.push(b'\t');
        write_g6(b, sum_fst / n as f64);
        b.push(b'\t');
        write_g6(b, sum_a / sum_all);
        b.push(b'\n');
        self.summary.row().map_err(io::Error::other)?;
        self.summary.commit()?;
        if let Some(w) = &mut self.windows {
            w.close_before(u64::MAX).map_err(io::Error::other)?;
            w.out.commit()?;
        }
        Ok(())
    }
    fn needs_genotypes(&self) -> bool {
        true
    }
    fn reads_genotypes(&self) -> bool {
        true
    }
    fn ext(&self) -> &'static str {
        "weir.fst"
    }
}
