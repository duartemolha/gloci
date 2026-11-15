//! `gloci vcf popgen`: vcftools' population genetics statistics in one pass
//! over the records that pass the shared selection stage, through the
//! `vcf stats` engine (decisions D158 to D161). Each statistic writes
//! `PREFIX.<extension>` with vcftools' name and columns; the files are
//! published together only when the whole run succeeds.

use crate::vcf_select::{SelectedInput, VcfSelect};
use crate::vcf_stats::{StatsRun, StatsSummary, run_stats, selection_fields};
use gloci_io::{IoError, UsageError};
use gloci_vcf::Stat;
use gloci_vcf::popgen::{
    FST_TWO_FILES, GenoR2, HapR2, InterchromR2, LdBudget, LdWindow, Populations, SitePi, TajimaD,
    WeirFst, WindowPi,
};
use std::num::NonZeroUsize;
use std::path::PathBuf;

const COMMAND: &str = "vcf popgen";

/// Which statistics to compute.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PopgenRequest {
    /// `PREFIX.sites.pi` (vcftools `--site-pi`).
    pub site_pi: bool,
    /// `PREFIX.windowed.pi` (vcftools `--window-pi`): window size and step.
    pub window_pi: Option<(u64, u64)>,
    /// `PREFIX.Tajima.D` (vcftools `--TajimaD`): the bin size.
    pub tajima_d: Option<u64>,
    /// `PREFIX.weir.fst` and `PREFIX.weir.fst.summary`, with a window also
    /// `PREFIX.windowed.weir.fst` (vcftools `--weir-fst-pop`).
    pub fst: Option<FstRequest>,
    /// Linkage disequilibrium.
    pub ld: LdRequest,
}

/// Weir and Cockerham's Fst between populations.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FstRequest {
    /// One file per population (`--fst-pop`), at least two.
    pub pop_files: Vec<PathBuf>,
    /// Window size and step (`--fst-window-size`, `--fst-window-step`).
    pub window: Option<(u64, u64)>,
}

/// Linkage disequilibrium between pairs of biallelic sites.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LdRequest {
    /// `PREFIX.geno.ld` (vcftools `--geno-r2`).
    pub geno: bool,
    /// `PREFIX.hap.ld` (vcftools `--hap-r2`).
    pub hap: bool,
    /// `PREFIX.interchrom.geno.ld` (vcftools `--interchrom-geno-r2`).
    pub interchrom_geno: bool,
    /// `PREFIX.interchrom.hap.ld` (vcftools `--interchrom-hap-r2`).
    pub interchrom_hap: bool,
    /// The pairs of `geno` and `hap`; its `min_r2` applies to all four.
    pub window: LdWindow,
}

impl LdRequest {
    /// True when any LD statistic is requested.
    pub fn any(&self) -> bool {
        self.geno || self.hap || self.interchrom_geno || self.interchrom_hap
    }
}

impl PopgenRequest {
    /// True when no statistic is requested.
    pub fn is_empty(&self) -> bool {
        !self.site_pi
            && self.window_pi.is_none()
            && self.tajima_d.is_none()
            && self.fst.is_none()
            && !self.ld.any()
    }
}

/// Options for `vcf popgen`.
#[derive(Clone, Debug)]
pub struct PopgenOptions {
    pub select: VcfSelect,
    pub request: PopgenRequest,
    /// Output files are `PREFIX.<extension>`.
    pub out_prefix: PathBuf,
    /// The bytes the LD statistics may hold sites in (`--max-mem`).
    pub max_mem: u64,
}

fn usage(msg: impl Into<String>) -> IoError {
    UsageError(msg.into()).into()
}

/// Check the request before any input is read.
pub fn check_request(o: &PopgenOptions) -> Result<(), IoError> {
    o.select.validate()?;
    let r = &o.request;
    if r.is_empty() {
        return Err(usage(
            "no statistic requested; add at least one of --site-pi, --window-pi SIZE, \
             --tajima-d BIN, --fst-pop FILE (twice or more), --geno-r2, --hap-r2, \
             --interchrom-geno-r2 or --interchrom-hap-r2",
        ));
    }
    if let Some(w) = r.window_pi {
        check_window(w, "--window-pi", "--window-pi-step")?;
    }
    if r.tajima_d == Some(0) {
        return Err(usage(
            "--tajima-d needs a bin size above 0; give the bin width in bases, such as \
             --tajima-d 10000",
        ));
    }
    if let Some(f) = &r.fst {
        if f.pop_files.len() < 2 {
            return Err(usage(FST_TWO_FILES));
        }
        if let Some(w) = f.window {
            check_window(w, "--fst-window-size", "--fst-window-step")?;
        }
    }
    check_ld(&r.ld)?;
    if o.out_prefix.as_os_str().is_empty() {
        return Err(usage(
            "--out-prefix is empty; give a prefix such as --out-prefix results/cohort",
        ));
    }
    Ok(())
}

/// The LD window: options only where a statistic uses them, and a window
/// that admits a pair.
fn check_ld(ld: &LdRequest) -> Result<(), IoError> {
    let w = &ld.window;
    let d = LdWindow::default();
    let windowed =
        w.snps != d.snps || w.bp != d.bp || w.min_snps != d.min_snps || w.min_bp != d.min_bp;
    if windowed && !ld.geno && !ld.hap {
        return Err(usage(
            "--ld-window, --ld-window-bp, --ld-window-min and --ld-window-bp-min choose the pairs \
             of --geno-r2 and --hap-r2; add one of them",
        ));
    }
    if w.min_r2.is_some() && !ld.any() {
        return Err(usage(
            "--min-r2 filters the pairs of --geno-r2, --hap-r2, --interchrom-geno-r2 and \
             --interchrom-hap-r2; add one of them",
        ));
    }
    if w.snps == 0 {
        return Err(usage(
            "--ld-window 0 admits no pair, since the two sites of a pair are at least 1 site \
             apart; give 1 or more",
        ));
    }
    if w.min_snps > w.snps {
        return Err(usage(format!(
            "--ld-window-min {} is above --ld-window {}, which admits no pair; give a minimum of \
             at most {}",
            w.min_snps, w.snps, w.snps
        )));
    }
    if w.min_bp > w.bp {
        return Err(usage(format!(
            "--ld-window-bp-min {} is above --ld-window-bp {}, which admits no pair; give a \
             minimum of at most {}",
            w.min_bp, w.bp, w.bp
        )));
    }
    Ok(())
}

/// A window's size and step, for the flags `size_flag` and `step_flag`:
/// both above 0, and the step at most the size (D158).
fn check_window((size, step): (u64, u64), size_flag: &str, step_flag: &str) -> Result<(), IoError> {
    if size == 0 {
        return Err(usage(format!(
            "{size_flag} needs a window size above 0; give the window width in bases, such as \
             {size_flag} 10000"
        )));
    }
    if step == 0 {
        return Err(usage(format!(
            "{step_flag} needs a step above 0; give it in bases, or leave it out for windows \
             that do not overlap"
        )));
    }
    if step > size {
        return Err(usage(format!(
            "{step_flag} {step} is larger than the window ({size}), which would leave positions \
             in no window; give a step of at most {size}"
        )));
    }
    Ok(())
}

/// The statistics a request asks for, over the selected samples
/// (`samples`: index in the record and name); the LD statistics hold sites
/// within `max_mem` bytes together.
fn build_stats(
    r: &PopgenRequest,
    samples: &[(usize, String)],
    max_mem: u64,
) -> Result<Vec<Box<dyn Stat>>, IoError> {
    let mut v: Vec<Box<dyn Stat>> = Vec::new();
    if r.site_pi {
        v.push(Box::new(SitePi::new()));
    }
    if let Some((size, step)) = r.window_pi {
        v.push(Box::new(WindowPi::new(size, step)));
    }
    if let Some(bin) = r.tajima_d {
        v.push(Box::new(TajimaD::new(bin, samples.len() as u64)));
    }
    if let Some(f) = &r.fst {
        let pops = Populations::load(&f.pop_files, samples)?;
        v.push(Box::new(WeirFst::new(pops, f.window)));
    }
    let ld = &r.ld;
    let indices = || samples.iter().map(|(i, _)| *i).collect::<Vec<usize>>();
    let budget = LdBudget::new(max_mem);
    if ld.geno {
        v.push(Box::new(GenoR2::new(indices(), ld.window, budget.clone())));
    }
    if ld.hap {
        v.push(Box::new(HapR2::new(indices(), ld.window, budget.clone())));
    }
    for (on, haplotypes) in [(ld.interchrom_geno, false), (ld.interchrom_hap, true)] {
        if on {
            v.push(Box::new(InterchromR2::new(
                indices(),
                haplotypes,
                ld.window.min_r2,
                budget.clone(),
            )));
        }
    }
    Ok(v)
}

/// Compute the requested statistics of the records of `input`, with the
/// selection on `threads` threads, and write `PREFIX.<extension>` for each.
pub fn vcf_popgen(
    input: SelectedInput,
    o: &PopgenOptions,
    threads: NonZeroUsize,
) -> Result<StatsSummary, IoError> {
    check_request(o)?;
    let name = input.reader().name().to_string();
    run_stats(
        input,
        &StatsRun {
            select: &o.select,
            prefix: &o.out_prefix,
            fields: selection_fields(&o.select),
            need_genotypes: true,
            threads,
            command: COMMAND,
        },
        &|samples, _| {
            if samples.is_empty() {
                return Err(usage(format!(
                    "vcf popgen computes its statistics over the selected samples, but {name} \
                     has no sample left to use: it has no sample columns, or --samples and \
                     --exclude-samples removed them all"
                )));
            }
            build_stats(&o.request, samples, o.max_mem)
        },
    )
}
