//! `gloci vcf popgen`: vcftools' population genetics statistics in one pass
//! over the records that pass the shared selection stage, through the
//! `vcf stats` engine (decisions D158, D159). Each statistic writes
//! `PREFIX.<extension>` with vcftools' name and columns; the files are
//! published together only when the whole run succeeds.

use crate::vcf_select::{SelectedInput, VcfSelect};
use crate::vcf_stats::{StatsSummary, run_stats, selection_fields};
use gloci_io::{IoError, UsageError};
use gloci_vcf::Stat;
use gloci_vcf::popgen::{SitePi, TajimaD, WindowPi};
use std::num::NonZeroUsize;
use std::path::PathBuf;

const COMMAND: &str = "vcf popgen";

/// Which statistics to compute.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PopgenRequest {
    /// `PREFIX.sites.pi` (vcftools `--site-pi`).
    pub site_pi: bool,
    /// `PREFIX.windowed.pi` (vcftools `--window-pi`): window size and step.
    pub window_pi: Option<(u64, u64)>,
    /// `PREFIX.Tajima.D` (vcftools `--TajimaD`): the bin size.
    pub tajima_d: Option<u64>,
}

impl PopgenRequest {
    /// True when no statistic is requested.
    pub fn is_empty(&self) -> bool {
        *self == PopgenRequest::default()
    }
}

/// Options for `vcf popgen`.
#[derive(Clone, Debug)]
pub struct PopgenOptions {
    pub select: VcfSelect,
    pub request: PopgenRequest,
    /// Output files are `PREFIX.<extension>`.
    pub out_prefix: PathBuf,
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
            "no statistic requested; add at least one of --site-pi, --window-pi SIZE or \
             --tajima-d BIN",
        ));
    }
    if let Some((size, step)) = r.window_pi {
        if size == 0 {
            return Err(usage(
                "--window-pi needs a window size above 0; give the window width in bases, such \
                 as --window-pi 10000",
            ));
        }
        if step == 0 {
            return Err(usage(
                "--window-pi-step needs a step above 0; give it in bases, or leave it out for \
                 windows that do not overlap",
            ));
        }
        if step > size {
            return Err(usage(format!(
                "--window-pi-step {step} is larger than the window ({size}), which would leave \
                 positions in no window; give a step of at most {size}"
            )));
        }
    }
    if r.tajima_d == Some(0) {
        return Err(usage(
            "--tajima-d needs a bin size above 0; give the bin width in bases, such as \
             --tajima-d 10000",
        ));
    }
    if o.out_prefix.as_os_str().is_empty() {
        return Err(usage(
            "--out-prefix is empty; give a prefix such as --out-prefix results/cohort",
        ));
    }
    Ok(())
}

/// The statistics a request asks for, over `n_samples` selected samples.
fn build_stats(r: &PopgenRequest, n_samples: u64) -> Vec<Box<dyn Stat>> {
    let mut v: Vec<Box<dyn Stat>> = Vec::new();
    if r.site_pi {
        v.push(Box::new(SitePi::new()));
    }
    if let Some((size, step)) = r.window_pi {
        v.push(Box::new(WindowPi::new(size, step)));
    }
    if let Some(bin) = r.tajima_d {
        v.push(Box::new(TajimaD::new(bin, n_samples)));
    }
    v
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
        &o.select,
        &o.out_prefix,
        selection_fields(&o.select),
        true,
        threads,
        COMMAND,
        &|samples, _| {
            if samples.is_empty() {
                return Err(usage(format!(
                    "vcf popgen computes its statistics over the selected samples, but {name} \
                     has no sample left to use: it has no sample columns, or --samples and \
                     --exclude-samples removed them all"
                )));
            }
            Ok(build_stats(&o.request, samples.len() as u64))
        },
    )
}
