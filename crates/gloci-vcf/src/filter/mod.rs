//! `vcf filter`: the site filters of spec 7.8 step 4 (decision D50).

mod lists;
mod rewrite;
mod samples;
mod site;

pub use lists::{PositionSet, read_ids, read_positions};
pub use rewrite::{RewritePlan, fmt_g6, rewrite_header, rewrite_text};
pub use samples::{SampleSelection, read_samples, select_samples};
pub use site::{BcfSite, SiteView, Thinner};

use gloci_core::RegionSet;
use std::collections::{HashMap, HashSet};

/// A record type for `--types` / `--exclude-types`. A record has a type
/// set, the union of its ALT kinds; `Indel` is `Ins` or `Del`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VcfType {
    Snv,
    Mnv,
    Ins,
    Del,
    Indel,
    Sv,
    Ref,
    Other,
}

/// The site filters, ANDed. Every field's default (`None`, empty, `false`)
/// applies no filter.
///
/// A statistic that cannot be computed (no called alleles for MAF/MAC, no
/// depth for mean DP) fails every threshold on it, like a missing QUAL
/// fails `min_qual`. The exception is `hwe`: sites the test does not
/// apply to (not biallelic diploid) pass.
#[derive(Clone, Debug, Default)]
pub struct SiteFilters {
    pub chroms: Option<HashSet<Vec<u8>>>,
    pub exclude_chroms: HashSet<Vec<u8>>,
    pub positions: Option<PositionSet>,
    pub exclude_positions: Option<PositionSet>,
    pub ids: Option<HashSet<Vec<u8>>>,
    pub exclude_ids: Option<HashSet<Vec<u8>>>,
    pub types: Option<Vec<VcfType>>,
    pub exclude_types: Option<Vec<VcfType>>,
    pub min_alleles: Option<usize>,
    pub max_alleles: Option<usize>,
    pub min_qual: Option<f64>,
    pub pass_only: bool,
    pub keep_filter: Option<Vec<Vec<u8>>>,
    pub remove_filter: Option<Vec<Vec<u8>>>,
    pub require_info: Vec<Vec<u8>>,
    pub exclude_info: Vec<Vec<u8>>,
    pub min_maf: Option<f64>,
    pub max_maf: Option<f64>,
    pub min_mac: Option<u64>,
    pub max_mac: Option<u64>,
    pub max_missing: Option<f64>,
    pub max_missing_count: Option<u64>,
    pub min_mean_dp: Option<f64>,
    pub max_mean_dp: Option<f64>,
    pub hwe: Option<f64>,
    pub thin: Option<u64>,
    pub phased_only: bool,
    pub exclude_regions: Option<HashMap<Vec<u8>, RegionSet>>,
    pub inclusive_points: bool,
}
