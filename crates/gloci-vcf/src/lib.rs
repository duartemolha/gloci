//! Genotype decoding (VCF text and BCF), per-site statistics and the exact
//! Hardy-Weinberg test, shared by the `vcf` commands.

mod bcf_fmt;
mod filter;
mod gt;
mod hwe;
mod stats;

pub use bcf_fmt::{BcfGenotypes, bcf_genotypes};
pub use filter::{
    BcfSite, PositionSet, SiteFilters, SiteView, Thinner, VcfType, read_ids, read_positions,
};
pub use gt::{Genotype, GenotypeSource, SampleFields, TextGenotypes, text_genotypes};
pub use hwe::hwe_p;
pub use stats::{MaskRules, SiteStats, site_stats};
