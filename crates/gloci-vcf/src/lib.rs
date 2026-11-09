//! Genotype decoding (VCF text and BCF), per-site statistics and the exact
//! Hardy-Weinberg test, shared by the `vcf` commands.

mod bcf_fmt;
mod expr;
mod filter;
mod gt;
mod hwe;
mod stats;

pub use bcf_fmt::{BcfGenotypes, bcf_genotypes};
pub use expr::{EvalScratch, Expr, MaskedGenotypes, parse as parse_expr};
pub use filter::{
    BcfSite, PositionSet, RewritePlan, SampleSelection, SiteFilters, SiteView, Thinner, VcfType,
    fmt_g6, read_ids, read_positions, read_samples, rewrite_header, rewrite_text, select_samples,
};
pub use gt::{Genotype, GenotypeSource, SampleFields, TextGenotypes, text_genotypes};
pub use hwe::hwe_p;
pub use stats::{MaskRules, SiteStats, site_stats};
