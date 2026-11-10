//! Genotype decoding (VCF text and BCF), per-site statistics and the exact
//! Hardy-Weinberg test, shared by the `vcf` commands.

mod bcf_fmt;
mod expr;
mod extract;
mod filter;
mod gt;
mod hwe;
mod stats;

pub use bcf_fmt::{BcfGenotypes, bcf_genotypes};
pub use expr::{EvalScratch, Expr, MaskedGenotypes, parse as parse_expr};
pub use extract::{
    InfoKey, InfoKind, format_cell, format_declared, format_index, genotype_012, info_cell,
    info_kind, split_columns,
};
pub use filter::{
    BcfSite, PositionSet, RewritePlan, SampleSelection, SiteFilters, SiteView, Thinner, VcfType,
    fmt_g6, read_ids, read_positions, read_samples, rewrite_header, rewrite_text, select_samples,
};
pub use gt::{Genotype, GenotypeSource, SampleFields, TextGenotypes, text_genotypes};
pub use hwe::{HweTails, hwe_p, hwe_tails};
pub use stats::{
    FilterSummary, Freq, GenoDepth, Hardy, IndelHist, MaskRules, MissingSite, Out, Singletons,
    SiteDepth, SiteMeanDepth, SiteQuality, SiteStats, SnpDensity, Stat, StatSink, StatWrite,
    StatsEngine, TsTv, site_stats, site_stats_into,
};
