//! Genotype decoding (VCF text and BCF), per-site statistics and the exact
//! Hardy-Weinberg test, shared by the `vcf` commands.

mod bcf_fmt;
pub mod concat;
mod expr;
mod extract;
mod filter;
mod gt;
mod hwe;
mod split;
mod stats;

pub use bcf_fmt::{BcfGenotypes, bcf_genotypes, bcf_genotypes_with};
pub use expr::{EvalScratch, Expr, MaskedGenotypes, parse as parse_expr};
pub use extract::{
    GtError, InfoKey, InfoKind, declared_keys, format_cell, format_declared, format_index,
    genotype_012, info_cell, info_kind, split_columns,
};
pub use filter::{
    BcfSite, PositionSet, RewritePlan, SampleSelection, SiteFilters, SiteView, Thinner, VcfType,
    fmt_g6, read_ids, read_positions, read_samples, rewrite_header, rewrite_text, select_samples,
};
pub use gt::{
    CacheReader, CachedGenotypes, Genotype, GenotypeSource, GtFields, SampleCache, SampleFields,
    TextGenotypes, text_genotypes, text_genotypes_with,
};
pub use hwe::{HweTails, hwe_p, hwe_tails};
pub use split::{Outcome as SplitOutcome, SplitTypes, Splitter};
pub use stats::{
    FilterSummary, Freq, GenoDepth, Hardy, Heterozygosity, IndelHist, MaskRules, MissingSample,
    MissingSite, Out, SUMMARY_SCHEMA_VERSION, SampleDepth, Singletons, SiteDepth, SiteMeanDepth,
    SiteQuality, SiteStats, SnpDensity, Stat, StatSink, StatWrite, StatsEngine, Summary, TsTv,
    site_stats, site_stats_into,
};
