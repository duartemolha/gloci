//! Genotype decoding (VCF text and BCF), per-site statistics and the exact
//! Hardy-Weinberg test, shared by the `vcf` commands.

pub mod annotate;
mod bcf_counts;
mod bcf_fmt;
pub mod compare;
pub mod concat;
mod expr;
mod extract;
mod filter;
mod gt;
mod hwe;
pub mod merge;
mod norm;
pub mod popgen;
pub mod query;
mod split;
pub mod srsort;
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
    CacheReader, CachedGenotypes, Genotype, GenotypeSource, GtFields, PackedGenotypes,
    PackedRecord, RecordSlot, SampleCache, SampleFields, TextGenotypes, text_genotypes,
    text_genotypes_with,
};
pub use hwe::{HweTails, hwe_p, hwe_tails};
pub use norm::{CheckRef, Dedup, Normalizer, Reference, RmDup, Verdict, clean_reference_bases};
pub use query::{Query, QueryFormat};
pub use split::{Outcome as SplitOutcome, SplitTypes, Splitter};
pub use stats::{
    FilterSummary, Freq, GenoDepth, Hardy, Heterozygosity, IndelHist, MaskRules, MissingSample,
    MissingSite, Out, SUMMARY_SCHEMA_VERSION, SampleDepth, SampleTally, Singletons, SiteDepth,
    SiteMeanDepth, SiteQuality, SiteStats, SnpDensity, Stat, StatSink, StatWrite, StatsEngine,
    Summary, TallyWants, TsTv, site_stats, site_stats_into,
};
