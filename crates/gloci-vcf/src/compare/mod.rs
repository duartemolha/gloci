//! `vcf compare`: the walk over two sorted VCF streams and the genotype
//! concordance of the records it pairs (decision D163).

pub mod genotypes;
pub mod walk;

pub use genotypes::{AlleleMap, DiscordanceMatrix, Gt, GtCompare, MatrixCell};
pub use walk::{MatchRule, Row, Side, Site, SiteClass, SiteSource, Walk, share_allele};
