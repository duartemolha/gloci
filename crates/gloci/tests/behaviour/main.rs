//! Cross-command behaviour: thread determinism, broken pipe, truncation,
//! atomic outputs, signals and the error-message format.

mod alignments_matrices;
mod bed_vcf_matrices;
mod bed_vcf_outputs;
mod chrom_map;
#[path = "../common/mod.rs"]
mod common;
mod fixtures_alignments;
mod fixtures_bed_vcf;
mod fixtures_samtools;
mod fixtures_stats;
mod helpers;
mod matrices;
mod process;
mod samtools_matrices;
mod stats_matrices;
mod stdin;
