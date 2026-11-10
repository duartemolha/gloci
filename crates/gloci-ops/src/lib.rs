//! gloci-ops: commands as library functions.

pub mod agg;
pub mod aln_cat;
pub mod bedsearch;
pub mod check;
pub mod closest;
pub mod cluster;
pub mod complement;
pub mod convert;
pub mod coverage;
pub mod genomecov;
pub mod getfasta;
pub mod groupby;
pub mod intersect;
pub mod map;
pub mod maskfasta;
pub mod merge;
pub mod nuc;
pub mod partition;
pub mod relate;
pub mod seq;
pub mod shift;
mod skipped;
pub mod sort;
pub mod spacing;
pub mod subtract;
pub mod transform;
pub mod vcf_concat;
pub mod vcf_extract;
pub mod vcf_filter;
pub mod vcf_select;
pub mod vcf_stats;
pub mod view;
pub mod windows;

pub use aln_cat::{cat_bams, cat_header};
pub use bedsearch::{ChromBlock, find_chrom_block, find_chrom_block_mapped};
pub use check::{CheckReport, Problems, check_file};
pub use gloci_engine::all_indexed;
pub use maskfasta::{Mask, MaskSet, MaskSummary, PastEnd, maskfasta, maskfasta_with};
pub use relate::{
    Acceptor, InputPaths, RelateInputs, RelateOptions, RelationFiles, ReportOutput, RunOptions,
    b_placeholder, open_report, rewrite_coords, run_relation, runs_per_chromosome,
};
pub use sort::{
    Dedup, MIN_MAX_MEM, SortConfig, SortSettings, SortStats, Sorter, Tie, check_definitions,
    check_definitions_pair, check_samples, check_samples_pair, merge_fan_in, merge_sorted,
    merge_sorted_batches, sort_input_options, sort_readers, sort_together,
};
pub use vcf_concat::{ConcatOptions, vcf_concat};
pub use vcf_extract::{
    ExtractMode, ExtractOptions, ExtractSummary, check_keys, check_tmp_dir, extract_matrix,
    extract_table,
};
pub use vcf_filter::{
    FilterSummary, VcfFilterInput, VcfFilterOptions, open_vcf_filter, vcf_filter, vcf_filter_file,
    vcf_filter_output,
};
pub use vcf_select::{
    FilterExpr, Hit, SelectStage, SelectedInput, VcfSelect, open_selected, run_selection,
    select_regions, validate_thresholds, vcf_header,
};
pub use vcf_stats::{StatsOptions, StatsRequest, StatsSummary, check_request, vcf_stats};
pub use view::{ViewOptions, open_view, read_regions_file, run_view, select_a};
