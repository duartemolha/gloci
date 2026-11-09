//! gloci-ops: commands as library functions.

pub mod agg;
pub mod bedsearch;
pub mod check;
pub mod closest;
pub mod complement;
pub mod convert;
pub mod coverage;
pub mod intersect;
pub mod map;
pub mod merge;
pub mod relate;
mod skipped;
pub mod sort;
pub mod subtract;
pub mod transform;
pub mod vcf_filter;
pub mod view;
pub mod windows;

pub use bedsearch::{ChromBlock, find_chrom_block, find_chrom_block_mapped};
pub use check::{CheckReport, Problems, check_file};
pub use gloci_engine::all_indexed;
pub use relate::{
    Acceptor, InputPaths, RelateInputs, RelateOptions, RelationFiles, ReportOutput, RunOptions,
    b_placeholder, open_report, rewrite_coords, run_relation, runs_per_chromosome,
};
pub use sort::{
    Dedup, MIN_MAX_MEM, SortConfig, SortSettings, SortStats, Sorter, Tie, check_definitions,
    check_definitions_pair, check_samples, check_samples_pair, merge_fan_in, merge_sorted,
    merge_sorted_batches, sort_input_options, sort_readers, sort_together,
};
pub use vcf_filter::{
    FilterExpr, FilterSummary, VcfFilterInput, VcfFilterOptions, open_vcf_filter,
    validate_thresholds, vcf_filter, vcf_filter_file, vcf_filter_output,
};
pub use view::{ViewOptions, open_view, read_regions_file, run_view, select_a};
