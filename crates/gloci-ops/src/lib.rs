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
pub mod sort;
pub mod subtract;
pub mod transform;
pub mod vcf_filter;
pub mod view;
pub mod windows;

pub use bedsearch::{ChromBlock, find_chrom_block};
pub use check::{CheckReport, Problems, check_file};
pub use relate::{
    Acceptor, RelateInputs, RelateOptions, ReportOutput, b_placeholder, open_report,
    rewrite_coords, run_relation,
};
pub use sort::{
    Dedup, MIN_MAX_MEM, SortConfig, SortSettings, SortStats, Sorter, Tie, check_definitions,
    check_samples, merge_sorted, sort_input_options, sort_readers, sort_together,
};
pub use vcf_filter::{FilterSummary, VcfFilterOptions, vcf_filter};
pub use view::{ViewOptions, open_view, read_regions_file, run_view};
