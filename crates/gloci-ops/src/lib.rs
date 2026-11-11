//! gloci-ops: commands as library functions.

pub mod agg;
pub mod aln_cat;
pub mod aln_flags;
pub mod aln_merge;
pub mod aln_sort;
pub mod bedsearch;
pub mod check;
pub mod closest;
pub mod cluster;
pub mod complement;
pub mod convert;
pub mod coverage;
pub mod dict;
pub mod faidx;
pub mod flagstat;
pub mod genomecov;
pub mod getfasta;
pub mod groupby;
pub mod idxstats;
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
pub use aln_merge::{
    merge_alignment_headers, merge_alignments, merge_alignments_lazy, threads_per_reader,
};
pub use aln_sort::{AlnSortOptions, check_alignments, sort_alignments};
pub use bedsearch::{ChromBlock, find_chrom_block, find_chrom_block_mapped};
pub use check::{CheckReport, Problems, check_file};
pub use gloci_engine::all_indexed;
pub use flagstat::{FlagStats, FlagstatFormat, flagstat, percent, render_flagstat};
pub use idxstats::{IdxstatsPath, idxstats};
pub use maskfasta::{Mask, MaskSet, MaskSummary, PastEnd, maskfasta, maskfasta_with};
pub use relate::{
    Acceptor, InputPaths, RelateInputs, RelateOptions, RelationFiles, ReportOutput, RunOptions,
    b_placeholder, open_report, rewrite_coords, run_relation, runs_per_chromosome,
};
pub use sort::{
    Dedup, MIN_MAX_MEM, Order, SortConfig, SortSettings, SortStats, Sorter, Tie, check_definitions,
    check_definitions_pair, check_samples, check_samples_pair, coordinate_key, merge_fan_in,
    merge_sorted, merge_sorted_batches, name_flag_rank, sort_input_options, sort_readers,
    sort_together, strnum_cmp,
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
pub use view::{ViewOptions, for_each_region, open_view, read_regions_file, run_view, select_a};
