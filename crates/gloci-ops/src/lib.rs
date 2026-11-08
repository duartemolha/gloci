//! gloci-ops: commands as library functions.

pub mod bedsearch;
pub mod check;
pub mod relate;
pub mod sort;
pub mod view;

pub use bedsearch::{ChromBlock, find_chrom_block};
pub use check::{CheckReport, Problems, check_file};
pub use relate::{
    RelateInputs, RelateOptions, ReportOutput, accepted, b_placeholder, open_report,
    rewrite_coords, run_relation,
};
pub use sort::{Dedup, SortConfig, SortStats, Sorter, Tie, merge_sorted};
pub use view::{ViewOptions, open_view, run_view};
