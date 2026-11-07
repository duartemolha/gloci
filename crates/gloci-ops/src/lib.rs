//! gloci-ops: commands as library functions.

pub mod bedsearch;
pub mod check;
pub mod sort;
pub mod view;

pub use bedsearch::{ChromBlock, find_chrom_block};
pub use check::{CheckReport, Problems, check_file};
pub use sort::{Dedup, SortConfig, SortStats, Sorter, Tie, merge_sorted};
pub use view::{ViewOptions, open_view, run_view};
