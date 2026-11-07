//! gloci-ops: commands as library functions.

pub mod sort;

pub use sort::{Dedup, SortConfig, SortStats, Sorter, Tie, merge_sorted};
