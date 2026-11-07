//! gloci-ops: commands as library functions.

pub mod sort;
pub mod view;

pub use sort::{Dedup, SortConfig, SortStats, Sorter, Tie, merge_sorted};
pub use view::{ViewOptions, run_view};
