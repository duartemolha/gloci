//! Interval statistics (`gloci stats`): `jaccard`, `fisher`; `reldist` next.

mod fisher;
mod jaccard;
mod merged;

pub use fisher::{
    FisherCounts, FisherOptions, FisherTable, fisher_counts, fisher_table, kt_fisher_exact,
    write_fisher,
};
pub use jaccard::{JaccardOptions, JaccardResult, jaccard, write_jaccard};
pub use merged::{DropUnknown, MergedReader, Tally, TallyCounts};
