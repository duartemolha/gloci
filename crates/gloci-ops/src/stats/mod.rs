//! Interval statistics (`gloci stats`): `jaccard` now, `fisher` and
//! `reldist` next.

mod jaccard;
mod merged;

pub use jaccard::{JaccardOptions, JaccardResult, jaccard, write_jaccard};
pub use merged::{MergedReader, Tally, TallyCounts};
