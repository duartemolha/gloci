//! Commands over BEDPE pairs (D146): `pair-intersect` (D154) and, later,
//! `pair-pair`. Each end of a pair is queried against an in-memory
//! [`IntervalIndex`](gloci_engine::IntervalIndex).

mod intersect;

pub use intersect::{PairIntersectOptions, PairSummary, PairType, pair_intersect, validate};
