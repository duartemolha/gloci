//! gloci-engine: the core intersection primitive (spec §8.2). For each A
//! record it yields every B record whose footprint overlaps A, tagged by
//! its source file. Every relation command is a reporting layer over it.

mod active;
pub mod hit;
mod nearest;
pub mod sweep;

pub use hit::{ARec, Hit, Neighbours, Visitor};
#[doc(hidden)]
pub use sweep::SweepStats;
pub use sweep::{Inputs, NearestConfig, SweepConfig, SweepSummary, sweep_sorted};
