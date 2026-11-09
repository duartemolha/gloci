//! gloci-engine: the core intersection primitive (spec §8.2). For each A
//! record it yields every B record whose footprint overlaps A, tagged by
//! its source file. Every relation command is a reporting layer over it.
//! The sorted sweep streams every input; unsorted mode (`--unsorted`)
//! loads B into per-chromosome indexes and streams A in any order.

mod active;
pub mod chunk;
pub mod hit;
pub mod kmerge;
mod nearest;
pub mod parallel;
pub mod sweep;
mod unsorted;

pub use chunk::ChunkSink;
pub use hit::{ARec, Hit, Neighbours, Visitor};
pub use kmerge::{kmerge, kmerge_with};
pub use parallel::{all_indexed, run_chunks};
#[doc(hidden)]
pub use sweep::SweepStats;
pub use sweep::{Inputs, NearestConfig, SweepConfig, SweepSummary, sweep_sorted};
pub use unsorted::sweep_unsorted;
