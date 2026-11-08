//! The `-c COL:OP` aggregation engine (D46), shared by `map` and `merge`.
//!
//! `parse_aggs` reads the typed specs, `resolve` binds them to a file's
//! columns, and `Accumulators` folds the B records overlapping each A record.

mod num;
mod ops;
mod spec;

pub use ops::{Accumulators, ResolvedAgg, resolve};
pub use spec::{AggSpec, Col, Op, parse_aggs};
