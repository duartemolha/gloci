//! What the sweep hands its visitor for each A record.

use gloci_core::{Interval, Strand};
use gloci_io::{IoError, RecordMeta};
use std::marker::PhantomData;

/// An A record: its position data and its raw record bytes.
#[derive(Debug, Clone, Copy)]
pub struct ARec<'a> {
    pub meta: RecordMeta,
    pub raw: &'a [u8],
}

/// A candidate hit (D36): a B record whose window-padded span overlaps A's
/// span under §4.3. Thresholds and strand rules are the reporting layer's.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Hit<'a> {
    /// The B input it comes from, by position in `Inputs::b`.
    pub file: u16,
    /// B's own footprint, without the window.
    pub span: Interval,
    pub strand: Strand,
    /// B's raw record bytes, valid only during the `visit` call.
    pub raw: &'a [u8],
    /// B's 1-based line (or record) number.
    pub loc: u64,
}

/// The nearest non-overlapping B on each side of A, per file. Task 6 fills
/// them when `SweepConfig::nearest` is set; until then (and by default)
/// they are empty.
#[derive(Debug, Default)]
pub struct Neighbours<'a> {
    _hits: PhantomData<&'a [Hit<'a>]>,
}

impl<'a> Neighbours<'a> {
    /// B of `file` ending at or before A's start, nearest first.
    pub fn left(&self, _file: u16) -> &[Hit<'a>] {
        &[]
    }

    /// B of `file` starting at or after A's end, nearest first.
    pub fn right(&self, _file: u16) -> &[Hit<'a>] {
        &[]
    }
}

/// Receives every A record, in A's input order, with its hits ordered by B
/// file, then B input order (D38).
pub trait Visitor {
    fn visit(
        &mut self,
        a: &ARec<'_>,
        hits: &[Hit<'_>],
        near: &Neighbours<'_>,
    ) -> Result<(), IoError>;
}
