//! What the sweep hands its visitor for each A record.

use gloci_core::{Interval, Strand};
use gloci_io::{IoError, RecordMeta};
use std::ops::Range;

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

/// The nearest B on each side of A, per file, filled when
/// `SweepConfig::nearest` is set (empty otherwise, or with `k = 0`).
///
/// - **Left** holds B that do not overlap A with `b.span.end <=
///   a.span.start`, **right** those with `b.span.start >= a.span.end`.
///   Each is sorted nearest first by gloci-core `distance()`, ties in input
///   order.
/// - Each side holds exactly the `k` nearest and every B tied with the
///   k-th, among the B at distance at most `max_distance`. So it contains
///   every B that could be among the `k` nearest on that side.
/// - When `strand` is not `Any`, only B whose strand meets the rule for
///   this A count. B with unknown strand never do, and an A with unknown
///   strand has no neighbours (D55).
///
/// The hits' raw bytes are valid only during the `visit` call.
#[derive(Debug, Default)]
pub struct Neighbours<'a> {
    pub(crate) left: Vec<Hit<'a>>,
    pub(crate) right: Vec<Hit<'a>>,
    /// Where each file's run ends in `left` and `right`, by file.
    pub(crate) ends: Vec<(usize, usize)>,
}

impl<'a> Neighbours<'a> {
    /// B of `file` ending at or before A's start, nearest first.
    pub fn left(&self, file: u16) -> &[Hit<'a>] {
        &self.left[self.run(file, |e| e.0)]
    }

    /// B of `file` starting at or after A's end, nearest first.
    pub fn right(&self, file: u16) -> &[Hit<'a>] {
        &self.right[self.run(file, |e| e.1)]
    }

    fn run(&self, file: u16, side: impl Fn((usize, usize)) -> usize) -> Range<usize> {
        let f = usize::from(file);
        match self.ends.get(f) {
            None => 0..0,
            Some(&end) => {
                let start = if f == 0 { 0 } else { side(self.ends[f - 1]) };
                start..side(end)
            }
        }
    }

    /// Close the current file's runs: what was pushed since the last call
    /// belongs to the next file.
    pub(crate) fn end_file(&mut self) {
        self.ends.push((self.left.len(), self.right.len()));
    }

    /// Empty, keeping the allocations, for neighbours that borrow
    /// something else.
    pub(crate) fn recycle<'b>(mut self) -> Neighbours<'b> {
        self.ends.clear();
        Neighbours {
            left: recycle(self.left),
            right: recycle(self.right),
            ends: self.ends,
        }
    }
}

/// Reuse the allocation of a hit vector for hits that borrow something
/// else: the vector is emptied, and collecting an empty `into_iter` of a
/// same-layout type reuses its buffer.
pub(crate) fn recycle<'x, 'y>(mut v: Vec<Hit<'x>>) -> Vec<Hit<'y>> {
    v.clear();
    v.into_iter()
        .map(|_| unreachable!("the vector is empty"))
        .collect()
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
