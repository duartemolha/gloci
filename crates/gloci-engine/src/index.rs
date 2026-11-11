//! The in-memory interval index (D61, D154): B records of any number of
//! chromosomes, queried by an interval of any chromosome, in any order.
//! `--unsorted` and the `pair-*` commands query it; its rules are the
//! sweep's.
//!
//! - **Layout.** One [`Lapper`] per chromosome, in a `Vec` indexed by
//!   [`ChromId`] (a chromosome with no entry holds `None`), so a query costs
//!   one array access however many contigs the reference has.
//! - **Storage.** An entry `[start, end)` is stored as `[start, end + 1)`
//!   (saturating): the closed range `[start, end]`, which a point never
//!   empties. An entry starting at `u64::MAX` has no such range and is kept
//!   in a side list. Each stored interval carries its exact `end`.
//! - **Query.** The closed range `[a.start - pad - 1, a.end + pad + 1]`
//!   (saturating) finds every entry whose padded span can overlap A, and
//!   each candidate is re-tested with spec section 4.3's exact rule
//!   ([`gloci_core::overlaps`]), so points and bookended intervals behave
//!   as they do in sorted mode. Ids come back ascending.

use gloci_core::{ChromId, Interval, Strand, Window, overlaps};
use rust_lapper::{Interval as LapInterval, Lapper};

/// One entry of the index: where it is, and the id the caller gets back.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IndexEntry {
    pub chrom: ChromId,
    pub span: Interval,
    pub id: u32,
}

/// A stored interval's value: its exact `end` and its id.
type Val = (u64, u32);

#[derive(Debug)]
struct PerChrom {
    lapper: Lapper<u64, Val>,
    /// Entries starting at `u64::MAX`, which `lapper` cannot hold:
    /// `(start, end, id)`.
    at_max: Vec<(u64, u64, u32)>,
}

/// Entries per chromosome, queried by interval (see the module docs).
#[derive(Debug, Default)]
pub struct IntervalIndex {
    chroms: Vec<Option<PerChrom>>,
}

impl IntervalIndex {
    /// Index `entries`, in any order.
    pub fn build(mut entries: Vec<IndexEntry>) -> IntervalIndex {
        entries.sort_by_key(|e| e.chrom);
        let n = entries.last().map_or(0, |e| e.chrom.0 as usize + 1);
        let mut chroms: Vec<Option<PerChrom>> = Vec::new();
        chroms.resize_with(n, || None);
        for run in entries.chunk_by(|x, y| x.chrom == y.chrom) {
            let (mut items, mut at_max) = (Vec::with_capacity(run.len()), Vec::new());
            for e in run {
                if e.span.start == u64::MAX {
                    at_max.push((e.span.start, e.span.end, e.id));
                } else {
                    items.push(LapInterval {
                        start: e.span.start,
                        stop: e.span.end.saturating_add(1),
                        val: (e.span.end, e.id),
                    });
                }
            }
            chroms[run[0].chrom.0 as usize] = Some(PerChrom {
                lapper: Lapper::new(items),
                at_max,
            });
        }
        IntervalIndex { chroms }
    }

    /// Replace `out` with the ids, ascending, of the entries on `chrom` that
    /// overlap `a` under section 4.3 (`inclusive_points` relaxes the point
    /// rule). With a `window`, each entry is padded first, as
    /// [`Window::apply`] pads it for a feature on `a_strand`.
    pub fn query(
        &self,
        chrom: ChromId,
        a: Interval,
        window: Option<Window>,
        a_strand: Strand,
        inclusive_points: bool,
        out: &mut Vec<u32>,
    ) {
        out.clear();
        let Some(Some(c)) = self.chroms.get(chrom.0 as usize) else {
            return;
        };
        let pad = window.map_or(0, |w| w.left.max(w.right));
        let lo = a.start.saturating_sub(pad).saturating_sub(1);
        let hi = a.end.saturating_add(pad).saturating_add(1);
        let hit = |b: Interval| match window {
            None => overlaps(a, b, inclusive_points),
            Some(w) => overlaps(a, w.apply(b, a_strand), inclusive_points),
        };
        // Lapper takes `[lo, hi + 1)`; `hi + 1` saturates only for entries
        // starting at `u64::MAX`, which are the side list.
        for iv in c.lapper.find(lo, hi.saturating_add(1)) {
            let (end, id) = iv.val;
            if hit(Interval {
                start: iv.start,
                end,
            }) {
                out.push(id);
            }
        }
        if hi == u64::MAX {
            for &(start, end, id) in &c.at_max {
                if hit(Interval { start, end }) {
                    out.push(id);
                }
            }
        }
        out.sort_unstable();
    }

    /// Whether the index holds nothing.
    pub fn is_empty(&self) -> bool {
        self.chroms.iter().all(Option::is_none)
    }
}
