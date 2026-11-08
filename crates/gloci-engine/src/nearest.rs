//! Nearest-neighbour candidates for `closest` (spec §4.4, D55): for each A
//! and B file, the `k` nearest B on each side that do not overlap A, plus
//! ties. [`Neighbours`](crate::Neighbours) states the contract.
//!
//! - **Active B.** The scan that finds A's hits also takes each live B
//!   that does not overlap A as a left (`end <= a.start`) or right
//!   (`start >= a.end`) candidate.
//! - **Left pool.** A B that leaves the active set (or is skipped on load)
//!   ends before A's key. So it is a left neighbour of this A and of every
//!   later one, whose spans all start at or after this key, and it
//!   overlaps none of them. The pool keeps such B per file, and per strand
//!   under a strand rule. It drops a B once `k` members end strictly later
//!   (they are strictly nearer for every later A), or once
//!   `end + max_distance < a.key`. It therefore holds the `k` latest ends
//!   and their ties.
//! - **Right look-ahead.** Before the scan, B are read past the load limit
//!   into the active set until no later B can be among the `k` nearest on
//!   the right. Every later B starts at or after its key, so once `k`
//!   candidates are known, a B whose key is beyond the k-th smallest start
//!   ends the look-ahead (it stays unread). A key beyond
//!   `a.end + max_distance` ends it too, as does the end of the chromosome
//!   or the input. A B with the same start as the k-th is a tie, so it is
//!   read. The test is on keys, not starts, because a VCF B's span can
//!   start after its key, so B in key order are not in start order.
//!
//! Each side is then sorted by `(distance, input order)` and cut to the
//! `k` nearest plus ties.
//!
//! With `split_by_b_strand` (and no strand rule), B form two groups, `-`
//! and the rest (`+` and unknown), each with its own pool and look-ahead
//! count, and the cut keeps the `k` nearest plus ties of each group.

use crate::hit::Hit;
use crate::sweep::NearestConfig;
use gloci_core::{Interval, Strand, StrandRule, distance, overlaps};
use std::collections::BinaryHeap;
use std::ops::Range;

/// Which B strands count as neighbours of the current A.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Want {
    /// None: a strand rule, and A's strand is unknown.
    Nothing,
    /// Every B, whatever its strand (no strand rule).
    Any,
    /// B on this (known) strand.
    Only(Strand),
}

impl Want {
    fn new(rule: StrandRule, a: Strand) -> Want {
        match (rule, a) {
            (StrandRule::Any, _) => Want::Any,
            (_, Strand::Unknown) => Want::Nothing,
            (StrandRule::Same, s) => Want::Only(s),
            (StrandRule::Opposite, Strand::Plus) => Want::Only(Strand::Minus),
            (StrandRule::Opposite, _) => Want::Only(Strand::Plus),
        }
    }

    fn meets(self, b: Strand) -> bool {
        match self {
            Want::Nothing => false,
            Want::Any => true,
            Want::Only(s) => b == s,
        }
    }

    /// The [`Buckets`] whose B count for this A.
    pub(crate) fn buckets(self) -> Range<usize> {
        match self {
            // Bucket 1 is empty unless the split is on.
            Want::Any => 0..2,
            Want::Only(Strand::Plus) => 0..1,
            Want::Only(Strand::Minus) => 1..2,
            Want::Nothing | Want::Only(Strand::Unknown) => 0..0,
        }
    }
}

/// Which of two buckets a B is kept in, by its strand (the left pools'
/// sub-pools, unsorted mode's sorted lists). Under a strand rule, `+` B go
/// to 0, `-` B to 1 and B of unknown strand nowhere (D55). Otherwise, with
/// `split_by_b_strand`, `-` B go to 1 and the rest to 0; with neither,
/// every B goes to 0.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Buckets {
    rule: bool,
    split: bool,
}

impl Buckets {
    pub(crate) fn new(cfg: &NearestConfig) -> Buckets {
        Buckets {
            rule: cfg.strand != StrandRule::Any,
            split: splits(cfg),
        }
    }

    pub(crate) fn of(self, b: Strand) -> Option<usize> {
        if self.rule && b == Strand::Unknown {
            // D55: under a strand rule a B of unknown strand never counts.
            return None;
        }
        Some(group(self.rule || self.split, b))
    }
}

/// The group of a B: with `split`, `-` B are group 1 and the rest group
/// 0; otherwise every B is group 0.
fn group(split: bool, b: Strand) -> usize {
    usize::from(split && b == Strand::Minus)
}

/// `split_by_b_strand` takes effect only without a strand rule.
fn splits(cfg: &NearestConfig) -> bool {
    cfg.split_by_b_strand && cfg.strand == StrandRule::Any
}

/// The current A, as the neighbour search sees it.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Target {
    pub(crate) span: Interval,
    pub(crate) want: Want,
    pub(crate) k: usize,
    max_distance: Option<u64>,
    inclusive_points: bool,
    /// `-` B are kept apart from the rest (`split_by_b_strand`).
    split: bool,
}

impl Target {
    /// `cfg.k` must be at least 1.
    pub(crate) fn new(
        cfg: &NearestConfig,
        span: Interval,
        strand: Strand,
        inclusive_points: bool,
    ) -> Target {
        debug_assert!(cfg.k > 0);
        Target {
            span,
            want: Want::new(cfg.strand, strand),
            k: cfg.k,
            max_distance: cfg.max_distance,
            inclusive_points,
            split: splits(cfg),
        }
    }

    /// How many groups of B are kept apart: 1, or 2 with the split.
    fn groups(&self) -> usize {
        1 + usize::from(self.split)
    }

    /// Does `b` overlap A (spec §4.3)?
    pub(crate) fn overlaps(&self, b: Interval) -> bool {
        overlaps(self.span, b, self.inclusive_points)
    }

    pub(crate) fn distance(&self, b: Interval) -> u64 {
        distance(self.span, b, self.inclusive_points)
    }

    /// Is `b` at most `max_distance` from A?
    pub(crate) fn within(&self, b: Interval) -> bool {
        self.max_distance.is_none_or(|m| self.distance(b) <= m)
    }
}

/// What `ActiveSet::scan` gathers for `closest`: evicted B go to the pool,
/// and live B that are neighbours of A to `left` or `right`.
pub(crate) struct NearScan<'s, 'h> {
    pub(crate) t: Target,
    pub(crate) a_key: u64,
    /// Past the load limit the scan looks at entries while
    /// `key <= right_bound` (from the look-ahead).
    pub(crate) right_bound: u64,
    pub(crate) pool: &'s mut LeftPool,
    pub(crate) left: &'s mut Vec<Hit<'h>>,
    pub(crate) right: &'s mut Vec<Hit<'h>>,
}

impl<'h> NearScan<'_, 'h> {
    /// A B left the active set.
    pub(crate) fn evicted(&mut self, span: Interval, strand: Strand, loc: u64, raw: &[u8]) {
        self.pool.offer(span, strand, loc, raw, self.a_key);
    }

    /// A live B: keep it if it is a neighbour of A.
    pub(crate) fn live(&mut self, h: Hit<'h>) {
        let t = &self.t;
        if !t.want.meets(h.strand)
            || overlaps(t.span, h.span, t.inclusive_points)
            || !t.within(h.span)
        {
            return;
        }
        if h.span.end <= t.span.start {
            self.left.push(h);
        } else if h.span.start >= t.span.end {
            self.right.push(h);
        }
    }
}

/// The B of one input that left the active set, kept while they can still
/// be among the `k` nearest left of a later A.
#[derive(Debug)]
pub(crate) struct LeftPool {
    k: usize,
    max_distance: Option<u64>,
    /// Which of `subs` each B goes to.
    buckets: Buckets,
    subs: [SubPool; 2],
    /// Scratch space for compacting an arena.
    spare: Vec<u8>,
    /// The most arena bytes (both sub-pools) held at any moment.
    peak_bytes: usize,
}

#[derive(Debug, Default)]
struct SubPool {
    /// Latest end first, then input order.
    members: Vec<Member>,
    /// Members' raw bytes, `arena[off..off + len]`.
    arena: Vec<u8>,
    /// Arena bytes of dropped members.
    dead: usize,
}

#[derive(Debug, Clone, Copy)]
struct Member {
    span: Interval,
    strand: Strand,
    loc: u64,
    off: usize,
    len: u32,
}

impl LeftPool {
    pub(crate) fn new(cfg: &NearestConfig) -> LeftPool {
        LeftPool {
            k: cfg.k,
            max_distance: cfg.max_distance,
            buckets: Buckets::new(cfg),
            subs: Default::default(),
            spare: Vec::new(),
            peak_bytes: 0,
        }
    }

    /// The most raw bytes this pool held at any moment.
    pub(crate) fn peak_bytes(&self) -> usize {
        self.peak_bytes
    }

    /// Members, over both strands.
    pub(crate) fn len(&self) -> usize {
        self.subs.iter().map(|s| s.members.len()).sum()
    }

    pub(crate) fn clear(&mut self) {
        for s in &mut self.subs {
            s.members.clear();
            s.arena.clear();
            s.dead = 0;
        }
    }

    /// Drop the members with `end + max_distance < a_key`: too far from
    /// this A and from every later one.
    pub(crate) fn expire(&mut self, a_key: u64) {
        let Some(max) = self.max_distance else {
            return;
        };
        for s in &mut self.subs {
            let keep = s
                .members
                .partition_point(|m| m.span.end.saturating_add(max) >= a_key);
            s.drop_from(keep, &mut self.spare);
        }
    }

    /// Take a B that left the active set while A's key was `a_key` (it
    /// ends before it), unless it cannot be a neighbour of a later A.
    pub(crate) fn offer(
        &mut self,
        span: Interval,
        strand: Strand,
        loc: u64,
        raw: &[u8],
        a_key: u64,
    ) {
        if self
            .max_distance
            .is_some_and(|m| span.end.saturating_add(m) < a_key)
        {
            return;
        }
        let Some(i) = self.buckets.of(strand) else {
            return;
        };
        if let Some(bytes) = self.subs[i].offer(span, strand, loc, raw, self.k, &mut self.spare) {
            let held = bytes + self.subs[1 - i].arena.len();
            self.peak_bytes = self.peak_bytes.max(held);
        }
    }

    /// Push the members of the strands `t` wants that are within its
    /// `max_distance`, as hits of `file`.
    pub(crate) fn extend<'h>(&'h self, file: u16, t: &Target, out: &mut Vec<Hit<'h>>) {
        for sub in &self.subs[t.want.buckets()] {
            for m in &sub.members {
                debug_assert!(
                    m.span.end <= t.span.start && !overlaps(t.span, m.span, t.inclusive_points),
                    "a pool member ends before A's key, so before A"
                );
                if !t.within(m.span) {
                    // Latest end first: the rest are farther.
                    break;
                }
                out.push(Hit {
                    file,
                    span: m.span,
                    strand: m.strand,
                    raw: &sub.arena[m.off..m.off + m.len as usize],
                    loc: m.loc,
                });
            }
        }
    }
}

impl SubPool {
    /// Insert in order unless `k` members end strictly later, then drop
    /// the members that now have `k` ending strictly later. Returns the
    /// arena's size just after the insertion, if any.
    fn offer(
        &mut self,
        span: Interval,
        strand: Strand,
        loc: u64,
        raw: &[u8],
        k: usize,
        spare: &mut Vec<u8>,
    ) -> Option<usize> {
        let end = span.end;
        let later = self.members.partition_point(|m| m.span.end > end);
        if later >= k {
            return None;
        }
        let at =
            later + self.members[later..].partition_point(|m| m.span.end == end && m.loc < loc);
        let off = self.arena.len();
        self.arena.extend_from_slice(raw);
        self.members.insert(
            at,
            Member {
                span,
                strand,
                loc,
                off,
                len: raw.len() as u32,
            },
        );
        let bytes = self.arena.len();
        if self.members.len() > k {
            let kth = self.members[k - 1].span.end;
            let keep = k + self.members[k..].partition_point(|m| m.span.end == kth);
            self.drop_from(keep, spare);
        }
        Some(bytes)
    }

    /// Drop the members from `keep` on. Once dropped bytes are more than
    /// half the arena, move the members' bytes to the front (through
    /// `spare`). Each compaction copies fewer bytes than were dropped
    /// since the last, so it costs O(1) per byte, and the arena never
    /// holds much more than twice its members, even while one A step
    /// skips a whole gap of B into the pool.
    fn drop_from(&mut self, keep: usize, spare: &mut Vec<u8>) {
        for m in &self.members[keep..] {
            self.dead += m.len as usize;
        }
        self.members.truncate(keep);
        if self.dead * 2 <= self.arena.len() {
            return;
        }
        spare.clear();
        for m in &mut self.members {
            let n = m.len as usize;
            let off = spare.len();
            spare.extend_from_slice(&self.arena[m.off..m.off + n]);
            m.off = off;
        }
        std::mem::swap(&mut self.arena, spare);
        self.dead = 0;
    }
}

/// While reading ahead, the `k` smallest starts of the right neighbours of
/// A found past the load limit, per group of B, which say when to stop.
pub(crate) struct LookAhead<'s> {
    t: Target,
    /// Per group, a max-heap of the `k` smallest starts seen.
    starts: &'s mut [BinaryHeap<u64>; 2],
    /// No B starting beyond this is within `max_distance`.
    cap: u64,
}

impl<'s> LookAhead<'s> {
    pub(crate) fn new(t: Target, starts: &'s mut [BinaryHeap<u64>; 2]) -> LookAhead<'s> {
        starts.iter_mut().for_each(BinaryHeap::clear);
        let cap = t
            .max_distance
            .map_or(u64::MAX, |m| t.span.end.saturating_add(m));
        LookAhead { t, starts, cap }
    }

    /// The largest key a B can have and still be among the `k` nearest
    /// right of A in its group. A B with a larger key starts beyond the
    /// k-th smallest start already seen in every group, or beyond
    /// `max_distance`.
    pub(crate) fn bound(&self) -> u64 {
        let kth = self.starts[..self.t.groups()]
            .iter()
            .map(|h| match h.peek() {
                Some(&s) if h.len() == self.t.k => s,
                _ => u64::MAX,
            })
            .max()
            .unwrap_or(u64::MAX);
        kth.min(self.cap)
    }

    /// A B past the load limit, which starts after A's end: a right
    /// neighbour if its strand counts.
    pub(crate) fn add(&mut self, start: u64, strand: Strand) {
        if !self.t.want.meets(strand) {
            return;
        }
        let starts = &mut self.starts[group(self.t.split, strand)];
        if starts.len() < self.t.k {
            starts.push(start);
        } else if let Some(mut top) = starts.peek_mut()
            && start < *top
        {
            *top = start;
        }
    }
}

/// Sort `v[from..]` nearest first, ties in input order, and keep, per
/// group, the `k` nearest and every B tied with the k-th.
pub(crate) fn keep_nearest(v: &mut Vec<Hit<'_>>, from: usize, t: &Target) {
    v[from..].sort_unstable_by_key(|h| (t.distance(h.span), h.loc));
    // Per group: how many are kept, and the k-th distance once reached.
    let mut kept = [0; 2];
    let mut kth = [None; 2];
    let mut w = from;
    for r in from..v.len() {
        let h = v[r];
        let (g, d) = (group(t.split, h.strand), t.distance(h.span));
        if kept[g] < t.k {
            kept[g] += 1;
            if kept[g] == t.k {
                kth[g] = Some(d);
            }
        } else if kth[g] != Some(d) {
            continue;
        }
        v[w] = h;
        w += 1;
    }
    v.truncate(w);
}
