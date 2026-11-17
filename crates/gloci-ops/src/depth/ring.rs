//! The per-input depth counts of `depth` and `bedcov` (D137): a ring of counters keyed by
//! absolute position, so memory follows the longest read and not the reference.

use gloci_core::Interval;
use gloci_io::AlnRef;

/// What counts toward the depth of a position.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CountRules {
    /// `D` bases count (`depth -J`), each gated by the quality of the next query base.
    pub deletions: bool,
    /// `N` bases count (`bedcov` without `-j`); `depth` never counts them.
    pub ref_skips: bool,
    /// A base counts when its quality is at least this (0 = every base). A missing quality
    /// (0xFF) passes, and so does every base of a read with SEQ `*`.
    pub min_baseq: u8,
}

/// Add 1 to the counter of each position in `[lo, hi)`, which the ring holds (`hi - lo` is at
/// most its length): the run is one or two slices of `counts`, counted without a mask.
fn count_range(counts: &mut [u32], lo: u64, hi: u64) {
    if lo >= hi {
        return;
    }
    let len = counts.len();
    let at = (lo as usize) & (len - 1);
    let n = ((hi - lo) as usize).min(len);
    let first = n.min(len - at);
    for c in &mut counts[at..at + first] {
        *c = c.saturating_add(1);
    }
    for c in &mut counts[..n - first] {
        *c = c.saturating_add(1);
    }
}

/// Call `f` with each run of `[lo, hi)` that lies outside `except` (sorted, disjoint).
fn for_each_gap(lo: u64, hi: u64, except: &[Interval], mut f: impl FnMut(u64, u64)) {
    let mut at = lo;
    for iv in &except[except.partition_point(|iv| iv.end <= lo)..] {
        if iv.start >= hi {
            break;
        }
        if at < iv.start {
            f(at, iv.start);
        }
        at = at.max(iv.end);
    }
    if at < hi {
        f(at, hi);
    }
}

/// The counts of one input's reads, from `base` to `covered_end`. A position at or past
/// `covered_end` has no read over it, and so depth 0.
#[derive(Debug, Default)]
pub struct DepthRing {
    counts: Vec<u32>,
    /// First position still held; every earlier one has been released.
    base: u64,
    /// One past the last position any read added so far covers.
    end: u64,
}

impl DepthRing {
    pub fn new() -> DepthRing {
        DepthRing::default()
    }

    /// One past the last position covered by a read added since the last [`reset`].
    ///
    /// [`reset`]: DepthRing::reset
    pub fn covered_end(&self) -> u64 {
        self.end
    }

    /// The slots held: `next_pow2(longest span + 1)` at most.
    pub fn capacity(&self) -> usize {
        self.counts.len()
    }

    /// Forget every read (a new reference starts).
    pub fn reset(&mut self) {
        self.base = 0;
        self.end = 0;
    }

    /// Release the positions before `pos`: reads arrive in order, so none of them is asked
    /// for again, and a read that starts at or past the covered end leaves nothing held.
    fn release_to(&mut self, pos: u64) {
        if pos >= self.end {
            self.base = pos;
            self.end = pos;
        } else if pos > self.base {
            self.base = pos;
        }
    }

    fn grow(&mut self, need: u64) {
        let want = usize::try_from(need.saturating_add(1))
            .unwrap_or(usize::MAX >> 1)
            .next_power_of_two()
            .max(64);
        if want <= self.counts.len() {
            return;
        }
        let mut fresh = vec![0u32; want];
        let (old_mask, new_mask) = (self.counts.len().wrapping_sub(1), want - 1);
        for p in self.base..self.end {
            fresh[(p as usize) & new_mask] = self.counts[(p as usize) & old_mask];
        }
        self.counts = fresh;
    }

    /// Extend the covered range over `span`, releasing what lies before it.
    fn cover(&mut self, span: Interval) {
        self.release_to(span.start);
        let new_end = span.end.max(self.end);
        self.grow(new_end - self.base);
        let mask = self.counts.len() - 1;
        // Positions beyond the covered end hold what released positions left: zero them.
        for p in self.end..new_end {
            self.counts[(p as usize) & mask] = 0;
        }
        self.end = new_end;
    }

    /// Cover `span` without counting anything (a placed unmapped read: its position prints,
    /// with depth 0).
    pub fn add_empty(&mut self, span: Interval) {
        self.cover(span);
    }

    /// Count the read `r`, whose D101 span is `span` (cut at the reference length), by walking
    /// its CIGAR from its position.
    ///
    /// The contract callers (`depth`, `--skip-mate-overlap`, `bedcov`) rely on:
    /// - Reads arrive in order, and the caller must have taken or printed every position below
    ///   `span.start` before calling: `add` releases them (a later `take(p)` with `p` below the
    ///   ring's base returns 0).
    /// - `clip` is an absolute position: bases at positions below it are not counted (0 =
    ///   none), as for the part of a read that overlaps its mate. It does not shrink the span:
    ///   `covered_end` grows to `span.end` and the positions in `span` still print.
    pub fn add(&mut self, r: &AlnRef<'_>, span: Interval, rules: &CountRules, clip: u64) {
        self.add_inner(r, span, rules, clip, &[]);
    }

    /// As `add`, but bases at positions inside any of `except` (sorted, disjoint) are not counted.
    pub fn add_except(
        &mut self,
        r: &AlnRef<'_>,
        span: Interval,
        rules: &CountRules,
        except: &[Interval],
    ) {
        self.add_inner(r, span, rules, 0, except);
    }

    fn add_inner(
        &mut self,
        r: &AlnRef<'_>,
        span: Interval,
        rules: &CountRules,
        clip: u64,
        except: &[Interval],
    ) {
        self.cover(span);
        let mask = self.counts.len() - 1;

        let quals = r.quals();
        let seq_len = r.seq_len();
        // Every base counts: no quality gate, or none to apply (SEQ `*`).
        let ungated = rules.min_baseq == 0 || seq_len == 0;
        let gate = |q: usize| -> bool {
            ungated || {
                let c = quals.get(q).copied().unwrap_or(0);
                c == 0xFF || c >= rules.min_baseq
            }
        };
        let counts = &mut self.counts;
        let mut rpos = span.start;
        let mut q = 0usize;
        r.for_each_op(&mut |op, n| {
            let n = u64::from(n);
            let lo = rpos.max(clip);
            let hi = rpos.saturating_add(n).min(span.end);
            match op {
                b'M' | b'=' | b'X' => {
                    for_each_gap(lo, hi, except, |a, b| {
                        if ungated {
                            count_range(counts, a, b);
                        } else {
                            for p in a..b {
                                if gate(q + (p - rpos) as usize) {
                                    let s = &mut counts[(p as usize) & mask];
                                    *s = s.saturating_add(1);
                                }
                            }
                        }
                    });
                    q += n as usize;
                    rpos = rpos.saturating_add(n);
                }
                b'D' => {
                    if rules.deletions && gate(q) {
                        for_each_gap(lo, hi, except, |a, b| count_range(counts, a, b));
                    }
                    rpos = rpos.saturating_add(n);
                }
                b'N' => {
                    if rules.ref_skips {
                        for_each_gap(lo, hi, except, |a, b| count_range(counts, a, b));
                    }
                    rpos = rpos.saturating_add(n);
                }
                b'I' | b'S' => q += n as usize,
                _ => {}
            }
        });
    }

    /// The depth at `pos`; the positions before it are released.
    pub fn take(&mut self, pos: u64) -> u32 {
        if pos >= self.end {
            return 0;
        }
        if pos < self.base {
            return 0;
        }
        self.base = pos;
        self.counts[(pos as usize) & (self.counts.len() - 1)]
    }
}

#[cfg(test)]
mod tests {
    use super::count_range;

    #[test]
    fn count_range_wraps_the_ring() {
        let mut counts = vec![0u32; 8];
        count_range(&mut counts, 6, 11);
        assert_eq!(counts, [1, 1, 1, 0, 0, 0, 1, 1]);
        count_range(&mut counts, 16, 24);
        assert_eq!(counts, [2, 2, 2, 1, 1, 1, 2, 2]);
        count_range(&mut counts, 5, 5);
        assert_eq!(counts, [2, 2, 2, 1, 1, 1, 2, 2]);
        counts[1] = u32::MAX;
        count_range(&mut counts, 9, 10);
        assert_eq!(counts[1], u32::MAX);
    }
}
