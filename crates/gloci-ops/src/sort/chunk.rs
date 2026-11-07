//! The in-memory chunk: record bytes in one arena plus a sort entry per
//! record, both held within the memory budget.

use super::Tie;
use rayon::slice::ParallelSliceMut;
use std::cmp::Ordering;

/// Sort entry of one record; its bytes are `arena[off..off + len]`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Entry {
    pub key: u64,
    pub end: u64,
    /// Global input sequence number (the stable tiebreak).
    pub seq: u64,
    pub off: usize,
    pub chrom: u32,
    pub len: u32,
}

pub(crate) const ENTRY_BYTES: usize = std::mem::size_of::<Entry>();

/// The size a chunk first grows to (or the budget, if smaller).
const FIRST_TARGET: usize = 64 << 10;

/// Records accumulated until the budget is reached, then sorted.
///
/// The capacity held by `arena` and `entries` together never exceeds the
/// budget, except for a single record larger than the budget, which gets a
/// chunk of its own. Growth goes through targets `budget / 2^j`, split
/// between the two vectors by the mix of record lengths seen so far; once
/// the next target would exceed `budget / 16` it jumps straight to the full
/// budget, so freed intermediate allocations stay small (at most about
/// `budget / 8` in total) instead of piling up in the allocator. A vector is
/// only moved to a larger allocation while its old copy, the copied part of
/// the new one and the other vector together fit in the budget. During the
/// jump the virtual size can transiently reach about `budget + budget / 16`
/// (the old copy plus the new allocation), while resident memory stays
/// within the budget because the new allocation is only touched as records
/// arrive.
///
/// Whether a record fits depends only on the budget and the records pushed
/// so far, never on the thread count.
pub(crate) struct Chunk {
    pub arena: Vec<u8>,
    pub entries: Vec<Entry>,
    budget: usize,
}

impl Chunk {
    pub fn new(budget: usize) -> Chunk {
        Chunk {
            arena: Vec::new(),
            entries: Vec::new(),
            budget,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Bytes allocated by the chunk.
    fn held(&self) -> usize {
        self.arena
            .capacity()
            .saturating_add(self.entries.capacity().saturating_mul(ENTRY_BYTES))
    }

    pub fn raw(&self, e: &Entry) -> &[u8] {
        &self.arena[e.off..e.off + e.len as usize]
    }

    /// Append a record; returns false (appending nothing) when the chunk is
    /// full. An empty chunk accepts any record.
    #[inline]
    pub fn try_push(&mut self, chrom: u32, key: u64, end: u64, seq: u64, raw: &[u8]) -> bool {
        let fits = self.arena.capacity() - self.arena.len() >= raw.len()
            && self.entries.len() < self.entries.capacity();
        if !fits && !self.grow(raw.len()) {
            return false;
        }
        let off = self.arena.len();
        self.arena.extend_from_slice(raw);
        self.entries.push(Entry {
            key,
            end,
            seq,
            off,
            chrom,
            len: raw.len() as u32,
        });
        true
    }

    /// Make room for one more record of `len` bytes, or return false.
    #[cold]
    fn grow(&mut self, len: usize) -> bool {
        if !self.entries.is_empty() {
            return self.plan(len);
        }
        if len.saturating_add(ENTRY_BYTES) > self.budget {
            // A record larger than the budget gets a chunk of its own.
            self.arena = Vec::with_capacity(len);
            self.entries = Vec::with_capacity(1);
            return true;
        }
        if self.plan(len) {
            return true;
        }
        // Capacity kept from an earlier chunk is in the way: start afresh.
        self.release();
        let planned = self.plan(len);
        debug_assert!(planned, "an empty chunk accepts a record within the budget");
        planned
    }

    fn plan(&mut self, len: usize) -> bool {
        let b = self.budget;
        let (cap_a, cap_e) = (self.arena.capacity(), self.entries.capacity());
        let need_a = self.arena.len() + len;
        let need_e = self.entries.len() + 1;
        let min_a = need_a.max(cap_a);
        let min_e = need_e.max(cap_e);
        let min_total = min_a.saturating_add(min_e.saturating_mul(ENTRY_BYTES));
        if min_total > b {
            return false;
        }
        let held = self.held();
        // The smallest `b / 2^j` that holds the record and at least doubles
        // what is held now.
        let floor = min_total
            .max(FIRST_TARGET)
            .max(held.saturating_mul(2))
            .min(b);
        let mut target = b;
        while target / 2 >= floor {
            target /= 2;
        }
        // Past 1/16 of the budget, go straight to the whole budget: the
        // freed intermediate allocations would otherwise pile up in the
        // allocator (they are not returned to the OS promptly) and push
        // the peak RSS well past the budget. The untouched part of the
        // allocation costs no resident memory.
        if target > b / 16 {
            target = b;
        }
        // Split the target by the mix of record bytes to entries so far.
        let per_rec = need_a / need_e;
        let recs = target / (per_rec + ENTRY_BYTES);
        let want_a = (recs * per_rec).max(min_a);
        let want_e = recs.max(min_e);
        let (new_a, new_e) = if want_a.saturating_add(want_e * ENTRY_BYTES) <= target {
            (want_a, want_e)
        } else {
            // Keep the minimums and share the rest of the target in proportion.
            let room = (target - min_total) as u128;
            let extra_a = (want_a - min_a) as u128;
            let extra_e = ((want_e - min_e) * ENTRY_BYTES) as u128;
            let give_a = (room * extra_a / (extra_a + extra_e)) as usize;
            let give_e = (room as usize - give_a) / ENTRY_BYTES;
            (min_a + give_a, min_e + give_e)
        };
        // While a vector is copied into its new allocation, its old copy is
        // held too.
        let grow_a = new_a > cap_a;
        let grow_e = new_e > cap_e;
        if (grow_a && held.saturating_add(cap_a) > b)
            || (grow_e && held.saturating_add(cap_e * ENTRY_BYTES) > b)
        {
            return false;
        }
        if grow_a {
            self.arena.reserve_exact(new_a - self.arena.len());
        }
        if grow_e {
            self.entries.reserve_exact(new_e - self.entries.len());
        }
        true
    }

    /// Empty the chunk for the next one, keeping its capacity unless it
    /// exceeds the budget (after an oversized record) or the chunk used
    /// less than three quarters of it (the mix of record lengths changed),
    /// so that the next chunk plans afresh.
    pub fn clear(&mut self) {
        let used = self.arena.len() + self.entries.len() * ENTRY_BYTES;
        let held = self.held();
        self.arena.clear();
        self.entries.clear();
        if held > self.budget || used < held / 4 * 3 {
            self.release();
        }
    }

    /// Free all memory.
    pub fn release(&mut self) {
        self.arena = Vec::new();
        self.entries = Vec::new();
    }

    /// Sort by `(rank[chrom], key, end, tiebreak)`. The order is total
    /// (sequence numbers are unique), so an unstable or parallel sort gives
    /// the same result for any thread count.
    pub fn sort(&mut self, ranks: &[u32], tie: Tie, parallel: bool) {
        let arena = &self.arena;
        let entries = &mut self.entries;
        match tie {
            Tie::InputOrder => sort_by(entries, parallel, |a, b| {
                cmp_pos(ranks, a, b).then(a.seq.cmp(&b.seq))
            }),
            Tie::FullLine => sort_by(entries, parallel, |a, b| {
                cmp_pos(ranks, a, b)
                    .then_with(|| {
                        let ra = &arena[a.off..a.off + a.len as usize];
                        let rb = &arena[b.off..b.off + b.len as usize];
                        ra.cmp(rb)
                    })
                    .then(a.seq.cmp(&b.seq))
            }),
        }
    }
}

#[inline]
fn cmp_pos(ranks: &[u32], a: &Entry, b: &Entry) -> Ordering {
    if a.chrom != b.chrom {
        ranks[a.chrom as usize].cmp(&ranks[b.chrom as usize])
    } else {
        (a.key, a.end).cmp(&(b.key, b.end))
    }
}

fn sort_by<F>(v: &mut [Entry], parallel: bool, cmp: F)
where
    F: Fn(&Entry, &Entry) -> Ordering + Sync,
{
    if parallel {
        v.par_sort_unstable_by(cmp);
    } else {
        v.sort_unstable_by(cmp);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fill(c: &mut Chunk, len: usize) -> usize {
        let raw = vec![b'a'; len];
        let mut n = 0;
        while c.try_push(0, n as u64, 0, n as u64, &raw) {
            n += 1;
            let oversized = n == 1 && len + ENTRY_BYTES > c.budget;
            assert!(
                oversized || c.held() <= c.budget,
                "held {} > budget {}",
                c.held(),
                c.budget
            );
        }
        n
    }

    #[test]
    fn capacity_stays_within_the_budget_and_is_well_used() {
        for budget in [1 << 10, 100_000, 1 << 20, 3 << 20] {
            for len in [1, 30, 200, 5_000] {
                let mut c = Chunk::new(budget);
                let n = fill(&mut c, len);
                let used = n * (len + ENTRY_BYTES);
                assert!(n >= 1);
                assert!(
                    used * 10 >= budget * 9 || used + len + ENTRY_BYTES > budget,
                    "budget {budget}, len {len}: only {used} bytes used"
                );
                // The next chunk with the same mix reuses the capacity
                // (released instead after an oversized record).
                c.clear();
                let held = c.held();
                assert_eq!(fill(&mut c, len), n);
                if len + ENTRY_BYTES <= budget {
                    assert_eq!(c.held(), held);
                } else {
                    assert_eq!(held, 0);
                }
            }
        }
    }

    #[test]
    fn changed_mix_plans_afresh() {
        let mut c = Chunk::new(1 << 20);
        fill(&mut c, 10);
        c.clear();
        let n = fill(&mut c, 2_000);
        c.clear();
        assert_eq!(fill(&mut c, 2_000), (1 << 20) / (2_000 + ENTRY_BYTES));
        assert!(n > 0);
    }

    #[test]
    fn oversized_record_gets_its_own_chunk() {
        let mut c = Chunk::new(1 << 10);
        assert!(c.try_push(0, 0, 0, 0, &[b'a'; 10]));
        assert!(!c.try_push(0, 0, 0, 1, &[b'b'; 4 << 10]));
        c.clear();
        assert!(c.try_push(0, 0, 0, 1, &[b'b'; 4 << 10]));
        assert!(!c.try_push(0, 0, 0, 2, &[b'c'; 1]));
        c.clear();
        assert!(c.held() <= 1 << 10);
        assert!(c.try_push(0, 0, 0, 2, &[b'c'; 1]));
    }

    #[test]
    #[cfg(target_pointer_width = "64")]
    fn entry_is_40_bytes() {
        assert_eq!(ENTRY_BYTES, 40);
    }

    #[test]
    fn zero_budget_holds_one_record_per_chunk() {
        let mut c = Chunk::new(0);
        assert!(c.try_push(0, 0, 0, 0, b"x"));
        assert!(!c.try_push(0, 0, 0, 1, b"y"));
    }
}
