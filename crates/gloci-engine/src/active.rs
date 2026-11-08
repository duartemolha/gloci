//! The active set of one B input: the records loaded for the current A
//! that can still overlap it or a later A, with their raw bytes in an arena.

use crate::hit::Hit;
use gloci_core::{Interval, Strand};
use gloci_io::RecordMeta;

/// A loaded B record; its raw bytes are `arena[off..off + len]`.
#[derive(Debug, Clone, Copy)]
struct Entry {
    key: u64,
    span: Interval,
    loc: u64,
    off: usize,
    len: u32,
    strand: Strand,
}

/// The live entries are `entries[head..]`, in load order, which is B's
/// input order and therefore key order. Evicting moves live entries
/// towards the end and leaves free slots before `head`. Arena bytes are
/// in the same order as the entries.
#[derive(Debug, Default)]
pub(crate) struct ActiveSet {
    entries: Vec<Entry>,
    head: usize,
    arena: Vec<u8>,
    /// Bytes of the arena that belong to evicted records.
    dead: usize,
}

impl ActiveSet {
    pub(crate) fn clear(&mut self) {
        self.entries.clear();
        self.head = 0;
        self.arena.clear();
        self.dead = 0;
    }

    /// Reclaim what evicted records left behind: the free entry slots once
    /// they outnumber the live ones, and the arena once more than half of
    /// it is dead. Both are amortised O(1) per record. Call it between A
    /// records, when no hit borrows the arena.
    pub(crate) fn compact(&mut self) {
        if self.head == self.entries.len() {
            self.clear();
            return;
        }
        if self.head * 2 > self.entries.len() {
            self.entries.drain(..self.head);
            self.head = 0;
        }
        if self.dead * 2 > self.arena.len() {
            let mut w = 0;
            for e in &mut self.entries[self.head..] {
                let n = e.len as usize;
                self.arena.copy_within(e.off..e.off + n, w);
                e.off = w;
                w += n;
            }
            self.arena.truncate(w);
            self.dead = 0;
        }
    }

    /// Add a record after every loaded one (B's input order).
    pub(crate) fn push(&mut self, m: &RecordMeta, raw: &[u8]) {
        let off = self.arena.len();
        self.arena.extend_from_slice(raw);
        self.entries.push(Entry {
            key: m.key,
            span: m.span,
            loc: m.loc,
            off,
            len: m.len,
            strand: m.strand,
        });
    }

    /// Process the current A, whose key is `a_key` and whose padded end is
    /// `lim`: evict every entry with `span.end + pad < a_key` (it cannot
    /// overlap this or any later A), and push every entry whose span
    /// passes `is_hit` to `hits`, in input order. Only entries with
    /// `key <= lim` are looked at: the rest cannot overlap A, and since
    /// `end >= key > lim >= a_key` none of them can be evicted yet.
    /// Returns the number of live entries.
    pub(crate) fn scan<'h>(
        &'h mut self,
        file: u16,
        lim: u64,
        a_key: u64,
        pad: u64,
        is_hit: impl Fn(Interval) -> bool,
        hits: &mut Vec<Hit<'h>>,
    ) -> usize {
        let ActiveSet {
            entries,
            head,
            arena,
            dead,
        } = self;
        let arena: &'h [u8] = arena;
        let evicted = |e: &Entry| e.span.end.saturating_add(pad) < a_key;
        let mut last_evicted = None;
        for (i, e) in entries.iter().enumerate().skip(*head) {
            if e.key > lim {
                break;
            }
            if evicted(e) {
                *dead += e.len as usize;
                last_evicted = Some(i);
            } else if is_hit(e.span) {
                hits.push(Hit {
                    file,
                    span: e.span,
                    strand: e.strand,
                    raw: &arena[e.off..e.off + e.len as usize],
                    loc: e.loc,
                });
            }
        }
        if let Some(last) = last_evicted {
            // Move the live entries before the last evicted one up against
            // the entries after it, keeping their order. Entries evicted at
            // the front cost nothing.
            let mut w = last + 1;
            for i in (*head..last).rev() {
                if !evicted(&entries[i]) {
                    w -= 1;
                    entries[w] = entries[i];
                }
            }
            *head = w;
        }
        entries.len() - *head
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gloci_core::ChromId;

    fn meta(key: u64, end: u64, len: u32) -> RecordMeta {
        RecordMeta {
            chrom: ChromId(0),
            key,
            span: Interval { start: key, end },
            strand: Strand::Unknown,
            loc: key,
            off: 0,
            len,
        }
    }

    /// The raw lines of the hits for an A at `[a_key, lim)` that every
    /// live entry overlaps, i.e. the live entries up to `lim`.
    fn live(set: &mut ActiveSet, a_key: u64, lim: u64) -> Vec<Vec<u8>> {
        let mut hits = Vec::new();
        set.scan(0, lim, a_key, 0, |_| true, &mut hits);
        hits.iter().map(|h| h.raw.to_vec()).collect()
    }

    #[test]
    fn eviction_keeps_input_order_and_compaction_keeps_bytes() {
        let mut set = ActiveSet::default();
        // Records `k` cover [k, ends[k]); raw bytes are their key repeated.
        let ends = [100, 3, 50, 4, 5, 60, 6, 7];
        for (k, &e) in ends.iter().enumerate() {
            let raw = vec![b'a' + k as u8; k + 1];
            set.push(&meta(k as u64, e, raw.len() as u32), &raw);
        }
        // An A at key 8 evicts the records ending before 8 (keys 1, 3, 4, 6, 7).
        let got = live(&mut set, 8, 8);
        assert_eq!(got, [b"a".to_vec(), b"ccc".to_vec(), b"ffffff".to_vec()]);
        assert_eq!(set.head, 5);
        assert!(set.dead * 2 > set.arena.len());
        set.compact();
        assert_eq!(set.head, 0);
        assert_eq!(set.dead, 0);
        assert_eq!(set.arena, b"acccffffff");
        assert_eq!(live(&mut set, 8, 8), got);
        // Entries past `lim` are neither reported nor evicted.
        set.push(&meta(200, 300, 1), b"z");
        assert_eq!(live(&mut set, 70, 150), [b"a".to_vec()]);
        assert_eq!(live(&mut set, 70, 250), [b"a".to_vec(), b"z".to_vec()]);
    }
}
