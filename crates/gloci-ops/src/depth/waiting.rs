//! The table behind both mate rules: reads waiting for their mate, each with a value, dropped
//! once no later read can be the mate.

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};

/// Name to (start of the mate, value). A min-heap on the mate's start drops entries nobody can
/// claim any more; an entry taken or replaced leaves a stale heap item, which is skipped when
/// it surfaces.
pub(super) struct Waiting<V> {
    pending: HashMap<Box<[u8]>, (u64, V)>,
    by_mpos: BinaryHeap<Reverse<(u64, Box<[u8]>)>>,
    peak: usize,
}

impl<V> Default for Waiting<V> {
    fn default() -> Self {
        Waiting {
            pending: HashMap::new(),
            by_mpos: BinaryHeap::new(),
            peak: 0,
        }
    }
}

impl<V> Waiting<V> {
    pub(super) fn contains(&self, name: &[u8]) -> bool {
        self.pending.contains_key(name)
    }

    /// Wait for the mate of `name` at `mpos`, replacing an entry of that name.
    pub(super) fn hold(&mut self, name: &[u8], mpos: u64, value: V) {
        let name: Box<[u8]> = name.into();
        self.by_mpos.push(Reverse((mpos, name.clone())));
        self.pending.insert(name, (mpos, value));
        self.peak = self.peak.max(self.pending.len());
    }

    /// Remove the entry of `name` and return its value.
    pub(super) fn take(&mut self, name: &[u8]) -> Option<V> {
        self.pending.remove(name).map(|(_, v)| v)
    }

    /// Drop the entries whose mate position is below `pos`.
    pub(super) fn advance(&mut self, pos: u64) {
        while let Some(Reverse((mpos, _))) = self.by_mpos.peek() {
            if *mpos >= pos {
                break;
            }
            let Some(Reverse((mpos, name))) = self.by_mpos.pop() else {
                break;
            };
            // A later read of the same name may have been held with its own mate position.
            if self.pending.get(&name).is_some_and(|&(m, _)| m == mpos) {
                self.pending.remove(&name);
            }
        }
        if self.pending.is_empty() {
            self.by_mpos.clear();
        }
    }

    pub(super) fn reset(&mut self) {
        self.pending.clear();
        self.by_mpos.clear();
    }

    pub(super) fn len(&self) -> usize {
        self.pending.len()
    }

    pub(super) fn is_empty(&self) -> bool {
        self.pending.is_empty()
    }

    pub(super) fn peak(&self) -> usize {
        self.peak
    }
}
