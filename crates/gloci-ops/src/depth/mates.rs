//! `--skip-mate-overlap` (D139): the table of reads whose mate is still ahead, so the part of
//! the second mate that overlaps the first is counted once.

use gloci_io::AlnRef;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};

const FLAG_PAIRED: u16 = 0x1;
const FLAG_MATE_UNMAPPED: u16 = 0x8;

/// Reads waiting for their mate: name to (end of the read, start of the mate), both 0-based and
/// absolute. A min-heap on the mate's start drops entries nobody can claim any more; an entry
/// the mate claimed first leaves a stale heap item, which is skipped when it surfaces.
#[derive(Default)]
pub struct MateTable {
    pending: HashMap<Box<[u8]>, (u64, u64)>,
    by_mpos: BinaryHeap<Reverse<(u64, Box<[u8]>)>>,
    peak: usize,
}

impl MateTable {
    pub fn new() -> MateTable {
        MateTable::default()
    }

    /// The position below which the read `r` (starting at `pos`, ending at `end`) is not
    /// counted: the end of its mate if that is waiting, else 0. A read whose mate is still
    /// ahead on this reference (`pos <= mpos <= end`) is stored for it.
    pub fn clip_for(&mut self, r: &AlnRef<'_>, pos: u64, end: u64) -> u64 {
        let flag = r.flag();
        if flag & FLAG_PAIRED == 0 || flag & FLAG_MATE_UNMAPPED != 0 {
            return 0;
        }
        if let Some((mate_end, _)) = self.pending.remove(r.name()) {
            return mate_end;
        }
        if let (Some(mref), Some(mpos)) = (r.mate_ref_id(), r.mate_pos())
            && Some(mref) == r.ref_id()
            && pos <= mpos
            && mpos <= end
        {
            let name: Box<[u8]> = r.name().into();
            self.by_mpos.push(Reverse((mpos, name.clone())));
            self.pending.insert(name, (end, mpos));
            self.peak = self.peak.max(self.pending.len());
        }
        0
    }

    /// Drop the entries whose mate would have started before `pos`: a read starting there
    /// cannot be it.
    pub fn advance(&mut self, pos: u64) {
        while let Some(Reverse((mpos, _))) = self.by_mpos.peek() {
            if *mpos >= pos {
                break;
            }
            let Some(Reverse((mpos, name))) = self.by_mpos.pop() else {
                break;
            };
            // A later read of the same name may have been stored with its own mate position.
            if self.pending.get(&name).is_some_and(|&(_, m)| m == mpos) {
                self.pending.remove(&name);
            }
        }
        if self.pending.is_empty() {
            self.by_mpos.clear();
        }
    }

    /// Forget everything: a new reference starts.
    pub fn reset(&mut self) {
        self.pending.clear();
        self.by_mpos.clear();
    }

    /// Reads waiting for a mate.
    pub fn len(&self) -> usize {
        self.pending.len()
    }

    pub fn is_empty(&self) -> bool {
        self.pending.is_empty()
    }

    /// The most reads that waited at once.
    pub fn peak(&self) -> usize {
        self.peak
    }
}
