//! `--skip-mate-overlap` (D139): the table of reads whose mate is still ahead, so the part of
//! the second mate that overlaps the first is counted once.

use super::waiting::Waiting;
use gloci_io::AlnRef;

const FLAG_PAIRED: u16 = 0x1;
const FLAG_MATE_UNMAPPED: u16 = 0x8;

/// Reads waiting for their mate: name to (end of the read, start of the mate), both 0-based and
/// absolute, in a [`Waiting`] table.
#[derive(Default)]
pub struct MateTable {
    waiting: Waiting<u64>,
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
        if let Some(mate_end) = self.waiting.take(r.name()) {
            return mate_end;
        }
        if let (Some(mref), Some(mpos)) = (r.mate_ref_id(), r.mate_pos())
            && Some(mref) == r.ref_id()
            && pos <= mpos
            && mpos <= end
        {
            self.waiting.hold(r.name(), mpos, end);
        }
        0
    }

    /// Drop the entries whose mate would have started before `pos`: a read starting there
    /// cannot be it.
    pub fn advance(&mut self, pos: u64) {
        self.waiting.advance(pos);
    }

    /// Forget everything: a new reference starts.
    pub fn reset(&mut self) {
        self.waiting.reset();
    }

    /// Reads waiting for a mate.
    pub fn len(&self) -> usize {
        self.waiting.len()
    }

    pub fn is_empty(&self) -> bool {
        self.waiting.is_empty()
    }

    /// The most reads that waited at once.
    pub fn peak(&self) -> usize {
        self.waiting.peak()
    }
}
