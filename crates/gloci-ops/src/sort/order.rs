//! Sort orders: the interval order, and samtools 1.21's alignment orders (D121).
//!
//! Alignment records are BAM bytes without `block_size`: `refID` (i32) at 0, `pos` (i32) at 4,
//! `l_read_name` (u8) at 8, the flag (u16) at 14, and the NUL-ended read name from byte 32.

use super::Tie;
use gloci_io::AlnSortOrder;
use std::cmp::Ordering;

/// How the sorter orders records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Order {
    /// `(chromosome rank, key, end, tiebreak)` for intervals.
    Interval(Tie),
    /// `samtools sort`'s coordinate, `-n` and `-N` orders for BAM records.
    Aln(AlnSortOrder),
}

/// The compile-time mode of a comparison: how heads and entries compare.
pub(crate) const M_INTERVAL: u8 = 0;
pub(crate) const M_FULL: u8 = 1;
pub(crate) const M_COORD: u8 = 2;
pub(crate) const M_NAME: u8 = 3;
pub(crate) const M_NAME_BYTES: u8 = 4;

impl Order {
    pub(crate) fn mode(self) -> u8 {
        match self {
            Order::Interval(Tie::InputOrder) => M_INTERVAL,
            Order::Interval(Tie::FullLine) => M_FULL,
            Order::Aln(AlnSortOrder::Coordinate) => M_COORD,
            Order::Aln(AlnSortOrder::Name) => M_NAME,
            Order::Aln(AlnSortOrder::NameBytes) => M_NAME_BYTES,
        }
    }
}

/// Whether a mode compares whole raw records (so merges keep their bytes).
pub(crate) const fn keeps_bytes(mode: u8) -> bool {
    mode != M_INTERVAL && mode != M_COORD
}

#[inline]
fn is_digit(c: u8) -> bool {
    c.is_ascii_digit()
}

/// samtools 1.21 `strnum_cmp` (natural order), ported line by line from `bam_sort.c`: unsigned
/// bytes, digit runs compared as numbers after their leading zeros (so `r01` equals `r1`).
/// The C strings end at a NUL, which is read here as the end of the slice.
pub fn strnum_cmp(a: &[u8], b: &[u8]) -> Ordering {
    let at = |s: &[u8], i: usize| s.get(i).copied().unwrap_or(0);
    let (mut pa, mut pb) = (0usize, 0usize);
    while at(a, pa) != 0 && at(b, pb) != 0 {
        if !is_digit(at(a, pa)) || !is_digit(at(b, pb)) {
            if at(a, pa) != at(b, pb) {
                return at(a, pa).cmp(&at(b, pb));
            }
            pa += 1;
            pb += 1;
        } else {
            // skip leading zeros
            while at(a, pa) == b'0' {
                pa += 1;
            }
            while at(b, pb) == b'0' {
                pb += 1;
            }
            // skip matching digits
            while is_digit(at(a, pa)) && at(a, pa) == at(b, pb) {
                pa += 1;
                pb += 1;
            }
            // Now mismatching, so see which ends the number sooner
            let diff = i32::from(at(a, pa)) - i32::from(at(b, pb));
            while is_digit(at(a, pa)) && is_digit(at(b, pb)) {
                pa += 1;
                pb += 1;
            }
            if is_digit(at(a, pa)) {
                return Ordering::Greater; // pa still going, so larger
            } else if is_digit(at(b, pb)) {
                return Ordering::Less; // pb still going, so larger
            } else if diff != 0 {
                return diff.cmp(&0); // same length, so earlier diff
            }
        }
    }
    if at(a, pa) != 0 {
        Ordering::Greater
    } else if at(b, pb) != 0 {
        Ordering::Less
    } else {
        Ordering::Equal
    }
}

/// The rank of a flag among records of one name (`-n`, `-N`): neither mate, then READ1, then
/// READ2; within each, primary, then supplementary, then secondary.
pub fn name_flag_rank(flag: u16) -> u32 {
    let f = u32::from(flag);
    ((f & 0xC0) << 8) | ((f & 0x100) << 3) | ((f & 0x800) >> 3)
}

/// The read name of a BAM record (without its NUL), empty when the record is too short.
pub fn qname_of(raw: &[u8]) -> &[u8] {
    let n = usize::from(raw.get(8).copied().unwrap_or(0));
    let name = raw.get(32..32 + n).unwrap_or(&[]);
    // htslib reads the name as a C string.
    let end = name.iter().position(|&b| b == 0).unwrap_or(name.len());
    &name[..end]
}

pub(crate) fn flag_of(raw: &[u8]) -> u16 {
    raw.get(14..16)
        .map_or(0, |b| u16::from_le_bytes([b[0], b[1]]))
}

/// The coordinate sort key of a BAM record: `(tid, (pos + 1) << 1 | reverse)`, with an unplaced
/// reference id (-1) as `n_ref`, so unplaced reads come last.
pub fn coordinate_key(raw: &[u8], n_ref: u32) -> (u32, u64) {
    let i32_at = |i: usize| {
        raw.get(i..i + 4)
            .map_or(-1, |b| i32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    };
    let tid = match i32_at(0) {
        t if t < 0 => n_ref,
        t => t as u32,
    };
    let pos = u64::from(i32_at(4).wrapping_add(1) as u32);
    let rev = u64::from(flag_of(raw) & 0x10 != 0);
    (tid, (pos << 1) | rev)
}

/// Compares two BAM records by read name and flag rank (`-n` with `natural`, `-N` without).
pub fn name_cmp_raw(a: &[u8], b: &[u8], natural: bool) -> Ordering {
    let (qa, qb) = (qname_of(a), qname_of(b));
    let by_name = if natural {
        strnum_cmp(qa, qb)
    } else {
        // strcmp: unsigned bytes, the end of the string smallest.
        qa.cmp(qb)
    };
    by_name.then_with(|| name_flag_rank(flag_of(a)).cmp(&name_flag_rank(flag_of(b))))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cmp::Ordering::{Equal, Greater, Less};

    #[test]
    fn strnum_cmp_matches_samtools() {
        let c = |a: &str, b: &str| strnum_cmp(a.as_bytes(), b.as_bytes());
        assert_eq!(c("r2", "r10"), Less);
        assert_eq!(c("a1b2", "a1b10"), Less);
        assert_eq!(c("r01", "r1"), Equal);
        assert_eq!(c("007", "7"), Equal);
        assert_eq!(c("r1", "r1a"), Less);
        assert_eq!(c("R1", "r1"), Less);
        assert_eq!(c("a10", "a9b"), Greater);
        assert_eq!(c("x", "x0"), Less);
        // bytes (-N): plain strcmp
        assert_eq!("r10".as_bytes().cmp("r2".as_bytes()), Less);
    }

    #[test]
    fn name_flag_rank_table() {
        let order = [0u16, 0x800, 0x100, 0x40, 0x40 | 0x800, 0x40 | 0x100, 0x80];
        let ranks: Vec<u32> = order.iter().map(|&f| name_flag_rank(f)).collect();
        assert!(ranks.windows(2).all(|w| w[0] < w[1]), "{ranks:?}");
    }

    #[test]
    fn coordinate_key_values() {
        let mut raw = vec![0u8; 40];
        raw[0..4].copy_from_slice(&(-1i32).to_le_bytes());
        raw[4..8].copy_from_slice(&(-1i32).to_le_bytes());
        assert_eq!(coordinate_key(&raw, 7), (7, 0));
        raw[0..4].copy_from_slice(&2i32.to_le_bytes());
        raw[4..8].copy_from_slice(&99i32.to_le_bytes());
        raw[14..16].copy_from_slice(&0x10u16.to_le_bytes());
        assert_eq!(coordinate_key(&raw, 7), (2, (100 << 1) | 1));
        assert_eq!(coordinate_key(&raw[..3], 7).0, 7);
    }
}
