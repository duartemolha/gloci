//! Naive reference implementations to test against.
//!
//! The §4.3 overlap rule is written out here independently of gloci-core,
//! so the oracle does not share code with the implementation under test.

use crate::strategies::TestRec;
use gloci_core::Interval;

fn inside(p: u64, start: u64, end: u64, inclusive: bool) -> bool {
    if inclusive {
        start <= p && p <= end
    } else {
        start < p && p < end
    }
}

/// Does `rec` overlap the region `chrom:iv` under §4.3? A zero-length
/// interval is a point: it must lie strictly inside a non-empty interval
/// (or on an edge with `inclusive`), and two points must be equal.
pub fn matches(rec: &TestRec, chrom: &[u8], iv: Interval, inclusive: bool) -> bool {
    if rec.chrom.as_bytes() != chrom {
        return false;
    }
    let (s, e) = (rec.start, rec.end);
    match (s == e, iv.start == iv.end) {
        (true, true) => s == iv.start,
        (true, false) => inside(s, iv.start, iv.end, inclusive),
        (false, true) => inside(iv.start, s, e, inclusive),
        (false, false) => s < iv.end && iv.start < e,
    }
}

/// The records of `recs` that overlap `chrom:iv`, in input order.
pub fn region_query(recs: &[TestRec], chrom: &[u8], iv: Interval, inclusive: bool) -> Vec<TestRec> {
    recs.iter()
        .filter(|r| matches(r, chrom, iv, inclusive))
        .cloned()
        .collect()
}

/// Stable sort by `(chrom_rank(chrom), start, end)`.
pub fn sort(mut recs: Vec<TestRec>, chrom_rank: impl Fn(&str) -> u32) -> Vec<TestRec> {
    recs.sort_by_key(|r| (chrom_rank(&r.chrom), r.start, r.end));
    recs
}

/// BED lines sorted by `(chrom_rank(chrom), start, end, whole line)` (the
/// `--full-key` order).
pub fn sort_full_line(mut recs: Vec<TestRec>, chrom_rank: impl Fn(&str) -> u32) -> Vec<String> {
    recs.sort_by_cached_key(|r| (chrom_rank(&r.chrom), r.start, r.end, r.bed_line()));
    recs.iter().map(TestRec::bed_line).collect()
}

/// `uniq`: one copy of each run of identical adjacent lines.
pub fn uniq(lines: &[String]) -> Vec<String> {
    let mut out = lines.to_vec();
    out.dedup();
    out
}

/// `uniq -d`: one copy of each run of identical adjacent lines that has at
/// least two members.
pub fn uniq_d(lines: &[String]) -> Vec<String> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < lines.len() {
        let mut j = i + 1;
        while j < lines.len() && lines[j] == lines[i] {
            j += 1;
        }
        if j - i >= 2 {
            out.push(lines[i].clone());
        }
        i = j;
    }
    out
}
