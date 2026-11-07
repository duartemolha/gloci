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

/// Rank of each name in karyotype order (spec §4.5 rule 3), written out
/// independently of gloci-core: sort by a key tuple, then number the
/// distinct names. Names not in `names` rank after all of them.
pub fn karyotype_rank(names: &[&str]) -> impl Fn(&str) -> u32 {
    fn key(name: &str) -> (u8, Vec<(u8, Vec<u8>)>, String) {
        let r = if name.len() > 3 && name[..3].eq_ignore_ascii_case("chr") {
            &name[3..]
        } else {
            name
        };
        // Digit runs as (0, zero-stripped digits padded by length) so that
        // byte order equals numeric order; other runs as (1, bytes).
        let runs = |s: &str| -> Vec<(u8, Vec<u8>)> {
            let mut out: Vec<(u8, Vec<u8>)> = Vec::new();
            let mut last_digit = None;
            for c in s.bytes() {
                let d = c.is_ascii_digit();
                if last_digit == Some(d) {
                    out.last_mut().unwrap().1.push(c);
                } else {
                    out.push((if d { 0 } else { 1 }, vec![c]));
                    last_digit = Some(d);
                }
            }
            for (k, v) in &mut out {
                if *k == 0 {
                    let t: Vec<u8> = v.iter().copied().skip_while(|&c| c == b'0').collect();
                    let mut padded = (t.len() as u64).to_be_bytes().to_vec();
                    padded.extend(t);
                    *v = padded;
                }
            }
            out
        };
        let digits = r.bytes().take_while(u8::is_ascii_digit).count();
        let numbered = digits > 0 && r[digits..].bytes().all(|c| c.is_ascii_alphabetic());
        let (class, sub) = if numbered {
            (0, runs(r))
        } else {
            match r {
                "X" => (1, vec![]),
                "Y" => (2, vec![]),
                "M" | "MT" => (3, vec![]),
                _ => (4, runs(r)),
            }
        };
        (class, sub, name.to_string())
    }
    let mut sorted: Vec<String> = names.iter().map(|s| s.to_string()).collect();
    sorted.sort_by_key(|n| key(n));
    sorted.dedup();
    move |n: &str| sorted.iter().position(|s| s == n).unwrap_or(sorted.len()) as u32
}
