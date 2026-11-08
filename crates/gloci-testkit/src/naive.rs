//! Naive reference implementations to test against.
//!
//! The §4.3 overlap rule is written out here independently of gloci-core,
//! so the oracle does not share code with the implementation under test.
//! The relation references (`naive_hits`, `naive_neighbours`) are the
//! exception: they apply gloci-core's `overlaps`, `distance` and `Window`
//! in a plain double loop, so their property tests check an engine's
//! bookkeeping (loading, eviction, look-ahead, chromosome changes, hit
//! order) rather than a second copy of §4.3 and §4.4, which gloci-core's
//! own table tests pin.

use crate::strategies::{Feature, TestRec};
use gloci_core::{Interval, Strand, StrandRule, Window, distance, overlaps};

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

/// Candidate hits of a relation (D36, D38), by brute force: for each A in
/// the order given, every B of `b[f]` on A's chromosome whose span, padded
/// by `window` relative to A's strand, overlaps A's span under §4.3. Each
/// hit is `(f, B's line)`, ordered by file, then by B's order in `b[f]`.
pub fn naive_hits<R: Feature>(
    a: &[R],
    b: &[Vec<R>],
    window: Option<Window>,
    inclusive_points: bool,
) -> Vec<Vec<(u16, String)>> {
    a.iter()
        .map(|ra| {
            let mut hits = Vec::new();
            for (f, file) in b.iter().enumerate() {
                for rb in file {
                    if rb.chrom() != ra.chrom() {
                        continue;
                    }
                    let span = match window {
                        Some(w) => w.apply(rb.span(), ra.strand()),
                        None => rb.span(),
                    };
                    if overlaps(ra.span(), span, inclusive_points) {
                        hits.push((f as u16, rb.line()));
                    }
                }
            }
            hits
        })
        .collect()
}

/// One file's neighbours of one A: `(left, right)`, as indices into the
/// file's records, nearest first.
pub type Sides = (Vec<usize>, Vec<usize>);

/// Nearest-neighbour candidates (Task 6), by brute force: for each A in
/// the order given and each `b[f]`, the B on A's chromosome that do not
/// overlap A under §4.3, split into **left** (`end <= A.start`) and
/// **right** (`start >= A.end`). Each side is sorted by gloci-core
/// `distance`, ties in B's order in `b[f]`, and cut to the `k` nearest
/// plus every B tied with the k-th. Only B at distance at most
/// `max_distance` count. Under `Same`/`Opposite`, only B whose strand meets
/// the rule against A's count: a B or A with unknown strand never does.
/// With `split_by_b_strand` and `Any`, the cut is made separately for `-`
/// B and for the rest, and each side is the union, sorted the same way.
pub fn naive_neighbours<R: Feature>(
    a: &[R],
    b: &[Vec<R>],
    k: usize,
    max_distance: Option<u64>,
    strand: StrandRule,
    split_by_b_strand: bool,
    inclusive_points: bool,
) -> Vec<Vec<Sides>> {
    let strand_ok = |x: Strand, y: Strand| match strand {
        StrandRule::Any => true,
        rule => {
            x != Strand::Unknown && y != Strand::Unknown && ((x == y) == (rule == StrandRule::Same))
        }
    };
    // One group, or `-` B (`Some(true)`) apart from the rest.
    let groups: &[Option<bool>] = if split_by_b_strand && strand == StrandRule::Any {
        &[Some(false), Some(true)]
    } else {
        &[None]
    };
    // `v` holds (distance, index, whether B is on the `-` strand).
    let nearest = |v: Vec<(u64, usize, bool)>| -> Vec<usize> {
        let mut out = Vec::new();
        for &g in groups {
            let mut group: Vec<(u64, usize)> = v
                .iter()
                .filter(|x| g.is_none_or(|minus| x.2 == minus))
                .map(|x| (x.0, x.1))
                .collect();
            group.sort();
            if group.len() > k {
                let kth = group[k - 1].0;
                group.retain(|&(d, _)| d <= kth);
            }
            out.extend(group);
        }
        out.sort();
        out.into_iter().map(|(_, i)| i).collect()
    };
    a.iter()
        .map(|ra| {
            let x = ra.span();
            b.iter()
                .map(|file| {
                    let (mut left, mut right) = (Vec::new(), Vec::new());
                    for (i, rb) in file.iter().enumerate() {
                        let y = rb.span();
                        if rb.chrom() != ra.chrom()
                            || !strand_ok(ra.strand(), rb.strand())
                            || overlaps(x, y, inclusive_points)
                        {
                            continue;
                        }
                        let d = distance(x, y, inclusive_points);
                        if max_distance.is_some_and(|m| d > m) {
                            continue;
                        }
                        let minus = rb.strand() == Strand::Minus;
                        if y.end <= x.start {
                            left.push((d, i, minus));
                        } else if y.start >= x.end {
                            right.push((d, i, minus));
                        }
                    }
                    (nearest(left), nearest(right))
                })
                .collect()
        })
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
