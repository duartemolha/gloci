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
use gloci_core::{Interval, OverlapCriteria, Strand, StrandRule, Window, distance, overlaps};

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

fn spans_overlap(a: Interval, b: Interval, inclusive: bool) -> bool {
    match (a.start == a.end, b.start == b.end) {
        (true, true) => a.start == b.start,
        (true, false) => inside(a.start, b.start, b.end, inclusive),
        (false, true) => inside(b.start, a.start, a.end, inclusive),
        (false, false) => a.start < b.end && b.start < a.end,
    }
}

fn shared_bases(a: Interval, b: Interval) -> u64 {
    a.end.min(b.end).saturating_sub(a.start.max(b.start))
}

fn strand_met(rule: StrandRule, a: Strand, b: Strand) -> bool {
    match rule {
        StrandRule::Any => true,
        _ if a == Strand::Unknown || b == Strand::Unknown => false,
        StrandRule::Same => a == b,
        StrandRule::Opposite => a != b,
    }
}

/// Spec 4.3 written out with plain arithmetic: overlapping spans, the
/// thresholds (each met when either side is a point; `-r` with one
/// fraction applies it to both sides; `-e` needs either of the given
/// fractions) and the strand rule. Only `Fraction::met` (exact decimal
/// comparison) comes from gloci-core.
fn naive_accepts(c: &OverlapCriteria, a: Interval, sa: Strand, b: Interval, sb: Strand) -> bool {
    if !spans_overlap(a, b, c.inclusive_points) || !strand_met(c.strand, sa, sb) {
        return false;
    }
    let shared = shared_bases(a, b);
    let point = a.start == a.end || b.start == b.end;
    if point {
        return true;
    }
    if c.min_bases.is_some_and(|n| shared < n) {
        return false;
    }
    let (mut fa, mut fb) = (c.min_frac_a, c.min_frac_b);
    if c.reciprocal {
        let one = fa.or(fb);
        fa = fa.or(one);
        fb = fb.or(one);
    }
    let ok_a = fa.map(|f| f.met(shared, a.end - a.start));
    let ok_b = fb.map(|f| f.met(shared, b.end - b.start));
    match (ok_a, ok_b) {
        (None, None) => true,
        (Some(x), None) | (None, Some(x)) => x,
        (Some(x), Some(y)) if c.either => x || y,
        (Some(x), Some(y)) => x && y,
    }
}

/// What [`naive_intersect`] reports (the `intersect --mode` choices that
/// write one line per record or pair).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NaiveMode {
    Any,
    None,
    Pairs {
        overlap_bases: bool,
        keep_unmatched: bool,
    },
    Count,
}

/// `intersect`, by brute force, as the lines gloci writes for BED inputs.
/// For each A in the order given: the B that are candidate hits (as
/// [`naive_hits`]) and pass `criteria` against B's own span (with a window,
/// only the strand rule applies, D37). `labels` is `Some` with several B
/// files: pairs rows then carry the file's label, and unmatched rows `.`.
/// An unmatched row's B columns follow D45 for BED6 B (the test records),
/// or 3 columns when the first B file is empty.
pub fn naive_intersect<R: Feature>(
    a: &[R],
    b: &[Vec<R>],
    criteria: &OverlapCriteria,
    window: Option<Window>,
    mode: NaiveMode,
    labels: Option<&[String]>,
) -> Vec<String> {
    let mut out = Vec::new();
    let placeholder = if b.first().is_some_and(|f| !f.is_empty()) {
        ".\t-1\t-1\t.\t-1\t."
    } else {
        ".\t-1\t-1"
    };
    for ra in a {
        let accepted = naive_accepted(ra, b, criteria, window);
        match mode {
            NaiveMode::Any if !accepted.is_empty() => out.push(ra.line()),
            NaiveMode::None if accepted.is_empty() => out.push(ra.line()),
            NaiveMode::Count => out.push(format!("{}\t{}", ra.line(), accepted.len())),
            NaiveMode::Pairs {
                overlap_bases,
                keep_unmatched,
            } => {
                let tail = |ov: u64| {
                    if overlap_bases {
                        format!("\t{ov}")
                    } else {
                        String::new()
                    }
                };
                for (f, rb) in &accepted {
                    let label = labels.map_or(String::new(), |l| format!("\t{}", l[*f]));
                    out.push(format!(
                        "{}{}\t{}{}",
                        ra.line(),
                        label,
                        rb.line(),
                        tail(shared_bases(ra.span(), rb.span()))
                    ));
                }
                if accepted.is_empty() && keep_unmatched {
                    let label = labels.map_or(String::new(), |_| "\t.".to_string());
                    out.push(format!(
                        "{}{}\t{}{}",
                        ra.line(),
                        label,
                        placeholder,
                        tail(0)
                    ));
                }
            }
            _ => {}
        }
    }
    out
}

/// The B records that count for `ra` (D36, D37): candidate hits (spans
/// overlapping after the window) that pass `criteria` against B's own span,
/// or only the strand rule when a window is given. As `(file, record)`, by
/// file and then input order.
fn naive_accepted<'b, R: Feature>(
    ra: &R,
    b: &'b [Vec<R>],
    criteria: &OverlapCriteria,
    window: Option<Window>,
) -> Vec<(usize, &'b R)> {
    let mut accepted = Vec::new();
    for (f, file) in b.iter().enumerate() {
        for rb in file {
            if rb.chrom() != ra.chrom() {
                continue;
            }
            let padded = match window {
                Some(w) => w.apply(rb.span(), ra.strand()),
                None => rb.span(),
            };
            if !spans_overlap(ra.span(), padded, criteria.inclusive_points) {
                continue;
            }
            let ok = match window {
                Some(_) => strand_met(criteria.strand, ra.strand(), rb.strand()),
                None => naive_accepts(criteria, ra.span(), ra.strand(), rb.span(), rb.strand()),
            };
            if ok {
                accepted.push((f, rb));
            }
        }
    }
    accepted
}

/// `subtract`, by brute force, as the lines gloci writes for BED inputs.
///
/// - **Default.** A point A is dropped when it has an accepted hit.
///   Otherwise A's bases minus the bases of every accepted hit's windowed
///   span, found by cutting A at every removal edge; each maximal run left is written as
///   A's line with new coordinates. A with nothing removed is written as is.
/// - **`whole`.** A is dropped when it has an accepted hit.
/// - **`sum_overlaps`.** Only the strand rule (and the window) decide which
///   hits take part; A is dropped when the shared bases of those hits sum
///   to at least `criteria.min_frac_a` of A (any such hit when there is no
///   fraction, or A is a point). Other thresholds in `criteria` are ignored.
pub fn naive_subtract(
    a: &[TestRec],
    b: &[Vec<TestRec>],
    criteria: &OverlapCriteria,
    window: Option<Window>,
    whole: bool,
    sum_overlaps: bool,
) -> Vec<String> {
    let strand_only = OverlapCriteria {
        min_frac_a: None,
        min_frac_b: None,
        reciprocal: false,
        either: false,
        min_bases: None,
        ..*criteria
    };
    let used = if sum_overlaps { &strand_only } else { criteria };
    let mut out = Vec::new();
    for ra in a {
        let accepted = naive_accepted(ra, b, used, window);
        let span = ra.span();
        let windowed = |rb: &TestRec| match window {
            Some(w) => w.apply(rb.span(), ra.strand()),
            None => rb.span(),
        };
        if sum_overlaps {
            let sum: u64 = accepted
                .iter()
                .map(|(_, rb)| shared_bases(span, windowed(rb)))
                .sum();
            let drop = match criteria.min_frac_a {
                Some(f) if span.start != span.end => f.met(sum, span.end - span.start),
                _ => !accepted.is_empty(),
            };
            if !drop {
                out.push(ra.line());
            }
            continue;
        }
        if whole || span.start == span.end {
            if accepted.is_empty() {
                out.push(ra.line());
            }
            continue;
        }
        // Cut A at every edge of a removal that falls inside it; a segment
        // between two cuts is removed when some removal spans all of it.
        let wins: Vec<Interval> = accepted.iter().map(|(_, rb)| windowed(rb)).collect();
        let mut cuts = vec![span.start, span.end];
        for w in &wins {
            cuts.extend(
                [w.start, w.end]
                    .into_iter()
                    .filter(|&c| span.start < c && c < span.end),
            );
        }
        cuts.sort_unstable();
        cuts.dedup();
        let kept: Vec<(u64, u64)> = cuts
            .windows(2)
            .map(|c| (c[0], c[1]))
            .filter(|&(from, to)| !wins.iter().any(|w| w.start <= from && to <= w.end))
            .collect();
        if kept == [(span.start, span.end)] {
            out.push(ra.line());
            continue;
        }
        // Join neighbouring kept segments into maximal runs.
        let mut runs: Vec<(u64, u64)> = Vec::new();
        for (from, to) in kept {
            match runs.last_mut() {
                Some(last) if last.1 == from => last.1 = to,
                _ => runs.push((from, to)),
            }
        }
        for (start, end) in runs {
            let piece = TestRec {
                start,
                end,
                ..ra.clone()
            };
            out.push(piece.bed_line());
        }
    }
    out
}

/// `(rank, start, end, strand slot, chrom, count, names)`.
type MergeRow = (u32, u64, u64, usize, String, usize, Vec<String>);

/// `merge` by brute force over BED records: all inputs are concatenated in
/// input order and stably sorted by `(chromosome rank, start)`, which is the
/// order a k-way merge yields; then clusters are swept per chromosome (and
/// per strand with `per_strand`). A record joins the open cluster when
/// `start <= cluster_end + max_gap`. Rows are `chrom start end count names`,
/// where `names` is the comma-joined column 4 in that order, sorted by
/// `(chromosome rank, start, end, strand)` with strands `+`, `-`, `.`.
pub fn naive_merge(
    inputs: &[Vec<TestRec>],
    chrom_rank: impl Fn(&str) -> u32,
    max_gap: i64,
    per_strand: bool,
) -> Vec<String> {
    let mut all: Vec<&TestRec> = inputs.iter().flatten().collect();
    all.sort_by_key(|r| (chrom_rank(&r.chrom), r.start));
    let slot = |r: &TestRec| {
        if !per_strand {
            0
        } else {
            match r.strand {
                '+' => 0,
                '-' => 1,
                _ => 2,
            }
        }
    };
    // (rank, start, end, slot, chrom, count, names)
    let mut rows: Vec<MergeRow> = Vec::new();
    let mut keys: Vec<(String, usize)> = all.iter().map(|r| (r.chrom.clone(), slot(r))).collect();
    keys.sort();
    keys.dedup();
    for (chrom, s) in keys {
        let mut cur: Option<(u64, u64, Vec<String>)> = None;
        let rank = chrom_rank(&chrom);
        let close = |c: Option<(u64, u64, Vec<String>)>, rows: &mut Vec<MergeRow>| {
            if let Some((st, en, names)) = c {
                rows.push((rank, st, en, s, chrom.clone(), names.len(), names));
            }
        };
        for r in all.iter().filter(|r| r.chrom == chrom && slot(r) == s) {
            match &mut cur {
                Some((_, end, names))
                    if (r.start as i128) <= (*end as i128) + (max_gap as i128) =>
                {
                    *end = (*end).max(r.end);
                    names.push(r.name.clone());
                }
                _ => {
                    close(cur.take(), &mut rows);
                    cur = Some((r.start, r.end, vec![r.name.clone()]));
                }
            }
        }
        close(cur.take(), &mut rows);
    }
    rows.sort_by(|a, b| (a.0, a.1, a.2, a.3).cmp(&(b.0, b.1, b.2, b.3)));
    rows.into_iter()
        .map(|(_, st, en, _, chrom, n, names)| {
            format!("{chrom}\t{st}\t{en}\t{n}\t{}", names.join(","))
        })
        .collect()
}

/// Which of several equally near B `naive_closest` reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NaiveTies {
    All,
    First,
    Last,
}

/// What `naive_closest` is asked for; the fields mirror `gloci closest`.
#[derive(Debug, Clone)]
pub struct NaiveClosest {
    pub k: usize,
    pub ties: NaiveTies,
    pub distance: bool,
    pub signed: Option<gloci_core::SignFrame>,
    pub ignore_overlaps: bool,
    pub ignore_upstream: bool,
    pub ignore_downstream: bool,
    pub max_distance: Option<u64>,
    pub each_side: bool,
    pub per_file: bool,
    pub strand: StrandRule,
    pub inclusive_points: bool,
    /// One label per B file, written when there are several.
    pub labels: Vec<String>,
}

/// `closest`, by brute force, as the TSV rows gloci writes for BED inputs.
///
/// Every B of A's chromosome that meets the strand rule is a candidate
/// (distance 0 when it overlaps A, else the gap plus 1, signed in the
/// frame). The ignore flags and `max_distance` drop candidates; the rest
/// are sorted by (distance, file, input order) and the `k` nearest are
/// kept with the ties at the k-th distance (all of them, or as many as `k`
/// still needs, first or last). `per_file` does this for each file with
/// candidates, in file order; `each_side` pairs the upstream and downstream
/// selections rank by rank (an overlap counts on both sides). No row at
/// all gives the placeholder row.
pub fn naive_closest(a: &[TestRec], b: &[Vec<TestRec>], q: &NaiveClosest) -> Vec<String> {
    use gloci_core::{SignFrame, signed_distance};
    let frame = q.signed.unwrap_or(SignFrame::Ref);
    let show_distance = q.distance || q.signed.is_some();
    let sign_shown = q.signed.is_some() || q.each_side;
    let strand_ok = |x: Strand, y: Strand| match q.strand {
        StrandRule::Any => true,
        rule => {
            x != Strand::Unknown && y != Strand::Unknown && ((x == y) == (rule == StrandRule::Same))
        }
    };
    // (distance, file, index in file, signed distance)
    type Cand = (u64, usize, usize, i64);
    let pick = |c: &[Cand]| -> Vec<Cand> {
        if c.len() <= q.k {
            return c.to_vec();
        }
        let d = c[q.k - 1].0;
        let mut out: Vec<Cand> = c.iter().filter(|x| x.0 < d).copied().collect();
        let need = q.k - out.len();
        let tied: Vec<Cand> = c.iter().filter(|x| x.0 == d).copied().collect();
        match q.ties {
            NaiveTies::All => out.extend(tied),
            NaiveTies::First => out.extend(&tied[..need]),
            NaiveTies::Last => out.extend(&tied[tied.len() - need..]),
        }
        out
    };
    let labelled = b.len() > 1;
    let cell = |c: Option<Cand>| -> String {
        let mut s = String::new();
        match c {
            Some((d, f, i, signed)) => {
                if labelled {
                    s.push_str(&format!("\t{}", q.labels[f]));
                }
                s.push_str(&format!("\t{}", b[f][i].bed_line()));
                if show_distance {
                    s.push_str(&format!("\t{}", if sign_shown { signed } else { d as i64 }));
                }
            }
            None => {
                if labelled {
                    s.push_str("\t.");
                }
                // as the implementation: B's first record sets the width, 3
                // columns when the first file is empty
                if b.first().is_some_and(|f| !f.is_empty()) {
                    s.push_str("\t.\t-1\t-1\t.\t-1\t.");
                } else {
                    s.push_str("\t.\t-1\t-1");
                }
                if show_distance {
                    s.push_str("\t-1");
                }
            }
        }
        s
    };
    let mut rows = Vec::new();
    for ra in a {
        let x = ra.span();
        let mut cands: Vec<Cand> = Vec::new();
        for (f, file) in b.iter().enumerate() {
            for (i, rb) in file.iter().enumerate() {
                if rb.chrom != ra.chrom || !strand_ok(ra.strand(), rb.strand()) {
                    continue;
                }
                let y = rb.span();
                let d = distance(x, y, q.inclusive_points);
                if (d == 0 && q.ignore_overlaps) || q.max_distance.is_some_and(|m| d > m) {
                    continue;
                }
                let signed =
                    signed_distance(x, ra.strand(), y, rb.strand(), frame, q.inclusive_points);
                if (signed < 0 && q.ignore_upstream) || (signed > 0 && q.ignore_downstream) {
                    continue;
                }
                cands.push((d, f, i, signed));
            }
        }
        cands.sort();
        // the candidate groups: all of them, or each file's
        let groups: Vec<Vec<Cand>> = if q.per_file {
            (0..b.len())
                .map(|f| cands.iter().filter(|c| c.1 == f).copied().collect())
                .collect()
        } else {
            vec![cands]
        };
        let start = rows.len();
        for g in groups {
            if q.each_side {
                let up = pick(&g.iter().filter(|c| c.3 <= 0).copied().collect::<Vec<_>>());
                let down = pick(&g.iter().filter(|c| c.3 >= 0).copied().collect::<Vec<_>>());
                for r in 0..up.len().max(down.len()) {
                    rows.push(format!(
                        "{}{}{}",
                        ra.bed_line(),
                        cell(up.get(r).copied()),
                        cell(down.get(r).copied())
                    ));
                }
            } else {
                for c in pick(&g) {
                    rows.push(format!("{}{}", ra.bed_line(), cell(Some(c))));
                }
            }
        }
        if rows.len() == start {
            let second = if q.each_side {
                cell(None)
            } else {
                String::new()
            };
            rows.push(format!("{}{}{}", ra.bed_line(), cell(None), second));
        }
    }
    rows
}

/// An aggregation `naive_map` knows. Columns are those of the BED6 test
/// record as `map` sees it: 2 is its start and 3 its end (integers), 4 its
/// name, and 5 its score, which `naive_map`'s `score` gives and may be
/// missing (`.`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NaiveAgg {
    Count,
    Bases,
    Indicator,
    Sum(u8),
    Min(u8),
    Max(u8),
    /// Values joined by `,` in hit order.
    Collapse(u8),
    /// Each value once, first-seen order.
    Distinct(u8),
    First(u8),
    Last(u8),
}

/// The text of a test record as `map` reads it: `bed_line` with `score` in
/// column 5 (`.` when `None`).
pub fn scored_line(r: &TestRec, score: Option<u64>) -> String {
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}",
        r.chrom,
        r.start,
        r.end,
        r.name,
        score.map_or(".".to_string(), |s| s.to_string()),
        r.strand
    )
}

/// `map`, by brute force, as the TSV rows gloci writes for BED inputs: for
/// each A in the order given, A's line (score 0) and one value per
/// aggregation over the B records that count (the same set as
/// `naive_intersect`'s: D36, D37), in hit order. Missing values (`None`
/// scores) are skipped by every op but `Count`. `Bases` adds the overlap of
/// A with each B's own span. With no values, `Count` and `Bases` give 0,
/// `Indicator` 0, and the others `.`. An A with no counting B is dropped
/// under `skip_unmatched`.
pub fn naive_map(
    a: &[TestRec],
    b: &[Vec<TestRec>],
    criteria: &OverlapCriteria,
    window: Option<Window>,
    aggs: &[NaiveAgg],
    skip_unmatched: bool,
    score: &dyn Fn(&TestRec) -> Option<u64>,
) -> Vec<String> {
    let mut out = Vec::new();
    for ra in a {
        let hits = naive_accepted(ra, b, criteria, window);
        if hits.is_empty() && skip_unmatched {
            continue;
        }
        let mut row = ra.bed_line();
        for agg in aggs {
            let value = |rb: &TestRec, c: u8| -> Option<String> {
                match c {
                    2 => Some(rb.start.to_string()),
                    3 => Some(rb.end.to_string()),
                    4 => Some(rb.name.clone()),
                    _ => score(rb).map(|s| s.to_string()),
                }
            };
            let values =
                |c: u8| -> Vec<String> { hits.iter().filter_map(|(_, rb)| value(rb, c)).collect() };
            let nums =
                |c: u8| -> Vec<i128> { values(c).iter().map(|v| v.parse().unwrap()).collect() };
            let join = |v: Vec<String>| {
                if v.is_empty() {
                    ".".to_string()
                } else {
                    v.join(",")
                }
            };
            let cell = match *agg {
                NaiveAgg::Count => hits.len().to_string(),
                NaiveAgg::Bases => hits
                    .iter()
                    .map(|(_, rb)| shared_bases(ra.span(), rb.span()))
                    .sum::<u64>()
                    .to_string(),
                NaiveAgg::Indicator => u8::from(!hits.is_empty()).to_string(),
                NaiveAgg::Sum(c) => {
                    let n = nums(c);
                    if n.is_empty() {
                        ".".into()
                    } else {
                        n.iter().sum::<i128>().to_string()
                    }
                }
                NaiveAgg::Min(c) => nums(c).iter().min().map_or(".".into(), |v| v.to_string()),
                NaiveAgg::Max(c) => nums(c).iter().max().map_or(".".into(), |v| v.to_string()),
                NaiveAgg::Collapse(c) => join(values(c)),
                NaiveAgg::Distinct(c) => {
                    let mut seen: Vec<String> = Vec::new();
                    for v in values(c) {
                        if !seen.contains(&v) {
                            seen.push(v);
                        }
                    }
                    join(seen)
                }
                NaiveAgg::First(c) => values(c).first().cloned().unwrap_or(".".into()),
                NaiveAgg::Last(c) => values(c).last().cloned().unwrap_or(".".into()),
            };
            row.push('\t');
            row.push_str(&cell);
        }
        out.push(row);
    }
    out
}
