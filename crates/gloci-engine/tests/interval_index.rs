//! The interval index (D154): spec section 4.3's point and bookend rules,
//! coordinates past 32 bits and near `u64::MAX`, windows, id order, and
//! agreement with a double loop over random entries.

use gloci_core::{ChromId, Interval, Strand, Window, overlaps};
use gloci_engine::{IndexEntry, IntervalIndex};
use proptest::prelude::*;

fn iv(start: u64, end: u64) -> Interval {
    Interval { start, end }
}

fn entry(chrom: u32, start: u64, end: u64, id: u32) -> IndexEntry {
    IndexEntry {
        chrom: ChromId(chrom),
        span: iv(start, end),
        id,
    }
}

fn query(
    idx: &IntervalIndex,
    chrom: u32,
    a: Interval,
    window: Option<Window>,
    strand: Strand,
    inclusive: bool,
) -> Vec<u32> {
    let mut out = vec![99];
    idx.query(ChromId(chrom), a, window, strand, inclusive, &mut out);
    out
}

fn plain(idx: &IntervalIndex, a: Interval, inclusive: bool) -> Vec<u32> {
    query(idx, 0, a, None, Strand::Unknown, inclusive)
}

#[test]
fn interval_index_rules() {
    // Section 4.3's table: (A, B, overlaps, overlaps with --inclusive-points).
    let table = [
        ((10, 10), (5, 15), true, true),
        ((5, 5), (5, 15), false, true),
        ((15, 15), (5, 15), false, true),
        ((10, 10), (10, 10), true, true),
        ((10, 10), (11, 11), false, false),
        ((5, 10), (10, 15), false, false),
        ((5, 11), (10, 15), true, true),
        // A interval, B point.
        ((5, 15), (10, 10), true, true),
        ((5, 15), (5, 5), false, true),
        ((5, 15), (15, 15), false, true),
        ((5, 15), (4, 4), false, false),
    ];
    for (a, b, strict, inclusive) in table {
        let idx = IntervalIndex::build(vec![entry(0, b.0, b.1, 7)]);
        let want = |yes: bool| if yes { vec![7] } else { vec![] };
        assert_eq!(
            plain(&idx, iv(a.0, a.1), false),
            want(strict),
            "{a:?} {b:?}"
        );
        assert_eq!(
            plain(&idx, iv(a.0, a.1), true),
            want(inclusive),
            "{a:?} {b:?} inclusive"
        );
    }
}

#[test]
fn interval_index_large_coordinates() {
    let big = 3_000_000_000u64;
    let max = u64::MAX;
    let idx = IntervalIndex::build(vec![
        entry(0, big, big + 100, 0),
        entry(0, max - 10, max - 5, 1),
        // Ends at the top of the range.
        entry(0, max - 4, max, 2),
        // Starts at the top of the range: a point, kept apart.
        entry(0, max, max, 3),
        entry(0, 5, 6, 4),
    ]);
    assert_eq!(plain(&idx, iv(big + 50, big + 60), false), vec![0]);
    assert_eq!(plain(&idx, iv(big - 1, big + 1), false), vec![0]);
    assert_eq!(plain(&idx, iv(max - 12, max - 8), false), vec![1]);
    assert_eq!(plain(&idx, iv(max - 6, max - 3), false), vec![1, 2]);
    assert_eq!(plain(&idx, iv(max - 1, max), false), vec![2]);
    // The point at u64::MAX meets a point there, and an interval ending
    // there only with --inclusive-points.
    assert_eq!(plain(&idx, iv(max, max), false), vec![3]);
    assert_eq!(plain(&idx, iv(max, max), true), vec![2, 3]);
    assert_eq!(plain(&idx, iv(max - 1, max), true), vec![2, 3]);
    assert!(plain(&idx, iv(10, 20), false).is_empty());
}

#[test]
fn interval_index_window_pads_b() {
    let idx = IntervalIndex::build(vec![entry(0, 100, 110, 0), entry(0, 300, 310, 1)]);
    let a = iv(120, 130);
    assert!(plain(&idx, a, false).is_empty());
    let w = |left, right, strand_aware| {
        Some(Window {
            left,
            right,
            strand_aware,
        })
    };
    // Padding B by 15 on the right reaches A; by 15 on the left does not.
    assert_eq!(query(&idx, 0, a, w(0, 15, false), Strand::Plus, false), [0]);
    assert!(query(&idx, 0, a, w(15, 0, false), Strand::Plus, false).is_empty());
    // On the minus strand a strand-aware window swaps the sides.
    assert!(query(&idx, 0, a, w(0, 15, true), Strand::Minus, false).is_empty());
    assert_eq!(query(&idx, 0, a, w(15, 0, true), Strand::Minus, false), [0]);
    // Both sides: the second entry is 170 away from A's end.
    assert_eq!(
        query(&idx, 0, a, w(200, 200, false), Strand::Plus, false),
        [0, 1]
    );
    // Bookended after padding is not an overlap.
    assert!(query(&idx, 0, a, w(0, 9, false), Strand::Plus, false).is_empty());
    assert_eq!(query(&idx, 0, a, w(0, 11, false), Strand::Plus, false), [0]);
}

#[test]
fn interval_index_ids_ascend_and_chromosomes_stay_apart() {
    // Ids are not in position order, and a chromosome with no entries or
    // beyond the last one finds nothing.
    let idx = IntervalIndex::build(vec![
        entry(2, 0, 50, 5),
        entry(0, 40, 90, 3),
        entry(0, 0, 100, 9),
        entry(0, 10, 20, 1),
        entry(0, 30, 60, 2),
        entry(2, 10, 20, 0),
    ]);
    assert_eq!(plain(&idx, iv(0, 100), false), vec![1, 2, 3, 9]);
    assert!(query(&idx, 1, iv(0, 100), None, Strand::Unknown, false).is_empty());
    assert_eq!(
        query(&idx, 2, iv(0, 100), None, Strand::Unknown, false),
        vec![0, 5]
    );
    assert!(query(&idx, 7, iv(0, 100), None, Strand::Unknown, false).is_empty());
    assert!(!idx.is_empty());
    assert!(IntervalIndex::build(Vec::new()).is_empty());
    assert!(plain(&IntervalIndex::default(), iv(0, 1), false).is_empty());
}

fn coord() -> impl Strategy<Value = u64> {
    prop_oneof![
        8 => 0u64..60,
        1 => (u64::MAX - 20)..=u64::MAX,
        1 => 3_000_000_000u64..3_000_000_030,
    ]
}

fn span() -> impl Strategy<Value = Interval> {
    (coord(), 0u64..12).prop_map(|(s, len)| iv(s, s.saturating_add(len)))
}

proptest! {
    #[test]
    fn interval_index_matches_double_loop(
        entries in proptest::collection::vec((0u32..3, span()), 0..40),
        queries in proptest::collection::vec((0u32..4, span()), 1..12),
        inclusive in any::<bool>(),
        pad in proptest::option::of((0u64..20, 0u64..20, any::<bool>())),
        minus in any::<bool>(),
    ) {
        let idx = IntervalIndex::build(
            entries
                .iter()
                .enumerate()
                .map(|(i, &(c, s))| entry(c, s.start, s.end, i as u32))
                .collect(),
        );
        let window = pad.map(|(left, right, strand_aware)| Window { left, right, strand_aware });
        let strand = if minus { Strand::Minus } else { Strand::Plus };
        for (c, a) in queries {
            let want: Vec<u32> = entries
                .iter()
                .enumerate()
                .filter(|(_, (ec, b))| {
                    *ec == c
                        && match window {
                            None => overlaps(a, *b, inclusive),
                            Some(w) => overlaps(a, w.apply(*b, strand), inclusive),
                        }
                })
                .map(|(i, _)| i as u32)
                .collect();
            prop_assert_eq!(query(&idx, c, a, window, strand, inclusive), want);
        }
    }
}

#[test]
fn interval_index_many_chromosomes() {
    // 200,000 contigs with one entry each: a build and a query cost no more
    // than the entries.
    let n = 200_000u32;
    let idx = IntervalIndex::build((0..n).rev().map(|c| entry(c, 10, 20, c)).collect());
    for c in [0, 1, 99_999, n - 1] {
        assert_eq!(
            query(&idx, c, iv(15, 16), None, Strand::Unknown, false),
            vec![c]
        );
    }
    assert!(query(&idx, n, iv(15, 16), None, Strand::Unknown, false).is_empty());
}
