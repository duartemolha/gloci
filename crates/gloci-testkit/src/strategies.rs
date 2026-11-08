//! proptest strategies for interval data.
//!
//! Records cover points, bookended features, duplicates, nesting, huge
//! intervals (up to 2^40), mixed strands, and genomes in which some
//! chromosomes have no records. Regions overlap, repeat and come in any
//! order, and some name a chromosome that never has records.

use gloci_core::Interval;
use proptest::prelude::*;
use proptest::sample::Index;

/// Chromosomes records are drawn from. Each case uses `chr1` plus a random
/// subset of the others, so some chromosomes of the pool have no records.
pub const CHROMS: &[&str] = &["chr1", "chr2", "chr10", "chrX", "chrM"];

/// A chromosome regions may name that never has records.
pub const ABSENT_CHROM: &str = "chrUn";

/// The largest coordinate the strategies generate.
pub const HUGE: u64 = 1 << 40;

/// Coordinates on and next to 16 kb index window (and bin) boundaries,
/// where the linear index prunes chunks.
fn near_window_edge() -> impl Strategy<Value = u64> {
    (prop::sample::select(vec![1u64, 2, 8, 64]), -2i64..=2)
        .prop_map(|(k, d)| ((k << 14) as i64 + d) as u64)
}

/// One BED6 record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TestRec {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub name: String,
}

impl TestRec {
    /// `chrom start end name 0 strand`, tab-separated, without a newline.
    pub fn bed_line(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t0\t{}",
            self.chrom, self.start, self.end, self.name, self.strand
        )
    }

    pub fn interval(&self) -> Interval {
        Interval {
            start: self.start,
            end: self.end,
        }
    }
}

/// Start positions. Without `huge`, every record ends below 2^29 (so a TBI
/// index suffices); with it, some cross 2^29 or reach 2^40 (CSI). The dense
/// low range makes records and regions share coordinates often, so points
/// land on edges (where `--inclusive-points` matters).
fn position(huge: bool) -> BoxedStrategy<u64> {
    if huge {
        prop_oneof![
            5 => 0u64..16,
            2 => 0u64..100_000,
            2 => near_window_edge(),
            2 => (1u64 << 29) - 100..(1u64 << 29) + 100,
            1 => 0u64..=HUGE,
        ]
        .boxed()
    } else {
        prop_oneof![
            6 => 0u64..16,
            3 => 0u64..100_000,
            2 => near_window_edge(),
            1 => (1u64 << 29) - 200_000..(1u64 << 29) - 100_000,
        ]
        .boxed()
    }
}

fn length(huge: bool) -> BoxedStrategy<u64> {
    if huge {
        prop_oneof![
            2 => Just(0u64),
            4 => 1u64..6,
            2 => 6u64..50_000,
            1 => 0u64..=HUGE,
        ]
        .boxed()
    } else {
        prop_oneof![
            2 => Just(0u64),
            5 => 1u64..6,
            2 => 6u64..50_000,
        ]
        .boxed()
    }
}

/// How a record relates to the one generated before it.
#[derive(Debug, Clone, Copy)]
enum Kind {
    Fresh,
    /// Starts where the previous record ends.
    Bookended,
    /// The same line again.
    Duplicate,
    /// The same interval with another name and strand.
    SameInterval,
    /// Inside the previous record.
    Nested,
    /// A point inside (or on an edge of) the previous record.
    PointInside,
}

fn kind() -> impl Strategy<Value = Kind> {
    prop_oneof![
        6 => Just(Kind::Fresh),
        2 => Just(Kind::Bookended),
        1 => Just(Kind::Duplicate),
        1 => Just(Kind::SameInterval),
        2 => Just(Kind::Nested),
        1 => Just(Kind::PointInside),
    ]
}

fn strand() -> impl Strategy<Value = char> {
    prop_oneof![Just('+'), Just('-'), Just('.')]
}

/// Unsorted BED records (sort them with `naive::sort` before writing).
pub fn bed_records() -> impl Strategy<Value = Vec<TestRec>> {
    let others = prop::sample::subsequence(CHROMS[1..].to_vec(), 0..CHROMS.len());
    (prop::bool::weighted(0.3), others)
        .prop_flat_map(|(huge, others)| {
            let mut genome = vec![CHROMS[0]];
            genome.extend(others);
            let item = (
                kind(),
                prop::sample::select(genome),
                position(huge),
                length(huge),
                strand(),
            );
            (Just(huge), prop::collection::vec(item, 0..40))
        })
        .prop_map(|(huge, items)| {
            let cap = if huge { HUGE } else { (1 << 29) - 1 };
            let mut out: Vec<TestRec> = Vec::with_capacity(items.len());
            for (i, (kind, chrom, pos, len, strand)) in items.into_iter().enumerate() {
                let start = pos.min(cap);
                let fresh = TestRec {
                    chrom: chrom.to_string(),
                    start,
                    end: start.saturating_add(len).min(cap),
                    strand,
                    name: format!("r{i}"),
                };
                let rec = match (kind, out.last()) {
                    (Kind::Fresh, _) | (_, None) => fresh,
                    (Kind::Bookended, Some(p)) => TestRec {
                        chrom: p.chrom.clone(),
                        start: p.end,
                        end: p.end.saturating_add(len % 1000).min(cap),
                        ..fresh
                    },
                    (Kind::Duplicate, Some(p)) => p.clone(),
                    (Kind::SameInterval, Some(p)) => TestRec {
                        chrom: p.chrom.clone(),
                        start: p.start,
                        end: p.end,
                        ..fresh
                    },
                    (Kind::Nested, Some(p)) => {
                        let s = p.start + pos % (p.end - p.start + 1);
                        TestRec {
                            chrom: p.chrom.clone(),
                            start: s,
                            end: s + len % (p.end - s + 1),
                            ..fresh
                        }
                    }
                    (Kind::PointInside, Some(p)) => {
                        let s = p.start + pos % (p.end - p.start + 1);
                        TestRec {
                            chrom: p.chrom.clone(),
                            start: s,
                            end: s,
                            ..fresh
                        }
                    }
                };
                out.push(rec);
            }
            out
        })
}

/// A and 1–3 B record sets for a relation property test, unsorted like
/// [`bed_records`] (sort each with `naive::sort` before writing).
#[derive(Debug, Clone)]
pub struct RelationInputs {
    pub a: Vec<TestRec>,
    pub b: Vec<Vec<TestRec>>,
}

/// How a B record derived from an A record relates to it.
#[derive(Debug, Clone, Copy)]
enum Derived {
    SameInterval,
    /// Ends where A starts.
    BookendLeft,
    /// Starts where A ends.
    BookendRight,
    PointOnStart,
    PointOnEnd,
    PointInside,
    Nested,
    Containing,
}

fn derived() -> impl Strategy<Value = Derived> {
    prop_oneof![
        Just(Derived::SameInterval),
        Just(Derived::BookendLeft),
        Just(Derived::BookendRight),
        Just(Derived::PointOnStart),
        Just(Derived::PointOnEnd),
        Just(Derived::PointInside),
        Just(Derived::Nested),
        Just(Derived::Containing),
    ]
}

/// Inputs for an A-against-B relation over [`CHROMS`]. A and each B are
/// drawn like [`bed_records`] (points, bookended, duplicate, nested and
/// huge features with mixed strands, on a random subset of chromosomes),
/// and each B also gets up to 12 records derived from A's: the same
/// interval, bookended on either side, a point on either edge or inside,
/// nested inside, or containing it. That makes the edge cases of §4.3
/// meet across files, not only within one.
pub fn relation_inputs() -> impl Strategy<Value = RelationInputs> {
    let picks = prop::collection::vec((any::<Index>(), derived(), 0u64..1000, strand()), 0..=12);
    let b_file = (bed_records(), picks);
    (bed_records(), prop::collection::vec(b_file, 1..=3)).prop_map(|(a, files)| {
        let b = files
            .into_iter()
            .map(|(mut recs, picks)| {
                if a.is_empty() {
                    return recs;
                }
                for (i, (idx, how, d, strand)) in picks.into_iter().enumerate() {
                    let p = &a[idx.index(a.len())];
                    let inside = p.start + d % (p.end - p.start + 1);
                    let (start, end) = match how {
                        Derived::SameInterval => (p.start, p.end),
                        Derived::BookendLeft => (p.start.saturating_sub(d + 1), p.start),
                        Derived::BookendRight => (p.end, p.end + d + 1),
                        Derived::PointOnStart => (p.start, p.start),
                        Derived::PointOnEnd => (p.end, p.end),
                        Derived::PointInside => (inside, inside),
                        Derived::Nested => (inside, inside + (d / 7) % (p.end - inside + 1)),
                        Derived::Containing => (p.start.saturating_sub(d), p.end + d),
                    };
                    recs.push(TestRec {
                        chrom: p.chrom.clone(),
                        start,
                        end,
                        strand,
                        name: format!("d{i}"),
                    });
                }
                recs
            })
            .collect();
        RelationInputs { a, b }
    })
}

fn region_start() -> impl Strategy<Value = u64> {
    prop_oneof![
        7 => 0u64..16,
        1 => 0u64..100_000,
        2 => near_window_edge(),
        1 => (1u64 << 29) - 100..(1u64 << 29) + 100,
        1 => 0u64..=HUGE,
    ]
}

fn region_len() -> impl Strategy<Value = u64> {
    prop_oneof![
        3 => Just(0u64),
        4 => 1u64..8,
        2 => 8u64..100_000,
        1 => 0u64..=2 * HUGE,
        1 => Just(u64::MAX), // to the end of the chromosome
    ]
}

/// 1..=`max` query regions, possibly overlapping, repeated, out of order,
/// zero-length, or on a chromosome with no records.
pub fn regions(max: usize) -> impl Strategy<Value = Vec<(Vec<u8>, Interval)>> {
    let chrom = prop_oneof![
        4 => Just(CHROMS[0]),
        4 => prop::sample::select(CHROMS),
        1 => Just(ABSENT_CHROM),
    ];
    let region = (chrom, region_start(), region_len()).prop_map(|(c, s, l)| {
        let iv = Interval {
            start: s,
            end: s.saturating_add(l),
        };
        (c.as_bytes().to_vec(), iv)
    });
    (
        prop::collection::vec(region, 1..=max.max(1)),
        any::<Index>(),
        any::<bool>(),
    )
        .prop_map(|(mut v, i, dup)| {
            // replace the last region with a copy of an earlier one
            if dup && v.len() >= 2 {
                let k = i.index(v.len() - 1);
                let last = v.len() - 1;
                v[last] = v[k].clone();
            }
            v
        })
}
