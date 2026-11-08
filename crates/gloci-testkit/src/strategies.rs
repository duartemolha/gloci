//! proptest strategies for interval data.
//!
//! Records cover points, bookended features, duplicates, nesting, huge
//! intervals (up to 2^40), mixed strands, and genomes in which some
//! chromosomes have no records. Regions overlap, repeat and come in any
//! order, and some name a chromosome that never has records.

use gloci_core::{Interval, Strand, Variant, VcfSpan, footprint};
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

/// What the relation references (`naive_hits`, `naive_neighbours`) need of
/// a test record, whatever its format.
pub trait Feature {
    fn chrom(&self) -> &str;
    /// The position its file is sorted by: BED start, VCF POS-1.
    fn key(&self) -> u64;
    /// Its footprint (spec §4.2). A VCF record's can start after its key.
    fn span(&self) -> Interval;
    fn strand(&self) -> Strand;
    /// The line as written to its file, without a newline.
    fn line(&self) -> String;
}

impl Feature for TestRec {
    fn chrom(&self) -> &str {
        &self.chrom
    }

    fn key(&self) -> u64 {
        self.start
    }

    fn span(&self) -> Interval {
        self.interval()
    }

    /// As gloci reads BED column 6 (D54).
    fn strand(&self) -> Strand {
        match self.strand {
            '+' => Strand::Plus,
            '-' => Strand::Minus,
            _ => Strand::Unknown,
        }
    }

    fn line(&self) -> String {
        self.bed_line()
    }
}

/// One VCF record as [`vcf_relation_inputs`] writes it:
/// `CHROM POS ID REF ALT . . INFO`, INFO being `END=…` or `.`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VcfRec {
    pub chrom: String,
    /// 1-based POS.
    pub pos: u64,
    pub id: String,
    pub reference: String,
    pub alt: String,
    /// INFO/END, for a `<DEL>`.
    pub end: Option<u64>,
}

impl Feature for VcfRec {
    fn chrom(&self) -> &str {
        &self.chrom
    }

    fn key(&self) -> u64 {
        self.pos - 1
    }

    /// The allele span (`--vcf-span allele`), by gloci-core's footprint
    /// rule, which the reader applies too.
    fn span(&self) -> Interval {
        let v = Variant {
            pos: self.pos,
            reference: self.reference.as_bytes(),
            alts: self.alt.as_bytes(),
            end: self.end,
            svlen: None,
        };
        footprint(&v, VcfSpan::Allele).expect("strategy VCF records have a footprint")
    }

    /// VCF records have no strand (D54).
    fn strand(&self) -> Strand {
        Strand::Unknown
    }

    fn line(&self) -> String {
        let info = self.end.map_or(".".to_string(), |e| format!("END={e}"));
        format!(
            "{}\t{}\t{}\t{}\t{}\t.\t.\t{info}",
            self.chrom, self.pos, self.id, self.reference, self.alt
        )
    }
}

/// A VCF record shape with a known allele span, for a record whose key
/// (POS-1) is `k`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VcfShape {
    /// `A → G`: `[k, k+1)`.
    Snv,
    /// `A^prefix T^len → A^prefix`, a deletion trimmed to the deleted
    /// bases: `[k+prefix, k+prefix+len)`.
    Del { prefix: u64, len: u64 },
    /// `A → A C^len`: a point after POS, `[k+1, k+1)`.
    Ins { len: u64 },
    /// `N → <DEL>` with `END=k+1+len`: `[k+1, k+1+len)`.
    Sv { len: u64 },
}

/// The VCF record of `shape` at key `key` (POS `key + 1`).
pub fn vcf_rec(chrom: &str, key: u64, shape: VcfShape, id: String) -> VcfRec {
    let bases = |b: &str, n: u64| b.repeat(n as usize);
    let (reference, alt, end) = match shape {
        VcfShape::Snv => ("A".to_string(), "G".to_string(), None),
        VcfShape::Del { prefix, len } => (
            bases("A", prefix) + &bases("T", len),
            bases("A", prefix),
            None,
        ),
        VcfShape::Ins { len } => ("A".to_string(), "A".to_string() + &bases("C", len), None),
        VcfShape::Sv { len } => ("N".to_string(), "<DEL>".to_string(), Some(key + 1 + len)),
    };
    VcfRec {
        chrom: chrom.to_string(),
        pos: key + 1,
        id,
        reference,
        alt,
        end,
    }
}

/// A VCF record whose allele span is `[start, end)`, in a shape chosen by
/// `sel`, if one exists: a point is an insertion after POS `start`; one
/// base may be an SNV; up to 64 bases may be a deletion trimmed to them,
/// with a 1–4 base prefix; anything else is a `<DEL>` at POS `start`.
/// Only an SNV can start at 0.
pub fn vcf_with_span(chrom: &str, start: u64, end: u64, sel: u64, id: String) -> Option<VcfRec> {
    let (key, shape) = if start == end {
        (start.checked_sub(1)?, VcfShape::Ins { len: 1 + sel % 3 })
    } else if end - start == 1 && (start == 0 || sel.is_multiple_of(3)) {
        (start, VcfShape::Snv)
    } else if start == 0 {
        return None;
    } else if end - start <= 64 && sel % 3 == 1 {
        let prefix = 1 + (sel / 3) % start.min(4);
        (
            start - prefix,
            VcfShape::Del {
                prefix,
                len: end - start,
            },
        )
    } else {
        (start - 1, VcfShape::Sv { len: end - start })
    };
    Some(vcf_rec(chrom, key, shape, id))
}

/// A relation test record in either format.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RelRec {
    Bed(TestRec),
    Vcf(VcfRec),
}

impl RelRec {
    fn inner(&self) -> &dyn Feature {
        match self {
            RelRec::Bed(r) => r,
            RelRec::Vcf(r) => r,
        }
    }
}

impl Feature for RelRec {
    fn chrom(&self) -> &str {
        self.inner().chrom()
    }

    fn key(&self) -> u64 {
        self.inner().key()
    }

    fn span(&self) -> Interval {
        self.inner().span()
    }

    fn strand(&self) -> Strand {
        self.inner().strand()
    }

    fn line(&self) -> String {
        self.inner().line()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TestFormat {
    Bed,
    Vcf,
}

/// The header [`TestFile::text`] writes before VCF records.
pub const VCF_HEADER: &str =
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// One input of a relation property test. Every record is in `format`.
#[derive(Debug, Clone)]
pub struct TestFile {
    pub format: TestFormat,
    pub recs: Vec<RelRec>,
}

impl TestFile {
    /// Stable sort by `(chrom_rank(chrom), key, span end)`, an order gloci
    /// accepts in either format.
    pub fn sort(&mut self, chrom_rank: impl Fn(&str) -> u32) {
        self.recs
            .sort_by_cached_key(|r| (chrom_rank(r.chrom()), r.key(), r.span().end));
    }

    /// The file's text, with [`VCF_HEADER`] first for VCF.
    pub fn text(&self) -> String {
        let mut out = match self.format {
            TestFormat::Bed => String::new(),
            TestFormat::Vcf => VCF_HEADER.to_string(),
        };
        for r in &self.recs {
            out.push_str(&r.line());
            out.push('\n');
        }
        out
    }

    /// The lines before the first record: record `i` is on line
    /// `header_lines() + i + 1`.
    pub fn header_lines(&self) -> u64 {
        match self.format {
            TestFormat::Bed => 0,
            TestFormat::Vcf => VCF_HEADER.lines().count() as u64,
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

/// The interval of a B record derived from A's span `p` by `how`, with
/// `d` (0..1000) choosing distances and offsets.
fn derived_span(p: Interval, how: Derived, d: u64) -> (u64, u64) {
    let inside = p.start + d % (p.end - p.start + 1);
    match how {
        Derived::SameInterval => (p.start, p.end),
        Derived::BookendLeft => (p.start.saturating_sub(d + 1), p.start),
        Derived::BookendRight => (p.end, p.end + d + 1),
        Derived::PointOnStart => (p.start, p.start),
        Derived::PointOnEnd => (p.end, p.end),
        Derived::PointInside => (inside, inside),
        Derived::Nested => (inside, inside + (d / 7) % (p.end - inside + 1)),
        Derived::Containing => (p.start.saturating_sub(d), p.end + d),
    }
}

/// Picks of records to derive from A's: which A, how, a distance, a strand.
fn derived_picks() -> impl Strategy<Value = Vec<(Index, Derived, u64, char)>> {
    prop::collection::vec((any::<Index>(), derived(), 0u64..1000, strand()), 0..=12)
}

/// Inputs for an A-against-B relation over [`CHROMS`]. A and each B are
/// drawn like [`bed_records`] (points, bookended, duplicate, nested and
/// huge features with mixed strands, on a random subset of chromosomes),
/// and each B also gets up to 12 records derived from A's: the same
/// interval, bookended on either side, a point on either edge or inside,
/// nested inside, or containing it. That makes the edge cases of §4.3
/// meet across files, not only within one.
pub fn relation_inputs() -> impl Strategy<Value = RelationInputs> {
    let b_file = (bed_records(), derived_picks());
    (bed_records(), prop::collection::vec(b_file, 1..=3)).prop_map(|(a, files)| {
        let b = files
            .into_iter()
            .map(|(mut recs, picks)| {
                if a.is_empty() {
                    return recs;
                }
                for (i, (idx, how, d, strand)) in picks.into_iter().enumerate() {
                    let p = &a[idx.index(a.len())];
                    let (start, end) = derived_span(p.interval(), how, d);
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

/// A and 1–3 B inputs of a relation, each in BED or VCF.
#[derive(Debug, Clone)]
pub struct MixedInputs {
    pub a: TestFile,
    pub b: Vec<TestFile>,
}

impl From<RelationInputs> for MixedInputs {
    /// The same records, every file BED.
    fn from(input: RelationInputs) -> MixedInputs {
        let bed = |recs: Vec<TestRec>| TestFile {
            format: TestFormat::Bed,
            recs: recs.into_iter().map(RelRec::Bed).collect(),
        };
        MixedInputs {
            a: bed(input.a),
            b: input.b.into_iter().map(bed).collect(),
        }
    }
}

fn vcf_shape(huge: bool) -> impl Strategy<Value = VcfShape> {
    prop_oneof![
        3 => Just(VcfShape::Snv),
        4 => (1u64..=12, 1u64..=8).prop_map(|(prefix, len)| VcfShape::Del { prefix, len }),
        2 => (1u64..=3).prop_map(|len| VcfShape::Ins { len }),
        2 => length(huge).prop_map(|len| VcfShape::Sv { len }),
    ]
}

/// Unsorted VCF records on `chr1` and a random subset of the other
/// [`CHROMS`] (sort them with [`TestFile::sort`]). Most spans start after
/// the key (POS-1): trimmed deletions up to 12 bases after it, insertions
/// (a point after POS) and `<DEL>`s with END (huge ones too). The rest are
/// SNVs, whose span starts at the key. Keys are dense, so they repeat and
/// fall inside other records' prefixes.
pub fn vcf_records() -> impl Strategy<Value = Vec<VcfRec>> {
    let others = prop::sample::subsequence(CHROMS[1..].to_vec(), 0..CHROMS.len());
    (prop::bool::weighted(0.3), others)
        .prop_flat_map(|(huge, others)| {
            let mut genome = vec![CHROMS[0]];
            genome.extend(others);
            let item = (
                prop::sample::select(genome),
                position(huge),
                vcf_shape(huge),
            );
            prop::collection::vec(item, 0..40)
        })
        .prop_map(|items| {
            items
                .into_iter()
                .enumerate()
                .map(|(i, (chrom, key, shape))| vcf_rec(chrom, key, shape, format!("v{i}")))
                .collect()
        })
}

/// A relation input in BED (weight `bed`) or VCF (weight `vcf`).
fn test_file(bed: u32, vcf: u32) -> impl Strategy<Value = TestFile> {
    prop_oneof![
        bed => bed_records().prop_map(|recs| TestFile {
            format: TestFormat::Bed,
            recs: recs.into_iter().map(RelRec::Bed).collect(),
        }),
        vcf => vcf_records().prop_map(|recs| TestFile {
            format: TestFormat::Vcf,
            recs: recs.into_iter().map(RelRec::Vcf).collect(),
        }),
    ]
}

/// Like [`relation_inputs`], but A (VCF 3 times in 5) and each B (VCF
/// every other time) are BED or VCF, with VCF records from
/// [`vcf_records`]: their span starts after their key, so key order is not
/// span order. Each B also gets up to 12 records derived from A's span as
/// in [`relation_inputs`], written in B's format ([`vcf_with_span`]; a
/// span no VCF record can have is skipped).
pub fn vcf_relation_inputs() -> impl Strategy<Value = MixedInputs> {
    let b_file = (test_file(1, 1), derived_picks());
    (test_file(2, 3), prop::collection::vec(b_file, 1..=3)).prop_map(|(a, files)| {
        let b = files
            .into_iter()
            .map(|(mut file, picks)| {
                if a.recs.is_empty() {
                    return file;
                }
                for (i, (idx, how, d, strand)) in picks.into_iter().enumerate() {
                    let p = &a.recs[idx.index(a.recs.len())];
                    let (start, end) = derived_span(p.span(), how, d);
                    let chrom = p.chrom().to_string();
                    let name = format!("d{i}");
                    let rec = match file.format {
                        TestFormat::Bed => Some(RelRec::Bed(TestRec {
                            chrom,
                            start,
                            end,
                            strand,
                            name,
                        })),
                        TestFormat::Vcf => {
                            vcf_with_span(&chrom, start, end, d, name).map(RelRec::Vcf)
                        }
                    };
                    file.recs.extend(rec);
                }
                file
            })
            .collect();
        MixedInputs { a, b }
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
