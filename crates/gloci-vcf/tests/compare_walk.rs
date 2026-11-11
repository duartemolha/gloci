//! The walk of `vcf compare` over two record lists (decision D163).

use gloci_core::ChromId;
use gloci_vcf::compare::{MatchRule, Side, Site, SiteClass, SiteSource, Walk};
use std::cmp::Ordering;
use std::collections::VecDeque;
use std::convert::Infallible;

/// A record: (chromosome id, POS, REF, ALT column).
type Rec = (u32, u64, &'static str, &'static str);

struct Lists {
    side: [VecDeque<Site>; 2],
    /// Rank of each chromosome id in the walk's order.
    rank: Vec<u32>,
    /// Records handed to the walk so far.
    pulled: std::cell::Cell<usize>,
}

impl SiteSource for Lists {
    type Err = Infallible;
    fn next(&mut self, side: Side) -> Result<Option<Site>, Infallible> {
        let s = self.side[usize::from(side == Side::Second)].pop_front();
        self.pulled
            .set(self.pulled.get() + usize::from(s.is_some()));
        Ok(s)
    }
    fn cmp_chrom(&self, a: ChromId, b: ChromId) -> Ordering {
        self.rank[a.0 as usize].cmp(&self.rank[b.0 as usize])
    }
}

fn sites(v: &[Rec]) -> VecDeque<Site> {
    v.iter()
        .map(|&(c, p, r, a)| Site::bare(ChromId(c), p, r, a))
        .collect()
}

/// Rows as `B 10=10`, `O 10=10`, `1 10`, `2 10` (`c:pos` when the chromosome is not 0).
fn walk(rule: MatchRule, rank: Vec<u32>, a: &[Rec], b: &[Rec]) -> Vec<String> {
    let src = Lists {
        side: [sites(a), sites(b)],
        rank,
        pulled: Default::default(),
    };
    let mut w = Walk::new(src, rule);
    let at = |s: &Site| format!("{}:{}", s.chrom.0, s.pos);
    let mut out = Vec::new();
    while let Some(r) = w.next_row().unwrap() {
        out.push(match (r.class, &r.first, &r.second) {
            (SiteClass::First, Some(x), _) => format!("1 {}", at(x)),
            (SiteClass::Second, _, Some(y)) => format!("2 {}", at(y)),
            (c, Some(x), Some(y)) => format!("{} {}={}", c.letter(), at(x), at(y)),
            _ => panic!("a row without its sites"),
        });
    }
    out
}

fn alleles(a: &[Rec], b: &[Rec]) -> Vec<String> {
    walk(MatchRule::Alleles, vec![0, 1, 2], a, b)
}

fn position(a: &[Rec], b: &[Rec]) -> Vec<String> {
    walk(MatchRule::Position, vec![0, 1, 2], a, b)
}

#[test]
fn compare_walk_table() {
    // identical
    assert_eq!(
        alleles(&[(0, 10, "A", "C")], &[(0, 10, "A", "C")]),
        ["B 0:10=0:10"]
    );
    // the same POS, a different ALT: overlapping, not the same site
    assert_eq!(
        alleles(&[(0, 10, "A", "C")], &[(0, 10, "A", "G")]),
        ["O 0:10=0:10"]
    );
    // A's second ALT is not needed
    assert_eq!(
        alleles(&[(0, 10, "A", "C,G")], &[(0, 10, "A", "C")]),
        ["B 0:10=0:10"]
    );
    // each B record is taken once: the first matches, the second is B-only
    assert_eq!(
        alleles(
            &[(0, 10, "A", "C,G")],
            &[(0, 10, "A", "C"), (0, 10, "A", "G")]
        ),
        ["B 0:10=0:10", "2 0:10"]
    );
    // A's second record takes the second B record
    assert_eq!(
        alleles(
            &[(0, 10, "A", "C"), (0, 10, "A", "G")],
            &[(0, 10, "A", "G"), (0, 10, "A", "C")]
        ),
        ["B 0:10=0:10", "B 0:10=0:10"]
    );
    // two different trimmed forms of one deletion
    assert_eq!(
        alleles(&[(0, 10, "AT", "A")], &[(0, 10, "ATG", "AG")]),
        ["B 0:10=0:10"]
    );
    // case is not significant
    assert_eq!(
        alleles(&[(0, 10, "a", "c")], &[(0, 10, "A", "C")]),
        ["B 0:10=0:10"]
    );
    // the same change at another POS is another site
    assert_eq!(
        alleles(&[(0, 10, "A", "C")], &[(0, 11, "A", "C")]),
        ["1 0:10", "2 0:11"]
    );
    // leftovers at one position pair as O, the rest are alone
    assert_eq!(
        alleles(
            &[(0, 10, "A", "C"), (0, 10, "A", "G"), (0, 10, "A", "T")],
            &[(0, 10, "A", "C"), (0, 10, "A", "<DEL>")]
        ),
        ["B 0:10=0:10", "O 0:10=0:10", "1 0:10"]
    );
    // a chromosome only in B
    assert_eq!(
        alleles(
            &[(0, 10, "A", "C")],
            &[(0, 10, "A", "C"), (1, 5, "A", "C"), (1, 9, "A", "C")]
        ),
        ["B 0:10=0:10", "2 1:5", "2 1:9"]
    );
    // a chromosome only in A, and in the middle
    assert_eq!(
        alleles(
            &[(0, 1, "A", "C"), (1, 5, "A", "C"), (2, 3, "A", "C")],
            &[(0, 1, "A", "C"), (2, 3, "A", "C")]
        ),
        ["B 0:1=0:1", "1 1:5", "B 2:3=2:3"]
    );
    // a REF-only record shares no allele
    assert_eq!(
        alleles(&[(0, 10, "A", ".")], &[(0, 10, "A", ".")]),
        ["O 0:10=0:10"]
    );
}

#[test]
fn compare_walk_position() {
    // same POS and REF: B whatever the ALT
    assert_eq!(
        position(&[(0, 10, "A", "C")], &[(0, 10, "A", "G")]),
        ["B 0:10=0:10"]
    );
    // different REF: O, both consumed
    assert_eq!(
        position(
            &[(0, 10, "AT", "A"), (0, 20, "A", "C")],
            &[(0, 10, "ATG", "AG"), (0, 20, "A", "C")]
        ),
        ["O 0:10=0:10", "B 0:20=0:20"]
    );
    // a REF of N or . takes the other's
    assert_eq!(
        position(&[(0, 10, "N", "C")], &[(0, 10, "ACG", "C")]),
        ["B 0:10=0:10"]
    );
    assert_eq!(
        position(&[(0, 10, "ACG", "C")], &[(0, 10, ".", "C")]),
        ["B 0:10=0:10"]
    );
    // NNN is not N
    assert_eq!(
        position(&[(0, 10, "NNN", "C")], &[(0, 10, "ACG", "C")]),
        ["O 0:10=0:10"]
    );
    // a REF span reaching the other file's POS changes nothing in this walk:
    // the earlier record is alone, and the records still pair by POS
    assert_eq!(
        position(
            &[(0, 70, "ACG", "T"), (0, 71, "C", "T")],
            &[(0, 71, "C", "T")]
        ),
        ["1 0:70", "B 0:71=0:71"]
    );
    assert_eq!(
        position(&[(0, 100, "A", "C")], &[(0, 99, "AT", "A")]),
        ["2 0:99", "1 0:100"]
    );
    // records are taken in order, one from each file: a second record at a
    // POS is alone
    assert_eq!(
        position(
            &[(0, 200, "A", "C"), (0, 200, "A", "G")],
            &[(0, 200, "A", "C")]
        ),
        ["B 0:200=0:200", "1 0:200"]
    );
}

#[test]
fn compare_chrom_tables() {
    // header order chr1, chr2, chr10 = ids 0, 1, 2 with ranks 0, 1, 2; the
    // walk follows the ranks, not the ids or the names
    let rank = vec![0, 1, 2];
    let a: Vec<Rec> = vec![(1, 5, "A", "C"), (2, 5, "A", "C")];
    let b: Vec<Rec> = vec![(0, 5, "A", "C"), (2, 5, "A", "C")];
    assert_eq!(
        walk(MatchRule::Alleles, rank.clone(), &a, &b),
        ["2 0:5", "1 1:5", "B 2:5=2:5"]
    );
    // ids in reverse of the rank order: id 1 comes first
    let rank = vec![2, 1, 0];
    assert_eq!(
        walk(
            MatchRule::Position,
            rank,
            &[(1, 5, "A", "C"), (0, 5, "A", "C")],
            &[(1, 5, "A", "C")]
        ),
        ["B 1:5=1:5", "1 0:5"]
    );
}

#[test]
fn compare_huge_contig_count() {
    // 200,000 contigs, one record each; B has every third: the walk holds
    // one position at a time and never indexes by contig count squared
    let n = 200_000u32;
    let a: VecDeque<Site> = (0..n)
        .map(|c| Site::bare(ChromId(c), 1, "A", "C"))
        .collect();
    let b: VecDeque<Site> = (0..n)
        .step_by(3)
        .map(|c| Site::bare(ChromId(c), 1, "A", "C"))
        .collect();
    let src = Lists {
        side: [a, b],
        rank: (0..n).collect(),
        pulled: Default::default(),
    };
    let mut w = Walk::new(src, MatchRule::Alleles);
    let (mut both, mut first) = (0u32, 0u32);
    while let Some(r) = w.next_row().unwrap() {
        match r.class {
            SiteClass::Both => both += 1,
            SiteClass::First => first += 1,
            c => panic!("unexpected {c:?}"),
        }
    }
    assert_eq!((both, first), (66_667, 133_333));
}

#[test]
fn compare_walk_holds_the_records_of_one_position() {
    // 50,000 positions with two records in each file: the walk holds the
    // records of one position (4) plus one read ahead per file, never more
    for rule in [MatchRule::Alleles, MatchRule::Position] {
        let mk = |alt: &str| -> VecDeque<Site> {
            (0..50_000u64)
                .flat_map(|p| {
                    [
                        Site::bare(ChromId(0), p + 1, "A", "C"),
                        Site::bare(ChromId(0), p + 1, "A", alt),
                    ]
                })
                .collect()
        };
        let src = Lists {
            side: [mk("G"), mk("T")],
            rank: vec![0],
            pulled: Default::default(),
        };
        let mut w = Walk::new(src, rule);
        let (mut yielded, mut peak) = (0usize, 0usize);
        while let Some(r) = w.next_row().unwrap() {
            yielded += usize::from(r.first.is_some()) + usize::from(r.second.is_some());
            peak = peak.max(w.source().pulled.get() - yielded);
        }
        assert_eq!(yielded, 200_000);
        assert!(peak <= 6, "{rule:?}: {peak} records held");
    }
}
