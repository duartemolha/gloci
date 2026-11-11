//! The walk over two sorted record streams (decision D163).
//!
//! `--match position` is the walk behind vcftools' `--diff-site-discordance`:
//! records are taken in order, one from each file; the same POS is `B` (same
//! REF; a REF of `N` or `.` takes the other's) or `O` (different REF, both
//! consumed); at different POS the earlier record is alone (`1` or `2`). A REF
//! span that reaches the other file's next POS makes no difference here
//! (vcftools' `--diff-site` walk calls that `O`; gloci does not write that
//! file, and the discordance walk keeps the two records apart). `--match alleles` groups the records of each file at one
//! CHROM and POS: each first-file record takes the first unmatched
//! second-file record that shares an allele (after §4.2 trimming) as `B`;
//! the leftovers pair up in order as `O` and the rest are `1` or `2`.

use super::genotypes::Gt;
use gloci_core::{ChromId, trimmed_alleles};
use std::cmp::Ordering;
use std::collections::VecDeque;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatchRule {
    Alleles,
    Position,
}

/// How a row of the walk relates the two files (`B 1 2 O` in vcftools).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SiteClass {
    /// The same site in both files.
    Both,
    /// Only in the first file.
    First,
    /// Only in the second file.
    Second,
    /// The two sites overlap without being the same.
    Overlap,
}

impl SiteClass {
    pub fn letter(self) -> char {
        match self {
            SiteClass::Both => 'B',
            SiteClass::First => '1',
            SiteClass::Second => '2',
            SiteClass::Overlap => 'O',
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Side {
    First,
    Second,
}

/// A record as the walk and the comparison see it.
#[derive(Clone, Debug)]
pub struct Site {
    /// The chromosome's name, as written in the file.
    pub chrom_name: Box<[u8]>,
    /// Its id in the walk's table; set by the [`SiteSource`].
    pub chrom: ChromId,
    /// 1-based.
    pub pos: u64,
    /// REF, upper-cased.
    pub ref_: Vec<u8>,
    /// The raw ALT column, upper-cased.
    pub alts: Vec<u8>,
    /// Record number in its file, for messages.
    pub rec_no: u64,
    /// The genotypes of the samples both files have, in one shared order.
    pub gts: Vec<Gt>,
}

impl Site {
    /// A site with no genotypes, for the walk alone.
    pub fn bare(chrom: ChromId, pos: u64, ref_: &str, alts: &str) -> Site {
        Site {
            chrom_name: Box::default(),
            chrom,
            pos,
            ref_: ref_.to_ascii_uppercase().into_bytes(),
            alts: alts.to_ascii_uppercase().into_bytes(),
            rec_no: 0,
            gts: Vec::new(),
        }
    }

    /// One past the last base of REF.
    pub fn end(&self) -> u64 {
        self.pos + self.ref_.len() as u64
    }
}

/// A REF vcftools lets the other file's REF replace.
pub fn is_unknown_ref(r: &[u8]) -> bool {
    r == b"N" || r == b"."
}

/// Do the two records share a trimmed (position, REF, ALT)? Symbolic,
/// breakend, `*` and missing ALTs have no sequence and never match.
pub fn share_allele(a: &Site, b: &Site) -> bool {
    fn usable(alt: &[u8]) -> bool {
        !alt.is_empty()
            && !alt
                .iter()
                .any(|c| matches!(c, b'<' | b'[' | b']' | b'*' | b'.'))
    }
    a.alts.split(|&c| c == b',').filter(|x| usable(x)).any(|x| {
        let (ap, ar, aa) = trimmed_alleles(a.pos.saturating_sub(1), &a.ref_, x);
        b.alts.split(|&c| c == b',').filter(|y| usable(y)).any(|y| {
            let (bp, br, ba) = trimmed_alleles(b.pos.saturating_sub(1), &b.ref_, y);
            ap == bp && ar.eq_ignore_ascii_case(br) && aa.eq_ignore_ascii_case(ba)
        })
    })
}

/// The two record streams a [`Walk`] reads, each already checked for order.
pub trait SiteSource {
    type Err;
    /// The next record of one file, `None` at its end.
    fn next(&mut self, side: Side) -> Result<Option<Site>, Self::Err>;
    /// The order of the walk's chromosomes.
    fn cmp_chrom(&self, a: ChromId, b: ChromId) -> Ordering;
}

/// One step of the walk.
#[derive(Debug)]
pub struct Row {
    pub class: SiteClass,
    pub first: Option<Site>,
    pub second: Option<Site>,
}

pub struct Walk<S: SiteSource> {
    src: S,
    rule: MatchRule,
    head: [Option<Site>; 2],
    eof: [bool; 2],
    out: VecDeque<Row>,
}

impl<S: SiteSource> Walk<S> {
    pub fn new(src: S, rule: MatchRule) -> Walk<S> {
        Walk {
            src,
            rule,
            head: [None, None],
            eof: [false, false],
            out: VecDeque::new(),
        }
    }

    pub fn source(&self) -> &S {
        &self.src
    }

    fn side(i: usize) -> Side {
        if i == 0 { Side::First } else { Side::Second }
    }

    fn fill(&mut self, i: usize) -> Result<(), S::Err> {
        if self.head[i].is_none() && !self.eof[i] {
            match self.src.next(Self::side(i))? {
                Some(s) => self.head[i] = Some(s),
                None => self.eof[i] = true,
            }
        }
        Ok(())
    }

    fn row(class: SiteClass, first: Option<Site>, second: Option<Site>) -> Row {
        Row {
            class,
            first,
            second,
        }
    }

    /// The next row, `None` when both files are read.
    pub fn next_row(&mut self) -> Result<Option<Row>, S::Err> {
        loop {
            if let Some(r) = self.out.pop_front() {
                return Ok(Some(r));
            }
            self.fill(0)?;
            self.fill(1)?;
            let (a, b) = match (&self.head[0], &self.head[1]) {
                (None, None) => return Ok(None),
                (Some(_), None) => {
                    return Ok(Some(Self::row(SiteClass::First, self.head[0].take(), None)));
                }
                (None, Some(_)) => {
                    return Ok(Some(Self::row(
                        SiteClass::Second,
                        None,
                        self.head[1].take(),
                    )));
                }
                (Some(a), Some(b)) => (a, b),
            };
            // The earlier file's record alone when the chromosome or the
            // position differs.
            let ord = if a.chrom == b.chrom {
                a.pos.cmp(&b.pos)
            } else {
                self.src.cmp_chrom(a.chrom, b.chrom)
            };
            if self.rule == MatchRule::Position && ord == Ordering::Equal && a.chrom == b.chrom {
                let class =
                    if a.ref_ == b.ref_ || is_unknown_ref(&a.ref_) || is_unknown_ref(&b.ref_) {
                        SiteClass::Both
                    } else {
                        SiteClass::Overlap
                    };
                return Ok(Some(Self::row(
                    class,
                    self.head[0].take(),
                    self.head[1].take(),
                )));
            }
            match ord {
                Ordering::Less => {
                    return Ok(Some(Self::row(SiteClass::First, self.head[0].take(), None)));
                }
                Ordering::Greater => {
                    return Ok(Some(Self::row(
                        SiteClass::Second,
                        None,
                        self.head[1].take(),
                    )));
                }
                Ordering::Equal => {}
            }
            // Alleles mode, one position in both files.
            let (chrom, pos) = (a.chrom, a.pos);
            let mut group: [Vec<Option<Site>>; 2] = [Vec::new(), Vec::new()];
            for (i, g) in group.iter_mut().enumerate() {
                loop {
                    self.fill(i)?;
                    match &self.head[i] {
                        Some(s) if s.chrom == chrom && s.pos == pos => g.push(self.head[i].take()),
                        _ => break,
                    }
                }
            }
            self.pair_group(group);
        }
    }

    fn pair_group(&mut self, [mut ga, mut gb]: [Vec<Option<Site>>; 2]) {
        let mut matched: Vec<Option<usize>> = vec![None; ga.len()];
        let mut used = vec![false; gb.len()];
        for (i, m) in matched.iter_mut().enumerate() {
            let a = ga[i].as_ref().expect("present");
            for (j, u) in used.iter_mut().enumerate() {
                if !*u && share_allele(a, gb[j].as_ref().expect("present")) {
                    *u = true;
                    *m = Some(j);
                    break;
                }
            }
        }
        // Leftovers at the position pair up in order as overlaps.
        let mut left_b = (0..gb.len()).filter(|&j| !used[j]);
        let mut overlap: Vec<Option<usize>> = vec![None; ga.len()];
        for (i, o) in overlap.iter_mut().enumerate() {
            if matched[i].is_none() {
                *o = left_b.next();
            }
        }
        let mut paired_b = vec![false; gb.len()];
        for i in 0..ga.len() {
            let (class, j) = match (matched[i], overlap[i]) {
                (Some(j), _) => (SiteClass::Both, Some(j)),
                (None, Some(j)) => (SiteClass::Overlap, Some(j)),
                (None, None) => (SiteClass::First, None),
            };
            let second = j.map(|j| {
                paired_b[j] = true;
                gb[j].take().expect("present")
            });
            self.out.push_back(Self::row(class, ga[i].take(), second));
        }
        for (j, done) in paired_b.iter().enumerate() {
            if !done {
                self.out
                    .push_back(Self::row(SiteClass::Second, None, gb[j].take()));
            }
        }
    }
}
