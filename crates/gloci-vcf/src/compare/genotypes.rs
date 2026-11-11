//! Genotype concordance between two paired records (decision D163).
//!
//! A genotype is compared as the unordered list of the alleles it names, an
//! allele being identified by what it is (its trimmed change, §4.2), not by
//! its index: `0/1` of `A>C,G` equals `0/2` of `A>G,C`. A genotype with any
//! missing allele is missing; vcftools counts `./1` as called (and stops on
//! it in the matrix), see `docs/differences.md#compare-half-missing`.

use super::walk::{MatchRule, Site, is_unknown_ref};
use gloci_core::trimmed_alleles;

/// One sample's genotype at one record, reduced to allele indices.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Gt {
    /// Not called: no GT, a masked genotype, or any `.` allele.
    Missing,
    Hap(u32),
    Dip(u32, u32),
    Poly(Box<[u32]>),
}

impl Gt {
    /// The genotype of `alleles` (`None` = `.`); no alleles at all is missing.
    pub fn from_alleles(alleles: &[Option<u32>]) -> Gt {
        match alleles {
            [] => Gt::Missing,
            [Some(a)] => Gt::Hap(*a),
            [Some(a), Some(b)] => Gt::Dip(*a, *b),
            _ if alleles.iter().any(Option::is_none) => Gt::Missing,
            _ => Gt::Poly(alleles.iter().flatten().copied().collect()),
        }
    }

    /// The largest allele index, `None` when missing.
    pub fn max_allele(&self) -> Option<u32> {
        match self {
            Gt::Missing => None,
            Gt::Hap(a) => Some(*a),
            Gt::Dip(a, b) => Some(*a.max(b)),
            Gt::Poly(v) => v.iter().copied().max(),
        }
    }
}

/// What comparing one sample's two genotypes gave.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GtCompare {
    BothMissing,
    FirstMissing,
    SecondMissing,
    Concordant,
    Discordant,
}

/// The matrix column or row of a genotype: `0/0`, `0/1`, `1/1`, `./.`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatrixCell {
    Cell(usize),
    /// Haploid or polyploid: not counted.
    Skip,
}

/// Counts of (file 2 genotype, file 1 genotype) over biallelic pairs with
/// the same REF and ALT: `cells[row][column]` with the row from file 2 and
/// the column from file 1, each in the order `0/0`, `0/1`, `1/1`, `./.`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DiscordanceMatrix {
    pub cells: [[u64; 4]; 4],
}

impl DiscordanceMatrix {
    pub fn add(&mut self, first: usize, second: usize) {
        self.cells[second][first] += 1;
    }
}

/// How two paired records' allele indices map to common allele ids (0 is the
/// reference allele). Reusable: [`AlleleMap::build`] overwrites it.
#[derive(Default)]
pub struct AlleleMap {
    identical: bool,
    a: Vec<u32>,
    b: Vec<u32>,
    keys: Vec<Vec<u8>>,
    /// REF and the ALT list are the same in both records (vcftools'
    /// MATCHING_ALLELES).
    pub matching: bool,
    /// Both records have one ALT and it is the same allele.
    pub biallelic_same: bool,
}

fn alt_iter(alts: &[u8]) -> impl Iterator<Item = &[u8]> {
    let none = alts.is_empty() || alts == b".";
    (!none)
        .then(|| alts.split(|&c| c == b','))
        .into_iter()
        .flatten()
}

/// Number of ALT alleles of a record.
pub fn alt_count(alts: &[u8]) -> usize {
    alt_iter(alts).count()
}

/// What an ALT is: its trimmed change when it is a sequence, else its text.
fn allele_key(out: &mut Vec<u8>, pos: u64, ref_: &[u8], alt: &[u8]) {
    out.clear();
    if alt
        .iter()
        .any(|c| matches!(c, b'<' | b'[' | b']' | b'*' | b'.'))
    {
        out.push(b'S');
        out.extend_from_slice(alt);
        return;
    }
    let (p, r, a) = trimmed_alleles(pos.saturating_sub(1), ref_, alt);
    out.push(b'Q');
    out.extend_from_slice(&p.to_le_bytes());
    out.extend_from_slice(r);
    out.push(0);
    out.extend_from_slice(a);
}

impl AlleleMap {
    /// Map the alleles of `a` and `b`. In `--match position` a REF of `N` or
    /// `.` takes the other record's REF, as vcftools does.
    pub fn build(&mut self, a: &Site, b: &Site, rule: MatchRule) {
        let (mut ra, mut rb) = (a.ref_.as_slice(), b.ref_.as_slice());
        if rule == MatchRule::Position && a.pos == b.pos {
            if is_unknown_ref(ra) && !is_unknown_ref(rb) {
                ra = rb;
            } else if is_unknown_ref(rb) && !is_unknown_ref(ra) {
                rb = ra;
            }
        }
        self.matching = ra == rb && a.alts == b.alts;
        self.identical = self.matching;
        let (na, nb) = (alt_count(&a.alts), alt_count(&b.alts));
        if self.identical {
            self.biallelic_same = na == 1;
            return;
        }
        self.a.clear();
        self.b.clear();
        self.keys.clear();
        let mut key = Vec::new();
        for (site, r, out) in [(a, ra, 0usize), (b, rb, 1)] {
            let map = if out == 0 { &mut self.a } else { &mut self.b };
            map.push(0);
            for alt in alt_iter(&site.alts) {
                allele_key(&mut key, site.pos, r, alt);
                let id = match self.keys.iter().position(|k| *k == key) {
                    Some(i) => i,
                    None => {
                        self.keys.push(key.clone());
                        self.keys.len() - 1
                    }
                };
                map.push(id as u32 + 1);
            }
        }
        self.biallelic_same = na == 1 && nb == 1 && self.a[1] == self.b[1];
    }

    #[inline]
    fn id_a(&self, i: u32) -> u32 {
        if self.identical {
            i
        } else {
            self.a[i as usize]
        }
    }

    #[inline]
    fn id_b(&self, i: u32) -> u32 {
        if self.identical {
            i
        } else {
            self.b[i as usize]
        }
    }

    /// Compare the genotypes of one sample at the paired records.
    pub fn compare(&self, a: &Gt, b: &Gt) -> GtCompare {
        match (a, b) {
            (Gt::Missing, Gt::Missing) => GtCompare::BothMissing,
            (Gt::Missing, _) => GtCompare::FirstMissing,
            (_, Gt::Missing) => GtCompare::SecondMissing,
            (Gt::Dip(x, y), Gt::Dip(u, v)) => {
                let (x, y) = (self.id_a(*x), self.id_a(*y));
                let (u, v) = (self.id_b(*u), self.id_b(*v));
                if (x.min(y), x.max(y)) == (u.min(v), u.max(v)) {
                    GtCompare::Concordant
                } else {
                    GtCompare::Discordant
                }
            }
            _ => {
                let mut ia = self.ids(a, true);
                let mut ib = self.ids(b, false);
                ia.sort_unstable();
                ib.sort_unstable();
                if ia == ib {
                    GtCompare::Concordant
                } else {
                    GtCompare::Discordant
                }
            }
        }
    }

    fn ids(&self, g: &Gt, first: bool) -> Vec<u32> {
        let f = |i: u32| if first { self.id_a(i) } else { self.id_b(i) };
        match g {
            Gt::Missing => Vec::new(),
            Gt::Hap(x) => vec![f(*x)],
            Gt::Dip(x, y) => vec![f(*x), f(*y)],
            Gt::Poly(v) => v.iter().map(|&i| f(i)).collect(),
        }
    }

    /// The matrix cell of a genotype of a record whose ALT is the one
    /// biallelic ALT (see [`AlleleMap::biallelic_same`]).
    pub fn cell(g: &Gt) -> MatrixCell {
        match g {
            Gt::Missing => MatrixCell::Cell(3),
            Gt::Dip(x, y) => MatrixCell::Cell((*x as usize).min(1) + (*y as usize).min(1)),
            Gt::Hap(_) | Gt::Poly(_) => MatrixCell::Skip,
        }
    }
}
