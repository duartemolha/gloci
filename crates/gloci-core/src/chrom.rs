//! Chromosome interning and ordering.

use crate::error::UnknownChrom;
use crate::genome::Genome;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::path::PathBuf;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ChromId(pub u32);

/// Where a chromosome order came from.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OrderSource {
    GenomeFile(PathBuf),
    VcfHeader(String),
    Lexicographic,
    Karyotype,
    FirstAppearance,
}

/// The order used when neither `-g` nor a VCF/BCF header sets one (spec §4.5).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ChromOrder {
    #[default]
    Karyotype,
    Lex,
}

impl std::str::FromStr for ChromOrder {
    type Err = String;
    fn from_str(s: &str) -> Result<ChromOrder, String> {
        match s {
            "karyotype" => Ok(ChromOrder::Karyotype),
            "lex" => Ok(ChromOrder::Lex),
            _ => Err(format!(
                "unknown chromosome order \"{s}\" (expected karyotype or lex)"
            )),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mode {
    Lexicographic,
    Karyotype,
    /// Order by id (first appearance, or fixed genome order).
    ById,
}

#[derive(Debug, Clone)]
pub struct ChromTable {
    map: HashMap<Box<[u8]>, ChromId>,
    names: Vec<Box<[u8]>>,
    mode: Mode,
    closed: bool,
    source: OrderSource,
}

impl ChromTable {
    fn empty(mode: Mode, closed: bool, source: OrderSource) -> ChromTable {
        ChromTable {
            map: HashMap::new(),
            names: Vec::new(),
            mode,
            closed,
            source,
        }
    }

    /// Open table ordered by name bytes.
    pub fn lexicographic() -> ChromTable {
        Self::empty(Mode::Lexicographic, false, OrderSource::Lexicographic)
    }

    /// Open table ordered by [`karyotype_cmp`].
    pub fn karyotype() -> ChromTable {
        Self::empty(Mode::Karyotype, false, OrderSource::Karyotype)
    }

    /// Open table ordered by first appearance.
    pub fn first_appearance() -> ChromTable {
        Self::empty(Mode::ById, false, OrderSource::FirstAppearance)
    }

    /// Closed table: genome names are ids `0..n` in genome order; unknown names are rejected.
    pub fn closed(genome: &Genome, source: OrderSource) -> ChromTable {
        let mut t = Self::empty(Mode::ById, true, source);
        for i in 0..genome.len() {
            t.insert(genome.name(i));
        }
        t
    }

    fn insert(&mut self, name: &[u8]) -> ChromId {
        let id = ChromId(self.names.len() as u32);
        let boxed: Box<[u8]> = name.into();
        self.names.push(boxed.clone());
        self.map.insert(boxed, id);
        id
    }

    pub fn intern(&mut self, name: &[u8]) -> Result<ChromId, UnknownChrom> {
        if let Some(&id) = self.map.get(name) {
            return Ok(id);
        }
        if self.closed {
            return Err(UnknownChrom(String::from_utf8_lossy(name).into_owned()));
        }
        Ok(self.insert(name))
    }

    pub fn get(&self, name: &[u8]) -> Option<ChromId> {
        self.map.get(name).copied()
    }

    pub fn name(&self, id: ChromId) -> &[u8] {
        &self.names[id.0 as usize]
    }

    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    pub fn cmp(&self, a: ChromId, b: ChromId) -> Ordering {
        match self.mode {
            Mode::Lexicographic => self.name(a).cmp(self.name(b)),
            Mode::Karyotype => karyotype_cmp(self.name(a), self.name(b)),
            Mode::ById => a.0.cmp(&b.0),
        }
    }

    /// Compare two names in the table's order without interning them.
    /// `None` when the order is by id and either name is not in the table
    /// (a genome or VCF header that lacks it, or a first-appearance table
    /// that has not seen it).
    pub fn cmp_names(&self, a: &[u8], b: &[u8]) -> Option<Ordering> {
        match self.mode {
            Mode::Lexicographic => Some(a.cmp(b)),
            Mode::Karyotype => Some(karyotype_cmp(a, b)),
            Mode::ById => Some(self.get(a)?.cmp(&self.get(b)?)),
        }
    }

    /// Dense rank of each id in the table's order.
    pub fn ranks(&self) -> Vec<u32> {
        let n = self.names.len();
        match self.mode {
            Mode::ById => (0..n as u32).collect(),
            Mode::Lexicographic | Mode::Karyotype => {
                let mut ids: Vec<u32> = (0..n as u32).collect();
                ids.sort_unstable_by(|&a, &b| self.cmp(ChromId(a), ChromId(b)));
                let mut ranks = vec![0u32; n];
                for (rank, id) in ids.into_iter().enumerate() {
                    ranks[id as usize] = rank as u32;
                }
                ranks
            }
        }
    }

    pub fn order_source(&self) -> &OrderSource {
        &self.source
    }
}

/// Compare chromosome names in karyotype order (spec §4.5 rule 3): a total
/// order that is a pure function of the two names, so the relative order of
/// two chromosomes never changes as other names are interned.
pub fn karyotype_cmp(a: &[u8], b: &[u8]) -> Ordering {
    let (ra, rb) = (strip_chr(a), strip_chr(b));
    let (ka, kb) = (classify(ra), classify(rb));
    ka.cmp(&kb)
        .then_with(|| match ka {
            Class::Numbered => cmp_numbered(ra, rb),
            Class::Other => natural_cmp(ra, rb),
            _ => Ordering::Equal,
        })
        .then_with(|| a.cmp(b))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Class {
    Numbered,
    X,
    Y,
    Mito,
    Other,
}

/// Strip one leading `chr` (any case) if something remains.
fn strip_chr(name: &[u8]) -> &[u8] {
    if name.len() > 3 && name[..3].eq_ignore_ascii_case(b"chr") {
        &name[3..]
    } else {
        name
    }
}

fn digit_len(s: &[u8]) -> usize {
    s.iter().take_while(|b| b.is_ascii_digit()).count()
}

fn classify(r: &[u8]) -> Class {
    let d = digit_len(r);
    if d > 0 && r[d..].iter().all(|b| b.is_ascii_alphabetic()) {
        return Class::Numbered;
    }
    match r {
        b"X" => Class::X,
        b"Y" => Class::Y,
        b"M" | b"MT" => Class::Mito,
        _ => Class::Other,
    }
}

/// Numeric comparison of two digit runs of any length.
fn cmp_digits(a: &[u8], b: &[u8]) -> Ordering {
    let strip = |s: &[u8]| -> usize { s.iter().take_while(|&&c| c == b'0').count() };
    let (a, b) = (&a[strip(a)..], &b[strip(b)..]);
    a.len().cmp(&b.len()).then_with(|| a.cmp(b))
}

fn cmp_numbered(a: &[u8], b: &[u8]) -> Ordering {
    let (da, db) = (digit_len(a), digit_len(b));
    cmp_digits(&a[..da], &b[..db]).then_with(|| a[da..].cmp(&b[db..]))
}

/// Split into maximal digit / non-digit runs and compare run by run. A digit
/// run sorts before a non-digit run; a sequence that is a prefix of the other
/// comes first.
fn natural_cmp(a: &[u8], b: &[u8]) -> Ordering {
    let (mut a, mut b) = (a, b);
    loop {
        match (run(&mut a), run(&mut b)) {
            (None, None) => return Ordering::Equal,
            (None, Some(_)) => return Ordering::Less,
            (Some(_), None) => return Ordering::Greater,
            (Some((da, ra)), Some((db, rb))) => {
                let o = match (da, db) {
                    (true, true) => cmp_digits(ra, rb),
                    (false, false) => ra.cmp(rb),
                    (true, false) => Ordering::Less,
                    (false, true) => Ordering::Greater,
                };
                if o != Ordering::Equal {
                    return o;
                }
            }
        }
    }
}

fn run<'a>(s: &mut &'a [u8]) -> Option<(bool, &'a [u8])> {
    let digit = s.first()?.is_ascii_digit();
    let n = s.iter().take_while(|c| c.is_ascii_digit() == digit).count();
    let (r, rest) = s.split_at(n);
    *s = rest;
    Some((digit, r))
}
