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
    FirstAppearance,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mode {
    Lexicographic,
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
            Mode::ById => a.0.cmp(&b.0),
        }
    }

    /// Dense rank of each id in the table's order.
    pub fn ranks(&self) -> Vec<u32> {
        let n = self.names.len();
        match self.mode {
            Mode::ById => (0..n as u32).collect(),
            Mode::Lexicographic => {
                let mut ids: Vec<u32> = (0..n as u32).collect();
                ids.sort_unstable_by(|&a, &b| self.names[a as usize].cmp(&self.names[b as usize]));
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
