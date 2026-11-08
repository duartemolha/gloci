//! Position and ID list files.

use gloci_io::{IoError, Loc};
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader};

/// `(chromosome, 1-based POS)` pairs, indexed by chromosome so a lookup
/// needs no allocation.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PositionSet(HashMap<Vec<u8>, HashSet<u64>>);

impl PositionSet {
    pub fn insert(&mut self, chrom: Vec<u8>, pos: u64) {
        self.0.entry(chrom).or_default().insert(pos);
    }

    pub fn contains(&self, chrom: &[u8], pos: u64) -> bool {
        self.0.get(chrom).is_some_and(|s| s.contains(&pos))
    }

    pub fn len(&self) -> usize {
        self.0.values().map(HashSet::len).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

impl FromIterator<(Vec<u8>, u64)> for PositionSet {
    fn from_iter<I: IntoIterator<Item = (Vec<u8>, u64)>>(iter: I) -> Self {
        let mut s = PositionSet::default();
        for (c, p) in iter {
            s.insert(c, p);
        }
        s
    }
}

/// Visit the data lines of a list file: `#` comments and blank lines are
/// skipped, a trailing `\r` is removed.
pub(crate) fn each_line(
    path: &str,
    mut f: impl FnMut(u64, &[u8]) -> Result<(), String>,
) -> Result<(), IoError> {
    let file = std::fs::File::open(path).map_err(|source| IoError::Open {
        path: path.into(),
        source,
    })?;
    let mut r = BufReader::new(file);
    let mut buf = Vec::new();
    let mut n = 0u64;
    loop {
        buf.clear();
        let read = r
            .read_until(b'\n', &mut buf)
            .map_err(|source| IoError::Read {
                path: path.into(),
                source,
            })?;
        if read == 0 {
            return Ok(());
        }
        n += 1;
        let mut line = &buf[..];
        while let [rest @ .., b'\n' | b'\r'] = line {
            line = rest;
        }
        if line.is_empty() || line[0] == b'#' {
            continue;
        }
        f(n, line).map_err(|msg| IoError::Parse {
            path: path.into(),
            loc: Loc::Line(n),
            msg,
        })?;
    }
}

/// Read `chrom<TAB>POS` lines (1-based POS); further columns are ignored.
pub fn read_positions(path: &str) -> Result<PositionSet, IoError> {
    let mut out = PositionSet::default();
    each_line(path, |_, line| {
        let mut cols = line.split(|&b| b == b'\t');
        let chrom = cols.next().unwrap_or_default();
        let Some(pos) = cols.next() else {
            return Err(format!(
                "expected \"chrom<TAB>POS\", found {:?}; separate the two columns with a tab",
                String::from_utf8_lossy(line)
            ));
        };
        let shown = || String::from_utf8_lossy(pos).into_owned();
        let pos = std::str::from_utf8(pos)
            .ok()
            .and_then(|s| s.parse::<u64>().ok())
            .ok_or_else(|| format!("POS {:?} is not an integer", shown()))?;
        if pos == 0 {
            return Err("POS 0 is not valid; positions are 1-based".into());
        }
        out.insert(chrom.to_vec(), pos);
        Ok(())
    })?;
    Ok(out)
}

/// Read one ID per line.
pub fn read_ids(path: &str) -> Result<HashSet<Vec<u8>>, IoError> {
    let mut out = HashSet::new();
    each_line(path, |_, line| {
        out.insert(line.to_vec());
        Ok(())
    })?;
    Ok(out)
}
