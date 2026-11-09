//! Genome (chromosome names and optional lengths) and its text parsers.

use crate::error::GenomeError;
use crate::interval::parse_coord;
use std::collections::HashSet;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Genome {
    entries: Vec<(Vec<u8>, Option<u64>)>,
}

impl Genome {
    pub fn new(entries: Vec<(Vec<u8>, Option<u64>)>) -> Result<Genome, GenomeError> {
        let mut seen = HashSet::new();
        for (name, _) in &entries {
            if !seen.insert(name.as_slice()) {
                return Err(GenomeError::Duplicate(
                    String::from_utf8_lossy(name).into_owned(),
                ));
            }
        }
        Ok(Genome { entries })
    }

    /// `.fai` / `.genome`: `name TAB length [TAB ...]` per line.
    pub fn parse_sizes(text: &[u8]) -> Result<Genome, GenomeError> {
        let mut entries = Vec::new();
        for (i, line) in lines(text) {
            if line.starts_with(b"#") {
                continue;
            }
            let mut f = line.split(|&b| b == b'\t');
            let name = f.next().unwrap_or_default();
            if name.is_empty() {
                return Err(GenomeError::Line {
                    line: i,
                    msg: "empty chromosome name".into(),
                });
            }
            let len = f.next().ok_or_else(|| GenomeError::Line {
                line: i,
                msg: "expected name and length separated by a TAB".into(),
            })?;
            let len = parse_coord(len).map_err(|e| GenomeError::Line {
                line: i,
                msg: e.to_string(),
            })?;
            entries.push((name.to_vec(), Some(len)));
        }
        Genome::new(entries)
    }

    /// SAM `.dict`: `@SQ SN:name LN:length` lines; everything else is ignored.
    pub fn parse_dict(text: &[u8]) -> Result<Genome, GenomeError> {
        let mut entries = Vec::new();
        for (i, line) in lines(text) {
            if !line.starts_with(b"@SQ\t") {
                continue;
            }
            let (mut name, mut len) = (None, None);
            for f in line.split(|&b| b == b'\t').skip(1) {
                if let Some(v) = f.strip_prefix(b"SN:") {
                    name = Some(v.to_vec());
                } else if let Some(v) = f.strip_prefix(b"LN:") {
                    len = Some(parse_coord(v).map_err(|e| GenomeError::Line {
                        line: i,
                        msg: e.to_string(),
                    })?);
                }
            }
            let name = name.ok_or_else(|| GenomeError::Line {
                line: i,
                msg: "@SQ line without SN:".into(),
            })?;
            entries.push((name, len));
        }
        Genome::new(entries)
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn name(&self, i: usize) -> &[u8] {
        &self.entries[i].0
    }

    pub fn length(&self, i: usize) -> Option<u64> {
        self.entries[i].1
    }
}

/// Non-empty lines with 1-based numbers, CR stripped.
fn lines(text: &[u8]) -> impl Iterator<Item = (usize, &[u8])> {
    text.split(|&b| b == b'\n')
        .enumerate()
        .map(|(i, l)| (i + 1, l.strip_suffix(b"\r").unwrap_or(l)))
        .filter(|(_, l)| !l.is_empty())
}
