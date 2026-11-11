//! `gloci maskfasta`: replace the bases under intervals, streaming the FASTA
//! and keeping its headers and line layout exactly (D113).

use gloci_core::{ChromTable, thousands};
use gloci_io::{FastaLine, FastaStream, IoError, RecordBatch, RecordReader, next_batch_size};
use std::collections::{HashMap, HashSet};
use std::io::Write;

/// How a masked base is written.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mask {
    /// Replace the base with this character.
    Hard(u8),
    /// Lowercase the base if it is a letter; other characters are kept.
    Soft,
}

/// The mask intervals of every chromosome, 0-based half-open, sorted and
/// merged so that a base is masked (and counted) once.
#[derive(Debug, Default)]
pub struct MaskSet {
    by_name: HashMap<Vec<u8>, Vec<(u64, u64)>>,
    /// Names in order of first appearance, for the warning about absent ones.
    order: Vec<Vec<u8>>,
}

impl MaskSet {
    /// Read every record of `readers` (any order, any supported interval
    /// format) and merge the intervals per chromosome. Zero-length intervals
    /// mask nothing and are dropped.
    pub fn load(
        readers: Vec<Box<dyn RecordReader>>,
        chroms: &mut ChromTable,
    ) -> Result<MaskSet, IoError> {
        let mut per: Vec<Vec<(u64, u64)>> = Vec::new();
        let mut seen: Vec<u32> = Vec::new();
        let mut per_seen: Vec<bool> = Vec::new();
        let mut batch = RecordBatch::default();
        for mut r in readers {
            let mut max = 1024;
            loop {
                batch.clear();
                let more = r.read_batch(chroms, &mut batch, max)?;
                for m in &batch.recs {
                    let id = m.chrom.0 as usize;
                    if per.len() <= id {
                        per.resize_with(id + 1, Vec::new);
                    }
                    if per_seen.len() <= id {
                        per_seen.resize(id + 1, false);
                    }
                    if !per_seen[id] {
                        per_seen[id] = true;
                        seen.push(m.chrom.0);
                    }
                    if !m.span.is_point() {
                        per[id].push((m.span.start, m.span.end));
                    }
                }
                max = next_batch_size(&batch, max);
                if !more {
                    break;
                }
            }
            crate::skipped::log_reader(r.as_ref());
        }
        let mut set = MaskSet::default();
        for id in seen {
            let mut v = std::mem::take(&mut per[id as usize]);
            v.sort_unstable();
            let mut merged: Vec<(u64, u64)> = Vec::with_capacity(v.len());
            for (s, e) in v {
                match merged.last_mut() {
                    Some(last) if s <= last.1 => last.1 = last.1.max(e),
                    _ => merged.push((s, e)),
                }
            }
            let name = chroms.name(gloci_core::ChromId(id)).to_vec();
            set.order.push(name.clone());
            set.by_name.entry(name).or_default().extend(merged);
        }
        Ok(set)
    }

    /// The merged intervals of chromosome `name`.
    pub fn get(&self, name: &[u8]) -> Option<&[(u64, u64)]> {
        self.by_name.get(name).map(|v| &v[..])
    }

    /// The bases the set covers, each counted once.
    pub fn total_bases(&self) -> u64 {
        self.by_name.values().flatten().map(|(s, e)| e - s).sum()
    }
}

/// What a run did.
#[derive(Debug, Default, PartialEq, Eq)]
pub struct MaskSummary {
    /// Bases masked, each counted once.
    pub masked_bases: u64,
    /// Chromosomes that have masks but are not in the FASTA, in the order
    /// the masks named them.
    pub missing_chroms: Vec<Vec<u8>>,
    /// With [`PastEnd::Clip`]: merged mask intervals that reached past the end
    /// of their sequence (and were clipped to it).
    pub past_end: u64,
    /// The first of them, described for a warning.
    pub first_past_end: Option<String>,
}

/// What a mask that reaches past the end of its sequence does.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PastEnd {
    /// A data error, reported once the sequence is done.
    Error,
    /// The mask is clipped to the sequence (one that starts at or past the
    /// end masks nothing) and counted in the summary.
    Clip,
}

/// The first word of a header: the name masks refer to.
fn first_word(header: &[u8]) -> &[u8] {
    let end = header
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .unwrap_or(header.len());
    &header[..end]
}

/// Mask the part of `buf` (the bases at `pos..pos + buf.len()` of its
/// sequence) that `masks` cover; `cursor` is the first mask that can still
/// reach later bases. Returns the number of bases masked.
fn apply(buf: &mut [u8], pos: u64, masks: &[(u64, u64)], cursor: &mut usize, mask: Mask) -> u64 {
    let end = pos + buf.len() as u64;
    while *cursor < masks.len() && masks[*cursor].1 <= pos {
        *cursor += 1;
    }
    let mut n = 0;
    let mut i = *cursor;
    while i < masks.len() && masks[i].0 < end {
        let s = masks[i].0.max(pos);
        let e = masks[i].1.min(end);
        let part = &mut buf[(s - pos) as usize..(e - pos) as usize];
        match mask {
            Mask::Hard(c) => part.fill(c),
            Mask::Soft => part.make_ascii_lowercase(),
        }
        n += e - s;
        i += 1;
    }
    n
}

/// The sequence being written: its name, masks and progress.
struct Current<'a> {
    name: Vec<u8>,
    masks: Option<&'a [(u64, u64)]>,
    cursor: usize,
    len: u64,
}

impl Current<'_> {
    /// The merged masks that reach past the end of the sequence, now that
    /// its length is known: how many, and the first.
    fn past_end(&self) -> (u64, Option<(u64, u64)>) {
        let Some(m) = self.masks else {
            return (0, None);
        };
        // Merged masks are disjoint and sorted, so their ends increase.
        let i = m.partition_point(|&(_, e)| e <= self.len);
        ((m.len() - i) as u64, m.get(i).copied())
    }

    /// The sequence is complete: settle masks that reach past its end.
    fn finish(
        &self,
        fasta: &str,
        policy: PastEnd,
        summary: &mut MaskSummary,
    ) -> Result<(), IoError> {
        let (n, first) = self.past_end();
        let Some((s, e)) = first else {
            return Ok(());
        };
        let name = String::from_utf8_lossy(&self.name);
        let what = format!(
            "a mask on {name} covers {}-{} (0-based, half-open) but {name} in {fasta} has {} bases",
            thousands(s),
            thousands(e),
            thousands(self.len),
        );
        match policy {
            PastEnd::Error => Err(IoError::Detect {
                path: fasta.to_string(),
                msg: format!(
                    "{what}; check that the masks belong to this FASTA (same assembly and \
chromosome names, or rename with --chrom-map), or pass --skip-invalid to clip such masks to the \
sequence end"
                ),
            }),
            PastEnd::Clip => {
                summary.past_end += n;
                summary.first_past_end.get_or_insert(what);
                Ok(())
            }
        }
    }
}

fn write_err(source: std::io::Error) -> IoError {
    IoError::Write {
        path: "output".into(),
        source,
    }
}

/// Copy `fasta` to `out`, masking the bases the intervals of `masks` cover.
/// Headers, line widths, CRLF line ends and a missing final line break are
/// kept; a piece of at most 1 MiB is held at a time. A mask that reaches past
/// the end of its sequence is an error, reported once that sequence is done.
pub fn maskfasta(
    fasta: FastaStream,
    masks: &MaskSet,
    mask: Mask,
    out: &mut dyn Write,
) -> Result<MaskSummary, IoError> {
    maskfasta_with(fasta, masks, mask, PastEnd::Error, out)
}

/// [`maskfasta`] with a choice of what a mask past its sequence's end does.
pub fn maskfasta_with(
    mut fasta: FastaStream,
    masks: &MaskSet,
    mask: Mask,
    past_end: PastEnd,
    out: &mut dyn Write,
) -> Result<MaskSummary, IoError> {
    let name = fasta.name().to_string();
    let mut summary = MaskSummary::default();
    let mut cur: Option<Current> = None;
    let mut present: HashSet<Vec<u8>> = HashSet::new();
    let mut piece: Vec<u8> = Vec::new();
    loop {
        let is_header = match fasta.next()? {
            None => break,
            Some(FastaLine::Header(h)) => {
                piece.clear();
                piece.extend_from_slice(h);
                true
            }
            Some(FastaLine::Seq(s)) => {
                piece.clear();
                piece.extend_from_slice(s);
                false
            }
        };
        let (crlf, ended) = (fasta.crlf(), fasta.line_ended());
        if is_header {
            if let Some(c) = &cur {
                c.finish(&name, past_end, &mut summary)?;
            }
            let key = first_word(&piece);
            present.insert(key.to_vec());
            cur = Some(Current {
                name: key.to_vec(),
                masks: masks.get(key),
                cursor: 0,
                len: 0,
            });
            out.write_all(b">").map_err(write_err)?;
        } else if let Some(c) = &mut cur {
            if let Some(m) = c.masks {
                summary.masked_bases += apply(&mut piece, c.len, m, &mut c.cursor, mask);
            }
            c.len += piece.len() as u64;
        }
        out.write_all(&piece).map_err(write_err)?;
        if crlf {
            out.write_all(b"\r").map_err(write_err)?;
        }
        if ended {
            out.write_all(b"\n").map_err(write_err)?;
        }
    }
    if let Some(c) = &cur {
        c.finish(&name, past_end, &mut summary)?;
    }
    summary.missing_chroms = masks
        .order
        .iter()
        .filter(|n| !present.contains(n.as_slice()))
        .cloned()
        .collect();
    Ok(summary)
}
