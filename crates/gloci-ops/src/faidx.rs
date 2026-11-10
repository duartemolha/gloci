//! `gloci faidx`: write the bases of regions of an indexed FASTA (replaces
//! `samtools faidx`, D131).
//!
//! Each region is streamed in the chunks [`IndexedFasta`] delivers (at most
//! [`gloci_io::FETCH_CHUNK`] bases), re-wrapped to the output line width with
//! the column carried across chunks, so memory does not depend on the region
//! length. A reverse-complemented region is read last chunk first and each
//! chunk is reverse-complemented as it goes.

use crate::seq::revcomp_in_place;
use gloci_core::Interval;
use gloci_io::{IndexedFasta, IoError};
use std::io::Write;
use std::str::FromStr;

/// The strand indicator appended to a header (`--mark-strand`).
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum MarkStrand {
    /// `/rc` on reverse-complemented regions, nothing otherwise.
    #[default]
    Rc,
    /// Nothing.
    No,
    /// `(+)` or `(-)`.
    Sign,
    /// The given text for the forward and the reverse strand.
    Custom(String, String),
}

impl FromStr for MarkStrand {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "rc" => Ok(MarkStrand::Rc),
            "no" => Ok(MarkStrand::No),
            "sign" => Ok(MarkStrand::Sign),
            _ => match s.strip_prefix("custom,").map(|r| r.split_once(',')) {
                Some(Some((pos, neg))) => Ok(MarkStrand::Custom(pos.into(), neg.into())),
                _ => Err(format!(
                    "unknown strand mark '{s}'; use rc, no, sign or custom,POS,NEG (for example custom,+,-)"
                )),
            },
        }
    }
}

impl MarkStrand {
    fn suffix(&self, reversed: bool) -> &str {
        match (self, reversed) {
            (MarkStrand::Rc, true) => "/rc",
            (MarkStrand::Sign, false) => "(+)",
            (MarkStrand::Sign, true) => "(-)",
            (MarkStrand::Custom(p, _), false) => p,
            (MarkStrand::Custom(_, n), true) => n,
            _ => "",
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct FaidxOptions {
    /// Bases per output line; `None` is each sequence's `.fai` line width and
    /// `Some(0)` is one line.
    pub line_width: Option<u64>,
    pub reverse_complement: bool,
    pub mark: MarkStrand,
}

/// One region to write. `iv` is 0-based half-open; an end of `u64::MAX` means
/// "to the end of the sequence". An end past the sequence is cut there with a
/// warning, and a start at or past the end gives the header alone.
#[derive(Debug, Clone)]
pub struct FaidxRegion {
    /// The header text: the region as the user typed it.
    pub label: String,
    pub chrom: Vec<u8>,
    pub iv: Interval,
}

/// Bases reverse-complemented at a time.
const RC_BLOCK: usize = 8192;

fn write_err(source: std::io::Error) -> IoError {
    IoError::Write {
        path: "output".into(),
        source,
    }
}

/// Writes bytes in lines of `width` (0: no wrapping), keeping the column
/// between calls.
struct Wrapper<'a> {
    out: &'a mut dyn Write,
    width: usize,
    col: usize,
}

impl Wrapper<'_> {
    fn put(&mut self, mut data: &[u8]) -> Result<(), IoError> {
        if self.width == 0 {
            self.col += data.len();
            return self.out.write_all(data).map_err(write_err);
        }
        while !data.is_empty() {
            let n = (self.width - self.col).min(data.len());
            self.out.write_all(&data[..n]).map_err(write_err)?;
            self.col += n;
            data = &data[n..];
            if self.col == self.width {
                self.out.write_all(b"\n").map_err(write_err)?;
                self.col = 0;
            }
        }
        Ok(())
    }

    fn finish(&mut self) -> Result<(), IoError> {
        if self.col > 0 {
            self.out.write_all(b"\n").map_err(write_err)?;
            self.col = 0;
        }
        Ok(())
    }
}

/// Write every region, in order.
pub fn faidx(
    fasta: &mut IndexedFasta,
    regions: &[FaidxRegion],
    o: &FaidxOptions,
    out: &mut dyn Write,
) -> Result<(), IoError> {
    let mut rc_buf: Vec<u8> = Vec::new();
    for r in regions {
        let shown = String::from_utf8_lossy(&r.chrom).into_owned();
        let Some(entry) = fasta.entry_of(&r.chrom) else {
            return Err(IoError::Detect {
                path: fasta.name().to_string(),
                msg: format!(
                    "{}: sequence {shown} (region {}) is not in the FASTA index; check the sequence name",
                    fasta.name(),
                    r.label
                ),
            });
        };
        let (len, line_bases) = (entry.len, entry.line_bases);
        let width = o.line_width.unwrap_or(line_bases);
        let width = usize::try_from(width).unwrap_or(usize::MAX);
        let mut iv = r.iv;
        if iv.start >= len || iv.start >= iv.end {
            if iv.start >= len {
                log::warn!(
                    "{} starts past the {len} bases of {shown}; writing its header with no bases",
                    r.label
                );
            }
            iv.end = iv.start;
        } else if iv.end > len {
            if iv.end != u64::MAX {
                log::warn!(
                    "{} ends past the {len} bases of {shown}; cut at the end of {shown}",
                    r.label
                );
            }
            iv.end = len;
        }
        writeln!(out, ">{}{}", r.label, o.mark.suffix(o.reverse_complement)).map_err(write_err)?;
        if iv.start >= iv.end {
            continue;
        }
        let mut w = Wrapper {
            out: &mut *out,
            width,
            col: 0,
        };
        if o.reverse_complement {
            // Each chunk is read back to front in small blocks, so the only
            // large buffer is the fetch chunk itself.
            fasta.fetch_rev(&r.chrom, iv, &mut |chunk| {
                for block in chunk.rchunks(RC_BLOCK) {
                    rc_buf.clear();
                    rc_buf.extend_from_slice(block);
                    revcomp_in_place(&mut rc_buf, false);
                    w.put(&rc_buf)?;
                }
                Ok(())
            })?;
        } else {
            fasta.fetch(&r.chrom, iv, &mut |chunk| w.put(chunk))?;
        }
        w.finish()?;
    }
    Ok(())
}
