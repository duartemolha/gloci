//! An alignment as BED text (D106, bedtools `bamtobed`): BED6, BED12, or one BED6 row per
//! block. The same BED6 line is what the commands that print an alignment's columns write
//! (D104), so it lives here and not in `convert`.

use super::{AlnHeader, AlnRef, aln_ref, blocks_within};
use crate::detect::Format;
use crate::error::{IoError, Loc};
use crate::record::{Header, RecordReader};
use gloci_core::{Interval, Strand};
use std::sync::Arc;

/// Which rows an alignment becomes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AlnBedLayout {
    /// One row: `RNAME start end QNAME score strand`.
    #[default]
    Bed6,
    /// One BED12 row, its blocks cut at `N` (and at `D` with `split_deletions`).
    Bed12,
    /// One BED6 row per block (bedtools `-split`).
    Blocks,
}

#[derive(Debug, Clone)]
pub struct AlnBedOptions {
    pub layout: AlnBedLayout,
    /// The integer aux tag the score column takes; `None` is the MAPQ.
    pub score_tag: Option<[u8; 2]>,
    /// Append the CIGAR as a seventh column (BED6 only).
    pub cigar: bool,
    /// Cut blocks at `D` too (D105).
    pub split_deletions: bool,
    /// The BED12 `itemRgb`.
    pub color: String,
}

impl Default for AlnBedOptions {
    fn default() -> Self {
        AlnBedOptions {
            layout: AlnBedLayout::Bed6,
            score_tag: None,
            cigar: false,
            split_deletions: false,
            color: "255,0,0".into(),
        }
    }
}

/// The record layout of an alignment reader's raw bytes, or `None` for a reader of anything else.
fn raw_layout(r: &dyn RecordReader) -> Option<(Format, Arc<AlnHeader>)> {
    let Header::Alignment(h) = r.header() else {
        return None;
    };
    let format = match r.format() {
        Format::Bam | Format::Cram => Format::Bam,
        Format::Sam => Format::Sam,
        _ => return None,
    };
    Some((format, h.clone()))
}

/// The reference blocks of the alignments of one input (D105), for the commands that compare
/// blocks (`--split`, `--split-deletions`): the segments between `N` operations (and `D` too
/// with `split_deletions`), cut at the end of the reference like the span, never empty.
pub struct AlnBlocks {
    format: Format,
    name: String,
    hdr: Arc<AlnHeader>,
    split_deletions: bool,
}

impl AlnBlocks {
    /// Blocks for the records of `r`, or `None` when `r` does not read alignments.
    pub fn for_reader(r: &dyn RecordReader, split_deletions: bool) -> Option<AlnBlocks> {
        let (format, hdr) = raw_layout(r)?;
        Some(AlnBlocks {
            format,
            name: r.name().to_string(),
            hdr,
            split_deletions,
        })
    }

    /// The blocks of the record `raw` (record or line number `loc`, for errors) whose span, as
    /// the reader gave it, is `span`.
    pub fn get(&self, raw: &[u8], loc: u64, span: Interval) -> Result<Vec<Interval>, IoError> {
        let a = aln_ref(self.format, raw, &self.hdr).map_err(|msg| IoError::Parse {
            path: self.name.clone(),
            loc: if self.format == Format::Sam {
                Loc::Line(loc)
            } else {
                Loc::Record(loc)
            },
            msg,
        })?;
        let ref_len = a.ref_id().map_or(u64::MAX, |i| self.hdr.refs[i].1);
        let mut out = Vec::new();
        blocks_within(&a, self.split_deletions, ref_len, &mut out);
        if out.is_empty() {
            out.push(span);
        }
        Ok(out)
    }
}

/// Writes alignment records as BED rows.
pub struct AlnBed {
    /// The layout of the records `raw` holds: BAM bytes (BAM and CRAM) or a SAM line.
    format: Format,
    name: String,
    hdr: Arc<AlnHeader>,
    o: AlnBedOptions,
    blocks: Vec<Interval>,
    row: Vec<u8>,
}

fn push_num(out: &mut Vec<u8>, n: u64) {
    out.extend_from_slice(n.to_string().as_bytes());
}

impl AlnBed {
    /// A converter for the records of `r`, or `None` when `r` does not read alignments.
    pub fn for_reader(r: &dyn RecordReader, o: AlnBedOptions) -> Option<AlnBed> {
        let (format, hdr) = raw_layout(r)?;
        Some(AlnBed {
            format,
            name: r.name().to_string(),
            hdr,
            o,
            blocks: Vec::new(),
            row: Vec::new(),
        })
    }

    fn fail(&self, loc: u64, msg: String) -> IoError {
        IoError::Parse {
            path: self.name.clone(),
            loc: if self.format == Format::Sam {
                Loc::Line(loc)
            } else {
                Loc::Record(loc)
            },
            msg,
        }
    }

    /// Append the rows of the record `raw` (record or line number `loc`, for errors) to
    /// `out`, each ending in a newline.
    pub fn write(&mut self, raw: &[u8], loc: u64, out: &mut Vec<u8>) -> Result<(), IoError> {
        let layout = self.o.layout;
        self.emit(raw, loc, layout, out)
    }

    /// The BED6 row of the record `raw`, without its newline, whatever the layout.
    pub fn bed6<'a>(&'a mut self, raw: &[u8], loc: u64) -> Result<&'a [u8], IoError> {
        let mut row = std::mem::take(&mut self.row);
        row.clear();
        let r = self.emit(raw, loc, AlnBedLayout::Bed6, &mut row);
        row.pop();
        self.row = row;
        r.map(|()| &self.row[..])
    }

    fn score(&self, a: &AlnRef<'_>, loc: u64) -> Result<i64, IoError> {
        let Some(tag) = self.o.score_tag else {
            return Ok(i64::from(a.mapq()));
        };
        let t = String::from_utf8_lossy(&tag).into_owned();
        let read = String::from_utf8_lossy(a.name()).into_owned();
        if let Some(v) = a.int_tag(tag) {
            return Ok(v);
        }
        Err(self.fail(
            loc,
            match a.tag_type(tag) {
                None => format!(
                    "read {read} has no {t} tag, which --score-tag {t} takes as the score; \
use a tag every read has (NM, AS, ...) or drop --score-tag"
                ),
                Some(ty) => format!(
                    "the {t} tag of read {read} has type {}, not an integer; --score-tag needs \
an integer (i) tag such as NM or AS",
                    ty as char
                ),
            },
        ))
    }

    fn emit(
        &mut self,
        raw: &[u8],
        loc: u64,
        layout: AlnBedLayout,
        out: &mut Vec<u8>,
    ) -> Result<(), IoError> {
        let hdr = self.hdr.clone();
        let a = aln_ref(self.format, raw, &hdr).map_err(|m| self.fail(loc, m))?;
        let (Some(rid), Some(_)) = (a.ref_id(), a.pos()) else {
            return Err(self.fail(
                loc,
                "the read has no reference or position, so it is not a mapped alignment".into(),
            ));
        };
        let (chrom, ref_len) = &hdr.refs[rid];
        let score = self.score(&a, loc)?;
        // The span every interval command uses (D101): the reference bases the CIGAR covers
        // (at least one), cut at the end of the reference.
        let Some((span, _)) = a.span(*ref_len) else {
            return Err(self.fail(loc, "the read is unmapped".into()));
        };
        let (start, end) = (span.start, span.end);
        let mut name = a.name().to_vec();
        let flag = a.flag();
        if flag & 0x40 != 0 {
            name.extend_from_slice(b"/1");
        }
        if flag & 0x80 != 0 {
            name.extend_from_slice(b"/2");
        }
        let strand = if a.strand() == Strand::Minus {
            b'-'
        } else {
            b'+'
        };
        let tail = |out: &mut Vec<u8>| {
            out.push(b'\t');
            out.extend_from_slice(&name);
            out.push(b'\t');
            out.extend_from_slice(score.to_string().as_bytes());
            out.push(b'\t');
            out.push(strand);
        };
        let head = |out: &mut Vec<u8>, s: u64, e: u64| {
            out.extend_from_slice(chrom);
            out.push(b'\t');
            push_num(out, s);
            out.push(b'\t');
            push_num(out, e);
        };
        if layout == AlnBedLayout::Bed6 {
            head(out, start, end);
            tail(out);
            if self.o.cigar {
                out.push(b'\t');
                let before = out.len();
                a.for_each_op(&mut |op, n| {
                    push_num(out, u64::from(n));
                    out.push(op);
                });
                if out.len() == before {
                    out.push(b'*');
                }
            }
            out.push(b'\n');
            return Ok(());
        }
        self.blocks.clear();
        blocks_within(&a, self.o.split_deletions, *ref_len, &mut self.blocks);
        if self.blocks.is_empty() {
            self.blocks.push(span);
        }
        if layout == AlnBedLayout::Blocks {
            for b in &self.blocks {
                head(out, b.start, b.end);
                tail(out);
                out.push(b'\n');
            }
            return Ok(());
        }
        // The row and its thick span open at the first block that is kept and close at the
        // last, so the first block starts at offset 0 and the last ends at chromEnd as BED12
        // requires (a read that starts or ends in N or, with --split-deletions, D, or is cut at
        // the reference end).
        let start = self.blocks.first().map_or(start, |b| b.start);
        let end = self.blocks.last().map_or(end, |b| b.end);
        head(out, start, end);
        tail(out);
        out.push(b'\t');
        push_num(out, start);
        out.push(b'\t');
        push_num(out, end);
        out.push(b'\t');
        out.extend_from_slice(self.o.color.as_bytes());
        out.push(b'\t');
        push_num(out, self.blocks.len() as u64);
        out.push(b'\t');
        for (i, b) in self.blocks.iter().enumerate() {
            if i > 0 {
                out.push(b',');
            }
            push_num(out, b.end - b.start);
        }
        out.push(b'\t');
        for (i, b) in self.blocks.iter().enumerate() {
            if i > 0 {
                out.push(b',');
            }
            push_num(out, b.start - start);
        }
        out.push(b'\n');
        Ok(())
    }
}
