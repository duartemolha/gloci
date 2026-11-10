//! SAM text as alignments (D100, D104): a borrowed view of one record line, the line reader
//! shared by the interval reader and [`super::source::SamSource`], and the interval reader.
//!
//! The raw form of a SAM record is its line, without the newline, so output copies it
//! unchanged. A record's `RNAME` must be one of the header's `@SQ` lines.

use super::header::AlnHeader;
use super::place::{IntervalRules, Verdict};
use super::record::footprint;
use super::{AlignOptions, AlnRef};
use crate::chrom_map::ChromMap;
use crate::detect::{Compression, Format};
use crate::diag::{Diagnostics, InvalidPolicy};
use crate::error::{IoError, Loc};
use crate::input::Input;
use crate::record::{Header, RecordBatch, RecordMeta, RecordReader, warn_undeclared};
use gloci_core::{ChromTable, Interval, Strand};
use memchr::memchr;
use std::collections::HashMap;
use std::io::BufRead;
use std::sync::Arc;

const FLAG_UNMAPPED: u16 = 0x4;
const FLAG_REVERSE: u16 = 0x10;

/// A view of one SAM record line, checked by [`SamAln::parse`]. It has the methods of
/// [`super::Aln`].
#[derive(Debug, Clone, Copy)]
pub struct SamAln<'a> {
    line: &'a [u8],
    ref_id: i32,
    pos: i32,
    flag: u16,
    mapq: u8,
    name: &'a [u8],
    /// The CIGAR text, empty for `*`.
    cigar: &'a [u8],
    /// The optional fields, `TAG:TYPE:VALUE` separated by tabs; empty when there are none.
    tags: &'a [u8],
}

fn number<T: std::str::FromStr>(field: &[u8], what: &str) -> Result<T, String> {
    std::str::from_utf8(field)
        .ok()
        .and_then(|s| s.parse::<T>().ok())
        .ok_or_else(|| {
            format!(
                "{what} '{}' is not a valid number in the range SAM allows",
                String::from_utf8_lossy(field)
            )
        })
}

/// The `(letter, length)` operations of CIGAR text, or `None` if it is malformed.
fn cigar_ops(text: &[u8]) -> impl Iterator<Item = Option<(u8, u32)>> + '_ {
    let mut rest = text;
    std::iter::from_fn(move || {
        if rest.is_empty() {
            return None;
        }
        let digits = rest.iter().take_while(|b| b.is_ascii_digit()).count();
        let op = match rest.get(digits) {
            Some(&c) if b"MIDNSHP=X".contains(&c) && digits > 0 && digits <= 9 => c,
            _ => {
                rest = &[];
                return Some(None);
            }
        };
        let n = std::str::from_utf8(&rest[..digits])
            .ok()
            .and_then(|s| s.parse::<u32>().ok());
        rest = &rest[digits + 1..];
        Some(n.map(|n| (op, n)))
    })
}

impl<'a> SamAln<'a> {
    /// Check and view the line `line` (no newline), finding its `RNAME` in `h` by a scan of the
    /// reference list. Readers that see every line use [`SamAln::parse_with`] and a map.
    pub fn parse(line: &'a [u8], h: &AlnHeader) -> Result<SamAln<'a>, String> {
        SamAln::parse_with(line, &|name| h.refs.iter().position(|(n, _)| n == name))
    }

    /// As [`SamAln::parse`], with `rid` finding a reference by name.
    pub fn parse_with(
        line: &'a [u8],
        rid: &dyn Fn(&[u8]) -> Option<usize>,
    ) -> Result<SamAln<'a>, String> {
        let mut f: [&[u8]; 11] = [&[]; 11];
        let mut rest = line;
        for (i, slot) in f.iter_mut().enumerate() {
            if i == 10 {
                // The last mandatory field ends at a TAB, or at the end of the line.
                let end = memchr(b'\t', rest).unwrap_or(rest.len());
                *slot = &rest[..end];
                rest = rest.get(end + 1..).unwrap_or(&[]);
                break;
            }
            let Some(end) = memchr(b'\t', rest) else {
                return Err(format!(
                    "the line has {} tab-separated fields; a SAM record has at least 11 \
(QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL)",
                    i + 1
                ));
            };
            *slot = &rest[..end];
            rest = &rest[end + 1..];
        }
        let [
            qname,
            flag,
            rname,
            pos,
            mapq,
            cigar,
            _rnext,
            _pnext,
            _tlen,
            seq,
            _qual,
        ] = f;
        if qname.is_empty() || qname.len() > 254 {
            return Err(format!(
                "the read name is {} characters long; a SAM read name has 1 to 254",
                qname.len()
            ));
        }
        let flag: u16 = number(flag, "FLAG")?;
        let mapq: u8 = number(mapq, "MAPQ")?;
        let pos1: u32 = number(pos, "POS")?;
        if pos1 > i32::MAX as u32 {
            return Err(format!("POS {pos1} is above 2147483647"));
        }
        let ref_id = if rname == b"*" {
            -1
        } else {
            match rid(rname) {
                Some(i) => i as i32,
                None => {
                    return Err(format!(
                        "RNAME '{}' is not among the @SQ lines of the header; add an \
@SQ line for it, or fix the name",
                        String::from_utf8_lossy(rname)
                    ));
                }
            }
        };
        let cigar = if cigar == b"*" { &[][..] } else { cigar };
        let mut query = 0u64;
        for op in cigar_ops(cigar) {
            let Some((c, n)) = op else {
                return Err(format!(
                    "CIGAR '{}' is not a list of operations such as 10M2D5M",
                    String::from_utf8_lossy(cigar)
                ));
            };
            if b"MI=XS".contains(&c) {
                query += u64::from(n);
            }
        }
        if flag & FLAG_UNMAPPED == 0 {
            if ref_id < 0 || pos1 == 0 {
                return Err(format!(
                    "a mapped read (flag {flag}) has RNAME '{}' and POS {pos1}; the flag \
should include 4 (unmapped), or the record needs a reference and position",
                    String::from_utf8_lossy(rname)
                ));
            }
            if seq != b"*" && !cigar.is_empty() && query != seq.len() as u64 {
                return Err(format!(
                    "the CIGAR covers {query} bases of the read but the sequence has {}",
                    seq.len()
                ));
            }
        }
        Ok(SamAln {
            line,
            ref_id,
            pos: pos1 as i32 - 1,
            flag,
            mapq,
            name: qname,
            cigar,
            tags: rest,
        })
    }

    pub fn ref_id(&self) -> Option<usize> {
        usize::try_from(self.ref_id).ok()
    }

    /// The 0-based position, or `None` when the record has none.
    pub fn pos(&self) -> Option<u64> {
        u64::try_from(self.pos).ok()
    }

    pub fn flag(&self) -> u16 {
        self.flag
    }

    pub fn mapq(&self) -> u8 {
        self.mapq
    }

    pub fn name(&self) -> &'a [u8] {
        self.name
    }

    /// The record line.
    pub fn raw(&self) -> &'a [u8] {
        self.line
    }

    /// Call `f(letter, length)` for each CIGAR operation.
    pub fn for_each_op(&self, f: &mut dyn FnMut(u8, u32)) {
        for (op, n) in cigar_ops(self.cigar).map_while(|o| o) {
            f(op, n);
        }
    }

    /// An integer (`i`) optional field.
    pub fn int_tag(&self, tag: [u8; 2]) -> Option<i64> {
        for field in self.tags.split(|&b| b == b'\t') {
            if field.len() > 5 && field[..2] == tag && &field[2..5] == b":i:" {
                return std::str::from_utf8(&field[5..]).ok()?.parse().ok();
            }
        }
        None
    }

    /// The footprint of a mapped read, cut at `ref_len`, and whether it was cut (D101).
    pub fn span(&self, ref_len: u64) -> Option<(Interval, bool)> {
        let covered: u64 = cigar_ops(self.cigar)
            .map_while(|o| o)
            .filter(|(op, _)| b"MDN=X".contains(op))
            .map(|(_, n)| u64::from(n))
            .sum();
        footprint(self.flag, self.pos(), covered.max(1), ref_len)
    }

    pub fn strand(&self) -> Strand {
        if self.flag & FLAG_REVERSE != 0 {
            Strand::Minus
        } else {
            Strand::Plus
        }
    }
}

/// The header and lines of a SAM file. Header lines come first; blank lines are dropped and
/// counted; the record numbers are the same as the BAM reader's.
pub(crate) struct SamLines {
    pub(crate) name: String,
    pub(crate) compression: Compression,
    pub(crate) aln: Arc<AlnHeader>,
    src: Box<dyn BufRead + Send>,
    /// The line last read, counted from 1 over the whole file.
    pub(crate) line: u64,
    /// Records returned so far.
    pub(crate) rec: u64,
    by_name: HashMap<Vec<u8>, usize>,
}

impl SamLines {
    pub(crate) fn open(input: Input, map: Option<&ChromMap>) -> Result<SamLines, IoError> {
        let Input {
            name,
            compression,
            mut reader,
            ..
        } = input;
        let read_err = |source| IoError::Read {
            path: name.clone(),
            source,
        };
        if map.is_some() {
            return Err(IoError::Detect {
                path: name.clone(),
                msg: format!(
                    "{name} is SAM, and --chrom-map cannot rename the reference names in SAM \
text; convert it to BAM first (gloci view {name} -O bam -o x.bam) and read that"
                ),
            });
        }
        let mut text = Vec::new();
        let mut line = 0u64;
        loop {
            let first = reader.fill_buf().map_err(&read_err)?.first().copied();
            if first != Some(b'@') {
                break;
            }
            let before = text.len();
            reader.read_until(b'\n', &mut text).map_err(&read_err)?;
            if text.last() != Some(&b'\n') {
                text.push(b'\n');
            }
            line += 1;
            if text[before..].contains(&0) {
                return Err(IoError::Detect {
                    path: name.clone(),
                    msg: format!("{name} line {line}: the header line contains a NUL byte"),
                });
            }
        }
        let h = AlnHeader::parse_sam(text).map_err(|msg| IoError::Detect {
            path: name.clone(),
            msg: format!("{name}: {msg}"),
        })?;
        let by_name = h
            .refs
            .iter()
            .enumerate()
            .map(|(i, (n, _))| (n.clone(), i))
            .collect();
        Ok(SamLines {
            aln: Arc::new(h),
            name,
            compression,
            src: reader,
            line,
            rec: 0,
            by_name,
        })
    }

    pub(crate) fn data_err(&self, msg: &str) -> IoError {
        IoError::Detect {
            path: self.name.clone(),
            msg: format!("{} line {}: {msg}", self.name, self.line),
        }
    }

    /// An invalid record: an error, or counted in `diag` under `--skip-invalid`.
    pub(crate) fn invalid(
        &self,
        policy: InvalidPolicy,
        diag: &mut Diagnostics,
        msg: &str,
    ) -> Result<(), IoError> {
        match policy {
            InvalidPolicy::Fail => Err(self.data_err(msg)),
            InvalidPolicy::Skip => {
                diag.record(self.line, msg);
                Ok(())
            }
        }
    }

    /// The record line `line` (from [`SamLines::next`]), checked.
    pub(crate) fn parse<'a>(&self, line: &'a [u8]) -> Result<SamAln<'a>, String> {
        SamAln::parse_with(line, &|n| self.by_name.get(n).copied())
    }

    /// Append the next record line, without its newline, to `buf`. `Ok(false)` at the end of
    /// the file. Blank lines are noted in `diag` and skipped.
    pub(crate) fn next(
        &mut self,
        buf: &mut Vec<u8>,
        diag: &mut Diagnostics,
    ) -> Result<bool, IoError> {
        loop {
            let start = buf.len();
            let got = self
                .src
                .read_until(b'\n', buf)
                .map_err(|source| IoError::Read {
                    path: self.name.clone(),
                    source,
                })?;
            if got == 0 {
                return Ok(false);
            }
            self.line += 1;
            if buf.last() == Some(&b'\n') {
                buf.pop();
            }
            match buf.get(start) {
                None => {
                    diag.record_blank(self.line);
                    continue;
                }
                Some(b'@') => {
                    let e = self.data_err(
                        "a header line after the first record; the @HD, @SQ, ... lines must \
come before the records",
                    );
                    buf.truncate(start);
                    return Err(e);
                }
                Some(_) => {}
            }
            self.rec += 1;
            return Ok(true);
        }
    }
}

/// Reads SAM records for interval commands.
pub struct SamReader {
    lines: SamLines,
    header: Header,
    diag: Diagnostics,
    policy: InvalidPolicy,
    opts: AlignOptions,
    rules: IntervalRules,
}

impl SamReader {
    pub(crate) fn open_sam(
        input: Input,
        policy: InvalidPolicy,
        map: Option<&ChromMap>,
        opts: AlignOptions,
    ) -> Result<SamReader, IoError> {
        let lines = SamLines::open(input, map)?;
        Ok(SamReader {
            header: Header::Alignment(lines.aln.clone()),
            rules: IntervalRules::new(lines.aln.refs.len()),
            lines,
            diag: Diagnostics::default(),
            policy,
            opts,
        })
    }

    /// An invalid record: an error, or counted under `--skip-invalid`.
    fn bad(&mut self, msg: &str) -> Result<(), IoError> {
        self.lines.invalid(self.policy, &mut self.diag, msg)
    }
}

impl RecordReader for SamReader {
    fn name(&self) -> &str {
        &self.lines.name
    }
    fn format(&self) -> Format {
        Format::Sam
    }
    fn compression(&self) -> Compression {
        self.lines.compression
    }
    fn header(&self) -> &Header {
        &self.header
    }
    fn diagnostics(&self) -> &Diagnostics {
        &self.diag
    }

    fn read_batch(
        &mut self,
        chroms: &mut ChromTable,
        batch: &mut RecordBatch,
        max: usize,
    ) -> Result<bool, IoError> {
        debug_assert!(max > 0, "read_batch max must be > 0");
        let hdr = self.lines.aln.clone();
        let mut n = 0;
        while n < max {
            let start = batch.buf.len();
            if !self.lines.next(&mut batch.buf, &mut self.diag)? {
                return Ok(false);
            }
            let a = match self.lines.parse(&batch.buf[start..]) {
                Ok(a) => a,
                Err(msg) => {
                    batch.buf.truncate(start);
                    self.bad(&msg)?;
                    continue;
                }
            };
            let verdict = self.rules.judge(
                &self.lines.name,
                &AlnRef::Sam(a),
                &hdr,
                &self.opts,
                chroms,
                &mut self.diag,
            );
            let p = match verdict {
                Verdict::Keep(p) => p,
                Verdict::Skip => {
                    batch.buf.truncate(start);
                    continue;
                }
                Verdict::Bad(msg) => {
                    batch.buf.truncate(start);
                    self.bad(&msg)?;
                    continue;
                }
            };
            if chroms.has_undeclared() {
                warn_undeclared(chroms, &self.lines.name, Some(Loc::Record(self.lines.rec)));
            }
            batch.recs.push(RecordMeta {
                chrom: p.chrom,
                key: p.pos,
                span: p.span,
                strand: p.strand,
                loc: self.lines.rec,
                off: start,
                len: (batch.buf.len() - start) as u32,
            });
            n += 1;
        }
        Ok(true)
    }
}
