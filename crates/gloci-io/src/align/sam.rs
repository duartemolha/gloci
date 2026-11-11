//! SAM text as alignments (D100, D104): a borrowed view of one record line, the line reader
//! shared by the interval reader and [`super::source::SamSource`], and the interval reader.
//!
//! The raw form of a SAM record is its line, without the newline, so output copies it
//! unchanged. A record's `RNAME` must be one of the header's `@SQ` lines.

use super::encode::{CIGAR_LETTERS, NO_CIGAR, ZERO_POS, encode_sam, normalise, parse_flag};
use super::header::AlnHeader;
use super::place::{IntervalRules, Verdict};
use super::record::footprint;

use super::{AlignOptions, Aln, AlnRef};
use crate::chrom_map::ChromMap;
use crate::detect::{Compression, Format};
use crate::diag::{Diagnostics, InvalidPolicy};
use crate::error::{IoError, Loc};
use crate::input::Input;
use crate::record::{Header, RecordBatch, RecordMeta, RecordReader, warn_undeclared};
use gloci_core::{ChromTable, Interval, Strand};
use memchr::memchr;
use std::borrow::Cow;
use std::io::BufRead;
use std::sync::Arc;

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
    /// The mate's reference as BAM would hold it (`=` resolved, -1 for none), and PNEXT - 1
    /// (-1 when the mate has no position).
    mate_ref: i32,
    mate_pos: i32,
    seq: &'a [u8],
    qual: &'a [u8],
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
            Some(&c) if CIGAR_LETTERS.contains(&c) && digits > 0 && digits <= 9 => c,
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
    /// reference list.
    pub fn parse(line: &'a [u8], h: &AlnHeader) -> Result<SamAln<'a>, String> {
        SamAln::parse_with(line, &|name| h.ref_id(name))
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
            rnext,
            pnext,
            _tlen,
            seq,
            qual,
        ] = f;
        if qname.is_empty() || qname.len() > 254 {
            return Err(format!(
                "the read name is {} characters long; a SAM read name has 1 to 254",
                qname.len()
            ));
        }
        let flag = parse_flag(flag)?;
        let mapq: u8 = number(mapq, "MAPQ")?;
        let pos1: u32 = number(pos, "POS")?;
        if pos1 > i32::MAX as u32 + 1 {
            return Err(format!(
                "POS {pos1} is above 2147483648; BAM positions are 32-bit"
            ));
        }
        let ref_id = if rname == b"*" {
            -1
        } else {
            match rid(rname) {
                Some(i) => i as i32,
                None => {
                    return Err(format!(
                        "RNAME '{}' is not among the @SQ lines of the header; add an \
@SQ line for it, fix the name, or pass --skip-invalid to skip such reads",
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
        if seq != b"*" && !cigar.is_empty() && query != seq.len() as u64 {
            return Err(format!(
                "the CIGAR covers {query} bases of the read but the sequence has {}",
                seq.len()
            ));
        }
        // As htslib reads it: no position means no reference, and no reference or CIGAR means
        // unmapped. The line itself is kept as it is.
        let (flag, ref_id, _) = normalise(flag, ref_id, i64::from(pos1) - 1, !cigar.is_empty());
        // The same checks of QUAL as `encode_sam` makes, with its messages, so that `quals()`
        // always agrees with `seq_len()`.
        if qual != b"*" {
            let l_seq = if seq == b"*" { 0 } else { seq.len() };
            if qual.len() != l_seq {
                return Err(format!(
                    "the quality string has {} characters but the sequence has {l_seq}",
                    qual.len()
                ));
            }
            if let Some(&q) = qual.iter().find(|q| !(33..=126).contains(*q)) {
                return Err(format!(
                    "the quality character {:?} is outside ! to ~",
                    q as char
                ));
            }
        }
        // As `encode_sam` reads the mate: a mate with no position has no reference. A PNEXT that
        // is not a number is the checker's to refuse, so it reads as no mate here.
        let mate_pos = number::<u32>(pnext, "PNEXT")
            .ok()
            .filter(|&p| p <= i32::MAX as u32 + 1)
            .map_or(-1, |p| (i64::from(p) - 1) as i32);
        let mate_ref = match rnext {
            _ if mate_pos < 0 => -1,
            b"*" => -1,
            b"=" => ref_id,
            n => rid(n).map_or(-1, |i| i as i32),
        };
        Ok(SamAln {
            line,
            ref_id,
            pos: (i64::from(pos1) - 1) as i32,
            flag,
            mapq,
            mate_ref,
            mate_pos,
            seq: if seq == b"*" { &[] } else { seq },
            qual,
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

    /// The reference of the mate, or `None` when the record has none.
    pub fn mate_ref_id(&self) -> Option<usize> {
        usize::try_from(self.mate_ref).ok()
    }

    /// The 0-based position of the mate, or `None` when the record has none.
    pub fn mate_pos(&self) -> Option<u64> {
        u64::try_from(self.mate_pos).ok()
    }

    /// The number of bases in the read (0 for SEQ `*`).
    pub fn seq_len(&self) -> usize {
        self.seq.len()
    }

    /// The raw Phred qualities: empty for SEQ `*`, and `0xFF` for every base when QUAL is `*`.
    /// [`SamAln::parse`] refuses a QUAL whose length differs from SEQ's or that holds a
    /// character outside `!` to `~` (as `encode_sam` does), so the result always has
    /// `seq_len()` values.
    pub fn quals(&self) -> Cow<'a, [u8]> {
        if self.seq.is_empty() {
            return Cow::Borrowed(&[]);
        }
        if self.qual == b"*" {
            return Cow::Owned(vec![0xFF; self.seq.len()]);
        }
        Cow::Owned(self.qual.iter().map(|q| q.wrapping_sub(33)).collect())
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

    /// The type character of an optional field (`i`, `f`, `A`, `Z`, `H` or `B`).
    pub fn tag_type(&self, tag: [u8; 2]) -> Option<u8> {
        self.tags
            .split(|&b| b == b'\t')
            .find(|f| f.len() > 4 && f[..2] == tag && f[2] == b':' && f[4] == b':')
            .map(|f| f[3])
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
    /// Where records are encoded to check them.
    scratch: Vec<u8>,
    /// The notes of `normalise` already warned about.
    warned: u8,
    /// BAM can hold the record in `scratch`.
    fits: bool,
}

/// What [`SamLines::check`] found out about a record line: the flag and MAPQ as htslib reads
/// them (what the filters judge), and whether the line says something else.
pub(crate) struct Checked {
    pub(crate) flag: u16,
    pub(crate) mapq: u8,
    normalised: bool,
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
            // CRLF line ends become LF, as htslib reads them.
            if text.ends_with(b"\r\n") {
                let n = text.len();
                text.remove(n - 2);
            }
            line += 1;
            if text[before..].contains(&0) {
                return Err(IoError::Detect {
                    path: name.clone(),
                    msg: format!("{name} line {line}: the header line contains a NUL byte"),
                });
            }
        }
        let mut h = AlnHeader::parse_sam(text).map_err(|msg| IoError::Detect {
            path: name.clone(),
            msg: format!("{name}: {msg}"),
        })?;
        if let Some(m) = map {
            m.check_file_names(&name, h.refs.iter().map(|r| &r.0[..]))?;
            h = h.renamed(m);
        }
        Ok(SamLines {
            aln: Arc::new(h),
            name,
            compression,
            src: reader,
            line,
            rec: 0,
            scratch: Vec::new(),
            warned: 0,
            fits: false,
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

    /// Check the record line `line` (from [`SamLines::next`]) as far as writing it as BAM
    /// would: every field and optional field, so a bad line is an error where it is read,
    /// whatever the output. The record is encoded into `self.scratch` on the way, and
    /// [`SamLines::bam`] gives the bytes. A 32-bit overflow is no error here, because SAM can
    /// hold what BAM cannot.
    pub(crate) fn check(&mut self, line: &[u8]) -> Result<Checked, String> {
        self.scratch.clear();
        let enc = encode_sam(line, &|n| self.aln.ref_id(n), &mut self.scratch, false)?;
        self.fits = enc.not_bam.is_none();
        // htslib treats such reads as unmapped and says so; once per input here.
        for (note, what) in [
            (ZERO_POS, "has a reference but no position"),
            (NO_CIGAR, "is mapped but has no CIGAR"),
        ] {
            if enc.notes & note != 0 && self.warned & note == 0 {
                self.warned |= note;
                let name = line.split(|&b| b == b'\t').next().unwrap_or(&[]);
                log::warn!(
                    "{} line {}: read {} {what}; reads like it are treated as unmapped",
                    self.name,
                    self.line,
                    String::from_utf8_lossy(name)
                );
            }
        }
        Ok(Checked {
            flag: u16::from_le_bytes([self.scratch[14], self.scratch[15]]),
            mapq: self.scratch[9],
            normalised: enc.normalised && self.fits,
        })
    }

    /// The BAM record (without `block_size`) of the line last checked, when BAM can hold it.
    pub(crate) fn bam(&self) -> Option<&[u8]> {
        self.fits.then_some(&self.scratch[..])
    }

    /// Where htslib reads the record differently from how the line is written (its flag,
    /// reference or mate changed), replace the line in `buf[start..]` by the SAM text of the
    /// record it reads, as samtools writes it. Other lines stay as written (D120).
    pub(crate) fn rewrite_normalised(&self, buf: &mut Vec<u8>, start: usize, c: &Checked) {
        if !c.normalised {
            return;
        }
        let mut text = Vec::new();
        if Aln::parse(&self.scratch, &self.aln)
            .is_ok_and(|a| a.write_sam(&self.aln, &mut text).is_ok())
        {
            buf.truncate(start);
            buf.extend_from_slice(&text);
        }
    }

    /// As [`SamLines::check`], and the view of the line for the interval readers.
    pub(crate) fn parse<'a>(&mut self, line: &'a [u8]) -> Result<(SamAln<'a>, Checked), String> {
        let a = SamAln::parse(line, &self.aln)?;
        let c = self.check(line)?;
        Ok((a, c))
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
            // CRLF line ends become LF.
            if buf.len() > start && buf.last() == Some(&b'\r') {
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
    /// The slim line being built (see [`slim_sam_line`]).
    scratch: Vec<u8>,
}

/// Replace the SAM line `buf[start..]` with its slim form for the interval commands that only
/// compare reads (D115), as BAM records are slimmed: the first nine fields, `*` for SEQ and
/// QUAL, and of the optional fields only `keep` (`--score-tag`). What reads a slim line (the
/// BED6 row, the blocks) uses QNAME, FLAG, RNAME, POS, MAPQ, CIGAR and that tag.
fn slim_sam_line(buf: &mut Vec<u8>, start: usize, keep: Option<[u8; 2]>, out: &mut Vec<u8>) {
    out.clear();
    for (i, field) in buf[start..].split(|&b| b == b'\t').enumerate() {
        match i {
            0..=8 => {
                if i > 0 {
                    out.push(b'\t');
                }
                out.extend_from_slice(field);
            }
            9 | 10 => out.extend_from_slice(b"\t*"),
            _ if keep.is_some_and(|t| field.starts_with(&t) && field.get(2) == Some(&b':')) => {
                out.push(b'\t');
                out.extend_from_slice(field);
            }
            _ => {}
        }
    }
    buf.truncate(start);
    buf.extend_from_slice(out);
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
            scratch: Vec::new(),
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
    fn whole_records(&self) -> bool {
        self.opts.full_records
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
            let (a, checked) = match self.lines.parse(&batch.buf[start..]) {
                Ok(r) => r,
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
                warn_undeclared(chroms, &self.lines.name, Some(Loc::Line(self.lines.line)));
            }
            self.lines
                .rewrite_normalised(&mut batch.buf, start, &checked);
            if !self.opts.full_records {
                slim_sam_line(
                    &mut batch.buf,
                    start,
                    self.opts.score_tag,
                    &mut self.scratch,
                );
            }
            batch.recs.push(RecordMeta {
                chrom: p.chrom,
                key: p.pos,
                span: p.span,
                strand: p.strand,
                loc: self.lines.line,
                off: start,
                len: (batch.buf.len() - start) as u32,
            });
            n += 1;
        }
        Ok(true)
    }
}
