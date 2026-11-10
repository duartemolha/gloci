//! BEDPE reader (D146): `chrom1 start1 end1 chrom2 start2 end2 [name score
//! strand1 strand2 extra...]`, one pair of ends per line.
//!
//! BEDPE is never detected (a BED3 file of pairs has the same shape); it is
//! opened by the commands that read it (`pair-intersect`, `pair-pair`,
//! `shuffle --format bedpe`), through [`BedpeReader`] and not through the
//! one-interval-per-record `RecordReader`.

use crate::bed::is_header_line;
use crate::chrom_map::ChromMap;
use crate::detect::{Compression, Format, detect_format};
use crate::diag::Diagnostics;
use crate::error::{GateKind, IoError};
use crate::fields::{column_count, field, split_fields};
use crate::input::{Input, InputOptions, open_input, open_input_from};
use crate::record::{
    Class, Parsed, ReaderOptions, RecordBatch, RecordReader, TextKind, TextReader, coord,
    intern_chrom,
};
use gloci_core::{ChromId, ChromTable, Interval, Strand};
use std::cell::RefCell;
use std::sync::Arc;

/// One end of a pair: where it is, and on which strand.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BedpeEnd {
    pub chrom: ChromId,
    pub span: Interval,
    pub strand: Strand,
}

/// One record's ends; its raw bytes live in `BedpeBatch::buf`. An end whose
/// chromosome is `.` is `None`: it overlaps nothing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BedpeMeta {
    pub ends: [Option<BedpeEnd>; 2],
    /// 1-based physical line number.
    pub loc: u64,
    /// Offset of the raw bytes in `BedpeBatch::buf`.
    pub off: usize,
    pub len: u32,
}

/// Records of one batch: raw bytes (no line terminators) plus metadata.
#[derive(Debug, Clone, Default)]
pub struct BedpeBatch {
    pub buf: Vec<u8>,
    pub recs: Vec<BedpeMeta>,
}

impl BedpeBatch {
    pub fn raw(&self, i: usize) -> &[u8] {
        let m = &self.recs[i];
        &self.buf[m.off..m.off + m.len as usize]
    }

    pub fn len(&self) -> usize {
        self.recs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.recs.is_empty()
    }

    pub fn clear(&mut self) {
        self.buf.clear();
        self.recs.clear();
    }
}

/// The name of a record: column 7, when the line has one.
pub fn bedpe_name(raw: &[u8]) -> Option<&[u8]> {
    field(raw, 6)
}

type Ends = [Option<BedpeEnd>; 2];

/// The line parser. It hands the text machinery the first present end, and
/// queues both ends for [`BedpeReader`], so a line is parsed once.
pub(crate) struct Bedpe {
    map: Option<Arc<ChromMap>>,
    queue: RefCell<Vec<Ends>>,
}

impl TextKind for Bedpe {
    fn classify(line: &[u8], _seen: bool) -> Class {
        if line.is_empty() {
            Class::Skip
        } else if is_header_line(line) {
            Class::Header
        } else {
            Class::Record
        }
    }

    fn parse(&self, line: &[u8], chroms: &mut ChromTable) -> Result<Parsed, String> {
        let ends = parse_ends(self.map.as_deref(), line, chroms)?;
        let first = ends.iter().flatten().next().copied();
        self.queue.borrow_mut().push(ends);
        // A record with no present end has no footprint; the chromosome and key
        // below are never read (the reader keeps only `ends`).
        Ok(match first {
            Some(e) => (e.chrom, e.span.start, e.span, e.strand),
            None => (ChromId(0), 0, Interval::point(0), Strand::Unknown),
        })
    }
}

fn strand(col: usize, name: &str, v: &[u8]) -> Result<Strand, String> {
    match v {
        b"+" => Ok(Strand::Plus),
        b"-" => Ok(Strand::Minus),
        b"." => Ok(Strand::Unknown),
        _ => Err(format!(
            "column {col} ({name}) is \"{}\"; a strand is +, - or .",
            String::from_utf8_lossy(v)
        )),
    }
}

/// A coordinate of an absent end: any integer (`-1` is the convention).
fn any_int(label: &str, v: &[u8]) -> Result<(), String> {
    let digits = v.strip_prefix(b"-").unwrap_or(v);
    if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
        return Err(format!(
            "{label} \"{}\" is not an integer",
            String::from_utf8_lossy(v)
        ));
    }
    Ok(())
}

fn parse_end(
    map: Option<&ChromMap>,
    chroms: &mut ChromTable,
    k: u8,
    f: [&[u8]; 3],
    strand: Strand,
) -> Result<Option<BedpeEnd>, String> {
    let [c, s, e] = f;
    let (ls, le) = (format!("start{k}"), format!("end{k}"));
    if c == b"." {
        any_int(&ls, s)?;
        any_int(&le, e)?;
        return Ok(None);
    }
    let start = coord(&ls, s)?;
    let end = coord(&le, e)?;
    if start > end {
        return Err(format!(
            "{ls} {start} is greater than {le} {end}; each end must have start <= end"
        ));
    }
    let name = map.map_or(c, |m| m.rename(c));
    let chrom = intern_chrom(chroms, name)?;
    Ok(Some(BedpeEnd {
        chrom,
        span: Interval { start, end },
        strand,
    }))
}

fn parse_ends(
    map: Option<&ChromMap>,
    line: &[u8],
    chroms: &mut ChromTable,
) -> Result<Ends, String> {
    let n = column_count(line);
    if !matches!(n, 6..=8) && n < 10 {
        let hint = if memchr::memchr(b' ', line).is_some() && n == 1 {
            " (the line contains spaces; BEDPE must be tab-separated)"
        } else {
            ""
        };
        return Err(format!(
            "expected 6, 7, 8 or at least 10 columns, found {n}{hint}; a BEDPE line is \
chrom1 start1 end1 chrom2 start2 end2 [name score strand1 strand2 ...]"
        ));
    }
    let mut f = [&b""[..]; 10];
    for (slot, v) in f.iter_mut().zip(split_fields(line)) {
        *slot = v;
    }
    let (s1, s2) = if n >= 10 {
        (strand(9, "strand1", f[8])?, strand(10, "strand2", f[9])?)
    } else {
        (Strand::Unknown, Strand::Unknown)
    };
    let a = parse_end(map, chroms, 1, [f[0], f[1], f[2]], s1)?;
    let b = parse_end(map, chroms, 2, [f[3], f[4], f[5]], s2)?;
    Ok([a, b])
}

/// Reads a BEDPE file (plain, gzip or BGZF; `-` is stdin) in batches.
pub struct BedpeReader {
    inner: TextReader<Bedpe>,
    scratch: RecordBatch,
    header: Vec<Vec<u8>>,
    map: Option<Arc<ChromMap>>,
}

/// Why a file that is plainly another format is refused: `None` when it is
/// not recognisably one. Only formats with a signature of their own count;
/// a headerless GFF falls through to the column check.
fn refused_format(head: &[u8]) -> Option<Format> {
    let f = detect_format(head).ok()?;
    let first = head
        .iter()
        .position(|b| !b.is_ascii_whitespace())
        .map_or(&head[..0], |i| &head[i..]);
    let sure = match f {
        Format::Vcf | Format::Bcf | Format::Bam | Format::Cram | Format::Fasta => true,
        Format::Gff3 | Format::Gtf => first.starts_with(b"##gff-version"),
        Format::Sam => first.starts_with(b"@"),
        _ => false,
    };
    sure.then_some(f)
}

fn bedpe_input(o: &ReaderOptions) -> InputOptions {
    InputOptions {
        threads: o.threads,
        format: Some(Format::Bedpe),
    }
}

impl BedpeReader {
    /// Open `path` as BEDPE whatever `o.format` says: BEDPE is never detected.
    /// A VCF, BCF, GFF or alignment file is refused (exit 2).
    pub fn open(path: &str, o: &ReaderOptions) -> Result<Self, IoError> {
        Self::from_input(open_input(path, &bedpe_input(o))?, o)
    }

    /// Like [`BedpeReader::open`], over a stream called `name`.
    pub fn open_from(
        name: &str,
        src: Box<dyn std::io::Read + Send>,
        o: &ReaderOptions,
    ) -> Result<Self, IoError> {
        Self::from_input(open_input_from(name, src, &bedpe_input(o))?, o)
    }

    fn from_input(mut input: Input, o: &ReaderOptions) -> Result<Self, IoError> {
        let head =
            std::io::BufRead::fill_buf(&mut input.reader).map_err(|source| IoError::Read {
                path: input.name.clone(),
                source,
            })?;
        if let Some(f) = refused_format(head) {
            let name = &input.name;
            let tail = if f.is_alignment() {
                format!("; to make pairs from alignments run: gloci convert {name} --bedpe")
            } else {
                String::new()
            };
            return Err(IoError::Gate {
                kind: GateKind::PairedOnly,
                msg: format!(
                    "{name} is {}, but this command reads BEDPE (chrom1 start1 end1 chrom2 start2 \
end2 ...); give it a BEDPE file{tail}",
                    f.name()
                ),
            });
        }
        let kind = Bedpe {
            map: o.chrom_map.clone(),
            queue: RefCell::new(Vec::new()),
        };
        let inner = TextReader::open(input, o.invalid, kind)?;
        let header = match &inner.header {
            crate::record::Header::Lines(l) => l
                .split(|&b| b == b'\n')
                .filter(|x| !x.is_empty())
                .map(<[u8]>::to_vec)
                .collect(),
            _ => Vec::new(),
        };
        Ok(BedpeReader {
            inner,
            scratch: RecordBatch::default(),
            header,
            map: o.chrom_map.clone(),
        })
    }

    pub fn name(&self) -> &str {
        self.inner.name()
    }

    pub fn compression(&self) -> Compression {
        self.inner.compression()
    }

    /// The leading `#`, `track` and `browser` lines, without terminators.
    pub fn header(&self) -> &[Vec<u8>] {
        &self.header
    }

    pub fn diagnostics(&self) -> &Diagnostics {
        self.inner.diagnostics()
    }

    /// Append up to `max` records to `batch` (the caller clears it); false
    /// once the input is exhausted (records appended by that call are valid).
    pub fn read_batch(
        &mut self,
        chroms: &mut ChromTable,
        batch: &mut BedpeBatch,
        max: usize,
    ) -> Result<bool, IoError> {
        self.scratch.clear();
        let more = self.inner.read_batch(chroms, &mut self.scratch, max);
        let ends = std::mem::take(&mut *self.inner.kind().queue.borrow_mut());
        let more = more?;
        if ends.len() != self.scratch.recs.len() {
            return Err(IoError::Internal(format!(
                "BEDPE reader parsed {} lines into {} records",
                ends.len(),
                self.scratch.recs.len()
            )));
        }
        let base = batch.buf.len();
        match &self.map {
            // Only the chromosome columns change; the rest of the line is kept.
            Some(m) => {
                for (r, e) in self.scratch.recs.iter().zip(ends) {
                    let raw = &self.scratch.buf[r.off..r.off + r.len as usize];
                    let off = batch.buf.len();
                    rewrite_chroms(m, raw, &mut batch.buf);
                    batch.recs.push(BedpeMeta {
                        ends: e,
                        loc: r.loc,
                        off,
                        len: (batch.buf.len() - off) as u32,
                    });
                }
            }
            None => {
                if base == 0 {
                    std::mem::swap(&mut batch.buf, &mut self.scratch.buf);
                } else {
                    batch.buf.extend_from_slice(&self.scratch.buf);
                }
                for (r, e) in self.scratch.recs.iter().zip(ends) {
                    batch.recs.push(BedpeMeta {
                        ends: e,
                        loc: r.loc,
                        off: base + r.off,
                        len: r.len,
                    });
                }
            }
        }
        Ok(more)
    }
}

/// Copy `raw` to `out`, renaming columns 1 and 4.
fn rewrite_chroms(map: &ChromMap, raw: &[u8], out: &mut Vec<u8>) {
    for (i, f) in split_fields(raw).enumerate() {
        if i > 0 {
            out.push(b'\t');
        }
        out.extend_from_slice(if i == 0 || i == 3 { map.rename(f) } else { f });
    }
}
