//! Record batches and the reader interface shared by every command.

use crate::detect::{Compression, Format};
use crate::diag::{Diagnostics, InvalidPolicy};
use crate::error::{IoError, Loc};
use crate::input::{Input, InputOptions, open_input, open_input_from};
use gloci_core::{ChromId, ChromTable, CoordError, Interval, OrderSource, VcfSpan};
use std::io::{BufRead, Read};
use std::num::NonZeroUsize;

/// One record's position data; its raw bytes live in `RecordBatch::buf`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RecordMeta {
    pub chrom: ChromId,
    /// The position the file is sorted by (BED start, GFF start-1, VCF POS-1).
    pub key: u64,
    /// The record's footprint.
    pub span: Interval,
    /// 1-based physical line number (or record number for binary formats).
    pub loc: u64,
    /// Offset of the raw bytes in `RecordBatch::buf`.
    pub off: usize,
    pub len: u32,
}

/// Records of one batch: raw bytes (no line terminators) plus metadata.
#[derive(Debug, Clone, Default)]
pub struct RecordBatch {
    pub buf: Vec<u8>,
    pub recs: Vec<RecordMeta>,
}

impl RecordBatch {
    pub fn raw(&self, i: usize) -> &[u8] {
        let m = &self.recs[i];
        &self.buf[m.off..m.off + m.len as usize]
    }

    pub fn clear(&mut self) {
        self.buf.clear();
        self.recs.clear();
    }

    pub fn len(&self) -> usize {
        self.recs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.recs.is_empty()
    }
}

/// Header lines read before the first record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Header {
    None,
    /// Header lines, each terminated by `\n`.
    Lines(Vec<u8>),
    Vcf(crate::vcf::VcfHeader),
}

/// Bytes of records to read per batch, and the bound on the record count.
const BATCH_BYTES: usize = 4 << 20;
const MAX_BATCH: usize = 16384;

/// Records to request next, so a batch is about `BATCH_BYTES` whatever the
/// record length (a VCF line can be megabytes).
pub fn next_batch_size(batch: &RecordBatch, current: usize) -> usize {
    if batch.is_empty() {
        return current;
    }
    let per_rec = batch.buf.len() / batch.len() + 1;
    (BATCH_BYTES / per_rec).clamp(1, MAX_BATCH)
}

pub trait RecordReader: Send {
    fn name(&self) -> &str;
    fn format(&self) -> Format;
    fn compression(&self) -> Compression;
    fn header(&self) -> &Header;
    /// Append up to `max` records to `batch` (the caller clears it).
    /// Returns `Ok(true)` while more input may remain and `Ok(false)` once the
    /// input is exhausted; records appended by that final call are valid.
    fn read_batch(
        &mut self,
        chroms: &mut ChromTable,
        batch: &mut RecordBatch,
        max: usize,
    ) -> Result<bool, IoError>;
    fn diagnostics(&self) -> &Diagnostics;
}

#[derive(Debug, Clone)]
pub struct ReaderOptions {
    pub threads: NonZeroUsize,
    pub format: Option<Format>,
    pub vcf_span: VcfSpan,
    pub invalid: InvalidPolicy,
}

pub fn open_reader(path: &str, o: &ReaderOptions) -> Result<Box<dyn RecordReader>, IoError> {
    dispatch(open_input(path, &input_options(o))?, o)
}

pub fn open_reader_from(
    name: &str,
    src: Box<dyn Read + Send>,
    o: &ReaderOptions,
) -> Result<Box<dyn RecordReader>, IoError> {
    dispatch(open_input_from(name, src, &input_options(o))?, o)
}

fn input_options(o: &ReaderOptions) -> InputOptions {
    InputOptions {
        threads: o.threads,
        format: o.format,
    }
}

fn dispatch(input: Input, o: &ReaderOptions) -> Result<Box<dyn RecordReader>, IoError> {
    match input.format {
        Format::Bed | Format::BedGraph => {
            Ok(Box::new(crate::bed::BedReader::open_bed(input, o.invalid)?))
        }
        Format::Gff3 | Format::Gtf => {
            Ok(Box::new(crate::gff::GffReader::open_gff(input, o.invalid)?))
        }
        Format::Vcf => Ok(Box::new(crate::vcf::VcfReader::open_vcf(
            input, o.invalid, o.vcf_span,
        )?)),
        Format::Bcf => Ok(Box::new(crate::bcf::BcfReader::open_bcf(
            input, o.invalid, o.vcf_span,
        )?)),
    }
}

// ---- shared machinery for line-oriented text formats ----

pub(crate) enum Class {
    /// Header line if before the first record, otherwise dropped.
    Header,
    Skip,
    Record,
    Fatal(String),
}

pub(crate) trait TextKind: Send + 'static {
    fn classify(line: &[u8], seen_record: bool) -> Class;
    fn parse(
        &self,
        line: &[u8],
        chroms: &mut ChromTable,
    ) -> Result<(ChromId, u64, Interval), String>;
}

/// Intern `name`, explaining a failure against a closed table.
pub(crate) fn intern_chrom(chroms: &mut ChromTable, name: &[u8]) -> Result<ChromId, String> {
    if name.is_empty() {
        return Err("empty chromosome name".into());
    }
    chroms.intern(name).map_err(|e| {
        let place = match chroms.order_source() {
            OrderSource::GenomeFile(p) => format!("the genome file {}", p.display()),
            OrderSource::BuiltinGenome(id) => {
                format!("the built-in genome {id} (gloci genomes show {id} lists its names)")
            }
            OrderSource::VcfHeader(n) => format!("the ##contig lines of {n}"),
            OrderSource::Lexicographic | OrderSource::Karyotype | OrderSource::FirstAppearance => {
                "the genome".to_string()
            }
        };
        let hint = match chroms.order_source() {
            OrderSource::VcfHeader(_) => format!(
                "; add ##contig=<ID={}> to the header, or pass -g with a genome that lists it",
                e.0
            ),
            _ => String::new(),
        };
        format!("chromosome \"{}\" is not in {place}{hint}", e.0)
    })
}

/// Warn once about each contig a record named that the header does not
/// declare (htslib accepts such records). `loc` is where the first record
/// on it was found, when known.
pub(crate) fn warn_undeclared(chroms: &mut ChromTable, name: &str, loc: Option<Loc>) {
    for id in chroms.take_undeclared() {
        let chrom = String::from_utf8_lossy(chroms.name(id)).into_owned();
        let at = loc.map_or(String::new(), |l| format!(" {l}"));
        log::warn!(
            "{name}{at}: chromosome \"{chrom}\" is not in the ##contig lines; \
add ##contig=<ID={chrom}> to the header"
        );
    }
}

/// Parse a coordinate field, naming it in the error.
pub(crate) fn coord(field: &str, v: &[u8]) -> Result<u64, String> {
    gloci_core::parse_coord(v).map_err(|e| match e {
        CoordError::NotInteger(_) | CoordError::Negative(_) | CoordError::TooLarge(_) => {
            format!("{field} {e}")
        }
        _ => e.to_string(),
    })
}

pub(crate) struct TextReader<K: TextKind> {
    name: String,
    format: Format,
    compression: Compression,
    pub(crate) header: Header,
    diag: Diagnostics,
    policy: InvalidPolicy,
    src: Box<dyn BufRead + Send>,
    line: u64,
    pending: Option<Vec<u8>>,
    seen_record: bool,
    eof: bool,
    kind: K,
}

impl<K: TextKind> TextReader<K> {
    pub(crate) fn open(input: Input, policy: InvalidPolicy, kind: K) -> Result<Self, IoError> {
        let mut r = TextReader {
            name: input.name,
            format: input.format,
            compression: input.compression,
            header: Header::None,
            diag: Diagnostics::default(),
            policy,
            src: input.reader,
            line: 0,
            pending: None,
            seen_record: false,
            eof: false,
            kind,
        };
        let mut head = Vec::new();
        let mut line = Vec::new();
        loop {
            line.clear();
            if !r.read_line(&mut line, 0)? {
                r.eof = true;
                break;
            }
            match K::classify(&line, false) {
                Class::Skip => {
                    if line.is_empty() {
                        r.diag.record_blank(r.line);
                    }
                }
                Class::Header => {
                    head.extend_from_slice(&line);
                    head.push(b'\n');
                }
                Class::Fatal(msg) => return Err(r.parse_err(msg)),
                Class::Record => {
                    r.line -= 1;
                    r.pending = Some(line);
                    r.seen_record = true;
                    break;
                }
            }
        }
        if !head.is_empty() {
            r.header = Header::Lines(head);
        }
        Ok(r)
    }

    fn parse_err(&self, msg: String) -> IoError {
        IoError::Parse {
            path: self.name.clone(),
            loc: Loc::Line(self.line),
            msg,
        }
    }

    /// Append the next line (terminator stripped) to `out`; `start` is where
    /// the new line begins in `out`. Returns false at end of input.
    fn read_line(&mut self, out: &mut Vec<u8>, start: usize) -> Result<bool, IoError> {
        if let Some(p) = self.pending.take() {
            out.extend_from_slice(&p);
            self.line += 1;
            return Ok(true);
        }
        let n = self
            .src
            .read_until(b'\n', out)
            .map_err(|source| IoError::Read {
                path: self.name.clone(),
                source,
            })?;
        if n == 0 {
            return Ok(false);
        }
        self.line += 1;
        if out.last() == Some(&b'\n') {
            out.pop();
            if out.len() > start && out.last() == Some(&b'\r') {
                out.pop();
            }
        }
        Ok(true)
    }
}

impl<K: TextKind> RecordReader for TextReader<K> {
    fn name(&self) -> &str {
        &self.name
    }
    fn format(&self) -> Format {
        self.format
    }
    fn compression(&self) -> Compression {
        self.compression
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
        let mut n = 0;
        while n < max {
            if self.eof {
                return Ok(false);
            }
            let start = batch.buf.len();
            if !self.read_line(&mut batch.buf, start)? {
                self.eof = true;
                return Ok(false);
            }
            let line = &batch.buf[start..];
            match K::classify(line, self.seen_record) {
                Class::Skip | Class::Header => {
                    if line.is_empty() {
                        self.diag.record_blank(self.line);
                    }
                    batch.buf.truncate(start);
                    continue;
                }
                Class::Fatal(msg) => return Err(self.parse_err(msg)),
                Class::Record => {}
            }
            let len = line.len();
            let parsed = if len > u32::MAX as usize {
                Err("line is longer than 4 GiB".to_string())
            } else {
                self.kind.parse(line, chroms)
            };
            match parsed {
                Ok((chrom, key, span)) => {
                    batch.recs.push(RecordMeta {
                        chrom,
                        key,
                        span,
                        loc: self.line,
                        off: start,
                        len: len as u32,
                    });
                    self.seen_record = true;
                    n += 1;
                    if chroms.has_undeclared() {
                        warn_undeclared(chroms, &self.name, Some(Loc::Line(self.line)));
                    }
                }
                Err(msg) => {
                    batch.buf.truncate(start);
                    match self.policy {
                        InvalidPolicy::Fail => return Err(self.parse_err(msg)),
                        InvalidPolicy::Skip => self.diag.record(self.line, &msg),
                    }
                }
            }
        }
        Ok(true)
    }
}
