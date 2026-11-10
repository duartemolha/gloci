//! Record batches and the reader interface shared by every command.

use crate::chrom_map::{ChromMap, Renamer};
use crate::detect::{Compression, Format};
use crate::diag::{Diagnostics, InvalidPolicy};
use crate::error::{IoError, Loc};
use crate::input::{Input, InputOptions, is_regular_file, open_input, open_input_from};
use gloci_core::{ChromId, ChromTable, CoordError, Interval, Strand, VcfSpan};
use std::io::{BufRead, Read};
use std::num::NonZeroUsize;
use std::sync::Arc;

/// One record's position data; its raw bytes live in `RecordBatch::buf`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RecordMeta {
    pub chrom: ChromId,
    /// The position the file is sorted by (BED start, GFF start-1, VCF POS-1).
    pub key: u64,
    /// The record's footprint.
    pub span: Interval,
    /// Strand (D54): BED column 6, GFF column 7; unknown otherwise.
    pub strand: Strand,
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
    /// The header of a BAM, SAM or CRAM input.
    Alignment(Arc<crate::align::AlnHeader>),
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
    ///
    /// # Contract
    /// `max` must be greater than 0 (`max > 0`).
    fn read_batch(
        &mut self,
        chroms: &mut ChromTable,
        batch: &mut RecordBatch,
        max: usize,
    ) -> Result<bool, IoError>;
    fn diagnostics(&self) -> &Diagnostics;
    /// Whether a sweep must read this input to its end, so that every record
    /// of it is validated whether or not the other inputs reach it: true of
    /// text read from the start of a file or stream, false of a reader that
    /// already limits itself to chosen chromosomes through an index.
    fn consume_to_end(&self) -> bool {
        false
    }
    /// Whether the records are whole, so that a writer may write them. False of an
    /// alignment reader that keeps only the fields comparisons need (D115): `--sort` then
    /// sorts slim records into its temporary file and reads them back slim.
    fn whole_records(&self) -> bool {
        true
    }
}

#[derive(Debug, Clone)]
pub struct ReaderOptions {
    pub threads: NonZeroUsize,
    pub format: Option<Format>,
    pub vcf_span: VcfSpan,
    pub invalid: InvalidPolicy,
    /// `--chrom-map` (D62): renames chromosomes as the records are read.
    pub chrom_map: Option<Arc<ChromMap>>,
    /// What the command reads beyond interval formats (D100): the input gate.
    pub accept: crate::align::Accept,
}

/// What [`open_reader_or_alignments`] found.
pub enum Opened {
    /// An interval format, read as records.
    Records(Box<dyn RecordReader>),
    /// BAM or SAM, read as every-record alignments (D120).
    Alignments(Box<dyn crate::align::AlnSource>),
}

/// Open `path` as [`open_reader`] does, except that an alignment file is opened as an
/// [`AlnSource`](crate::align::AlnSource) when `o.accept.alignments` is set: for a command
/// (`view`) that handles alignments as whole records and other formats as intervals, and
/// must not open a stream twice to find out which it has.
pub fn open_reader_or_alignments(path: &str, o: &ReaderOptions) -> Result<Opened, IoError> {
    let input = open_input(path, &input_options(o))?;
    if input.format.is_alignment() && o.accept.alignments.is_some() {
        return crate::align::source_from_input(input, o).map(Opened::Alignments);
    }
    finish_reader(path, dispatch(input, o)?).map(Opened::Records)
}

pub fn open_reader(path: &str, o: &ReaderOptions) -> Result<Box<dyn RecordReader>, IoError> {
    finish_reader(path, dispatch(open_input(path, &input_options(o))?, o)?)
}

fn finish_reader(path: &str, r: Box<dyn RecordReader>) -> Result<Box<dyn RecordReader>, IoError> {
    // A sweep reads a plain file to its end (see `RecordReader::consume_to_end`).
    if r.compression() == Compression::None
        && crate::seekable::range_format(r.format())
        && is_regular_file(path)
    {
        return Ok(Box::new(crate::seekable::DrainedToEnd(r)));
    }
    Ok(r)
}

/// With `--chrom-map`, refuse a BGZF text file that has an index when two of
/// its chromosomes would get one name, as a run through the index refuses
/// it: a relation command calls this on each input it opens, so that its
/// answer does not depend on how many threads are used. An index that cannot
/// be read is not this check's business.
pub fn check_map_against_index(
    path: &str,
    r: &dyn RecordReader,
    o: &ReaderOptions,
) -> Result<(), IoError> {
    let Some(map) = &o.chrom_map else {
        return Ok(());
    };
    if r.compression() != Compression::Bgzf
        || r.format() == Format::Bcf
        || !crate::is_regular_file(path)
    {
        return Ok(());
    }
    let data = std::path::Path::new(path);
    let Some(info) = crate::index::find_index(data) else {
        return Ok(());
    };
    if let Ok(names) = crate::index::index_reference_names(data, &info) {
        map.check_file_names(path, names.iter().map(|n| &n[..]))?;
    }
    Ok(())
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
    crate::align::gate(&input.name, input.format, &o.accept)?;
    match input.format {
        Format::Bed | Format::BedGraph => Ok(Box::new(crate::bed::BedReader::open_bed(
            input,
            o.invalid,
            o.chrom_map.clone(),
        )?)),
        Format::Wig => Ok(Box::new(crate::wig::WigReader::open(
            input,
            o.invalid,
            o.chrom_map.clone(),
        )?)),
        Format::Gff3 | Format::Gtf => Ok(Box::new(crate::gff::GffReader::open_gff(
            input,
            o.invalid,
            o.chrom_map.clone(),
        )?)),
        Format::Vcf => Ok(Box::new(crate::vcf::VcfReader::open_vcf(
            input,
            o.invalid,
            o.vcf_span,
            o.chrom_map.clone(),
        )?)),
        Format::Bcf => Ok(Box::new(crate::bcf::BcfReader::open_bcf(
            input,
            o.invalid,
            o.vcf_span,
            o.chrom_map.as_deref(),
        )?)),
        Format::Bam if o.accept.alignments.is_some() => {
            Ok(Box::new(crate::align::BamReader::open_bam(
                input,
                o.invalid,
                o.chrom_map.as_deref(),
                o.accept.alignments.clone().unwrap_or_default(),
            )?))
        }
        Format::Psl => Ok(Box::new(crate::psl::PslReader::open_psl(
            input,
            o.invalid,
            o.chrom_map.clone(),
        )?)),
        Format::Rmsk => Ok(Box::new(crate::rmsk::RmskReader::open_rmsk(
            input,
            o.invalid,
            o.chrom_map.clone(),
        )?)),
        Format::Sam if o.accept.alignments.is_some() => {
            Ok(Box::new(crate::align::SamReader::open_sam(
                input,
                o.invalid,
                o.chrom_map.as_deref(),
                o.accept.alignments.clone().unwrap_or_default(),
            )?))
        }
        // `gate` has returned an error for each of these.
        Format::Bam | Format::Cram | Format::Sam | Format::Fasta | Format::Bedpe => {
            Err(IoError::Detect {
                path: input.name.clone(),
                msg: format!(
                    "{} is {}, which is not readable here",
                    input.name,
                    input.format.name()
                ),
            })
        }
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
    fn parse(&self, line: &[u8], chroms: &mut ChromTable) -> Result<Parsed, String>;
}

/// A parsed text record: chromosome, sort key, footprint, strand.
pub(crate) type Parsed = (ChromId, u64, Interval, Strand);

/// Intern `name`, explaining a failure against a closed table.
pub(crate) fn intern_chrom(chroms: &mut ChromTable, name: &[u8]) -> Result<ChromId, String> {
    if name.is_empty() {
        return Err("empty chromosome name".into());
    }
    chroms
        .intern(name)
        .map_err(|e| crate::naming::unknown_chrom_err(chroms.order_source(), &e.0))
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
    /// An error found after records of the same batch were read: the batch
    /// returns those records first, and the next call raises it, so the
    /// earliest problem in the file is the one reported whatever the batch
    /// boundaries.
    deferred: Option<IoError>,
    kind: K,
    /// `--chrom-map`: rewrites each record's chromosome column as it is read.
    renamer: Option<Renamer>,
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
            deferred: None,
            kind,
            renamer: None,
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

    /// A reader of the lines at `src`, which holds whole records starting at
    /// line `first_line` of the file called `name`, after `header`: the
    /// reader a byte range of a plain file gets (see `seekable`).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn over_range(
        name: &str,
        format: Format,
        header: Header,
        policy: InvalidPolicy,
        src: Box<dyn BufRead + Send>,
        first_line: u64,
        kind: K,
        map: Option<Arc<ChromMap>>,
    ) -> Self {
        TextReader {
            name: name.to_string(),
            format,
            compression: Compression::None,
            header,
            diag: Diagnostics::default(),
            policy,
            src,
            line: first_line - 1,
            pending: None,
            seen_record: true,
            eof: false,
            deferred: None,
            kind,
            renamer: map.map(Renamer::new),
        }
    }

    /// The line parser, for a reader that keeps state in it.
    pub(crate) fn kind(&self) -> &K {
        &self.kind
    }

    /// Rename chromosomes with `map` from now on (D62).
    pub(crate) fn set_chrom_map(&mut self, map: Option<Arc<ChromMap>>) {
        self.renamer = map.map(Renamer::new);
    }

    fn parse_err(&self, msg: String) -> IoError {
        IoError::Parse {
            path: self.name.clone(),
            loc: Loc::Line(self.line),
            msg,
        }
    }

    /// Report `e`, found after `n` records of this batch: at once when there
    /// are none, else on the next call (see `deferred`).
    fn fail(&mut self, n: usize, e: IoError) -> Result<bool, IoError> {
        if n == 0 {
            return Err(e);
        }
        self.deferred = Some(e);
        Ok(true)
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
        }
        if out.len() > start && out.last() == Some(&b'\r') {
            out.pop();
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
        debug_assert!(max > 0, "read_batch max must be > 0");
        if let Some(e) = self.deferred.take() {
            return Err(e);
        }
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
                Class::Fatal(msg) => {
                    batch.buf.truncate(start);
                    return self.fail(n, self.parse_err(msg));
                }
                Class::Record => {}
            }
            if let Some(r) = &mut self.renamer {
                r.apply(&mut batch.buf, start);
            }
            let line = &batch.buf[start..];
            let len = line.len();
            let parsed = if len > u32::MAX as usize {
                Err("line is longer than 4 GiB".to_string())
            } else {
                self.kind.parse(line, chroms)
            };
            match parsed {
                Ok((chrom, key, span, strand)) => {
                    batch.recs.push(RecordMeta {
                        chrom,
                        key,
                        span,
                        strand,
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
                        InvalidPolicy::Fail => return self.fail(n, self.parse_err(msg)),
                        InvalidPolicy::Skip => self.diag.record(self.line, &msg),
                    }
                }
            }
        }
        Ok(true)
    }
}
