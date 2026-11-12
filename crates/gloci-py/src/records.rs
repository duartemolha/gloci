//! Records (D89): a pipeline's output read back with gloci's own readers, so that iteration,
//! `count()` and the Arrow tables ([`crate::arrow`]) share one notion of a record (R1).
//!
//! The last step writes into a pipe that a host stage reads ([`reread`]), on its own thread and
//! never holding the GIL. How it reads depends on the step ([`embed::read_back`]):
//! - **A table** ([`Shape::Table`], [`Shape::HeaderRow`]): each line that is neither blank nor a
//!   `#` line is a row with no coordinates. A step that prints a header row (`stats jaccard`)
//!   has it name the columns; it is no row.
//! - **Records**, by their format as gloci's readers tell it from the content:
//!   - BED, bedGraph, GFF, GTF and VCF, through the format's reader: the coordinates are the
//!     record's footprint (a VCF record's by the step's `--vcf-span`). A line the reader does
//!     not take as a record, nor as a header, is a row with no coordinates (a BEDPE line whose
//!     first end is absent, `.  -1  -1`), so no line of the output is lost.
//!   - BCF records as their VCF text, BAM and SAM records as their SAM text: every record, an
//!     unmapped one with no footprint.
//!   - FASTA: one record per sequence, its `>` line and its sequence lines.
//!   - Anything else: a table with no header row.
//!
//! Blank lines, and `#`, `track` and `browser` lines, are never records.
//!
//! Iteration hands records to Python in batches of about [`BATCH_BYTES`] ([`IterSink`], read
//! by [`RecordIterator`]), and each [`Record`] is built under the GIL when it is asked for.

use crate::run::{Driver, spec_of};
use crossbeam_channel::{Receiver, SendTimeoutError, Sender, TryRecvError};
use gloci::embed::{CommandSpec, ReadBack, Shape};
use gloci_core::ChromTable;
use gloci_io::align::{Accept, AlignOptions, AlnRef, aln_view};
use gloci_io::{
    BcfText, Format, InputOptions, InvalidPolicy, IoError, Opened, ReaderOptions, RecordBatch,
    detect_format, next_batch_size, open_input_from, open_reader_from,
    open_reader_or_alignments_from,
};
use pyo3::basic::CompareOp;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::{PyIterator, PyString, PyTuple};
use std::io::{self, BufRead, Read};
use std::num::NonZeroUsize;
use std::sync::{Arc, Mutex, MutexGuard, PoisonError};
use std::time::Duration;

/// The text of the records a batch hands to Python.
pub(crate) const BATCH_BYTES: usize = 4 << 20;
/// The most records in a batch, whatever their size.
const BATCH_ROWS: usize = 1 << 16;
/// How long a full batch channel is waited on before the run's cancel token is looked at.
const SEND_POLL: Duration = Duration::from_millis(50);

/// How a record's fields are laid out, for the Arrow columns.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Layout {
    /// BED, bedGraph and other text with BED's first three columns.
    Bed,
    /// VCF, and BCF read as VCF text.
    Vcf,
    /// GFF3 and GTF.
    Gff,
    /// SAM text, and BAM read as SAM text.
    Sam,
    /// One record per sequence: its name and its sequence.
    Fasta,
    /// Rows of values with no coordinates.
    Table,
}

/// One record of the output, as a sink sees it.
pub(crate) struct Row<'a> {
    /// The record's text, without its line break (FASTA: its `>` line and sequence lines).
    pub line: &'a [u8],
    /// A FASTA record's name and sequence (its lines joined), its fields; else its fields are
    /// the line's tab-separated fields.
    pub fasta: Option<(&'a [u8], &'a [u8])>,
    /// The field that names the chromosome, when the record has one.
    pub chrom: Option<usize>,
    /// The footprint: `(start, end)`, 0-based, half-open.
    pub span: Option<(u64, u64)>,
}

impl Row<'_> {
    /// Call `f` with each field and its index.
    pub fn for_each_field(&self, mut f: impl FnMut(usize, &[u8])) {
        match self.fasta {
            Some((name, seq)) => {
                f(0, name);
                f(1, seq);
            }
            None => self
                .line
                .split(|&b| b == b'\t')
                .enumerate()
                .for_each(|(i, v)| f(i, v)),
        }
    }
}

/// What the output's records go to: a count, Python ([`IterSink`]) or an Arrow table.
pub(crate) trait RowSink {
    /// The layout of the records, said once before the first row.
    fn layout(&mut self, _layout: Layout) {}
    /// A line before the first row that may name the columns: a `#` line, or (with
    /// `header_row`) a table's header row.
    fn names(&mut self, _line: &[u8], _header_row: bool) {}
    /// A record. An error ends the read (the consumer is gone: a broken pipe).
    fn row(&mut self, row: &Row<'_>) -> io::Result<()>;
    /// Every record has been given.
    fn finish(&mut self) -> anyhow::Result<()> {
        Ok(())
    }
}

/// Counts the records (`count()`, R5): no Python object is made.
#[derive(Default)]
pub(crate) struct CountSink(pub u64);

impl RowSink for CountSink {
    fn row(&mut self, _row: &Row<'_>) -> io::Result<()> {
        self.0 += 1;
        Ok(())
    }
}

fn read_error(display: &str) -> impl Fn(io::Error) -> IoError + '_ {
    move |source| IoError::Read {
        path: display.to_string(),
        source,
    }
}

/// Read the output `src` (called `display` in messages), which the step `back` describes, into
/// `sink`. Runs on a host stage's thread: the readers' warnings were the step's to give, so
/// this thread's log records are dropped.
pub(crate) fn reread(
    src: Box<dyn Read + Send>,
    display: &str,
    back: ReadBack,
    sink: &mut dyn RowSink,
) -> anyhow::Result<()> {
    crate::logs::mute_this_thread();
    // Decompressed, and its format not detected yet: a table is read as lines whatever it is.
    let input = open_input_from(
        display,
        src,
        &InputOptions {
            threads: NonZeroUsize::MIN,
            format: Some(Format::Bed),
        },
    )?;
    let mut reader = input.reader;
    match back.shape {
        Shape::Table => read_table(&mut *reader, display, false, sink)?,
        Shape::HeaderRow => read_table(&mut *reader, display, true, sink)?,
        Shape::Records => {
            let head = reader.fill_buf().map_err(read_error(display))?;
            match detect_format(head) {
                Ok(f @ (Format::Bed | Format::BedGraph | Format::Gff3 | Format::Gtf)) => {
                    read_text(reader, display, f, back, sink)?
                }
                Ok(Format::Vcf) => read_text(reader, display, Format::Vcf, back, sink)?,
                Ok(Format::Bcf) => read_bcf(reader, display, back, sink)?,
                Ok(f @ (Format::Bam | Format::Sam | Format::Cram)) => {
                    read_alignments(reader, display, f, sink)?
                }
                Ok(Format::Fasta) => read_fasta(&mut *reader, display, sink)?,
                // Text gloci does not read as records (`groupby -g 1`, `getfasta --tab`).
                _ => read_table(&mut *reader, display, false, sink)?,
            }
        }
    }
    sink.finish()
}

/// The next line of `r` into `line`, without its line break or a CR before it; false at the
/// end.
fn next_line(r: &mut dyn BufRead, line: &mut Vec<u8>, display: &str) -> Result<bool, IoError> {
    line.clear();
    if r.read_until(b'\n', line).map_err(read_error(display))? == 0 {
        return Ok(false);
    }
    if line.last() == Some(&b'\n') {
        line.pop();
    }
    if line.last() == Some(&b'\r') {
        line.pop();
    }
    Ok(true)
}

fn is_blank(line: &[u8]) -> bool {
    line.iter().all(u8::is_ascii_whitespace)
}

/// A header line of a text format: `#`, `track` or `browser` (D14).
fn is_header_line(line: &[u8]) -> bool {
    let word =
        |w: &[u8]| line.starts_with(w) && matches!(line.get(w.len()), None | Some(b' ' | b'\t'));
    line.first() == Some(&b'#') || word(b"track") || word(b"browser")
}

/// A table: every line but blank and `#` lines is a row; with `header_row`, the first is the
/// header row.
fn read_table(
    r: &mut dyn BufRead,
    display: &str,
    header_row: bool,
    sink: &mut dyn RowSink,
) -> anyhow::Result<()> {
    sink.layout(Layout::Table);
    let mut line = Vec::new();
    let mut header_due = header_row;
    let mut rows = 0u64;
    while next_line(r, &mut line, display)? {
        if is_blank(&line) {
            continue;
        }
        if line[0] == b'#' {
            if rows == 0 {
                sink.names(&line, false);
            }
            continue;
        }
        if header_due {
            header_due = false;
            sink.names(&line, true);
            continue;
        }
        sink.row(&Row {
            line: &line,
            fasta: None,
            chrom: None,
            span: None,
        })?;
        rows += 1;
    }
    Ok(())
}

/// A FASTA output: one record per `>` line and the sequence lines after it.
fn read_fasta(r: &mut dyn BufRead, display: &str, sink: &mut dyn RowSink) -> anyhow::Result<()> {
    sink.layout(Layout::Fasta);
    let mut line = Vec::new();
    // The record being read: its text, where its name ends in it, and its sequence.
    let mut text: Vec<u8> = Vec::new();
    let mut name_end = 0;
    let mut seq: Vec<u8> = Vec::new();
    let mut open = false;
    let mut emit = |text: &[u8], name_end: usize, seq: &[u8]| {
        sink.row(&Row {
            line: text,
            fasta: Some((&text[1..name_end], seq)),
            chrom: None,
            span: None,
        })
    };
    while next_line(r, &mut line, display)? {
        if line.first() == Some(&b'>') {
            if open {
                emit(&text, name_end, &seq)?;
            }
            text.clear();
            seq.clear();
            text.extend_from_slice(&line);
            name_end = line.len();
            open = true;
        } else if open && !is_blank(&line) {
            text.push(b'\n');
            text.extend_from_slice(&line);
            seq.extend_from_slice(&line);
        }
    }
    if open {
        emit(&text, name_end, &seq)?;
    }
    Ok(())
}

fn lock<T>(m: &Mutex<T>) -> MutexGuard<'_, T> {
    m.lock().unwrap_or_else(PoisonError::into_inner)
}

/// A reader that keeps a copy of what it reads, so that the lines a record reader passes over
/// (headers, blank lines, lines it cannot read as records) are still in hand.
struct Tap {
    inner: Box<dyn BufRead + Send>,
    seen: Arc<Mutex<Vec<u8>>>,
}

impl Read for Tap {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        lock(&self.seen).extend_from_slice(&buf[..n]);
        Ok(n)
    }
}

/// The lines of the tapped bytes, numbered from 1 as the readers number them.
struct TapLines {
    /// Where the next line starts in the tapped bytes.
    pos: usize,
    /// The lines taken so far.
    line: u64,
}

impl TapLines {
    /// The next line of `buf`, without its line break or a CR before it: up to its line
    /// break, else the rest, which the caller knows to be a whole line (the reader has read
    /// past it, or the input has ended).
    fn next<'b>(&mut self, buf: &'b [u8]) -> Option<&'b [u8]> {
        let rest = buf.get(self.pos..).filter(|r| !r.is_empty())?;
        let (line, used) = match memchr::memchr(b'\n', rest) {
            Some(i) => (&rest[..i], i + 1),
            None => (rest, rest.len()),
        };
        self.pos += used;
        self.line += 1;
        Some(line.strip_suffix(b"\r").unwrap_or(line))
    }
}

/// A line a text reader passed over: nothing when blank or a header line (a `#` line before
/// the first record may name the columns), else a row with no coordinates.
fn passed_over(line: &[u8], rows: &mut u64, sink: &mut dyn RowSink) -> io::Result<()> {
    if is_blank(line) {
        return Ok(());
    }
    if is_header_line(line) {
        if *rows == 0 && line[0] == b'#' {
            sink.names(line, false);
        }
        return Ok(());
    }
    *rows += 1;
    sink.row(&Row {
        line,
        fasta: None,
        chrom: None,
        span: None,
    })
}

/// BED, bedGraph, GFF, GTF or VCF text through its reader, which skips the lines it cannot
/// read as records; those come from the tapped bytes, by their line numbers.
fn read_text(
    reader: Box<dyn BufRead + Send>,
    display: &str,
    format: Format,
    back: ReadBack,
    sink: &mut dyn RowSink,
) -> anyhow::Result<()> {
    sink.layout(match format {
        Format::Vcf => Layout::Vcf,
        Format::Gff3 | Format::Gtf => Layout::Gff,
        _ => Layout::Bed,
    });
    let seen = Arc::new(Mutex::new(Vec::new()));
    let tap = Tap {
        inner: reader,
        seen: seen.clone(),
    };
    let opts = ReaderOptions {
        threads: NonZeroUsize::MIN,
        format: Some(format),
        vcf_span: back.vcf_span,
        invalid: InvalidPolicy::Skip,
        chrom_map: None,
        accept: Accept::default(),
    };
    let mut r = open_reader_from(display, Box::new(tap), &opts)?;
    let mut chroms = ChromTable::first_appearance();
    let mut batch = RecordBatch::default();
    let mut lines = TapLines { pos: 0, line: 0 };
    let lost = |n: u64| IoError::Internal(format!("line {n} of {display} was not seen"));
    let mut rows = 0u64;
    let mut max = 1;
    loop {
        batch.clear();
        let more = r.read_batch(&mut chroms, &mut batch, max)?;
        {
            let buf = lock(&seen);
            for m in &batch.recs {
                while lines.line + 1 < m.loc {
                    let line = lines.next(&buf).ok_or_else(|| lost(lines.line + 1))?;
                    passed_over(line, &mut rows, sink)?;
                }
                let line = lines.next(&buf).ok_or_else(|| lost(m.loc))?;
                rows += 1;
                sink.row(&Row {
                    line,
                    fasta: None,
                    chrom: Some(0),
                    span: Some((m.span.start, m.span.end)),
                })?;
            }
        }
        // Drop the bytes of the lines taken.
        let mut buf = lock(&seen);
        buf.drain(..lines.pos);
        lines.pos = 0;
        drop(buf);
        max = next_batch_size(&batch, max);
        if !more {
            break;
        }
    }
    // The input has ended: the lines after the last record.
    let buf = lock(&seen);
    while let Some(line) = lines.next(&buf) {
        passed_over(line, &mut rows, sink)?;
    }
    Ok(())
}

/// BCF, as the VCF text of its records.
fn read_bcf(
    reader: Box<dyn BufRead + Send>,
    display: &str,
    back: ReadBack,
    sink: &mut dyn RowSink,
) -> anyhow::Result<()> {
    sink.layout(Layout::Vcf);
    let opts = ReaderOptions {
        threads: NonZeroUsize::MIN,
        format: Some(Format::Bcf),
        vcf_span: back.vcf_span,
        invalid: InvalidPolicy::Fail,
        chrom_map: None,
        accept: Accept::default(),
    };
    let mut r = open_reader_from(display, Box::new(reader), &opts)?;
    if let gloci_io::Header::Vcf(h) = r.header()
        && let Some(chrom_line) = h
            .raw
            .split(|&b| b == b'\n')
            .find(|l| l.starts_with(b"#CHROM"))
    {
        sink.names(chrom_line, false);
    }
    let mut text = BcfText::for_reader(&*r)?
        .ok_or_else(|| IoError::Internal(format!("{display} is not BCF")))?;
    let mut chroms = ChromTable::first_appearance();
    let mut batch = RecordBatch::default();
    let mut max = 1;
    loop {
        batch.clear();
        let more = r.read_batch(&mut chroms, &mut batch, max)?;
        for (i, m) in batch.recs.iter().enumerate() {
            let line = text.render(batch.raw(i), m.loc)?;
            sink.row(&Row {
                line,
                fasta: None,
                chrom: Some(0),
                span: Some((m.span.start, m.span.end)),
            })?;
        }
        max = next_batch_size(&batch, max);
        if !more {
            break;
        }
    }
    Ok(())
}

/// BAM, SAM (or CRAM), as the SAM text of every record.
fn read_alignments(
    reader: Box<dyn BufRead + Send>,
    display: &str,
    format: Format,
    sink: &mut dyn RowSink,
) -> anyhow::Result<()> {
    sink.layout(Layout::Sam);
    let opts = ReaderOptions {
        threads: NonZeroUsize::MIN,
        format: Some(format),
        vcf_span: Default::default(),
        invalid: InvalidPolicy::Fail,
        chrom_map: None,
        accept: Accept {
            alignments: Some(AlignOptions {
                full_records: true,
                ..AlignOptions::default()
            }),
            convert_only: false,
        },
    };
    let Opened::Alignments(mut src) = open_reader_or_alignments_from(display, reader, &opts)?
    else {
        return Err(IoError::Internal(format!("{display} is not read as alignments")).into());
    };
    let h = src.header().clone();
    let layout = src.layout();
    let mut raw = Vec::new();
    let mut line = Vec::new();
    while let Some(n) = src.next_record(&mut raw)? {
        let bad = |msg: String| IoError::Parse {
            path: display.to_string(),
            loc: gloci_io::Loc::Record(n),
            msg,
        };
        let a = aln_view(layout, &raw, &h).map_err(bad)?;
        line.clear();
        match a {
            AlnRef::Bam(b) => b.write_sam(&h, &mut line).map_err(bad)?,
            AlnRef::Sam(_) => line.extend_from_slice(&raw),
        }
        let ref_len = |i: usize| h.refs.get(i).map_or(u64::MAX, |r| r.1);
        let span = a
            .ref_id()
            .and_then(|i| a.span(ref_len(i)))
            .map(|(iv, _)| (iv.start, iv.end));
        sink.row(&Row {
            line: &line,
            fasta: None,
            chrom: a.ref_id().map(|_| 2),
            span,
        })?;
        raw.clear();
    }
    Ok(())
}

// ---- Iteration -------------------------------------------------------------------------------

/// Where a record's parts are in its batch's text.
struct RowMeta {
    line: (usize, usize),
    /// A FASTA record's name and sequence.
    fasta: Option<[(usize, usize); 2]>,
    chrom: Option<usize>,
    span: Option<(u64, u64)>,
}

/// Records handed to Python at once.
#[derive(Default)]
pub(crate) struct Batch {
    text: Vec<u8>,
    rows: Vec<RowMeta>,
}

/// Sends the records to the iterator in batches of about [`BATCH_BYTES`].
pub(crate) struct IterSink {
    tx: Sender<Batch>,
    batch: Batch,
}

impl IterSink {
    pub(crate) fn new(tx: Sender<Batch>) -> IterSink {
        IterSink {
            tx,
            batch: Batch::default(),
        }
    }

    /// Send the batch, waiting while the iterator has one it has not taken; a cancelled run,
    /// or an iterator that is gone, ends the read.
    fn flush(&mut self) -> io::Result<()> {
        if self.batch.rows.is_empty() {
            return Ok(());
        }
        let mut batch = std::mem::take(&mut self.batch);
        loop {
            match self.tx.send_timeout(batch, SEND_POLL) {
                Ok(()) => return Ok(()),
                Err(SendTimeoutError::Timeout(b)) => {
                    if gloci_io::cancel::is_cancelled() {
                        return Err(gloci_io::cancel::io_error());
                    }
                    batch = b;
                }
                Err(SendTimeoutError::Disconnected(_)) => {
                    return Err(io::ErrorKind::BrokenPipe.into());
                }
            }
        }
    }
}

impl RowSink for IterSink {
    fn row(&mut self, row: &Row<'_>) -> io::Result<()> {
        let b = &mut self.batch;
        let mut put = |bytes: &[u8]| {
            let at = b.text.len();
            b.text.extend_from_slice(bytes);
            (at, b.text.len())
        };
        let line = put(row.line);
        let fasta = row.fasta.map(|(name, seq)| [put(name), put(seq)]);
        b.rows.push(RowMeta {
            line,
            fasta,
            chrom: row.chrom,
            span: row.span,
        });
        if b.text.len() >= BATCH_BYTES || b.rows.len() >= BATCH_ROWS {
            self.flush()?;
        }
        Ok(())
    }

    fn finish(&mut self) -> anyhow::Result<()> {
        Ok(self.flush()?)
    }
}

/// A Python `str` of `bytes`, which gloci's output has as UTF-8 (anything else is replaced).
fn py_str<'py>(py: Python<'py>, bytes: &[u8]) -> Bound<'py, PyString> {
    match std::str::from_utf8(bytes) {
        Ok(s) => PyString::new(py, s),
        Err(_) => PyString::new(py, &String::from_utf8_lossy(bytes)),
    }
}

/// One record of a pipeline's output (D89): its coordinates as gloci's readers read them back,
/// `None` where the output has none, its fields, and its text. A sequence over its fields;
/// `str()` is its line.
#[pyclass(frozen, sequence, module = "gloci", name = "Record")]
pub(crate) struct Record {
    /// The chromosome, or None.
    #[pyo3(get)]
    chrom: Option<Py<PyString>>,
    /// The 0-based start of its footprint, or None.
    #[pyo3(get)]
    start: Option<u64>,
    /// The end of its footprint (half-open), or None.
    #[pyo3(get)]
    end: Option<u64>,
    /// Its fields: the line's tab-separated values (FASTA: the name and the sequence).
    #[pyo3(get)]
    fields: Py<PyTuple>,
    /// Its text, as the output has it, without the line break.
    #[pyo3(get)]
    line: Py<PyString>,
}

impl Record {
    fn of(py: Python<'_>, text: &[u8], m: &RowMeta) -> PyResult<Record> {
        let part = |(a, b): (usize, usize)| &text[a..b];
        let line = part(m.line);
        let fields: Vec<Bound<'_, PyString>> = match m.fasta {
            Some([name, seq]) => vec![py_str(py, part(name)), py_str(py, part(seq))],
            None => line.split(|&b| b == b'\t').map(|f| py_str(py, f)).collect(),
        };
        let chrom = m
            .chrom
            .and_then(|i| fields.get(i))
            .map(|s| s.clone().unbind());
        Ok(Record {
            chrom,
            start: m.span.map(|s| s.0),
            end: m.span.map(|s| s.1),
            fields: PyTuple::new(py, &fields)?.unbind(),
            line: py_str(py, line).unbind(),
        })
    }

    /// The line, for a source that reads records (D88).
    pub(crate) fn line_str<'py>(&self, py: Python<'py>) -> &Bound<'py, PyString> {
        self.line.bind(py)
    }

    fn key<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        PyTuple::new(
            py,
            [
                self.chrom
                    .as_ref()
                    .map_or_else(|| py.None(), |c| c.clone_ref(py).into_any()),
                self.start.into_pyobject(py)?.into_any().unbind(),
                self.end.into_pyobject(py)?.into_any().unbind(),
                self.fields.clone_ref(py).into_any(),
                self.line.clone_ref(py).into_any(),
            ],
        )
    }
}

#[pymethods]
impl Record {
    fn __len__(&self, py: Python<'_>) -> usize {
        self.fields.bind(py).len()
    }

    fn __getitem__<'py>(
        &self,
        py: Python<'py>,
        index: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        self.fields.bind(py).as_any().get_item(index)
    }

    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyIterator>> {
        self.fields.bind(py).as_any().try_iter()
    }

    fn __contains__(&self, py: Python<'_>, value: &Bound<'_, PyAny>) -> PyResult<bool> {
        self.fields.bind(py).contains(value)
    }

    fn __str__(&self, py: Python<'_>) -> Py<PyString> {
        self.line.clone_ref(py)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let chrom = match &self.chrom {
            Some(c) => c.bind(py).repr()?.to_string(),
            None => "None".into(),
        };
        let num = |v: Option<u64>| v.map_or_else(|| "None".into(), |v| v.to_string());
        Ok(format!(
            "Record(chrom={chrom}, start={}, end={}, line={})",
            num(self.start),
            num(self.end),
            self.line.bind(py).repr()?
        ))
    }

    fn __richcmp__(
        &self,
        py: Python<'_>,
        other: &Bound<'_, PyAny>,
        op: CompareOp,
    ) -> PyResult<Py<PyAny>> {
        let Ok(other) = other.cast::<Record>() else {
            return Ok(py.NotImplemented());
        };
        let same = self.key(py)?.eq(other.get().key(py)?)?;
        Ok(match op {
            CompareOp::Eq => same.into_pyobject(py)?.to_owned().into_any().unbind(),
            CompareOp::Ne => (!same).into_pyobject(py)?.to_owned().into_any().unbind(),
            _ => py.NotImplemented(),
        })
    }

    fn __hash__(&self, py: Python<'_>) -> PyResult<isize> {
        self.key(py)?.hash()
    }
}

/// What a [`RecordIterator`] holds while the run goes on.
struct Live {
    /// The batches the host stage sends; None once it has ended, or when the run is being
    /// cancelled.
    rx: Option<Receiver<Batch>>,
    batch: Batch,
    next: usize,
    driver: Driver,
    /// Each stage's command, for usage messages.
    specs: Vec<Option<&'static CommandSpec>>,
}

/// The records of a pipeline's output, as the run makes them (D89, R2): `iter(pipeline)`.
///
/// The steps run while it is read. Breaking out of a loop over it, or dropping it, cancels the
/// steps and joins their threads (with the GIL released), and their temporary files are
/// removed. An error of the run is raised once the records before it have been read (D90).
#[pyclass(module = "gloci._native", name = "RecordIterator")]
pub(crate) struct RecordIterator {
    live: Mutex<Option<Live>>,
}

impl RecordIterator {
    pub(crate) fn new(
        rx: Receiver<Batch>,
        driver: Driver,
        stages: &[Vec<String>],
    ) -> RecordIterator {
        RecordIterator {
            live: Mutex::new(Some(Live {
                rx: Some(rx),
                batch: Batch::default(),
                next: 0,
                driver,
                specs: stages.iter().map(|a| spec_of(a)).collect(),
            })),
        }
    }
}

impl Live {
    /// The next record, driving the run (D91) until one comes; None at the end. The run's
    /// error is raised once every record before it was given.
    fn next(&mut self, py: Python<'_>) -> PyResult<Option<Record>> {
        loop {
            if let Some(m) = self.batch.rows.get(self.next) {
                self.next += 1;
                return Record::of(py, &self.batch.text, m).map(Some);
            }
            if self.driver.interrupted() {
                // A signal or a log handler stopped the run: no more records.
                self.rx = None;
            }
            if let Some(rx) = &self.rx {
                match rx.try_recv() {
                    Ok(batch) => {
                        self.batch = batch;
                        self.next = 0;
                        continue;
                    }
                    Err(TryRecvError::Empty) => {}
                    Err(TryRecvError::Disconnected) => self.rx = None,
                }
            }
            if self.rx.is_none() && self.driver.ended() {
                let specs = &self.specs;
                self.driver
                    .settle(py, |i| specs.get(i).copied().flatten())?;
                return Ok(None);
            }
            self.driver.turn(py, self.rx.as_ref());
        }
    }
}

impl Drop for Live {
    fn drop(&mut self) {
        // The host stage may wait to send a batch: let it go before the run is cancelled and
        // joined (the driver's drop).
        self.rx = None;
    }
}

#[pymethods]
impl RecordIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&self, py: Python<'_>) -> PyResult<Option<Record>> {
        crate::guarded(py, || {
            let mut slot = match self.live.try_lock() {
                Ok(slot) => slot,
                Err(std::sync::TryLockError::Poisoned(p)) => p.into_inner(),
                Err(std::sync::TryLockError::WouldBlock) => {
                    return Err(PyRuntimeError::new_err(
                        "this gloci iterator is being read by another thread",
                    ));
                }
            };
            let Some(live) = slot.as_mut() else {
                return Ok(None);
            };
            let next = live.next(py);
            if !matches!(next, Ok(Some(_))) {
                // The end, or the run's error: the run is over, and so is the iteration.
                let done = slot.take();
                drop(slot);
                drop(done);
            }
            next
        })
    }

    /// Stop the run now, as dropping the iterator does: the steps are cancelled and joined.
    fn close(&self) -> PyResult<()> {
        let mut slot = match self.live.try_lock() {
            Ok(slot) => slot,
            Err(std::sync::TryLockError::Poisoned(p)) => p.into_inner(),
            Err(std::sync::TryLockError::WouldBlock) => {
                return Err(PyRuntimeError::new_err(
                    "this gloci iterator is being read by another thread",
                ));
            }
        };
        let done = slot.take();
        drop(slot);
        drop(done);
        Ok(())
    }
}
