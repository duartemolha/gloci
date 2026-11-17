//! The bigWig and bigBed writer (D171), on bigtools.
//!
//! Every record is checked here before bigtools sees it, because bigtools checks little: only
//! the start of a bigBed record against its chromosome, nothing about a bigWig's overlaps
//! reaching the file, and a chromosome count it stores in 16 bits. A refused record names
//! why; the caller names where it came from. The file is written under a temporary name and
//! renamed into place by [`BbiWriter::finish`], once [`super::check`] has read it back.
//!
//! bigtools runs on its own thread, fed batches of records through a bounded channel, so
//! memory stays flat in the input size. Its sections reach the file in input order whatever
//! the thread count, so the bytes do not depend on `--threads`.
//!
//! bigtools spools each chromosome to temporary files in `$TMPDIR`, about 10 open at once per
//! chromosome with records, and loses some of their errors. So `create` raises the soft
//! open-file limit to the hard one and checks that `$TMPDIR` takes a file, a stopped or
//! panicking bigtools is reported with those likely causes, and the file is checked before it
//! is renamed into place: a short or inconsistent file is never kept (D171).

use super::autosql::{
    BedLayout, autosql_field_count, autosql_standard_columns, bed_autosql, bed_layout,
    check_bed_line, standard_size,
};
use super::check;
use crate::atomic::{AtomicFile, panic_note, set_thread_panic_note};
use crate::detect::{BIGBED_MAGIC, BIGWIG_MAGIC};
use crate::error::{IoError, UsageError};
use crate::output::{OutFormat, OutputSpec};
use bigtools::bed::bedparser::BedValueError;
use bigtools::beddata::BedParserStreamingIterator;
use bigtools::{BBIProcessError, BedEntry, BigBedWrite, BigWigWrite, InputSortType, Value};
use crossbeam_channel::{Receiver, Sender, bounded};
use gloci_core::Genome;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::thread::JoinHandle;
use std::time::Instant;

/// Records per message to the writer thread.
const BATCH: usize = 4096;
/// Messages the channel holds before `push` waits for the writer.
const QUEUE: usize = 4;
/// A bigWig or bigBed stores its chromosome count in 16 bits.
const MAX_CHROMS: usize = 65_535;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BbiKind {
    BigWig,
    BigBed,
}

impl BbiKind {
    pub(crate) fn name(self) -> &'static str {
        match self {
            BbiKind::BigWig => "bigWig",
            BbiKind::BigBed => "bigBed",
        }
    }
}

/// Chromosome sizes and how to name their source in errors ("-g hg38", "the header of reads.bam").
pub struct BbiSizes {
    pub genome: Genome,
    pub source: String,
}

pub struct BbiWriteOptions {
    pub kind: BbiKind,
    pub sizes: BbiSizes,
    pub threads: NonZeroUsize,
    /// bigBed: the autoSql text to store, as given.
    pub autosql: Option<String>,
    /// Names the autoSql in errors ("bad.as", or the option that gave it).
    pub autosql_name: String,
    /// bigBed: the definedFieldCount to keep, a bigBed input's, instead of the one the autoSql
    /// or the first record gives; capped at the standard columns the first record has (D171).
    pub defined_field_count: Option<usize>,
}

#[derive(Debug)]
pub enum PushError {
    /// The record is refused; the caller names where it came from.
    Record(String),
    /// The writer stopped; its error.
    Writer(IoError),
}

impl std::fmt::Display for PushError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PushError::Record(m) => f.write_str(m),
            PushError::Writer(e) => e.fmt(f),
        }
    }
}

impl std::error::Error for PushError {}

/// What `--output-format`/`-o` of a command that cannot write a bigWig or bigBed is told.
pub fn bbi_output_refusal(f: OutFormat) -> String {
    match f {
        OutFormat::BigBed => "bigBed output is written only by gloci convert (from BED or \
bigBed); write BED here (-o x.bed) and convert it with gloci convert x.bed -g GENOME -o x.bb"
            .to_string(),
        OutFormat::BigWig => {
            "bigWig output is written only by gloci convert (from bedGraph, a 4-column \
BED, WIG or bigWig) and gloci genomecov --bedgraph or --bga; write bedGraph here (-o x.bedgraph) \
and convert it"
                .to_string()
        }
        other => unreachable!(
            "bbi_output_refusal is for bigWig and bigBed, not {}",
            other.name()
        ),
    }
}

enum Msg<T> {
    Batch(Vec<(String, T)>),
    Abort,
}

/// The records a writer thread reads, batch by batch.
struct Feed<T> {
    rx: Receiver<Msg<T>>,
    cur: std::vec::IntoIter<(String, T)>,
    aborted: bool,
}

impl<T> Iterator for Feed<T> {
    type Item = Result<(String, T), BedValueError>;
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(item) = self.cur.next() {
                return Some(Ok(item));
            }
            if self.aborted {
                return None;
            }
            match self.rx.recv() {
                Ok(Msg::Batch(b)) => self.cur = b.into_iter(),
                Ok(Msg::Abort) => {
                    self.aborted = true;
                    return Some(Err(BedValueError::InvalidInput("write aborted".into())));
                }
                Err(_) => return None,
            }
        }
    }
}

type WorkerResult = Result<(), BBIProcessError<BedValueError>>;

/// The writer thread and the way to feed it.
enum Worker {
    Wig(Sender<Msg<Value>>, JoinHandle<WorkerResult>),
    Bed(Sender<Msg<BedEntry>>, JoinHandle<WorkerResult>),
}

/// Where a bigBed's standard-column count comes from, when not from its first record.
#[derive(Clone, Copy, Debug)]
enum Standard {
    /// A bigBed input's definedFieldCount: at most the first record's standard columns.
    Kept(usize),
    /// The given autoSql's: the first record's columns must be valid up to it.
    Declared(usize),
}

/// The checks and bookkeeping of the chromosome being written.
struct Cursor {
    chrom: String,
    start: u64,
    end: u64,
}

pub struct BbiWriter {
    kind: BbiKind,
    dest: String,
    af: Option<AtomicFile>,
    sizes: HashMap<Vec<u8>, Option<u64>>,
    source: String,
    threads: NonZeroUsize,
    autosql: Option<String>,
    autosql_name: String,
    /// The field count of the given autoSql, parsed by `create`.
    autosql_fields: Option<usize>,
    /// How many leading columns are standard BED columns, when not the first record's own:
    /// kept from a bigBed input, or declared by the autoSql.
    standard: Option<Standard>,
    layout: Option<BedLayout>,
    worker: Option<Worker>,
    /// The writer thread stopped with an error: nothing more is accepted.
    failed: bool,
    wig: Vec<(String, Value)>,
    bed: Vec<(String, BedEntry)>,
    records: u64,
    cursor: Option<Cursor>,
    done: HashSet<String>,
    /// A record that covers at least one base has been accepted.
    covers_a_base: bool,
}

impl BbiWriter {
    /// Usage error for stdout, a device or pipe (`spec.direct`), an in-memory sink, or a spec that is not bigWig/bigBed.
    pub fn create(spec: &OutputSpec, o: BbiWriteOptions) -> Result<BbiWriter, IoError> {
        let want = match o.kind {
            BbiKind::BigWig => OutFormat::BigWig,
            BbiKind::BigBed => OutFormat::BigBed,
        };
        if spec.format != Some(want) {
            return Err(IoError::Internal(format!(
                "a {} writer was given the output format {}",
                o.kind.name(),
                spec.format.map_or("none", OutFormat::name)
            )));
        }
        let path: PathBuf = match &spec.path {
            // An in-memory sink is always `direct`.
            Some(p) if !spec.direct => p.clone(),
            _ => {
                return Err(UsageError(format!(
                    "{} output must be a regular file (-o FILE.{}): it is written under a \
temporary name and renamed, so standard output, a device, a pipe or an in-memory sink cannot \
take it",
                    o.kind.name(),
                    if o.kind == BbiKind::BigWig {
                        "bw"
                    } else {
                        "bb"
                    }
                ))
                .into());
            }
        };
        // Parse the autoSql before the output exists: a bad one is the user's usage error.
        let usage = |m| UsageError(format!("{}: {m}", o.autosql_name));
        let (autosql_fields, declared) = match &o.autosql {
            Some(text) => (
                Some(autosql_field_count(text).map_err(usage)?),
                Some(autosql_standard_columns(text).map_err(usage)?),
            ),
            None => (None, None),
        };
        let standard = match (o.defined_field_count, declared) {
            (Some(n), _) => Some(Standard::Kept(n)),
            (None, Some(n)) => Some(Standard::Declared(n)),
            (None, None) => None,
        };
        let dest = path.display().to_string();
        raise_open_file_limit();
        check_temp_dir(&dest)?;
        let af = AtomicFile::create(&path).map_err(|source| IoError::Create {
            path: dest.clone(),
            source,
        })?;
        let sizes = (0..o.sizes.genome.len())
            .map(|i| (o.sizes.genome.name(i).to_vec(), o.sizes.genome.length(i)))
            .collect();
        Ok(BbiWriter {
            kind: o.kind,
            dest,
            af: Some(af),
            sizes,
            source: o.sizes.source,
            threads: o.threads,
            autosql: o.autosql,
            autosql_name: o.autosql_name,
            autosql_fields,
            standard,
            layout: None,
            worker: None,
            failed: false,
            wig: Vec::new(),
            bed: Vec::new(),
            records: 0,
            cursor: None,
            done: HashSet::new(),
            covers_a_base: false,
        })
    }

    /// One bedGraph (BigWig) or BED (BigBed) line, no newline. After a `Writer` error nothing
    /// more is accepted: every later call returns an error and starts no writer.
    pub fn push(&mut self, raw: &[u8]) -> Result<(), PushError> {
        if self.failed {
            return Err(PushError::Writer(already_failed()));
        }
        let r = self.push_checked(raw);
        if matches!(r, Err(PushError::Writer(_))) {
            self.failed = true;
        }
        r
    }

    fn push_checked(&mut self, raw: &[u8]) -> Result<(), PushError> {
        let line = std::str::from_utf8(raw).map_err(|e| {
            PushError::Record(format!(
                "the record is not valid UTF-8 (at byte {}); a bigWig or bigBed holds text in \
UTF-8",
                e.valid_up_to() + 1
            ))
        })?;
        let rec = PushError::Record;
        let f: Vec<&str> = line.split('\t').collect();
        let mut value = 0.0f32;
        let mut first_layout = None;
        match self.kind {
            BbiKind::BigWig => {
                if f.len() != 4 {
                    return Err(rec(format!(
                        "a bigWig record needs 4 columns (chrom, start, end, value), this one \
has {}",
                        f.len()
                    )));
                }
                value = parse_value(f[3]).map_err(rec)?;
            }
            BbiKind::BigBed => {
                if let Some(layout) = &self.layout {
                    check_bed_line(raw, layout).map_err(rec)?;
                } else {
                    let mut layout = bed_layout(raw).map_err(rec)?;
                    if let Some(n) = self.autosql_fields
                        && n != layout.columns
                    {
                        return Err(rec(format!(
                            "the autoSql declares {n} fields, but the records have {} \
columns; they must agree",
                            layout.columns
                        )));
                    }
                    match self.standard {
                        None => {}
                        Some(Standard::Kept(n)) => {
                            layout.standard = standard_size(n.min(layout.standard));
                        }
                        Some(Standard::Declared(n)) => {
                            layout.standard = standard_size(n.min(layout.columns));
                            check_bed_line(raw, &layout).map_err(|e| {
                                rec(format!(
                                    "{e}, but {} declares the first {} columns as BED's \
standard ones",
                                    self.autosql_name, layout.standard
                                ))
                            })?;
                        }
                    }
                    // Kept only once the whole record is accepted.
                    first_layout = Some(layout);
                }
            }
        }
        let (chrom, start, end) = (f[0], pos("start", f[1])?, pos("end", f[2])?);
        if end < start {
            return Err(rec(format!("the end {end} is before the start {start}")));
        }
        if start == end && self.kind == BbiKind::BigWig {
            return Err(rec(
                "a bigWig cannot hold a zero-length record (start equals end)".into(),
            ));
        }
        if start == end && start == 0 {
            return Err(rec(
                "a bigBed cannot hold a zero-length record at position 0 (bigtools' reader \
rejects it)"
                    .into(),
            ));
        }
        match self.sizes.get(chrom.as_bytes()) {
            None => {
                return Err(rec(format!(
                    "chromosome {chrom} is not in {}; every chromosome of a bigWig or bigBed \
needs a length",
                    self.source
                )));
            }
            Some(None) => {
                return Err(rec(format!(
                    "{} gives no length for chromosome {chrom}, which a bigWig or bigBed needs",
                    self.source
                )));
            }
            Some(Some(len)) if *len > u64::from(u32::MAX) => {
                return Err(rec(format!(
                    "chromosome {chrom} is {len} bases long, more than the 32-bit positions of \
a bigWig or bigBed can hold"
                )));
            }
            Some(Some(len)) if end > *len => {
                return Err(rec(format!(
                    "the end {end} is past the end of chromosome {chrom} ({len} bases, from {})",
                    self.source
                )));
            }
            // bigtools refuses a record that starts at the chromosome's end, after the fact.
            Some(Some(len)) if start == *len => {
                return Err(rec(format!(
                    "a zero-length record at {start} is at the very end of chromosome {chrom} \
({len} bases, from {}); a bigBed record must start inside its chromosome",
                    self.source
                )));
            }
            Some(Some(_)) => {}
        }
        match &mut self.cursor {
            Some(c) if c.chrom == chrom => {
                if start < c.start {
                    return Err(rec(format!(
                        "the start {start} is before the previous record's start {} on \
{chrom}; a bigWig or bigBed needs each chromosome sorted by start (gloci sort)",
                        c.start
                    )));
                }
                if self.kind == BbiKind::BigWig && start < c.end {
                    return Err(rec(format!(
                        "{chrom}:{start}-{end} overlaps the previous record {chrom}:{}-{}; a \
bigWig cannot hold overlapping records",
                        c.start, c.end
                    )));
                }
                c.start = start;
                c.end = end;
            }
            cur => {
                if self.done.contains(chrom) {
                    return Err(rec(format!(
                        "the records of {chrom} are not together: they resume after another \
chromosome; a bigWig or bigBed needs each chromosome's records in one run (gloci sort)"
                    )));
                }
                let seen = self.done.len() + usize::from(cur.is_some());
                if seen >= MAX_CHROMS {
                    return Err(rec(
                        "more than 65,535 chromosomes have records; a bigWig or bigBed holds at \
most 65,535"
                            .into(),
                    ));
                }
                if let Some(old) = cur.take() {
                    self.done.insert(old.chrom);
                }
                *cur = Some(Cursor {
                    chrom: chrom.to_string(),
                    start,
                    end,
                });
            }
        }
        // Checked above: both positions fit 32 bits, as the chromosome's length does.
        let (start, end) = (start as u32, end as u32);
        if first_layout.is_some() {
            self.layout = first_layout;
        }
        if self.worker.is_none() {
            self.start_worker()?;
        }
        match self.kind {
            BbiKind::BigWig => {
                self.wig
                    .push((chrom.to_string(), Value { start, end, value }));
                if self.wig.len() >= BATCH {
                    self.send_batch()?;
                }
            }
            BbiKind::BigBed => {
                // Everything after the end column, tab separated, as bigtools stores it.
                let rest = f[3..].join("\t");
                self.bed
                    .push((chrom.to_string(), BedEntry { start, end, rest }));
                if self.bed.len() >= BATCH {
                    self.send_batch()?;
                }
            }
        }
        self.covers_a_base |= end > start;
        self.records += 1;
        Ok(())
    }

    fn start_worker(&mut self) -> Result<(), PushError> {
        let werr = |e: std::io::Error| {
            PushError::Writer(IoError::Write {
                path: self.dest.clone(),
                source: e,
            })
        };
        let file = self
            .af
            .as_ref()
            .expect("file present until finish")
            .try_clone_file()
            .map_err(werr)?;
        // A panic in bigtools (its temporary files failing, mostly) is reported with the likely
        // causes, on its runtime's threads as on this one, instead of as gloci's internal error.
        let note = unfinished(
            &self.dest,
            &format!("the {} writer stopped", self.kind.name()),
        )
        .to_string();
        let runtime = if self.threads.get() == 1 {
            tokio::runtime::Builder::new_current_thread().build()
        } else {
            let note = note.clone();
            tokio::runtime::Builder::new_multi_thread()
                .worker_threads(self.threads.get())
                .on_thread_start(move || set_thread_panic_note(note.clone()))
                .build()
        }
        .map_err(werr)?;
        let chrom_sizes: HashMap<String, u32> = self
            .sizes
            .iter()
            .filter_map(|(n, l)| {
                let name = String::from_utf8(n.clone()).ok()?;
                Some((name, u32::try_from((*l)?).ok()?))
            })
            .collect();
        let spawn = |body: Box<dyn FnOnce() -> WorkerResult + Send>| {
            std::thread::Builder::new()
                .name("gloci-bbi-writer".into())
                .stack_size(crate::WORKER_STACK)
                .spawn(body)
                .map_err(werr)
        };
        self.worker = Some(match self.kind {
            BbiKind::BigWig => {
                let (tx, rx) = bounded(QUEUE);
                let handle = spawn(Box::new(move || {
                    let _note = panic_note(note);
                    let feed = Feed {
                        rx,
                        cur: Vec::new().into_iter(),
                        aborted: false,
                    };
                    let mut w = BigWigWrite::new(file, chrom_sizes);
                    w.options.input_sort_type = InputSortType::START;
                    w.write(BedParserStreamingIterator::wrap_iter(feed, true), runtime)
                }))?;
                Worker::Wig(tx, handle)
            }
            BbiKind::BigBed => {
                let layout = self
                    .layout
                    .expect("the layout is set with the first record");
                let text = self.autosql.clone().unwrap_or_else(|| bed_autosql(&layout));
                let (tx, rx) = bounded(QUEUE);
                let handle = spawn(Box::new(move || {
                    let _note = panic_note(note);
                    let feed = Feed {
                        rx,
                        cur: Vec::new().into_iter(),
                        aborted: false,
                    };
                    let mut w = BigBedWrite::new(file, chrom_sizes);
                    w.options.input_sort_type = InputSortType::START;
                    w.autosql = Some(text);
                    w.write(BedParserStreamingIterator::wrap_iter(feed, true), runtime)
                }))?;
                Worker::Bed(tx, handle)
            }
        });
        Ok(())
    }

    /// Sends the pending batch. A send that fails means the writer stopped: its error.
    fn send_batch(&mut self) -> Result<(), PushError> {
        let sent = match self.worker.as_ref().expect("worker started") {
            Worker::Wig(tx, _) => tx.send(Msg::Batch(std::mem::take(&mut self.wig))).is_ok(),
            Worker::Bed(tx, _) => tx.send(Msg::Batch(std::mem::take(&mut self.bed))).is_ok(),
        };
        if sent {
            return Ok(());
        }
        let worker = self.worker.take().expect("worker started");
        self.failed = true;
        Err(PushError::Writer(
            match join(worker, self.kind, &self.dest) {
                Err(e) => e,
                Ok(()) => unfinished(
                    &self.dest,
                    &format!(
                        "the {} writer stopped before its input ended",
                        self.kind.name()
                    ),
                ),
            },
        ))
    }

    /// Ends the input, joins the writer, checks the file and fixes definedFieldCount (bigBed),
    /// renames into place.
    pub fn finish(mut self) -> Result<(), IoError> {
        if self.failed {
            return Err(already_failed());
        }
        if self.records == 0 {
            return Err(IoError::Detect {
                path: self.dest.clone(),
                msg: format!(
                    "no records to write; a {} needs at least one record",
                    self.kind.name()
                ),
            });
        }
        if !self.covers_a_base {
            // bigtools' writer never finishes when no record covers a base, as it has no
            // zoom data to write: stop it (the drop aborts) instead of waiting for it.
            return Err(IoError::Detect {
                path: self.dest.clone(),
                msg: format!(
                    "a {} cannot be written when every record has zero length (start equals \
end); the writer needs at least one record that covers a base",
                    self.kind.name()
                ),
            });
        }
        if let Err(e) = self.send_batch() {
            return Err(match e {
                PushError::Writer(e) => e,
                PushError::Record(m) => IoError::Internal(m),
            });
        }
        let worker = self.worker.take().expect("worker started");
        join(worker, self.kind, &self.dest)?;
        let af = self.af.take().expect("file present until finish");
        finish_file(
            &af,
            self.kind,
            self.layout.map(|l| l.standard),
            self.records,
            &self.dest,
        )?;
        af.commit().map_err(|source| IoError::Write {
            path: self.dest.clone(),
            source,
        })
    }
}

/// Ends the input of a worker and waits for it. Its I/O errors and its panics (when a panic
/// unwinds rather than ending the process) name the likely causes.
fn join(worker: Worker, kind: BbiKind, dest: &str) -> Result<(), IoError> {
    let handle = match worker {
        Worker::Wig(tx, h) => {
            drop(tx);
            h
        }
        Worker::Bed(tx, h) => {
            drop(tx);
            h
        }
    };
    match handle.join() {
        Ok(Ok(())) => Ok(()),
        Ok(Err(BBIProcessError::IoError(source))) => {
            Err(unfinished(dest, &format!("bigtools failed: {source}")))
        }
        // gloci checks every record first, so bigtools refusing one is gloci's bug.
        Ok(Err(e)) => Err(IoError::Internal(format!("bigtools: {e}"))),
        Err(payload) => {
            let why = payload
                .downcast_ref::<&str>()
                .map(|s| (*s).to_string())
                .or_else(|| payload.downcast_ref::<String>().cloned())
                .unwrap_or_else(|| "a panic with no message".into());
            Err(unfinished(
                dest,
                &format!("the {} writer stopped ({why})", kind.name()),
            ))
        }
    }
}

impl Drop for BbiWriter {
    /// Dropping without `finish` aborts the write; the temp file goes with the `AtomicFile`.
    fn drop(&mut self) {
        let Some(worker) = self.worker.take() else {
            return;
        };
        match &worker {
            Worker::Wig(tx, _) => {
                let _ = tx.send(Msg::Abort);
            }
            Worker::Bed(tx, _) => {
                let _ = tx.send(Msg::Abort);
            }
        }
        let _ = join(worker, self.kind, &self.dest);
    }
}

impl crate::writer::RecordWriter for BbiWriter {
    /// `io::ErrorKind::InvalidData` means a refused record, a data error (exit 1) naming it:
    /// callers rely on that kind. Any other kind carries the writer's own `IoError`.
    fn write_record(&mut self, raw: &[u8]) -> std::io::Result<()> {
        self.push(raw).map_err(|e| match e {
            PushError::Record(m) => std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("record {} of the output: {m}", self.records + 1),
            ),
            PushError::Writer(e) => std::io::Error::other(e),
        })
    }
    fn finish(self: Box<Self>) -> Result<(), IoError> {
        BbiWriter::finish(*self)
    }
}

fn pos(name: &str, v: &str) -> Result<u64, PushError> {
    crate::record::coord(name, v.as_bytes()).map_err(PushError::Record)
}

/// A bedGraph value as the 32-bit float a bigWig stores.
fn parse_value(s: &str) -> Result<f32, String> {
    let wide: f64 = s
        .parse()
        .map_err(|_| format!("the value '{s}' is not a number"))?;
    if !wide.is_finite() {
        return Err(format!(
            "the value '{s}' is not a finite number; a bigWig cannot hold nan or infinity"
        ));
    }
    match s.parse::<f32>() {
        Ok(v) if v.is_finite() => Ok(v),
        _ => Err(format!(
            "the value '{s}' is outside what a 32-bit float can hold (about 3.4e38)"
        )),
    }
}

/// The error of a write bigtools could not finish: `what` went wrong, then the likely causes.
fn unfinished(dest: &str, what: &str) -> IoError {
    IoError::Detect {
        path: dest.to_string(),
        msg: format!("cannot write {dest}: {what}; {}", resources_hint()),
    }
}

/// The likely causes of a write bigtools could not finish: its temporary files.
fn resources_hint() -> String {
    #[cfg(unix)]
    let limit = match open_file_limits() {
        Some((soft, hard)) => format!(
            "ulimit -n: soft {}, hard {}",
            limit_text(soft),
            limit_text(hard)
        ),
        None => "ulimit -n".to_string(),
    };
    #[cfg(not(unix))]
    let limit = "the system's limit".to_string();
    format!(
        "bigtools keeps about 10 temporary files open per chromosome, in $TMPDIR ({}): check the \
open-file limit ({limit}) and free disk space in $TMPDIR and for the output",
        std::env::temp_dir().display()
    )
}

/// bigtools creates its temporary files in `$TMPDIR`: one that cannot take a file fails here,
/// naming it, before anything is written.
fn check_temp_dir(dest: &str) -> Result<(), IoError> {
    tempfile::tempfile().map(drop).map_err(|e| IoError::Detect {
        path: dest.to_string(),
        msg: format!(
            "cannot write {dest}: bigtools writes its temporary files to $TMPDIR ({}), and none \
can be created there ({e}); set TMPDIR to a writable directory with free disk space",
            std::env::temp_dir().display()
        ),
    })
}

/// The soft and hard open-file limits.
#[cfg(unix)]
fn open_file_limits() -> Option<(libc::rlim_t, libc::rlim_t)> {
    let mut lim = libc::rlimit {
        rlim_cur: 0,
        rlim_max: 0,
    };
    // SAFETY: getrlimit only writes the struct it is given.
    let ok = unsafe { libc::getrlimit(libc::RLIMIT_NOFILE, &mut lim) } == 0;
    ok.then_some((lim.rlim_cur, lim.rlim_max))
}

#[cfg(unix)]
fn limit_text(n: libc::rlim_t) -> String {
    if n == libc::RLIM_INFINITY {
        "unlimited".into()
    } else {
        n.to_string()
    }
}

/// Raises the soft open-file limit to the hard one, once per process: bigtools keeps about 10
/// temporary files open per chromosome with records while it writes (D171).
fn raise_open_file_limit() {
    #[cfg(unix)]
    {
        static ONCE: std::sync::Once = std::sync::Once::new();
        ONCE.call_once(|| {
            let Some((soft, hard)) = open_file_limits() else {
                return;
            };
            if soft >= hard {
                return;
            }
            let set = |n| {
                let lim = libc::rlimit {
                    rlim_cur: n,
                    rlim_max: hard,
                };
                // SAFETY: setrlimit only reads the struct it is given.
                unsafe { libc::setrlimit(libc::RLIMIT_NOFILE, &lim) == 0 }
            };
            match raise_soft_limit(soft, hard, set) {
                Some(n) => log::debug!(
                    "raised the open-file limit from {soft} to {n} (hard limit {}) for the \
bigWig/bigBed writer",
                    limit_text(hard)
                ),
                None => log::debug!(
                    "cannot raise the open-file limit above {soft} (hard limit {})",
                    limit_text(hard)
                ),
            }
        });
    }
}

/// The largest soft limit above `soft`, at most `hard`, that `set` accepts, and the one in
/// force when it returns; `None` when it accepts none. macOS refuses an unlimited hard limit
/// as the soft one, so `OPEN_MAX` is tried next, then a search between the largest value
/// accepted and the smallest refused.
#[cfg(unix)]
fn raise_soft_limit(
    soft: libc::rlim_t,
    hard: libc::rlim_t,
    mut set: impl FnMut(libc::rlim_t) -> bool,
) -> Option<libc::rlim_t> {
    /// macOS's `OPEN_MAX`.
    const OPEN_MAX: libc::rlim_t = 10_240;
    if set(hard) {
        return Some(hard);
    }
    let (mut lo, mut hi) = (soft, hard);
    if lo < OPEN_MAX && OPEN_MAX < hi {
        if set(OPEN_MAX) {
            lo = OPEN_MAX;
        } else {
            hi = OPEN_MAX;
        }
    }
    while hi - lo > 1 {
        let mid = lo + (hi - lo) / 2;
        if set(mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (lo > soft).then_some(lo)
}

/// Checks the finished temp file before the rename: the magic number at both ends, then its
/// structure and record count ([`check`]). For a bigBed, then writes the number of standard
/// columns into the header's definedFieldCount.
fn finish_file(
    af: &AtomicFile,
    kind: BbiKind,
    standard: Option<usize>,
    records: u64,
    dest: &str,
) -> Result<(), IoError> {
    let write_err = |source| IoError::Write {
        path: dest.to_string(),
        source,
    };
    // The temp file is open for writing only: reopen it to read.
    let mut f = std::fs::OpenOptions::new()
        .read(true)
        .write(true)
        .open(af.temp_path())
        .map_err(write_err)?;
    // bigtools writes the trailing magic last, inside a buffer whose flush errors on drop are
    // swallowed: a file without it is incomplete.
    let little = byte_order(&mut f, kind).ok_or_else(|| {
        unfinished(
            dest,
            &format!(
                "{dest} is incomplete (it does not end with the {} magic number)",
                kind.name()
            ),
        )
    })?;
    let t = Instant::now();
    let checked = check::check_written(af.temp_path(), kind, records).map_err(|why| {
        unfinished(
            dest,
            &format!(
                "the {} bigtools wrote is incomplete or inconsistent, so it was not kept: {why}",
                kind.name()
            ),
        )
    })?;
    log::debug!(
        "checked {dest} before renaming it into place: {} sections, {} zoom levels, in {:.2?}",
        checked.sections,
        checked.zoom_levels,
        t.elapsed()
    );
    if let Some(standard) = standard {
        write_defined_field_count(&mut f, little, standard).map_err(write_err)?;
    }
    f.flush().map_err(write_err)
}

/// Whether `f` is little-endian, by the magic number at its start, which must also end it;
/// `None` when either is missing.
fn byte_order(f: &mut File, kind: BbiKind) -> Option<bool> {
    let magic = match kind {
        BbiKind::BigWig => BIGWIG_MAGIC,
        BbiKind::BigBed => BIGBED_MAGIC,
    };
    let mut head = [0u8; 4];
    f.seek(SeekFrom::Start(0)).ok()?;
    f.read_exact(&mut head).ok()?;
    let little = u32::from_le_bytes(head) == magic;
    if !little && u32::from_be_bytes(head) != magic {
        return None;
    }
    let len = f.seek(SeekFrom::End(0)).ok()?;
    let mut tail = [0u8; 4];
    if len < 8 {
        return None;
    }
    f.seek(SeekFrom::Start(len - 4)).ok()?;
    f.read_exact(&mut tail).ok()?;
    (tail == head).then_some(little)
}

/// Writes `standard` into the header's definedFieldCount (bytes 34-35) in the file's byte order.
fn write_defined_field_count(f: &mut File, little: bool, standard: usize) -> std::io::Result<()> {
    let n = u16::try_from(standard).unwrap_or(12);
    f.seek(SeekFrom::Start(34))?;
    f.write_all(&if little {
        n.to_le_bytes()
    } else {
        n.to_be_bytes()
    })
}

fn already_failed() -> IoError {
    IoError::Internal(
        "the bigWig/bigBed writer already failed; no more records are accepted".into(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use gloci_core::Genome;

    fn writer(dir: &std::path::Path) -> BbiWriter {
        let genome = Genome::new(vec![(b"chr1".to_vec(), Some(1_000_000))]).unwrap();
        let spec = crate::output::parse_track_output_spec(Some(&dir.join("p.bw")), None).unwrap();
        BbiWriter::create(
            &spec,
            BbiWriteOptions {
                kind: BbiKind::BigWig,
                sizes: BbiSizes {
                    genome,
                    source: "-g test".into(),
                },
                threads: NonZeroUsize::new(1).unwrap(),
                autosql: None,
                autosql_name: "autosql".into(),
                defined_field_count: None,
            },
        )
        .unwrap()
    }

    /// After a writer error, `push` and `finish` return an error and never start a second
    /// writer on the temp file or panic.
    #[test]
    fn a_failed_writer_is_poisoned() {
        let dir = tempfile::tempdir().unwrap();
        let mut w = writer(dir.path());
        w.push(b"chr1\t0\t1\t1").unwrap();
        // Replace the running writer with one that has already stopped with an error.
        let real = w.worker.take().unwrap();
        if let Worker::Wig(tx, _) = &real {
            let _ = tx.send(Msg::Abort);
        }
        let _ = join(real, w.kind, &w.dest);
        let (tx, rx) = bounded(QUEUE);
        drop(rx);
        let handle =
            std::thread::spawn(|| Err(BBIProcessError::IoError(std::io::Error::other("boom"))));
        w.worker = Some(Worker::Wig(tx, handle));
        let mut got = None;
        for i in 1..=BATCH as u64 + 1 {
            if let Err(e) = w.push(format!("chr1\t{}\t{}\t1", i * 2, i * 2 + 1).as_bytes()) {
                got = Some(e);
                break;
            }
        }
        match got {
            Some(PushError::Writer(IoError::Detect { msg, .. })) => {
                assert!(msg.contains("bigtools failed: boom"), "{msg}");
                assert!(msg.contains("ulimit -n"), "{msg}");
            }
            other => panic!("expected the writer's error, got {other:?}"),
        }
        assert!(w.worker.is_none());
        assert!(matches!(
            w.push(b"chr1\t900000\t900001\t1"),
            Err(PushError::Writer(IoError::Internal(_)))
        ));
        assert!(w.worker.is_none(), "no second writer was started");
        assert!(matches!(w.finish(), Err(IoError::Internal(_))));
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
    }

    /// A crafted temp file: `head` magic, 40 bytes of body, `tail` magic.
    fn crafted(head: [u8; 4], tail: &[u8]) -> (tempfile::TempDir, AtomicFile) {
        let dir = tempfile::tempdir().unwrap();
        let af = AtomicFile::create(&dir.path().join("c.out")).unwrap();
        let mut bytes = head.to_vec();
        bytes.extend_from_slice(&[0u8; 40]);
        bytes.extend_from_slice(tail);
        std::fs::write(af.temp_path(), &bytes).unwrap();
        (dir, af)
    }

    fn order_of(kind: BbiKind, head: [u8; 4], tail: &[u8]) -> Option<bool> {
        let (_dir, af) = crafted(head, tail);
        byte_order(&mut File::open(af.temp_path()).unwrap(), kind)
    }

    #[test]
    fn a_complete_little_endian_bigwig_is_accepted() {
        let m = BIGWIG_MAGIC.to_le_bytes();
        assert_eq!(order_of(BbiKind::BigWig, m, &m), Some(true));
    }

    #[test]
    fn a_tail_cut_by_one_byte_is_incomplete() {
        let m = BIGWIG_MAGIC.to_le_bytes();
        assert_eq!(order_of(BbiKind::BigWig, m, &m[..3]), None);
        let (_dir, af) = crafted(m, &m[..3]);
        match finish_file(&af, BbiKind::BigWig, None, 1, "c.out") {
            Err(IoError::Detect { msg, .. }) => {
                assert!(msg.contains("c.out is incomplete"), "{msg}");
                assert!(msg.contains("free disk space"), "{msg}");
                assert!(msg.contains("ulimit -n"), "{msg}");
                assert!(msg.contains("$TMPDIR"), "{msg}");
            }
            other => panic!("expected the incomplete error, got {other:?}"),
        }
    }

    #[test]
    fn a_big_endian_bigbed_gets_a_big_endian_field_count() {
        let m = BIGBED_MAGIC.to_be_bytes();
        assert_eq!(order_of(BbiKind::BigBed, m, &m), Some(false));
        let (_dir, af) = crafted(m, &m);
        let mut f = std::fs::OpenOptions::new()
            .read(true)
            .write(true)
            .open(af.temp_path())
            .unwrap();
        write_defined_field_count(&mut f, false, 6).unwrap();
        let after = std::fs::read(af.temp_path()).unwrap();
        assert_eq!(u16::from_be_bytes([after[34], after[35]]), 6);
    }

    #[test]
    fn a_bigwig_magic_under_a_bigbed_writer_is_refused() {
        let m = BIGWIG_MAGIC.to_le_bytes();
        assert_eq!(order_of(BbiKind::BigBed, m, &m), None);
    }

    /// A crafted file passes the magic check but not the structure check, and is not kept.
    #[test]
    fn a_file_with_only_its_magic_numbers_is_not_kept() {
        let m = BIGWIG_MAGIC.to_le_bytes();
        let (dir, af) = crafted(m, &m);
        match finish_file(&af, BbiKind::BigWig, None, 1, "c.out") {
            Err(IoError::Detect { msg, .. }) => {
                assert!(msg.contains("incomplete or inconsistent"), "{msg}");
                assert!(msg.contains("ulimit -n"), "{msg}");
            }
            other => panic!("expected the check's error, got {other:?}"),
        }
        drop(af);
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
    }

    /// The soft limit `raise_soft_limit` settles on, and the one `set` left in force.
    #[cfg(unix)]
    fn raised(
        soft: libc::rlim_t,
        hard: libc::rlim_t,
        cap: libc::rlim_t,
    ) -> (Option<libc::rlim_t>, libc::rlim_t) {
        let mut now = soft;
        let got = raise_soft_limit(soft, hard, |n| {
            let ok = n <= cap;
            if ok {
                now = n;
            }
            ok
        });
        (got, now)
    }

    #[cfg(unix)]
    #[test]
    fn the_soft_open_file_limit_is_raised_as_far_as_allowed() {
        // Linux: the hard limit is accepted.
        assert_eq!(raised(1024, 4096, 4096), (Some(4096), 4096));
        // macOS: an unlimited hard limit is refused; OPEN_MAX, then the largest accepted.
        assert_eq!(
            raised(256, libc::RLIM_INFINITY, 10_240),
            (Some(10_240), 10_240)
        );
        assert_eq!(
            raised(256, libc::RLIM_INFINITY, 245_760),
            (Some(245_760), 245_760)
        );
        assert_eq!(raised(256, libc::RLIM_INFINITY, 5000), (Some(5000), 5000));
        // Nothing above the soft limit is accepted.
        assert_eq!(raised(256, 4096, 256), (None, 256));
    }

    #[test]
    fn a_panicking_writer_names_the_likely_causes() {
        let (tx, rx) = bounded(QUEUE);
        drop(rx);
        let handle = std::thread::spawn(|| -> WorkerResult { panic!("Couldn't send section.") });
        match join(Worker::Wig(tx, handle), BbiKind::BigWig, "x.bw") {
            Err(IoError::Detect { msg, .. }) => {
                assert!(
                    msg.starts_with("cannot write x.bw: the bigWig writer stopped"),
                    "{msg}"
                );
                assert!(msg.contains("Couldn't send section."), "{msg}");
                assert!(msg.contains("ulimit -n"), "{msg}");
                assert!(msg.contains("$TMPDIR"), "{msg}");
                assert!(!msg.contains("internal error"), "{msg}");
            }
            other => panic!("expected the writer's error, got {other:?}"),
        }
    }
}
