//! The bigWig and bigBed writer (D171), on bigtools.
//!
//! Every record is checked here before bigtools sees it, because bigtools checks little: only
//! the start of a bigBed record against its chromosome, nothing about a bigWig's overlaps
//! reaching the file, and a chromosome count it stores in 16 bits. A refused record names
//! why; the caller names where it came from. The file is written under a temporary name and
//! renamed into place by [`BbiWriter::finish`].
//!
//! bigtools runs on its own thread, fed batches of records through a bounded channel, so
//! memory stays flat in the input size. Its sections reach the file in input order whatever
//! the thread count, so the bytes do not depend on `--threads`.

use super::autosql::{BedLayout, autosql_field_count, bed_autosql, bed_layout, check_bed_line};
use crate::atomic::AtomicFile;
use crate::detect::{BIGBED_MAGIC, BIGWIG_MAGIC};
use crate::error::{IoError, UsageError};
use crate::output::{OutFormat, OutputSpec};
use bigtools::bed::bedparser::BedValueError;
use bigtools::beddata::BedParserStreamingIterator;
use bigtools::{BBIProcessError, BedEntry, BigBedWrite, BigWigWrite, InputSortType, Value};
use crossbeam_channel::{Receiver, Sender, bounded};
use gloci_core::Genome;
use std::collections::{HashMap, HashSet};
use std::io::{Read, Seek, SeekFrom, Write};
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::thread::JoinHandle;

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
    fn name(self) -> &'static str {
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
            "bigWig output is written only by gloci convert (from bedGraph, WIG or \
bigWig) and gloci genomecov --bedgraph or --bga; write bedGraph here (-o x.bedgraph) and convert it"
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
    /// The field count of the given autoSql, parsed by `create`.
    autosql_fields: Option<usize>,
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
        let autosql_fields = match &o.autosql {
            Some(text) => Some(
                autosql_field_count(text)
                    .map_err(|m| UsageError(format!("{}: {m}", o.autosql_name)))?,
            ),
            None => None,
        };
        let dest = path.display().to_string();
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
            autosql_fields,
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
                    let layout = bed_layout(raw).map_err(rec)?;
                    if let Some(n) = self.autosql_fields
                        && n != layout.columns
                    {
                        return Err(rec(format!(
                            "the autoSql declares {n} fields, but the records have {} \
columns; they must agree",
                            layout.columns
                        )));
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
        let runtime = if self.threads.get() == 1 {
            tokio::runtime::Builder::new_current_thread().build()
        } else {
            tokio::runtime::Builder::new_multi_thread()
                .worker_threads(self.threads.get())
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
        Err(PushError::Writer(match join(worker, &self.dest) {
            Err(e) => e,
            Ok(()) => IoError::Internal("the bigWig/bigBed writer stopped early".into()),
        }))
    }

    /// Ends the input, joins the writer, fixes definedFieldCount (bigBed), renames into place.
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
        join(worker, &self.dest)?;
        let af = self.af.take().expect("file present until finish");
        // bigtools writes the trailing magic last, inside a buffer whose flush errors on drop
        // are swallowed: a file without it is incomplete.
        finish_file(&af, self.kind, self.layout.map(|l| l.standard), &self.dest)?;
        af.commit().map_err(|source| IoError::Write {
            path: self.dest.clone(),
            source,
        })
    }
}

/// Ends the input of a worker and waits for it.
fn join(worker: Worker, dest: &str) -> Result<(), IoError> {
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
        Ok(Err(BBIProcessError::IoError(source))) => Err(IoError::Write {
            path: dest.to_string(),
            source,
        }),
        Ok(Err(e)) => Err(IoError::Internal(format!("bigtools: {e}"))),
        Err(_) => Err(IoError::Internal(
            "the bigWig/bigBed writer panicked".into(),
        )),
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
        let _ = join(worker, &self.dest);
    }
}

impl crate::writer::RecordWriter for BbiWriter {
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

/// Checks the finished temp file ends with its magic number, and for a bigBed writes the
/// number of standard columns into the header's definedFieldCount (bytes 34-35), in the byte
/// order the magic at the start of the file shows.
fn finish_file(
    af: &AtomicFile,
    kind: BbiKind,
    standard: Option<usize>,
    dest: &str,
) -> Result<(), IoError> {
    let write_err = |source| IoError::Write {
        path: dest.to_string(),
        source,
    };
    let incomplete = || IoError::Detect {
        path: dest.to_string(),
        msg: format!(
            "{dest} is incomplete (it does not end with the {} magic number); check free disk \
space",
            kind.name()
        ),
    };
    // The temp file is open for writing only: reopen it to read.
    let mut f = std::fs::OpenOptions::new()
        .read(true)
        .write(true)
        .open(af.temp_path())
        .map_err(write_err)?;
    let magic = match kind {
        BbiKind::BigWig => BIGWIG_MAGIC,
        BbiKind::BigBed => BIGBED_MAGIC,
    };
    let mut head = [0u8; 4];
    if f.read_exact(&mut head).is_err() {
        return Err(incomplete());
    }
    let little = u32::from_le_bytes(head) == magic;
    if !little && u32::from_be_bytes(head) != magic {
        return Err(incomplete());
    }
    let len = f.seek(SeekFrom::End(0)).map_err(write_err)?;
    let mut tail = [0u8; 4];
    if len < 8
        || f.seek(SeekFrom::Start(len - 4)).is_err()
        || f.read_exact(&mut tail).is_err()
        || tail != head
    {
        return Err(incomplete());
    }
    if let Some(standard) = standard {
        let n = u16::try_from(standard).unwrap_or(12);
        f.seek(SeekFrom::Start(34)).map_err(write_err)?;
        f.write_all(&if little {
            n.to_le_bytes()
        } else {
            n.to_be_bytes()
        })
        .map_err(write_err)?;
    }
    f.flush().map_err(write_err)
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
        let _ = join(real, &w.dest);
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
            Some(PushError::Writer(IoError::Write { source, .. })) => {
                assert_eq!(source.to_string(), "boom")
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
}
