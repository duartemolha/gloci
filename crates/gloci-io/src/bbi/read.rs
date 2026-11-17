//! The bigWig and bigBed reader (D171).
//!
//! A bbi file carries its own chromosome list and an index over its data, so the reader does
//! not scan from the start: it orders the file's chromosomes by the command's `ChromTable`
//! and walks each one through the index, records as stored (starts ascending). A bigWig
//! record is a bedGraph row; a bigBed record is a BED line.

use crate::atomic::panic_note;
use crate::chrom_map::ChromMap;
use crate::detect::{Compression, Format, bbi_format};
use crate::diag::Diagnostics;
use crate::error::{IoError, Loc};
use crate::record::{Header, ReaderOptions, RecordBatch, RecordMeta, RecordReader, intern_chrom};
use bigtools::{
    BBIReadError, BedEntry, BigBedIntervalIter, BigBedRead, BigWigIntervalIter, BigWigRead, Value,
};
use gloci_core::{ChromId, ChromTable, Interval, Strand};
use std::cmp::Ordering;
use std::collections::VecDeque;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::sync::Arc;

/// A bbi file that is open and between chromosomes.
enum Source {
    Wig(BigWigRead<File>),
    Bed(BigBedRead<File>),
}

/// A chromosome being walked.
enum Cursor {
    Wig(BigWigIntervalIter<File, BigWigRead<File>>),
    Bed(BigBedIntervalIter<File, BigBedRead<File>>),
}

enum Item {
    Wig(Value),
    Bed(BedEntry),
}

impl Cursor {
    fn next(&mut self) -> Option<Result<Item, BBIReadError>> {
        match self {
            Cursor::Wig(it) => it.next().map(|r| r.map(Item::Wig)),
            Cursor::Bed(it) => it.next().map(|r| r.map(Item::Bed)),
        }
    }

    fn finish(self) -> Source {
        match self {
            Cursor::Wig(it) => Source::Wig(it.into()),
            Cursor::Bed(it) => Source::Bed(it.into()),
        }
    }
}

/// The chromosome being walked: the file's name for it, the name records carry, and its id
/// once its first record has been seen (a chromosome with no records is never interned).
struct Walk {
    out_name: Vec<u8>,
    id: Option<ChromId>,
}

pub(crate) struct BbiReader {
    path: String,
    format: Format,
    renamer: Option<Arc<ChromMap>>,
    header: Header,
    diag: Diagnostics,
    /// The file's chromosomes, in its own list's order, until the first batch orders them.
    chroms: Vec<String>,
    /// The chromosomes still to read, in the command's order.
    queue: Option<VecDeque<String>>,
    source: Option<Source>,
    cursor: Option<(Cursor, Walk)>,
    /// Records yielded so far.
    count: u64,
}

fn what(format: Format) -> &'static str {
    format.name()
}

/// The error for a file that is not (or no longer) a valid bbi file.
fn invalid(path: &str, format: Format, detail: &str) -> IoError {
    IoError::Detect {
        path: path.to_string(),
        msg: format!(
            "{path} is not a valid {} file (it may be truncated): {detail}",
            what(format)
        ),
    }
}

fn not_this_format(path: &str, format: Format) -> IoError {
    IoError::Detect {
        path: path.to_string(),
        msg: format!(
            "{path} is not a {} file (it does not start with the {} magic number); \
check the file or drop --format",
            what(format),
            what(format)
        ),
    }
}

/// Open `path` as `format`, checking the magic number first.
fn open_source(path: &str, format: Format) -> Result<Source, IoError> {
    let _note = panic_note(format!(
        "{path} is not a valid {} file (it may be truncated)",
        what(format)
    ));
    let mut f = File::open(path).map_err(|source| IoError::Open {
        path: path.to_string(),
        source,
    })?;
    let mut head = [0u8; 4];
    let got = {
        let mut n = 0;
        while n < 4 {
            match f.read(&mut head[n..]) {
                Ok(0) => break,
                Ok(k) => n += k,
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                Err(source) => {
                    return Err(IoError::Read {
                        path: path.to_string(),
                        source,
                    });
                }
            }
        }
        n
    };
    if bbi_format(&head[..got]) != Some(format) {
        return Err(not_this_format(path, format));
    }
    // The readers read the header from the current position.
    f.seek(SeekFrom::Start(0)).map_err(|source| IoError::Read {
        path: path.to_string(),
        source,
    })?;
    match format {
        Format::BigWig => BigWigRead::open(f)
            .map(Source::Wig)
            .map_err(|e| invalid(path, format, &e.to_string())),
        Format::BigBed => BigBedRead::open(f)
            .map(Source::Bed)
            .map_err(|e| invalid(path, format, &e.to_string())),
        _ => unreachable!("open_source is called for bbi formats"),
    }
}

fn source_chroms(source: &Source) -> Vec<(String, u32)> {
    let list = match source {
        Source::Wig(r) => r.chroms(),
        Source::Bed(r) => r.chroms(),
    };
    list.iter().map(|c| (c.name.clone(), c.length)).collect()
}

/// The chromosomes a bigWig or bigBed declares, as (file's own name, length), in its
/// chromosome list's order.
pub fn bbi_chrom_sizes(path: &str) -> Result<Vec<(Vec<u8>, u64)>, IoError> {
    let mut head = [0u8; 4];
    let got = {
        let mut f = File::open(path).map_err(|source| IoError::Open {
            path: path.to_string(),
            source,
        })?;
        f.read(&mut head).map_err(|source| IoError::Read {
            path: path.to_string(),
            source,
        })?
    };
    let format = bbi_format(&head[..got]).ok_or_else(|| IoError::Detect {
        path: path.to_string(),
        msg: format!("{path} is not a bigWig or bigBed file"),
    })?;
    let source = open_source(path, format)?;
    Ok(source_chroms(&source)
        .into_iter()
        .map(|(n, len)| (n.into_bytes(), u64::from(len)))
        .collect())
}

pub(crate) fn open_bbi(
    path: &str,
    format: Format,
    o: &ReaderOptions,
) -> Result<Box<dyn RecordReader>, IoError> {
    let source = open_source(path, format)?;
    let chroms: Vec<String> = source_chroms(&source).into_iter().map(|(n, _)| n).collect();
    if let Some(map) = &o.chrom_map {
        map.check_file_names(path, chroms.iter().map(|c| c.as_bytes()))?;
    }
    Ok(Box::new(BbiReader {
        path: path.to_string(),
        format,
        renamer: o.chrom_map.clone(),
        header: Header::None,
        diag: Diagnostics::default(),
        chroms,
        queue: None,
        source: Some(source),
        cursor: None,
        count: 0,
    }))
}

impl BbiReader {
    fn renamed(&self, file_name: &str) -> Vec<u8> {
        match &self.renamer {
            Some(m) => m.rename(file_name.as_bytes()).to_vec(),
            None => file_name.as_bytes().to_vec(),
        }
    }

    /// Order the file's chromosomes by the command's table: those it can place in its order,
    /// the rest after them in the file's own order.
    fn order(&self, chroms: &ChromTable) -> VecDeque<String> {
        let named: Vec<(Vec<u8>, &String)> =
            self.chroms.iter().map(|c| (self.renamed(c), c)).collect();
        let placed = |n: &[u8]| chroms.cmp_names(n, n).is_some();
        let (mut known, rest): (Vec<_>, Vec<_>) = named.into_iter().partition(|(n, _)| placed(n));
        known.sort_by(|(a, _), (b, _)| chroms.cmp_names(a, b).unwrap_or(Ordering::Equal));
        known
            .into_iter()
            .chain(rest)
            .map(|(_, c)| c.clone())
            .collect()
    }

    /// Start walking `file_name`.
    fn start(&mut self, file_name: &str) -> Result<(), IoError> {
        let path = self.path.clone();
        let format = self.format;
        let cursor = match self
            .source
            .take()
            .expect("source is held between chromosomes")
        {
            Source::Wig(r) => r.get_interval_move(file_name, 0, u32::MAX).map(Cursor::Wig),
            Source::Bed(r) => r.get_interval_move(file_name, 0, u32::MAX).map(Cursor::Bed),
        }
        .map_err(|e| invalid(&path, format, &e.to_string()))?;
        let out_name = self.renamed(file_name);
        self.cursor = Some((cursor, Walk { out_name, id: None }));
        Ok(())
    }
}

impl RecordReader for BbiReader {
    fn name(&self) -> &str {
        &self.path
    }
    fn format(&self) -> Format {
        self.format
    }
    fn compression(&self) -> Compression {
        Compression::None
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
        let _note = panic_note(format!(
            "{} is not a valid {} file (it may be truncated)",
            self.path,
            what(self.format)
        ));
        if self.queue.is_none() {
            self.queue = Some(self.order(chroms));
        }
        let mut n = 0;
        while n < max {
            if self.cursor.is_none() {
                let Some(next) = self.queue.as_mut().and_then(VecDeque::pop_front) else {
                    return Ok(false);
                };
                self.start(&next)?;
            }
            let (cursor, walk) = self.cursor.as_mut().expect("cursor was just set");
            let Some(item) = cursor.next() else {
                let (cursor, _) = self.cursor.take().expect("cursor is set");
                self.source = Some(cursor.finish());
                continue;
            };
            let item = item.map_err(|e| invalid(&self.path, self.format, &e.to_string()))?;
            let (start, end) = match &item {
                Item::Wig(v) => (v.start, v.end),
                Item::Bed(e) => (e.start, e.end),
            };
            if end < start {
                return Err(invalid(
                    &self.path,
                    self.format,
                    &format!("a record ends ({end}) before it starts ({start})"),
                ));
            }
            let id = match walk.id {
                Some(id) => id,
                None => {
                    let out_name = walk.out_name.clone();
                    let id = intern_chrom(chroms, &out_name).map_err(|m| IoError::Parse {
                        path: self.path.clone(),
                        loc: Loc::Record(self.count + 1),
                        msg: m,
                    })?;
                    walk.id = Some(id);
                    id
                }
            };
            let off = batch.buf.len();
            let mut strand = Strand::Unknown;
            batch.buf.extend_from_slice(&walk.out_name);
            match &item {
                Item::Wig(v) => {
                    use std::io::Write;
                    write!(batch.buf, "\t{start}\t{end}\t{}", v.value).expect("write to Vec");
                }
                Item::Bed(e) => {
                    use std::io::Write;
                    write!(batch.buf, "\t{start}\t{end}").expect("write to Vec");
                    if !e.rest.is_empty() {
                        batch.buf.push(b'\t');
                        batch.buf.extend_from_slice(e.rest.as_bytes());
                        strand = match e.rest.split('\t').nth(2) {
                            Some("+") => Strand::Plus,
                            Some("-") => Strand::Minus,
                            _ => Strand::Unknown,
                        };
                    }
                }
            }
            self.count += 1;
            batch.recs.push(RecordMeta {
                chrom: id,
                key: u64::from(start),
                span: Interval {
                    start: u64::from(start),
                    end: u64::from(end),
                },
                strand,
                loc: self.count,
                off,
                len: (batch.buf.len() - off) as u32,
            });
            n += 1;
        }
        Ok(true)
    }
}
