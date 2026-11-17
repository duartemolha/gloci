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
use gloci_core::{ChromId, ChromTable, Interval, Strand, overlaps};
use std::cmp::Ordering;
use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;
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
    file_name: String,
    out_name: Vec<u8>,
    id: Option<ChromId>,
    /// The merged regions to read, ascending; empty for the whole chromosome.
    ranges: Vec<Interval>,
    /// The range being read.
    at: usize,
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
    /// Region reads: the merged regions by the file's chromosome name.
    plan: Option<HashMap<String, Vec<Interval>>>,
    inclusive: bool,
}

/// The error for a file that is not (or no longer) a valid bbi file.
fn invalid(path: &str, format: Format, detail: &str) -> IoError {
    IoError::Detect {
        path: path.to_string(),
        msg: format!(
            "{path} is not a valid {} file (it may be truncated): {detail}",
            format.name()
        ),
    }
}

fn not_this_format(path: &str, format: Format) -> IoError {
    IoError::Detect {
        path: path.to_string(),
        msg: format!(
            "{path} is not a {} file (it does not start with the {} magic number); \
check the file or drop --format",
            format.name(),
            format.name()
        ),
    }
}

/// Up to the first 4 bytes of `f`.
fn read_magic(f: &mut File, path: &str) -> Result<Vec<u8>, IoError> {
    let mut head = [0u8; 4];
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
    Ok(head[..n].to_vec())
}

/// Open `path` as `format`, checking the magic number first.
fn open_source(path: &str, format: Format) -> Result<Source, IoError> {
    let _note = panic_note(format!(
        "{path} is not a valid {} file (it may be truncated)",
        format.name()
    ));
    let mut f = File::open(path).map_err(|source| IoError::Open {
        path: path.to_string(),
        source,
    })?;
    let head = read_magic(&mut f, path)?;
    if bbi_format(&head) != Some(format) {
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
    let mut f = File::open(path).map_err(|source| IoError::Open {
        path: path.to_string(),
        source,
    })?;
    let head = read_magic(&mut f, path)?;
    let format = bbi_format(&head).ok_or_else(|| IoError::Detect {
        path: path.to_string(),
        msg: format!("{path} is not a bigWig or bigBed file"),
    })?;
    let source = open_source(path, format)?;
    Ok(source_chroms(&source)
        .into_iter()
        .map(|(n, len)| (n.into_bytes(), u64::from(len)))
        .collect())
}

/// What the bigBed at `path` says about its columns: its autoSql text as written (`None` when
/// it has none) and its header's definedFieldCount. Read as the reader reads, so a bad file is
/// the reader's error, never a message from bigtools or a panic.
pub fn bbi_bed_columns(path: &str) -> Result<(Option<String>, u16), IoError> {
    let Source::Bed(mut r) = open_source(path, Format::BigBed)? else {
        unreachable!("open_source opens a bigBed as one");
    };
    let _note = panic_note(format!(
        "{path} is not a valid bigBed file (it may be truncated)"
    ));
    let autosql = r
        .autosql()
        .map_err(|e| invalid(path, Format::BigBed, &e.to_string()))?;
    Ok((autosql, r.info().header.defined_field_count))
}

pub(crate) fn open_bbi(
    path: &str,
    format: Format,
    o: &ReaderOptions,
) -> Result<Box<dyn RecordReader>, IoError> {
    Ok(Box::new(new_reader(path, format, o)?))
}

fn new_reader(path: &str, format: Format, o: &ReaderOptions) -> Result<BbiReader, IoError> {
    let source = open_source(path, format)?;
    let chroms: Vec<String> = source_chroms(&source).into_iter().map(|(n, _)| n).collect();
    if let Some(map) = &o.chrom_map {
        map.check_file_names(path, chroms.iter().map(|c| c.as_bytes()))?;
    }
    Ok(BbiReader {
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
        plan: None,
        inclusive: false,
    })
}

/// Sort `ivs` and merge those that overlap or touch (D10).
fn merge_ranges(mut ivs: Vec<Interval>) -> Vec<Interval> {
    ivs.sort_unstable();
    let mut out: Vec<Interval> = Vec::new();
    for iv in ivs {
        match out.last_mut() {
            Some(last) if iv.start <= last.end => last.end = last.end.max(iv.end),
            _ => out.push(iv),
        }
    }
    out
}

/// The records of the bigWig or bigBed at `path` within `regions` (named as records are after
/// --chrom-map), read through its index: bigWig rows clipped to each region, bigBed rows whole
/// and once, by spec §4.3 with `inclusive_points`; chromosomes in the order of the table
/// read_batch gets (D171).
pub fn open_bbi_regions(
    path: &Path,
    regions: &[(Vec<u8>, Interval)],
    o: &ReaderOptions,
    inclusive_points: bool,
) -> Result<Box<dyn RecordReader>, IoError> {
    let name = path.display().to_string();
    let mut f = File::open(path).map_err(|source| IoError::Open {
        path: name.clone(),
        source,
    })?;
    let head = read_magic(&mut f, &name)?;
    let format = o
        .format
        .filter(|f| f.is_bbi())
        .or_else(|| bbi_format(&head))
        .ok_or_else(|| IoError::Detect {
            path: name.clone(),
            msg: format!("{name} is not a bigWig or bigBed file"),
        })?;
    drop(f);
    let mut reader = new_reader(&name, format, o)?;
    let mut plan: HashMap<String, Vec<Interval>> = HashMap::new();
    for (chrom, iv) in regions {
        let file_name = match &o.chrom_map {
            Some(m) => m.file_name(chrom),
            None => Some(&chrom[..]),
        };
        let Some(file_name) = file_name.and_then(|n| std::str::from_utf8(n).ok()) else {
            continue;
        };
        if reader.chroms.iter().any(|c| c == file_name) {
            plan.entry(file_name.to_string()).or_default().push(*iv);
        }
    }
    reader.plan = Some(
        plan.into_iter()
            .map(|(c, ivs)| (c, merge_ranges(ivs)))
            .collect(),
    );
    reader.inclusive = inclusive_points;
    Ok(Box::new(reader))
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
            .filter(|c| self.plan.as_ref().is_none_or(|p| p.contains_key(c)))
            .collect()
    }

    /// Start walking `file_name`: the whole chromosome, or its merged regions.
    fn start(&mut self, file_name: &str) -> Result<(), IoError> {
        let ranges = self
            .plan
            .as_ref()
            .and_then(|p| p.get(file_name))
            .cloned()
            .unwrap_or_default();
        let walk = Walk {
            file_name: file_name.to_string(),
            out_name: self.renamed(file_name),
            id: None,
            ranges,
            at: 0,
        };
        self.begin(walk)
    }

    /// Open the query of `walk`'s current range. bigtools returns the records that touch it,
    /// both ends included (a point at either end too), so bigBed records are filtered by §4.3
    /// afterwards.
    fn begin(&mut self, walk: Walk) -> Result<(), IoError> {
        let (lo, hi) = match walk.ranges.get(walk.at) {
            None => (0, u32::MAX),
            Some(r) => (
                u32::try_from(r.start).unwrap_or(u32::MAX),
                u32::try_from(r.end).unwrap_or(u32::MAX),
            ),
        };
        let (path, format) = (self.path.clone(), self.format);
        let source = self
            .source
            .take()
            .expect("source is held between chromosomes");
        let cursor = match source {
            Source::Wig(r) => r
                .get_interval_move(&walk.file_name, lo, hi)
                .map(Cursor::Wig),
            Source::Bed(r) => r
                .get_interval_move(&walk.file_name, lo, hi)
                .map(Cursor::Bed),
        }
        .map_err(|e| invalid(&path, format, &e.to_string()))?;
        self.cursor = Some((cursor, walk));
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
            self.format.name()
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
                let (cursor, mut walk) = self.cursor.take().expect("cursor is set");
                self.source = Some(cursor.finish());
                walk.at += 1;
                if walk.at < walk.ranges.len() {
                    self.begin(walk)?;
                }
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
            if let Some(range) = walk.ranges.get(walk.at) {
                let span = Interval {
                    start: u64::from(start),
                    end: u64::from(end),
                };
                // A bigBed record is whole: kept when it overlaps this region and did not
                // already overlap the one before, which gave it.
                let keep = match &item {
                    Item::Wig(_) => start < end,
                    Item::Bed(_) => {
                        overlaps(span, *range, self.inclusive)
                            && !walk
                                .at
                                .checked_sub(1)
                                .is_some_and(|i| overlaps(span, walk.ranges[i], self.inclusive))
                    }
                };
                if !keep {
                    continue;
                }
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
