//! Tabix (TBI) and CSI indexes: building one from a finished BGZF file
//! (D20), finding the index next to a data file (D26), and reading the
//! records that overlap a set of regions through it (§4.3, D9, D10).
//!
//! Index spans are a superset of htslib's, so htslib queries against gloci
//! indexes stay correct:
//! - BED `[start, end)`, with a zero-length record indexed as `[s, s+1)`;
//! - GFF/GTF `[start-1, end)`, a zero-length feature as for BED;
//! - VCF/BCF `[POS-1, max(ref end, INFO/END))`. The ref end is gloci's `ref`
//!   footprint (so SVLEN-only SVs are covered), INFO/END counts whatever the
//!   ALT alleles are, and BCF also counts `rlen`.

use crate::atomic::AtomicFile;
use crate::bcf::{self, BcfDicts, Decoder, read_bcf_header};
use crate::bed::Bed;
use crate::detect::{Compression, Format, detect_compression};
use crate::diag::{Diagnostics, InvalidPolicy};
use crate::error::{IoError, Loc, UsageError};
use crate::gff::Gff;
use crate::input::{BGZF_EOF, read_up_to};
use crate::record::{
    Class, Header, ReaderOptions, RecordBatch, RecordMeta, RecordReader, TextKind, open_reader_from, warn_undeclared,
};
use crate::sorted::OrderGuard;
use crate::vcf::{Vcf, VcfHeader};
use gloci_core::{ChromId, ChromTable, Interval, RegionSet, VcfSpan, parse_coord, thousands};
use noodles_bgzf::VirtualPosition;
use noodles_bgzf::io::MultithreadedReader;
use noodles_core::Position;
use noodles_csi::BinningIndex;
use noodles_csi::binning_index::index::ReferenceSequence;
use noodles_csi::binning_index::index::header::{
    Builder as HeaderBuilder, Header as CsiHeader, ReferenceSequenceNames,
};
use noodles_csi::binning_index::index::reference_sequence::Index as BinIndexKind;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_csi::binning_index::index::reference_sequence::index::{BinnedIndex, LinearIndex};
use noodles_csi::binning_index::indexer::AddRecordError;
use noodles_csi::binning_index::{
    Index as BinningIndexOf, Indexer, ReferenceSequence as _, optimize_chunks,
};
use std::collections::BTreeMap;
use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Seek, SeekFrom, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::time::SystemTime;

/// Width of the smallest bin and of a linear-index window: 2^14 = 16 kb.
const MIN_SHIFT: u8 = 14;
/// TBI depth; 2^(14 + 3*5) = 2^29 is the largest end a TBI index can hold.
const TBI_DEPTH: u8 = 5;
/// The deepest CSI index noodles reads and writes covers 2^(14 + 3*9).
const MAX_INDEX_END: u64 = 1 << 41;
/// Records checked for order and indexed at a time ...
const BATCH: usize = 4096;
/// ... or fewer, once their raw bytes reach this, so memory stays flat with
/// multi-megabyte VCF/BCF records (100k+ samples).
const BATCH_BYTES: usize = 16 << 20;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IndexKind {
    Tbi,
    Csi,
}

impl IndexKind {
    pub fn extension(self) -> &'static str {
        match self {
            IndexKind::Tbi => "tbi",
            IndexKind::Csi => "csi",
        }
    }

    fn other(self) -> IndexKind {
        match self {
            IndexKind::Tbi => IndexKind::Csi,
            IndexKind::Csi => IndexKind::Tbi,
        }
    }
}

/// An index found next to a data file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndexInfo {
    pub kind: IndexKind,
    pub path: PathBuf,
    /// The index is older than the data file (D26).
    pub stale: bool,
}

/// `DATA.tbi` or `DATA.csi`.
pub(crate) fn index_path(data: &Path, kind: IndexKind) -> PathBuf {
    let mut s: OsString = data.as_os_str().to_owned();
    s.push(".");
    s.push(kind.extension());
    PathBuf::from(s)
}

/// Look for `DATA.tbi` and `DATA.csi`, preferring the newer one (the CSI
/// when both are equally new, since it can describe any file).
pub fn find_index(data: &Path) -> Option<IndexInfo> {
    let data_mtime = fs::metadata(data).and_then(|m| m.modified()).ok();
    let mut best: Option<(Option<SystemTime>, IndexInfo)> = None;
    for kind in [IndexKind::Tbi, IndexKind::Csi] {
        let path = index_path(data, kind);
        let Ok(meta) = fs::metadata(&path) else {
            continue;
        };
        if !meta.is_file() {
            continue;
        }
        let mtime = meta.modified().ok();
        if best.as_ref().is_some_and(|(t, _)| *t > mtime) {
            continue;
        }
        let stale = matches!((mtime, data_mtime), (Some(i), Some(d)) if i < d);
        best = Some((mtime, IndexInfo { kind, path, stale }));
    }
    best.map(|(_, info)| info)
}

/// A BGZF file with its end-of-file marker; anything else is refused.
fn check_bgzf(data: &Path, name: &str) -> Result<(), IoError> {
    let open_err = |source| IoError::Open {
        path: name.to_string(),
        source,
    };
    let read_err = |source| IoError::Read {
        path: name.to_string(),
        source,
    };
    let mut f = File::open(data).map_err(open_err)?;
    let len = f.metadata().map_err(open_err)?.len();
    let head = read_up_to(&mut f, 18).map_err(read_err)?;
    if detect_compression(&head) != Compression::Bgzf {
        return Err(UsageError(format!(
            "{name} is not BGZF-compressed; recompress with: gloci view {name} -o {name}.gz"
        ))
        .into());
    }
    let terminated = len >= 28 && {
        f.seek(SeekFrom::Start(len - 28)).map_err(read_err)?;
        read_up_to(&mut f, 28).map_err(read_err)? == BGZF_EOF
    };
    if !terminated {
        return Err(IoError::Truncated {
            path: name.to_string(),
        });
    }
    Ok(())
}

fn loc(binary: bool, n: u64) -> Loc {
    if binary { Loc::Record(n) } else { Loc::Line(n) }
}

/// A noodles position; index positions are within `1..=2^41`.
fn position(v: u64) -> Position {
    usize::try_from(v)
        .ok()
        .and_then(Position::new)
        .expect("index positions are within 1..=2^41")
}

/// Remove a trailing `\n` (and a `\r` before it) from the line at `buf[off..]`.
fn strip_eol(buf: &mut Vec<u8>, off: usize) {
    if buf.last() == Some(&b'\n') {
        buf.pop();
        if buf.len() > off && buf.last() == Some(&b'\r') {
            buf.pop();
        }
    }
}

/// INFO/END of a VCF line, if present and an integer.
fn vcf_info_end(line: &[u8]) -> Option<u64> {
    let info = line.split(|&b| b == b'\t').nth(7)?;
    info.split(|&b| b == b';')
        .find_map(|kv| kv.strip_prefix(b"END="))
        .and_then(|v| parse_coord(v).ok())
}

fn le_u32(b: &[u8], at: usize) -> u32 {
    u32::from_le_bytes([b[at], b[at + 1], b[at + 2], b[at + 3]])
}

// ---- per-format record parsing ----

/// Record parsing shared by index building and region reading.
enum Parser {
    Bed(Bed),
    Gff(Gff),
    Vcf(Vcf),
    Bcf(Box<(Decoder, BcfDicts)>),
}

impl Parser {
    fn text(format: Format, span: VcfSpan) -> Parser {
        match format {
            Format::Bed | Format::BedGraph => Parser::Bed(Bed),
            Format::Gff3 | Format::Gtf => Parser::Gff(Gff),
            Format::Vcf => Parser::Vcf(Vcf::new(span)),
            Format::Bcf => unreachable!("BCF is not a text format"),
        }
    }

    fn bcf(span: VcfSpan, dicts: BcfDicts) -> Parser {
        Parser::Bcf(Box::new((Decoder::new(span, &dicts), dicts)))
    }

    fn classify(&self, line: &[u8], seen_record: bool) -> Class {
        match self {
            Parser::Bed(_) => <Bed as TextKind>::classify(line, seen_record),
            Parser::Gff(_) => <Gff as TextKind>::classify(line, seen_record),
            Parser::Vcf(_) => <Vcf as TextKind>::classify(line, seen_record),
            Parser::Bcf(_) => Class::Record,
        }
    }

    fn parse(
        &mut self,
        raw: &[u8],
        chroms: &mut ChromTable,
    ) -> Result<(ChromId, u64, Interval), String> {
        match self {
            Parser::Bed(k) => k.parse(raw, chroms),
            Parser::Gff(k) => k.parse(raw, chroms),
            Parser::Vcf(k) => k.parse(raw, chroms),
            Parser::Bcf(b) => {
                let (dec, dicts) = &mut **b;
                dec.decode(raw, dicts, chroms)
            }
        }
    }

    /// The record's index span `(chrom, begin, end)` (D20). VCF/BCF
    /// parsers must have been made with the `ref` span.
    fn index_span(
        &mut self,
        raw: &[u8],
        chroms: &mut ChromTable,
    ) -> Result<(ChromId, u64, u64), String> {
        let (chrom, key, span) = self.parse(raw, chroms)?;
        let end = match self {
            Parser::Bed(_) | Parser::Gff(_) if span.is_point() => span.start.saturating_add(1),
            Parser::Bed(_) | Parser::Gff(_) => span.end,
            Parser::Vcf(_) => span.end.max(vcf_info_end(raw).unwrap_or(0)),
            Parser::Bcf(b) => {
                // shared section: CHROM, POS, rlen, ...
                let rlen = le_u32(raw, 16) as i32;
                let by_rlen = if rlen > 0 { key + rlen as u64 } else { 0 };
                let info_end = b.0.info_end(raw)?.unwrap_or(0);
                span.end.max(info_end).max(by_rlen)
            }
        };
        Ok((chrom, key, end))
    }

    /// Is the raw record on reference `rid` (BCF) or chromosome `name`?
    /// A record too short to tell is let through, so parsing reports it.
    fn on_reference(&self, raw: &[u8], rid: usize, name: &[u8]) -> bool {
        match self {
            Parser::Bcf(_) => {
                raw.len() < 12
                    || le_u32(raw, 0) < 24
                    || usize::try_from(le_u32(raw, 8) as i32).ok() == Some(rid)
            }
            _ => {
                let end = memchr::memchr(b'\t', raw).unwrap_or(raw.len());
                &raw[..end] == name
            }
        }
    }
}

// ---- index building ----

/// One record to index.
struct Entry {
    rid: usize,
    begin: u64,
    end: u64,
    chunk: Chunk,
    loc: u64,
}

/// What a scan learned besides the records.
struct Scanned {
    /// Reference names in id order (text formats only).
    names: Vec<Vec<u8>>,
    ref_count: usize,
    /// Leading lines htslib must skip: up to the last leading BED
    /// `track`/`browser` line (`#` lines are skipped as comments anyway).
    skip: u32,
    max_end: u64,
    /// The most raw record bytes held at once.
    #[cfg(test)]
    peak_bytes: usize,
}

/// Reads records with their virtual positions, checks order in
/// first-appearance mode and hands them to `sink` in batches.
struct Scanner<'a> {
    name: &'a str,
    binary: bool,
    chroms: ChromTable,
    guard: OrderGuard,
    batch: RecordBatch,
    entries: Vec<Entry>,
    max_end: u64,
    #[cfg(test)]
    peak_bytes: usize,
    sink: &'a mut dyn FnMut(&Entry) -> Result<(), IoError>,
}

impl Scanner<'_> {
    fn err(&self, n: u64, msg: impl Into<String>) -> IoError {
        IoError::Parse {
            path: self.name.to_string(),
            loc: loc(self.binary, n),
            msg: msg.into(),
        }
    }

    fn read_err(&self, source: io::Error) -> IoError {
        IoError::Read {
            path: self.name.to_string(),
            source,
        }
    }

    /// Queue the record whose raw bytes are `batch.buf[off..]`.
    #[allow(clippy::too_many_arguments)]
    fn push(
        &mut self,
        chrom: ChromId,
        rid: usize,
        begin: u64,
        end: u64,
        n: u64,
        off: usize,
        chunk: Chunk,
    ) -> Result<(), IoError> {
        if end > MAX_INDEX_END {
            return Err(self.err(
                n,
                format!(
                    "the record ends at {}, beyond {}, the largest position an index can hold",
                    thousands(end),
                    thousands(MAX_INDEX_END)
                ),
            ));
        }
        let len = u32::try_from(self.batch.buf.len() - off)
            .map_err(|_| self.err(n, "record is longer than 4 GiB"))?;
        self.batch.recs.push(RecordMeta {
            chrom,
            key: begin,
            span: Interval { start: begin, end },
            loc: n,
            off,
            len,
        });
        self.entries.push(Entry {
            rid,
            begin,
            end,
            chunk,
            loc: n,
        });
        self.max_end = self.max_end.max(end);
        #[cfg(test)]
        {
            self.peak_bytes = self.peak_bytes.max(self.batch.buf.len());
        }
        if self.entries.len() >= BATCH || self.batch.buf.len() >= BATCH_BYTES {
            self.flush()?;
        }
        Ok(())
    }

    fn flush(&mut self) -> Result<(), IoError> {
        self.guard.check(&self.chroms, &self.batch)?;
        for e in &self.entries {
            (self.sink)(e)?;
        }
        self.batch.clear();
        self.entries.clear();
        Ok(())
    }

    /// Text formats; returns the line skip count.
    fn text(&mut self, r: &mut MultithreadedReader<File>, format: Format) -> Result<u32, IoError> {
        let mut parser = Parser::text(format, VcfSpan::Ref);
        let bed = matches!(format, Format::Bed | Format::BedGraph);
        let (mut n, mut skip, mut seen) = (0u64, 0u32, false);
        loop {
            let start = r.virtual_position();
            let off = self.batch.buf.len();
            let got = r
                .read_until(b'\n', &mut self.batch.buf)
                .map_err(|e| self.read_err(e))?;
            if got == 0 {
                break;
            }
            n += 1;
            let end = r.virtual_position();
            strip_eol(&mut self.batch.buf, off);
            let line = &self.batch.buf[off..];
            match parser.classify(line, seen) {
                Class::Record => {}
                Class::Header if !seen && bed && line.first() != Some(&b'#') => {
                    skip = u32::try_from(n).unwrap_or(u32::MAX);
                    self.batch.buf.truncate(off);
                    continue;
                }
                Class::Header | Class::Skip => {
                    self.batch.buf.truncate(off);
                    continue;
                }
                Class::Fatal(msg) => return Err(self.err(n, msg)),
            }
            seen = true;
            let (chrom, begin, end_pos) = parser
                .index_span(line, &mut self.chroms)
                .map_err(|msg| self.err(n, msg))?;
            self.push(
                chrom,
                chrom.0 as usize,
                begin,
                end_pos,
                n,
                off,
                Chunk::new(start, end),
            )?;
        }
        Ok(skip)
    }

    /// BCF; returns the reference count (the contig dictionary size).
    fn bcf(&mut self, r: &mut MultithreadedReader<File>) -> Result<usize, IoError> {
        let h = read_bcf_header(r).map_err(|msg| IoError::Detect {
            msg: format!("{}: {msg}", self.name),
            path: self.name.to_string(),
        })?;
        let dicts = h.bcf.expect("read_bcf_header builds the dictionaries");
        // Chromosome ids follow the dictionary, so the order check enforces
        // the dictionary order BCF indexes need.
        for c in dicts.contigs.iter().filter(|c| !c.is_empty()) {
            self.chroms
                .intern(c)
                .expect("a first-appearance table is open");
        }
        let ref_count = dicts.contigs.len();
        let mut parser = Parser::bcf(VcfSpan::Ref, dicts);
        let mut n = 0u64;
        loop {
            let start = r.virtual_position();
            let off = self.batch.buf.len();
            let got = bcf::append(r, &mut self.batch.buf, 8).map_err(|e| self.read_err(e))?;
            if got == 0 {
                break;
            }
            n += 1;
            if got < 8 {
                return Err(self.err(n, bcf::RECORD_TRUNCATED));
            }
            let b = &self.batch.buf[off..];
            let body = u64::from(le_u32(b, 0)) + u64::from(le_u32(b, 4));
            if body + 8 > u64::from(u32::MAX) {
                return Err(self.err(n, "record is larger than 4 GiB"));
            }
            let body = body as usize;
            let got = bcf::append(r, &mut self.batch.buf, body).map_err(|e| self.read_err(e))?;
            if got < body {
                return Err(self.err(n, bcf::RECORD_TRUNCATED));
            }
            let end = r.virtual_position();
            let raw = &self.batch.buf[off..];
            let (chrom, begin, end_pos) = parser
                .index_span(raw, &mut self.chroms)
                .map_err(|msg| self.err(n, msg))?;
            // decode() checked that CHROM is a dictionary index
            let rid = le_u32(raw, 8) as usize;
            self.push(chrom, rid, begin, end_pos, n, off, Chunk::new(start, end))?;
        }
        Ok(ref_count)
    }
}

/// Read every record of `data` (BGZF) with its chunk and pass it to `sink`.
fn scan(
    data: &Path,
    name: &str,
    format: Format,
    threads: NonZeroUsize,
    sink: &mut dyn FnMut(&Entry) -> Result<(), IoError>,
) -> Result<Scanned, IoError> {
    let file = File::open(data).map_err(|source| IoError::Open {
        path: name.to_string(),
        source,
    })?;
    let mut r = MultithreadedReader::with_worker_count(threads, file);
    let binary = format == Format::Bcf;
    let mut sc = Scanner {
        name,
        binary,
        chroms: ChromTable::first_appearance(),
        guard: OrderGuard::new(name, format),
        batch: RecordBatch::default(),
        entries: Vec::with_capacity(BATCH),
        max_end: 0,
        #[cfg(test)]
        peak_bytes: 0,
        sink,
    };
    let (skip, ref_count) = if binary {
        (0, sc.bcf(&mut r)?)
    } else {
        let skip = sc.text(&mut r, format)?;
        (skip, sc.chroms.len())
    };
    sc.flush()?;
    let names = if binary {
        Vec::new()
    } else {
        (0..sc.chroms.len() as u32)
            .map(|i| sc.chroms.name(ChromId(i)).to_vec())
            .collect()
    };
    Ok(Scanned {
        names,
        ref_count,
        skip,
        max_end: sc.max_end,
        #[cfg(test)]
        peak_bytes: sc.peak_bytes,
    })
}

/// A sparse htslib-style linear index of one reference: runs of 16 kb
/// windows, each starting at `(first window, offset of the first record
/// that reaches it)`.
#[derive(Debug, Default)]
struct Linear {
    runs: Vec<(u64, VirtualPosition)>,
    covered: u64,
}

impl Linear {
    fn add(&mut self, end: u64, start: VirtualPosition) {
        let last = (end - 1) >> MIN_SHIFT;
        if last >= self.covered {
            self.runs.push((self.covered, start));
            self.covered = last + 1;
        }
    }

    fn at(&self, window: u64) -> VirtualPosition {
        match self.runs.partition_point(|&(w, _)| w <= window) {
            0 => VirtualPosition::default(),
            i => self.runs[i - 1].1,
        }
    }
}

/// The result of one indexing pass; `index` is `None` when a record ended
/// beyond the range of the chosen depth.
struct Pass<I> {
    index: Option<BinningIndexOf<I>>,
    scanned: Scanned,
    linear: Vec<Linear>,
}

fn text_header(format: Format, s: &Scanned) -> Option<CsiHeader> {
    let builder = match format {
        Format::Bed | Format::BedGraph => HeaderBuilder::bed().set_line_skip_count(s.skip),
        Format::Gff3 | Format::Gtf => HeaderBuilder::gff(),
        Format::Vcf => HeaderBuilder::vcf(),
        Format::Bcf => return None,
    };
    let mut names = ReferenceSequenceNames::default();
    for n in &s.names {
        names.insert(n.clone().into());
    }
    Some(builder.set_reference_sequence_names(names).build())
}

fn pass<I: BinIndexKind + Default>(
    data: &Path,
    name: &str,
    format: Format,
    threads: NonZeroUsize,
    depth: u8,
) -> Result<Pass<I>, IoError> {
    let mut indexer = Indexer::<I>::new(MIN_SHIFT, depth).expect("depth is within 5..=9");
    let mut linear: Vec<Linear> = Vec::new();
    let mut full = false;
    let binary = format == Format::Bcf;
    let mut sink = |e: &Entry| -> Result<(), IoError> {
        if full {
            return Ok(());
        }
        let ctx = Some((e.rid, position(e.begin + 1), position(e.end), true));
        match indexer.add_record(ctx, e.chunk) {
            Ok(()) => {}
            Err(AddRecordError::EndExceedsMaxPosition { .. }) => {
                full = true;
                return Ok(());
            }
            Err(err) => {
                return Err(IoError::Parse {
                    path: name.to_string(),
                    loc: loc(binary, e.loc),
                    msg: format!("cannot index the record: {err}"),
                });
            }
        }
        if linear.len() <= e.rid {
            linear.resize_with(e.rid + 1, Linear::default);
        }
        linear[e.rid].add(e.end, e.chunk.start());
        Ok(())
    };
    let scanned = scan(data, name, format, threads, &mut sink)?;
    let index = if full {
        None
    } else {
        let indexer = match text_header(format, &scanned) {
            Some(h) => indexer.set_header(h),
            None => indexer,
        };
        Some(indexer.build(scanned.ref_count))
    };
    Ok(Pass {
        index,
        scanned,
        linear,
    })
}

/// The smallest depth (at least 5) whose range covers `max_end`.
fn csi_depth(max_end: u64) -> u8 {
    let mut d = TBI_DEPTH;
    while (1u64 << (u32::from(MIN_SHIFT) + 3 * u32::from(d))) < max_end {
        d += 1;
    }
    d
}

/// First position covered by bin `id`.
fn bin_start(id: usize, min_shift: u8, depth: u8) -> u64 {
    let mut first = 0usize;
    for l in 0..=u32::from(depth) {
        let next = first + (1usize << (3 * l));
        if id < next {
            let shift = u32::from(min_shift) + 3 * (u32::from(depth) - l);
            return ((id - first) as u64) << shift;
        }
        first = next;
    }
    0
}

/// Give every bin the loffset htslib would: the offset of the first record
/// that reaches the bin's first 16 kb window. noodles records the first
/// record *in* the bin's subtree instead, which lets htslib prune a long
/// record that sits in an ancestor bin.
fn with_htslib_loffsets(index: noodles_csi::Index, linear: &[Linear]) -> noodles_csi::Index {
    let (min_shift, depth) = (index.min_shift(), index.depth());
    let refs = index
        .reference_sequences()
        .iter()
        .enumerate()
        .map(|(rid, rs)| {
            let mut loffsets = BinnedIndex::default();
            for &id in rs.bins().keys() {
                let window = bin_start(id, min_shift, depth) >> MIN_SHIFT;
                let off = linear
                    .get(rid)
                    .map_or_else(VirtualPosition::default, |l| l.at(window));
                loffsets.insert(id, off);
            }
            ReferenceSequence::new(rs.bins().clone(), loffsets, rs.metadata().cloned())
        })
        .collect();
    let mut b = noodles_csi::Index::builder()
        .set_min_shift(min_shift)
        .set_depth(depth)
        .set_reference_sequences(refs);
    if let Some(h) = index.header() {
        b = b.set_header(h.clone());
    }
    if let Some(n) = index.unplaced_unmapped_record_count() {
        b = b.set_unplaced_unmapped_record_count(n);
    }
    b.build()
}

fn tbi_bytes(index: &noodles_tabix::Index) -> io::Result<Vec<u8>> {
    let mut w = noodles_tabix::io::Writer::new(Vec::new());
    w.write_index(index)?;
    w.into_inner().finish()
}

fn csi_bytes(index: &noodles_csi::Index) -> io::Result<Vec<u8>> {
    let mut w = noodles_csi::io::Writer::new(Vec::new());
    w.write_index(index)?;
    w.into_inner().finish()
}

/// Build the TBI or CSI index (D20) of the finished BGZF file `data` and
/// return its kind and bytes, without writing anything. `data` may be a
/// temp file that becomes `final_data_path`; messages name
/// `final_data_path`.
///
/// Order is checked in first-appearance mode (BCF: dictionary order). TBI is
/// used unless an index span ends beyond 2^29; then CSI with min_shift 14
/// and the smallest depth (at least 5) that covers the largest end. BCF
/// always gets CSI, keyed by the contig dictionary.
pub(crate) fn build_index(
    data: &Path,
    final_data_path: &Path,
    format: Format,
    threads: NonZeroUsize,
) -> Result<(IndexKind, Vec<u8>), IoError> {
    let name = final_data_path.display().to_string();
    check_bgzf(data, &name)?;
    let changed = || IoError::Detect {
        path: name.clone(),
        msg: format!("{name} changed while it was being indexed; run the command again"),
    };
    let csi = |p: Pass<BinnedIndex>| match p.index {
        Some(index) => Ok(with_htslib_loffsets(index, &p.linear)),
        None => Err(changed()),
    };
    let (kind, bytes) = if format == Format::Bcf {
        let p = pass::<BinnedIndex>(data, &name, format, threads, TBI_DEPTH)?;
        let p = match p.index {
            Some(_) => p,
            None => {
                let depth = csi_depth(p.scanned.max_end);
                pass::<BinnedIndex>(data, &name, format, threads, depth)?
            }
        };
        (IndexKind::Csi, csi_bytes(&csi(p)?))
    } else {
        let p = pass::<LinearIndex>(data, &name, format, threads, TBI_DEPTH)?;
        match p.index {
            Some(index) => (IndexKind::Tbi, tbi_bytes(&index)),
            None => {
                let depth = csi_depth(p.scanned.max_end);
                let q = pass::<BinnedIndex>(data, &name, format, threads, depth)?;
                (IndexKind::Csi, csi_bytes(&csi(q)?))
            }
        }
    };
    let bytes = bytes.map_err(|source| IoError::Write {
        path: index_path(final_data_path, kind).display().to_string(),
        source,
    })?;
    Ok((kind, bytes))
}

/// Write index `bytes` to a temp file next to `path`; publish it with
/// [`AtomicFile::commit`].
pub(crate) fn stage_index(path: &Path, bytes: &[u8]) -> Result<AtomicFile, IoError> {
    let shown = || path.display().to_string();
    let mut af = AtomicFile::create(path).map_err(|source| IoError::Open {
        path: shown(),
        source,
    })?;
    af.write_all(bytes).map_err(|source| IoError::Write {
        path: shown(),
        source,
    })?;
    Ok(af)
}

/// Remove `DATA.tbi` and `DATA.csi` if they exist (D19).
pub(crate) fn remove_indexes(data: &Path) -> Result<(), IoError> {
    for kind in [IndexKind::Tbi, IndexKind::Csi] {
        let p = index_path(data, kind);
        match fs::remove_file(&p) {
            Ok(()) => {}
            Err(e) if e.kind() == io::ErrorKind::NotFound => {}
            Err(source) => {
                return Err(IoError::Write {
                    path: p.display().to_string(),
                    source,
                });
            }
        }
    }
    Ok(())
}

/// Build a TBI or CSI index for the BGZF file `data` (see `build_index`),
/// write it atomically to `final_data_path` + `.tbi`/`.csi`, remove an index
/// of the other kind, and return the index path. Messages name
/// `final_data_path`.
pub fn write_index(
    data: &Path,
    final_data_path: &Path,
    format: Format,
    threads: NonZeroUsize,
) -> Result<PathBuf, IoError> {
    let (kind, bytes) = build_index(data, final_data_path, format, threads)?;
    let path = index_path(final_data_path, kind);
    let af = stage_index(&path, &bytes)?;
    af.commit().map_err(|source| IoError::Write {
        path: path.display().to_string(),
        source,
    })?;
    let other = index_path(final_data_path, kind.other());
    match fs::remove_file(&other) {
        Ok(()) => {}
        Err(e) if e.kind() == io::ErrorKind::NotFound => {}
        Err(source) => {
            return Err(IoError::Write {
                path: other.display().to_string(),
                source,
            });
        }
    }
    Ok(path)
}

// ---- index-backed region reader ----

enum Loaded {
    Tbi(noodles_tabix::Index),
    Csi(noodles_csi::Index),
}

impl Loaded {
    fn read(info: &IndexInfo) -> Result<Loaded, IoError> {
        let err = |source| IoError::Read {
            path: info.path.display().to_string(),
            source,
        };
        Ok(match info.kind {
            IndexKind::Tbi => Loaded::Tbi(noodles_tabix::fs::read(&info.path).map_err(err)?),
            IndexKind::Csi => Loaded::Csi(noodles_csi::fs::read(&info.path).map_err(err)?),
        })
    }

    fn names(&self) -> Option<&ReferenceSequenceNames> {
        let h = match self {
            Loaded::Tbi(i) => i.header(),
            Loaded::Csi(i) => i.header(),
        };
        h.map(CsiHeader::reference_sequence_names)
    }

    /// Where the records of reference `rid` begin in the file: the smallest
    /// chunk start among its bins and its metadata pseudo-bin. `None` when
    /// the index holds no records for it.
    fn first_offset(&self, rid: usize) -> Option<VirtualPosition> {
        fn first<I: BinIndexKind>(rs: &ReferenceSequence<I>) -> Option<VirtualPosition> {
            rs.bins()
                .values()
                .flat_map(|b| b.chunks())
                .map(|c| c.start())
                .chain(rs.metadata().map(|m| m.start_position()))
                .min()
        }
        match self {
            Loaded::Tbi(i) => i.reference_sequences().get(rid).and_then(first),
            Loaded::Csi(i) => i.reference_sequences().get(rid).and_then(first),
        }
    }

    /// One past the largest position the index can hold.
    fn max_position(&self) -> u64 {
        let (min_shift, depth) = match self {
            Loaded::Tbi(i) => (i.min_shift(), i.depth()),
            Loaded::Csi(i) => (i.min_shift(), i.depth()),
        };
        1u64.checked_shl(u32::from(min_shift) + 3 * u32::from(depth))
            .unwrap_or(u64::MAX)
    }

    /// Chunks, in file order, that hold every record whose index span
    /// overlaps the 0-based window `[lo, hi)`.
    fn chunks(&self, rid: usize, lo: u64, hi: u64) -> Vec<Chunk> {
        match self {
            Loaded::Tbi(i) => {
                let Some(rs) = i.reference_sequences().get(rid) else {
                    return Vec::new();
                };
                // TBI linear index: the first record reaching lo's window
                let min_off = rs.min_offset(i.min_shift(), i.depth(), position(lo + 1));
                bin_chunks(rs, i.min_shift(), i.depth(), lo, hi, min_off)
            }
            Loaded::Csi(i) => {
                let Some(rs) = i.reference_sequences().get(rid) else {
                    return Vec::new();
                };
                let min_off = csi_min_offset(rs.index(), i.min_shift(), i.depth(), lo);
                bin_chunks(rs, i.min_shift(), i.depth(), lo, hi, min_off)
            }
        }
    }
}

/// htslib's CSI lower bound for a query starting at `lo`: the loffset of
/// lo's leaf bin or, if it has none, of the nearest bin to its left or above.
fn csi_min_offset(loffsets: &BinnedIndex, min_shift: u8, depth: u8, lo: u64) -> VirtualPosition {
    let first_leaf = ((1usize << (3 * u32::from(depth))) - 1) / 7;
    let mut bin = first_leaf + (lo >> min_shift) as usize;
    loop {
        if let Some(&off) = loffsets.get(&bin) {
            return off;
        }
        if bin == 0 {
            return VirtualPosition::default();
        }
        let parent = (bin - 1) >> 3;
        bin = if bin > (parent << 3) + 1 {
            bin - 1
        } else {
            parent
        };
    }
}

/// The chunks of every bin that overlaps `[lo, hi)`, merged, without the
/// parts before `min_off` (which hold only records that end before `lo`).
fn bin_chunks<I: BinIndexKind>(
    rs: &ReferenceSequence<I>,
    min_shift: u8,
    depth: u8,
    lo: u64,
    hi: u64,
    min_off: VirtualPosition,
) -> Vec<Chunk> {
    let bins = rs.bins();
    let mut chunks = Vec::new();
    let mut first = 0u64;
    for l in 0..=u32::from(depth) {
        let shift = u32::from(min_shift) + 3 * (u32::from(depth) - l);
        let (b, e) = (first + (lo >> shift), first + ((hi - 1) >> shift));
        if e - b < bins.len() as u64 {
            for id in b..=e {
                if let Some(bin) = bins.get(&(id as usize)) {
                    chunks.extend_from_slice(bin.chunks());
                }
            }
        } else {
            for (&id, bin) in bins {
                if (b..=e).contains(&(id as u64)) {
                    chunks.extend_from_slice(bin.chunks());
                }
            }
        }
        first += 1 << (3 * l);
    }
    optimize_chunks(&chunks, min_off)
        .into_iter()
        .filter_map(|c| {
            let start = c.start().max(min_off);
            (start < c.end()).then(|| Chunk::new(start, c.end()))
        })
        .collect()
}

/// The requested regions on one reference.
struct Plan {
    rid: usize,
    name: Vec<u8>,
    /// Records are tested against each region, not against their union: a
    /// point between two bookended regions is in neither (§4.3).
    set: RegionSet,
    /// Disjoint query windows in coordinate order: each region widened by
    /// one base on both sides (so points on region edges, and zero-length
    /// records indexed as `[s, s+1)`, are found), merged, and clamped to
    /// the index's range. 0-based, half-open.
    windows: Vec<(u64, u64)>,
}

impl Plan {
    fn new(rid: usize, name: Vec<u8>, regions: Vec<Interval>, max_pos: u64) -> Plan {
        let set = RegionSet::new(regions);
        let mut windows: Vec<(u64, u64)> = Vec::new();
        for r in set.regions() {
            let lo = r.start.saturating_sub(1);
            let hi = r.end.saturating_add(1).min(max_pos);
            if lo >= hi {
                continue; // beyond anything the index can hold
            }
            match windows.last_mut() {
                Some(w) if lo <= w.1 => w.1 = w.1.max(hi),
                _ => windows.push((lo, hi)),
            }
        }
        Plan {
            rid,
            name,
            set,
            windows,
        }
    }
}

/// Yields the records whose span overlaps a region, each once, in file
/// order. Windows are visited in coordinate order; each window's chunks are
/// read in file order until a record starts at or after the window's end,
/// and every record read is tested against all regions of its reference.
/// A record is emitted only if it lies after the last emitted one, which
/// drops later sightings of a record already emitted.
///
/// Why nothing is lost: a record Q matching region r is read while r's
/// window W is processed (Q's index span overlaps W, so Q is in W's chunks,
/// and Q starts before W's end). If an earlier window V emitted some E
/// after Q, then Q starts before V's end and reaches W, so Q's index span
/// overlaps V too and Q was read, and emitted, before E while V was
/// processed. Chunks are therefore always read from their start: skipping
/// what an earlier window read could skip a record that the earlier window
/// jumped over between two of its chunks.
struct RegionReader {
    name: String,
    format: Format,
    header: Header,
    diag: Diagnostics,
    policy: InvalidPolicy,
    inclusive: bool,
    src: noodles_bgzf::io::Reader<BufReader<File>>,
    parser: Parser,
    index: Loaded,
    plans: Vec<Plan>,
    /// Current plan and window.
    p: usize,
    w: usize,
    /// The current window's chunks (valid when `loaded`) and the next one.
    chunks: Vec<Chunk>,
    c: usize,
    loaded: bool,
    /// The end of the chunk being read.
    chunk_end: Option<VirtualPosition>,
    last_emitted: Option<VirtualPosition>,
}

impl RegionReader {
    fn read_err(&self, source: io::Error) -> IoError {
        IoError::Read {
            path: self.name.clone(),
            source,
        }
    }

    /// Append the next record (a line without its terminator, or a framed
    /// BCF record) to `buf`; false at end of input.
    fn read_raw(&mut self, buf: &mut Vec<u8>) -> Result<bool, IoError> {
        let off = buf.len();
        if let Parser::Bcf(_) = self.parser {
            let got = bcf::append(&mut self.src, buf, 8).map_err(|e| self.read_err(e))?;
            if got == 0 {
                return Ok(false);
            }
            if got < 8 {
                return self.bad_frame(buf, off);
            }
            let body = u64::from(le_u32(buf, off)) + u64::from(le_u32(buf, off + 4));
            if body + 8 > u64::from(u32::MAX) {
                return self.bad_frame(buf, off);
            }
            let body = body as usize;
            let got = bcf::append(&mut self.src, buf, body).map_err(|e| self.read_err(e))?;
            if got < body {
                return self.bad_frame(buf, off);
            }
        } else {
            let got = self
                .src
                .read_until(b'\n', buf)
                .map_err(|e| self.read_err(e))?;
            if got == 0 {
                return Ok(false);
            }
            strip_eol(buf, off);
        }
        Ok(true)
    }

    fn bad_frame(&self, buf: &mut Vec<u8>, off: usize) -> Result<bool, IoError> {
        buf.truncate(off);
        Err(IoError::Detect {
            path: self.name.clone(),
            msg: format!(
                "{}: a BCF record in an indexed region is truncated or corrupt; check the file with: gloci check {}",
                self.name, self.name
            ),
        })
    }

    /// An invalid record at `at`: an error, or counted under `--skip-invalid`.
    /// Index-backed reads have no line numbers, so the place is the record's
    /// virtual offset.
    fn bad(&mut self, at: VirtualPosition, msg: &str) -> Result<(), IoError> {
        let (c, u): (u64, u16) = at.into();
        let place = format!("record at BGZF offset {c}:{u}");
        match self.policy {
            InvalidPolicy::Fail => Err(IoError::Detect {
                path: self.name.clone(),
                msg: format!("{} {place}: {msg}", self.name),
            }),
            InvalidPolicy::Skip => {
                self.diag.record_place(&place, msg);
                Ok(())
            }
        }
    }

    /// Append the next matching record to `batch`; false when done.
    fn next(&mut self, chroms: &mut ChromTable, batch: &mut RecordBatch) -> Result<bool, IoError> {
        loop {
            let Some(plan) = self.plans.get(self.p) else {
                return Ok(false);
            };
            let rid = plan.rid;
            let Some(&(lo, hi)) = plan.windows.get(self.w) else {
                self.p += 1;
                self.w = 0;
                self.loaded = false;
                continue;
            };
            if !self.loaded {
                self.chunks = self.index.chunks(rid, lo, hi);
                self.c = 0;
                self.chunk_end = None;
                self.loaded = true;
            }
            let chunk_end = match self.chunk_end {
                Some(e) => e,
                None => {
                    let Some(&ch) = self.chunks.get(self.c) else {
                        self.w += 1;
                        self.loaded = false;
                        continue;
                    };
                    if self.src.virtual_position() != ch.start() {
                        self.src.seek(ch.start()).map_err(|e| self.read_err(e))?;
                    }
                    self.chunk_end = Some(ch.end());
                    ch.end()
                }
            };
            let start = self.src.virtual_position();
            if start >= chunk_end {
                self.chunk_end = None;
                self.c += 1;
                continue;
            }
            let off = batch.buf.len();
            if !self.read_raw(&mut batch.buf)? {
                self.chunk_end = None;
                self.c += 1;
                continue;
            }
            let raw = &batch.buf[off..];
            let plan = &self.plans[self.p];
            let class = self.parser.classify(raw, true);
            if matches!(class, Class::Header | Class::Skip)
                || !self.parser.on_reference(raw, rid, &plan.name)
            {
                batch.buf.truncate(off);
                continue;
            }
            let parsed = match class {
                Class::Fatal(msg) => Err(msg),
                _ if raw.len() > u32::MAX as usize => Err("line is longer than 4 GiB".into()),
                _ => self.parser.parse(raw, chroms),
            };
            let (chrom, key, span) = match parsed {
                Ok(p) => p,
                Err(msg) => {
                    batch.buf.truncate(off);
                    self.bad(start, &msg)?;
                    continue;
                }
            };
            if chroms.has_undeclared() {
                warn_undeclared(chroms, &self.name, None);
            }
            if key >= hi {
                // past this window; later windows may still want it
                batch.buf.truncate(off);
                self.chunk_end = None;
                self.w += 1;
                self.loaded = false;
                continue;
            }
            if plan.set.matches(span, self.inclusive) && self.last_emitted.is_none_or(|l| start > l)
            {
                let len = (batch.buf.len() - off) as u32;
                batch.recs.push(RecordMeta {
                    chrom,
                    key,
                    span,
                    loc: 0,
                    off,
                    len,
                });
                self.last_emitted = Some(start);
                return Ok(true);
            }
            batch.buf.truncate(off);
        }
    }
}

impl RecordReader for RegionReader {
    fn name(&self) -> &str {
        &self.name
    }
    fn format(&self) -> Format {
        self.format
    }
    fn compression(&self) -> Compression {
        Compression::Bgzf
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
        for _ in 0..max {
            if !self.next(chroms, batch)? {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

/// Read the records of the indexed BGZF file `path` whose span (per
/// `o.vcf_span`) overlaps one of `regions` under §4.3 (D9). Each record is
/// yielded once, in file order: the regions are merged per reference (D10)
/// and the references are visited in the order their records appear in the
/// file, which need not be their id order (htslib indexes a BCF whose
/// contigs are contiguous but not in dictionary order). Regions on
/// chromosomes the index does not know yield nothing. The header is read
/// from the start of the file.
///
/// Records carry `loc` 0: their line numbers are unknown. Parse errors name
/// the record's BGZF virtual offset instead.
pub fn open_region_reader(
    path: &Path,
    regions: &[(Vec<u8>, Interval)],
    o: &ReaderOptions,
    inclusive_points: bool,
) -> Result<Box<dyn RecordReader>, IoError> {
    let name = path.display().to_string();
    check_bgzf(path, &name)?;
    let info = find_index(path).ok_or_else(|| {
        UsageError(format!(
            "{name} has no index; create one with: gloci index {name}"
        ))
    })?;
    let open_err = |source| IoError::Open {
        path: name.clone(),
        source,
    };
    let head_opts = ReaderOptions {
        threads: NonZeroUsize::MIN,
        ..o.clone()
    };
    let head = open_reader_from(
        &name,
        Box::new(File::open(path).map_err(open_err)?),
        &head_opts,
    )?;
    let (format, header) = (head.format(), head.header().clone());
    drop(head);
    let index = Loaded::read(&info)?;

    // region chromosome -> index reference id
    let parser = match (&header, format) {
        (
            Header::Vcf(VcfHeader {
                bcf: Some(dicts), ..
            }),
            Format::Bcf,
        ) => Parser::bcf(o.vcf_span, dicts.clone()),
        _ => Parser::text(format, o.vcf_span),
    };
    let rid_of = |chrom: &[u8]| -> Option<usize> {
        match &parser {
            Parser::Bcf(b) => b.1.contigs.iter().position(|c| !c.is_empty() && c == chrom),
            _ => index.names().and_then(|n| n.get_index_of(chrom)),
        }
    };
    let mut by_rid: BTreeMap<usize, (Vec<u8>, Vec<Interval>)> = BTreeMap::new();
    for (chrom, iv) in regions {
        if let Some(rid) = rid_of(chrom) {
            by_rid
                .entry(rid)
                .or_insert_with(|| (chrom.clone(), Vec::new()))
                .1
                .push(*iv);
        }
    }
    let max_pos = index.max_position();
    let mut plans: Vec<Plan> = by_rid
        .into_iter()
        .map(|(rid, (name, ivs))| Plan::new(rid, name, ivs, max_pos))
        .collect();
    // File order, so that `last_emitted` only ever drops a record already
    // emitted. References without records yield nothing; their place is
    // irrelevant.
    plans.sort_by_key(|p| index.first_offset(p.rid).unwrap_or(VirtualPosition::MAX));

    let file = File::open(path).map_err(open_err)?;
    Ok(Box::new(RegionReader {
        name,
        format,
        header,
        diag: Diagnostics::default(),
        policy: o.invalid,
        inclusive: inclusive_points,
        src: noodles_bgzf::io::Reader::new(BufReader::new(file)),
        parser,
        index,
        plans,
        p: 0,
        w: 0,
        chunks: Vec::new(),
        c: 0,
        loaded: false,
        chunk_end: None,
        last_emitted: None,
    }))
}

/// The reference names the index `info` of `data` knows, in reference id
/// order: the names stored in a TBI/CSI header, or, for a BCF's CSI (which
/// has none), the contig dictionary of the BCF header.
pub fn index_reference_names(data: &Path, info: &IndexInfo) -> Result<Vec<Vec<u8>>, IoError> {
    if let Some(names) = Loaded::read(info)?.names() {
        return Ok(names.iter().map(|n| n.to_vec()).collect());
    }
    let name = data.display().to_string();
    let file = File::open(data).map_err(|source| IoError::Open {
        path: name.clone(),
        source,
    })?;
    let mut src = noodles_bgzf::io::Reader::new(BufReader::new(file));
    let h = read_bcf_header(&mut src).map_err(|msg| IoError::Detect {
        msg: format!(
            "{name}: the index {} lists no chromosome names, as only a BCF index may, \
but the BCF header cannot be read: {msg}",
            info.path.display()
        ),
        path: name.clone(),
    })?;
    let dicts = h.bcf.expect("read_bcf_header builds the dictionaries");
    Ok(dicts
        .contigs
        .into_iter()
        .filter(|c| !c.is_empty())
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Records over 1 MB each (a 300k-sample VCF): the scan holds at most
    /// `BATCH_BYTES` plus one record, not a whole 4096-record batch, and
    /// queries still find the right records.
    #[test]
    fn scan_memory_is_bounded_by_bytes() {
        let samples = "\t0|1".repeat(300_000);
        let mut vcf = b"##fileformat=VCFv4.3\n##contig=<ID=chr1>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n"
            .to_vec();
        let mut longest = 0;
        let mut records = 0;
        for i in 0..24u64 {
            let line = format!(
                "chr1\t{}\t.\tA\tG\t.\tPASS\t.\tGT{samples}\n",
                1000 + i * 100
            );
            longest = longest.max(line.len());
            records += line.len();
            vcf.extend_from_slice(line.as_bytes());
        }
        assert!(records > BATCH_BYTES + 8 * longest);
        let mut w = noodles_bgzf::io::Writer::new(Vec::new());
        w.write_all(&vcf).unwrap();
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("big.vcf.gz");
        fs::write(&p, w.finish().unwrap()).unwrap();

        let mut n = 0;
        let mut sink = |_: &Entry| -> Result<(), IoError> {
            n += 1;
            Ok(())
        };
        let threads = NonZeroUsize::new(2).unwrap();
        let s = scan(&p, "big.vcf.gz", Format::Vcf, threads, &mut sink).unwrap();
        assert_eq!(n, 24);
        assert!(
            s.peak_bytes <= BATCH_BYTES + longest,
            "held {} bytes, budget {} + one record of {longest}",
            s.peak_bytes,
            BATCH_BYTES
        );

        write_index(&p, &p, Format::Vcf, threads).unwrap();
        let o = ReaderOptions {
            threads,
            format: None,
            vcf_span: VcfSpan::Allele,
            invalid: InvalidPolicy::Fail,
        };
        let mut chroms = ChromTable::first_appearance();
        let mut batch = RecordBatch::default();
        let mut query = |iv: Interval| -> Vec<(u64, u32)> {
            let regions = [(b"chr1".to_vec(), iv)];
            let mut r = open_region_reader(&p, &regions, &o, false).unwrap();
            let mut out = Vec::new();
            loop {
                batch.clear();
                let more = r.read_batch(&mut chroms, &mut batch, 4).unwrap();
                out.extend(batch.recs.iter().map(|m| (m.key, m.len)));
                if !more {
                    break;
                }
            }
            out
        };
        let len = |i: u64| {
            format!("chr1\t{}\t.\tA\tG\t.\tPASS\t.\tGT{samples}", 1000 + i * 100).len() as u32
        };
        assert_eq!(
            query(Interval {
                start: 2099,
                end: 2100
            }),
            [(2099, len(11))]
        );
        let all = query(Interval {
            start: 0,
            end: 1 << 20,
        });
        assert_eq!(all.len(), 24);
        assert!(
            all.iter()
                .enumerate()
                .all(|(i, &(k, l))| k == 999 + 100 * i as u64 && l == len(i as u64))
        );
    }

    #[test]
    fn bin_starts_follow_the_csi_layout() {
        assert_eq!(bin_start(0, 14, 5), 0);
        assert_eq!(bin_start(1, 14, 5), 0);
        assert_eq!(bin_start(2, 14, 5), 1 << 26);
        assert_eq!(bin_start(4681, 14, 5), 0); // first leaf
        assert_eq!(bin_start(4682, 14, 5), 1 << 14);
        assert_eq!(bin_start(585 + 3, 14, 5), 3 << 17);
    }

    #[test]
    fn csi_depth_covers_the_largest_end() {
        assert_eq!(csi_depth(1), 5);
        assert_eq!(csi_depth(1 << 29), 5);
        assert_eq!(csi_depth((1 << 29) + 1), 6);
        assert_eq!(csi_depth(600_000_010), 6);
        assert_eq!(csi_depth(1 << 40), 9);
        assert_eq!(csi_depth(1 << 41), 9);
    }

    #[test]
    fn sparse_linear_index() {
        let vp = |n: u64| VirtualPosition::from(n);
        let mut l = Linear::default();
        l.add(10, vp(1)); // window 0
        l.add(100_000, vp(2)); // windows 1..=6
        l.add(50_000, vp(3)); // already covered
        l.add(200_000, vp(4)); // windows 7..=12
        assert_eq!(l.at(0), vp(1));
        assert_eq!(l.at(1), vp(2));
        assert_eq!(l.at(6), vp(2));
        assert_eq!(l.at(7), vp(4));
        assert_eq!(l.at(99), vp(4));
    }

    #[test]
    fn plan_windows_widen_merge_and_clamp() {
        let iv = |s, e| Interval { start: s, end: e };
        let p = Plan::new(
            0,
            b"c".to_vec(),
            vec![iv(10, 20), iv(0, 0), iv(21, 30), iv(100, 100), iv(10, 20)],
            1 << 29,
        );
        assert_eq!(p.set.regions().len(), 4);
        // [0,0) -> (0,1); [10,20) and [21,30) widen into touching windows
        assert_eq!(p.windows, [(0, 1), (9, 31), (99, 101)]);
        let p = Plan::new(
            0,
            b"c".to_vec(),
            vec![iv(5, u64::MAX), iv(1 << 30, 1 << 31)],
            1 << 29,
        );
        assert_eq!(p.windows, [(4, 1 << 29)]);
        assert!(p.set.matches(iv(7, 7), false));
        assert!(!p.set.matches(iv(5, 5), false));
        assert!(p.set.matches(iv(5, 5), true));
    }
}
