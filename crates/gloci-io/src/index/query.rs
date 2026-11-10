//! The index-backed region reader (§4.3, D9, D10).

use super::*;
use crate::align::{AlnHeader, AlnSource};
use crate::chrom_map::{ChromMap, Renamer};
use std::collections::HashMap;
use std::sync::Arc;

/// The largest decompressed CSI read into memory to check it.
const MAX_INDEX_BYTES: u64 = 4 << 30;

enum Loaded {
    Tbi(noodles_tabix::Index),
    Csi(noodles_csi::Index),
    Bai(noodles_bam::bai::Index),
}

impl Loaded {
    /// How many records precede reference `rid`'s in a coordinate-sorted BAM, from the
    /// index's per-reference counts: the number the streaming reader gives as "record M of
    /// the file". `None` when the index lacks the counts of an earlier reference.
    fn records_before(&self, rid: usize) -> Option<u64> {
        fn count<I: BinIndexKind>(refs: &[ReferenceSequence<I>], rid: usize) -> Option<u64> {
            refs.get(..rid)?.iter().try_fold(0u64, |n, r| {
                let m = r.metadata()?;
                Some(n + m.mapped_record_count() + m.unmapped_record_count())
            })
        }
        match self {
            Loaded::Bai(i) => count(i.reference_sequences(), rid),
            Loaded::Csi(i) => count(i.reference_sequences(), rid),
            Loaded::Tbi(_) => None,
        }
    }

    /// Load the index of the data file `data`, after checking its structure against the data
    /// file's size and, when it is known, the number of references its header lists.
    fn read(info: &IndexInfo, data: &Path, n_ref: Option<usize>) -> Result<Loaded, IoError> {
        let data_len = fs::metadata(data).map_or(u64::MAX, |m| m.len());
        let err = |source| IoError::Read {
            path: info.path.display().to_string(),
            source,
        };
        let corrupt = |msg: String| IoError::Detect {
            path: info.path.display().to_string(),
            msg: format!(
                "{}: {msg}; rebuild it with: gloci index {}",
                info.path.display(),
                info.path.with_extension("").display()
            ),
        };
        Ok(match info.kind {
            IndexKind::Tbi => {
                let bytes = validate::read_bgzf_capped(
                    File::open(&info.path).map_err(err)?,
                    MAX_INDEX_BYTES,
                )
                .map_err(err)?;
                validate::check_tbi(&bytes, data_len).map_err(corrupt)?;
                Loaded::Tbi(noodles_tabix::fs::read(&info.path).map_err(err)?)
            }
            IndexKind::Csi => {
                // noodles sizes its tables from the file's counts: check them first.
                let bytes = validate::read_bgzf_capped(
                    File::open(&info.path).map_err(err)?,
                    MAX_INDEX_BYTES,
                )
                .map_err(err)?;
                validate::check_csi(&bytes, data_len, n_ref).map_err(corrupt)?;
                Loaded::Csi(noodles_csi::fs::read(&info.path).map_err(err)?)
            }
            IndexKind::Bai => {
                let bytes = fs::read(&info.path).map_err(err)?;
                validate::check_bai(&bytes, data_len, n_ref).map_err(corrupt)?;
                Loaded::Bai(
                    noodles_bam::bai::io::Reader::new(&bytes[..])
                        .read_index()
                        .map_err(err)?,
                )
            }
            IndexKind::Crai => {
                return Err(IoError::Detect {
                    path: info.path.display().to_string(),
                    msg: format!(
                        "{} is a CRAM index, which this build cannot read yet; convert the CRAM \
file with: gloci view FILE.cram -o x.bam --write-index",
                        info.path.display()
                    ),
                });
            }
        })
    }

    fn names(&self) -> Option<&ReferenceSequenceNames> {
        let h = match self {
            Loaded::Tbi(i) => i.header(),
            Loaded::Csi(i) => i.header(),
            Loaded::Bai(_) => None,
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
            Loaded::Bai(i) => i.reference_sequences().get(rid).and_then(first),
        }
    }

    /// One past the largest position the index can hold.
    fn max_position(&self) -> u64 {
        let (min_shift, depth) = match self {
            Loaded::Tbi(i) => (i.min_shift(), i.depth()),
            Loaded::Csi(i) => (i.min_shift(), i.depth()),
            Loaded::Bai(i) => (i.min_shift(), i.depth()),
        };
        1u64.checked_shl(u32::from(min_shift) + 3 * u32::from(depth))
            .unwrap_or(u64::MAX)
    }

    /// Chunks, in file order, that hold every record whose index span
    /// overlaps the 0-based window `[lo, hi)`.
    fn chunks(&self, rid: usize, lo: u64, hi: u64) -> Vec<Chunk> {
        match self {
            Loaded::Bai(i) => {
                let Some(rs) = i.reference_sequences().get(rid) else {
                    return Vec::new();
                };
                let min_off = rs.min_offset(i.min_shift(), i.depth(), position(lo + 1));
                bin_chunks(rs, i.min_shift(), i.depth(), lo, hi, min_off)
            }
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
pub(super) struct Plan {
    rid: usize,
    name: Vec<u8>,
    /// Records are tested against each region, not against their union: a
    /// point between two bookended regions is in neither (§4.3).
    pub(super) set: RegionSet,
    /// Disjoint query windows in coordinate order: each region widened by
    /// one base on both sides (so points on region edges, and zero-length
    /// records indexed as `[s, s+1)`, are found), merged, and clamped to
    /// the index's range. 0-based, half-open.
    pub(super) windows: Vec<(u64, u64)>,
}

impl Plan {
    pub(super) fn new(rid: usize, name: Vec<u8>, regions: Vec<Interval>, max_pos: u64) -> Plan {
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
    /// `--chrom-map` for text records (a BCF's contigs are renamed in its
    /// header, which `Parser::Bcf` already reads).
    renamer: Option<Renamer>,
    index: Arc<Loaded>,
    plans: Vec<Plan>,
    /// Whole-chromosome reads: every record of each plan's reference is
    /// emitted, with no footprint test, and plans are visited in the order
    /// given (each reference's records lie in one stretch of the file).
    whole: bool,
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
    /// Records emitted so far: each record's `loc` is its number in this
    /// reader's output, which keeps file order for tie-breaking (D38).
    emitted: u64,
    /// Data records of the current reference read so far (each once, even
    /// when chunks overlap) and the position of the last counted one. Only
    /// whole-reference reads start at the reference's first record, so only
    /// there is it the record's number within the chromosome.
    seen: u64,
    seen_last: Option<VirtualPosition>,
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
        if let Parser::Bam(_) = self.parser {
            return self.read_bam_raw(buf, off);
        }
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

    /// Append the next BAM record, without its `block_size`.
    fn read_bam_raw(&mut self, buf: &mut Vec<u8>, off: usize) -> Result<bool, IoError> {
        let got = bcf::append(&mut self.src, buf, 4).map_err(|e| self.read_err(e))?;
        if got == 0 {
            return Ok(false);
        }
        if got < 4 {
            return self.bad_bam_frame(
                buf,
                off,
                "record is truncated; the file was probably cut off",
            );
        }
        let block = i32::from_le_bytes([buf[off], buf[off + 1], buf[off + 2], buf[off + 3]]);
        buf.truncate(off);
        let Ok(size) = usize::try_from(block) else {
            return self.bad_bam_frame(
                buf,
                off,
                &format!("block_size {block} is negative, so the next record cannot be found"),
            );
        };
        let got = bcf::append(&mut self.src, buf, size).map_err(|e| self.read_err(e))?;
        if got < size {
            return self.bad_bam_frame(
                buf,
                off,
                "record is truncated; the file was probably cut off",
            );
        }
        Ok(true)
    }

    /// A BAM record whose frame cannot be read: its place is its BGZF offset, or (in a
    /// whole-reference read) its number in the reference, as for other record errors.
    fn bad_bam_frame(&mut self, buf: &mut Vec<u8>, off: usize, msg: &str) -> Result<bool, IoError> {
        buf.truncate(off);
        let at = self.src.virtual_position();
        self.seen += 1;
        self.bad(at, &format!("{msg}; the file is corrupt"))?;
        // --skip-invalid cannot resynchronise a broken frame: stop reading this file.
        Ok(false)
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
    /// Index-backed reads have no line numbers. A whole-reference read names
    /// the record by its number within its chromosome ("record N of CHROM");
    /// a region read, which starts mid-file, by its BGZF virtual offset.
    fn bad(&mut self, at: VirtualPosition, msg: &str) -> Result<(), IoError> {
        let place = match self.plans.get(self.p) {
            Some(plan) if self.whole => {
                // plans name the file's own chromosome; show the new name
                let shown = self
                    .renamer
                    .as_ref()
                    .map_or(&plan.name[..], |r| r.map().rename(&plan.name));
                let shown = String::from_utf8_lossy(shown);
                // A BAM's streaming reader words it so, with the number in the whole file.
                match (&self.parser, self.index.records_before(plan.rid)) {
                    (Parser::Bam(_), Some(before)) => format!(
                        "record {} of {shown} (record {} of the file)",
                        gloci_core::thousands(self.seen),
                        gloci_core::thousands(before + self.seen)
                    ),
                    _ => format!("record {} of {shown}", self.seen),
                }
            }
            _ => {
                let (c, u): (u64, u16) = at.into();
                format!("record at BGZF offset {c}:{u}")
            }
        };
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
                self.seen = 0;
                self.seen_last = None;
                if self.whole {
                    // A new reference: nothing of it has been emitted, and its
                    // records may lie before the previous reference's.
                    self.last_emitted = None;
                }
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
            if let Some(r) = &mut self.renamer {
                r.apply(&mut batch.buf, off);
            }
            if self.seen_last.is_none_or(|l| start > l) {
                self.seen += 1;
                self.seen_last = Some(start);
            }
            let raw_len = batch.buf.len() - off;
            let parsed = match class {
                Class::Fatal(msg) => Err(msg),
                _ if raw_len > u32::MAX as usize => Err("line is longer than 4 GiB".into()),
                _ => match &mut self.parser {
                    Parser::Bam(b) => {
                        match b.parse(&self.name, &mut batch.buf, off, chroms, &mut self.diag) {
                            Ok(Some(p)) => Ok(p),
                            Ok(None) => {
                                batch.buf.truncate(off);
                                continue;
                            }
                            Err(msg) => Err(msg),
                        }
                    }
                    p => p.parse(&batch.buf[off..], chroms),
                },
            };
            let (chrom, key, span, strand) = match parsed {
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
            if (self.whole || plan.set.matches(span, self.inclusive))
                && self.last_emitted.is_none_or(|l| start > l)
            {
                let len = (batch.buf.len() - off) as u32;
                self.emitted += 1;
                batch.recs.push(RecordMeta {
                    chrom,
                    key,
                    span,
                    strand,
                    loc: self.emitted,
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

/// Index-backed [`AlnSource`]: the BAM records of some regions.
struct IndexedAlnSource {
    aln: Arc<AlnHeader>,
    reader: RegionReader,
    chroms: ChromTable,
    batch: RecordBatch,
}

impl AlnSource for IndexedAlnSource {
    fn name(&self) -> &str {
        &self.reader.name
    }
    fn layout(&self) -> Format {
        Format::Bam
    }
    fn header(&self) -> &Arc<AlnHeader> {
        &self.aln
    }
    fn diagnostics(&self) -> &Diagnostics {
        &self.reader.diag
    }

    fn next_record(&mut self, buf: &mut Vec<u8>) -> Result<Option<u64>, IoError> {
        self.batch.clear();
        if !self.reader.next(&mut self.chroms, &mut self.batch)? {
            return Ok(None);
        }
        buf.extend_from_slice(self.batch.raw(0));
        Ok(Some(self.reader.emitted))
    }
}

/// The name the file itself gives the chromosome that `--chrom-map` (when
/// given) calls `requested`, which is what an index lists and the file's
/// text records say; `None` when no chromosome of the file is called
/// that. This is how a region's chromosome, named as the records are
/// after renaming, is turned back into the file's own name before the
/// index is queried.
pub fn file_chrom_name<'a>(map: Option<&'a ChromMap>, requested: &'a [u8]) -> Option<&'a [u8]> {
    match map {
        Some(m) => m.file_name(requested),
        None => Some(requested),
    }
}

/// A renamer for the text records of a `format` file, when `o` has a map.
fn text_renamer(format: Format, o: &ReaderOptions) -> Option<Renamer> {
    match (&o.chrom_map, format) {
        (Some(m), Format::Bed | Format::BedGraph | Format::Gff3 | Format::Gtf | Format::Vcf) => {
            Some(Renamer::new(m.clone()))
        }
        _ => None,
    }
}

/// WIG lines are not BED, so an index next to a WIG file cannot be used.
fn wig_unindexed(name: &str) -> IoError {
    UsageError(format!(
        "{name} is WIG, which cannot be read through an index; convert it with \
`gloci convert {name}`, which writes bedGraph, and index that"
    ))
    .into()
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
/// Line numbers are unknown, so each record's `loc` is its number among the
/// records this reader yields (file order, which ties are broken by); parse
/// errors name the record's BGZF virtual offset instead (whole-chromosome
/// readers name "record N of CHROM").
pub fn open_region_reader(
    path: &Path,
    regions: &[(Vec<u8>, Interval)],
    o: &ReaderOptions,
    inclusive_points: bool,
) -> Result<Box<dyn RecordReader>, IoError> {
    let name = path.display().to_string();
    check_bgzf(path, &name)?;
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
    IndexedFile::open(path, format, header, o)?.regions_reader(regions, inclusive_points)
}

/// The reference names of a file whose index lists none, by reference id: a BCF's contig
/// dictionary (unused ids are empty) or a BAM's reference list.
fn dictionary_names(data: &Path, info: &IndexInfo) -> Result<(Vec<Vec<u8>>, bool), IoError> {
    let name = data.display().to_string();
    let file = File::open(data).map_err(|source| IoError::Open {
        path: name.clone(),
        source,
    })?;
    let mut src = noodles_bgzf::io::Reader::new(BufReader::new(file));
    let is_bam = src
        .fill_buf()
        .map_err(|source| IoError::Read {
            path: name.clone(),
            source,
        })?
        .starts_with(b"BAM\x01");
    if is_bam {
        let h = AlnHeader::read_bam(&mut src).map_err(|msg| IoError::Detect {
            msg: format!("{name}: {msg}"),
            path: name.clone(),
        })?;
        return Ok((h.refs.into_iter().map(|(n, _)| n).collect(), true));
    }
    let h = read_bcf_header(&mut src).map_err(|msg| IoError::Detect {
        msg: format!(
            "{name}: the index {} lists no chromosome names, as only a BCF or BAM index may, \
but the header cannot be read: {msg}",
            info.path.display()
        ),
        path: name.clone(),
    })?;
    let dicts = h.bcf.expect("read_bcf_header builds the dictionaries");
    Ok((dicts.contigs, false))
}

/// The reference names the index `info` of `data` knows, in reference id
/// order: the names stored in a TBI/CSI header, or, for a BCF's or BAM's
/// index (which has none), the contig dictionary of the BCF header or the
/// reference list of the BAM header.
pub fn index_reference_names(data: &Path, info: &IndexInfo) -> Result<Vec<Vec<u8>>, IoError> {
    let loaded = Loaded::read(info, data, None)?;
    if let Some(names) = loaded.names() {
        return Ok(names.iter().map(|n| n.to_vec()).collect());
    }
    let (names, is_bam) = dictionary_names(data, info)?;
    Ok(names
        .into_iter()
        .enumerate()
        .filter(|(rid, c)| !c.is_empty() && (is_bam || loaded.first_offset(*rid).is_some()))
        .map(|(_, c)| c)
        .collect())
}

/// The chromosome names declared in an index that actually hold records
/// (at least one chunk offset is present).
pub fn index_non_empty_reference_names(
    data: &Path,
    info: &IndexInfo,
) -> Result<Vec<Vec<u8>>, IoError> {
    let loaded = Loaded::read(info, data, None)?;
    if let Some(names) = loaded.names() {
        return Ok(names
            .iter()
            .enumerate()
            .filter(|(rid, _)| loaded.first_offset(*rid).is_some())
            .map(|(_, n)| n.to_vec())
            .collect());
    }
    Ok(dictionary_names(data, info)?
        .0
        .into_iter()
        .enumerate()
        .filter(|(rid, c)| !c.is_empty() && loaded.first_offset(*rid).is_some())
        .map(|(_, c)| c)
        .collect())
}

/// An indexed BGZF file whose header and index are loaded once, from which
/// whole-chromosome readers are made cheaply (per-chromosome parallelism,
/// D41). The readers share the index.
pub struct IndexedFile {
    path: PathBuf,
    name: String,
    format: Format,
    header: Header,
    info: IndexInfo,
    index: Arc<Loaded>,
    opts: ReaderOptions,
    /// Reference ids by name as the records have it after `--chrom-map`,
    /// for text formats (built forward: the map is not reversed).
    by_new_name: HashMap<Vec<u8>, usize>,
    /// The readers made from this file keep their one-time warnings in their diagnostics.
    defer_warnings: bool,
}

impl IndexedFile {
    /// Load the index of `path` (with `format` and `header` already read
    /// from it, for example by the reader the command opened). The readers
    /// made from it parse records with `o`.
    pub fn open(
        path: &Path,
        format: Format,
        header: Header,
        o: &ReaderOptions,
    ) -> Result<IndexedFile, IoError> {
        let name = path.display().to_string();
        check_bgzf(path, &name)?;
        let info = find_index(path).ok_or_else(|| {
            UsageError(format!(
                "{name} has no index; create one with: gloci index {name}"
            ))
        })?;
        if format == Format::Wig {
            return Err(wig_unindexed(&name));
        }
        // A tabix/CSI index describes interval text; the other alignment, sequence and
        // convert-only formats have no text parser here (`Parser::text`). BAM has its own.
        if (format.is_alignment() && format != Format::Bam)
            || format.is_convert_only()
            || format == Format::Fasta
        {
            return Err(UsageError(format!(
                "{name} is {}, which cannot be read through a tabix or CSI index as intervals; \
convert it to BED first: gloci convert {name} -o x.bed.gz",
                format.name()
            ))
            .into());
        }
        let n_ref = match (&header, format) {
            (Header::Alignment(a), Format::Bam) => Some(a.refs.len()),
            _ => None,
        };
        let index = Loaded::read(&info, path, n_ref)?;
        let by_new_name = if format == Format::Bam {
            let Header::Alignment(aln) = &header else {
                return Err(IoError::Internal(
                    "IndexedFile::open of a BAM needs the alignment header".into(),
                ));
            };
            // A map renames the header's references, which `ref_id` does not know by
            // their new names; without one, `ref_id` answers.
            if o.chrom_map.is_some() {
                aln.refs
                    .iter()
                    .enumerate()
                    .map(|(rid, (n, _))| (n.clone(), rid))
                    .collect()
            } else {
                HashMap::new()
            }
        } else {
            if let (Some(m), Some(names)) = (&o.chrom_map, index.names()) {
                m.check_file_names(&name, names.iter().map(|n| &n[..]))?;
            }
            let map = o.chrom_map.as_deref();
            index
                .names()
                .map(|names| {
                    names
                        .iter()
                        .enumerate()
                        .map(|(rid, n)| (map.map_or(&n[..], |m| m.rename(n)).to_vec(), rid))
                        .collect()
                })
                .unwrap_or_default()
        };
        Ok(IndexedFile {
            defer_warnings: false,
            by_new_name,
            path: path.to_path_buf(),
            name,
            format,
            header,
            info,
            index: Arc::new(index),
            opts: o.clone(),
        })
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn info(&self) -> &IndexInfo {
        &self.info
    }

    /// Have the readers made from here keep the one-time warnings of an alignment reader (a
    /// read past its reference) in their diagnostics, for the caller to print once for all
    /// chunks of a run (`Diagnostics::log_deferred_warnings`).
    pub fn defer_warnings(&mut self) {
        self.defer_warnings = true;
    }

    /// The reads with no place the index counts (BAM): the unmapped reads at the end of a
    /// coordinate-sorted file, which no chromosome's reader reaches. 0 for other indexes.
    pub fn unplaced_reads(&self) -> u64 {
        if self.format != Format::Bam {
            return 0;
        }
        match &*self.index {
            Loaded::Bai(i) => i.unplaced_unmapped_record_count().unwrap_or(0),
            Loaded::Csi(i) => i.unplaced_unmapped_record_count().unwrap_or(0),
            Loaded::Tbi(_) => 0,
        }
    }

    /// Reference ids of the chromosomes, by name (as renamed by
    /// `--chrom-map`): the BCF contig
    /// dictionary for a BCF's CSI (which stores no names), else the index's
    /// own names.
    fn rid_of(&self, chrom: &[u8]) -> Option<usize> {
        match &self.header {
            Header::Vcf(VcfHeader {
                bcf: Some(dicts), ..
            }) if self.format == Format::Bcf => dicts
                .contigs
                .iter()
                .position(|c| !c.is_empty() && c == chrom),
            // Without `--chrom-map` the header's own O(1) name map answers.
            Header::Alignment(aln) if self.opts.chrom_map.is_none() => aln.ref_id(chrom),
            _ => self.by_new_name.get(chrom).copied(),
        }
    }

    /// The name the file itself gives reference `rid` (before `--chrom-map`);
    /// BCF's dictionary is not renamed by the map.
    fn own_name(&self, rid: usize, renamed: &[u8]) -> Vec<u8> {
        if self.format == Format::Bcf || self.format == Format::Bam {
            return renamed.to_vec();
        }
        self.index
            .names()
            .and_then(|n| n.get_index(rid))
            .map_or_else(|| renamed.to_vec(), |n| n.to_vec())
    }

    /// The reference names the index knows, in reference id order, as
    /// [`index_reference_names`] gives them.
    pub fn reference_names(&self) -> Vec<Vec<u8>> {
        if let Header::Alignment(aln) = &self.header {
            return aln.refs.iter().map(|(n, _)| n.clone()).collect();
        }
        if let Some(names) = self.index.names() {
            let map = self.opts.chrom_map.as_deref();
            return names
                .iter()
                .map(|n| map.map_or(&n[..], |m| m.rename(n)).to_vec())
                .collect();
        }
        match &self.header {
            Header::Vcf(VcfHeader {
                bcf: Some(dicts), ..
            }) => dicts
                .contigs
                .iter()
                .enumerate()
                .filter(|(rid, c)| !c.is_empty() && self.index.first_offset(*rid).is_some())
                .map(|(_, c)| c.clone())
                .collect(),
            _ => Vec::new(),
        }
    }

    /// The reference names that hold records, in the order their records
    /// appear in the file (by first chunk offset), which is the index's own
    /// order for a TBI but the file's order, not the dictionary's, for a BCF.
    pub fn file_order_reference_names(&self) -> Vec<Vec<u8>> {
        let mut named: Vec<(VirtualPosition, Vec<u8>)> = self
            .reference_names()
            .into_iter()
            .filter_map(|n| {
                let at = self.index.first_offset(self.rid_of(&n)?)?;
                Some((at, n))
            })
            .collect();
        named.sort_by_key(|(at, _)| *at);
        named.into_iter().map(|(_, n)| n).collect()
    }

    /// The reference names that hold records.
    pub fn non_empty_reference_names(&self) -> Vec<Vec<u8>> {
        self.reference_names()
            .into_iter()
            .filter(|n| self.weight(n) > 0)
            .collect()
    }

    /// An estimate of how much of the file `chrom`'s records take: the
    /// compressed bytes between its first and last chunk offsets, plus one;
    /// 0 when the index holds no records for it.
    pub fn weight(&self, chrom: &[u8]) -> u64 {
        let Some(rid) = self.rid_of(chrom) else {
            return 0;
        };
        fn span<I: BinIndexKind>(rs: &ReferenceSequence<I>) -> Option<(u64, u64)> {
            let chunks = rs.bins().values().flat_map(|b| b.chunks());
            let lo = chunks.clone().map(|c| c.start().compressed()).min()?;
            let hi = chunks.map(|c| c.end().compressed()).max()?;
            Some((lo, hi))
        }
        let s = match &*self.index {
            Loaded::Tbi(i) => i.reference_sequences().get(rid).and_then(span),
            Loaded::Csi(i) => i.reference_sequences().get(rid).and_then(span),
            Loaded::Bai(i) => i.reference_sequences().get(rid).and_then(span),
        };
        s.map_or(0, |(lo, hi)| hi.saturating_sub(lo) + 1)
    }

    /// The records that overlap `regions` (named as records are after
    /// `--chrom-map`), each once, in file order, as [`open_region_reader`]
    /// reads them. Regions on chromosomes the index does not know yield
    /// nothing.
    pub fn regions_reader(
        &self,
        regions: &[(Vec<u8>, Interval)],
        inclusive_points: bool,
    ) -> Result<Box<dyn RecordReader>, IoError> {
        Ok(Box::new(self.region_reader(
            regions,
            inclusive_points,
            false,
        )?))
    }

    /// The parser of this file's records. `all` is for BAM: keep every record, placed
    /// unmapped reads too, tested by the span samtools finds them by (`view`, D116).
    fn parser(&self, all: bool) -> Result<Parser, IoError> {
        Ok(match (&self.header, self.format) {
            (
                Header::Vcf(VcfHeader {
                    bcf: Some(dicts), ..
                }),
                Format::Bcf,
            ) => Parser::bcf(self.opts.vcf_span, dicts.clone()),
            (Header::Alignment(aln), Format::Bam) => {
                Parser::Bam(Box::new(crate::align::BamRegionParse::new(
                    aln.clone(),
                    self.opts.accept.alignments.clone().unwrap_or_default(),
                    all,
                    self.defer_warnings,
                )))
            }
            (_, Format::Bam) => {
                return Err(IoError::Internal(
                    "a BAM index reader needs the alignment header".into(),
                ));
            }
            _ => Parser::text(self.format, self.opts.vcf_span),
        })
    }

    /// The records of the BAM file that overlap `regions`, each once, in file order, as
    /// every-record [`AlnSource`] reads them: unmapped reads placed on a region included,
    /// and the span a region is tested against is samtools' (D116). The read filters apply;
    /// record numbers count the records yielded.
    pub fn aln_regions_source(
        &self,
        regions: &[(Vec<u8>, Interval)],
        inclusive_points: bool,
    ) -> Result<Box<dyn AlnSource>, IoError> {
        let Header::Alignment(aln) = &self.header else {
            return Err(IoError::Internal(
                "an alignment source needs a BAM file".into(),
            ));
        };
        Ok(Box::new(IndexedAlnSource {
            aln: aln.clone(),
            reader: self.region_reader(regions, inclusive_points, true)?,
            chroms: ChromTable::first_appearance(),
            batch: RecordBatch::default(),
        }))
    }

    fn region_reader(
        &self,
        regions: &[(Vec<u8>, Interval)],
        inclusive_points: bool,
        all: bool,
    ) -> Result<RegionReader, IoError> {
        log::debug!(
            "reading {} through its index {} for {} region(s)",
            self.name,
            self.info.path.display(),
            regions.len()
        );
        let mut by_rid: BTreeMap<usize, (Vec<u8>, Vec<Interval>)> = BTreeMap::new();
        for (chrom, iv) in regions {
            let Some(rid) = self.rid_of(chrom) else {
                continue;
            };
            // BCF's dictionary is already renamed; text uses the file's name.
            by_rid
                .entry(rid)
                .or_insert_with(|| (self.own_name(rid, chrom), Vec::new()))
                .1
                .push(*iv);
        }
        let max_pos = self.index.max_position();
        let mut plans: Vec<Plan> = by_rid
            .into_iter()
            .map(|(rid, (name, ivs))| Plan::new(rid, name, ivs, max_pos))
            .collect();
        // File order, so that `last_emitted` only ever drops a record already
        // emitted. References without records yield nothing; their place is
        // irrelevant.
        plans.sort_by_key(|p| {
            self.index
                .first_offset(p.rid)
                .unwrap_or(VirtualPosition::MAX)
        });
        let file = File::open(&self.path).map_err(|source| IoError::Open {
            path: self.name.clone(),
            source,
        })?;
        Ok(RegionReader {
            name: self.name.clone(),
            format: self.format,
            header: self.header.clone(),
            diag: Diagnostics::default(),
            policy: self.opts.invalid,
            inclusive: inclusive_points,
            src: noodles_bgzf::io::Reader::new(BufReader::new(file)),
            renamer: text_renamer(self.format, &self.opts),
            parser: self.parser(all)?,
            index: self.index.clone(),
            plans,
            whole: false,
            p: 0,
            w: 0,
            chunks: Vec::new(),
            c: 0,
            loaded: false,
            chunk_end: None,
            last_emitted: None,
            emitted: 0,
            seen: 0,
            seen_last: None,
        })
    }

    /// Every record on `chroms`, chromosome by chromosome in the order
    /// given, each chromosome's records in file order. There is no footprint
    /// test, so zero-length records at position 0 are read too. Chromosomes
    /// the index does not know yield nothing.
    pub fn chromosomes_reader(&self, chroms: &[Vec<u8>]) -> Result<Box<dyn RecordReader>, IoError> {
        let max_pos = self.index.max_position();
        let plans: Vec<Plan> = chroms
            .iter()
            .filter_map(|c| {
                self.rid_of(c).map(|rid| {
                    Plan::new(
                        rid,
                        self.own_name(rid, c),
                        vec![Interval {
                            start: 0,
                            end: max_pos,
                        }],
                        max_pos,
                    )
                })
            })
            .collect();
        let file = File::open(&self.path).map_err(|source| IoError::Open {
            path: self.name.clone(),
            source,
        })?;
        Ok(Box::new(RegionReader {
            name: self.name.clone(),
            format: self.format,
            header: self.header.clone(),
            diag: Diagnostics::default(),
            policy: self.opts.invalid,
            inclusive: true,
            src: noodles_bgzf::io::Reader::new(BufReader::new(file)),
            renamer: text_renamer(self.format, &self.opts),
            parser: self.parser(false)?,
            index: self.index.clone(),
            plans,
            whole: true,
            p: 0,
            w: 0,
            chunks: Vec::new(),
            c: 0,
            loaded: false,
            chunk_end: None,
            last_emitted: None,
            emitted: 0,
            seen: 0,
            seen_last: None,
        }))
    }
}
