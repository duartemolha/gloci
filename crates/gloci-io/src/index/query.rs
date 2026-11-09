//! The index-backed region reader (§4.3, D9, D10).

use super::*;
use std::sync::Arc;

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
            let parsed = match class {
                Class::Fatal(msg) => Err(msg),
                _ if raw.len() > u32::MAX as usize => Err("line is longer than 4 GiB".into()),
                _ => self.parser.parse(raw, chroms),
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
/// errors name the record's BGZF virtual offset instead.
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
        index: Arc::new(index),
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
    }))
}

/// The reference names the index `info` of `data` knows, in reference id
/// order: the names stored in a TBI/CSI header, or, for a BCF's CSI (which
/// has none), the contig dictionary of the BCF header.
pub fn index_reference_names(data: &Path, info: &IndexInfo) -> Result<Vec<Vec<u8>>, IoError> {
    let loaded = Loaded::read(info)?;
    if let Some(names) = loaded.names() {
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
        .enumerate()
        .filter(|(rid, c)| !c.is_empty() && loaded.first_offset(*rid).is_some())
        .map(|(_, c)| c)
        .collect())
}

/// The chromosome names declared in an index that actually hold records
/// (at least one chunk offset is present).
pub fn index_non_empty_reference_names(
    data: &Path,
    info: &IndexInfo,
) -> Result<Vec<Vec<u8>>, IoError> {
    let loaded = Loaded::read(info)?;
    if let Some(names) = loaded.names() {
        return Ok(names
            .iter()
            .enumerate()
            .filter(|(rid, _)| loaded.first_offset(*rid).is_some())
            .map(|(_, n)| n.to_vec())
            .collect());
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
        let index = Loaded::read(&info)?;
        Ok(IndexedFile {
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

    /// Reference ids of the chromosomes, by name: the BCF contig
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
            _ => self.index.names().and_then(|n| n.get_index_of(chrom)),
        }
    }

    /// The reference names the index knows, in reference id order, as
    /// [`index_reference_names`] gives them.
    pub fn reference_names(&self) -> Vec<Vec<u8>> {
        if let Some(names) = self.index.names() {
            return names.iter().map(|n| n.to_vec()).collect();
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
        };
        s.map_or(0, |(lo, hi)| hi.saturating_sub(lo) + 1)
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
                        c.clone(),
                        vec![Interval {
                            start: 0,
                            end: max_pos,
                        }],
                        max_pos,
                    )
                })
            })
            .collect();
        let parser = match (&self.header, self.format) {
            (
                Header::Vcf(VcfHeader {
                    bcf: Some(dicts), ..
                }),
                Format::Bcf,
            ) => Parser::bcf(self.opts.vcf_span, dicts.clone()),
            _ => Parser::text(self.format, self.opts.vcf_span),
        };
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
            parser,
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
        }))
    }
}
