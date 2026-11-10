//! Index building (D20): scanning, TBI/CSI construction, staged publish.

use super::*;

/// A BGZF file with its end-of-file marker; anything else is refused.
pub(super) fn check_bgzf(data: &Path, name: &str) -> Result<(), IoError> {
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
        // A .gz name is rewritten in place (plain gzip is not seekable, so
        // it cannot be indexed); anything else gets a .gz name.
        let out = if name.ends_with(".gz") {
            name.to_string()
        } else {
            format!("{name}.gz")
        };
        return Err(UsageError(format!(
            "{name} is not BGZF-compressed; recompress with: gloci view {name} -o {out}"
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

/// One record to index.
pub(super) struct Entry {
    rid: usize,
    begin: u64,
    end: u64,
    chunk: Chunk,
    loc: u64,
}

/// What a scan learned besides the records.
pub(super) struct Scanned {
    /// Reference names in id order (text formats only).
    names: Vec<Vec<u8>>,
    ref_count: usize,
    /// Leading lines htslib must skip: up to the last leading BED
    /// `track`/`browser` line (`#` lines are skipped as comments anyway).
    skip: u32,
    max_end: u64,
    /// The most raw record bytes held at once.
    #[cfg(test)]
    pub(super) peak_bytes: usize,
}

/// Reads records with their virtual positions, checks order in
/// first-appearance mode and hands them to `sink` in batches.
struct Scanner<'a> {
    name: &'a str,
    binary: bool,
    /// Columns of a generic TSV being indexed (`format` is then nominal).
    tsv: Option<&'a TsvColumns>,
    chroms: ChromTable,
    guard: OrderGuard,
    batch: RecordBatch,
    entries: Vec<Entry>,
    max_end: u64,
    #[cfg(test)]
    pub(super) peak_bytes: usize,
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
            strand: Strand::Unknown, // index scans never read strand
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
        self.guard.check(&self.chroms, &self.batch).map_err(|e| {
            let e: IoError = e;
            match e {
                // The records of a BCF can be contiguous per contig and still
                // fail here: the index needs the contigs in dictionary order.
                IoError::Order(mut o) if self.tsv.is_some() => {
                    o.hint = tsv_sort_hint(&o.path, self.tsv.expect("checked above"));
                    IoError::Order(o)
                }
                IoError::Order(mut o) if self.binary => {
                    o.hint.push_str(
                        "; an index of a BCF needs its contigs in the order of the header's \
contig dictionary, even when each contig is contiguous (bcftools index is less strict)",
                    );
                    IoError::Order(o)
                }
                e => e,
            }
        })?;
        for e in &self.entries {
            (self.sink)(e)?;
        }
        self.batch.clear();
        self.entries.clear();
        Ok(())
    }

    /// Text formats; returns the line skip count.
    fn text(&mut self, r: &mut MultithreadedReader<File>, format: Format) -> Result<u32, IoError> {
        let mut parser = match self.tsv {
            Some(cols) => Parser::Tsv(TsvKind::new(cols)),
            None => Parser::text(format, VcfSpan::Ref),
        };
        let skip_lines = self.tsv.map_or(0, |c| c.skip_lines as u64);
        let bed = self.tsv.is_none() && matches!(format, Format::Bed | Format::BedGraph);
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
            if n <= skip_lines {
                skip = u32::try_from(n).unwrap_or(u32::MAX);
                self.batch.buf.truncate(off);
                continue;
            }
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
        // dictionary order. That is stricter than htslib: the noodles indexer
        // needs increasing reference ids, where `bcftools index` accepts
        // contigs that are contiguous in any order.
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
pub(super) fn scan(
    data: &Path,
    name: &str,
    format: Format,
    tsv: Option<&TsvColumns>,
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
        tsv,
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
pub(super) struct Linear {
    runs: Vec<(u64, VirtualPosition)>,
    covered: u64,
}

impl Linear {
    pub(super) fn add(&mut self, end: u64, start: VirtualPosition) {
        let last = (end - 1) >> MIN_SHIFT;
        if last >= self.covered {
            self.runs.push((self.covered, start));
            self.covered = last + 1;
        }
    }

    pub(super) fn at(&self, window: u64) -> VirtualPosition {
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

fn text_header(format: Format, tsv: Option<&TsvColumns>, s: &Scanned) -> Option<CsiHeader> {
    let builder = match (format, tsv) {
        (_, Some(cols)) => cols.header_builder(),
        (Format::Bed | Format::BedGraph, _) => HeaderBuilder::bed().set_line_skip_count(s.skip),
        // WIG is rejected before indexing (`build_index_with`).
        (Format::Wig, _) => return None,
        (Format::Gff3 | Format::Gtf, _) => HeaderBuilder::gff(),
        (Format::Vcf, _) => HeaderBuilder::vcf(),
        (Format::Bcf, _) => return None,
        // Refused before indexing (`build_index_with`).
        (Format::Bam | Format::Cram | Format::Sam | Format::Fasta, _) => return None,
        (Format::Psl | Format::Rmsk, _) => return None,
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
    tsv: Option<&TsvColumns>,
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
    let scanned = scan(data, name, format, tsv, threads, &mut sink)?;
    let index = if full {
        None
    } else {
        let indexer = match text_header(format, tsv, &scanned) {
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
pub(super) fn csi_depth(max_end: u64) -> u8 {
    let mut d = TBI_DEPTH;
    while (1u64 << (u32::from(MIN_SHIFT) + 3 * u32::from(d))) < max_end {
        d += 1;
    }
    d
}

/// First position covered by bin `id`.
pub(super) fn bin_start(id: usize, min_shift: u8, depth: u8) -> u64 {
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
pub(super) fn with_htslib_loffsets(
    index: noodles_csi::Index,
    linear: &[Linear],
) -> noodles_csi::Index {
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

pub(super) fn csi_bytes(index: &noodles_csi::Index) -> io::Result<Vec<u8>> {
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
    build_index_with(data, final_data_path, format, None, threads)
}

/// `build_index` for a generic TSV when `tsv` is given (`format` is then
/// only nominal and must be a text format).
pub(super) fn build_index_with(
    data: &Path,
    final_data_path: &Path,
    format: Format,
    tsv: Option<&TsvColumns>,
    threads: NonZeroUsize,
) -> Result<(IndexKind, Vec<u8>), IoError> {
    let name = final_data_path.display().to_string();
    if format == Format::Wig && tsv.is_none() {
        return Err(UsageError(format!(
            "{name} is WIG, which cannot be indexed; convert it first with \
`gloci convert {name}`, which writes bedGraph, and index that"
        ))
        .into());
    }
    if tsv.is_none() && format == Format::Bam {
        return super::bam::build_bam_index(data, final_data_path, false, threads);
    }
    if tsv.is_none() && format == Format::Sam {
        return Err(UsageError(format!(
            "{name} is SAM text, which has no index; write BAM instead, with: \
gloci view {name} -o x.bam --write-index"
        ))
        .into());
    }
    if tsv.is_none() && (format.is_alignment() || format == Format::Fasta) {
        // CRAM (CRAI) and FASTA (FAI) are indexed by their own commands.
        return Err(UsageError(format!(
            "{name} is {}, and indexing it is not available in this build yet",
            format.name()
        ))
        .into());
    }
    if tsv.is_none() && format.is_convert_only() {
        return Err(UsageError(format!(
            "{name} is {}, which `gloci index` does not index; convert it to BED first with \
`gloci convert {name} -o out.bed.gz` and index that",
            format.name()
        ))
        .into());
    }
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
        let p = pass::<BinnedIndex>(data, &name, format, tsv, threads, TBI_DEPTH)?;
        let p = match p.index {
            Some(_) => p,
            None => {
                let depth = csi_depth(p.scanned.max_end);
                pass::<BinnedIndex>(data, &name, format, tsv, threads, depth)?
            }
        };
        (IndexKind::Csi, csi_bytes(&csi(p)?))
    } else {
        let p = pass::<LinearIndex>(data, &name, format, tsv, threads, TBI_DEPTH)?;
        match p.index {
            Some(index) => (IndexKind::Tbi, tbi_bytes(&index)),
            None => {
                let depth = csi_depth(p.scanned.max_end);
                let q = pass::<BinnedIndex>(data, &name, format, tsv, threads, depth)?;
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
    let mut af = AtomicFile::create(path).map_err(|source| IoError::Create {
        path: shown(),
        source,
    })?;
    af.write_all(bytes).map_err(|source| IoError::Write {
        path: shown(),
        source,
    })?;
    Ok(af)
}

/// Remove the indexes `data` can have (`DATA.tbi`, `DATA.csi`, and for a `.bam` also
/// `DATA.bai` and `x.bai`) if they exist (D19).
pub(crate) fn remove_indexes(data: &Path) -> Result<(), IoError> {
    remove_paths(all_index_paths(data))
}

fn remove_paths(paths: Vec<PathBuf>) -> Result<(), IoError> {
    for p in paths {
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
/// write it atomically to `final_data_path` + `.tbi`/`.csi`/`.bai`, remove an index
/// of the other kind, and return the index path. Messages name
/// `final_data_path`.
pub fn write_index(
    data: &Path,
    final_data_path: &Path,
    format: Format,
    threads: NonZeroUsize,
) -> Result<PathBuf, IoError> {
    let (kind, bytes) = build_index(data, final_data_path, format, threads)?;
    publish_index(final_data_path, kind, &bytes)
}

/// Write `bytes` atomically as the index of `final_data_path` and remove an
/// index of the other kind; returns the index path.
pub(super) fn publish_index(
    final_data_path: &Path,
    kind: IndexKind,
    bytes: &[u8],
) -> Result<PathBuf, IoError> {
    let path = index_path(final_data_path, kind);
    let af = stage_index(&path, bytes)?;
    af.commit().map_err(|source| IoError::Write {
        path: path.display().to_string(),
        source,
    })?;
    // Any other index of this file describes data that no longer is there, or is shadowed.
    remove_paths(
        all_index_paths(final_data_path)
            .into_iter()
            .filter(|p| *p != path)
            .collect(),
    )?;
    Ok(path)
}

// ---- index-backed region reader ----

/// How to sort a TSV by its index columns while keeping the skipped leading
/// lines and the comment lines at the top, where the index expects them.
fn tsv_sort_hint(path: &str, c: &TsvColumns) -> String {
    let comment = match c.comment {
        b'\'' => "\"'\"".to_string(),
        b => format!("'{}'", b as char),
    };
    let body = format!("zcat {path} | tail -n +{}", c.skip_lines + 1);
    let head = if c.skip_lines > 0 {
        format!("zcat {path} | head -n {}; ", c.skip_lines)
    } else {
        String::new()
    };
    format!(
        "sort by chromosome (each one in a single block) and then by the begin column, keeping \
the first lines and the comment lines on top: ({head}{body} | grep '^'{comment}; {body} | \
grep -v '^'{comment} | LC_ALL=C sort -t $'\\t' -k{s},{s} -k{b},{b}n) | bgzip > {path}.sorted",
        s = c.seq,
        b = c.begin,
    )
}
