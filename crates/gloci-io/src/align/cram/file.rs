//! A CRAM read through its `.crai`: whole chromosomes for per-chromosome runs (D41), and
//! regions for `view` and `--region` (D107, D116), the CRAM counterparts of the BAM index
//! readers.
//!
//! - **The index** lists, per slice, the reference it covers (a slice of several references has
//!   a line per run of reads on one reference). `CramFile` keeps those lines by reference id,
//!   in file order. A region selects the slices whose lines overlap its window (widened by one
//!   base each side, as the BAM reader's), and a chromosome selects every slice it has a line
//!   for.
//! - **Reading** a slice means reading its container from the file (one `read` of its bytes,
//!   CRC-checked like a streamed one) and decoding that slice alone ([`Decoder::decode_where`]).
//!   A reader keeps the container and the slice it decoded last, so the next reference of a
//!   slice shared by several (slices of several references) does not decode it again.
//! - **Records** of other references than the one being read are dropped, uncounted. The ones
//!   kept go through the BAM region parser ([`BamRegionParse`]), so spans, filters, skips and
//!   one-time warnings are the BAM's (D101, D102).
//! - **Reference memory** is one [`RefCache`] for every reader of the file, with a slot per
//!   worker: workers on different chromosomes keep their reference each, and a chromosome is
//!   read from the FASTA once while a worker is on it (D109).
//! - **Errors** name "record N of CHROM": N is counted from the chromosome's first record, as
//!   the streaming reader counts it, so every `-t` reports one place. A `.crai` cannot count
//!   records before a chromosome, so a CRAM never adds "(record M of the file)", streamed or not.

use super::decode::Decoder;
use super::frame::{CramFrames, RawContainer, read_container_header, slice_heads};
use super::reference::{Counters, RefCache};
use super::{CramInfo, ExtentDecoder, check_eof_container, parse_header, reference_of};
use crate::align::region::BamRegionParse;
use crate::align::{AlnHeader, AlnSource};
use crate::detect::{Compression, Format};
use crate::diag::{Diagnostics, InvalidPolicy};
use crate::error::IoError;
use crate::index::{CraiEntry, CramIndex, IndexInfo};
use crate::record::{
    Header, ReaderOptions, RecordBatch, RecordMeta, RecordReader, warn_undeclared,
};
use gloci_core::{ChromTable, Interval, RegionSet, thousands};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// What every reader made from one [`CramFile`] shares.
struct Shared {
    name: String,
    path: PathBuf,
    /// The `.crai` the file is read through.
    index: PathBuf,
    /// The header as records are named after `--chrom-map`.
    aln: Arc<AlnHeader>,
    info: Arc<CramInfo>,
    /// The `.crai` lines of each reference id, in file order.
    by_ref: Vec<Vec<CraiEntry>>,
    /// The start of every data container with lines, ascending: a container's number is its
    /// rank here (as the streaming reader counts them, unless some hold no reads).
    containers: Vec<u64>,
    cache: Arc<RefCache>,
    counters: Arc<Counters>,
    opts: crate::align::AlignOptions,
    policy: InvalidPolicy,
}

/// A CRAM file whose `.crai` is loaded once, from which chromosome and region readers are made
/// cheaply (the CRAM counterpart of an indexed BAM, see `IndexedFile`).
pub(crate) struct CramFile {
    shared: Arc<Shared>,
    /// The first slice (container, offset in it) of each reference, which fixes their order in
    /// the file.
    first: Vec<Option<(u64, u64)>>,
    /// How much of the file each reference's slices take, plus one; 0 for none.
    weights: Vec<u64>,
    /// Reference ids by name as records are named after `--chrom-map`.
    by_name: HashMap<Vec<u8>, usize>,
    /// Containers with reads that have no reference.
    unplaced: Vec<u64>,
    defer_warnings: bool,
}

/// One slice to read: its container, the slice's place in the container's data, and its size.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct Pick {
    container: u64,
    slice: u64,
}

/// A reference to read, and the slices that hold what is wanted of it.
struct Plan {
    rid: usize,
    set: Option<RegionSet>,
    picks: Vec<Pick>,
}

fn stale_or_bad(name: &str, index: &Path, msg: String) -> IoError {
    IoError::Detect {
        path: index.display().to_string(),
        msg: format!(
            "{}: {msg}; rebuild it with: gloci index {name}",
            index.display()
        ),
    }
}

impl CramFile {
    /// Open the CRAM `path` with its `.crai` `index`: the header and reference are checked as
    /// for a streamed read, the index is loaded and checked against the header. `aln` is the
    /// header as records are named after `--chrom-map`.
    pub(crate) fn open(
        path: &Path,
        aln: Arc<AlnHeader>,
        index: &IndexInfo,
        o: &ReaderOptions,
    ) -> Result<CramFile, IoError> {
        let name = path.display().to_string();
        let file = File::open(path).map_err(|source| IoError::Open {
            path: name.clone(),
            source,
        })?;
        let mut reader = BufReader::new(file);
        let head = super::read_head(&mut reader, &name)?;
        check_eof_container(&name)?;
        let raw = parse_header(&name, &head)?;
        let info = Arc::new(CramInfo::new(&name, raw, reference_of(o))?);
        let n_refs = info.aln.refs.len();
        if aln.refs.len() != n_refs {
            return Err(IoError::Internal(format!(
                "{name}: the header given has {} references, the CRAM's has {n_refs}",
                aln.refs.len()
            )));
        }
        let crai = CramIndex::read(&index.path)?;
        let mut entries = crai.entries;
        entries.sort_by_key(|e| (e.container_offset, e.slice_offset, e.start));
        let mut by_ref: Vec<Vec<CraiEntry>> = vec![Vec::new(); n_refs];
        let mut containers: Vec<u64> = Vec::new();
        let mut unplaced: Vec<u64> = Vec::new();
        for e in entries {
            if containers.last() != Some(&e.container_offset) {
                containers.push(e.container_offset);
            }
            match usize::try_from(e.ref_id) {
                Ok(rid) if rid < n_refs => by_ref[rid].push(e),
                Ok(rid) => {
                    return Err(stale_or_bad(
                        &name,
                        &index.path,
                        format!(
                            "a line names reference id {rid}, but {name} has {n_refs} references"
                        ),
                    ));
                }
                Err(_) if e.ref_id == -1 => {
                    if unplaced.last() != Some(&e.container_offset) {
                        unplaced.push(e.container_offset);
                    }
                }
                Err(_) => {
                    return Err(stale_or_bad(
                        &name,
                        &index.path,
                        format!("a line names reference id {}", e.ref_id),
                    ));
                }
            }
        }
        let first = by_ref
            .iter()
            .map(|v| v.first().map(|e| (e.container_offset, e.slice_offset)))
            .collect();
        let weights = by_ref
            .iter()
            .map(|v| {
                if v.is_empty() {
                    return 0;
                }
                let mut seen: Option<(u64, u64)> = None;
                let mut total = 1u64;
                for e in v {
                    let key = (e.container_offset, e.slice_offset);
                    if seen != Some(key) {
                        total = total.saturating_add(e.slice_size);
                        seen = Some(key);
                    }
                }
                total
            })
            .collect();
        let by_name = aln
            .refs
            .iter()
            .enumerate()
            .map(|(rid, (n, _))| (n.clone(), rid))
            .collect();
        let counters = Arc::new(Counters::default());
        // A slot per worker: workers on different chromosomes keep their reference each.
        let cache = Arc::new(RefCache::with_slots(
            info.clone(),
            counters.clone(),
            o.threads.get(),
        ));
        Ok(CramFile {
            shared: Arc::new(Shared {
                name,
                path: path.to_path_buf(),
                index: index.path.clone(),
                aln,
                info,
                by_ref,
                containers,
                cache,
                counters,
                opts: o.accept.alignments.clone().unwrap_or_default(),
                policy: o.invalid,
            }),
            first,
            weights,
            by_name,
            unplaced,
            defer_warnings: false,
        })
    }

    /// Have the readers made from here keep the one-time warnings in their diagnostics, for
    /// the caller to print once for all chunks of a run.
    pub(crate) fn defer_warnings(&mut self) {
        self.defer_warnings = true;
    }

    /// What the readers made from here read from the FASTA, and the most they held at once.
    pub(crate) fn reference_counts(&self) -> super::ReferenceCounts {
        self.shared.counters.snapshot()
    }

    /// The reads with no place, which no chromosome's reader reaches: the unmapped reads at
    /// the end of a coordinate-sorted file. Counted from the slices that hold them.
    pub(crate) fn unplaced_reads(&self) -> u64 {
        let count = || -> Result<u64, IoError> {
            let mut total = 0u64;
            let mut dec: Option<ExtentDecoder> = None;
            for &offset in &self.unplaced {
                let c = read_container(&self.shared, offset)?;
                let heads = slice_heads(&c, &self.shared.name, self.shared.info.ref_count())?;
                let extents = if heads.iter().any(|h| h.ref_id == -2) {
                    let d = dec.get_or_insert_with(|| ExtentDecoder::new(self.shared.info.clone()));
                    d.extents(&c, &heads)?
                } else {
                    vec![None; heads.len()]
                };
                for (h, ext) in heads.iter().zip(extents) {
                    match (h.ref_id, ext) {
                        (-1, _) => total += h.records as u64,
                        (-2, Some(v)) => {
                            total += v.iter().filter(|(r, _, _)| *r < 0).count() as u64
                        }
                        _ => {}
                    }
                }
            }
            Ok(total)
        };
        count().unwrap_or_else(|e| {
            log::warn!(
                "cannot count the reads with no place in {}: {e}",
                self.shared.name
            );
            0
        })
    }

    fn rid_of(&self, chrom: &[u8]) -> Option<usize> {
        self.by_name.get(chrom).copied()
    }

    /// The references of the header, in id order.
    pub(crate) fn reference_names(&self) -> Vec<Vec<u8>> {
        self.shared
            .aln
            .refs
            .iter()
            .map(|(n, _)| n.clone())
            .collect()
    }

    /// The references that hold reads, in the order the file shows them.
    pub(crate) fn file_order_reference_names(&self) -> Vec<Vec<u8>> {
        let mut named: Vec<((u64, u64), usize)> = self
            .first
            .iter()
            .enumerate()
            .filter_map(|(rid, f)| f.map(|f| (f, rid)))
            .collect();
        named.sort();
        named
            .into_iter()
            .map(|(_, rid)| self.shared.aln.refs[rid].0.clone())
            .collect()
    }

    pub(crate) fn non_empty_reference_names(&self) -> Vec<Vec<u8>> {
        self.weights
            .iter()
            .enumerate()
            .filter(|(_, w)| **w > 0)
            .map(|(rid, _)| self.shared.aln.refs[rid].0.clone())
            .collect()
    }

    /// How much of the file `chrom`'s slices take (plus one), or 0 when it has none.
    pub(crate) fn weight(&self, chrom: &[u8]) -> u64 {
        self.rid_of(chrom).map_or(0, |rid| self.weights[rid])
    }

    /// The slices of `rid` whose lines overlap a window of `set` (all of them without one),
    /// each once, in file order.
    fn picks(&self, rid: usize, set: Option<&RegionSet>) -> Vec<Pick> {
        // The regions as windows, 0-based half-open, merged. The `.crai` gives a slice the
        // span of its reads clipped at the reference's end, and every read covers a base, so a
        // slice is wanted exactly when its span overlaps a region, as samtools finds them (a
        // region that starts past the end of the reference finds nothing there, though reads
        // run past it). An empty region (a point, with `inclusive_points`) takes the base
        // before it too.
        let windows: Option<Vec<(u64, u64)>> = set.map(|s| {
            let mut w: Vec<(u64, u64)> = Vec::new();
            for r in s.regions() {
                let (lo, hi) = if r.start >= r.end {
                    (r.start.saturating_sub(1), r.start.saturating_add(1))
                } else {
                    (r.start, r.end)
                };
                match w.last_mut() {
                    Some(l) if lo <= l.1 => l.1 = l.1.max(hi),
                    _ => w.push((lo, hi)),
                }
            }
            w
        });
        let mut picks: Vec<Pick> = Vec::new();
        for e in &self.shared.by_ref[rid] {
            if let Some(w) = &windows {
                // The line's bases, 0-based half-open.
                let lo = u64::try_from(e.start.saturating_sub(1)).unwrap_or(0);
                let hi = lo.saturating_add(u64::try_from(e.span.max(1)).unwrap_or(1));
                let at = w.partition_point(|&(_, whi)| whi <= lo);
                if w.get(at).is_none_or(|&(wlo, _)| wlo >= hi) {
                    continue;
                }
            }
            let pick = Pick {
                container: e.container_offset,
                slice: e.slice_offset,
            };
            if picks.last() != Some(&pick) {
                picks.push(pick);
            }
        }
        picks
    }

    fn reader(&self, plans: Vec<Plan>, whole: bool, inclusive: bool, all: bool) -> Reader {
        let s = &self.shared;
        Reader {
            parse: BamRegionParse::new(s.aln.clone(), s.opts.clone(), all, self.defer_warnings),
            header: Header::Alignment(s.aln.clone()),
            shared: s.clone(),
            diag: Diagnostics::default(),
            inclusive,
            whole,
            plans,
            p: 0,
            s: 0,
            decoder: None,
            container: None,
            cur: Vec::new(),
            cur_key: None,
            at: 0,
            in_slice: 0,
            emitted: 0,
            seen: 0,
        }
    }

    /// The records that overlap `regions` (named as records are after `--chrom-map`), each
    /// once, in file order within a reference.
    pub(crate) fn regions_reader(
        &self,
        regions: &[(Vec<u8>, Interval)],
        inclusive_points: bool,
        all: bool,
    ) -> Reader {
        log::debug!(
            "reading {} through its index for {} region(s)",
            self.shared.name,
            regions.len()
        );
        let mut by_rid: std::collections::BTreeMap<usize, Vec<Interval>> = Default::default();
        for (chrom, iv) in regions {
            if let Some(rid) = self.rid_of(chrom) {
                by_rid.entry(rid).or_default().push(*iv);
            }
        }
        let mut plans: Vec<Plan> = by_rid
            .into_iter()
            .map(|(rid, ivs)| {
                let set = RegionSet::new(ivs);
                let picks = self.picks(rid, Some(&set));
                Plan {
                    rid,
                    set: Some(set),
                    picks,
                }
            })
            .collect();
        // File order, as the BAM reader visits them.
        plans.sort_by_key(|p| self.first[p.rid].unwrap_or((u64::MAX, u64::MAX)));
        self.reader(plans, false, inclusive_points, all)
    }

    /// The records of every slice of `chroms`, chromosome by chromosome in the order given.
    /// There is no footprint test.
    pub(crate) fn chromosomes_reader(&self, chroms: &[Vec<u8>]) -> Reader {
        let plans = chroms
            .iter()
            .filter_map(|c| self.rid_of(c))
            .map(|rid| Plan {
                rid,
                set: None,
                picks: self.picks(rid, None),
            })
            .collect();
        self.reader(plans, true, true, false)
    }

    /// The records that overlap `regions` as every-record [`AlnSource`] reads them (D116).
    pub(crate) fn aln_regions_source(
        &self,
        regions: &[(Vec<u8>, Interval)],
        inclusive_points: bool,
    ) -> Box<dyn AlnSource> {
        Box::new(RegionSource {
            aln: self.shared.aln.clone(),
            reader: self.regions_reader(regions, inclusive_points, true),
            chroms: ChromTable::first_appearance(),
            batch: RecordBatch::default(),
        })
    }
}

/// The container that starts at byte `offset`, read from the file as the streaming reader
/// reads it.
fn read_container(s: &Shared, offset: u64) -> Result<RawContainer, IoError> {
    let open_err = |source| IoError::Open {
        path: s.name.clone(),
        source,
    };
    let mut f = File::open(&s.path).map_err(open_err)?;
    f.seek(SeekFrom::Start(offset))
        .map_err(|source| IoError::Read {
            path: s.name.clone(),
            source,
        })?;
    let number = match s.containers.binary_search(&offset) {
        Ok(i) => i as u64 + 1,
        Err(i) => i as u64 + 1,
    };
    // The index says a container starts here: if no container header does, the index is at
    // fault (stale, or made for another file), not the CRAM.
    let mut probe = Vec::new();
    let mut reader = BufReader::new(f);
    let ok = read_container_header(&mut reader, &mut probe).is_ok();
    if !ok {
        return Err(IoError::Detect {
            path: s.index.display().to_string(),
            msg: format!(
                "{}: the index names a container of {} at byte {}, where the file has none; \
the index looks out of date or belongs to another file; rebuild it with: gloci index {}",
                s.index.display(),
                s.name,
                thousands(offset),
                s.name
            ),
        });
    }
    reader
        .seek(SeekFrom::Start(offset))
        .map_err(|source| IoError::Read {
            path: s.name.clone(),
            source,
        })?;
    let mut frames = CramFrames::at(&s.name, Box::new(reader), offset, number);
    match frames.next()? {
        Some(c) => Ok(c),
        None => Err(IoError::Detect {
            path: s.name.clone(),
            msg: format!(
                "{}: the index names a container at byte {}, where the file has its end; \
rebuild the index with: gloci index {}",
                s.name,
                thousands(offset),
                s.name
            ),
        }),
    }
}

/// Reads the records of some references of a CRAM (see the module notes).
pub(crate) struct Reader {
    shared: Arc<Shared>,
    header: Header,
    diag: Diagnostics,
    parse: BamRegionParse,
    inclusive: bool,
    /// Every record of each plan's reference is emitted, with no footprint test.
    whole: bool,
    plans: Vec<Plan>,
    /// The plan being read and the next of its slices.
    p: usize,
    s: usize,
    decoder: Option<Decoder>,
    /// The container last read, and the slice decoded last into `cur`.
    container: Option<RawContainer>,
    cur: Vec<u8>,
    cur_key: Option<Pick>,
    at: usize,
    /// Records of the current slice passed so far (for places of a region read).
    in_slice: u64,
    /// Records emitted so far: each one's `loc`, which keeps file order for ties.
    emitted: u64,
    /// Records of the current reference read so far, the number errors name.
    seen: u64,
}

impl Reader {
    /// An invalid record: an error, or counted under `--skip-invalid`.
    fn bad(&mut self, msg: &str) -> Result<(), IoError> {
        let place = match self.plans.get(self.p) {
            Some(plan) if self.whole => format!(
                "record {} of {}",
                thousands(self.seen),
                String::from_utf8_lossy(&self.shared.aln.refs[plan.rid].0)
            ),
            _ => {
                let c = self.cur_key.map_or(0, |k| k.container);
                let n = self
                    .shared
                    .containers
                    .binary_search(&c)
                    .map_or(0, |i| i as u64 + 1);
                format!(
                    "record {} of a slice of container {} (at byte {})",
                    thousands(self.in_slice),
                    thousands(n),
                    thousands(c)
                )
            }
        };
        match self.shared.policy {
            InvalidPolicy::Fail => Err(IoError::Detect {
                path: self.shared.name.clone(),
                msg: format!("{} {place}: {msg}", self.shared.name),
            }),
            InvalidPolicy::Skip => {
                self.diag.record_place(&place, msg);
                Ok(())
            }
        }
    }

    /// Decode the slice `pick` into `self.cur`, reading its container unless it is the one
    /// held.
    fn load(&mut self, pick: Pick) -> Result<(), IoError> {
        self.at = 0;
        self.in_slice = 0;
        if self.cur_key == Some(pick) {
            return Ok(());
        }
        let s = &self.shared;
        if self.container.as_ref().map(|c| c.offset) != Some(pick.container) {
            self.container = Some(read_container(s, pick.container)?);
        }
        let Some(c) = self.container.as_ref() else {
            return Err(IoError::Internal("a container was just read".into()));
        };
        if self.decoder.is_none() {
            self.decoder = Some(Decoder::new(s.info.clone(), s.cache.clone())?.pinning());
        }
        let Some(decoder) = self.decoder.as_mut() else {
            return Err(IoError::Internal("a decoder was just made".into()));
        };
        self.cur.clear();
        self.cur_key = None;
        let n = decoder.decode_where(c, &|h| h.landmark as u64 == pick.slice, &mut self.cur)?;
        if n == 0 {
            return Err(IoError::Detect {
                path: s.name.clone(),
                msg: format!(
                    "{}: its index names a slice at offset {} of {}, which holds no such slice; \
the index is out of date or damaged, rebuild it with: gloci index {}",
                    s.name,
                    thousands(pick.slice),
                    c.place(),
                    s.name
                ),
            });
        }
        self.cur_key = Some(pick);
        Ok(())
    }

    /// Append the next matching record to `batch`; false when done.
    fn next(&mut self, chroms: &mut ChromTable, batch: &mut RecordBatch) -> Result<bool, IoError> {
        loop {
            if self.at >= self.cur.len() {
                let Some(plan) = self.plans.get(self.p) else {
                    return Ok(false);
                };
                let Some(&pick) = plan.picks.get(self.s) else {
                    // A new reference: its records are numbered from its first.
                    self.p += 1;
                    self.s = 0;
                    self.seen = 0;
                    continue;
                };
                self.s += 1;
                self.load(pick)?;
                continue;
            }
            let Some(plan) = self.plans.get(self.p) else {
                return Ok(false);
            };
            let rid = plan.rid;
            // The decoder wrote a whole frame: its `block_size`, then the record.
            let frame = &self.cur[self.at..];
            let size = frame
                .get(..4)
                .map(|b| u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as usize);
            let Some(size) = size.filter(|&n| frame.len() >= 4 + n) else {
                return Err(IoError::Internal(format!(
                    "{}: a decoded slice ends inside a record",
                    self.shared.name
                )));
            };
            let (from, to) = (self.at + 4, self.at + 4 + size);
            self.at = to;
            self.in_slice += 1;
            let raw = &self.cur[from..to];
            // (A record too short to tell is let through, so parsing reports it.)
            let on_ref = raw.len() < 4
                || usize::try_from(i32::from_le_bytes([raw[0], raw[1], raw[2], raw[3]])).ok()
                    == Some(rid);
            if !on_ref {
                continue;
            }
            self.seen += 1;
            let off = batch.buf.len();
            batch.buf.extend_from_slice(raw);
            let parsed = self.parse.parse(
                &self.shared.name,
                &mut batch.buf,
                off,
                chroms,
                &mut self.diag,
            );
            let (chrom, key, span, strand) = match parsed {
                Ok(Some(p)) => p,
                Ok(None) => {
                    batch.buf.truncate(off);
                    continue;
                }
                Err(msg) => {
                    batch.buf.truncate(off);
                    self.bad(&msg)?;
                    continue;
                }
            };
            if chroms.has_undeclared() {
                warn_undeclared(chroms, &self.shared.name, None);
            }
            let matched = self.whole
                || self.plans[self.p]
                    .set
                    .as_ref()
                    .is_some_and(|set| set.matches(span, self.inclusive));
            if !matched {
                batch.buf.truncate(off);
                continue;
            }
            self.emitted += 1;
            let len = (batch.buf.len() - off) as u32;
            batch.recs.push(RecordMeta {
                chrom,
                key,
                span,
                strand,
                loc: self.emitted,
                off,
                len,
            });
            return Ok(true);
        }
    }
}

impl RecordReader for Reader {
    fn name(&self) -> &str {
        &self.shared.name
    }
    fn format(&self) -> Format {
        Format::Bam
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
    fn whole_records(&self) -> bool {
        self.shared.opts.full_records
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

/// Index-backed [`AlnSource`]: the BAM records of some regions of a CRAM.
struct RegionSource {
    aln: Arc<AlnHeader>,
    reader: Reader,
    chroms: ChromTable,
    batch: RecordBatch,
}

impl AlnSource for RegionSource {
    fn name(&self) -> &str {
        &self.reader.shared.name
    }
    fn layout(&self) -> Format {
        Format::Bam
    }
    fn format(&self) -> Format {
        Format::Cram
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
