//! Decoding CRAM containers to BAM records: the slice checks before noodles decodes, the reads'
//! extents against their reference bases (padded with `N` where they run past), the
//! transcoding with htslib's derived fields, and the position-only decoding the `.crai` needs.

use super::super::encode::reg2bin;
use super::frame::{
    CompressionInfo, RawContainer, SliceHead, check_compression_block, check_slice_blocks,
    compression, detect_err, itf8_values, slice_heads, unsupported_codec,
};
use super::reference::{Counters, CramInfo, RefAdapter, RefCache, RefState};
use crate::error::IoError;
use gloci_core::thousands;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam::alignment::Record as _;
use noodles_sam::alignment::io::Write as _;
use std::borrow::Cow;
use std::io;
use std::sync::Arc;

/// Records a slice may claim beyond one per byte of its blocks: noodles allocates a record for
/// each before reading any, so a corrupt count must not ask for gigabytes.
const SLACK_RECORDS: usize = 4096;
/// The least overhang past a reference's end gloci decodes against `N`.
const MIN_OVERHANG: i64 = 1 << 20;

// ---- decoding ----

/// Where a slice's reads find their reference bases.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Bases {
    /// A repository over `--reference`.
    Repository,
    /// The slice's embedded reference block.
    Embedded,
    /// Every base is stored, but noodles slices a placeholder sequence at each mapped read (a
    /// slice of several references without `RR`).
    Placeholder,
    /// Nowhere: every base is stored.
    Stored,
}

/// What the reads of a decoded slice need beyond what they were decoded with.
pub(crate) enum Extent {
    Fits,
    /// These references need sequences this long (`N` after their end).
    Pad(Vec<(usize, usize)>),
    /// The embedded bases need to run to this reference position (1-based).
    EmbeddedTo(i64),
    /// The placeholders of these references need this length.
    Placeholder(Vec<(usize, usize)>),
}

/// While container `c` of `name` is decoded, a panic in noodles-cram means malformed data the
/// checks before it let through (a crafted file): say so, and where (`#cram-limits`).
fn decoder_note(name: &str, c: &RawContainer) -> crate::atomic::PanicNote {
    crate::atomic::panic_note(format!(
        "{name}: the CRAM decoder failed on malformed data in {}; the file is corrupt or was \
crafted (see #cram-limits in docs/differences.md)",
        c.place()
    ))
}

/// The references of the mapped reads of `records`.
fn mapped_refs(records: &[cram::Record<'_>], info: &CramInfo) -> Vec<usize> {
    let mut out = Vec::new();
    for r in records {
        if r.flags().is_ok_and(|f| !f.is_unmapped())
            && let Some(Ok(id)) = r.reference_sequence_id(&info.sam)
            && id < info.aln.refs.len()
            && !out.contains(&id)
        {
            out.push(id);
        }
    }
    out
}

/// The block with content id `id` among a slice's external blocks.
fn block<'a>(ext: &'a [(i32, Cow<'_, [u8]>)], id: i32) -> Option<&'a [u8]> {
    ext.iter().find(|(i, _)| *i == id).map(|(_, b)| &b[..])
}

/// The values of data series `key` when it has an external block of its own, or one value for
/// every one of the slice's `n` records when it is coded as a constant.
fn series(
    comp: &CompressionInfo,
    ext: &[(i32, Cow<'_, [u8]>)],
    key: &[u8; 2],
    n: usize,
) -> Option<Vec<i32>> {
    if let Some(&v) = comp.constant.get(key) {
        return Some(vec![v; n]);
    }
    itf8_values(block(ext, comp.own_block(key)?)?)
}

/// Check the decoded blocks of slice `head` before noodles reads its records: it allocates a
/// record per claimed record and reserves each read's claimed features, and it panics on a
/// mapped read with no reference in a slice of several references and on a mate past the
/// slice's end. Where a data series has an external block of its own its values are checked;
/// elsewhere noodles is trusted (crafted files only: `#cram-crafted-input`). Returns, for a
/// slice of several references whose reads' references are known, those of its mapped reads.
fn precheck(
    head: &SliceHead,
    n_refs: usize,
    comp: &CompressionInfo,
    core: &[u8],
    ext: &[(i32, Cow<'_, [u8]>)],
) -> Result<Option<Vec<usize>>, String> {
    let bytes = core.len() + ext.iter().map(|(_, b)| b.len()).sum::<usize>();
    if head.records > SLACK_RECORDS + bytes {
        return Err(format!(
            "a slice claims {} records in {} bytes",
            thousands(head.records as u64),
            thousands(bytes as u64)
        ));
    }
    // Each feature reads one feature code: their count cannot pass the codes stored. (A constant
    // count is checked as one value: no read may need more.)
    let fns = match comp.constant.get(b"FN") {
        Some(&v) => Some(vec![v]),
        None => series(comp, ext, b"FN", head.records),
    };
    if let Some(fns) = fns {
        let codes = comp
            .external
            .get(b"FC")
            .and_then(|&id| block(ext, id))
            .map_or(bytes, <[u8]>::len);
        let total: i64 = fns.iter().map(|&n| i64::from(n)).sum();
        if fns.iter().any(|&n| n < 0) || total > codes as i64 {
            return Err(format!(
                "its reads claim {} features, and the slice holds {} feature codes",
                thousands(total.max(0) as u64),
                thousands(codes as u64)
            ));
        }
    }
    let n = head.records;
    // A series read once per record, in a block of its own, holds one value per record: a
    // count that disagrees is a corrupt slice (the checks below rely on it).
    let per_record = |key: &[u8; 2]| -> Result<Option<Vec<i32>>, String> {
        match series(comp, ext, key, n) {
            Some(v) if v.len() != n => Err(format!(
                "a slice claims {} records, and its {} series holds {} values",
                thousands(n as u64),
                String::from_utf8_lossy(key),
                thousands(v.len() as u64)
            )),
            v => Ok(v),
        }
    };
    let flags = per_record(b"BF")?;
    let cram_flags = per_record(b"CF")?;
    // A mate given as a distance must be a later record of the slice.
    if let (Some(cf), Some(nf)) = (&cram_flags, series(comp, ext, b"NF", n)) {
        let constant = comp.constant.get(b"NF").copied();
        let mut next = nf.iter();
        for (i, &f) in cf.iter().enumerate() {
            if f & 0x2 == 0 && f & 0x4 != 0 {
                let Some(d) = constant.or_else(|| next.next().copied()) else {
                    break;
                };
                if d < 0 || i as i64 + i64::from(d) + 1 >= head.records as i64 {
                    return Err(format!(
                        "read {} of a slice of {} names its mate {d} reads further on",
                        thousands(i as u64 + 1),
                        thousands(head.records as u64)
                    ));
                }
            }
        }
    }
    if head.ref_id != -2 {
        return Ok(None);
    }
    let (Some(bf), Some(ri)) = (flags, per_record(b"RI")?) else {
        return Ok(None);
    };
    let mut needed = Vec::new();
    for (i, (&f, &r)) in bf.iter().zip(&ri).enumerate() {
        if f & 0x4 != 0 {
            continue;
        }
        if r < 0 {
            return Err(format!(
                "mapped read {} of a slice of several references has no reference",
                thousands(i as u64 + 1)
            ));
        }
        if (r as usize) < n_refs && !needed.contains(&(r as usize)) {
            needed.push(r as usize);
        }
    }
    Ok(Some(needed))
}

/// Decodes containers to BAM records. Each worker of the pool has one; the reference memory is
/// the cache they share.
pub(crate) struct Decoder {
    pub(crate) info: Arc<CramInfo>,
    /// The reference memory every decoder shares (slots and generations).
    pub(crate) cache: Arc<RefCache>,
    /// A repository of empty placeholders: for slices that read none (noodles still takes
    /// one), and to find, by decoding the records once without their bases, which references
    /// a slice's reads need when its RI series cannot tell.
    pub(crate) probe: fasta::Repository,
    pub(crate) probe_state: Arc<RefState>,
    pub(crate) bam: noodles_bam::io::Writer<Vec<u8>>,
}

/// Where a slice's reads find the repository they need.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Route {
    /// A slot of the shared cache: one reference.
    Slot(usize),
    /// A shared generation of reference bases (`true`) or placeholders.
    Generation(bool),
    /// None: the slice reads no repository.
    Idle,
}

impl Decoder {
    pub(crate) fn new(info: Arc<CramInfo>, cache: Arc<RefCache>) -> Result<Decoder, IoError> {
        // Placeholders of no length (`RefState::default()`): the probe never builds bases.
        let probe_state = Arc::new(RefState::default());
        let adapter = RefAdapter {
            info: info.clone(),
            fasta: None,
            state: probe_state.clone(),
            counters: cache.counters.clone(),
        };
        Ok(Decoder {
            info,
            cache,
            probe: fasta::Repository::new(adapter),
            probe_state,
            bam: noodles_bam::io::Writer::from(Vec::new()),
        })
    }

    pub(crate) fn corrupt(&self, c: &RawContainer, e: &dyn std::fmt::Display) -> IoError {
        let name = &self.info.name;
        detect_err(
            name,
            format!(
                "{name}: {} cannot be decoded: {e}; the file is corrupt",
                c.place()
            ),
        )
    }

    /// The gloci error for a noodles failure while decoding `head` of `c` with a repository
    /// whose adapter shares `state`.
    pub(crate) fn decode_error(
        &self,
        c: &RawContainer,
        head: &SliceHead,
        state: &RefState,
        e: io::Error,
    ) -> IoError {
        if let Some(f) = state.failure.lock().ok().and_then(|mut s| s.take()) {
            return f;
        }
        if e.kind() == io::ErrorKind::InvalidData
            && e.to_string()
                .starts_with("reference sequence checksum mismatch")
            && head.ref_id >= 0
        {
            return self
                .info
                .mismatch(head.ref_id as usize, i64::from(head.start), head.end());
        }
        self.corrupt(c, &e)
    }

    /// Where the reads of `head` find their bases, after checking that they can (D109).
    pub(crate) fn bases(
        &self,
        c: &RawContainer,
        head: &SliceHead,
        required: bool,
    ) -> Result<Bases, IoError> {
        let bases = match head.ref_id {
            -1 => Bases::Stored,
            _ if required => Bases::Repository,
            -2 => Bases::Placeholder,
            _ if head.embedded != -1 => Bases::Embedded,
            _ => Bases::Stored,
        };
        if head.ref_id >= 0 && required {
            let rid = head.ref_id as usize;
            let (name, len) = &self.info.aln.refs[rid];
            if !self.info.in_fasta[rid] {
                return Err(self.info.needs_reference(name));
            }
            if head.end() > *len as i64 {
                return Err(self.corrupt(
                    c,
                    &format!(
                        "a slice of {} runs to {}, past its end at {}",
                        self.info.chrom(rid),
                        thousands(head.end() as u64),
                        thousands(*len)
                    ),
                ));
            }
        }
        Ok(bases)
    }

    /// Append the BAM records (each with its `block_size`) of container `c` to `out`.
    pub fn decode(&mut self, c: &RawContainer, out: &mut Vec<u8>) -> Result<(), IoError> {
        let _note = decoder_note(&self.info.name, c);
        let heads = slice_heads(c, &self.info.name, self.info.aln.refs.len())?;
        let comp = compression(c);
        if let Some(what) = &comp.unsupported {
            return Err(unsupported_codec(&self.info.name, c, what));
        }
        let required = comp.required;
        check_compression_block(c).map_err(|m| self.corrupt(c, &m))?;
        let mut container = cram::io::reader::Container::default();
        cram::io::Reader::new(&c.bytes[..])
            .read_container(&mut container)
            .map_err(|e| self.corrupt(c, &e))?;
        let ch = container
            .compression_header()
            .map_err(|e| self.corrupt(c, &e))?;
        // The records borrow the header and the leases the cache; the decoder itself stays free
        // to transcode them.
        let info = self.info.clone();
        let cache = self.cache.clone();
        for (head, slice) in heads.iter().zip(container.slices()) {
            let slice = slice.map_err(|e| self.corrupt(c, &e))?;
            let bases = self.bases(c, head, required)?;
            check_slice_blocks(c, head, comp.external.get(b"RN").copied())
                .map_err(|m| self.corrupt(c, &m))?;
            let (core, ext) = slice.decode_blocks().map_err(|e| self.corrupt(c, &e))?;
            if bases == Bases::Embedded && block(&ext, head.embedded).is_none() {
                // noodles finds the embedded block by its id, and `expect`s it.
                return Err(self.corrupt(
                    c,
                    &"a slice names an embedded reference block it does not have",
                ));
            }
            let needed = precheck(head, info.aln.refs.len(), &comp, &core, &ext)
                .map_err(|m| self.corrupt(c, &m))?;
            // The references the slice's reads need, when known before decoding.
            let mut refs: Option<Vec<usize>> = match head.ref_id {
                r if r >= 0 => Some(vec![r as usize]),
                -2 => needed,
                _ => Some(Vec::new()),
            };
            if refs.is_none() && bases != Bases::Stored && bases != Bases::Embedded {
                // Which references the reads need, from the records decoded without bases.
                self.probe.clear();
                if let Ok(mut l) = self.probe_state.loaded.lock() {
                    l.clear();
                }
                let records = slice
                    .records(self.probe.clone(), &info.sam, &ch, &core, &ext)
                    .map_err(|e| self.decode_error(c, head, &self.probe_state, e))?;
                refs = Some(mapped_refs(&records, &info));
            }
            // The length each needs (beyond its own, for padding or a placeholder's reach).
            let mut lengths: Vec<(usize, usize)> = Vec::new();
            let mut route = match (bases, &refs) {
                (Bases::Repository, Some(r)) if r.len() == 1 => Route::Slot(r[0]),
                (Bases::Repository, _) => Route::Generation(true),
                (Bases::Placeholder, _) => Route::Generation(false),
                _ => Route::Idle,
            };
            // Reads past the end of their reference bases are decoded again against `N`.
            let mut padded: Option<Vec<(i32, Cow<'_, [u8]>)>> = None;
            let mut done = false;
            for _ in 0..6 {
                let need: Option<Vec<(usize, usize)>> = refs.as_ref().map(|rs| {
                    rs.iter()
                        .map(|&r| {
                            let own = usize::try_from(info.aln.refs[r].1).unwrap_or(0);
                            let len = lengths
                                .iter()
                                .find(|(x, _)| *x == r)
                                .map_or(own, |&(_, l)| l.max(own));
                            (r, len)
                        })
                        .collect()
                });
                let lease = match route {
                    Route::Slot(rid) => Some(cache.acquire(rid, c.number)?),
                    Route::Generation(rr) => Some(cache.lease(rr, need.as_deref())?),
                    Route::Idle => None,
                };
                let (repo, state) = match &lease {
                    Some(l) => (l.repo.clone(), l.state.clone()),
                    None => (self.probe.clone(), self.probe_state.clone()),
                };
                let srcs = padded.as_deref().unwrap_or(&ext);
                let records = slice
                    .records(repo, &info.sam, &ch, &core, srcs)
                    .map_err(|e| self.decode_error(c, head, &state, e))?;
                cache.note_live();
                if refs.is_none() {
                    refs = Some(mapped_refs(&records, &info));
                }
                let available = |rid: usize| {
                    let (n, len) = &info.aln.refs[rid];
                    state.length_of(n, usize::try_from(*len).unwrap_or(0)) as i64
                };
                match self.extent(c, head, bases, &available, &records, srcs)? {
                    Extent::Fits => {
                        self.transcode(c, &records, out)?;
                        done = true;
                        break;
                    }
                    Extent::Pad(more) | Extent::Placeholder(more) => {
                        drop(records);
                        drop(lease);
                        for (r, len) in more {
                            match lengths.iter_mut().find(|(x, _)| *x == r) {
                                Some((_, l)) => *l = (*l).max(len),
                                None => lengths.push((r, len)),
                            }
                            if let Some(rs) = refs.as_mut()
                                && !rs.contains(&r)
                            {
                                rs.push(r);
                            }
                        }
                        if let Route::Slot(_) = route {
                            route = Route::Generation(true);
                        }
                    }
                    Extent::EmbeddedTo(end) => {
                        drop(records);
                        let want = usize::try_from(end - i64::from(head.start) + 1).unwrap_or(0);
                        let mut v = ext.to_vec();
                        for (id, src) in &mut v {
                            if *id == head.embedded && src.len() < want {
                                let mut b = src.to_vec();
                                b.resize(want, b'N');
                                *src = Cow::Owned(b);
                            }
                        }
                        padded = Some(v);
                    }
                }
            }
            if !done {
                return Err(IoError::Internal(format!(
                    "{}: the reads of {} still run past their reference bases after padding",
                    self.info.name,
                    c.place()
                )));
            }
        }
        Ok(())
    }

    /// Whether every mapped read of `records` lies within the reference bases it decodes with
    /// (`available(rid)` long for reference `rid`). A read may run past its reference's end, and
    /// is decoded against `N` there, by at most 1 Mb or four times the slice's longest read;
    /// further is gloci's limit, an error.
    pub(crate) fn extent(
        &self,
        c: &RawContainer,
        head: &SliceHead,
        bases: Bases,
        available: &dyn Fn(usize) -> i64,
        records: &[cram::Record<'_>],
        srcs: &[(i32, Cow<'_, [u8]>)],
    ) -> Result<Extent, IoError> {
        if bases == Bases::Stored {
            return Ok(Extent::Fits);
        }
        let embedded_len = block(srcs, head.embedded).map_or(0, |b| b.len() as i64);
        // (reference, how far past its end) of the read that runs furthest past one.
        let mut over: Option<(usize, i64)> = None;
        let mut longest = 0i64;
        let mut need_end = 0i64;
        let mut need: Vec<(usize, usize)> = Vec::new();
        let mut want = |rid: usize, len: i64| match need.iter_mut().find(|(r, _)| *r == rid) {
            Some((_, l)) => *l = (*l).max(len as usize),
            None => need.push((rid, len as usize)),
        };
        for r in records {
            let flags = r.flags().map_err(|e| self.corrupt(c, &e))?;
            if flags.is_unmapped() {
                continue;
            }
            let Some(Ok(start)) = r.alignment_start() else {
                return Err(self.corrupt(c, &"a mapped read has no position"));
            };
            let start = usize::from(start) as i64;
            let end = match r.alignment_end() {
                Some(Ok(e)) => usize::from(e) as i64,
                _ => start - 1,
            };
            let rid = match r.reference_sequence_id(&self.info.sam) {
                Some(Ok(id)) if id < self.info.aln.refs.len() => id,
                _ => return Err(self.corrupt(c, &"a mapped read has no reference")),
            };
            let len = self.info.aln.refs[rid].1 as i64;
            if end > len && over.is_none_or(|(_, o)| end - len > o) {
                over = Some((rid, end - len));
            }
            match bases {
                Bases::Repository | Bases::Placeholder => {
                    let reach = start.max(end);
                    if reach > available(rid) {
                        want(rid, reach);
                    }
                }
                Bases::Embedded => {
                    let first = i64::from(head.start);
                    if start < first {
                        return Err(self.corrupt(
                            c,
                            &"a read starts before the reference bases embedded in its slice",
                        ));
                    }
                    if end > first + embedded_len - 1 {
                        need_end = need_end.max(end);
                    }
                }
                Bases::Stored => {}
            }
            // The read's length (its bases, stored or not); the checks above make `sequence()`
            // safe to ask for it.
            longest = longest.max(r.sequence().len() as i64);
        }
        let cap = MIN_OVERHANG.max(4 * longest);
        // gloci's limit, not a fault of the file: samtools decodes such reads.
        let limit = |what: String| {
            let name = &self.info.name;
            detect_err(
                name,
                format!(
                    "{name}: in {}, {what}, more than gloci's limit of {} (1 Mb, or four times \
the slice's longest read; see #cram-limits in docs/differences.md)",
                    c.place(),
                    thousands(cap as u64)
                ),
            )
        };
        if let Some((rid, o)) = over
            && o > cap
        {
            return Err(limit(format!(
                "a read runs {} bases past the end of {} ({} bases)",
                thousands(o as u64),
                self.info.chrom(rid),
                thousands(self.info.aln.refs[rid].1)
            )));
        }
        if bases == Bases::Embedded && need_end > 0 {
            let past = need_end - (i64::from(head.start) + embedded_len - 1);
            if past > cap {
                return Err(limit(format!(
                    "a read runs {} bases past the reference bases embedded in its slice",
                    thousands(past as u64)
                )));
            }
        }
        Ok(if need_end > 0 {
            Extent::EmbeddedTo(need_end)
        } else if need.is_empty() {
            Extent::Fits
        } else if bases == Bases::Placeholder {
            Extent::Placeholder(need)
        } else {
            Extent::Pad(need)
        })
    }

    pub(crate) fn transcode(
        &mut self,
        c: &RawContainer,
        records: &[cram::Record<'_>],
        out: &mut Vec<u8>,
    ) -> Result<(), IoError> {
        for r in records {
            self.bam.get_mut().clear();
            if let Err(e) = self.bam.write_alignment_record(&self.info.sam, r) {
                return Err(self.corrupt(c, &format!("a record cannot be written as BAM: {e}")));
            }
            let rec = self.bam.get_mut();
            htslib_fields(&mut rec[4..]);
            out.extend_from_slice(rec);
        }
        Ok(())
    }
}

/// Set the fields htslib derives when it decodes a CRAM record: MAPQ 0 for an unmapped read
/// (CRAM does not store it) and `bin` from the footprint, as the SAM encoder computes it.
pub(crate) fn htslib_fields(r: &mut [u8]) {
    if r.len() < 32 {
        return;
    }
    let pos = i64::from(i32::from_le_bytes([r[4], r[5], r[6], r[7]]));
    let l_name = usize::from(r[8]);
    let n_cigar = usize::from(u16::from_le_bytes([r[12], r[13]]));
    let flag = u16::from_le_bytes([r[14], r[15]]);
    let unmapped = flag & 0x4 != 0;
    if unmapped && r[9] == 255 {
        r[9] = 0;
    }
    let at = 32 + l_name;
    let rlen: i64 = r
        .get(at..at + 4 * n_cigar)
        .unwrap_or(&[])
        .chunks_exact(4)
        .map(|w| u32::from_le_bytes([w[0], w[1], w[2], w[3]]))
        .filter(|v| matches!(v & 0xf, 0 | 2 | 3 | 7 | 8))
        .map(|v| i64::from(v >> 4))
        .sum();
    let end = if unmapped { pos + 1 } else { pos + rlen.max(1) };
    r[10..12].copy_from_slice(&reg2bin(pos, end).to_le_bytes());
}

// ---- for the index ----

/// `(ref id, start, end)` of a record, as htslib's index builder sees it.
pub(crate) type Extents = Vec<(i32, i64, i64)>;

/// The records of slices that hold several references, for the `.crai`: each record's
/// reference id (-1 for none), start and end, as htslib's index builder sees them. No reference
/// bases are read: the adapter is told the container needs none (so it gives empty
/// placeholders), and the reads' bases are never put together.
pub(crate) struct ExtentDecoder {
    info: Arc<CramInfo>,
    repo: fasta::Repository,
    state: Arc<RefState>,
}

impl ExtentDecoder {
    pub fn new(info: Arc<CramInfo>) -> ExtentDecoder {
        // Placeholders of no length: the reads' bases are never put together.
        let state = Arc::new(RefState::default());
        let adapter = RefAdapter {
            info: info.clone(),
            fasta: None,
            state: state.clone(),
            counters: Arc::new(Counters::default()),
        };
        ExtentDecoder {
            info,
            repo: fasta::Repository::new(adapter),
            state,
        }
    }

    pub fn info(&self) -> &CramInfo {
        &self.info
    }

    /// For each slice of `c` (whose headers are `heads`), the extents of its records when it
    /// holds several references, else `None`. The container is parsed once. A read covers at
    /// least its first base, and an unmapped read only that; an unplaced read is at 0.
    pub fn extents(
        &mut self,
        c: &RawContainer,
        heads: &[SliceHead],
    ) -> Result<Vec<Option<Extents>>, IoError> {
        let mut out: Vec<Option<Extents>> = vec![None; heads.len()];
        if heads.iter().all(|h| h.ref_id != -2) {
            return Ok(out);
        }
        let _note = decoder_note(&self.info.name, c);
        let name = &self.info.name;
        let bad = |e: &dyn std::fmt::Display| {
            detect_err(
                name,
                format!(
                    "{name}: {} cannot be decoded: {e}; the file is corrupt",
                    c.place()
                ),
            )
        };
        let comp = compression(c);
        if let Some(what) = &comp.unsupported {
            return Err(unsupported_codec(name, c, what));
        }
        check_compression_block(c).map_err(|m| bad(&m))?;
        let mut container = cram::io::reader::Container::default();
        cram::io::Reader::new(&c.bytes[..])
            .read_container(&mut container)
            .map_err(|e| bad(&e))?;
        let ch = container.compression_header().map_err(|e| bad(&e))?;
        for ((head, slice), slot) in heads.iter().zip(container.slices()).zip(&mut out) {
            if head.ref_id != -2 {
                continue;
            }
            let slice = slice.map_err(|e| bad(&e))?;
            check_slice_blocks(c, head, comp.external.get(b"RN").copied()).map_err(|m| bad(&m))?;
            let (core, ext) = slice.decode_blocks().map_err(|e| bad(&e))?;
            precheck(head, self.info.aln.refs.len(), &comp, &core, &ext).map_err(|m| bad(&m))?;
            self.repo.clear();
            if let Ok(mut l) = self.state.loaded.lock() {
                l.clear();
            }
            let records = slice
                .records(self.repo.clone(), &self.info.sam, &ch, &core, &ext)
                .map_err(|e| bad(&e))?;
            let mut v = Vec::with_capacity(records.len());
            for r in &records {
                let flags = r.flags().map_err(|e| bad(&e))?;
                let rid = match r.reference_sequence_id(&self.info.sam) {
                    Some(Ok(id)) => id as i32,
                    Some(Err(e)) => return Err(bad(&e)),
                    None => -1,
                };
                let start = match r.alignment_start() {
                    Some(Ok(p)) => usize::from(p) as i64,
                    Some(Err(e)) => return Err(bad(&e)),
                    None => 0,
                };
                // htslib: a read covers at least its first base, and an unmapped read only that.
                let span = if flags.is_unmapped() {
                    1
                } else {
                    match r.alignment_span() {
                        Some(Ok(s)) => (s as i64).max(1),
                        _ => 1,
                    }
                };
                v.push((rid, start, start + span - 1));
            }
            *slot = Some(v);
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn head(ref_id: i32, records: usize) -> SliceHead {
        SliceHead {
            ref_id,
            start: 1,
            span: 1,
            records,
            embedded: -1,
            landmark: 0,
            size: 0,
        }
    }

    /// A slice's external blocks.
    type Blocks = Vec<(i32, Cow<'static, [u8]>)>;

    /// Data series in external blocks of their own: `(key, block id, bytes)`.
    fn blocks(series: &[(&[u8; 2], i32, Vec<u8>)]) -> (CompressionInfo, Blocks) {
        let mut comp = CompressionInfo {
            required: true,
            ..CompressionInfo::default()
        };
        let mut ext = Vec::new();
        for (key, id, bytes) in series {
            comp.external.insert(**key, *id);
            *comp.uses.entry(*id).or_default() += 1;
            ext.push((*id, Cow::Owned(bytes.clone())));
        }
        (comp, ext)
    }

    const MINUS_ONE: [u8; 5] = [0xff, 0xff, 0xff, 0xff, 0x0f];

    /// The checks that keep noodles from allocating or indexing by a corrupt count.
    #[test]
    fn prechecks_refuse_counts_the_slice_cannot_hold() {
        // A record count far past the slice's bytes.
        let (comp, ext) = blocks(&[]);
        let e = precheck(&head(0, 5000), 1, &comp, &[], &ext).unwrap_err();
        assert!(e.contains("claims 5,000 records in 0 bytes"), "{e}");
        assert!(precheck(&head(0, 4096), 1, &comp, &[], &ext).is_ok());
        // Features: 5 + 6 claimed, 10 feature codes stored.
        let (comp, ext) = blocks(&[(b"FN", 1, vec![5, 6]), (b"FC", 2, vec![b'X'; 10])]);
        let e = precheck(&head(0, 2), 1, &comp, &[], &ext).unwrap_err();
        assert!(
            e.contains("claim 11 features, and the slice holds 10"),
            "{e}"
        );
        let (comp, ext) = blocks(&[(b"FN", 1, vec![5, 6]), (b"FC", 2, vec![b'X'; 11])]);
        assert!(precheck(&head(0, 2), 1, &comp, &[], &ext).is_ok());
        // A constant feature count no read could hold (one-symbol Huffman, no bytes).
        let (mut comp, ext) = blocks(&[(b"FC", 2, vec![b'X'; 11])]);
        comp.constant.insert(*b"FN", 1 << 30);
        assert!(precheck(&head(0, 2), 1, &comp, &[], &ext).is_err());
        // A mate five reads further on, in a slice of three.
        let (comp, ext) = blocks(&[(b"CF", 1, vec![0x4, 0, 0]), (b"NF", 2, vec![5])]);
        let e = precheck(&head(0, 3), 1, &comp, &[], &ext).unwrap_err();
        assert!(
            e.contains("read 1 of a slice of 3 names its mate 5 reads further on"),
            "{e}"
        );
        let (comp, ext) = blocks(&[(b"CF", 1, vec![0x4, 0, 0]), (b"NF", 2, vec![1])]);
        assert!(precheck(&head(0, 3), 1, &comp, &[], &ext).is_ok());
        // A series of its own whose count disagrees with the slice's.
        let (comp, ext) = blocks(&[(b"BF", 1, vec![0, 0, 0])]);
        let e = precheck(&head(0, 4), 1, &comp, &[], &ext).unwrap_err();
        assert!(
            e.contains("a slice claims 4 records, and its BF series holds 3 values"),
            "{e}"
        );
        // A detached read's mate is not a distance.
        let (comp, ext) = blocks(&[(b"CF", 1, vec![0x6, 0, 0]), (b"NF", 2, vec![5])]);
        assert!(precheck(&head(0, 3), 1, &comp, &[], &ext).is_ok());
    }

    /// In a slice of several references: a mapped read with no reference is refused (noodles
    /// `expect`s one), and the references of the mapped reads are found.
    #[test]
    fn prechecks_find_the_references_of_a_multi_reference_slice() {
        let mut ri = MINUS_ONE.to_vec();
        ri.push(1);
        let (comp, ext) = blocks(&[(b"BF", 1, vec![0, 4]), (b"RI", 2, ri.clone())]);
        let e = precheck(&head(-2, 2), 3, &comp, &[], &ext).unwrap_err();
        assert!(
            e.contains("mapped read 1 of a slice of several references"),
            "{e}"
        );
        let (comp, ext) = blocks(&[(b"BF", 1, vec![4, 0]), (b"RI", 2, ri)]);
        assert_eq!(
            precheck(&head(-2, 2), 3, &comp, &[], &ext),
            Ok(Some(vec![1]))
        );
        // RI coded as a constant (one-symbol Huffman), as htslib writes it when all agree.
        let (mut comp, ext) = blocks(&[(b"BF", 1, vec![0, 0])]);
        comp.constant.insert(*b"RI", 2);
        assert_eq!(
            precheck(&head(-2, 2), 3, &comp, &[], &ext),
            Ok(Some(vec![2]))
        );
        // A block that another series reads too is not trusted for values.
        let (mut comp, ext) = blocks(&[(b"BF", 1, vec![0, 0]), (b"RI", 2, vec![0, 1])]);
        comp.uses.insert(2, 2);
        assert_eq!(precheck(&head(-2, 2), 3, &comp, &[], &ext), Ok(None));
    }
}
