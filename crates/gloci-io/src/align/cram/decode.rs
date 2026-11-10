//! Decoding CRAM containers to BAM records: the slice checks before noodles decodes, the reads'
//! extents against their reference bases (padded with `N` where they run past), the
//! transcoding with htslib's derived fields, and the position-only decoding the `.crai` needs.

use super::super::encode::reg2bin;
use super::frame::{
    CompressionInfo, RawContainer, SliceHead, compression, detect_err, itf8_values, slice_heads,
    unsupported_codec,
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
use std::sync::atomic::Ordering;

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
    /// The repository's sequences need this many `N` after their end.
    Pad(usize),
    /// The embedded bases need to run to this reference position (1-based).
    EmbeddedTo(i64),
    /// The placeholder of each of these references needs this length.
    Placeholder(Vec<(usize, usize)>),
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
    let flags = series(comp, ext, b"BF", n).filter(|v| v.len() == n);
    let cram_flags = series(comp, ext, b"CF", n).filter(|v| v.len() == n);
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
    let (Some(bf), Some(ri)) = (flags, series(comp, ext, b"RI", n)) else {
        return Ok(None);
    };
    if ri.len() != head.records {
        return Ok(None);
    }
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

/// Decodes containers to BAM records. Each worker of the pool has one.
pub(crate) struct Decoder {
    pub(crate) info: Arc<CramInfo>,
    /// The references shared by every decoder (slices that need one reference).
    pub(crate) cache: Arc<RefCache>,
    /// This decoder's own repository: slices that need several references, are stored without
    /// the reference, or have reads past a reference's end.
    pub(crate) repo: fasta::Repository,
    pub(crate) state: Arc<RefState>,
    /// What `repo` held when the shared count last heard of it.
    pub(crate) private_len: usize,
    /// The reference id of the last slice decoded on `repo`.
    pub(crate) last_ref: Option<i32>,
    /// The `RR` of the last container decoded.
    pub(crate) last_required: Option<bool>,
    pub(crate) bam: noodles_bam::io::Writer<Vec<u8>>,
}

impl Decoder {
    pub(crate) fn new(info: Arc<CramInfo>, cache: Arc<RefCache>) -> Result<Decoder, IoError> {
        let state = Arc::new(RefState::default());
        let adapter = RefAdapter::new(info.clone(), state.clone(), cache.counters.clone())?;
        Ok(Decoder {
            info,
            cache,
            repo: fasta::Repository::new(adapter),
            state,
            private_len: 0,
            last_ref: None,
            last_required: None,
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

    /// Tell the shared counts what this decoder's own repository holds now.
    fn sync_private(&mut self) {
        let now = self.repo.len();
        self.cache.note_live(self.private_len, now);
        self.private_len = now;
    }

    /// Empty this decoder's repository (and with `reset`, forget its padding and placeholders).
    pub(crate) fn clear(&mut self, reset: bool) {
        self.repo.clear();
        if let Ok(mut c) = self.state.cached.lock() {
            c.clear();
        }
        self.state.placeholder_bytes.store(0, Ordering::Relaxed);
        if reset {
            self.state.pad.store(0, Ordering::Relaxed);
            if let Ok(mut p) = self.state.placeholders.lock() {
                p.clear();
            }
        }
        self.sync_private();
    }

    /// Where the reads of `head` find their bases, after checking that they can (D109).
    pub(crate) fn bases(
        &mut self,
        c: &RawContainer,
        head: &SliceHead,
        required: bool,
    ) -> Result<Bases, IoError> {
        // This decoder's repository starts afresh on another reference; a slice of several
        // references keeps what it holds, and the adapter asks for a clear if the slice needs
        // anything else.
        if head.ref_id >= 0 && self.last_ref != Some(head.ref_id) {
            let name = &self.info.aln.refs[head.ref_id as usize].0;
            let keep = self
                .state
                .cached
                .lock()
                .is_ok_and(|c| c.iter().all(|n| n == name));
            if !keep {
                self.clear(true);
            }
        }
        self.last_ref = Some(head.ref_id);
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
        let heads = slice_heads(c, &self.info.name, self.info.aln.refs.len())?;
        let comp = compression(c);
        if let Some(what) = &comp.unsupported {
            return Err(unsupported_codec(&self.info.name, c, what));
        }
        let required = comp.required;
        self.state.required.store(required, Ordering::Relaxed);
        if self.last_required != Some(required) {
            // Placeholders never serve a container that reads the reference, nor the reverse.
            self.clear(true);
            self.last_ref = None;
            self.last_required = Some(required);
        }
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
            // A slice whose reads need one reference reads it from the shared cache.
            let one = match (bases, head.ref_id, &needed) {
                (Bases::Repository, r, _) if r >= 0 => Some(r as usize),
                (Bases::Repository, _, Some(v)) if v.len() == 1 => Some(v[0]),
                _ => None,
            };
            let lease = match one {
                Some(rid) => Some(cache.acquire(rid, c.number)?),
                None => None,
            };
            let mut private = lease.is_none();
            // A slice that needs a sequence this decoder's repository lacks, while it holds
            // another, is decoded again after a clear; reads past the end of their reference
            // bases are decoded again against `N`, on this decoder's repository.
            let mut padded: Option<Vec<(i32, Cow<'_, [u8]>)>> = None;
            let mut done = false;
            self.state.grow.store(false, Ordering::Relaxed);
            self.state.evict.store(false, Ordering::Relaxed);
            for _ in 0..6 {
                let srcs = padded.as_deref().unwrap_or(&ext);
                let (repo, state) = match (&lease, private) {
                    (Some(l), false) => (l.repo.clone(), l.state.clone()),
                    _ => (self.repo.clone(), self.state.clone()),
                };
                let records = match slice.records(repo, &info.sam, &ch, &core, srcs) {
                    Ok(r) => r,
                    Err(_) if private && self.state.evict.swap(false, Ordering::Relaxed) => {
                        self.clear(false);
                        self.state.grow.store(true, Ordering::Relaxed);
                        continue;
                    }
                    Err(e) => return Err(self.decode_error(c, head, &state, e)),
                };
                self.sync_private();
                let pad = if private {
                    self.state.pad.load(Ordering::Relaxed)
                } else {
                    0
                };
                match self.extent(c, head, bases, pad, &records, srcs)? {
                    Extent::Fits => {
                        self.transcode(c, &records, out)?;
                        done = true;
                        break;
                    }
                    Extent::Pad(pad) => {
                        drop(records);
                        private = true;
                        self.clear(false);
                        self.state.grow.store(true, Ordering::Relaxed);
                        self.state.pad.store(pad, Ordering::Relaxed);
                    }
                    Extent::Placeholder(lens) => {
                        drop(records);
                        self.clear(false);
                        self.state.grow.store(true, Ordering::Relaxed);
                        if let Ok(mut p) = self.state.placeholders.lock() {
                            for (rid, len) in lens {
                                let e = p.entry(info.aln.refs[rid].0.clone()).or_default();
                                *e = (*e).max(len);
                            }
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
            drop(lease);
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
    /// (sequences `pad` bases longer than their reference). A read may run past its
    /// reference's end, and is decoded against `N` there, by at most 1 Mb or four times the
    /// slice's longest read; further is an error.
    pub(crate) fn extent(
        &self,
        c: &RawContainer,
        head: &SliceHead,
        bases: Bases,
        pad: usize,
        records: &[cram::Record<'_>],
        srcs: &[(i32, Cow<'_, [u8]>)],
    ) -> Result<Extent, IoError> {
        if bases == Bases::Stored {
            return Ok(Extent::Fits);
        }
        let pad = pad as i64;
        let embedded_len = block(srcs, head.embedded).map_or(0, |b| b.len() as i64);
        let placeholders = self.state.placeholders.lock().ok();
        let placeholder = |rid: usize| {
            placeholders
                .as_ref()
                .and_then(|p| p.get(&self.info.aln.refs[rid].0).copied())
                .unwrap_or(0) as i64
        };
        // (reference, how far past its end) of the read that runs furthest past one.
        let mut over: Option<(usize, i64)> = None;
        let mut longest = 0i64;
        let mut need_pad = 0i64;
        let mut need_end = 0i64;
        let mut need_len: Vec<(usize, usize)> = Vec::new();
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
            longest = longest.max(r.quality_scores().len() as i64);
            let rid = match r.reference_sequence_id(&self.info.sam) {
                Some(Ok(id)) if id < self.info.aln.refs.len() => id,
                _ => return Err(self.corrupt(c, &"a mapped read has no reference")),
            };
            let len = self.info.aln.refs[rid].1 as i64;
            if end > len && over.is_none_or(|(_, o)| end - len > o) {
                over = Some((rid, end - len));
            }
            match bases {
                Bases::Repository => {
                    if end > len + pad {
                        need_pad = need_pad.max(end - len);
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
                Bases::Placeholder => {
                    let want = start.max(end);
                    if want > placeholder(rid) {
                        match need_len.iter_mut().find(|(r, _)| *r == rid) {
                            Some((_, l)) => *l = (*l).max(want as usize),
                            None => need_len.push((rid, want as usize)),
                        }
                    }
                }
                Bases::Stored => {}
            }
        }
        let cap = MIN_OVERHANG.max(4 * longest);
        if let Some((rid, o)) = over
            && o > cap
        {
            return Err(self.corrupt(
                c,
                &format!(
                    "a read runs {} bases past the end of {} ({} bases); gloci decodes reads at most {} bases past the end of their reference",
                    thousands(o as u64),
                    self.info.chrom(rid),
                    thousands(self.info.aln.refs[rid].1),
                    thousands(cap as u64)
                ),
            ));
        }
        if bases == Bases::Embedded && need_end > 0 {
            let past = need_end - (i64::from(head.start) + embedded_len - 1);
            if past > cap {
                return Err(self.corrupt(
                    c,
                    &format!(
                        "a read runs {} bases past the reference bases embedded in its slice; gloci decodes reads at most {} bases past them",
                        thousands(past as u64),
                        thousands(cap as u64)
                    ),
                ));
            }
        }
        Ok(if need_pad > 0 {
            Extent::Pad(need_pad as usize)
        } else if need_end > 0 {
            Extent::EmbeddedTo(need_end)
        } else if !need_len.is_empty() {
            Extent::Placeholder(need_len)
        } else {
            Extent::Fits
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
        let state = Arc::new(RefState::default());
        state.required.store(false, Ordering::Relaxed);
        state.grow.store(true, Ordering::Relaxed);
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
            let (core, ext) = slice.decode_blocks().map_err(|e| bad(&e))?;
            precheck(head, self.info.aln.refs.len(), &comp, &core, &ext).map_err(|m| bad(&m))?;
            self.repo.clear();
            if let Ok(mut c) = self.state.cached.lock() {
                c.clear();
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

    /// Data series in external blocks of their own: `(key, block id, bytes)`.
    fn blocks(
        series: &[(&[u8; 2], i32, Vec<u8>)],
    ) -> (CompressionInfo, Vec<(i32, Cow<'static, [u8]>)>) {
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
