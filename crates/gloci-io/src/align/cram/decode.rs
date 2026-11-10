//! Decoding CRAM containers to BAM records: the slice checks before noodles decodes, the reads'
//! extents against their reference bases (padded with `N` where they run past), the
//! transcoding with htslib's derived fields, and the position-only decoding the `.crai` needs.

use super::super::encode::reg2bin;
use super::frame::{
    CompressionInfo, RawContainer, SliceHead, compression, detect_err, slice_heads,
    unsupported_codec,
};
use super::reference::{CramInfo, RefAdapter, RefState};
use crate::error::IoError;
use gloci_core::thousands;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam::alignment::Record as _;
use noodles_sam::alignment::io::Write as _;
use std::borrow::Cow;
use std::io;
use std::sync::atomic::Ordering;
use std::sync::{Arc, Mutex};

// ---- decoding ----

/// Where a slice's reads find their reference bases.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Bases {
    /// The repository (`--reference`).
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
    /// The placeholder sequences need this length.
    Placeholder(usize),
}

/// Decodes containers to BAM records. Each worker of the pool has one.
pub(crate) struct Decoder {
    pub(crate) info: Arc<CramInfo>,
    pub(crate) repo: fasta::Repository,
    pub(crate) state: Arc<RefState>,
    /// The reference id of the last slice decoded.
    pub(crate) last_ref: Option<i32>,
    /// The `RR` of the last container decoded.
    pub(crate) last_required: Option<bool>,
    pub(crate) bam: noodles_bam::io::Writer<Vec<u8>>,
}

impl Decoder {
    pub(crate) fn new(info: Arc<CramInfo>) -> Result<Decoder, IoError> {
        let fasta = match &info.fasta {
            None => None,
            Some((f, _)) => {
                let f = f
                    .lock()
                    .map_err(|_| IoError::Internal("the reference lock is poisoned".into()))?;
                Some(Mutex::new(f.reopen()?))
            }
        };
        let state = Arc::new(RefState::default());
        let adapter = RefAdapter {
            info: info.clone(),
            fasta,
            state: state.clone(),
        };
        Ok(Decoder {
            info,
            repo: fasta::Repository::new(adapter),
            state,
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

    /// The gloci error for a noodles failure while decoding `head` of `c`.
    pub(crate) fn decode_error(&self, c: &RawContainer, head: &SliceHead, e: io::Error) -> IoError {
        if let Some(f) = self.state.failure.lock().ok().and_then(|mut s| s.take()) {
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

    /// Empty the repository (and with `reset_pad`, forget the padding).
    pub(crate) fn clear(&mut self, reset_pad: bool) {
        self.repo.clear();
        if let Ok(mut c) = self.state.cached.lock() {
            c.clear();
        }
        if reset_pad {
            self.state.pad.store(0, Ordering::Relaxed);
            self.state.placeholder.store(0, Ordering::Relaxed);
        }
    }

    /// Where the reads of `head` find their bases, after checking that they can (D109).
    pub(crate) fn bases(
        &mut self,
        c: &RawContainer,
        head: &SliceHead,
        required: bool,
    ) -> Result<Bases, IoError> {
        // One sequence at a time (D109): a slice on another reference starts afresh. A slice of
        // several references keeps what the repository holds; the adapter asks for a clear if
        // the slice needs anything else.
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
        let CompressionInfo {
            required,
            unsupported,
        } = compression(c);
        if let Some(what) = unsupported {
            return Err(unsupported_codec(&self.info.name, c, &what));
        }
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
        // The records borrow the header; the decoder itself stays free to transcode them.
        let info = self.info.clone();
        for (head, slice) in heads.iter().zip(container.slices()) {
            let slice = slice.map_err(|e| self.corrupt(c, &e))?;
            let bases = self.bases(c, head, required)?;
            let (core, ext) = slice.decode_blocks().map_err(|e| self.corrupt(c, &e))?;
            if bases == Bases::Embedded && !ext.iter().any(|(id, _)| *id == head.embedded) {
                // noodles finds the embedded block by its id, and `expect`s it.
                return Err(self.corrupt(
                    c,
                    &"a slice names an embedded reference block it does not have",
                ));
            }
            // A slice that needs a sequence the repository lacks, while it holds another, is
            // decoded again after a clear; reads past the end of their reference bases are
            // decoded again against `N`.
            let mut padded: Option<Vec<(i32, Cow<'_, [u8]>)>> = None;
            let mut done = false;
            self.state.grow.store(false, Ordering::Relaxed);
            self.state.evict.store(false, Ordering::Relaxed);
            for _ in 0..4 {
                let srcs = padded.as_deref().unwrap_or(&ext);
                let records = match slice.records(self.repo.clone(), &info.sam, &ch, &core, srcs) {
                    Ok(r) => r,
                    Err(_) if self.state.evict.swap(false, Ordering::Relaxed) => {
                        self.clear(false);
                        self.state.grow.store(true, Ordering::Relaxed);
                        continue;
                    }
                    Err(e) => return Err(self.decode_error(c, head, e)),
                };
                self.state
                    .max_live
                    .fetch_max(self.repo.len(), Ordering::Relaxed);
                match self.extent(c, head, bases, &records, srcs)? {
                    Extent::Fits => {
                        self.transcode(c, &records, out)?;
                        done = true;
                        break;
                    }
                    Extent::Pad(pad) => {
                        drop(records);
                        self.clear(false);
                        self.state.grow.store(true, Ordering::Relaxed);
                        self.state.pad.store(pad, Ordering::Relaxed);
                    }
                    Extent::Placeholder(len) => {
                        drop(records);
                        self.clear(false);
                        self.state.grow.store(true, Ordering::Relaxed);
                        self.state.placeholder.store(len, Ordering::Relaxed);
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

    /// Whether every mapped read of `records` lies within the reference bases it decodes with.
    pub(crate) fn extent(
        &self,
        c: &RawContainer,
        head: &SliceHead,
        bases: Bases,
        records: &[cram::Record<'_>],
        srcs: &[(i32, Cow<'_, [u8]>)],
    ) -> Result<Extent, IoError> {
        if bases == Bases::Stored {
            return Ok(Extent::Fits);
        }
        let pad = self.state.pad.load(Ordering::Relaxed) as i64;
        let embedded_len = srcs
            .iter()
            .find(|(id, _)| *id == head.embedded)
            .map_or(0, |(_, s)| s.len() as i64);
        let placeholder = self.state.placeholder.load(Ordering::Relaxed) as i64;
        let mut need_pad = 0i64;
        let mut need_end = 0i64;
        let mut need_len = 0i64;
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
            match bases {
                Bases::Repository => {
                    let len = match r.reference_sequence_id(&self.info.sam) {
                        Some(Ok(id)) if id < self.info.aln.refs.len() => {
                            self.info.aln.refs[id].1 as i64
                        }
                        _ => return Err(self.corrupt(c, &"a mapped read has no reference")),
                    };
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
                    if start.max(end) > placeholder {
                        need_len = need_len.max(start.max(end));
                    }
                }
                Bases::Stored => {}
            }
        }
        Ok(if need_pad > 0 {
            Extent::Pad(need_pad as usize)
        } else if need_end > 0 {
            Extent::EmbeddedTo(need_end)
        } else if need_len > 0 {
            Extent::Placeholder(need_len as usize)
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

/// The records of a slice that holds several references, for the `.crai`: each record's
/// reference id (-1 for none), start and end, as htslib's index builder sees them. No reference
/// bases are read: the adapter is told the container needs none (so it gives empty
/// placeholders), and the reads' bases are never put together.
pub(crate) struct ExtentDecoder {
    info: Arc<CramInfo>,
    repo: fasta::Repository,
}

impl ExtentDecoder {
    pub fn new(info: Arc<CramInfo>) -> ExtentDecoder {
        let state = Arc::new(RefState::default());
        state.required.store(false, Ordering::Relaxed);
        let adapter = RefAdapter {
            info: info.clone(),
            fasta: None,
            state,
        };
        ExtentDecoder {
            info,
            repo: fasta::Repository::new(adapter),
        }
    }

    pub fn info(&self) -> &CramInfo {
        &self.info
    }

    /// `(ref id, start, end)` of each record of slice `index` (0-based) of `c`, a slice of
    /// several references. A read covers at least its first base, and an unmapped read only
    /// that; an unplaced read is at 0.
    pub fn extents(
        &mut self,
        c: &RawContainer,
        index: usize,
    ) -> Result<Vec<(i32, i64, i64)>, IoError> {
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
        if let Some(what) = compression(c).unsupported {
            return Err(unsupported_codec(name, c, &what));
        }
        let mut container = cram::io::reader::Container::default();
        cram::io::Reader::new(&c.bytes[..])
            .read_container(&mut container)
            .map_err(|e| bad(&e))?;
        let ch = container.compression_header().map_err(|e| bad(&e))?;
        let Some(slice) = container.slices().nth(index) else {
            return Err(bad(&"a slice is missing"));
        };
        let slice = slice.map_err(|e| bad(&e))?;
        let (core, ext) = slice.decode_blocks().map_err(|e| bad(&e))?;
        self.repo.clear();
        let records = slice
            .records(self.repo.clone(), &self.info.sam, &ch, &core, &ext)
            .map_err(|e| bad(&e))?;
        let mut out = Vec::with_capacity(records.len());
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
            out.push((rid, start, start + span - 1));
        }
        Ok(out)
    }
}
