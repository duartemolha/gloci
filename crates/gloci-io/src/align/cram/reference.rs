//! The CRAM reference (D109): what every decoder of one CRAM shares ([`CramInfo`]: the header,
//! the `--reference` FASTA checked against it, the D109 messages), and the noodles repository
//! adapter over the FASTA ([`RefAdapter`]), which never answers `None`.

use super::super::header::AlnHeader;
use super::frame::detect_err;
use crate::error::IoError;
use crate::fasta::IndexedFasta;
use gloci_core::{Interval, thousands};
use noodles_fasta as fasta;
use noodles_sam as sam;
use std::io;
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

// ---- the CRAM header and the reference ----

/// What every decoder of one CRAM shares.
pub(crate) struct CramInfo {
    /// The input's name, for messages.
    pub(crate) name: String,
    /// The header for noodles: the references and read groups, in file order.
    pub(crate) sam: sam::Header,
    /// The header as read (references named as in the file).
    pub(crate) aln: AlnHeader,
    /// Each reference's `M5` and `UR`, for messages.
    m5: Vec<Option<String>>,
    pub(crate) ur: Vec<Option<String>>,
    /// The `--reference` FASTA: a handle each decoder reopens, and its name.
    pub(crate) fasta: Option<(Mutex<IndexedFasta>, String)>,
    /// Whether the FASTA has each reference.
    pub(crate) in_fasta: Vec<bool>,
}

pub(crate) fn field<'a>(line: &'a [u8], key: &[u8]) -> Option<&'a [u8]> {
    line.split(|&b| b == b'\t')
        .skip(1)
        .find_map(|f| f.strip_prefix(key))
}

impl CramInfo {
    /// Build from the header `h` of the CRAM `name`, opening and checking `reference` (D109).
    pub(crate) fn new(
        name: &str,
        h: AlnHeader,
        reference: Option<&Path>,
    ) -> Result<CramInfo, IoError> {
        let n = h.refs.len();
        let mut m5 = Vec::with_capacity(n);
        let mut ur = Vec::with_capacity(n);
        let mut groups = Vec::new();
        for line in h.text.split(|&b| b == b'\n') {
            let line = line.strip_suffix(b"\r").unwrap_or(line);
            let text = |v: &[u8]| String::from_utf8_lossy(v).into_owned();
            if line.starts_with(b"@SQ\t") {
                m5.push(field(line, b"M5:").map(text));
                ur.push(field(line, b"UR:").map(text));
            } else if line.starts_with(b"@RG\t")
                && let Some(id) = field(line, b"ID:")
            {
                groups.push(id.to_vec());
            }
        }
        m5.resize(n, None);
        ur.resize(n, None);
        let mut b = sam::Header::builder();
        for (r, len) in &h.refs {
            let len = usize::try_from(*len)
                .ok()
                .and_then(NonZeroUsize::new)
                .ok_or_else(|| {
                    detect_err(
                        name,
                        format!(
                            "{name}: the @SQ line of {} has LN:{len}; a CRAM reference needs a \
length of at least 1",
                            String::from_utf8_lossy(r)
                        ),
                    )
                })?;
            b =
                b.add_reference_sequence(
                    r.clone(),
                    sam::header::record::value::Map::<
                        sam::header::record::value::map::ReferenceSequence,
                    >::new(len),
                );
        }
        for g in groups {
            b = b.add_read_group(g, sam::header::record::value::Map::default());
        }
        let mut in_fasta = vec![false; n];
        let fasta = match reference {
            None => None,
            Some(p) => {
                let f = IndexedFasta::open(p)?;
                let shown = p.display().to_string();
                for (i, (r, len)) in h.refs.iter().enumerate() {
                    let Some(flen) = f.len_of(r) else {
                        continue;
                    };
                    if flen != *len {
                        let made = ur[i]
                            .as_deref()
                            .map(|u| format!(" (its header names {u})"))
                            .unwrap_or_default();
                        return Err(detect_err(
                            name,
                            format!(
                                "{name} was not made with {shown}: {} is {} bases long in its \
header and {} in {shown}; pass the FASTA it was made with{made}: --reference FASTA",
                                String::from_utf8_lossy(r),
                                thousands(*len),
                                thousands(flen)
                            ),
                        ));
                    }
                    in_fasta[i] = true;
                }
                Some((Mutex::new(f), shown))
            }
        };
        Ok(CramInfo {
            name: name.to_string(),
            sam: b.build(),
            aln: h,
            m5,
            ur,
            fasta,
            in_fasta,
        })
    }

    pub(crate) fn chrom(&self, rid: usize) -> String {
        String::from_utf8_lossy(&self.aln.refs[rid].0).into_owned()
    }

    /// The D109 error for a reference the decoder needs and does not have.
    pub(crate) fn needs_reference(&self, name: &[u8]) -> IoError {
        let chrom = String::from_utf8_lossy(name);
        let rid = self.aln.ref_id(name);
        let mut ids = Vec::new();
        if let Some(m) = rid.and_then(|i| self.m5[i].as_deref()) {
            ids.push(format!("M5 {m}"));
        }
        if let Some(u) = rid.and_then(|i| self.ur[i].as_deref()) {
            ids.push(format!("UR {u}"));
        }
        let ids = if ids.is_empty() {
            String::new()
        } else {
            format!(" ({})", ids.join(", "))
        };
        let cram = &self.name;
        let msg = match &self.fasta {
            None => format!(
                "{cram} needs the reference sequence of {chrom}{ids} to decode its reads; pass \
it with --reference FASTA"
            ),
            Some((_, fa)) => format!(
                "{cram} needs the reference sequence of {chrom}{ids} to decode its reads, and \
{fa} has no sequence {chrom}; pass the FASTA it was made with: --reference FASTA"
            ),
        };
        detect_err(cram, msg)
    }

    /// The D109 error for a slice whose reference bases do not have the checksum it stores.
    pub(crate) fn mismatch(&self, rid: usize, start: i64, end: i64) -> IoError {
        let fa = self
            .fasta
            .as_ref()
            .map_or("the reference", |(_, n)| n.as_str());
        let chrom = self.chrom(rid);
        detect_err(
            &self.name,
            format!(
                "the sequence {chrom} in {fa} does not match the one {} was made with \
({chrom}:{}-{}); pass the FASTA it was made with: --reference FASTA",
                self.name,
                thousands(start.max(0) as u64),
                thousands(end.max(0) as u64)
            ),
        )
    }

    pub(crate) fn ref_count(&self) -> usize {
        self.aln.refs.len()
    }

    pub(crate) fn name(&self) -> &str {
        &self.name
    }
}

/// What a decoder and its [`RefAdapter`] share.
#[derive(Default)]
pub(crate) struct RefState {
    /// `N` bases added after each sequence, for reads that run past the end of their reference
    /// (htslib reads them against `N`).
    pub(crate) pad: AtomicUsize,
    /// The `RR` of the container being decoded: without it no bases are fetched.
    pub(crate) required: AtomicBool,
    /// Without `RR`, the length of the `N` placeholder given for each name. noodles still asks
    /// for the sequence of each mapped read of a slice of several references, and slices it at
    /// the read's positions though it reads no base, so the placeholder must reach them.
    pub(crate) placeholder: AtomicUsize,
    /// The adapter's last failure, with gloci's message (noodles passes on only an `io::Error`).
    pub(crate) failure: Mutex<Option<IoError>>,
    /// The names fetched since the repository was last cleared.
    pub(crate) cached: Mutex<Vec<Vec<u8>>>,
    /// The slice being decoded may add sequences to the repository. Without it, a fetch while
    /// the repository holds another sequence fails with `evict` set, and the decoder clears the
    /// repository and decodes the slice again: a decoder holds only what one slice needs.
    pub(crate) grow: AtomicBool,
    pub(crate) evict: AtomicBool,
    /// Sequences read from the FASTA.
    pub(crate) fetches: AtomicUsize,
    /// The most sequences the decoder's repository held at once.
    pub(crate) max_live: AtomicUsize,
}

/// The repository adapter over the `--reference` FASTA (D109). It never answers `None`.
pub(crate) struct RefAdapter {
    pub(crate) info: Arc<CramInfo>,
    pub(crate) fasta: Option<Mutex<IndexedFasta>>,
    pub(crate) state: Arc<RefState>,
}

pub(crate) const ADAPTER_FAILED: &str = "the reference adapter failed";
pub(crate) const EVICT: &str = "the repository holds another sequence";

impl RefAdapter {
    pub(crate) fn fetch(&mut self, name: &[u8]) -> Result<fasta::Record, IoError> {
        let record = |seq: Vec<u8>| {
            fasta::Record::new(
                fasta::record::Definition::new(name.to_vec(), None),
                fasta::record::Sequence::from(seq),
            )
        };
        if !self.state.required.load(Ordering::Relaxed) {
            // A container stored without the reference reads no base of it.
            if let Ok(mut c) = self.state.cached.lock() {
                c.push(name.to_vec());
            }
            return Ok(record(vec![
                b'N';
                self.state.placeholder.load(Ordering::Relaxed)
            ]));
        }
        let Some(f) = self.fasta.as_mut().and_then(|m| m.get_mut().ok()) else {
            return Err(self.info.needs_reference(name));
        };
        let Some(len) = f.len_of(name) else {
            return Err(self.info.needs_reference(name));
        };
        let mut cached = self
            .state
            .cached
            .lock()
            .map_err(|_| IoError::Internal("the reference cache lock is poisoned".into()))?;
        if !cached.is_empty() && !self.state.grow.load(Ordering::Relaxed) {
            self.state.evict.store(true, Ordering::Relaxed);
            return Err(IoError::Internal(EVICT.into()));
        }
        cached.push(name.to_vec());
        drop(cached);
        let pad = self.state.pad.load(Ordering::Relaxed);
        let mut seq = Vec::with_capacity(usize::try_from(len).unwrap_or(0).saturating_add(pad));
        f.fetch(
            name,
            Interval { start: 0, end: len },
            &mut |chunk: &[u8]| {
                seq.extend(chunk.iter().map(u8::to_ascii_uppercase));
                Ok(())
            },
        )?;
        seq.resize(seq.len() + pad, b'N');
        self.state.fetches.fetch_add(1, Ordering::Relaxed);
        Ok(record(seq))
    }
}

impl fasta::repository::Adapter for RefAdapter {
    fn get(&mut self, name: &[u8]) -> Option<io::Result<fasta::Record>> {
        Some(self.fetch(name).map_err(|e| {
            if !self.state.evict.load(Ordering::Relaxed)
                && let Ok(mut slot) = self.state.failure.lock()
                && slot.is_none()
            {
                *slot = Some(e);
            }
            io::Error::other(ADAPTER_FAILED)
        }))
    }
}
