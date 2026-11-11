//! The CRAM reference (D109): what every decoder of one CRAM shares ([`CramInfo`]: the header,
//! the `--reference` FASTA checked against it, the D109 messages), and the noodles repository
//! adapter over the FASTA ([`RefAdapter`]), which never answers `None`.

use super::super::header::{AlnHeader, SortOrder};
use super::frame::detect_err;
use crate::error::IoError;
use crate::fasta::IndexedFasta;
use gloci_core::{Interval, thousands};
use noodles_fasta as fasta;
use noodles_sam as sam;
use std::collections::{HashMap, HashSet};
use std::io;
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex, Weak};

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

    /// The name of reference `rid`, or "unplaced reads" for -1.
    pub(crate) fn chrom_or_unplaced(&self, rid: i32) -> String {
        usize::try_from(rid)
            .ok()
            .filter(|&r| r < self.aln.refs.len())
            .map_or_else(|| "unplaced reads".to_string(), |r| self.chrom(r))
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

/// Counts shared by every repository of one CRAM, for the memory tests (D109).
#[derive(Debug, Default)]
pub(crate) struct Counters {
    /// Sequences read from the FASTA.
    pub(crate) fetches: AtomicUsize,
    /// The most sequences held at once, by all decoders together.
    pub(crate) max_live: AtomicUsize,
    /// The most placeholder bytes held at once, by all decoders together.
    pub(crate) max_placeholder: AtomicUsize,
    /// Generations started.
    pub(crate) generations: AtomicUsize,
}

/// What the decoders of one CRAM read from the FASTA, and the most they held at once: the
/// bound D109 sets on reference memory.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReferenceCounts {
    /// Sequences read from the FASTA.
    pub fetches: usize,
    /// The most sequences (and placeholders) held at once, by all decoders together.
    pub max_live: usize,
    /// The most bytes of placeholders held at once, by all decoders together (CRAMs stored
    /// without the reference, in slices of several references).
    pub max_placeholder_bytes: usize,
    /// Shared generations started (each one loads again what its slices need).
    pub generations: usize,
}

impl Counters {
    pub(crate) fn snapshot(&self) -> ReferenceCounts {
        ReferenceCounts {
            fetches: self.fetches.load(Ordering::Relaxed),
            max_live: self.max_live.load(Ordering::Relaxed),
            max_placeholder_bytes: self.max_placeholder.load(Ordering::Relaxed),
            generations: self.generations.load(Ordering::Relaxed),
        }
    }
}

/// What a repository and its [`RefAdapter`] share.
#[derive(Default)]
pub(crate) struct RefState {
    /// Fetch reference bases (a container's `RR`); without it, placeholders.
    pub(crate) required: AtomicBool,
    /// Size each placeholder to its reference's length (else empty: the `.crai`'s decoder
    /// never puts bases together).
    pub(crate) full_placeholders: AtomicBool,
    /// `N` bases added after a sequence, for reads that run past the end of their reference
    /// (htslib reads them against `N`), by name.
    pub(crate) pads: Mutex<HashMap<Vec<u8>, usize>>,
    /// Placeholder lengths beyond a reference's own, by name. noodles asks for the sequence of
    /// each mapped read of a slice of several references even when no base is read, and slices
    /// it at the read's positions, so the placeholder must reach them; it is zero bytes, so
    /// only noodles' copy of it takes memory.
    pub(crate) placeholders: Mutex<HashMap<Vec<u8>, usize>>,
    /// Bytes of placeholders the repository holds.
    pub(crate) placeholder_bytes: AtomicUsize,
    /// What the repository holds: each name and the length served for it.
    pub(crate) loaded: Mutex<HashMap<Vec<u8>, usize>>,
    /// The adapter's last failure, with gloci's message (noodles passes on only an `io::Error`).
    pub(crate) failure: Mutex<Option<IoError>>,
}

impl RefState {
    /// A state for a repository of reference bases (`required`) or of placeholders.
    pub(crate) fn new(required: bool) -> Arc<RefState> {
        let s = RefState::default();
        s.required.store(required, Ordering::Relaxed);
        s.full_placeholders.store(true, Ordering::Relaxed);
        Arc::new(s)
    }

    /// The length a repository with this state serves (or would serve) for `name`, whose
    /// reference is `len` long.
    pub(crate) fn length_of(&self, name: &[u8], len: usize) -> usize {
        if let Some(&l) = self.loaded.lock().ok().as_ref().and_then(|m| m.get(name)) {
            return l;
        }
        let extra = |m: &Mutex<HashMap<Vec<u8>, usize>>| {
            m.lock()
                .ok()
                .and_then(|m| m.get(name).copied())
                .unwrap_or(0)
        };
        if self.required.load(Ordering::Relaxed) {
            len + extra(&self.pads)
        } else if self.full_placeholders.load(Ordering::Relaxed) {
            len.max(extra(&self.placeholders))
        } else {
            0
        }
    }
}

/// The repository adapter over the `--reference` FASTA (D109). It never answers `None`.
pub(crate) struct RefAdapter {
    pub(crate) info: Arc<CramInfo>,
    pub(crate) fasta: Option<Mutex<IndexedFasta>>,
    pub(crate) state: Arc<RefState>,
    pub(crate) counters: Arc<Counters>,
}

pub(crate) const ADAPTER_FAILED: &str = "the reference adapter failed";

impl RefAdapter {
    /// An adapter with its own handle on the `--reference` FASTA.
    pub(crate) fn new(
        info: Arc<CramInfo>,
        state: Arc<RefState>,
        counters: Arc<Counters>,
    ) -> Result<RefAdapter, IoError> {
        let fasta = match &info.fasta {
            None => None,
            Some((f, _)) => {
                let f = f
                    .lock()
                    .map_err(|_| IoError::Internal("the reference lock is poisoned".into()))?;
                Some(Mutex::new(f.reopen()?))
            }
        };
        Ok(RefAdapter {
            info,
            fasta,
            state,
            counters,
        })
    }

    pub(crate) fn fetch(&mut self, name: &[u8]) -> Result<fasta::Record, IoError> {
        let record = |seq: Vec<u8>| {
            fasta::Record::new(
                fasta::record::Definition::new(name.to_vec(), None),
                fasta::record::Sequence::from(seq),
            )
        };
        let note = |state: &RefState, len: usize| {
            if let Ok(mut l) = state.loaded.lock() {
                l.insert(name.to_vec(), len);
            }
        };
        let ln = self
            .info
            .aln
            .ref_id(name)
            .map_or(0, |r| usize::try_from(self.info.aln.refs[r].1).unwrap_or(0));
        if !self.state.required.load(Ordering::Relaxed) {
            // A container stored without the reference reads no base of it.
            let len = self.state.length_of(name, ln);
            self.state
                .placeholder_bytes
                .fetch_add(len, Ordering::Relaxed);
            note(&self.state, len);
            return Ok(record(vec![0; len]));
        }
        let Some(f) = self.fasta.as_mut().and_then(|m| m.get_mut().ok()) else {
            return Err(self.info.needs_reference(name));
        };
        let Some(len) = f.len_of(name) else {
            return Err(self.info.needs_reference(name));
        };
        let pad = self
            .state
            .pads
            .lock()
            .ok()
            .and_then(|p| p.get(name).copied())
            .unwrap_or(0);
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
        self.counters.fetches.fetch_add(1, Ordering::Relaxed);
        note(&self.state, seq.len());
        Ok(record(seq))
    }
}

impl fasta::repository::Adapter for RefAdapter {
    fn get(&mut self, name: &[u8]) -> Option<io::Result<fasta::Record>> {
        Some(self.fetch(name).map_err(|e| {
            if let Ok(mut slot) = self.state.failure.lock()
                && slot.is_none()
            {
                *slot = Some(e);
            }
            io::Error::other(ADAPTER_FAILED)
        }))
    }
}

/// Load `name` into `repo` (whose adapter shares `state`) now, returning gloci's error.
fn preload(repo: &fasta::Repository, state: &RefState, name: &[u8]) -> Result<(), IoError> {
    match repo.get(name) {
        Some(Ok(_)) => Ok(()),
        Some(Err(e)) => Err(state
            .failure
            .lock()
            .ok()
            .and_then(|mut f| f.take())
            .unwrap_or_else(|| IoError::Internal(format!("loading a reference failed: {e}")))),
        None => Err(IoError::Internal(
            "the reference repository answered nothing".into(),
        )),
    }
}

/// How many references the slots hold at most.
const SLOTS: usize = 2;

/// One shared repository holding one reference.
struct Slot {
    rid: Option<usize>,
    repo: Option<fasta::Repository>,
    state: Arc<RefState>,
    /// Decoders using it now.
    users: usize,
    /// The container that used it last (least recently used goes first).
    last_use: u64,
}

/// One repository shared by every decoder, for slices that need several references, padding
/// past a reference's end, or placeholders: references are added as slices need them, each
/// loaded once; leases count which are in use.
pub(crate) struct Generation {
    pub(crate) repo: fasta::Repository,
    pub(crate) state: Arc<RefState>,
    /// Active leases per reference name.
    leased: Mutex<HashMap<Vec<u8>, usize>>,
    /// The names a decoder is loading now (changed under the cache's `generations` lock).
    loading: Mutex<HashSet<Vec<u8>>>,
}

/// A decoder's use of some references of a generation (or of a slot), given back on drop.
pub(crate) struct Lease<'a> {
    cache: &'a RefCache,
    /// The slot, or the generation and the names leased in it.
    slot: Option<usize>,
    generation: Option<(Arc<Generation>, Vec<Vec<u8>>)>,
    pub(crate) repo: fasta::Repository,
    pub(crate) state: Arc<RefState>,
}

impl Drop for Lease<'_> {
    fn drop(&mut self) {
        if let Some(i) = self.slot {
            self.cache.release_slot(i);
        }
        if let Some((g, names)) = &self.generation
            && let Ok(mut l) = g.leased.lock()
        {
            for n in names {
                if let Some(c) = l.get_mut(n) {
                    *c = c.saturating_sub(1);
                }
            }
        }
    }
}

impl Lease<'_> {
    /// Keep this lease's slot in use after the lease ends, until the pin drops: a reader of
    /// one chromosome keeps its reference between slices, so another worker needing a slot
    /// cannot evict it. `None` for a lease on a shared generation.
    pub(crate) fn pin(&self, cache: &Arc<RefCache>, rid: usize) -> Option<Pin> {
        let i = self.slot?;
        let mut slots = cache.slots.lock().ok()?;
        slots[i].users += 1;
        Some(Pin {
            cache: cache.clone(),
            slot: i,
            rid,
        })
    }
}

/// A slot kept in use by a reader (see [`Lease::pin`]).
pub(crate) struct Pin {
    cache: Arc<RefCache>,
    slot: usize,
    pub(crate) rid: usize,
}

impl Drop for Pin {
    fn drop(&mut self) {
        self.cache.release_slot(self.slot);
        // The reader is done with its chromosome, and no other reader of this file reads it
        // (a per-chromosome run gives each chromosome to one): an idle slot frees it now,
        // rather than holding it until another reference needs the slot.
        if let Ok(mut slots) = self.cache.slots.lock()
            && slots[self.slot].users == 0
            && slots[self.slot].rid == Some(self.rid)
        {
            slots[self.slot].rid = None;
            slots[self.slot].repo = None;
        }
        self.cache.freed.notify_all();
    }
}

/// The reference memory every decoder of one CRAM shares (D109).
///
/// - **Slots**: up to [`SLOTS`] repositories, each holding one reference, for slices whose
///   reads need one. A slot is reused only when no decoder is using it (the least recently used
///   first); a decoder that finds none free waits. So sorted input reads each reference from
///   the FASTA once and holds at most two at a time, whatever the thread count.
/// - **Generations**: one repository of reference bases and one of placeholders, shared by all
///   decoders, for slices that need several references (unsorted or name-sorted input), reads
///   past a reference's end, or a CRAM stored without the reference. A slice leases every
///   reference it needs at once; each is loaded once (marked under the cache's lock, loaded
///   outside it; a decoder that needs one being loaded waits for it, as noodles' repository
///   would let decoders that miss at once load it twice). When a slice needs a reference the
///   current generation lacks (or needs it longer) and the generation holds more bytes that
///   no slice needs or uses than the longest reference, a new generation starts; the old one is
///   freed when its last lease ends. So small references touched by one slice and not the next
///   stay, as in samtools, and the large ones are not read again.
pub(crate) struct RefCache {
    info: Arc<CramInfo>,
    slots: Mutex<Vec<Slot>>,
    freed: Condvar,
    /// The current generations of reference bases and of placeholders.
    generations: Mutex<[Option<Arc<Generation>>; 2]>,
    /// Every generation still alive (for the counts).
    alive: Mutex<Vec<Weak<Generation>>>,
    /// Signalled when a decoder has loaded what it marked as loading.
    loaded_cv: Condvar,
    /// The longest reference: the most bytes a generation keeps that no slice needs or uses
    /// before a slice that needs a reference it lacks starts a new one.
    longest: usize,
    pub(crate) counters: Arc<Counters>,
}

impl RefCache {
    pub(crate) fn new(info: Arc<CramInfo>, counters: Arc<Counters>) -> RefCache {
        RefCache::with_slots(info, counters, SLOTS)
    }

    /// A cache of `slots` references (at least [`SLOTS`]): a per-chromosome run keeps one per
    /// worker, so workers on different chromosomes do not evict each other's.
    pub(crate) fn with_slots(
        info: Arc<CramInfo>,
        counters: Arc<Counters>,
        slots: usize,
    ) -> RefCache {
        let longest = info
            .aln
            .refs
            .iter()
            .map(|(_, l)| usize::try_from(*l).unwrap_or(0))
            .max()
            .unwrap_or(0);
        let slots = (0..slots.max(SLOTS))
            .map(|_| Slot {
                rid: None,
                repo: None,
                state: RefState::new(true),
                users: 0,
                last_use: 0,
            })
            .collect();
        RefCache {
            info,
            slots: Mutex::new(slots),
            freed: Condvar::new(),
            generations: Mutex::new([None, None]),
            alive: Mutex::new(Vec::new()),
            loaded_cv: Condvar::new(),
            longest,
            counters,
        }
    }

    /// A decoder is done with slot `i`: when none uses it, an idle slot whose reference the
    /// shared generation also holds gives its copy up.
    fn release_slot(&self, i: usize) {
        let idle = self.slots.lock().ok().and_then(|mut slots| {
            slots[i].users -= 1;
            (slots[i].users == 0).then_some(slots[i].rid).flatten()
        });
        if let Some(rid) = idle {
            let name = &self.info.aln.refs[rid].0;
            let held = self.generations.lock().ok().is_some_and(|g| {
                g[0].as_ref()
                    .is_some_and(|g| g.state.loaded.lock().is_ok_and(|l| l.contains_key(name)))
            });
            if held
                && let Ok(mut slots) = self.slots.lock()
                && slots[i].rid == Some(rid)
                && slots[i].users == 0
            {
                slots[i].rid = None;
                slots[i].repo = None;
            }
        }
        self.freed.notify_all();
    }

    fn poisoned() -> IoError {
        IoError::Internal("the reference cache lock is poisoned".into())
    }

    /// The slot holding reference `rid`, for a slice of container `tick`; or the shared
    /// generation of reference bases when it holds `rid` already (no second copy).
    pub(crate) fn acquire(&self, rid: usize, tick: u64) -> Result<Lease<'_>, IoError> {
        let in_generation = self
            .generations
            .lock()
            .map_err(|_| Self::poisoned())?
            .first()
            .and_then(Option::as_ref)
            .and_then(|g| {
                g.state
                    .loaded
                    .lock()
                    .ok()
                    .map(|l| l.contains_key(&self.info.aln.refs[rid].0))
            })
            .unwrap_or(false);
        // A file not sorted by coordinate meets its references in any order: one shared
        // generation serves all its slices, so no reference is held twice.
        let sorted = self.info.aln.sort_order == SortOrder::Coordinate;
        if in_generation || !sorted {
            let len = usize::try_from(self.info.aln.refs[rid].1).unwrap_or(0);
            return self.lease(true, Some(&[(rid, len)]));
        }
        let mut slots = self.slots.lock().map_err(|_| Self::poisoned())?;
        loop {
            let found = slots.iter().position(|s| s.rid == Some(rid));
            let index = match found {
                Some(i) => Some(i),
                None => {
                    // An idle slot: an empty one first, then the least recently used.
                    let idle = (0..slots.len())
                        .filter(|&i| slots[i].users == 0)
                        .min_by_key(|&i| (slots[i].rid.is_some(), slots[i].last_use));
                    if let Some(i) = idle {
                        let state = RefState::new(true);
                        let adapter = RefAdapter::new(
                            self.info.clone(),
                            state.clone(),
                            self.counters.clone(),
                        )?;
                        let repo = fasta::Repository::new(adapter);
                        // Loaded now, under the lock: noodles' repository would have every
                        // decoder that misses at once load it again.
                        slots[i] = Slot {
                            rid: None,
                            repo: None,
                            state: RefState::new(true),
                            users: 0,
                            last_use: 0,
                        };
                        preload(&repo, &state, &self.info.aln.refs[rid].0)?;
                        slots[i] = Slot {
                            rid: Some(rid),
                            repo: Some(repo),
                            state,
                            users: 0,
                            last_use: tick,
                        };
                    }
                    idle
                }
            };
            if let Some(i) = index
                && let Some(repo) = slots[i].repo.clone()
            {
                slots[i].users += 1;
                slots[i].last_use = slots[i].last_use.max(tick);
                let state = slots[i].state.clone();
                drop(slots);
                self.note_live();
                return Ok(Lease {
                    cache: self,
                    slot: Some(i),
                    generation: None,
                    repo,
                    state,
                });
            }
            slots = self.freed.wait(slots).map_err(|_| Self::poisoned())?;
        }
    }

    fn new_generation(&self, required: bool) -> Result<Arc<Generation>, IoError> {
        self.counters.generations.fetch_add(1, Ordering::Relaxed);
        let state = RefState::new(required);
        let adapter = RefAdapter::new(self.info.clone(), state.clone(), self.counters.clone())?;
        let g = Arc::new(Generation {
            repo: fasta::Repository::new(adapter),
            state,
            leased: Mutex::new(HashMap::new()),
            loading: Mutex::new(HashSet::new()),
        });
        if let Ok(mut a) = self.alive.lock() {
            a.retain(|w| w.strong_count() > 0);
            a.push(Arc::downgrade(&g));
        }
        Ok(g)
    }

    /// The shared repository of reference bases (`required`) or placeholders for a slice that
    /// needs the references `need`, each at least as long as given (its length, plus padding;
    /// or a placeholder's reach). `None`: the references are not known before decoding, and
    /// noodles loads them as it meets them.
    pub(crate) fn lease(
        &self,
        required: bool,
        need: Option<&[(usize, usize)]>,
    ) -> Result<Lease<'_>, IoError> {
        let name = |rid: usize| &self.info.aln.refs[rid].0;
        let ln = |rid: usize| usize::try_from(self.info.aln.refs[rid].1).unwrap_or(0);
        let need = need.unwrap_or(&[]);
        let which = usize::from(!required);
        let names: Vec<Vec<u8>> = need.iter().map(|&(r, _)| name(r).clone()).collect();
        let mut gens = self.generations.lock().map_err(|_| Self::poisoned())?;
        // Under the lock: choose the generation, mark what this slice loads, lease it all.
        let (generation, to_load) = loop {
            let current = match &gens[which] {
                Some(g) => g.clone(),
                None => {
                    let g = self.new_generation(required)?;
                    gens[which] = Some(g.clone());
                    g
                }
            };
            // What the current generation lacks (and no decoder is loading), holds too short,
            // and the bytes it holds that no slice needs or uses.
            let (missing, short, idle) = {
                let loaded = current.state.loaded.lock().map_err(|_| Self::poisoned())?;
                let leased = current.leased.lock().map_err(|_| Self::poisoned())?;
                let loading = current.loading.lock().map_err(|_| Self::poisoned())?;
                let missing = need
                    .iter()
                    .any(|&(r, _)| !loaded.contains_key(name(r)) && !loading.contains(name(r)));
                let short = need
                    .iter()
                    .any(|&(r, len)| loaded.get(name(r)).is_some_and(|&have| have < len));
                let idle: usize = loaded
                    .iter()
                    .filter(|(n, _)| {
                        !names.contains(n) && leased.get(*n).copied().unwrap_or(0) == 0
                    })
                    .map(|(_, &len)| len)
                    .sum();
                (missing, short, idle)
            };
            let generation = if short || (missing && idle > self.longest) {
                // The old generation goes when its last lease ends (now, if none is left).
                drop(current);
                gens[which] = None;
                let g = self.new_generation(required)?;
                gens[which] = Some(g.clone());
                g
            } else {
                current
            };
            let loaded = generation
                .state
                .loaded
                .lock()
                .map_err(|_| Self::poisoned())?;
            let mut loading = generation.loading.lock().map_err(|_| Self::poisoned())?;
            if names.iter().any(|n| loading.contains(n)) {
                // Another decoder is loading one: wait for it, then look again.
                drop((loaded, loading));
                gens = self.loaded_cv.wait(gens).map_err(|_| Self::poisoned())?;
                continue;
            }
            let to_load: Vec<(usize, usize)> = need
                .iter()
                .copied()
                .filter(|&(r, _)| !loaded.contains_key(name(r)))
                .collect();
            for &(r, _) in &to_load {
                loading.insert(name(r).clone());
            }
            drop((loaded, loading));
            if let Ok(mut l) = generation.leased.lock() {
                for n in &names {
                    *l.entry(n.clone()).or_default() += 1;
                }
            }
            break (generation, to_load);
        };
        drop(gens);
        // Load outside the lock (other decoders go on); each reference once, as marked.
        let loaded = (|| -> Result<(), IoError> {
            for &(r, len) in &to_load {
                let n = name(r);
                let extra = |m: &Mutex<HashMap<Vec<u8>, usize>>, v: usize| {
                    if let Ok(mut m) = m.lock() {
                        m.insert(n.clone(), v);
                    }
                };
                if required {
                    extra(&generation.state.pads, len.saturating_sub(ln(r)));
                    // An idle slot holding it gives its copy up (the generation has its own).
                    if let Ok(mut slots) = self.slots.lock() {
                        for slot in slots.iter_mut() {
                            if slot.rid == Some(r) && slot.users == 0 {
                                slot.rid = None;
                                slot.repo = None;
                            }
                        }
                    }
                } else {
                    extra(&generation.state.placeholders, len);
                }
                preload(&generation.repo, &generation.state, n)?;
            }
            Ok(())
        })();
        if let Ok(_gens) = self.generations.lock()
            && let Ok(mut loading) = generation.loading.lock()
        {
            for &(r, _) in &to_load {
                loading.remove(name(r));
            }
        }
        self.loaded_cv.notify_all();
        let lease = Lease {
            cache: self,
            slot: None,
            repo: generation.repo.clone(),
            state: generation.state.clone(),
            generation: Some((generation, names)),
        };
        // On a failure the lease is given back as it drops.
        loaded?;
        self.note_live();
        Ok(lease)
    }

    /// Count what the slots and every generation alive hold, and keep the most.
    pub(crate) fn note_live(&self) {
        let slots: usize = self
            .slots
            .lock()
            .map(|s| {
                s.iter()
                    .filter_map(|s| s.repo.as_ref().map(|r| r.len()))
                    .sum()
            })
            .unwrap_or(0);
        let (mut gens, mut bytes) = (0, 0);
        if let Ok(alive) = self.alive.lock() {
            for g in alive.iter().filter_map(Weak::upgrade) {
                gens += g.repo.len();
                bytes += g.state.placeholder_bytes.load(Ordering::Relaxed);
            }
        }
        let c = &self.counters;
        c.max_live.fetch_max(slots + gens, Ordering::Relaxed);
        c.max_placeholder.fetch_max(bytes, Ordering::Relaxed);
    }
}
