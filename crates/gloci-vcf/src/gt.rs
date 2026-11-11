//! Genotype types and the VCF text decoder.

use memchr::{memchr, memchr_iter};
use std::cell::RefCell;

/// One sample's genotype. The buffer is reused across samples.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Genotype {
    /// One entry per allele (so its length is the ploidy); `None` = `.`.
    pub alleles: Vec<Option<u32>>,
    /// True when the genotype has at least two alleles and every separator
    /// is `|`.
    pub phased: bool,
}

/// The FORMAT fields `vcf filter` needs, for one sample.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SampleFields {
    /// `None` when GT is absent from FORMAT or from a short sample column.
    /// Read-only for callers: decoders recycle the allele buffer through a
    /// private spare, so build values with `Default` and assign fields only
    /// to test fakes.
    pub gt: Option<Genotype>,
    /// Allele buffer parked while `gt` is `None`, so it is never reallocated.
    spare: Vec<Option<u32>>,
    pub dp: Option<u32>,
    pub gq: Option<u32>,
    /// FT is `PASS`, `.` or absent (D50).
    pub ft_pass: bool,
    /// Set by [`crate::MaskedGenotypes`] when masking removed this sample's
    /// genotype, so a statistic can tell a filtered genotype from a missing
    /// one. Decoders leave it false.
    pub masked: bool,
}

impl SampleFields {
    pub(crate) fn reset(&mut self) {
        // Keep the genotype buffer: `gt` is cleared in place when GT is
        // present and dropped only by `clear_gt`.
        self.dp = None;
        self.gq = None;
        self.ft_pass = true;
        self.masked = false;
    }

    /// Overwrite every field with `o`'s, reusing the allele buffer.
    pub(crate) fn copy_from(&mut self, o: &SampleFields) {
        self.dp = o.dp;
        self.gq = o.gq;
        self.ft_pass = o.ft_pass;
        self.masked = o.masked;
        match &o.gt {
            Some(g) => {
                let b = self.gt_buf();
                b.alleles.extend_from_slice(&g.alleles);
                b.phased = g.phased;
            }
            None => self.clear_gt(),
        }
    }

    /// Mark GT absent, parking the allele buffer.
    pub(crate) fn clear_gt(&mut self) {
        if let Some(g) = self.gt.take() {
            self.spare = g.alleles;
        }
    }

    /// The genotype buffer, cleared, ready to be filled.
    pub(crate) fn gt_buf(&mut self) -> &mut Genotype {
        if self.gt.is_none() {
            self.gt = Some(Genotype {
                alleles: std::mem::take(&mut self.spare),
                phased: false,
            });
        }
        let g = self.gt.as_mut().expect("just set");
        g.alleles.clear();
        g.phased = false;
        g
    }
}

/// Random access to the samples of one record.
pub trait GenotypeSource {
    fn n_samples(&self) -> usize;
    /// Decode sample `i` into `out`, overwriting every field. Panics when
    /// `i >= n_samples()`; a malformed sample is an `Err`.
    fn sample(&self, i: usize, out: &mut SampleFields) -> Result<(), String>;
}

/// Decoded samples of the current record, shared between the stages that
/// read them: the first stage to ask for a sample decodes it, the others copy
/// the result. Buffers are reused from record to record.
#[derive(Default)]
pub struct SampleCache {
    slots: RefCell<Vec<SampleFields>>,
    valid: RefCell<Vec<bool>>,
}

impl SampleCache {
    /// Forget the previous record's samples; `n` is the new record's sample count.
    pub fn start_record(&self, n: usize) {
        let mut valid = self.valid.borrow_mut();
        valid.clear();
        valid.resize(n, false);
        let mut slots = self.slots.borrow_mut();
        if slots.len() < n {
            slots.resize_with(n, SampleFields::default);
        }
    }

    /// The cached samples alone, for a consumer that reads only samples the
    /// record's first stage already decoded (a sample not decoded is an error).
    pub fn reader(&self) -> CacheReader<'_> {
        CacheReader { cache: self }
    }

    /// `inner` read through the cache.
    pub fn wrap<'a>(&'a self, inner: &'a dyn GenotypeSource) -> CachedGenotypes<'a> {
        CachedGenotypes { inner, cache: self }
    }
}

/// The samples [`SampleCache`] holds for the current record.
pub struct CacheReader<'a> {
    cache: &'a SampleCache,
}

impl GenotypeSource for CacheReader<'_> {
    fn n_samples(&self) -> usize {
        self.cache.valid.borrow().len()
    }

    fn sample(&self, i: usize, out: &mut SampleFields) -> Result<(), String> {
        if !self.cache.valid.borrow().get(i).copied().unwrap_or(false) {
            return Err(format!(
                "internal error: sample {} was not decoded by the selection stage",
                i + 1
            ));
        }
        out.copy_from(&self.cache.slots.borrow()[i]);
        Ok(())
    }
}

/// A [`GenotypeSource`] that decodes each sample at most once per record
/// (since the last [`SampleCache::start_record`]).
pub struct CachedGenotypes<'a> {
    inner: &'a dyn GenotypeSource,
    cache: &'a SampleCache,
}

impl GenotypeSource for CachedGenotypes<'_> {
    fn n_samples(&self) -> usize {
        self.inner.n_samples()
    }

    fn sample(&self, i: usize, out: &mut SampleFields) -> Result<(), String> {
        let mut valid = self.cache.valid.borrow_mut();
        if i >= valid.len() {
            // The record was never started in the cache: decode directly.
            return self.inner.sample(i, out);
        }
        let mut slots = self.cache.slots.borrow_mut();
        if !valid[i] {
            self.inner.sample(i, &mut slots[i])?;
            valid[i] = true;
        }
        out.copy_from(&slots[i]);
        Ok(())
    }
}

pub(crate) fn ft_passes(v: &[u8]) -> bool {
    v.is_empty() || v == b"." || v == b"PASS"
}

pub(crate) fn parse_uint(v: &[u8], what: &str) -> Result<Option<u32>, String> {
    if v == b"." || v.is_empty() {
        return Ok(None);
    }
    let mut n: u32 = 0;
    for &b in v {
        if !b.is_ascii_digit() {
            return Err(bad(what, v));
        }
        n = n
            .checked_mul(10)
            .and_then(|n| n.checked_add(u32::from(b - b'0')))
            .ok_or_else(|| bad(what, v))?;
    }
    Ok(Some(n))
}

fn bad(what: &str, v: &[u8]) -> String {
    format!(
        "{what} value {:?} is not a non-negative integer",
        String::from_utf8_lossy(v)
    )
}

/// Parse a text GT value into `g` (already cleared). A VCF 4.4 leading
/// phase indicator (`|0|1`, `/0|1`) states the first allele's phasing and
/// counts as a separator for `phased`.
fn parse_gt(v: &[u8], g: &mut Genotype) -> Result<(), String> {
    let mut seps = 0usize;
    let mut phased_seps = 0usize;
    let mut start = 0;
    if let Some(&p @ (b'|' | b'/')) = v.first()
        && v.len() > 1
    {
        seps += 1;
        phased_seps += usize::from(p == b'|');
        start = 1;
    }
    let mut i = start;
    loop {
        let at_end = i == v.len();
        if at_end || v[i] == b'/' || v[i] == b'|' {
            g.alleles.push(parse_allele(&v[start..i], v)?);
            if at_end {
                break;
            }
            seps += 1;
            phased_seps += usize::from(v[i] == b'|');
            start = i + 1;
        }
        i += 1;
    }
    g.phased = g.alleles.len() > 1 && seps == phased_seps;
    Ok(())
}

fn parse_allele(a: &[u8], whole: &[u8]) -> Result<Option<u32>, String> {
    if a.is_empty() {
        return Err(format!(
            "GT value {:?} has an empty allele",
            String::from_utf8_lossy(whole)
        ));
    }
    parse_uint(a, "GT allele")
}

/// Which of the optional FORMAT fields a decoder reads, besides GT (always
/// read). A field left out stays `None` (FT: passing) in every [`SampleFields`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GtFields {
    pub dp: bool,
    pub gq: bool,
    pub ft: bool,
}

impl GtFields {
    /// DP, GQ and FT.
    pub const ALL: GtFields = GtFields {
        dp: true,
        gq: true,
        ft: true,
    };
    /// GT only.
    pub const GT_ONLY: GtFields = GtFields {
        dp: false,
        gq: false,
        ft: false,
    };
}

/// A VCF data line with the FORMAT key positions located once.
pub struct TextGenotypes<'a> {
    line: &'a [u8],
    gt: Option<usize>,
    dp: Option<usize>,
    gq: Option<usize>,
    ft: Option<usize>,
    /// Highest subfield index any wanted key has.
    last: Option<usize>,
    /// Start and end of every sample column within `line`.
    samples: Vec<(u32, u32)>,
}

/// Index a VCF data line (no trailing newline needed). Columns 1-9 must be
/// present; a line with no FORMAT column has zero samples.
pub fn text_genotypes(line: &[u8]) -> Result<TextGenotypes<'_>, String> {
    text_genotypes_with(line, GtFields::ALL)
}

/// [`text_genotypes`] reading only GT and the fields `fields` names.
pub fn text_genotypes_with(line: &[u8], fields: GtFields) -> Result<TextGenotypes<'_>, String> {
    let line = line.strip_suffix(b"\n").unwrap_or(line);
    let line = line.strip_suffix(b"\r").unwrap_or(line);
    if u32::try_from(line.len()).is_err() {
        return Err("a VCF line over 4 GiB is not supported".into());
    }
    let mut tabs = memchr_iter(b'\t', line);
    // Seven tabs end columns 1-7; the eighth ends INFO; the ninth (if any) ends FORMAT.
    for k in 0..7 {
        if tabs.next().is_none() {
            return Err(format!(
                "the line has {} columns; a VCF line needs at least 8",
                k + 1
            ));
        }
    }
    let info_end = tabs.next();
    let mut out = TextGenotypes {
        line,
        gt: None,
        dp: None,
        gq: None,
        ft: None,
        last: None,
        samples: Vec::new(),
    };
    let Some(info_end) = info_end else {
        return Ok(out); // 8 columns, no FORMAT, no samples
    };
    let fmt_start = info_end + 1;
    let fmt_end = tabs.next().unwrap_or(line.len());
    for (idx, key) in line[fmt_start..fmt_end].split(|&b| b == b':').enumerate() {
        let slot = match key {
            b"GT" => &mut out.gt,
            b"DP" if fields.dp => &mut out.dp,
            b"GQ" if fields.gq => &mut out.gq,
            b"FT" if fields.ft => &mut out.ft,
            _ => continue,
        };
        slot.get_or_insert(idx);
    }
    out.last = [out.gt, out.dp, out.gq, out.ft].into_iter().flatten().max();
    if fmt_end < line.len() {
        let mut start = fmt_end + 1;
        for t in tabs {
            out.samples.push((start as u32, t as u32));
            start = t + 1;
        }
        out.samples.push((start as u32, line.len() as u32));
    }
    Ok(out)
}

impl GenotypeSource for TextGenotypes<'_> {
    fn n_samples(&self) -> usize {
        self.samples.len()
    }

    fn sample(&self, i: usize, out: &mut SampleFields) -> Result<(), String> {
        let (s, e) = self.samples[i];
        let col = &self.line[s as usize..e as usize];
        out.reset();
        out.clear_gt();
        // Walk the colon-separated subfields in place, stopping once the
        // highest wanted index has been seen.
        let Some(last) = self.last else {
            return Ok(());
        };
        let mut rest = col;
        let mut idx = 0;
        loop {
            let (field, more) = match memchr(b':', rest) {
                Some(p) => (&rest[..p], Some(&rest[p + 1..])),
                None => (rest, None),
            };
            if Some(idx) == self.gt {
                if !field.is_empty() {
                    // A bare `.` is a missing sample: one missing allele.
                    parse_gt(field, out.gt_buf()).map_err(|e| format!("sample {}: {e}", i + 1))?;
                }
            } else if Some(idx) == self.dp {
                out.dp = parse_uint(field, "DP").map_err(|e| format!("sample {}: {e}", i + 1))?;
            } else if Some(idx) == self.gq {
                out.gq = parse_uint(field, "GQ").map_err(|e| format!("sample {}: {e}", i + 1))?;
            } else if Some(idx) == self.ft {
                out.ft_pass = ft_passes(field);
            }
            match more {
                Some(m) if idx < last => {
                    rest = m;
                    idx += 1;
                }
                _ => break,
            }
        }
        Ok(())
    }
}

/// Where a record's samples sit in a [`PackedGenotypes`] arena.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RecordSlot {
    start: usize,
    n: usize,
}

const F_GT: u16 = 1;
const F_PHASED: u16 = 2;
const F_DP: u16 = 4;
const F_GQ: u16 = 8;
const F_FT_FAIL: u16 = 16;
const F_MASKED: u16 = 32;

/// The allele value that stands for a missing allele in the arena.
const NO_ALLELE: u32 = u32::MAX;

/// One packed sample: 16 bytes, plus 4 bytes per allele in the arena pool.
#[derive(Debug, Clone, Copy)]
struct Packed {
    dp: u32,
    gq: u32,
    /// Offset of the sample's alleles in the pool.
    off: u32,
    /// Number of alleles.
    len: u16,
    flags: u16,
}

/// An arena holding the decoded samples of many records, so a worker thread
/// can hand them to the thread that consumes the records. A diploid sample
/// costs 24 bytes. Masked samples are stored as the masking leaves them (no
/// GT, DP or GQ, `masked` set); samples a `keep` list leaves out are stored
/// as missing.
#[derive(Debug, Default)]
pub struct PackedGenotypes {
    samples: Vec<Packed>,
    alleles: Vec<u32>,
}

impl PackedGenotypes {
    /// Forget every record, keeping the buffers.
    pub fn clear(&mut self) {
        self.samples.clear();
        self.alleles.clear();
    }

    /// Samples held (across all records).
    pub fn n_samples_held(&self) -> usize {
        self.samples.len()
    }

    /// Decode every sample of `src` into the arena as one record. `masked`
    /// (empty, or one flag per sample) marks the samples the mask rules
    /// removed; `fields` names the optional fields to keep.
    pub fn push_record(
        &mut self,
        src: &dyn GenotypeSource,
        fields: GtFields,
        masked: &[bool],
    ) -> Result<RecordSlot, String> {
        self.push_record_selected(src, fields, masked, None)
    }

    /// [`Self::push_record`] decoding only the samples `keep` marks (`None` =
    /// all); the others are stored as missing and never read from `src`.
    pub fn push_record_selected(
        &mut self,
        src: &dyn GenotypeSource,
        fields: GtFields,
        masked: &[bool],
        keep: Option<&[bool]>,
    ) -> Result<RecordSlot, String> {
        let n = src.n_samples();
        let start = self.samples.len();
        let mark = (start, self.alleles.len());
        let mut f = SampleFields::default();
        for i in 0..n {
            let r = if keep.is_some_and(|k| !k.get(i).copied().unwrap_or(false)) {
                self.samples.push(Packed {
                    dp: 0,
                    gq: 0,
                    off: 0,
                    len: 0,
                    flags: 0,
                });
                Ok(())
            } else {
                src.sample(i, &mut f).and_then(|()| {
                    self.push_sample(&f, fields, masked.get(i).copied().unwrap_or(false))
                })
            };
            if let Err(e) = r {
                self.samples.truncate(mark.0);
                self.alleles.truncate(mark.1);
                return Err(e);
            }
        }
        Ok(RecordSlot { start, n })
    }

    fn push_sample(
        &mut self,
        f: &SampleFields,
        fields: GtFields,
        masked: bool,
    ) -> Result<(), String> {
        let off = u32::try_from(self.alleles.len())
            .map_err(|_| "a batch of genotypes is too large to pack".to_string())?;
        let mut p = Packed {
            dp: 0,
            gq: 0,
            off,
            len: 0,
            flags: 0,
        };
        if masked {
            p.flags |= F_MASKED;
        } else {
            if let Some(g) = &f.gt {
                p.len = u16::try_from(g.alleles.len()).map_err(|_| {
                    format!(
                        "a genotype with {} alleles is not supported",
                        g.alleles.len()
                    )
                })?;
                p.flags |= F_GT;
                if g.phased {
                    p.flags |= F_PHASED;
                }
                for a in &g.alleles {
                    self.alleles.push(match a {
                        Some(NO_ALLELE) => {
                            return Err(format!("genotype allele {NO_ALLELE} is too large"));
                        }
                        Some(a) => *a,
                        None => NO_ALLELE,
                    });
                }
            }
            if fields.dp
                && let Some(dp) = f.dp
            {
                p.dp = dp;
                p.flags |= F_DP;
            }
            if fields.gq
                && let Some(gq) = f.gq
            {
                p.gq = gq;
                p.flags |= F_GQ;
            }
        }
        if fields.ft && !f.ft_pass {
            p.flags |= F_FT_FAIL;
        }
        self.samples.push(p);
        Ok(())
    }

    /// The record at `slot`, which must come from this arena (and not predate
    /// a [`Self::clear`]).
    pub fn record(&self, slot: RecordSlot) -> PackedRecord<'_> {
        PackedRecord { arena: self, slot }
    }
}

/// One record of a [`PackedGenotypes`] arena.
pub struct PackedRecord<'a> {
    arena: &'a PackedGenotypes,
    slot: RecordSlot,
}

impl GenotypeSource for PackedRecord<'_> {
    fn n_samples(&self) -> usize {
        self.slot.n
    }

    fn sample(&self, i: usize, out: &mut SampleFields) -> Result<(), String> {
        let Some(p) = (i < self.slot.n)
            .then(|| self.arena.samples.get(self.slot.start + i))
            .flatten()
        else {
            return Err(format!(
                "internal error: sample {} is outside the packed record",
                i + 1
            ));
        };
        out.reset();
        out.clear_gt();
        out.masked = p.flags & F_MASKED != 0;
        out.ft_pass = p.flags & F_FT_FAIL == 0;
        if p.flags & F_DP != 0 {
            out.dp = Some(p.dp);
        }
        if p.flags & F_GQ != 0 {
            out.gq = Some(p.gq);
        }
        if p.flags & F_GT != 0 {
            let off = p.off as usize;
            let alleles = self
                .arena
                .alleles
                .get(off..off + usize::from(p.len))
                .ok_or_else(|| "internal error: a packed genotype is truncated".to_string())?;
            let g = out.gt_buf();
            g.alleles
                .extend(alleles.iter().map(|&a| (a != NO_ALLELE).then_some(a)));
            g.phased = p.flags & F_PHASED != 0;
        }
        Ok(())
    }
}
