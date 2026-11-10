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

    /// `inner` read through the cache.
    pub fn wrap<'a>(&'a self, inner: &'a dyn GenotypeSource) -> CachedGenotypes<'a> {
        CachedGenotypes { inner, cache: self }
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
            b"DP" => &mut out.dp,
            b"GQ" => &mut out.gq,
            b"FT" => &mut out.ft,
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
