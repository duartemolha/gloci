//! Site views, the static and statistics filters, and `--thin`.

use super::{SiteFilters, VcfType};
use crate::stats::SiteStats;
use gloci_core::{ChromId, Interval, VariantKind, alt_kind};
use gloci_io::BcfDicts;
use gloci_io::bcf::typed::{self, Int};
use memchr::memchr;
use std::fmt::Write as _;

/// The columns of one record the site filters read. For a BCF record,
/// `info` holds the INFO keys only (`DB;DP`), `id` and `filter` are `.`
/// when empty; everything else matches the text record.
#[derive(Clone, Copy, Debug)]
pub struct SiteView<'a> {
    pub chrom: &'a [u8],
    /// 1-based.
    pub pos: u64,
    pub id: &'a [u8],
    pub ref_: &'a [u8],
    /// The raw ALT column (comma-separated).
    pub alts: &'a [u8],
    pub qual: Option<f64>,
    pub filter: &'a [u8],
    pub info: &'a [u8],
    /// The record's footprint (`RecordMeta.span`), used by
    /// `exclude_regions`.
    pub span: Interval,
}

fn parse_qual(f: &[u8]) -> Result<Option<f64>, String> {
    if f == b"." || f.is_empty() {
        return Ok(None);
    }
    std::str::from_utf8(f)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .map(Some)
        .ok_or_else(|| format!("QUAL {:?} is not a number", String::from_utf8_lossy(f)))
}

impl<'a> SiteView<'a> {
    /// The view of a VCF text line (without the newline). Borrows `line`.
    pub fn from_text(line: &'a [u8], span: Interval) -> Result<SiteView<'a>, String> {
        let mut cols: [&[u8]; 8] = [b""; 8];
        let mut rest = line;
        for (i, col) in cols.iter_mut().enumerate() {
            match memchr(b'\t', rest) {
                Some(p) => {
                    *col = &rest[..p];
                    rest = &rest[p + 1..];
                }
                None if i == 7 => *col = rest,
                None => {
                    return Err(format!(
                        "the line has {} columns; a VCF record needs at least 8",
                        i + 1
                    ));
                }
            }
        }
        let pos = std::str::from_utf8(cols[1])
            .ok()
            .and_then(|s| s.parse::<u64>().ok())
            .ok_or_else(|| {
                format!(
                    "POS {:?} is not an integer",
                    String::from_utf8_lossy(cols[1])
                )
            })?;
        Ok(SiteView {
            chrom: cols[0],
            pos,
            id: cols[2],
            ref_: cols[3],
            alts: cols[4],
            qual: parse_qual(cols[5])?,
            filter: cols[6],
            info: cols[7],
            span,
        })
    }
}

/// Reusable buffers that decode the shared part of a framed BCF record
/// (`l_shared`, `l_indiv`, shared, individual) into a [`SiteView`].
#[derive(Default)]
pub struct BcfSite {
    alts: Vec<u8>,
    filter: Vec<u8>,
    info: Vec<u8>,
    qual: String,
}

const SHORT: &str = "the record is shorter than its header says";

impl BcfSite {
    pub fn view<'a>(
        &'a mut self,
        raw: &'a [u8],
        dicts: &'a BcfDicts,
        span: Interval,
    ) -> Result<SiteView<'a>, String> {
        if raw.len() < 8 {
            return Err(SHORT.into());
        }
        let l_shared = u32::from_le_bytes([raw[0], raw[1], raw[2], raw[3]]) as usize;
        let shared = raw.get(8..8 + l_shared).ok_or(SHORT)?;
        if shared.len() < 24 {
            return Err(format!(
                "the shared section is {} bytes; a BCF record needs at least 24",
                shared.len()
            ));
        }
        let word = |o: usize| [shared[o], shared[o + 1], shared[o + 2], shared[o + 3]];
        let contig = i32::from_le_bytes(word(0));
        let pos = match i32::from_le_bytes(word(4)) {
            p if p < 0 => return Err(format!("POS {} is not valid", i64::from(p) + 1)),
            p => p as u64 + 1,
        };
        let chrom = usize::try_from(contig)
            .ok()
            .and_then(|i| dicts.contigs.get(i))
            .ok_or_else(|| format!("contig index {contig} is not in the header"))?;
        // Go through the shortest decimal form so that 0.7f32 compares like
        // the text value 0.7.
        let qual = match u32::from_le_bytes(word(12)) {
            0x7f80_0001 => None,
            bits => {
                self.qual.clear();
                write!(self.qual, "{}", f32::from_bits(bits)).expect("write to String");
                self.qual.parse::<f64>().ok()
            }
        };
        let n_allele_info = u32::from_le_bytes(word(16));
        let (n_allele, n_info) = (n_allele_info >> 16, n_allele_info & 0xffff);

        let mut cur = &shared[24..];
        let id = typed::read(&mut cur)?
            .chars()
            .ok_or("the ID is not a string")?;
        let id = if id.is_empty() { &b"."[..] } else { id };
        let mut ref_: &[u8] = b"";
        self.alts.clear();
        for i in 0..n_allele {
            let a = typed::read(&mut cur)?
                .chars()
                .ok_or("an allele is not a character string")?;
            match i {
                0 => ref_ = a,
                1 => self.alts.extend_from_slice(a),
                _ => {
                    self.alts.push(b',');
                    self.alts.extend_from_slice(a);
                }
            }
        }
        if n_allele <= 1 {
            self.alts.push(b'.');
        }

        let name = |key: i32| -> Result<&'a [u8], String> {
            usize::try_from(key)
                .ok()
                .and_then(|i| dicts.strings.get(i))
                .map(Vec::as_slice)
                .ok_or_else(|| format!("dictionary index {key} is not in the header"))
        };
        self.filter.clear();
        let f = typed::read(&mut cur)?;
        for i in 0..f.n {
            if let Some(Int::Value(k)) = f.int(i) {
                if !self.filter.is_empty() {
                    self.filter.push(b';');
                }
                self.filter.extend_from_slice(name(k)?);
            }
        }
        if self.filter.is_empty() {
            self.filter.push(b'.');
        }
        self.info.clear();
        for _ in 0..n_info {
            let k = typed::read(&mut cur)?;
            let key = match k.int(0) {
                Some(Int::Value(key)) if k.n == 1 => key,
                _ => return Err("an INFO key is not a single integer".into()),
            };
            typed::read(&mut cur)?;
            if !self.info.is_empty() {
                self.info.push(b';');
            }
            self.info.extend_from_slice(name(key)?);
        }
        if self.info.is_empty() {
            self.info.push(b'.');
        }
        Ok(SiteView {
            chrom,
            pos,
            id,
            ref_,
            alts: &self.alts,
            qual,
            filter: &self.filter,
            info: &self.info,
            span,
        })
    }
}

const SNV: u8 = 1;
const MNV: u8 = 2;
const INS: u8 = 4;
const DEL: u8 = 8;
const SV: u8 = 16;
const REF: u8 = 32;
const OTHER: u8 = 64;

impl VcfType {
    fn mask(self) -> u8 {
        match self {
            VcfType::Snv => SNV,
            VcfType::Mnv => MNV,
            VcfType::Ins => INS,
            VcfType::Del => DEL,
            VcfType::Indel => INS | DEL,
            VcfType::Sv => SV,
            VcfType::Ref => REF,
            VcfType::Other => OTHER,
        }
    }
}

fn mask_of(list: &[VcfType]) -> u8 {
    list.iter().fold(0, |m, t| m | t.mask())
}

fn kind_bit(k: VariantKind) -> u8 {
    match k {
        VariantKind::Snv => SNV,
        VariantKind::Mnv => MNV,
        VariantKind::Ins => INS,
        VariantKind::Del => DEL,
        VariantKind::Sv => SV,
        VariantKind::Ref => REF,
        VariantKind::Other => OTHER,
    }
}

fn type_set(s: &SiteView) -> u8 {
    s.alts.split(|&b| b == b',').fold(0, |m, alt| {
        m | kind_bit(alt_kind(s.pos.saturating_sub(1), s.ref_, alt))
    })
}

fn n_alleles(alts: &[u8]) -> usize {
    if alts.is_empty() || alts == b"." {
        1
    } else {
        2 + alts.iter().filter(|&&b| b == b',').count()
    }
}

fn info_has(info: &[u8], key: &[u8]) -> bool {
    info.split(|&b| b == b';')
        .any(|item| item.split(|&b| b == b'=').next() == Some(key))
}

fn tags(filter: &[u8]) -> impl Iterator<Item = &[u8]> {
    filter.split(|&b| b == b';')
}

fn has_any(filter: &[u8], list: &[Vec<u8>]) -> bool {
    tags(filter).any(|t| list.iter().any(|l| l == t))
}

impl SiteFilters {
    /// Do the filters need the genotypes (the statistics)?
    pub fn needs_genotypes(&self) -> bool {
        self.min_maf.is_some()
            || self.max_maf.is_some()
            || self.min_mac.is_some()
            || self.max_mac.is_some()
            || self.max_missing.is_some()
            || self.max_missing_count.is_some()
            || self.min_mean_dp.is_some()
            || self.max_mean_dp.is_some()
            || self.hwe.is_some()
            || self.phased_only
    }

    /// The `--thin` state, or `None` without `--thin`. Run it after every
    /// other filter has passed.
    pub fn thinner(&self) -> Option<Thinner> {
        self.thin.map(Thinner::new)
    }

    /// The filters that need only the record's columns. Allocation free.
    /// `pass_only` and `keep_filter` both apply, so `--pass-only --keep-filter
    /// PASS` drops FILTER `.` records, which `--pass-only` alone keeps.
    pub fn static_pass(&self, s: &SiteView) -> bool {
        if let Some(c) = &self.chroms
            && !c.contains(s.chrom)
        {
            return false;
        }
        if self.exclude_chroms.contains(s.chrom) {
            return false;
        }
        if let Some(p) = &self.positions
            && !p.contains(s.chrom, s.pos)
        {
            return false;
        }
        if let Some(p) = &self.exclude_positions
            && p.contains(s.chrom, s.pos)
        {
            return false;
        }
        if self.ids.is_some() || self.exclude_ids.is_some() {
            let has = |set: &std::collections::HashSet<Vec<u8>>| {
                s.id.split(|&b| b == b';')
                    .any(|i| i != b"." && set.contains(i))
            };
            if self.ids.as_ref().is_some_and(|set| !has(set)) {
                return false;
            }
            if self.exclude_ids.as_ref().is_some_and(has) {
                return false;
            }
        }
        if self.types.is_some() || self.exclude_types.is_some() {
            let set = type_set(s);
            if let Some(t) = &self.types
                && set & mask_of(t) == 0
            {
                return false;
            }
            if let Some(t) = &self.exclude_types
                && set & mask_of(t) != 0
            {
                return false;
            }
        }
        if self.min_alleles.is_some() || self.max_alleles.is_some() {
            let n = n_alleles(s.alts);
            if self.min_alleles.is_some_and(|m| n < m) || self.max_alleles.is_some_and(|m| n > m) {
                return false;
            }
        }
        if let Some(min) = self.min_qual
            && !s.qual.is_some_and(|q| q >= min)
        {
            return false;
        }
        if self.pass_only && !matches!(s.filter, b"PASS" | b"." | b"") {
            return false;
        }
        if let Some(k) = &self.keep_filter
            && !has_any(s.filter, k)
        {
            return false;
        }
        if let Some(r) = &self.remove_filter
            && has_any(s.filter, r)
        {
            return false;
        }
        if self.require_info.iter().any(|k| !info_has(s.info, k))
            || self.exclude_info.iter().any(|k| info_has(s.info, k))
        {
            return false;
        }
        if let Some(regions) = &self.exclude_regions
            && regions
                .get(s.chrom)
                .is_some_and(|r| r.matches(s.span, self.inclusive_points))
        {
            return false;
        }
        true
    }

    /// The filters over the recomputed statistics.
    pub fn stats_pass(&self, st: &SiteStats) -> bool {
        let ok = |v: Option<f64>, f: &dyn Fn(f64) -> bool| v.is_some_and(f);
        if let Some(m) = self.min_maf
            && !ok(st.maf(), &|v| v >= m)
        {
            return false;
        }
        if let Some(m) = self.max_maf
            && !ok(st.maf(), &|v| v <= m)
        {
            return false;
        }
        if let Some(m) = self.min_mac
            && st.mac().is_none_or(|v| v < m)
        {
            return false;
        }
        if let Some(m) = self.max_mac
            && st.mac().is_none_or(|v| v > m)
        {
            return false;
        }
        if self.max_missing.is_some_and(|m| st.missing_fraction() > m)
            || self.max_missing_count.is_some_and(|m| st.missing > m)
        {
            return false;
        }
        if let Some(m) = self.min_mean_dp
            && !ok(st.mean_dp(), &|v| v >= m)
        {
            return false;
        }
        if let Some(m) = self.max_mean_dp
            && !ok(st.mean_dp(), &|v| v <= m)
        {
            return false;
        }
        if let Some(p) = self.hwe
            && st.hwe_p().is_some_and(|v| v < p)
        {
            return false;
        }
        !self.phased_only || st.all_phased
    }
}

/// `--thin N`: keep a site when it is the first on its chromosome, or at
/// least N bases after the last kept one. The input must be sorted (the
/// caller checks).
#[derive(Clone, Debug, Default)]
pub struct Thinner {
    n: u64,
    last: Vec<Option<u64>>,
}

impl Thinner {
    pub fn new(n: u64) -> Thinner {
        Thinner {
            n,
            last: Vec::new(),
        }
    }

    pub fn keep(&mut self, chrom: ChromId, pos: u64) -> bool {
        let i = chrom.0 as usize;
        if self.last.len() <= i {
            self.last.resize(i + 1, None);
        }
        let slot = &mut self.last[i];
        if slot.is_some_and(|last| pos.saturating_sub(last) < self.n) {
            return false;
        }
        *slot = Some(pos);
        true
    }
}
