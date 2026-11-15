//! The whole-file statistics. Each is bounded by the number of distinct keys,
//! never by the number of records: `FilterSummary` and `IndelHist` keep a
//! map keyed by FILTER string or indel length and write in `finish`;
//! `SnpDensity` streams its rows and keeps only the current bin.

use super::engine::{Out, Stat, StatSink};
use super::sorted::SortedInput;
use crate::filter::SiteView;
use crate::gt::GenotypeSource;
use crate::stats::SiteStats;
use gloci_core::cfmt::{write_g6, write_i64, write_u64};
use std::collections::{BTreeMap, HashMap};
use std::io;

/// Substitution models in output order; the index is `pair_index`.
const MODELS: [&str; 6] = ["AC", "AG", "AT", "CG", "CT", "GT"];

fn code(b: u8) -> Option<u8> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// For a single-base REF and ALT over ACGT (either case): the model index
/// and whether it is a transition.
fn substitution(ref_: &[u8], alt: &[u8]) -> Option<(usize, bool)> {
    let (&[r], &[a]) = (ref_, alt) else {
        return None;
    };
    let (r, a) = (code(r)?, code(a)?);
    if r == a {
        return None;
    }
    let (lo, hi) = (r.min(a), r.max(a));
    let idx = match (lo, hi) {
        (0, 1) => 0,
        (0, 2) => 1,
        (0, 3) => 2,
        (1, 2) => 3,
        (1, 3) => 4,
        _ => 5,
    };
    Some((idx, matches!(idx, 1 | 4)))
}

fn first_alt(alts: &[u8]) -> &[u8] {
    alts.split(|&c| c == b',').next().unwrap_or(b"")
}

/// The row key of a FILTER: its `;`-separated tags in ascending byte order
/// (vcftools sorts them, so `b;a` and `a;b` are one row; BCF joins tags in
/// dictionary order, which this normalises too). Allocation-free once the
/// scratch buffers have grown.
fn normalise<'a>(ranges: &mut Vec<(u32, u32)>, key: &'a mut Vec<u8>, filter: &'a [u8]) -> &'a [u8] {
    if !filter.contains(&b';') {
        return filter;
    }
    ranges.clear();
    let mut start = 0usize;
    for (i, &c) in filter.iter().enumerate() {
        if c == b';' {
            ranges.push((start as u32, i as u32));
            start = i + 1;
        }
    }
    ranges.push((start as u32, filter.len() as u32));
    let tag = |r: &(u32, u32)| &filter[r.0 as usize..r.1 as usize];
    ranges.sort_unstable_by(|a, b| tag(a).cmp(tag(b)));
    key.clear();
    for (k, r) in ranges.iter().enumerate() {
        if k > 0 {
            key.push(b';');
        }
        key.extend_from_slice(tag(r));
    }
    &*key
}

/// `--FILTER-summary` (`.FILTER.summary`): per FILTER string, the records
/// and the transitions/transversions among those whose first ALT is a
/// single-base substitution. Rows run from most to fewest records, ties by
/// name in descending order, as vcftools does.
pub struct FilterSummary {
    by_filter: HashMap<Vec<u8>, [u64; 3]>,
    /// Scratch for the sorted key of a multi-tag FILTER.
    key: Vec<u8>,
    ranges: Vec<(u32, u32)>,
    out: Out,
}

impl FilterSummary {
    pub fn new() -> FilterSummary {
        FilterSummary {
            by_filter: HashMap::new(),
            key: Vec::new(),
            ranges: Vec::new(),
            out: Out::closed(),
        }
    }
}

impl Default for FilterSummary {
    fn default() -> Self {
        Self::new()
    }
}

impl Stat for FilterSummary {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, "FILTER.summary")
    }
    fn site(
        &mut self,
        s: &SiteView,
        _gt: Option<&dyn GenotypeSource>,
        _st: Option<&SiteStats>,
    ) -> Result<(), String> {
        let sub = substitution(s.ref_, first_alt(s.alts));
        let key = normalise(&mut self.ranges, &mut self.key, s.filter);
        let e = match self.by_filter.get_mut(key) {
            Some(e) => e,
            None => self.by_filter.entry(key.to_vec()).or_default(),
        };
        e[0] += 1;
        if let Some((_, ts)) = sub {
            e[if ts { 1 } else { 2 }] += 1;
        }
        Ok(())
    }
    fn finish(&mut self) -> io::Result<()> {
        let out = &mut self.out;
        out.buf
            .extend_from_slice(b"FILTER\tN_VARIANTS\tN_Ts\tN_Tv\tTs/Tv\n");
        let mut rows: Vec<(&Vec<u8>, &[u64; 3])> = self.by_filter.iter().collect();
        rows.sort_by(|a, b| b.1[0].cmp(&a.1[0]).then_with(|| b.0.cmp(a.0)));
        for (name, c) in rows {
            let b = &mut out.buf;
            b.extend_from_slice(name);
            for v in c {
                b.push(b'\t');
                write_u64(b, *v);
            }
            b.push(b'\t');
            write_g6(b, c[1] as f64 / c[2] as f64);
            b.push(b'\n');
        }
        out.row().map_err(io::Error::other)?;
        out.commit()
    }
    fn needs_genotypes(&self) -> bool {
        false
    }
    fn ext(&self) -> &'static str {
        "FILTER.summary"
    }
}

/// `--TsTv-summary` (`.TsTv.summary`): substitution counts over biallelic
/// single-base SNPs, then the transition and transversion totals.
pub struct TsTv {
    models: [u64; 6],
    ts: u64,
    tv: u64,
    out: Out,
}

impl TsTv {
    pub fn new() -> TsTv {
        TsTv {
            models: [0; 6],
            ts: 0,
            tv: 0,
            out: Out::closed(),
        }
    }
}

impl Default for TsTv {
    fn default() -> Self {
        Self::new()
    }
}

impl Stat for TsTv {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, "TsTv.summary")
    }
    fn site(
        &mut self,
        s: &SiteView,
        _gt: Option<&dyn GenotypeSource>,
        _st: Option<&SiteStats>,
    ) -> Result<(), String> {
        if s.alts.contains(&b',') {
            return Ok(());
        }
        if let Some((m, ts)) = substitution(s.ref_, s.alts) {
            self.models[m] += 1;
            if ts {
                self.ts += 1;
            } else {
                self.tv += 1;
            }
        }
        Ok(())
    }
    fn finish(&mut self) -> io::Result<()> {
        let b = &mut self.out.buf;
        b.extend_from_slice(b"MODEL\tCOUNT\n");
        let rows = MODELS
            .into_iter()
            .zip(self.models)
            .chain([("Ts", self.ts), ("Tv", self.tv)]);
        for (name, n) in rows {
            b.extend_from_slice(name.as_bytes());
            b.push(b'\t');
            write_u64(b, n);
            b.push(b'\n');
        }
        self.out.row().map_err(io::Error::other)?;
        self.out.commit()
    }
    fn needs_genotypes(&self) -> bool {
        false
    }
    fn ext(&self) -> &'static str {
        "TsTv.summary"
    }
}

/// `--indel-hist` (`.indel.hist`): a histogram of ALT-minus-REF lengths.
/// A site whose REF and ALTs are all one base counts once at length 0; any
/// other site counts each ALT made of A/C/G/T/N whose length differs from
/// REF's (symbolic, breakend, `*` and same-length ALTs do not count).
pub struct IndelHist {
    hist: BTreeMap<i64, u64>,
    out: Out,
}

impl IndelHist {
    pub fn new() -> IndelHist {
        IndelHist {
            hist: BTreeMap::new(),
            out: Out::closed(),
        }
    }
}

impl Default for IndelHist {
    fn default() -> Self {
        Self::new()
    }
}

impl Stat for IndelHist {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, "indel.hist")
    }
    fn site(
        &mut self,
        s: &SiteView,
        _gt: Option<&dyn GenotypeSource>,
        _st: Option<&SiteStats>,
    ) -> Result<(), String> {
        if s.alts == b"." || s.alts.is_empty() {
            return Ok(());
        }
        let ref_len = s.ref_.len() as i64;
        let mut alts = s.alts.split(|&c| c == b',');
        if ref_len == 1 && alts.clone().all(|a| a.len() == 1) {
            *self.hist.entry(0).or_default() += 1;
            return Ok(());
        }
        for a in alts.by_ref() {
            let bases = a
                .iter()
                .all(|c| matches!(c.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'));
            let d = a.len() as i64 - ref_len;
            if bases && d != 0 {
                *self.hist.entry(d).or_default() += 1;
            }
        }
        Ok(())
    }
    fn finish(&mut self) -> io::Result<()> {
        let total: u64 = self.hist.values().sum();
        let b = &mut self.out.buf;
        b.extend_from_slice(b"LENGTH\tCOUNT\tPRCT\n");
        for (&len, &n) in &self.hist {
            write_i64(b, len);
            b.push(b'\t');
            write_u64(b, n);
            b.push(b'\t');
            write_g6(b, n as f64 / total as f64 * 100.0);
            b.push(b'\n');
        }
        self.out.row().map_err(io::Error::other)?;
        self.out.commit()
    }
    fn needs_genotypes(&self) -> bool {
        false
    }
    fn ext(&self) -> &'static str {
        "indel.hist"
    }
}

/// `--snp-density BIN` (`.snpden`): positions per bin of `BIN` bases
/// (records at one position count once, as in vcftools), bins starting at
/// multiples of `BIN` on the 1-based POS (POS 1000 is in bin 1000). Each
/// chromosome, in file order, runs from the bin of its first record to the
/// bin of its last, with empty bins written as zeros.
///
/// Rows stream: only the current bin and its count are held, plus the names
/// of the chromosomes already finished, so memory is bounded by the number
/// of chromosomes, not by records or bins (BIN = 1 included). Input must be
/// sorted by chromosome and position: a chromosome that reappears, or a
/// position in an earlier bin, is an error, because rows already written
/// cannot be merged.
pub struct SnpDensity {
    bin: u64,
    /// The previous record's place: vcftools counts one SNP per position.
    last: (Vec<u8>, u64),
    /// The order check, on the bin: positions may go back within one.
    sorted: SortedInput,
    /// The chromosome being counted, its current bin and that bin's count.
    cur: Option<(Vec<u8>, u64, u64)>,
    out: Out,
}

impl SnpDensity {
    /// `bin` must be positive.
    pub fn new(bin: u64) -> SnpDensity {
        assert!(bin > 0, "the SNP density bin size must be positive");
        SnpDensity {
            bin,
            last: (Vec::new(), u64::MAX),
            sorted: SortedInput::new("--snp-density", "density"),
            cur: None,
            out: Out::closed(),
        }
    }

    fn emit(&mut self, start: u64, n: u64) -> Result<(), String> {
        let b = &mut self.out.buf;
        let chrom = &self.cur.as_ref().expect("a current chromosome").0;
        b.extend_from_slice(chrom);
        b.push(b'\t');
        write_u64(b, start);
        b.push(b'\t');
        write_u64(b, n);
        b.push(b'\t');
        write_g6(b, n as f64 / self.bin as f64 * 1000.0);
        b.push(b'\n');
        self.out.row()
    }
}

impl Stat for SnpDensity {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        self.out.open(sink, "snpden")?;
        self.out
            .buf
            .extend_from_slice(b"CHROM\tBIN_START\tSNP_COUNT\tVARIANTS/KB\n");
        self.out.row().map_err(io::Error::other)
    }
    fn site(
        &mut self,
        s: &SiteView,
        _gt: Option<&dyn GenotypeSource>,
        _st: Option<&SiteStats>,
    ) -> Result<(), String> {
        if self.last.1 == s.pos && self.last.0 == s.chrom {
            return Ok(());
        }
        self.last.0.clear();
        self.last.0.extend_from_slice(s.chrom);
        self.last.1 = s.pos;
        let bin = s.pos / self.bin * self.bin;
        if self.sorted.check_key(s, bin)? {
            if let Some((_, cur_bin, count)) = self.cur.as_ref() {
                let (cur_bin, count) = (*cur_bin, *count);
                self.emit(cur_bin, count)?;
            }
            self.cur = Some((s.chrom.to_vec(), bin, 1));
            return Ok(());
        }
        let (_, cur_bin, count) = self.cur.as_mut().expect("same chromosome");
        if bin == *cur_bin {
            *count += 1;
            return Ok(());
        }
        let (cur_bin, count) = (*cur_bin, *count);
        self.emit(cur_bin, count)?;
        let mut gap = cur_bin + self.bin;
        while gap < bin {
            self.emit(gap, 0)?;
            gap += self.bin;
        }
        let c = self.cur.as_mut().expect("same chromosome");
        c.1 = bin;
        c.2 = 1;
        Ok(())
    }
    fn finish(&mut self) -> io::Result<()> {
        if let Some((_, cur_bin, count)) = self.cur.as_ref() {
            let (cur_bin, count) = (*cur_bin, *count);
            self.emit(cur_bin, count).map_err(io::Error::other)?;
        }
        self.out.commit()
    }
    fn needs_genotypes(&self) -> bool {
        false
    }
    fn ext(&self) -> &'static str {
        "snpden"
    }
}
