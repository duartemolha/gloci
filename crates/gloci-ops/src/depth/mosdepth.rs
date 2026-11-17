//! `gloci depth --out-prefix`: the files mosdepth 0.3.10 writes (D173), computed as the depth
//! core streams the positions. [`Tables`] is the sink: per reference it keeps the depth
//! histogram and the region accumulators, writes the reference's rows when the reference is
//! done, and drops them, so memory is flat in the genome's size.

use super::rows::{DepthSink, ReadFilter, ReadLength};
use super::{AllPositions, CountRules, DepthOptions, ModeRules, run_sources};
use gloci_core::ChromTable;
use gloci_core::cfmt::{write_fixed, write_u64};
use gloci_io::{
    AlnHeader, AlnSource, Format, InvalidPolicy, IoError, Loc, OutFormat, OutputOptions,
    OutputSpec, ReaderOptions, RecordBatch, StagedSet, field, next_batch_size, open_output,
    open_reader, plain_text_output_spec,
};
use std::collections::HashMap;
use std::ffi::OsString;
use std::io::{self, Write};
use std::path::PathBuf;
use std::sync::Arc;

/// Where the regions of `--by` come from.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum By {
    Nothing,
    /// Windows of this many bases on every reference.
    Windows(u64),
    /// The regions of this BED file.
    Bed(String),
}

#[derive(Debug, Clone)]
pub struct MosdepthOptions {
    /// Output path stem: files are `<prefix>.mosdepth.global.dist.txt` and so on.
    pub prefix: PathBuf,
    pub by: By,
    /// Ascending, duplicates kept, each at least 1 (the caller checks).
    pub thresholds: Vec<u32>,
    pub exclude_flags: u16,
    pub min_mapq: u8,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MosdepthSummary {
    pub max_ring: usize,
    pub max_pending_pairs: usize,
    /// Windows or BED regions written.
    pub regions: u64,
    pub files: Vec<PathBuf>,
}

/// Depths above this count at [`CAPPED_SLOT`], as in mosdepth.
const HIST_CAP: u32 = 400_000;
const CAPPED_SLOT: usize = 399_990;
/// A window's mean, rounded, counts in at most this bin of the region distribution.
const WINDOW_MEAN_CAP: usize = 511;
/// Bytes of a table gathered before they are written.
const SPILL_AT: usize = 1 << 16;

pub(crate) fn hist_slot(v: u32) -> usize {
    if v > HIST_CAP {
        CAPPED_SLOT
    } else {
        v as usize
    }
}

/// `hist[slot] += n`, growing the histogram to hold the slot.
fn bump(hist: &mut Vec<u64>, slot: usize, n: u64) {
    if hist.len() <= slot {
        hist.resize(slot + 1, 0);
    }
    hist[slot] += n;
}

fn merge_hist(total: &mut Vec<u64>, h: &[u64]) {
    if total.len() < h.len() {
        total.resize(h.len(), 0);
    }
    for (t, v) in total.iter_mut().zip(h) {
        *t += v;
    }
}

/// mosdepth's distribution lines: for each depth from the highest down, the fraction of bases at
/// that depth or more. Depths above 300 with no bases and fractions below 8e-5 are skipped. The
/// order the fractions are summed in fixes the rounding.
pub(crate) fn write_distribution(name: &[u8], hist: &[u64], out: &mut Vec<u8>) {
    let sum: u64 = hist.iter().sum();
    if sum < 1 {
        return;
    }
    let mut cum = 0f64;
    for d in (0..hist.len()).rev() {
        let v = hist[d];
        if d > 300 && v == 0 {
            continue;
        }
        cum += v as f64 / sum as f64;
        if cum < 8e-5 {
            continue;
        }
        out.extend_from_slice(name);
        out.push(b'\t');
        write_u64(out, d as u64);
        out.push(b'\t');
        write_fixed(out, cum, 2);
        out.push(b'\n');
    }
}

/// The length, depth sum, minimum and maximum over some bases: a summary line.
#[derive(Clone, Copy)]
pub(crate) struct DepthStat {
    cum_length: u64,
    cum_depth: u64,
    /// `u32::MAX` when there are no bases.
    min: u32,
    max: u32,
}

impl Default for DepthStat {
    fn default() -> Self {
        DepthStat {
            cum_length: 0,
            cum_depth: 0,
            min: u32::MAX,
            max: 0,
        }
    }
}

impl DepthStat {
    /// `n` bases at depth `v`.
    pub(crate) fn add(&mut self, v: u32, n: u64) {
        if n == 0 {
            return;
        }
        self.cum_length += n;
        self.cum_depth += u64::from(v) * n;
        self.min = self.min.min(v);
        self.max = self.max.max(v);
    }

    pub(crate) fn merge(&mut self, o: &DepthStat) {
        self.cum_length += o.cum_length;
        self.cum_depth += o.cum_depth;
        self.min = self.min.min(o.min);
        self.max = self.max.max(o.max);
    }

    /// `name length bases mean min max`; a minimum over no bases prints 0, and so does the mean.
    pub(crate) fn write_line(&self, name: &[u8], out: &mut Vec<u8>) {
        let mean = if self.cum_length > 0 {
            self.cum_depth as f64 / self.cum_length as f64
        } else {
            0.0
        };
        out.extend_from_slice(name);
        for v in [self.cum_length, self.cum_depth] {
            out.push(b'\t');
            write_u64(out, v);
        }
        out.push(b'\t');
        write_fixed(out, mean, 2);
        let min = if self.min == u32::MAX { 0 } else { self.min };
        for v in [min, self.max] {
            out.push(b'\t');
            write_u64(out, u64::from(v));
        }
        out.push(b'\n');
    }
}

/// One region of a BED file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct BedRegion {
    pub start: u64,
    pub end: u64,
    pub name: Option<Box<[u8]>>,
}

/// The regions of `bed` on each reference of `hdr` (by reference id), stable-sorted by start. A
/// chromosome the header lacks is skipped with one warning naming `input`; a region that ends
/// past its reference is an error naming the line.
pub(crate) fn load_bed_regions(
    bed: &str,
    hdr: &AlnHeader,
    input: &str,
    ro: &ReaderOptions,
) -> Result<Vec<Vec<BedRegion>>, IoError> {
    let o = ReaderOptions {
        format: Some(Format::Bed),
        invalid: InvalidPolicy::Fail,
        chrom_map: None,
        ..ro.clone()
    };
    let mut r = open_reader(bed, &o)?;
    let by_name: HashMap<&[u8], usize> = hdr
        .refs
        .iter()
        .enumerate()
        .map(|(i, (n, _))| (&n[..], i))
        .collect();
    let mut per: Vec<Vec<BedRegion>> = vec![Vec::new(); hdr.refs.len()];
    let mut skipped: Vec<(Vec<u8>, u64)> = Vec::new();
    let mut skipped_at: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut chroms = ChromTable::first_appearance();
    let mut batch = RecordBatch::default();
    let mut max = 1;
    loop {
        batch.clear();
        let more = r.read_batch(&mut chroms, &mut batch, max)?;
        for (i, m) in batch.recs.iter().enumerate() {
            let chrom = chroms.name(m.chrom);
            let Some(&rid) = by_name.get(chrom) else {
                let k = *skipped_at.entry(chrom.to_vec()).or_insert_with(|| {
                    skipped.push((chrom.to_vec(), 0));
                    skipped.len() - 1
                });
                skipped[k].1 += 1;
                continue;
            };
            let ln = hdr.refs[rid].1;
            if m.span.end > ln {
                let c = String::from_utf8_lossy(chrom);
                return Err(IoError::Parse {
                    path: bed.to_string(),
                    loc: Loc::Line(m.loc),
                    msg: format!(
                        "{c} {} {} ends past the end of {c} (LN {ln})",
                        m.span.start, m.span.end
                    ),
                });
            }
            let name = field(batch.raw(i), 3)
                .filter(|n| !n.is_empty())
                .map(|n| n.to_vec().into_boxed_slice());
            per[rid].push(BedRegion {
                start: m.span.start,
                end: m.span.end,
                name,
            });
        }
        max = next_batch_size(&batch, max);
        if !more {
            break;
        }
    }
    for (chrom, n) in skipped {
        log::warn!(
            "{bed}: {n} region(s) on {} skipped: {input} has no such reference",
            String::from_utf8_lossy(&chrom)
        );
    }
    for regions in &mut per {
        regions.sort_by_key(|r| r.start);
    }
    Ok(per)
}

/// A table being written: lines gather in `buf` and are handed to `w` in chunks.
struct Buf<'w> {
    w: &'w mut dyn Write,
    buf: Vec<u8>,
}

impl<'w> Buf<'w> {
    fn new(w: &'w mut dyn Write) -> Buf<'w> {
        Buf { w, buf: Vec::new() }
    }

    fn spill(&mut self) -> io::Result<()> {
        if self.buf.len() >= SPILL_AT {
            self.w.write_all(&self.buf)?;
            self.buf.clear();
        }
        Ok(())
    }

    fn finish(&mut self) -> io::Result<()> {
        self.w.write_all(&self.buf)?;
        self.buf.clear();
        self.w.flush()
    }
}

/// The targets of [`Tables`]: the two files every run writes, and those `--by` and
/// `--thresholds` add.
pub(crate) struct Targets<'w> {
    pub global_dist: &'w mut dyn Write,
    pub summary: &'w mut dyn Write,
    pub region_dist: Option<&'w mut dyn Write>,
    pub regions: Option<&'w mut dyn Write>,
    pub thresholds: Option<&'w mut dyn Write>,
}

/// The regions a run is over.
pub(crate) enum Regions {
    None,
    Windows(u64),
    /// The BED regions of each reference, sorted by start.
    Bed(Vec<Vec<BedRegion>>),
}

/// The window being filled.
struct Win {
    start: u64,
    end: u64,
    mean: f64,
    counts: Vec<u64>,
}

/// The accumulators of the BED regions of one reference.
struct BedState {
    /// The first region not yet in `active`.
    next: usize,
    /// Regions that have started and not ended, by index.
    active: Vec<usize>,
    mean: Vec<f64>,
    /// `thresholds.len()` counts per region.
    counts: Vec<u64>,
}

/// What is kept for the reference whose positions are arriving.
struct Cur {
    rid: usize,
    hist: Vec<u64>,
    stat: DepthStat,
    region_hist: Vec<u64>,
    region_stat: DepthStat,
    win: Option<Win>,
    bed: Option<BedState>,
}

fn region_row(
    out: &mut Vec<u8>,
    chrom: &[u8],
    start: u64,
    end: u64,
    name: Option<&[u8]>,
    mean: f64,
) {
    out.extend_from_slice(chrom);
    for v in [start, end] {
        out.push(b'\t');
        write_u64(out, v);
    }
    if let Some(n) = name {
        out.push(b'\t');
        out.extend_from_slice(n);
    }
    out.push(b'\t');
    write_fixed(out, mean, 2);
    out.push(b'\n');
}

fn threshold_row(
    out: &mut Vec<u8>,
    chrom: &[u8],
    start: u64,
    end: u64,
    name: Option<&[u8]>,
    counts: &[u64],
) {
    out.extend_from_slice(chrom);
    for v in [start, end] {
        out.push(b'\t');
        write_u64(out, v);
    }
    out.push(b'\t');
    out.extend_from_slice(name.unwrap_or(b"unknown"));
    for &c in counts {
        out.push(b'\t');
        write_u64(out, c);
    }
    out.push(b'\n');
}

/// mosdepth's tables as a [`DepthSink`]. The positions of every reference must arrive
/// (`AllPositions::AllRefs`); `finish` writes what comes last.
pub(crate) struct Tables<'w> {
    hdr: Arc<AlnHeader>,
    regions: Regions,
    thresholds: Vec<u32>,
    global: Buf<'w>,
    summary: Buf<'w>,
    region_dist: Option<Buf<'w>>,
    region_rows: Option<Buf<'w>>,
    threshold_rows: Option<Buf<'w>>,
    placed: Vec<bool>,
    cur: Option<Cur>,
    total_hist: Vec<u64>,
    total_stat: DepthStat,
    total_region_hist: Vec<u64>,
    total_region_stat: DepthStat,
    summary_header: bool,
    /// Windows or BED regions written.
    pub(crate) region_count: u64,
}

impl<'w> Tables<'w> {
    pub(crate) fn new(
        hdr: &Arc<AlnHeader>,
        regions: Regions,
        thresholds: Vec<u32>,
        t: Targets<'w>,
    ) -> io::Result<Tables<'w>> {
        let mut threshold_rows = t.thresholds.map(Buf::new);
        if let Some(b) = &mut threshold_rows {
            b.buf.extend_from_slice(b"#chrom\tstart\tend\tregion");
            for t in &thresholds {
                b.buf.push(b'\t');
                write_u64(&mut b.buf, u64::from(*t));
                b.buf.push(b'X');
            }
            b.buf.push(b'\n');
        }
        Ok(Tables {
            hdr: hdr.clone(),
            regions,
            thresholds,
            global: Buf::new(t.global_dist),
            summary: Buf::new(t.summary),
            region_dist: t.region_dist.map(Buf::new),
            region_rows: t.regions.map(Buf::new),
            threshold_rows,
            placed: vec![false; hdr.refs.len()],
            cur: None,
            total_hist: Vec::new(),
            total_stat: DepthStat::default(),
            total_region_hist: Vec::new(),
            total_region_stat: DepthStat::default(),
            summary_header: false,
            region_count: 0,
        })
    }

    fn has_regions(&self) -> bool {
        !matches!(self.regions, Regions::None)
    }

    fn start_ref(&self, rid: usize) -> Cur {
        let ln = self.hdr.refs[rid].1;
        let n = self.thresholds.len();
        let win = match self.regions {
            Regions::Windows(size) => Some(Win {
                start: 0,
                end: size.min(ln),
                mean: 0.0,
                counts: vec![0; n],
            }),
            _ => None,
        };
        let bed = match &self.regions {
            Regions::Bed(per) => Some(BedState {
                next: 0,
                active: Vec::new(),
                mean: vec![0.0; per[rid].len()],
                counts: vec![0; per[rid].len() * n],
            }),
            _ => None,
        };
        Cur {
            rid,
            hist: Vec::new(),
            stat: DepthStat::default(),
            region_hist: Vec::new(),
            region_stat: DepthStat::default(),
            win,
            bed,
        }
    }

    /// The reference `rid`'s state, started if it is not the current one.
    fn cur_for(&mut self, rid: usize) -> Cur {
        match self.cur.take() {
            Some(c) if c.rid == rid => c,
            _ => self.start_ref(rid),
        }
    }

    /// Positions `[a, b)` of `rid` at depth `v`.
    fn add_run(&mut self, rid: usize, a: u64, b: u64, v: u32) -> io::Result<()> {
        if a >= b {
            return Ok(());
        }
        let mut cur = self.cur_for(rid);
        let n = b - a;
        bump(&mut cur.hist, hist_slot(v), n);
        cur.stat.add(v, n);
        let r = self.add_to_regions(&mut cur, a, b, v);
        self.cur = Some(cur);
        r
    }

    fn add_to_regions(&mut self, cur: &mut Cur, a: u64, b: u64, v: u32) -> io::Result<()> {
        let ln = self.hdr.refs[cur.rid].1;
        if let (Regions::Windows(size), Some(win)) = (&self.regions, &mut cur.win) {
            cur.region_stat.add(v, b - a);
            let mut p = a;
            while p < b {
                while p >= win.end {
                    finish_window(
                        win,
                        &self.hdr.refs[cur.rid].0,
                        &self.thresholds,
                        &mut cur.region_hist,
                        (&mut self.region_rows, &mut self.threshold_rows),
                        &mut self.region_count,
                    )?;
                    win.start = win.end;
                    win.end = (win.start + size).min(ln);
                }
                let end = b.min(win.end);
                if v > 0 {
                    let q = f64::from(v) / (win.end - win.start) as f64;
                    for _ in p..end {
                        win.mean += q;
                    }
                    for (c, t) in win.counts.iter_mut().zip(&self.thresholds) {
                        if v < *t {
                            break;
                        }
                        *c += end - p;
                    }
                }
                p = end;
            }
        }
        if let (Regions::Bed(per), Some(st)) = (&self.regions, &mut cur.bed) {
            let regs = &per[cur.rid];
            while st.next < regs.len() && regs[st.next].start < b {
                if regs[st.next].end > regs[st.next].start {
                    st.active.push(st.next);
                }
                st.next += 1;
            }
            let nth = self.thresholds.len();
            let mut keep = 0;
            for k in 0..st.active.len() {
                let i = st.active[k];
                let r = &regs[i];
                let (lo, hi) = (a.max(r.start), b.min(r.end));
                if hi > lo {
                    let n = hi - lo;
                    cur.region_stat.add(v, n);
                    bump(&mut cur.region_hist, hist_slot(v), n);
                    if v > 0 {
                        let q = f64::from(v) / (r.end - r.start) as f64;
                        for _ in 0..n {
                            st.mean[i] += q;
                        }
                        for (j, t) in self.thresholds.iter().enumerate() {
                            if v < *t {
                                break;
                            }
                            st.counts[i * nth + j] += n;
                        }
                    }
                }
                if r.end > b {
                    st.active[keep] = i;
                    keep += 1;
                }
            }
            st.active.truncate(keep);
        }
        if let Some(b) = &mut self.region_rows {
            b.spill()?;
        }
        if let Some(b) = &mut self.threshold_rows {
            b.spill()?;
        }
        Ok(())
    }

    fn write_summary_header(&mut self) {
        if !self.summary_header {
            self.summary_header = true;
            self.summary
                .buf
                .extend_from_slice(b"chrom\tlength\tbases\tmean\tmin\tmax\n");
        }
    }

    /// What comes last: the `total` distribution lines, then the summary's `total` and
    /// `total_region` lines. Flushes every table.
    pub(crate) fn finish(&mut self) -> io::Result<()> {
        write_distribution(b"total", &self.total_hist, &mut self.global.buf);
        if let Some(b) = &mut self.region_dist {
            write_distribution(b"total", &self.total_region_hist, &mut b.buf);
        }
        self.write_summary_header();
        self.total_stat.write_line(b"total", &mut self.summary.buf);
        if self.has_regions() {
            self.total_region_stat
                .write_line(b"total_region", &mut self.summary.buf);
        }
        self.global.finish()?;
        self.summary.finish()?;
        for b in [
            &mut self.region_dist,
            &mut self.region_rows,
            &mut self.threshold_rows,
        ]
        .into_iter()
        .flatten()
        {
            b.finish()?;
        }
        Ok(())
    }
}

/// Write the finished window's rows, count it in the region distribution and reset it.
fn finish_window(
    win: &mut Win,
    chrom: &[u8],
    thresholds: &[u32],
    region_hist: &mut Vec<u64>,
    rows: (&mut Option<Buf<'_>>, &mut Option<Buf<'_>>),
    count: &mut u64,
) -> io::Result<()> {
    let mean = win.mean;
    if let Some(b) = rows.0 {
        region_row(&mut b.buf, chrom, win.start, win.end, None, mean);
    }
    if let Some(b) = rows.1 {
        threshold_row(&mut b.buf, chrom, win.start, win.end, None, &win.counts);
    }
    bump(region_hist, (mean.round() as usize).min(WINDOW_MEAN_CAP), 1);
    *count += 1;
    win.mean = 0.0;
    win.counts.iter_mut().for_each(|c| *c = 0);
    debug_assert_eq!(win.counts.len(), thresholds.len());
    Ok(())
}

impl DepthSink for Tables<'_> {
    fn position(&mut self, rid: usize, pos: u64, depths: &[u32]) -> io::Result<()> {
        let v = depths.iter().fold(0u32, |s, &d| s.saturating_add(d));
        self.add_run(rid, pos, pos + 1, v)
    }

    fn zeros(&mut self, rid: usize, a: u64, b: u64) -> io::Result<()> {
        self.add_run(rid, a, b, 0)
    }

    fn placed(&mut self, rid: usize) {
        self.placed[rid] = true;
    }

    fn reference_done(&mut self, rid: usize) -> io::Result<()> {
        let mut cur = self.cur_for(rid);
        let ln = self.hdr.refs[rid].1;
        let name = self.hdr.refs[rid].0.clone();
        // The windows the positions did not reach (none, when every position arrived).
        if let (Regions::Windows(size), Some(win)) = (&self.regions, &mut cur.win)
            && ln > 0
        {
            loop {
                finish_window(
                    win,
                    &name,
                    &self.thresholds,
                    &mut cur.region_hist,
                    (&mut self.region_rows, &mut self.threshold_rows),
                    &mut self.region_count,
                )?;
                if win.end >= ln {
                    break;
                }
                win.start = win.end;
                win.end = (win.start + size).min(ln);
            }
        }
        if let (Regions::Bed(per), Some(st)) = (&mut self.regions, &cur.bed) {
            let regs = std::mem::take(&mut per[rid]);
            let nth = self.thresholds.len();
            for (i, r) in regs.iter().enumerate() {
                let n = r.name.as_deref();
                if let Some(b) = &mut self.region_rows {
                    region_row(&mut b.buf, &name, r.start, r.end, n, st.mean[i]);
                }
                if let Some(b) = &mut self.threshold_rows {
                    let counts = &st.counts[i * nth..(i + 1) * nth];
                    threshold_row(&mut b.buf, &name, r.start, r.end, n, counts);
                }
                self.region_count += 1;
            }
        }
        if self.placed[rid] {
            write_distribution(&name, &cur.hist, &mut self.global.buf);
            if let Some(b) = &mut self.region_dist {
                write_distribution(&name, &cur.region_hist, &mut b.buf);
            }
            self.write_summary_header();
            cur.stat.write_line(&name, &mut self.summary.buf);
            if self.has_regions() {
                let mut region_name = name.clone();
                region_name.extend_from_slice(b"_region");
                cur.region_stat
                    .write_line(&region_name, &mut self.summary.buf);
            }
            merge_hist(&mut self.total_hist, &cur.hist);
            self.total_stat.merge(&cur.stat);
            merge_hist(&mut self.total_region_hist, &cur.region_hist);
            self.total_region_stat.merge(&cur.region_stat);
        }
        self.global.spill()?;
        self.summary.spill()?;
        for b in [
            &mut self.region_dist,
            &mut self.region_rows,
            &mut self.threshold_rows,
        ]
        .into_iter()
        .flatten()
        {
            b.spill()?;
        }
        Ok(())
    }
}

fn with_suffix(prefix: &std::path::Path, suffix: &str) -> PathBuf {
    let mut s: OsString = prefix.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}

/// mosdepth's files for `src` (opened from `path` with no reader filters), as D173 part B says.
/// Every file is staged and published together; on an error none is.
pub fn mosdepth_files(
    src: Box<dyn AlnSource>,
    path: &str,
    o: &MosdepthOptions,
    ro: &ReaderOptions,
    oo: &OutputOptions,
) -> Result<MosdepthSummary, IoError> {
    let hdr = src.header().clone();
    let regions = match &o.by {
        By::Nothing => Regions::None,
        By::Windows(size) => Regions::Windows(*size),
        By::Bed(bed) => Regions::Bed(load_bed_regions(bed, &hdr, path, ro)?),
    };
    let has_regions = !matches!(regions, Regions::None);
    // The tables are text and the index is the CSI of the BED files, never `--write-index`.
    let oo = OutputOptions {
        write_index: false,
        ..*oo
    };
    let text = |suffix: &str| -> Result<_, IoError> {
        let spec = plain_text_output_spec(Some(&with_suffix(&o.prefix, suffix)));
        open_output(&spec, OutFormat::Tsv, &oo)
    };
    let bgzf = |suffix: &str| -> Result<_, IoError> {
        let spec = OutputSpec {
            path: Some(with_suffix(&o.prefix, suffix)),
            format: Some(OutFormat::Bed),
            bgzf: true,
            direct: false,
        };
        open_output(&spec, OutFormat::Bed, &oo)
    };
    let mut global = text(".mosdepth.global.dist.txt")?;
    let mut summary = text(".mosdepth.summary.txt")?;
    let mut region_dist = has_regions
        .then(|| text(".mosdepth.region.dist.txt"))
        .transpose()?;
    let mut region_rows = has_regions.then(|| bgzf(".regions.bed.gz")).transpose()?;
    let mut threshold_rows = (has_regions && !o.thresholds.is_empty())
        .then(|| bgzf(".thresholds.bed.gz"))
        .transpose()?;

    let d = DepthOptions {
        rules: CountRules {
            deletions: false,
            ref_skips: false,
            min_baseq: 0,
        },
        min_read_length: 0,
        all: AllPositions::AllRefs,
        regions: None,
        skip_mate_overlap: false,
        column_names: false,
    };
    let mode = ModeRules {
        read_length: ReadLength::SeqLessClips,
        filter: Some(ReadFilter {
            exclude_flags: o.exclude_flags,
            min_mapq: o.min_mapq,
        }),
        shared_bases: true,
    };
    let (run, mut tables) = run_sources(vec![src], &[path.to_string()], &d, mode, ro, |hdr| {
        let t = Targets {
            global_dist: &mut global,
            summary: &mut summary,
            region_dist: region_dist.as_mut().map(|w| w as &mut dyn Write),
            regions: region_rows.as_mut().map(|w| w as &mut dyn Write),
            thresholds: threshold_rows.as_mut().map(|w| w as &mut dyn Write),
        };
        Tables::new(hdr, regions, o.thresholds.clone(), t).map_err(|source| IoError::Write {
            path: "output".into(),
            source,
        })
    })?;
    tables.finish().map_err(|source| IoError::Write {
        path: "output".into(),
        source,
    })?;
    let region_count = tables.region_count;
    drop(tables);

    let mut set = StagedSet::default();
    global.finish_staged(&mut set)?;
    summary.finish_staged(&mut set)?;
    if let Some(w) = region_dist {
        w.finish_staged(&mut set)?;
    }
    if let Some(w) = region_rows {
        w.finish_staged_with_csi(Format::Bed, &mut set)?;
    }
    if let Some(w) = threshold_rows {
        w.finish_staged_with_csi(Format::Bed, &mut set)?;
    }
    let files = set.publish()?;
    Ok(MosdepthSummary {
        max_ring: run.max_ring,
        max_pending_pairs: run.max_pending_mates,
        regions: region_count,
        files,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use gloci_core::VcfSpan;
    use gloci_io::{Accept, AlignOptions, open_aln_source};
    use std::num::NonZeroUsize;
    use std::sync::Mutex;

    fn ro() -> ReaderOptions {
        ReaderOptions {
            accept: Accept {
                alignments: Some(AlignOptions {
                    exclude_flags: 0,
                    ..AlignOptions::default()
                }),
                ..Accept::default()
            },
            chrom_map: None,
            threads: NonZeroUsize::MIN,
            format: None,
            vcf_span: VcfSpan::Ref,
            invalid: InvalidPolicy::Fail,
        }
    }

    fn input(file: &str) -> String {
        format!(
            "{}/../../xtask/diff/inputs/samtools/{file}",
            env!("CARGO_MANIFEST_DIR")
        )
    }

    fn header(refs: &[(&str, u64)]) -> Arc<AlnHeader> {
        let mut text = String::from("@HD\tVN:1.6\n");
        for (n, l) in refs {
            text.push_str(&format!("@SQ\tSN:{n}\tLN:{l}\n"));
        }
        let src = crate::depth::test_sink::sam_source(&text, &ro());
        src.header().clone()
    }

    fn dist(name: &str, pairs: &[(usize, u64)]) -> String {
        let len = pairs.iter().map(|p| p.0).max().map_or(0, |m| m + 1);
        let mut hist = vec![0; len];
        for &(d, n) in pairs {
            hist[d] = n;
        }
        let mut out = Vec::new();
        write_distribution(name.as_bytes(), &hist, &mut out);
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn distribution_lines_follow_mosdepth() {
        assert_eq!(
            dist("chr4", &[(0, 19969), (1, 20), (2, 10), (320, 1)]),
            "chr4\t2\t0.00\nchr4\t1\t0.00\nchr4\t0\t1.00\n"
        );
        let all = dist("chr4", &[(320, 1)]);
        let lines: Vec<&str> = all.lines().collect();
        assert_eq!(lines.len(), 302);
        assert_eq!(lines[0], "chr4\t320\t1.00");
        assert_eq!(lines[1], "chr4\t300\t1.00");
        assert_eq!(lines[301], "chr4\t0\t1.00");
        assert_eq!(dist("chr4", &[(5, 0)]), "");
        let mut out = Vec::new();
        write_distribution(b"x", &[], &mut out);
        assert!(out.is_empty());
    }

    #[test]
    fn depth_stat_lines() {
        let line = |s: &DepthStat, name: &str| {
            let mut out = Vec::new();
            s.write_line(name.as_bytes(), &mut out);
            String::from_utf8(out).unwrap()
        };
        let mut s = DepthStat::default();
        s.add(0, 2);
        s.add(3, 1);
        assert_eq!(line(&s, "x"), "x\t3\t3\t1.00\t0\t3\n");
        assert_eq!(
            line(&DepthStat::default(), "x_region"),
            "x_region\t0\t0\t0.00\t0\t0\n"
        );
        let mut t = DepthStat::default();
        t.add(7, 1);
        t.merge(&s);
        assert_eq!(line(&t, "t"), "t\t4\t10\t2.50\t0\t7\n");
    }

    #[test]
    fn hist_caps_at_400000() {
        assert_eq!(hist_slot(400_001), 399_990);
        assert_eq!(hist_slot(400_000), 400_000);
        assert_eq!(hist_slot(7), 7);
    }

    /// The tables of one run over `refs`, fed `runs` (rid, start, end, depth) and `placed`.
    struct Out {
        global: String,
        summary: String,
        region_dist: String,
        regions: String,
        thresholds: String,
    }

    fn tables(
        refs: &[(&str, u64)],
        regions: Regions,
        thresholds: &[u32],
        placed: &[usize],
        runs: &[(usize, u64, u64, u32)],
    ) -> Out {
        let hdr = header(refs);
        let (mut g, mut s, mut rd, mut r, mut t) =
            (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
        {
            let has = !matches!(regions, Regions::None);
            let targets = Targets {
                global_dist: &mut g,
                summary: &mut s,
                region_dist: has.then_some(&mut rd as &mut dyn Write),
                regions: has.then_some(&mut r as &mut dyn Write),
                thresholds: (!thresholds.is_empty()).then_some(&mut t as &mut dyn Write),
            };
            let mut tb = Tables::new(&hdr, regions, thresholds.to_vec(), targets).unwrap();
            for &p in placed {
                tb.placed(p);
            }
            for rid in 0..refs.len() {
                for &(r, a, b, v) in runs.iter().filter(|x| x.0 == rid) {
                    if v == 0 {
                        tb.zeros(r, a, b).unwrap();
                    } else {
                        for p in a..b {
                            tb.position(r, p, &[v]).unwrap();
                        }
                    }
                }
                tb.reference_done(rid).unwrap();
            }
            tb.finish().unwrap();
        }
        let st = |v: Vec<u8>| String::from_utf8(v).unwrap();
        Out {
            global: st(g),
            summary: st(s),
            region_dist: st(rd),
            regions: st(r),
            thresholds: st(t),
        }
    }

    #[test]
    fn region_mean_accumulates_base_by_base() {
        let o = tables(
            &[("c", 400)],
            Regions::Windows(200),
            &[],
            &[0],
            &[
                (0, 0, 11, 1),
                (0, 11, 200, 0),
                (0, 200, 225, 1),
                (0, 225, 400, 0),
            ],
        );
        // The exact quotients would print 0.06 and 0.12; mosdepth 0.3.10 printed these.
        assert_eq!(o.regions, "c\t0\t200\t0.05\nc\t200\t400\t0.13\n");
    }

    #[test]
    fn window_dist_rounds_half_away_from_zero() {
        // Windows of 100 on c: 50 bases at depth 1 (mean .5), 49 at depth 1 (mean .49), and a
        // window at depth 600.
        let o = tables(
            &[("c", 300)],
            Regions::Windows(100),
            &[],
            &[0],
            &[
                (0, 0, 50, 1),
                (0, 50, 100, 0),
                (0, 100, 149, 1),
                (0, 149, 200, 0),
                (0, 200, 300, 600),
            ],
        );
        // Three windows: bin 1 (mean .5), bin 0 (.49), bin 511 (600). The distribution keeps
        // each depth up to 300 (a third of the windows reach each of 2..300 through bin 511).
        let line = |d: usize| format!("c\t{d}\t");
        let find = |d: usize| {
            o.region_dist
                .lines()
                .find(|l| l.starts_with(&line(d)))
                .unwrap()
        };
        assert_eq!(find(511), "c\t511\t0.33");
        assert_eq!(find(2), "c\t2\t0.33");
        assert_eq!(find(1), "c\t1\t0.67");
        assert_eq!(find(0), "c\t0\t1.00");
        assert!(o.region_dist.contains("total\t1\t0.67\ntotal\t0\t1.00\n"));
    }

    #[test]
    fn reference_without_records_is_left_out() {
        let o = tables(
            &[("c1", 10), ("c2", 10)],
            Regions::Windows(10),
            &[],
            &[0],
            &[(0, 0, 10, 0), (1, 0, 10, 0)],
        );
        assert_eq!(o.global, "c1\t0\t1.00\ntotal\t0\t1.00\n");
        assert!(o.summary.contains("c1\t10\t0\t0.00\t0\t0\n"));
        assert!(o.summary.contains("c1_region\t10\t0\t0.00\t0\t0\n"));
        assert!(!o.summary.contains("c2"));
        assert!(
            o.summary
                .starts_with("chrom\tlength\tbases\tmean\tmin\tmax\n")
        );
        assert!(o.summary.ends_with("total_region\t10\t0\t0.00\t0\t0\n"));
        assert_eq!(o.regions, "c1\t0\t10\t0.00\nc2\t0\t10\t0.00\n");
    }

    #[test]
    fn thresholds_header_and_rows() {
        let o = tables(
            &[("c", 10)],
            Regions::Windows(10),
            &[1, 1, 2],
            &[0],
            &[(0, 0, 4, 1), (0, 4, 6, 2), (0, 6, 10, 0)],
        );
        assert_eq!(
            o.thresholds,
            "#chrom\tstart\tend\tregion\t1X\t1X\t2X\nc\t0\t10\tunknown\t6\t6\t2\n"
        );
    }

    #[test]
    fn bed_regions_count_each_base_in_each_region() {
        // Two overlapping regions and a named one.
        let regions = Regions::Bed(vec![vec![
            BedRegion {
                start: 0,
                end: 4,
                name: Some(b"a".to_vec().into_boxed_slice()),
            },
            BedRegion {
                start: 2,
                end: 6,
                name: None,
            },
        ]]);
        let o = tables(
            &[("c", 8)],
            regions,
            &[2],
            &[0],
            &[(0, 0, 3, 1), (0, 3, 5, 2), (0, 5, 8, 0)],
        );
        assert_eq!(o.regions, "c\t0\t4\ta\t1.25\nc\t2\t6\t1.25\n");
        assert_eq!(
            o.thresholds,
            "#chrom\tstart\tend\tregion\t2X\nc\t0\t4\ta\t1\nc\t2\t6\tunknown\t2\n"
        );
        // Region bases: 8, depth sum 10, a base in two regions counted twice.
        assert!(
            o.summary.contains("c_region\t8\t10\t1.25\t0\t2\n"),
            "{}",
            o.summary
        );
    }

    fn bed_regions(file: &str) -> Result<Vec<Vec<BedRegion>>, IoError> {
        let bam = input("mosdepth.bam");
        let src = open_aln_source(&bam, &ro()).unwrap();
        load_bed_regions(&input(file), src.header(), &bam, &ro())
    }

    fn brief(r: &[BedRegion]) -> Vec<(u64, u64, String)> {
        r.iter()
            .map(|r| {
                (
                    r.start,
                    r.end,
                    r.name
                        .as_deref()
                        .map_or("-".to_string(), |n| String::from_utf8_lossy(n).into_owned()),
                )
            })
            .collect()
    }

    struct Capture;
    static WARNINGS: Mutex<Vec<String>> = Mutex::new(Vec::new());
    impl log::Log for Capture {
        fn enabled(&self, _: &log::Metadata) -> bool {
            true
        }
        fn log(&self, r: &log::Record) {
            WARNINGS.lock().unwrap().push(r.args().to_string());
        }
        fn flush(&self) {}
    }

    #[test]
    fn bed_regions_sorted_stably_and_named() {
        let _ = log::set_logger(&Capture);
        log::set_max_level(log::LevelFilter::Warn);
        let per = bed_regions("mosdepth.bed").unwrap();
        let s = |a: u64, b: u64, n: &str| (a, b, n.to_string());
        assert_eq!(
            brief(&per[0]),
            [
                s(100, 180, "p"),
                s(200, 260, "q"),
                s(200, 210, "q-head"),
                s(500, 550, "-"),
                s(650, 650, "point"),
                s(700, 800, "v")
            ]
        );
        assert_eq!(brief(&per[1]), [s(0, 300, "dups2")]);
        assert_eq!(brief(&per[2]), [s(0, 200, "empty3")]);
        assert_eq!(brief(&per[3]), [s(5000, 5001, "deep")]);
        let warned: Vec<String> = WARNINGS
            .lock()
            .unwrap()
            .iter()
            .filter(|w| w.contains("chrZ") && w.contains("mosdepth.bed"))
            .cloned()
            .collect();
        assert_eq!(warned.len(), 1, "{warned:?}");
        assert!(
            warned[0].contains("1 region(s) on chrZ skipped"),
            "{warned:?}"
        );
    }

    #[test]
    fn bed_region_past_the_end_is_an_error() {
        let dir = std::env::temp_dir().join(format!("gloci-mosdepth-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let f = dir.join("past.bed");
        std::fs::write(&f, "chr1\t480\t1001\n").unwrap();
        let bam = input("mosdepth.bam");
        let src = open_aln_source(&bam, &ro()).unwrap();
        let e = load_bed_regions(f.to_str().unwrap(), src.header(), &bam, &ro())
            .unwrap_err()
            .to_string();
        std::fs::remove_dir_all(&dir).unwrap();
        assert!(e.contains("line 1") && e.contains("LN 1000"), "{e}");
        assert!(e.contains("past.bed"), "{e}");
    }

    #[test]
    fn blank_bed_lines_are_skipped() {
        let per = bed_regions("mosdepth-blank.bed").unwrap();
        let names: Vec<_> = brief(&per[0]).into_iter().map(|r| r.2).collect();
        assert_eq!(names, ["a", "b"]);
    }

    #[test]
    fn files_are_staged_and_published_together() {
        let dir = std::env::temp_dir().join(format!("gloci-mosdepth-files-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bam = input("mosdepth.bam");
        // A copy, so that the warning about chrZ is this test's own.
        let bed = dir.join("targets.bed");
        std::fs::copy(input("mosdepth.bed"), &bed).unwrap();
        let o = MosdepthOptions {
            prefix: dir.join("p"),
            by: By::Bed(bed.to_string_lossy().into_owned()),
            thresholds: vec![1, 5],
            exclude_flags: 1796,
            min_mapq: 0,
        };
        let oo = OutputOptions {
            threads: NonZeroUsize::MIN,
            level: 6,
            write_index: true,
        };
        let src = open_aln_source(&bam, &ro()).unwrap();
        let s = mosdepth_files(src, &bam, &o, &ro(), &oo).unwrap();
        let mut names: Vec<String> = s
            .files
            .iter()
            .map(|f| f.file_name().unwrap().to_string_lossy().into_owned())
            .collect();
        names.sort();
        assert_eq!(
            names,
            [
                "p.mosdepth.global.dist.txt",
                "p.mosdepth.region.dist.txt",
                "p.mosdepth.summary.txt",
                "p.regions.bed.gz",
                "p.regions.bed.gz.csi",
                "p.thresholds.bed.gz",
                "p.thresholds.bed.gz.csi"
            ]
        );
        assert_eq!(s.regions, 9);
        let summary = std::fs::read_to_string(dir.join("p.mosdepth.summary.txt")).unwrap();
        assert!(summary.starts_with("chrom\tlength\tbases\tmean\tmin\tmax\nchr1\t1000\t"));
        assert!(summary.contains("\ntotal_region\t"), "{summary}");
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
