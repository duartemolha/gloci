//! The rows of `depth`: which positions print, in which order, for one run of references
//! (D138). One run reads one record stream per input, merged by `(reference, position)`.

use super::mates::MateTable;
use super::pairs::PairTable;
use super::ring::{CountRules, DepthRing};
use super::{AllPositions, DepthOptions, ModeRules};
use gloci_core::Interval;
use gloci_io::{AlnHeader, AlnRef, AlnSource, Diagnostics, IoError, Loc, OrderError, aln_view};
use std::io::{self, Write};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

const FLAG_UNMAPPED: u16 = 0x4;
/// Bytes of rows gathered before they are handed on: as large as the output's own buffer, so
/// a chunk passes through it without being copied again.
const FLUSH_AT: usize = 1 << 20;

/// Where a depth run's positions go: `depth`'s rows (`Emitter`), the `--summary` table, or
/// mosdepth's files. Positions arrive in `@SQ` order, then position order, each at most once.
///
/// `read` and `placed` run one record ahead of the positions, so a record on reference N+1 can
/// arrive before `reference_done(N)` and before N's last positions (and, with several inputs,
/// before another input's positions): key per-reference state by the record's reference id,
/// never by "the current reference".
pub(crate) trait DepthSink {
    /// The depths at `pos` (0-based) of reference `rid`, one per input in input order.
    fn position(&mut self, rid: usize, pos: u64, depths: &[u32]) -> io::Result<()>;
    /// Positions `[a, b)` of `rid`, where every input's depth is 0.
    fn zeros(&mut self, rid: usize, a: u64, b: u64) -> io::Result<()>;
    /// Every position of `rid` that this run gives has been given.
    fn reference_done(&mut self, _rid: usize) -> io::Result<()> {
        Ok(())
    }
    /// A record of input `input` on `rid` that passed every filter (reader, plan filter, read
    /// length), whether or not it covers a base. `span` is its D101 span cut at `LN`, or `None`
    /// when it is unmapped or the cut span is empty.
    fn read(&mut self, _input: usize, _rid: usize, _a: &AlnRef<'_>, _span: Option<Interval>) {}
    /// A record placed on `rid` (reference and position set), before the plan's filter.
    fn placed(&mut self, _rid: usize) {}
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// How a read's length is taken for `--min-read-length`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ReadLength {
    /// D137, `depth -l`: `l_seq` less the soft clips at the ends.
    SeqLessClips,
    /// The CIGAR's query length (`M I S = X`), as `samtools coverage -l`.
    CigarQuery,
}

/// What is done about the bases of a pair's two reads that overlap.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum MateRule {
    Off,
    /// D139: the second read is clipped at the end of its mate.
    ClipAtMateEnd,
    /// D173, mosdepth's: the second read of a proper pair does not count the positions the
    /// first read's aligned segments cover.
    SharedBases,
}

/// Records dropped before anything else sees them, once the sink has been told they were placed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ReadFilter {
    pub exclude_flags: u16,
    pub min_mapq: u8,
}

impl ReadFilter {
    /// No excluded bit set and MAPQ at least the minimum.
    pub(crate) fn passes(&self, a: &AlnRef<'_>) -> bool {
        a.flag() & self.exclude_flags == 0 && a.mapq() >= self.min_mapq
    }
}

/// What a run of rows is computed from: the options as `depth` normalises them, the
/// references, and the regions merged per reference.
pub(crate) struct Plan<'a> {
    pub rules: CountRules,
    pub min_read_length: u64,
    pub all: AllPositions,
    pub read_length: ReadLength,
    /// What to do about an overlapping mate pair (D139).
    pub mates: MateRule,
    pub filter: Option<ReadFilter>,
    pub refs: &'a [(Vec<u8>, u64)],
    /// The union of the regions on each reference, sorted, disjoint and cut at `LN`.
    pub regions: Option<Vec<Vec<Interval>>>,
    /// The past-`LN` warning has been given.
    pub warned_past: &'a AtomicBool,
}

impl<'a> Plan<'a> {
    pub(crate) fn new(
        o: &DepthOptions,
        mode: ModeRules,
        refs: &'a [(Vec<u8>, u64)],
        warned_past: &'a AtomicBool,
    ) -> Plan<'a> {
        let regions = o.regions.as_ref().map(|list| {
            let by_name: std::collections::HashMap<&[u8], usize> = refs
                .iter()
                .enumerate()
                .map(|(i, (n, _))| (&n[..], i))
                .collect();
            let mut per: Vec<Vec<Interval>> = vec![Vec::new(); refs.len()];
            for (name, iv) in list {
                let Some(&rid) = by_name.get(&name[..]) else {
                    continue;
                };
                let end = iv.end.min(refs[rid].1);
                if iv.start < end {
                    per[rid].push(Interval {
                        start: iv.start,
                        end,
                    });
                }
            }
            for ivs in &mut per {
                ivs.sort_by_key(|iv| (iv.start, iv.end));
                let mut merged: Vec<Interval> = Vec::with_capacity(ivs.len());
                for iv in ivs.drain(..) {
                    match merged.last_mut() {
                        Some(last) if iv.start <= last.end => last.end = last.end.max(iv.end),
                        _ => merged.push(iv),
                    }
                }
                *ivs = merged;
            }
            per
        });
        // With regions, every position in them prints, whatever level of `--all-*` is asked.
        let all = if regions.is_some() && o.all != AllPositions::Covered {
            AllPositions::AllRefs
        } else {
            o.all
        };
        Plan {
            rules: o.rules,
            min_read_length: o.min_read_length,
            all,
            read_length: mode.read_length,
            mates: if mode.shared_bases {
                MateRule::SharedBases
            } else if o.skip_mate_overlap {
                MateRule::ClipAtMateEnd
            } else {
                MateRule::Off
            },
            filter: mode.filter,
            refs,
            regions,
            warned_past,
        }
    }
}

/// Formats rows into a buffer and hands it on every [`FLUSH_AT`] bytes.
pub(crate) struct Emitter<'a> {
    out: &'a mut dyn Write,
    hdr: Arc<AlnHeader>,
    /// Depth columns per row: one per input.
    columns: usize,
    buf: Vec<u8>,
    pub rows: u64,
    /// `CHROM\tPOS` of the row last started.
    head: RowHead,
}

impl<'a> Emitter<'a> {
    pub(crate) fn new(out: &'a mut dyn Write, hdr: &Arc<AlnHeader>, columns: usize) -> Emitter<'a> {
        Emitter {
            out,
            hdr: hdr.clone(),
            columns,
            buf: Vec::with_capacity(FLUSH_AT + 4096),
            rows: 0,
            head: RowHead::default(),
        }
    }

    /// Text that comes before the first row.
    pub(crate) fn preface(&mut self, text: &[u8]) {
        self.buf.extend_from_slice(text);
    }

    fn hand_on(&mut self) -> io::Result<()> {
        if !self.buf.is_empty() {
            self.out.write_all(&self.buf)?;
            self.buf.clear();
        }
        Ok(())
    }

    /// Start the row of the 1-based position `pos1` on the reference `rid`.
    fn row_start(&mut self, rid: usize, pos1: u64) {
        self.head.to(rid, &self.hdr.refs[rid].0, pos1);
        self.buf.extend_from_slice(&self.head.text);
    }

    fn column(&mut self, depth: u32) {
        if depth < 10 {
            self.buf.extend_from_slice(&[b'\t', b'0' + depth as u8]);
        } else if depth < 100 {
            self.buf.extend_from_slice(&[
                b'\t',
                b'0' + (depth / 10) as u8,
                b'0' + (depth % 10) as u8,
            ]);
        } else {
            self.buf.push(b'\t');
            push_u64(&mut self.buf, u64::from(depth));
        }
    }

    fn row_end(&mut self) -> io::Result<()> {
        self.buf.push(b'\n');
        self.rows += 1;
        if self.buf.len() >= FLUSH_AT {
            self.hand_on()?;
        }
        Ok(())
    }
}

impl DepthSink for Emitter<'_> {
    fn position(&mut self, rid: usize, pos: u64, depths: &[u32]) -> io::Result<()> {
        self.row_start(rid, pos + 1);
        for &d in depths {
            self.column(d);
        }
        self.row_end()
    }

    /// The rows of the positions `[a, b)` of `rid`, where every column is 0.
    fn zeros(&mut self, rid: usize, a: u64, b: u64) -> io::Result<()> {
        for p in a..b {
            self.head.to(rid, &self.hdr.refs[rid].0, p + 1);
            self.buf.extend_from_slice(&self.head.text);
            for _ in 0..self.columns {
                self.buf.extend_from_slice(b"\t0");
            }
            self.buf.push(b'\n');
            if self.buf.len() >= FLUSH_AT {
                self.hand_on()?;
            }
        }
        self.rows += b.saturating_sub(a);
        Ok(())
    }

    fn flush(&mut self) -> io::Result<()> {
        self.hand_on()
    }
}

/// The start of a row, `CHROM\tPOS`, kept from row to row: consecutive rows count the
/// position up in place, without writing the name or dividing.
struct RowHead {
    text: Vec<u8>,
    /// Where the digits of the position start in `text`.
    digits_at: usize,
    /// The reference of `text`; `usize::MAX` before the first row.
    rid: usize,
    value: u64,
}

impl Default for RowHead {
    fn default() -> RowHead {
        RowHead {
            text: Vec::new(),
            digits_at: 0,
            rid: usize::MAX,
            value: 0,
        }
    }
}

impl RowHead {
    /// Hold the position `v` (1-based) of the reference `rid`, called `name`.
    fn to(&mut self, rid: usize, name: &[u8], v: u64) {
        if rid != self.rid {
            self.text.clear();
            self.text.extend_from_slice(name);
            self.text.push(b'\t');
            self.digits_at = self.text.len();
            self.rid = rid;
        } else if v == self.value {
            return;
        } else if v == self.value.wrapping_add(1) {
            self.step();
            return;
        }
        self.text.truncate(self.digits_at);
        push_u64(&mut self.text, v);
        self.value = v;
    }

    /// Hold the next position.
    fn step(&mut self) {
        self.value += 1;
        let mut i = self.text.len() - 1;
        loop {
            if self.text[i] < b'9' {
                self.text[i] += 1;
                return;
            }
            self.text[i] = b'0';
            if i == self.digits_at {
                self.text.insert(i, b'1');
                return;
            }
            i -= 1;
        }
    }
}

/// Append `v` in decimal.
pub(crate) fn push_u64(out: &mut Vec<u8>, mut v: u64) {
    let mut tmp = [0u8; 20];
    let mut i = tmp.len();
    loop {
        i -= 1;
        tmp[i] = b'0' + (v % 10) as u8;
        v /= 10;
        if v == 0 {
            break;
        }
    }
    out.extend_from_slice(&tmp[i..]);
}

/// The next record to count, as one input offers it.
#[derive(Clone, Copy)]
struct Cur {
    rid: usize,
    pos: u64,
    span: Interval,
}

/// One input of a run: its records, its counts, and the record it offers next.
pub(crate) struct Input {
    /// `None`: no record of this input belongs to the run.
    src: Option<Box<dyn AlnSource>>,
    pub(crate) name: String,
    layout: gloci_io::Format,
    hdr: Arc<AlnHeader>,
    buf: Vec<u8>,
    pub(crate) ring: DepthRing,
    /// Reads waiting for their mate (`--skip-mate-overlap`).
    mates: MateTable,
    /// Proper pairs waiting for their mate (mosdepth's rule).
    pairs: PairTable,
    /// The segments the current read shares with its mate.
    shared: Vec<Interval>,
    cur: Option<Cur>,
    /// Reference, position and record number of the last record read, for the order check.
    last: Option<(usize, u64, u64)>,
}

impl Input {
    pub(crate) fn new(src: Option<Box<dyn AlnSource>>, name: &str, hdr: Arc<AlnHeader>) -> Input {
        Input {
            layout: src.as_ref().map_or(gloci_io::Format::Bam, |s| s.layout()),
            src,
            name: name.to_string(),
            hdr,
            buf: Vec::new(),
            ring: DepthRing::new(),
            mates: MateTable::new(),
            pairs: PairTable::new(),
            shared: Vec::new(),
            cur: None,
            last: None,
        }
    }

    /// What the reader noted about this input (counts of reads left out, invalid records).
    pub(crate) fn diagnostics(&self) -> Diagnostics {
        self.src
            .as_ref()
            .map(|s| s.diagnostics().clone())
            .unwrap_or_default()
    }

    fn at(&self, rid: usize, pos: u64) -> String {
        format!(
            "{}:{}",
            String::from_utf8_lossy(&self.hdr.refs[rid].0),
            pos + 1
        )
    }

    /// Read to the next record that counts.
    fn advance<S: DepthSink>(
        &mut self,
        k: usize,
        plan: &Plan<'_>,
        sink: &mut S,
    ) -> Result<(), IoError> {
        self.cur = None;
        let Some(src) = self.src.as_mut() else {
            return Ok(());
        };
        loop {
            self.buf.clear();
            let Some(rec) = src.next_record(&mut self.buf)? else {
                return Ok(());
            };
            // The source checked the record.
            let a = aln_view(self.layout, &self.buf, &self.hdr).map_err(|msg| IoError::Parse {
                path: self.name.clone(),
                loc: Loc::Record(rec),
                msg,
            })?;
            let (Some(rid), Some(pos)) = (a.ref_id(), a.pos()) else {
                continue;
            };
            if let Some((prid, ppos, prec)) = self.last
                && (rid, pos) < (prid, ppos)
            {
                return Err(OrderError {
                    path: self.name.clone(),
                    prev_loc: Loc::Record(prec),
                    prev: self.at(prid, ppos),
                    loc: Loc::Record(rec),
                    rec: self.at(rid, pos),
                    hint: format!(
                        "run: gloci sort {} -o sorted.bam, and use that file",
                        self.name
                    ),
                }
                .into());
            }
            self.last = Some((rid, pos, rec));
            sink.placed(rid);
            if let Some(f) = &plan.filter
                && !f.passes(&a)
            {
                continue;
            }
            if plan.min_read_length > 0 {
                let len = match plan.read_length {
                    ReadLength::SeqLessClips => read_length(&a),
                    ReadLength::CigarQuery => query_length(&a),
                };
                if len < plan.min_read_length {
                    continue;
                }
            }
            let ref_len = self.hdr.refs[rid].1;
            let unmapped = a.flag() & FLAG_UNMAPPED != 0;
            let cut = if unmapped {
                None
            } else {
                match a.span(ref_len) {
                    Some((iv, past)) => {
                        if past && !plan.warned_past.swap(true, Ordering::Relaxed) {
                            log::warn!(
                                "{}: a read ends past the end of its reference ({} has LN {}); \
depth counts only the bases up to LN, where samtools prints rows beyond it (this warning is given once)",
                                self.name,
                                String::from_utf8_lossy(&self.hdr.refs[rid].0),
                                ref_len
                            );
                        }
                        (iv.start < iv.end).then_some(iv)
                    }
                    None => None,
                }
            };
            sink.read(k, rid, &a, cut);
            let span = if unmapped {
                Interval {
                    start: pos,
                    end: pos + 1,
                }
            } else {
                let Some(iv) = cut else { continue };
                iv
            };
            self.cur = Some(Cur { rid, pos, span });
            return Ok(());
        }
    }

    fn add_current(&mut self, plan: &Plan<'_>) -> Result<(), IoError> {
        let Some(c) = self.cur else { return Ok(()) };
        // `advance` parsed this record already; the buffer has not changed since.
        let a = aln_view(self.layout, &self.buf, &self.hdr).map_err(|msg| {
            IoError::Internal(format!(
                "{}: a record that parsed once does not parse again: {msg}",
                self.name
            ))
        })?;
        if a.flag() & FLAG_UNMAPPED != 0 {
            // A placed unmapped read covers its one position with nothing counted.
            self.ring.add_empty(c.span);
        } else {
            // Consume or store first, then drop what no later read can claim, so a mate at
            // exactly `mpos` finds its entry.
            match plan.mates {
                MateRule::SharedBases => {
                    self.pairs
                        .shared_with_mate(&a, c.pos, c.span.end, &mut self.shared);
                    self.pairs.advance(c.pos);
                    self.ring.add_except(&a, c.span, &plan.rules, &self.shared);
                }
                MateRule::ClipAtMateEnd => {
                    let clip = self.mates.clip_for(&a, c.pos, c.span.end);
                    self.mates.advance(c.pos);
                    self.ring.add(&a, c.span, &plan.rules, clip);
                }
                MateRule::Off => self.ring.add(&a, c.span, &plan.rules, 0),
            }
        }
        Ok(())
    }
}

/// Say at `-v` what the reader left out of `name`, and warn of invalid records it skipped.
pub(crate) fn log_diagnostics(name: &str, d: &Diagnostics) {
    d.log_alignment_counts(name);
    if let Some(m) = d.summary(name) {
        log::warn!("{m}");
    }
}

/// The read length `--min-read-length` compares: `l_seq` less the soft clips at its two
/// ends, or, with SEQ `*`, the bases the CIGAR puts in the read (`M I = X`).
fn read_length(a: &AlnRef<'_>) -> u64 {
    let seq = a.seq_len() as u64;
    let (mut first, mut last, mut n) = (None, None, 0usize);
    let mut aligned = 0u64;
    a.for_each_op(&mut |op, len| {
        if first.is_none() {
            first = Some((op, len));
        }
        last = Some((op, len));
        n += 1;
        if matches!(op, b'M' | b'I' | b'=' | b'X') {
            aligned += u64::from(len);
        }
    });
    if seq == 0 {
        return aligned;
    }
    let mut len = seq;
    if let Some((b'S', l)) = first {
        len = len.saturating_sub(u64::from(l));
    }
    // samtools takes the last operation's soft clip off only when it is not also the first
    // (a CIGAR of one `S` is clipped once).
    if n > 1
        && let Some((b'S', l)) = last
    {
        len = len.saturating_sub(u64::from(l));
    }
    len
}

/// The query length `samtools coverage -l` compares: the CIGAR's `M I S = X`.
fn query_length(a: &AlnRef<'_>) -> u64 {
    let mut n = 0u64;
    a.for_each_op(&mut |op, len| {
        if matches!(op, b'M' | b'I' | b'S' | b'=' | b'X') {
            n += u64::from(len);
        }
    });
    n
}

/// The positions of the references `first..end`, from `inputs`' records (which belong to those
/// references), given to `sink`. Returns the largest ring any input held and the most reads
/// any input's mate table held at once.
pub(crate) fn run_group<S: DepthSink>(
    inputs: &mut [Input],
    first: usize,
    end: usize,
    plan: &Plan<'_>,
    sink: &mut S,
) -> Result<(usize, usize), IoError> {
    let write = |e: io::Error| IoError::Write {
        path: "output".into(),
        source: e,
    };
    let mut max_ring = 0usize;
    let mut max_mates = 0usize;
    let mut printer = Printer {
        plan,
        sink,
        depths: vec![0; inputs.len()],
    };
    for (k, i) in inputs.iter_mut().enumerate() {
        i.advance(k, plan, printer.sink)?;
    }
    let mut cur_ref: Option<usize> = None;
    // The first position of `cur_ref` not yet printed.
    let mut next = 0u64;
    loop {
        // The input with the smallest (reference, position); the lower index wins ties.
        let mut pick: Option<(usize, Cur)> = None;
        for (k, i) in inputs.iter().enumerate() {
            if let Some(c) = i.cur
                && pick.is_none_or(|(_, p)| (c.rid, c.pos) < (p.rid, p.pos))
            {
                pick = Some((k, c));
            }
        }
        let Some((k, c)) = pick else { break };
        if cur_ref != Some(c.rid) {
            if let Some(r) = cur_ref {
                printer
                    .pending(inputs, r, next, plan.refs[r].1)
                    .map_err(write)?;
                printer.sink.reference_done(r).map_err(write)?;
            }
            if plan.all == AllPositions::AllRefs {
                for r in cur_ref.map_or(first, |r| r + 1)..c.rid {
                    printer.whole(inputs, r).map_err(write)?;
                    printer.sink.reference_done(r).map_err(write)?;
                }
            }
            cur_ref = Some(c.rid);
            next = 0;
            for i in inputs.iter_mut() {
                i.ring.reset();
                i.mates.reset();
                i.pairs.reset();
            }
        }
        printer.pending(inputs, c.rid, next, c.pos).map_err(write)?;
        next = next.max(c.pos);
        inputs[k].add_current(plan)?;
        max_ring = max_ring.max(inputs[k].ring.capacity());
        max_mates = max_mates
            .max(inputs[k].mates.peak())
            .max(inputs[k].pairs.peak());
        inputs[k].advance(k, plan, printer.sink)?;
    }
    if let Some(r) = cur_ref {
        printer
            .pending(inputs, r, next, plan.refs[r].1)
            .map_err(write)?;
        printer.sink.reference_done(r).map_err(write)?;
    }
    if plan.all == AllPositions::AllRefs {
        for r in cur_ref.map_or(first, |r| r + 1)..end {
            printer.whole(inputs, r).map_err(write)?;
            printer.sink.reference_done(r).map_err(write)?;
        }
    }
    printer.sink.flush().map_err(write)?;
    Ok((max_ring, max_mates))
}

struct Printer<'p, 's, S: DepthSink> {
    plan: &'p Plan<'p>,
    sink: &'s mut S,
    /// One depth per input, for the position being given.
    depths: Vec<u32>,
}

impl<S: DepthSink> Printer<'_, '_, S> {
    /// A whole reference no read of the run is on.
    fn whole(&mut self, inputs: &mut [Input], rid: usize) -> io::Result<()> {
        for i in inputs.iter_mut() {
            i.ring.reset();
        }
        self.range(inputs, rid, 0, self.plan.refs[rid].1)
    }

    /// The positions `[a, b)` of `rid` that print now: all of them under `--all-*`, else those
    /// some input's reads cover.
    fn pending(&mut self, inputs: &mut [Input], rid: usize, a: u64, b: u64) -> io::Result<()> {
        let hi = match self.plan.all {
            AllPositions::Covered => inputs
                .iter()
                .map(|i| i.ring.covered_end())
                .max()
                .unwrap_or(0)
                .min(b),
            _ => b,
        };
        self.range(inputs, rid, a, hi)
    }

    /// Print `[a, b)` of `rid` where the regions allow.
    fn range(&mut self, inputs: &mut [Input], rid: usize, a: u64, b: u64) -> io::Result<()> {
        let b = b.min(self.plan.refs[rid].1);
        if a >= b {
            return Ok(());
        }
        let Some(regions) = &self.plan.regions else {
            return self.positions(inputs, rid, a, b);
        };
        let ivs = &regions[rid];
        let from = ivs.partition_point(|iv| iv.end <= a);
        for iv in &ivs[from..] {
            if iv.start >= b {
                break;
            }
            self.positions(inputs, rid, a.max(iv.start), b.min(iv.end))?;
        }
        Ok(())
    }

    fn positions(&mut self, inputs: &mut [Input], rid: usize, a: u64, b: u64) -> io::Result<()> {
        // At and past every input's covered end, every depth is 0.
        let covered = inputs
            .iter()
            .map(|i| i.ring.covered_end())
            .max()
            .unwrap_or(0)
            .clamp(a, b);
        for p in a..covered {
            for (d, i) in self.depths.iter_mut().zip(inputs.iter_mut()) {
                *d = i.ring.take(p);
            }
            self.sink.position(rid, p, &self.depths)?;
        }
        if covered < b {
            self.sink.zeros(rid, covered, b)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::RowHead;

    #[test]
    fn row_head_counts_up_and_jumps() {
        let mut h = RowHead::default();
        let steps = (1..=12)
            .chain(98..=101)
            .chain([999, 1000, 7, 8, 8, 1_000_000_000_000]);
        for v in steps {
            h.to(0, b"chr1", v);
            assert_eq!(h.text, format!("chr1\t{v}").as_bytes(), "after to({v})");
        }
        h.to(3, b"chrX", 8);
        assert_eq!(h.text, b"chrX\t8");
        h.to(3, b"chrX", 9);
        assert_eq!(h.text, b"chrX\t9");
        h.to(0, b"chr1", 10);
        assert_eq!(h.text, b"chr1\t10");
        h.to(0, b"chr1", u64::MAX - 1);
        h.to(0, b"chr1", u64::MAX);
        assert_eq!(h.text, format!("chr1\t{}", u64::MAX).as_bytes());
    }
}

#[cfg(test)]
mod sink_tests {
    use super::super::test_sink::{Event, Recorder, sam_source};
    use super::super::{AllPositions, DepthOptions, ModeRules, run_sources};
    use super::{ReadFilter, ReadLength};
    use crate::depth::CountRules;
    use gloci_core::{Interval, VcfSpan};
    use gloci_io::{Accept, AlignOptions, InvalidPolicy, ReaderOptions};
    use std::num::NonZeroUsize;

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

    fn opts(all: AllPositions, min_read_length: u64) -> DepthOptions {
        DepthOptions {
            rules: CountRules::default(),
            min_read_length,
            all,
            regions: None,
            skip_mate_overlap: false,
            column_names: false,
        }
    }

    /// A SAM of `refs` (name, LN) and `reads` (name flag ref pos mapq cigar).
    fn sam(refs: &[(&str, u64)], reads: &[&str]) -> String {
        let mut t = String::from("@HD\tVN:1.6\tSO:coordinate\n");
        for (n, l) in refs {
            t.push_str(&format!("@SQ\tSN:{n}\tLN:{l}\n"));
        }
        for r in reads {
            let f: Vec<&str> = r.split(' ').collect();
            t.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t*\t*\n",
                f[0], f[1], f[2], f[3], f[4], f[5]
            ));
        }
        t
    }

    fn run(text: &str, o: &DepthOptions, mode: ModeRules) -> Recorder {
        let src = sam_source(text, &ro());
        let (_, rec) = run_sources(vec![src], &["test.sam".to_string()], o, mode, &ro(), |_| {
            Ok(Recorder::default())
        })
        .unwrap();
        rec
    }

    #[test]
    fn sink_gets_every_position_once_in_order() {
        let text = sam(
            &[("c1", 20), ("c2", 5), ("c3", 8)],
            &["r1 0 c1 3 60 3M2D2M", "r2 0 c1 4 60 4M", "r3 0 c3 2 60 2M"],
        );
        let rec = run(&text, &opts(AllPositions::AllRefs, 0), ModeRules::depth());
        // Every position once, in order: the runs tile each reference from 0 to LN.
        let mut expect_next = (0usize, 0u64);
        for &(rid, a, b, _) in &rec.runs {
            assert_eq!((rid, a), expect_next, "runs: {:?}", rec.runs);
            expect_next = if b == [20, 5, 8][rid] {
                (rid + 1, 0)
            } else {
                (rid, b)
            };
        }
        assert_eq!(expect_next, (3, 0), "runs: {:?}", rec.runs);
        let got: Vec<u32> = (2..9).map(|p| rec.depth_at(0, p)).collect();
        assert_eq!(got, [1, 2, 2, 1, 1, 1, 1]);
        let done: Vec<usize> = rec
            .events
            .iter()
            .filter_map(|e| {
                if let Event::Done(r) = e {
                    Some(*r)
                } else {
                    None
                }
            })
            .collect();
        assert_eq!(done, [0, 1, 2]);
        // Each reference's done comes after its last run.
        for (rid, &n) in rec.runs_at_done.iter().enumerate() {
            assert_eq!(rec.runs[n - 1].0, rid);
            assert_eq!(rec.runs[n - 1].2, [20, 5, 8][rid]);
        }
    }

    #[test]
    fn read_hook_sees_each_passing_record() {
        let refs = [("c1", 100), ("c2", 50)];
        let reads = |text: &str, len: u64, mode: ModeRules| {
            let rec = run(text, &opts(AllPositions::Covered, len), mode);
            rec.events
                .into_iter()
                .filter_map(|e| match e {
                    Event::Read { name, span, .. } => Some((name, span)),
                    _ => None,
                })
                .collect::<Vec<_>>()
        };
        let iv = |start, end| Some(Interval { start, end });
        let t = sam(&refs, &["m 0 c1 3 60 45M", "u 4 c1 7 0 *"]);
        assert_eq!(
            reads(&t, 0, ModeRules::depth()),
            [("m".to_string(), iv(2, 47)), ("u".to_string(), None)]
        );
        let t = sam(&refs, &["long 0 c1 1 60 30S30M30S", "noseq 0 c1 1 60 50M"]);
        let coverage = ModeRules {
            read_length: ReadLength::CigarQuery,
            filter: None,
            shared_bases: false,
        };
        assert_eq!(reads(&t, 60, coverage), [("long".to_string(), iv(0, 30))]);
        assert_eq!(reads(&t, 60, ModeRules::depth()), []);
    }

    #[test]
    fn read_hook_sees_mapped_reads_with_an_empty_span() {
        // `a` is mapped with CIGAR `*`; `b` starts at LN + 1 (0-based LN), so the cut at LN
        // leaves nothing. Neither covers a base, and both are still reads.
        let text = sam(&[("c1", 100)], &["a 0 c1 3 60 *", "b 0 c1 101 60 5M"]);
        let rec = run(&text, &opts(AllPositions::Covered, 0), ModeRules::depth());
        let reads: Vec<_> = rec
            .events
            .iter()
            .filter_map(|e| match e {
                Event::Read { name, span, .. } => Some((name.as_str(), *span)),
                _ => None,
            })
            .collect();
        assert_eq!(reads, [("a", None), ("b", None)]);
    }

    #[test]
    fn placed_hook_runs_before_the_plan_filter() {
        let text = sam(
            &[("c1", 100), ("c2", 50)],
            &["d 1024 c2 1 60 5M", "q 0 c2 3 5 5M"],
        );
        let mode = ModeRules {
            read_length: ReadLength::SeqLessClips,
            filter: Some(ReadFilter {
                exclude_flags: 0x400,
                min_mapq: 20,
            }),
            shared_bases: false,
        };
        let rec = run(&text, &opts(AllPositions::Covered, 0), mode);
        let placed = rec
            .events
            .iter()
            .filter(|e| **e == Event::Placed(1))
            .count();
        assert_eq!(placed, 2);
        assert!(!rec.events.iter().any(|e| matches!(e, Event::Read { .. })));
    }
}
