//! The rows of `depth`: which positions print, in which order, for one run of references
//! (D138). One run reads one record stream per input, merged by `(reference, position)`.

use super::ring::{CountRules, DepthRing};
use super::{AllPositions, DepthOptions};
use gloci_core::Interval;
use gloci_io::{AlnHeader, AlnRef, AlnSource, Diagnostics, IoError, Loc, OrderError, aln_ref};
use std::io;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

const FLAG_UNMAPPED: u16 = 0x4;
/// Bytes of rows gathered before they are handed on.
const FLUSH_AT: usize = 1 << 16;

/// What a run of rows is computed from: the options as `depth` normalises them, the
/// references, and the regions merged per reference.
pub(crate) struct Plan<'a> {
    pub rules: CountRules,
    pub min_read_length: u64,
    pub all: AllPositions,
    pub refs: &'a [(Vec<u8>, u64)],
    /// The union of the regions on each reference, sorted, disjoint and cut at `LN`.
    pub regions: Option<Vec<Vec<Interval>>>,
    /// The past-`LN` warning has been given.
    pub warned_past: &'a AtomicBool,
}

impl<'a> Plan<'a> {
    pub(crate) fn new(
        o: &DepthOptions,
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
            refs,
            regions,
            warned_past,
        }
    }
}

/// Formats rows into a buffer and hands it on every 64 KiB.
pub(crate) struct Emitter<'a> {
    sink: &'a mut dyn FnMut(&[u8]) -> io::Result<()>,
    buf: Vec<u8>,
    pub rows: u64,
}

impl<'a> Emitter<'a> {
    pub(crate) fn new(sink: &'a mut dyn FnMut(&[u8]) -> io::Result<()>) -> Emitter<'a> {
        Emitter {
            sink,
            buf: Vec::with_capacity(FLUSH_AT + 4096),
            rows: 0,
        }
    }

    pub(crate) fn flush(&mut self) -> io::Result<()> {
        if !self.buf.is_empty() {
            (self.sink)(&self.buf)?;
            self.buf.clear();
        }
        Ok(())
    }

    fn row_start(&mut self, name: &[u8], pos1: u64) {
        self.buf.extend_from_slice(name);
        self.buf.push(b'\t');
        push_u64(&mut self.buf, pos1);
    }

    fn column(&mut self, depth: u32) {
        self.buf.push(b'\t');
        push_u64(&mut self.buf, u64::from(depth));
    }

    fn row_end(&mut self) -> io::Result<()> {
        self.buf.push(b'\n');
        self.rows += 1;
        if self.buf.len() >= FLUSH_AT {
            self.flush()?;
        }
        Ok(())
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
    fn advance(&mut self, plan: &Plan<'_>) -> Result<(), IoError> {
        self.cur = None;
        let Some(src) = self.src.as_mut() else {
            return Ok(());
        };
        loop {
            self.buf.clear();
            let Some(rec) = src.next_record(&mut self.buf)? else {
                return Ok(());
            };
            let a = aln_ref(self.layout, &self.buf, &self.hdr).map_err(|msg| IoError::Parse {
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
            if plan.min_read_length > 0 && read_length(&a) < plan.min_read_length {
                continue;
            }
            let ref_len = self.hdr.refs[rid].1;
            let span = if a.flag() & FLAG_UNMAPPED != 0 {
                Interval {
                    start: pos,
                    end: pos + 1,
                }
            } else {
                match a.span(ref_len) {
                    Some((iv, false)) => iv,
                    Some((iv, true)) => {
                        if !plan.warned_past.swap(true, Ordering::Relaxed) {
                            log::warn!(
                                "{}: a read ends past the end of its reference ({} has LN {}); \
depth counts only the bases up to LN, where samtools prints rows beyond it (this warning is given once)",
                                self.name,
                                String::from_utf8_lossy(&self.hdr.refs[rid].0),
                                ref_len
                            );
                        }
                        iv
                    }
                    None => continue,
                }
            };
            if span.start >= span.end {
                continue;
            }
            self.cur = Some(Cur { rid, pos, span });
            return Ok(());
        }
    }

    fn add_current(&mut self, plan: &Plan<'_>) -> Result<(), IoError> {
        let Some(c) = self.cur else { return Ok(()) };
        // `advance` parsed this record already; the buffer has not changed since.
        let a = aln_ref(self.layout, &self.buf, &self.hdr).map_err(|msg| {
            IoError::Internal(format!(
                "{}: a record that parsed once does not parse again: {msg}",
                self.name
            ))
        })?;
        if a.flag() & FLAG_UNMAPPED != 0 {
            // A placed unmapped read covers its one position with nothing counted.
            self.ring.add_empty(c.span);
        } else {
            self.ring.add(&a, c.span, &plan.rules, 0);
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
    if n > 1
        && let Some((b'S', l)) = last
    {
        len = len.saturating_sub(u64::from(l));
    }
    len
}

/// The rows of the references `first..end`, from `inputs`' records (which belong to those
/// references), written to `em`. Returns the largest ring any input held.
pub(crate) fn run_group(
    inputs: &mut [Input],
    first: usize,
    end: usize,
    plan: &Plan<'_>,
    em: &mut Emitter<'_>,
) -> Result<usize, IoError> {
    let write = |e: io::Error| IoError::Write {
        path: "output".into(),
        source: e,
    };
    let mut max_ring = 0usize;
    for i in inputs.iter_mut() {
        i.advance(plan)?;
    }
    let mut cur_ref: Option<usize> = None;
    // The first position of `cur_ref` not yet printed.
    let mut next = 0u64;
    let mut printer = Printer { plan, em };
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
            }
            if plan.all == AllPositions::AllRefs {
                for r in cur_ref.map_or(first, |r| r + 1)..c.rid {
                    printer.whole(inputs, r).map_err(write)?;
                }
            }
            cur_ref = Some(c.rid);
            next = 0;
            for i in inputs.iter_mut() {
                i.ring.reset();
            }
        }
        printer.pending(inputs, c.rid, next, c.pos).map_err(write)?;
        next = next.max(c.pos);
        inputs[k].add_current(plan)?;
        max_ring = max_ring.max(inputs[k].ring.capacity());
        inputs[k].advance(plan)?;
    }
    if let Some(r) = cur_ref {
        printer
            .pending(inputs, r, next, plan.refs[r].1)
            .map_err(write)?;
    }
    if plan.all == AllPositions::AllRefs {
        for r in cur_ref.map_or(first, |r| r + 1)..end {
            printer.whole(inputs, r).map_err(write)?;
        }
    }
    printer.em.flush().map_err(write)?;
    Ok(max_ring)
}

struct Printer<'p, 'a, 'e> {
    plan: &'p Plan<'p>,
    em: &'e mut Emitter<'a>,
}

impl Printer<'_, '_, '_> {
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
        let name = &self.plan.refs[rid].0;
        for p in a..b {
            self.em.row_start(name, p + 1);
            for i in inputs.iter_mut() {
                self.em.column(i.ring.take(p));
            }
            self.em.row_end()?;
        }
        Ok(())
    }
}
