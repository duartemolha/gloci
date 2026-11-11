//! `pair-intersect` (D154, bedtools pairtobed): the pairs of A, by how their
//! ends overlap the features of B.
//!
//! B goes into an [`IntervalIndex`]; A streams in batches, each batch is
//! worked on by one of `threads` workers, and the finished batches are
//! written in A's order, so the output does not depend on the thread count.
//! At most `2 * threads` batches are in flight.

use crate::intersect::Text;
use crate::skipped::log_reader;
use gloci_core::{ChromId, ChromTable, Interval, OverlapCriteria, Strand, StrandRule, thousands};
use gloci_engine::{IndexEntry, IntervalIndex};
use gloci_io::{
    BedpeBatch, BedpeEnd, BedpeMeta, BedpeReader, Format, Header, IoError, RecordBatch,
    RecordReader, UsageError, column_count,
};
use std::collections::BTreeMap;
use std::io::Write;
use std::num::NonZeroUsize;
use std::sync::mpsc::{channel, sync_channel};
use std::sync::{Arc, Mutex};

/// What `pair-intersect` reports (bedtools `-type`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairType {
    /// A with each B overlapping end 1, then each overlapping end 2.
    Either,
    /// A, alone, when neither end overlaps B.
    Neither,
    /// As `Either`, when both ends overlap B.
    Both,
    /// As `Either`, when exactly one end overlaps B.
    Xor,
    /// A alone (padded with a placeholder B) when no end overlaps B; as
    /// `Either` when exactly one does.
    NotBoth,
    /// A with each B overlapping the span between the ends.
    InnerSpan,
    /// A with each B overlapping the span from the first start to the last end.
    OuterSpan,
    /// A, alone, when nothing overlaps its inner span.
    NotInnerSpan,
    /// A, alone, when nothing overlaps its outer span.
    NotOuterSpan,
}

impl PairType {
    /// The spans of a pair are the unit, not its ends.
    pub fn is_span(self) -> bool {
        matches!(
            self,
            PairType::InnerSpan
                | PairType::OuterSpan
                | PairType::NotInnerSpan
                | PairType::NotOuterSpan
        )
    }

    /// Whether the report is A's pairs alone, which are BEDPE lines.
    pub fn writes_a_alone(self) -> bool {
        matches!(
            self,
            PairType::Neither | PairType::NotInnerSpan | PairType::NotOuterSpan
        )
    }

    pub fn name(self) -> &'static str {
        match self {
            PairType::Either => "either",
            PairType::Neither => "neither",
            PairType::Both => "both",
            PairType::Xor => "xor",
            PairType::NotBoth => "notboth",
            PairType::InnerSpan => "ispan",
            PairType::OuterSpan => "ospan",
            PairType::NotInnerSpan => "notispan",
            PairType::NotOuterSpan => "notospan",
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct PairIntersectOptions {
    pub kind: PairType,
    pub criteria: OverlapCriteria,
}

/// What a run saw.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PairSummary {
    /// Pairs read from A.
    pub a_pairs: u64,
    /// Pairs a span type left out because an end is absent (one or both) or
    /// the ends are on different chromosomes.
    pub span_skipped: u64,
}

/// B's records after loading: the index, and each record's text and shape.
struct BStore {
    index: IntervalIndex,
    /// Per id: span and strand.
    meta: Vec<(Interval, Strand)>,
    /// Per id: offset and length of its text in `text`.
    at: Vec<(usize, u32)>,
    /// B as text, one record after another (a BCF record as VCF).
    text: Vec<u8>,
    /// The row written for a missing B (`notboth`).
    placeholder: Vec<u8>,
    unknown_strand: bool,
}

fn usage(msg: impl Into<String>) -> IoError {
    UsageError(msg.into()).into()
}

/// Columns of a missing B's placeholder, from B file's format and first record.
fn width(format: Format, header: &Header, first: Option<&[u8]>) -> usize {
    match format {
        Format::Gff3 | Format::Gtf => 9,
        Format::Vcf | Format::Bcf => match header {
            Header::Vcf(h) if !h.samples.is_empty() => 9 + h.samples.len(),
            _ => 8,
        },
        _ => first.map_or(3, column_count),
    }
}

fn load_b(mut b: Box<dyn RecordReader>, chroms: &mut ChromTable) -> Result<BStore, IoError> {
    let mut text = Text::new(b.as_ref())?;
    let mut store = BStore {
        index: IntervalIndex::default(),
        meta: Vec::new(),
        at: Vec::new(),
        text: Vec::new(),
        placeholder: Vec::new(),
        unknown_strand: false,
    };
    let mut entries = Vec::new();
    let mut first: Option<Vec<u8>> = None;
    let (mut batch, mut max) = (RecordBatch::default(), 1);
    loop {
        batch.clear();
        let more = b.read_batch(chroms, &mut batch, max)?;
        max = gloci_io::next_batch_size(&batch, max);
        for (i, m) in batch.recs.iter().enumerate() {
            let id = u32::try_from(store.meta.len())
                .map_err(|_| usage("B has more than 4 billion records; split it"))?;
            let line = text.of(batch.raw(i), m.loc)?;
            if first.is_none() {
                first = Some(line.to_vec());
            }
            store.at.push((store.text.len(), line.len() as u32));
            store.text.extend_from_slice(line);
            store.meta.push((m.span, m.strand));
            store.unknown_strand |= m.strand == Strand::Unknown;
            entries.push(IndexEntry {
                chrom: m.chrom,
                span: m.span,
                id,
            });
        }
        if !more {
            break;
        }
    }
    log_reader(b.as_ref());
    store.index = IntervalIndex::build(entries);
    let w = width(b.format(), b.header(), first.as_deref());
    crate::relate::b_placeholder(b.format(), w, &mut store.placeholder);
    Ok(store)
}

#[derive(Default)]
struct Counts {
    pairs: u64,
    span_skipped: u64,
    /// An end of unknown strand met a strand rule.
    a_unknown: bool,
}

#[derive(Default)]
struct Scratch {
    cand: Vec<u32>,
    h1: Vec<u32>,
    h2: Vec<u32>,
}

struct Worker<'a> {
    b: &'a BStore,
    kind: PairType,
    criteria: &'a OverlapCriteria,
}

impl Worker<'_> {
    /// The ids of the B that meet `span` of `chrom` under the criteria.
    fn hits(&self, chrom: ChromId, span: Interval, strand: Strand, s: &mut Scratch, out: Select) {
        let c = self.criteria;
        self.b
            .index
            .query(chrom, span, None, strand, c.inclusive_points, &mut s.cand);
        let dst = match out {
            Select::One => &mut s.h1,
            Select::Two => &mut s.h2,
        };
        dst.clear();
        dst.extend(s.cand.iter().copied().filter(|&id| {
            let (bs, bstrand) = self.b.meta[id as usize];
            c.accepts(span, strand, bs, bstrand)
        }));
    }

    fn ab(&self, a: &[u8], id: u32, out: &mut Vec<u8>) {
        let (off, len) = self.b.at[id as usize];
        out.extend_from_slice(a);
        out.push(b'\t');
        out.extend_from_slice(&self.b.text[off..off + len as usize]);
        out.push(b'\n');
    }

    fn alone(a: &[u8], out: &mut Vec<u8>) {
        out.extend_from_slice(a);
        out.push(b'\n');
    }

    fn pair(&self, raw: &[u8], m: &BedpeMeta, s: &mut Scratch, out: &mut Vec<u8>, n: &mut Counts) {
        if self.kind.is_span() {
            return self.span(raw, m, s, out, n);
        }
        let strand_rule = self.criteria.strand != StrandRule::Any;
        for (k, end) in m.ends.iter().enumerate() {
            let dst = if k == 0 { Select::One } else { Select::Two };
            match end {
                Some(e) => {
                    n.a_unknown |= strand_rule && e.strand == Strand::Unknown;
                    self.hits(e.chrom, e.span, e.strand, s, dst);
                }
                None => match dst {
                    Select::One => s.h1.clear(),
                    Select::Two => s.h2.clear(),
                },
            }
        }
        let (h1, h2) = (&s.h1, &s.h2);
        let rows = |ids: &[&Vec<u32>], out: &mut Vec<u8>| {
            for ids in ids {
                for &id in *ids {
                    self.ab(raw, id, out);
                }
            }
        };
        match self.kind {
            PairType::Either => rows(&[h1, h2], out),
            PairType::Both if !h1.is_empty() && !h2.is_empty() => rows(&[h1, h2], out),
            PairType::Xor if h1.is_empty() != h2.is_empty() => rows(&[h1, h2], out),
            PairType::Neither if h1.is_empty() && h2.is_empty() => Self::alone(raw, out),
            PairType::NotBoth => {
                if h1.is_empty() && h2.is_empty() {
                    out.extend_from_slice(raw);
                    out.push(b'\t');
                    out.extend_from_slice(&self.b.placeholder);
                    out.push(b'\n');
                } else if h1.is_empty() != h2.is_empty() {
                    rows(&[h1, h2], out);
                }
            }
            _ => {}
        }
    }

    /// The inner and outer spans of a pair on one chromosome. The inner span
    /// is the gap between the ends, empty when they touch or overlap.
    fn spans(e1: &BedpeEnd, e2: &BedpeEnd) -> (Option<Interval>, Interval) {
        let inner = if e1.span.end <= e2.span.start {
            Some((e1.span.end, e2.span.start))
        } else if e2.span.end <= e1.span.start {
            Some((e2.span.end, e1.span.start))
        } else {
            None
        };
        let inner = inner
            .filter(|(s, e)| s < e)
            .map(|(start, end)| Interval { start, end });
        let outer = Interval {
            start: e1.span.start.min(e2.span.start),
            end: e1.span.end.max(e2.span.end),
        };
        (inner, outer)
    }

    fn span(&self, raw: &[u8], m: &BedpeMeta, s: &mut Scratch, out: &mut Vec<u8>, n: &mut Counts) {
        let [Some(e1), Some(e2)] = &m.ends else {
            // A span type needs both ends on one chromosome: an absent end (chromosome `.`)
            // leaves the pair out, whether one end is absent or both.
            n.span_skipped += 1;
            return;
        };
        if e1.chrom != e2.chrom {
            n.span_skipped += 1;
            return;
        }
        let (inner, outer) = Self::spans(e1, e2);
        let (inside, span) = match self.kind {
            PairType::InnerSpan | PairType::NotInnerSpan => (inner.is_some(), inner),
            _ => (true, Some(outer)),
        };
        s.h1.clear();
        if let (true, Some(span)) = (inside, span) {
            self.hits(e1.chrom, span, Strand::Unknown, s, Select::One);
        }
        match self.kind {
            PairType::InnerSpan | PairType::OuterSpan => {
                for &id in &s.h1 {
                    self.ab(raw, id, out);
                }
            }
            _ if s.h1.is_empty() => Self::alone(raw, out),
            _ => {}
        }
    }

    fn batch(&self, batch: &BedpeBatch, out: &mut Vec<u8>) -> Counts {
        let mut n = Counts::default();
        let mut s = Scratch::default();
        for (i, m) in batch.recs.iter().enumerate() {
            n.pairs += 1;
            self.pair(batch.raw(i), m, &mut s, out, &mut n);
        }
        n
    }
}

#[derive(Clone, Copy)]
enum Select {
    One,
    Two,
}

/// The result of one batch, by sequence number.
struct Done {
    out: Vec<u8>,
    counts: Counts,
}

fn werr(source: std::io::Error) -> IoError {
    IoError::Write {
        path: "output".to_string(),
        source,
    }
}

/// The option checks that need no input, so a command can run them before it
/// creates its output.
pub fn validate(o: &PairIntersectOptions) -> Result<(), IoError> {
    if o.kind.is_span() && o.criteria.strand != StrandRule::Any {
        return Err(usage(format!(
            "-s and -S cannot be used with --type {}: a span is not on a strand; drop -s/-S",
            o.kind.name()
        )));
    }
    Ok(())
}

/// Run `pair-intersect`: A's pairs against B, written to `out` (TSV rows for
/// A+B, BEDPE lines for A alone). `chroms` interns both inputs' chromosomes.
/// With `threads > 1` the output is byte-identical to `threads == 1`.
pub fn pair_intersect(
    mut a: BedpeReader,
    b: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    o: &PairIntersectOptions,
    threads: NonZeroUsize,
    out: &mut dyn Write,
) -> Result<PairSummary, IoError> {
    validate(o)?;
    let store = load_b(b, chroms)?;
    let worker = Worker {
        b: &store,
        kind: o.kind,
        criteria: &o.criteria,
    };
    let mut summary = PairSummary::default();
    let mut a_unknown = false;
    let mut take = |c: &Counts, summary: &mut PairSummary| {
        summary.a_pairs += c.pairs;
        summary.span_skipped += c.span_skipped;
        a_unknown |= c.a_unknown;
    };
    let (mut batch, mut max) = (BedpeBatch::default(), 1usize);
    let mut next = |a: &mut BedpeReader,
                    chroms: &mut ChromTable,
                    batch: &mut BedpeBatch|
     -> Result<bool, IoError> {
        batch.clear();
        let more = a.read_batch(chroms, batch, max)?;
        if !batch.is_empty() {
            let per_rec = batch.buf.len() / batch.len() + 1;
            max = ((4usize << 20) / per_rec).clamp(1, 16384);
        }
        Ok(more)
    };
    if threads.get() == 1 {
        let mut buf = Vec::new();
        loop {
            let more = next(&mut a, chroms, &mut batch)?;
            buf.clear();
            let c = worker.batch(&batch, &mut buf);
            take(&c, &mut summary);
            out.write_all(&buf).map_err(werr)?;
            if !more {
                break;
            }
        }
    } else {
        parallel(&mut a, chroms, &worker, threads, out, &mut next, &mut |c| {
            take(c, &mut summary)
        })?;
    }
    log_reader_bedpe(&a);
    warn(o, &store, a.name(), a_unknown, &summary);
    Ok(summary)
}

fn log_reader_bedpe(a: &BedpeReader) {
    if let Some(s) = a.diagnostics().summary(a.name()) {
        log::warn!("{s}");
    }
}

fn warn(o: &PairIntersectOptions, store: &BStore, a_name: &str, a_unknown: bool, s: &PairSummary) {
    if o.criteria.strand != StrandRule::Any {
        if a_unknown {
            log::warn!("{a_name} has pair ends of unknown strand; they never match under -s or -S");
        }
        if store.unknown_strand {
            log::warn!("B has records with unknown strand; they never match under -s or -S");
        }
    }
    if s.span_skipped > 0 {
        log::warn!(
            "{} pairs were left out of --type {}: an end is absent (chromosome .) or the ends \
are on different chromosomes",
            thousands(s.span_skipped),
            o.kind.name()
        );
    }
}

type NextBatch<'f> =
    dyn FnMut(&mut BedpeReader, &mut ChromTable, &mut BedpeBatch) -> Result<bool, IoError> + 'f;

/// Batches of A across workers, written in order; at most `2 * threads`
/// batches in flight.
fn parallel(
    a: &mut BedpeReader,
    chroms: &mut ChromTable,
    worker: &Worker<'_>,
    threads: NonZeroUsize,
    out: &mut dyn Write,
    next: &mut NextBatch<'_>,
    take: &mut dyn FnMut(&Counts),
) -> Result<(), IoError> {
    let n = threads.get();
    let limit = 2 * n;
    std::thread::scope(|scope| {
        let (job_tx, job_rx) = sync_channel::<(u64, BedpeBatch)>(limit);
        let job_rx = Arc::new(Mutex::new(job_rx));
        let (done_tx, done_rx) = channel::<(u64, Done)>();
        for _ in 0..n {
            let (job_rx, done_tx) = (job_rx.clone(), done_tx.clone());
            scope.spawn(move || {
                loop {
                    let job = job_rx.lock().unwrap_or_else(|e| e.into_inner()).recv();
                    let Ok((seq, batch)) = job else { break };
                    let mut buf = Vec::new();
                    let counts = worker.batch(&batch, &mut buf);
                    if done_tx.send((seq, Done { out: buf, counts })).is_err() {
                        break;
                    }
                }
            });
        }
        drop(done_tx);
        let mut pending: BTreeMap<u64, Done> = BTreeMap::new();
        let (mut sent, mut written) = (0u64, 0u64);
        let mut flush = |pending: &mut BTreeMap<u64, Done>, written: &mut u64| {
            while let Some(d) = pending.remove(written) {
                take(&d.counts);
                out.write_all(&d.out).map_err(werr)?;
                *written += 1;
            }
            Ok::<(), IoError>(())
        };
        let mut result = Ok(());
        loop {
            let mut batch = BedpeBatch::default();
            let more = match next(a, chroms, &mut batch) {
                Ok(m) => m,
                Err(e) => {
                    result = Err(e);
                    break;
                }
            };
            while sent - written >= limit as u64 && result.is_ok() {
                match done_rx.recv() {
                    Ok((seq, d)) => {
                        pending.insert(seq, d);
                        result = flush(&mut pending, &mut written);
                    }
                    Err(_) => {
                        result = Err(IoError::Internal("a pair-intersect worker died".into()));
                    }
                }
            }
            if result.is_err() {
                break;
            }
            if !batch.is_empty() {
                if job_tx.send((sent, batch)).is_err() {
                    result = Err(IoError::Internal("a pair-intersect worker died".into()));
                    break;
                }
                sent += 1;
            }
            if !more {
                break;
            }
        }
        drop(job_tx);
        while result.is_ok() && written < sent {
            match done_rx.recv() {
                Ok((seq, d)) => {
                    pending.insert(seq, d);
                    result = flush(&mut pending, &mut written);
                }
                Err(_) => {
                    result = Err(IoError::Internal("a pair-intersect worker died".into()));
                }
            }
        }
        result
    })
}
