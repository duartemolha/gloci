//! `gloci depth --summary`: one row per reference, as `samtools coverage` 1.21 prints it (D173).
//! A [`DepthSink`] over the depth core: the positions give covered bases and the depth sum, the
//! `read` hook gives the reads, their MAPQ and the base qualities.

use super::rows::{DepthSink, ReadLength};
use super::{AllPositions, CountRules, DepthOptions, DepthSummary, ModeRules, run_sources};
use gloci_core::Interval;
use gloci_core::cfmt::{write_g, write_u64};
use gloci_io::{AlnHeader, AlnRef, AlnSource, IoError, ReaderOptions, UsageError};
use std::io::{self, Write};
use std::sync::Arc;

/// What `--summary` takes beyond the depth rules.
#[derive(Debug, Clone)]
pub struct CoverageOptions {
    pub min_baseq: u8,
    /// Skip reads whose CIGAR query length (M I S = X) is below this (samtools coverage -l).
    pub min_read_length: u64,
    /// The one region: reference name and window (0-based start; end as given, or LN when open).
    pub region: Option<(Vec<u8>, Interval)>,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CoverageSummary {
    pub rows: u64,
    /// The most counters (not bytes) any input's ring held at once.
    pub max_ring: usize,
}

/// The counts of one reference.
#[derive(Clone, Copy, Default)]
struct RefStats {
    reads: u64,
    mapq_sum: u64,
    qual_sum: u64,
    qual_bases: u64,
    covered_bases: u64,
    depth_sum: u64,
    /// A counted read with a span reached this reference.
    covered: bool,
}

struct SummarySink {
    stats: Vec<RefStats>,
    min_baseq: u8,
    /// The reference and window of `--region`; reads and positions outside it never arrive.
    window: Option<(usize, Interval)>,
    hdr: Arc<AlnHeader>,
}

/// The sum and count of the qualities of the M, = and X bases of `a` at positions inside `span`
/// and `window` that pass the quality gate (`q >= min_baseq`, or a missing quality, `0xFF`, as
/// `DepthRing` gates them). A read with SEQ `*` has none.
fn quality_in(a: &AlnRef<'_>, span: Interval, window: Interval, min_baseq: u8) -> (u64, u64) {
    if a.seq_len() == 0 {
        return (0, 0);
    }
    let quals = a.quals();
    let from = span.start.max(window.start);
    let to = span.end.min(window.end);
    let (mut sum, mut count) = (0u64, 0u64);
    let (mut rpos, mut q) = (span.start, 0usize);
    a.for_each_op(&mut |op, n| {
        let n = u64::from(n);
        match op {
            b'M' | b'=' | b'X' => {
                for p in rpos.max(from)..rpos.saturating_add(n).min(to) {
                    let c = quals.get(q + (p - rpos) as usize).copied().unwrap_or(0);
                    if c == 0xFF || c >= min_baseq {
                        sum += u64::from(c);
                        count += 1;
                    }
                }
                q += n as usize;
                rpos = rpos.saturating_add(n);
            }
            b'D' | b'N' => rpos = rpos.saturating_add(n),
            b'I' | b'S' => q += n as usize,
            _ => {}
        }
    });
    (sum, count)
}

impl DepthSink for SummarySink {
    fn position(&mut self, rid: usize, _pos: u64, depths: &[u32]) -> io::Result<()> {
        let d: u64 = depths.iter().map(|&d| u64::from(d)).sum();
        let s = &mut self.stats[rid];
        s.depth_sum += d;
        s.covered_bases += u64::from(d > 0);
        Ok(())
    }

    fn zeros(&mut self, _rid: usize, _a: u64, _b: u64) -> io::Result<()> {
        Ok(())
    }

    fn read(&mut self, _input: usize, rid: usize, a: &AlnRef<'_>, span: Option<Interval>) {
        let ref_len = self.hdr.refs[rid].1;
        let s = &mut self.stats[rid];
        s.reads += 1;
        s.mapq_sum += u64::from(a.mapq());
        if let Some(span) = span {
            s.covered = true;
            let window = self.window.map_or(
                Interval {
                    start: 0,
                    end: ref_len,
                },
                |(_, w)| w,
            );
            let (sum, n) = quality_in(a, span, window, self.min_baseq);
            s.qual_sum += sum;
            s.qual_bases += n;
        }
    }
}

fn write_row(buf: &mut Vec<u8>, name: &[u8], start: u64, end: u64, s: &RefStats) {
    let len = end.saturating_sub(start) as f64;
    // A zero divisor prints 0.
    let ratio = |num: u64, den: f64| if den > 0.0 { num as f64 / den } else { 0.0 };
    buf.extend_from_slice(name);
    for v in [start + 1, end, s.reads, s.covered_bases] {
        buf.push(b'\t');
        write_u64(buf, v);
    }
    buf.push(b'\t');
    write_g(buf, 100.0 * ratio(s.covered_bases, len), 6);
    buf.push(b'\t');
    write_g(buf, ratio(s.depth_sum, len), 6);
    buf.push(b'\t');
    write_g(buf, ratio(s.qual_sum, s.qual_bases as f64), 3);
    buf.push(b'\t');
    write_g(buf, ratio(s.mapq_sum, s.reads as f64), 3);
    buf.push(b'\n');
}

fn write_err(source: io::Error) -> IoError {
    IoError::Write {
        path: "output".into(),
        source,
    }
}

/// `samtools coverage` over `opened` (opened from `paths`), written to `out`: the header line,
/// then a row per reference, those a counted, mapped read reaches first. With a region, its one
/// row.
pub fn coverage_summary(
    opened: Vec<Box<dyn AlnSource>>,
    paths: &[String],
    o: &CoverageOptions,
    ro: &ReaderOptions,
    out: &mut dyn Write,
) -> Result<CoverageSummary, IoError> {
    let d = DepthOptions {
        rules: CountRules {
            deletions: false,
            ref_skips: false,
            min_baseq: o.min_baseq,
        },
        min_read_length: o.min_read_length,
        all: AllPositions::Covered,
        regions: o.region.clone().map(|r| vec![r]),
        skip_mate_overlap: false,
        column_names: false,
    };
    let mode = ModeRules {
        read_length: ReadLength::CigarQuery,
        filter: None,
        shared_bases: false,
    };
    let (DepthSummary { max_ring, .. }, sink) = run_sources(opened, paths, &d, mode, ro, |hdr| {
        let window = match &o.region {
            None => None,
            Some((name, iv)) => {
                let rid = hdr
                    .refs
                    .iter()
                    .position(|(n, _)| n == name)
                    .ok_or_else(|| {
                        UsageError(format!(
                            "region reference {} is not in the header",
                            String::from_utf8_lossy(name)
                        ))
                    })?;
                Some((rid, *iv))
            }
        };
        Ok(SummarySink {
            stats: vec![RefStats::default(); hdr.refs.len()],
            min_baseq: o.min_baseq,
            window,
            hdr: hdr.clone(),
        })
    })?;
    let mut buf: Vec<u8> =
        b"#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
            .to_vec();
    let mut rows = 0u64;
    let mut flush = |buf: &mut Vec<u8>, force: bool| -> Result<(), IoError> {
        if force || buf.len() >= 1 << 20 {
            out.write_all(buf).map_err(write_err)?;
            buf.clear();
        }
        Ok(())
    };
    if let Some((rid, w)) = sink.window {
        write_row(
            &mut buf,
            &sink.hdr.refs[rid].0,
            w.start,
            w.end,
            &sink.stats[rid],
        );
        rows += 1;
    } else {
        for covered in [true, false] {
            for (rid, (name, len)) in sink.hdr.refs.iter().enumerate() {
                let s = &sink.stats[rid];
                if s.covered == covered {
                    write_row(&mut buf, name, 0, *len, s);
                    rows += 1;
                    flush(&mut buf, false)?;
                }
            }
        }
    }
    flush(&mut buf, true)?;
    out.flush().map_err(write_err)?;
    Ok(CoverageSummary { rows, max_ring })
}
