//! `gloci depth`: the depth of coverage of alignments at each position, as `samtools depth`
//! (D137, D138). The counting rules live in [`DepthRing`], one per input; the rows in
//! `rows.rs`. The inputs are one merged stream at every thread count; threads decode.

mod mates;
mod ring;
mod rows;

pub use mates::MateTable;
pub use ring::{CountRules, DepthRing};

use crate::aln_merge::threads_per_reader;
use gloci_core::{Interval, Region};
use gloci_io::{
    AlnHeader, AlnSource, IoError, ReaderOptions, UsageError, open_aln_source, regions_source_from,
};
pub(crate) use rows::log_diagnostics;
use rows::{Emitter, Input, Plan, run_group};
use std::collections::HashSet;
use std::io::{self, Write};
use std::num::NonZeroUsize;
use std::sync::Arc;
use std::sync::atomic::AtomicBool;

/// Which positions print besides those a read covers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AllPositions {
    /// Only positions inside a passing read's span in some input.
    Covered,
    /// `-a`: every position of each reference that has a passing read.
    RefsWithReads,
    /// `-aa`: every position of every reference.
    AllRefs,
}

/// Regions as `(reference, interval)` pairs.
pub type WantedRegions = Vec<(Vec<u8>, Interval)>;

#[derive(Debug, Clone)]
pub struct DepthOptions {
    pub rules: CountRules,
    /// Skip reads shorter than this (`--min-read-length`).
    pub min_read_length: u64,
    pub all: AllPositions,
    /// Print only positions in the union of these regions (0-based, half-open, by reference
    /// name as the header lists it).
    pub regions: Option<Vec<(Vec<u8>, Interval)>>,
    /// Count each base of an overlapping mate pair once.
    pub skip_mate_overlap: bool,
    /// Write `#CHROM POS <path>...` first.
    pub column_names: bool,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct DepthSummary {
    pub rows: u64,
    /// The most counters (not bytes) any input's ring held at once.
    pub max_ring: usize,
    /// The most reads any input's mate table held at once (`--skip-mate-overlap`).
    pub max_pending_mates: usize,
}

fn write_err(source: io::Error) -> IoError {
    IoError::Write {
        path: "output".into(),
        source,
    }
}

/// The regions of `--region` and `--regions-file` as `(reference, interval)` pairs, checked
/// against `header`, the header of the first input (named `first`), which every input shares:
/// the input is opened once, so a pipe or stdin works. A region on a reference the file lacks
/// is dropped with a warning. `None` when none is given.
pub fn depth_regions(
    header: &AlnHeader,
    first: &str,
    regions: &[String],
    regions_file: Option<&str>,
    ro: &ReaderOptions,
) -> Result<Option<WantedRegions>, IoError> {
    if regions.is_empty() && regions_file.is_none() {
        return Ok(None);
    }
    let known: HashSet<&[u8]> = header.refs.iter().map(|(n, _)| &n[..]).collect();
    let mut wanted: Vec<(Vec<u8>, Interval)> = Vec::new();
    for s in regions {
        let r = Region::parse(s, |n| known.contains(n)).map_err(|e| UsageError(e.to_string()))?;
        let iv = r.interval(None);
        wanted.push((r.chrom, iv));
    }
    if let Some(f) = regions_file {
        wanted.extend(crate::view::read_regions_file(f, ro)?);
    }
    let mut warned = HashSet::new();
    wanted.retain(|(c, _)| {
        if known.contains(&c[..]) {
            return true;
        }
        if warned.insert(c.clone()) {
            log::warn!(
                "region on \"{}\" skipped: {first} has no such reference",
                String::from_utf8_lossy(c)
            );
        }
        false
    });
    Ok(Some(wanted))
}

/// Refuse inputs whose `@SQ` lists differ, naming the first difference (D126).
pub(crate) fn check_same_references(
    names: &[String],
    headers: &[Arc<AlnHeader>],
) -> Result<(), IoError> {
    let first = &headers[0];
    for (k, h) in headers.iter().enumerate().skip(1) {
        let n = first.refs.len().max(h.refs.len());
        for i in 0..n {
            let (a, b) = (first.refs.get(i), h.refs.get(i));
            if a == b {
                continue;
            }
            let show = |r: Option<&(Vec<u8>, u64)>| match r {
                Some((name, len)) => format!("{} (LN {len})", String::from_utf8_lossy(name)),
                None => "missing".to_string(),
            };
            return Err(IoError::Detect {
                path: names[k].clone(),
                msg: format!(
                    "{} and {} have different @SQ lists: @SQ {} is {} in the first and {} in the \
second; every input needs the same references in the same order (rebuild one with the other's \
header, for example gloci reheader)",
                    names[0],
                    names[k],
                    i + 1,
                    show(a),
                    show(b)
                ),
            });
        }
    }
    Ok(())
}

/// The inputs of `depth`, each opened once (`-` is stdin), sharing `threads` to decode them:
/// each gets `threads / inputs`, at least 1, so the threads do not grow with the inputs.
pub fn open_depth_inputs(
    paths: &[String],
    ro: &ReaderOptions,
    threads: NonZeroUsize,
) -> Result<Vec<Box<dyn AlnSource>>, IoError> {
    if paths.is_empty() {
        return Err(UsageError("depth needs at least one input file".into()).into());
    }
    let mut ro = ro.clone();
    ro.threads = threads_per_reader(threads, paths.len());
    paths.iter().map(|p| open_aln_source(p, &ro)).collect()
}

/// `samtools depth` over `paths`, written to `out`: [`depth_sources`] of
/// [`open_depth_inputs`].
pub fn depth(
    paths: &[String],
    o: &DepthOptions,
    ro: &ReaderOptions,
    threads: NonZeroUsize,
    out: &mut dyn Write,
) -> Result<DepthSummary, IoError> {
    let opened = open_depth_inputs(paths, ro, threads)?;
    let mut ro = ro.clone();
    ro.threads = threads;
    depth_sources(opened, paths, o, &ro, out)
}

/// `samtools depth` over `opened`, the inputs [`open_depth_inputs`] opened from `paths`, written
/// to `out`. Rows come in `@SQ` order. The inputs are read as one merged stream at every thread
/// count (D138): threads go to BGZF and CRAM decoding, so the output and the diagnostics never
/// depend on them. With regions, a regular file with an index is read through it; anything
/// else is scanned from where it stands.
pub fn depth_sources(
    mut opened: Vec<Box<dyn AlnSource>>,
    paths: &[String],
    o: &DepthOptions,
    ro: &ReaderOptions,
    out: &mut dyn Write,
) -> Result<DepthSummary, IoError> {
    if opened.is_empty() || opened.len() != paths.len() {
        return Err(UsageError("depth needs at least one input file".into()).into());
    }
    let headers: Vec<Arc<AlnHeader>> = opened.iter().map(|s| s.header().clone()).collect();
    let names: Vec<String> = opened.iter().map(|s| s.name().to_string()).collect();
    check_same_references(&names, &headers)?;
    let hdr = headers[0].clone();
    let warned_past = AtomicBool::new(false);
    let plan = Plan::new(o, &hdr.refs, &warned_past);
    // One stream. Regions come through the indexes where there are some, else by scanning;
    // a reopened input gets its share of `ro.threads`, as `open_depth_inputs` gives.
    if let Some(regions) = &o.regions {
        let mut ro = ro.clone();
        ro.threads = threads_per_reader(ro.threads, paths.len());
        opened = opened
            .into_iter()
            .zip(paths)
            .map(|(s, p)| regions_source_from(s, p, regions, &ro))
            .collect::<Result<_, _>>()?;
    }

    if o.column_names {
        let mut line = b"#CHROM\tPOS".to_vec();
        for p in paths {
            line.push(b'\t');
            line.extend_from_slice(p.as_bytes());
        }
        line.push(b'\n');
        out.write_all(&line).map_err(write_err)?;
    }

    let mut inputs: Vec<Input> = opened
        .into_iter()
        .zip(&names)
        .map(|(s, n)| Input::new(Some(s), n, hdr.clone()))
        .collect();
    let mut summary = DepthSummary::default();
    let mut write = |b: &[u8]| out.write_all(b);
    let mut em = Emitter::new(&mut write);
    (summary.max_ring, summary.max_pending_mates) =
        run_group(&mut inputs, 0, hdr.refs.len(), &plan, &mut em)?;
    summary.rows = em.rows;
    for i in &inputs {
        log_diagnostics(&i.name, &i.diagnostics());
    }
    out.flush().map_err(write_err)?;
    Ok(summary)
}
