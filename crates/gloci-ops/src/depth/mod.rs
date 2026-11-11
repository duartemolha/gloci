//! `gloci depth`: the depth of coverage of alignments at each position, as `samtools depth`
//! (D137, D138). The counting rules live in [`DepthRing`], one per input; the rows in
//! `rows.rs`. With indexes and several threads each run of references is its own chunk.

mod ring;
mod rows;

pub use ring::{CountRules, DepthRing};

use gloci_core::{Interval, Region};
use gloci_engine::{ChunkSink, run_chunks};
use gloci_io::{
    AlnHeader, AlnSource, Format, Header, IndexedFile, IoError, ReaderOptions, RecordWriter,
    UsageError, find_index, open_aln_source, open_aln_source_regions,
};
use rows::{Emitter, Input, Plan, run_group};
use std::collections::{HashMap, HashSet};
use std::io::{self, Write};
use std::num::NonZeroUsize;
use std::path::Path;
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

#[derive(Debug, Clone)]
pub struct DepthOptions {
    pub rules: CountRules,
    /// Skip reads shorter than this (`--min-read-length`).
    pub min_read_length: u64,
    pub all: AllPositions,
    /// Print only positions in the union of these regions (0-based, half-open, by reference
    /// name as the header lists it).
    pub regions: Option<Vec<(Vec<u8>, Interval)>>,
    /// Count each base of an overlapping mate pair once (Task 7).
    pub skip_mate_overlap: bool,
    /// Write `#CHROM POS <path>...` first.
    pub column_names: bool,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct DepthSummary {
    pub rows: u64,
    /// The most counters (not bytes) any input's ring held at once.
    pub max_ring: usize,
    pub max_pending_mates: usize,
}

/// Writes each record as a line, for the chunks drained in order.
struct Lines<'a>(&'a mut dyn Write);

impl RecordWriter for Lines<'_> {
    fn write_record(&mut self, raw: &[u8]) -> io::Result<()> {
        self.0.write_all(raw)?;
        self.0.write_all(b"\n")
    }
    fn finish(self: Box<Self>) -> Result<(), IoError> {
        Ok(())
    }
}

fn write_err(source: io::Error) -> IoError {
    IoError::Write {
        path: "output".into(),
        source,
    }
}

/// The regions of `--region` and `--regions-file` as `(reference, interval)` pairs, checked
/// against the references of `first` (the first input, whose header every input shares). A
/// region on a reference the file lacks is dropped with a warning. `None` when none is given.
pub fn depth_regions(
    first: &str,
    regions: &[String],
    regions_file: Option<&str>,
    ro: &ReaderOptions,
) -> Result<Option<Vec<(Vec<u8>, Interval)>>, IoError> {
    if regions.is_empty() && regions_file.is_none() {
        return Ok(None);
    }
    if first == "-" {
        return Err(UsageError(
            "--region and --regions-file need a file to read the references from, not standard \
input; write the input to a file first"
                .into(),
        )
        .into());
    }
    let header = open_aln_source(first, ro)?.header().clone();
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
fn check_same_references(names: &[String], headers: &[Arc<AlnHeader>]) -> Result<(), IoError> {
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

/// [`depth_in`] with spill files, if any, in the system temp directory.
pub fn depth(
    paths: &[String],
    o: &DepthOptions,
    ro: &ReaderOptions,
    threads: NonZeroUsize,
    out: &mut dyn Write,
) -> Result<DepthSummary, IoError> {
    depth_in(paths, o, ro, threads, &std::env::temp_dir(), out)
}

/// One chunk of a per-chromosome run: the references `first..end` and the regions to read.
#[derive(Clone)]
struct Job {
    first: usize,
    end: usize,
    regions: Vec<(Vec<u8>, Interval)>,
}

/// `samtools depth` over `paths`, written to `out`. Rows come in `@SQ` order. With every input
/// indexed and `threads` above 1, runs of references are computed in parallel and written in
/// order, with the same output as one stream; their spill files go to `tmp_dir`.
pub fn depth_in(
    paths: &[String],
    o: &DepthOptions,
    ro: &ReaderOptions,
    threads: NonZeroUsize,
    tmp_dir: &Path,
    out: &mut dyn Write,
) -> Result<DepthSummary, IoError> {
    if paths.is_empty() {
        return Err(UsageError("depth needs at least one input file".into()).into());
    }
    let mut opened: Vec<Box<dyn AlnSource>> = paths
        .iter()
        .map(|p| open_aln_source(p, ro))
        .collect::<Result<_, _>>()?;
    let headers: Vec<Arc<AlnHeader>> = opened.iter().map(|s| s.header().clone()).collect();
    let names: Vec<String> = opened.iter().map(|s| s.name().to_string()).collect();
    check_same_references(&names, &headers)?;
    let hdr = headers[0].clone();
    let warned_past = AtomicBool::new(false);
    let plan = Plan::new(o, &hdr.refs, &warned_past);

    if o.column_names {
        let mut line = b"#CHROM\tPOS".to_vec();
        for p in paths {
            line.push(b'\t');
            line.extend_from_slice(p.as_bytes());
        }
        line.push(b'\n');
        out.write_all(&line).map_err(write_err)?;
    }

    // The files that can be read per chromosome: every one has an index.
    let indexed = threads.get() > 1
        && hdr.refs.len() > 1
        && opened.iter().zip(paths).all(|(s, p)| {
            matches!(s.format(), Format::Bam | Format::Cram)
                && p != "-"
                && Path::new(p).is_file()
                && find_index(Path::new(p)).is_some()
        });
    if indexed {
        let files: Vec<IndexedFile> = opened
            .iter()
            .zip(paths)
            .map(|(s, p)| {
                IndexedFile::open(Path::new(p), s.format(), Header::Alignment(hdr.clone()), ro)
            })
            .collect::<Result<_, _>>()?;
        drop(opened);
        return chunked(paths, &names, &files, &hdr, o, &plan, threads, tmp_dir, out);
    }

    // One stream. Regions come through the indexes where there are some, else by scanning.
    if let Some(regions) = &o.regions {
        for (k, p) in paths.iter().enumerate() {
            opened[k] = open_aln_source_regions(p, regions, ro)?;
        }
    }
    let mut inputs: Vec<Input> = opened
        .into_iter()
        .zip(&names)
        .map(|(s, n)| Input::new(Some(s), n, hdr.clone()))
        .collect();
    let mut summary = DepthSummary::default();
    let mut write = |b: &[u8]| out.write_all(b);
    let mut em = Emitter::new(&mut write);
    summary.max_ring = run_group(&mut inputs, 0, hdr.refs.len(), &plan, &mut em)?;
    summary.rows = em.rows;
    for i in inputs {
        i.log_diagnostics();
    }
    out.flush().map_err(write_err)?;
    Ok(summary)
}

#[allow(clippy::too_many_arguments)]
fn chunked(
    paths: &[String],
    names: &[String],
    files: &[IndexedFile],
    hdr: &Arc<AlnHeader>,
    o: &DepthOptions,
    plan: &Plan<'_>,
    threads: NonZeroUsize,
    tmp_dir: &Path,
    out: &mut dyn Write,
) -> Result<DepthSummary, IoError> {
    let _ = paths;
    let refs = &hdr.refs;
    let weights: Vec<u64> = refs
        .iter()
        .map(|(n, _)| files.iter().map(|f| f.weight(n)).sum::<u64>() + 1)
        .collect();
    let mut by_name: HashMap<&[u8], Vec<Interval>> = HashMap::new();
    if let Some(list) = &o.regions {
        for (n, iv) in list {
            by_name.entry(&n[..]).or_default().push(*iv);
        }
    }
    let mut first = 0;
    let jobs: Vec<Job> = crate::relate::group_chromosomes(
        refs.iter().map(|(n, _)| n.clone()).collect(),
        &weights,
        4 * threads.get(),
    )
    .into_iter()
    .map(|g| {
        let end = first + g.len();
        let mut regions = Vec::new();
        for r in first..end {
            let name = &refs[r].0;
            if o.regions.is_some() {
                for iv in by_name.get(&name[..]).into_iter().flatten() {
                    regions.push((name.clone(), *iv));
                }
            } else if weights[r] > 1 {
                regions.push((
                    name.clone(),
                    Interval {
                        start: 0,
                        end: refs[r].1,
                    },
                ));
            }
        }
        let job = Job {
            first,
            end,
            regions,
        };
        first = end;
        job
    })
    .collect();
    log::debug!("depth running per chromosome: {} chunk(s)", jobs.len());
    let step = |job: &Job, sink: &mut ChunkSink| -> Result<DepthSummary, IoError> {
        let mut inputs = Vec::with_capacity(files.len());
        for (f, n) in files.iter().zip(names) {
            let src = if job.regions.is_empty() {
                None
            } else {
                Some(f.aln_regions_source(&job.regions, false)?)
            };
            inputs.push(Input::new(src, n, hdr.clone()));
        }
        // A chunk's rows are one record per batch, without the final newline `Lines` adds.
        let mut write = |b: &[u8]| sink.write_record(&b[..b.len() - 1]);
        let mut em = Emitter::new(&mut write);
        let max_ring = run_group(&mut inputs, job.first, job.end, plan, &mut em)?;
        let rows = em.rows;
        for i in inputs {
            i.log_diagnostics();
        }
        Ok(DepthSummary {
            rows,
            max_ring,
            max_pending_mates: 0,
        })
    };
    let mut lines = Lines(out);
    let done = run_chunks(jobs, threads, tmp_dir, &step, &mut lines)?;
    let mut total = DepthSummary::default();
    for s in done {
        total.rows += s.rows;
        total.max_ring = total.max_ring.max(s.max_ring);
        total.max_pending_mates = total.max_pending_mates.max(s.max_pending_mates);
    }
    out.flush().map_err(write_err)?;
    Ok(total)
}
