//! Per-chromosome parallel relations on indexed inputs (D41) and plain
//! seekable files (D76).

use super::files::{RelationFiles, RelationStep, RunOptions, Source, relation_chroms};
use super::{ChunkScope, RelateInputs, RelateOptions, StrandWatch, has_threshold};
use gloci_core::{ChromTable, StrandRule};
use gloci_engine::{SweepSummary, run_chunks};
use gloci_io::{IndexedFile, IoError, RecordReader, RecordWriter, SeekableFile, UsageError};
use std::collections::HashSet;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

/// A chunk of a per-chromosome run: A's chromosomes, and for each B the
/// chromosomes its reader is given.
#[derive(Clone)]
struct Job {
    a: Vec<Vec<u8>>,
    b: Vec<Vec<Vec<u8>>>,
}

/// One input, ready to be read a chromosome at a time.
enum Src {
    Indexed(IndexedFile),
    Plain(SeekableFile),
}

impl Src {
    fn name(&self) -> &str {
        match self {
            Src::Indexed(f) => f.name(),
            Src::Plain(f) => f.name(),
        }
    }

    fn reference_names(&self) -> Vec<Vec<u8>> {
        match self {
            Src::Indexed(f) => f.reference_names(),
            Src::Plain(f) => f.reference_names(),
        }
    }

    /// The chromosomes in the order the file's records show them.
    fn file_order_reference_names(&self) -> Vec<Vec<u8>> {
        match self {
            Src::Indexed(f) => f.file_order_reference_names(),
            Src::Plain(f) => f.reference_names(),
        }
    }

    fn non_empty_reference_names(&self) -> Vec<Vec<u8>> {
        match self {
            Src::Indexed(f) => f.non_empty_reference_names(),
            Src::Plain(f) => f.reference_names(),
        }
    }

    fn weight(&self, chrom: &[u8]) -> u64 {
        match self {
            Src::Indexed(f) => f.weight(chrom),
            Src::Plain(f) => f.weight(chrom),
        }
    }

    fn chromosomes_reader(&self, chroms: &[Vec<u8>]) -> Result<Box<dyn RecordReader>, IoError> {
        match self {
            Src::Indexed(f) => f.chromosomes_reader(chroms),
            Src::Plain(f) => f.chromosomes_reader(chroms),
        }
    }
}

/// Load each input's index once, or scan each plain input for where its
/// chromosomes lie (the headers come from the readers already open), then
/// run `step` on chunks of whole chromosomes in parallel, writing the
/// chunks in chromosome-table order. `sources` says how each input (A, then
/// each B) is read.
///
/// What streaming would refuse, streaming refuses: an input whose
/// chromosomes are not contiguous, are in an order the table rejects, or
/// are unknown to a closed table is read whole by `step` (as an indexed
/// input out of table order is), which raises the error with its line
/// numbers exactly as at `-t 1`.
pub(super) fn run_per_chromosome<R: Send>(
    files: RelationFiles,
    sources: &[Source],
    relate: &RelateOptions,
    run: &RunOptions,
    out: &mut dyn RecordWriter,
    step: &RelationStep<R>,
) -> Result<(SweepSummary, Vec<R>), IoError> {
    if relate.window.is_some() && has_threshold(&relate.criteria) {
        return Err(UsageError(
            "-w/--left/--right cannot be combined with -f, -F, -r, -e or --min-bases (D37)".into(),
        )
        .into());
    }
    let RelationFiles {
        inputs,
        paths,
        chroms,
        a_regions,
        ..
    } = files;
    let mut chroms = match chroms {
        Some(c) => c,
        None => relation_chroms(&inputs, relate)?,
    };
    // Streaming scans all of a plain A for its regions, and reads B only as
    // far as the A records it selects take it: a run per chromosome would not
    // read (or validate) the same lines, so it streams.
    if a_regions.is_some() && sources.contains(&Source::Plain) {
        log::debug!("--region with a plain input; streaming");
        let (summary, r) = step(inputs, Some(chroms), out)?;
        return Ok((summary, vec![r]));
    }
    let open =
        |p: &str, r: &dyn RecordReader, ro, source: Source| -> Result<Option<Src>, IoError> {
            Ok(match source {
                Source::Indexed => Some(Src::Indexed(IndexedFile::open(
                    Path::new(p),
                    r.format(),
                    r.header().clone(),
                    ro,
                )?)),
                Source::Plain => {
                    SeekableFile::open(Path::new(p), r.format(), r.header().clone(), ro)?
                        .map(Src::Plain)
                }
            })
        };
    let mut opened = Vec::with_capacity(sources.len());
    let readers = std::iter::once(&inputs.a).chain(&inputs.b);
    let paths_in =
        std::iter::once((&paths.a, &paths.ro_a)).chain(paths.b.iter().map(|p| (p, &paths.ro_b)));
    for ((r, (p, ro)), &source) in readers.zip(paths_in).zip(sources) {
        match open(p, r.as_ref(), ro, source)? {
            Some(f) => opened.push(f),
            None => {
                log::debug!("{p} cannot be read per chromosome; streaming");
                drop(opened);
                let (summary, r) = step(inputs, Some(chroms), out)?;
                return Ok((summary, vec![r]));
            }
        }
    }
    let is_alignment: Vec<bool> = std::iter::once(&inputs.a)
        .chain(&inputs.b)
        .map(|r| r.format().is_alignment())
        .collect();
    // An alignment reader's one-time warnings are printed once, by this run.
    for f in &mut opened {
        if let Src::Indexed(f) = f {
            f.defer_warnings();
        }
    }
    let mut opened = opened.into_iter();
    let file_a = opened.next().expect("A is always an input");
    let files_b: Vec<Src> = opened.collect();
    let all = || std::iter::once(&file_a).chain(&files_b);

    // The table with every plain input's chromosomes in it, which knows
    // them for the order check; the real one stays as it is for a fallback
    // (interning on a declared table is reported on the next read).
    let mut probe = chroms.clone();
    let mut refused = false;
    // (An indexed B's names are interned only against a closed table, where
    // an unknown one is refused; an open table does not need them.)
    for (n, f) in all().enumerate() {
        if matches!(f, Src::Plain(_)) || n == 0 || probe.is_closed() {
            for name in f.reference_names() {
                refused |= probe.intern(&name).is_err();
            }
        }
    }
    // With no A records no chunk would read B, which streaming reads and
    // validates to its end.
    let a_empty = file_a
        .reference_names()
        .iter()
        .all(|c| file_a.weight(c) == 0);
    if refused
        || (a_empty && files_b.iter().any(|f| matches!(f, Src::Plain(_))))
        || all().any(|f| out_of_table_order(&probe, f))
    {
        log::debug!(
            "an input is out of table order or has a chromosome the table lacks; streaming"
        );
        drop((file_a, files_b));
        let (summary, r) = step(inputs, Some(chroms), out)?;
        return Ok((summary, vec![r]));
    }
    chroms = probe;
    drop(inputs);

    for f in all() {
        if let Src::Indexed(f) = f {
            f.info().warn_if_stale(f.name());
        }
    }
    let names_a = file_a.reference_names();
    let mut ordered = names_a;
    // With regions, A is empty off the chromosomes that have one, and so is
    // everything computed from it: skip them (B included).
    if let Some(r) = &a_regions {
        let with_regions: HashSet<&[u8]> = r.wanted.iter().map(|(c, _)| c.as_slice()).collect();
        ordered.retain(|c| with_regions.contains(c.as_slice()));
    }
    ordered.sort_by(|x, y| {
        let (x, y) = (chroms.get(x), chroms.get(y));
        chroms.cmp(x.expect("interned"), y.expect("interned"))
    });
    // A chromosome with no A records is never entered by a sweep.
    ordered.retain(|c| file_a.weight(c) > 0);
    let by_rank = |x: &Vec<u8>, y: &Vec<u8>| {
        let (x, y) = (chroms.get(x), chroms.get(y));
        chroms.cmp(x.expect("interned"), y.expect("interned"))
    };
    // Runs of consecutive chromosomes, balanced by the index's estimate of
    // their size in every input.
    let weights: Vec<u64> = ordered
        .iter()
        .map(|c| {
            std::iter::once(&file_a)
                .chain(&files_b)
                .map(|f| f.weight(c))
                .sum::<u64>()
                + 1
        })
        .collect();
    // A plain input is read to its end whatever A does (streaming drains it
    // after A ends), so every chromosome of a plain B is read by exactly one
    // chunk. Those A lacks go to the chunk holding the next A chromosome,
    // whose sweep skips past them as A moves on, and those after A's last
    // chromosome to the last chunk, which drains them.
    let a_chroms: HashSet<&Vec<u8>> = ordered.iter().collect();
    let groups = group_chromosomes(
        ordered.clone(),
        &weights,
        CHUNKS_PER_THREAD * run.threads.get(),
    );
    let mut chunks: Vec<Job> = groups
        .into_iter()
        .map(|g| Job {
            b: files_b.iter().map(|_| g.clone()).collect(),
            a: g,
        })
        .collect();
    // (No chunks, with all inputs indexed, leaves nothing to add to.)
    let last_a = chunks.last().and_then(|c| c.a.last()).cloned();
    // Past A's last chromosome: chunks of their own with no A, which read
    // (and validate) those chromosomes, spread over the workers by weight.
    let mut trailing: Vec<Vec<u8>> = Vec::new();
    for (f, src) in files_b.iter().enumerate() {
        let (Src::Plain(_), Some(last_a)) = (src, &last_a) else {
            continue;
        };
        let mut extras: Vec<Vec<u8>> = src
            .reference_names()
            .into_iter()
            .filter(|c| !a_chroms.contains(c))
            .collect();
        extras.sort_by(by_rank);
        // One walk over the extras and the chunks, both in table order.
        let mut at = 0;
        for x in extras {
            if by_rank(&x, last_a).is_gt() {
                trailing.push(x);
                continue;
            }
            while at + 1 < chunks.len()
                && by_rank(&x, chunks[at].a.last().expect("chunks hold chromosomes")).is_gt()
            {
                at += 1;
            }
            chunks[at].b[f].push(x);
        }
        for c in &mut chunks {
            c.b[f].sort_by(by_rank);
        }
    }
    trailing.sort_by(by_rank);
    trailing.dedup();
    if !trailing.is_empty() {
        let weights: Vec<u64> = trailing
            .iter()
            .map(|c| files_b.iter().map(|f| f.weight(c)).sum::<u64>() + 1)
            .collect();
        let groups = group_chromosomes(trailing, &weights, CHUNKS_PER_THREAD * run.threads.get());
        for g in groups {
            chunks.push(Job {
                a: Vec::new(),
                b: files_b
                    .iter()
                    .map(|f| match f {
                        Src::Plain(_) => g.clone(),
                        Src::Indexed(_) => Vec::new(),
                    })
                    .collect(),
            });
        }
    }

    let watching = relate.criteria.strand != StrandRule::Any;
    let flags: Vec<Arc<AtomicBool>> = (0..=files_b.len())
        .map(|_| Arc::new(AtomicBool::new(false)))
        .collect();
    let watch = |r: Box<dyn RecordReader>, flag: &Arc<AtomicBool>| -> Box<dyn RecordReader> {
        if watching {
            Box::new(StrandWatch {
                inner: r,
                unknown: flag.clone(),
            })
        } else {
            r
        }
    };
    let chunk_step =
        |job: &Job, sink: &mut gloci_engine::ChunkSink| -> Result<(SweepSummary, R), IoError> {
            let chunk = &job.a;
            let a = match (&a_regions, &file_a) {
                (Some(r), Src::Indexed(f)) => {
                    let in_chunk: HashSet<&[u8]> = chunk.iter().map(|c| c.as_slice()).collect();
                    let wanted: Vec<_> = r
                        .wanted
                        .iter()
                        .filter(|(c, _)| in_chunk.contains(c.as_slice()))
                        .cloned()
                        .collect();
                    f.regions_reader(&wanted, r.inclusive_points)?
                }
                _ => file_a.chromosomes_reader(chunk)?,
            };
            let a = watch(a, &flags[0]);
            let b = files_b
                .iter()
                .zip(&flags[1..])
                .zip(&job.b)
                .map(|((f, flag), list)| Ok(watch(f.chromosomes_reader(list)?, flag)))
                .collect::<Result<Vec<_>, IoError>>()?;
            let _scope = ChunkScope::enter();
            step(RelateInputs { a, b }, None, sink)
        };
    log::debug!(
        "running per chromosome: {} chunk(s) of whole chromosomes, {} plain input(s) scanned",
        chunks.len(),
        all().filter(|f| matches!(f, Src::Plain(_))).count()
    );
    let results = run_chunks(chunks, run.threads, &run.tmp_dir, &chunk_step, out)?;

    let mut summary = SweepSummary::default();
    let mut rs = Vec::with_capacity(results.len());
    for (s, r) in results {
        summary.merge(&s);
        rs.push(r);
    }
    // What streaming reports once for a whole alignment input, from the chunks together:
    // the counts of what it left out (D101, D102) and its one-time warnings. The reads with
    // no place at the end of A are in no chromosome's chunk: the index counts them.
    if let (Src::Indexed(fa), true) = (&file_a, is_alignment[0])
        && a_regions.is_none()
        && let Some((_, d)) = summary.skipped.first_mut()
    {
        d.unmapped += fa.unplaced_reads();
    }
    for (i, ((name, d), _)) in summary
        .skipped
        .iter()
        .zip(&is_alignment)
        .enumerate()
        .filter(|(_, (_, aln))| **aln)
    {
        d.log_deferred_warnings();
        // An indexed -b is read only where A needs it, and with regions only the regions
        // are: those counts are of the reads read, not of the whole file.
        if i > 0 || a_regions.is_some() {
            d.log_alignment_counts_among(
                name,
                " (among the reads read: an indexed -b, or --region, reads only where needed)",
            );
        } else {
            d.log_alignment_counts(name);
        }
    }
    if watching {
        for (f, flag) in std::iter::once(&file_a).chain(&files_b).zip(&flags) {
            if flag.load(Ordering::Relaxed) {
                log::warn!(
                    "{} has records with unknown strand; they never match under -s or -S",
                    f.name()
                );
            }
        }
    }
    if summary.a_records > 0 {
        let a_refs = file_a.non_empty_reference_names();
        if !a_refs.is_empty() {
            let a_set: HashSet<&Vec<u8>> = a_refs.iter().collect();
            for f in &files_b {
                let b_refs = f.non_empty_reference_names();
                if !b_refs.is_empty() && !b_refs.iter().any(|b| a_set.contains(b)) {
                    log::warn!(
                        "{}",
                        gloci_io::naming::format_disjoint_chroms_warning(
                            file_a.name(),
                            &a_refs,
                            f.name(),
                            &b_refs
                        )
                    );
                }
            }
        }
    }
    Ok((summary, rs))
}

/// Chunks per thread: enough to balance uneven chromosomes across threads,
/// few enough that per-chunk costs stay small on genomes with many contigs.
const CHUNKS_PER_THREAD: usize = 4;

/// Split `chroms` (in table order) into at most `target` runs of consecutive
/// chromosomes of about equal total weight: a run ends once the weight so
/// far reaches the next multiple of `total / target`.
pub(crate) fn group_chromosomes(
    chroms: Vec<Vec<u8>>,
    weights: &[u64],
    target: usize,
) -> Vec<Vec<Vec<u8>>> {
    let target = target.max(1) as u128;
    let total: u128 = weights.iter().map(|&w| u128::from(w)).sum::<u128>().max(1);
    let mut groups = Vec::new();
    let mut current = Vec::new();
    let (mut acc, mut k) = (0u128, 1u128);
    for (c, &w) in chroms.into_iter().zip(weights) {
        current.push(c);
        acc += u128::from(w);
        if acc * target >= k * total {
            groups.push(std::mem::take(&mut current));
            while k <= target && acc * target >= k * total {
                k += 1;
            }
        }
    }
    if !current.is_empty() {
        groups.push(current);
    }
    groups
}

/// Whether the chromosomes of `f`, in the order its records show them, go
/// backwards in `chroms`' order. Names the table cannot order are skipped
/// (an unknown name is reported when the table is closed over them).
fn out_of_table_order(chroms: &ChromTable, f: &Src) -> bool {
    // Names the table cannot order are dropped first, so a backwards step
    // across one of them is still seen.
    let known: Vec<Vec<u8>> = f
        .file_order_reference_names()
        .into_iter()
        .filter(|n| chroms.cmp_names(n, n).is_some())
        .collect();
    known
        .windows(2)
        .any(|w| chroms.cmp_names(&w[0], &w[1]) == Some(std::cmp::Ordering::Greater))
}

#[cfg(test)]
mod tests {
    use super::group_chromosomes;

    fn names(n: usize) -> Vec<Vec<u8>> {
        (0..n).map(|i| format!("c{i}").into_bytes()).collect()
    }

    #[test]
    fn groups_are_bounded_ordered_and_complete() {
        for (n, target) in [(5000, 16), (3, 16), (16, 16), (100, 1), (1, 4)] {
            let chroms = names(n);
            let weights: Vec<u64> = (0..n as u64).map(|i| 1 + (i * 7919) % 1000).collect();
            let groups = group_chromosomes(chroms.clone(), &weights, target);
            assert!(
                groups.len() <= target.max(1),
                "{n} {target}: {}",
                groups.len()
            );
            assert!(groups.iter().all(|g| !g.is_empty()));
            assert_eq!(groups.concat(), chroms);
        }
        // A heavy chromosome ends its run; the light ones after it share one.
        let g = group_chromosomes(names(4), &[1, 1000, 1, 1], 4);
        assert_eq!(g, vec![names(2), vec![b"c2".to_vec(), b"c3".to_vec()]]);
    }
}
