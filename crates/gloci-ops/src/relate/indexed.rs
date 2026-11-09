//! Per-chromosome parallel relations on indexed inputs (D41) and plain
//! seekable files (D76).

use super::files::{RelationFiles, RelationStep, RunOptions, Source, relation_chroms};
use super::regions::RegionFilter;
use super::{ChunkScope, RelateInputs, RelateOptions, StrandWatch, has_threshold};
use crate::view::not_found;
use gloci_core::{ChromTable, Interval, StrandRule};
use gloci_engine::{SweepSummary, run_chunks};
use gloci_io::{IndexedFile, IoError, RecordReader, RecordWriter, SeekableFile, UsageError};
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

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

    /// The records that overlap `wanted` (all on `chunk`'s chromosomes).
    fn regions_reader(
        &self,
        chunk: &[Vec<u8>],
        wanted: &[(Vec<u8>, Interval)],
        inclusive_points: bool,
    ) -> Result<Box<dyn RecordReader>, IoError> {
        match self {
            Src::Indexed(f) => f.regions_reader(wanted, inclusive_points),
            Src::Plain(f) => Ok(Box::new(RegionFilter::new(
                f.chromosomes_reader(chunk)?,
                wanted,
                inclusive_points,
            ))),
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
    let mut opened = opened.into_iter();
    let file_a = opened.next().expect("A is always an input");
    let files_b: Vec<Src> = opened.collect();
    let all = || std::iter::once(&file_a).chain(&files_b);

    // The table with every plain input's chromosomes in it, which knows
    // them for the order check; the real one stays as it is for a fallback
    // (interning on a declared table is reported on the next read).
    let mut probe = chroms.clone();
    let mut refused = false;
    for f in all().filter(|f| matches!(f, Src::Plain(_))) {
        for name in f.reference_names() {
            refused |= probe.intern(&name).is_err();
        }
    }
    if refused || all().any(|f| out_of_table_order(&probe, f)) {
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
        if let Src::Indexed(f) = f
            && f.info().stale
        {
            log::warn!(
                "the index {} is older than {name}; it may be out of date (rebuild it with: gloci index {name})",
                f.info().path.display(),
                name = f.name()
            );
        }
    }
    let names_a = file_a.reference_names();
    if let Src::Indexed(f) = &file_a {
        intern_index_names(&mut chroms, f, &names_a)?;
    }
    if chroms.is_closed() {
        for f in &files_b {
            if let Src::Indexed(f) = f {
                intern_index_names(&mut chroms, f, &f.reference_names())?;
            }
        }
    }
    if let (Some(r), Src::Plain(_)) = (&a_regions, &file_a)
        && r.warn_unseen
    {
        // Streaming A warns when it has read the file to its end.
        let have: std::collections::HashSet<Vec<u8>> = names_a.iter().cloned().collect();
        let wanted: std::collections::BTreeSet<&Vec<u8>> =
            r.wanted.iter().map(|(c, _)| c).collect();
        for c in wanted.into_iter().filter(|c| !have.contains(*c)) {
            log::warn!("{}", not_found(c, file_a.name(), |n| have.contains(n)));
        }
    }
    let mut ordered = names_a;
    // With regions, A is empty off the chromosomes that have one, and so is
    // everything computed from it: skip them (B included).
    if let Some(r) = &a_regions {
        let with_regions: std::collections::HashSet<&[u8]> =
            r.wanted.iter().map(|(c, _)| c.as_slice()).collect();
        ordered.retain(|c| with_regions.contains(c.as_slice()));
    }
    ordered.sort_by(|x, y| {
        let (x, y) = (chroms.get(x), chroms.get(y));
        chroms.cmp(x.expect("interned"), y.expect("interned"))
    });
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
    let chunks = group_chromosomes(ordered, &weights, CHUNKS_PER_THREAD * run.threads.get());

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
    let chunk_step = |chunk: &Vec<Vec<u8>>,
                      sink: &mut gloci_engine::ChunkSink|
     -> Result<(SweepSummary, R), IoError> {
        let a = match &a_regions {
            None => file_a.chromosomes_reader(chunk)?,
            Some(r) => {
                let wanted: Vec<_> = r
                    .wanted
                    .iter()
                    .filter(|(c, _)| chunk.contains(c))
                    .cloned()
                    .collect();
                file_a.regions_reader(chunk, &wanted, r.inclusive_points)?
            }
        };
        let a = watch(a, &flags[0]);
        let b = files_b
            .iter()
            .zip(&flags[1..])
            .map(|(f, flag)| Ok(watch(f.chromosomes_reader(chunk)?, flag)))
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
            for (i, f) in files_b.iter().enumerate() {
                if summary.b_records.get(i).copied().unwrap_or(0) > 0 {
                    continue;
                }
                let b_refs = f.non_empty_reference_names();
                if !b_refs.is_empty() && !a_refs.iter().any(|a| b_refs.contains(a)) {
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
fn group_chromosomes(chroms: Vec<Vec<u8>>, weights: &[u64], target: usize) -> Vec<Vec<Vec<u8>>> {
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

/// Intern the chromosomes `f`'s index lists; one missing from a closed
/// table is an error naming the index file and the chromosome.
fn intern_index_names(
    chroms: &mut ChromTable,
    f: &IndexedFile,
    names: &[Vec<u8>],
) -> Result<(), IoError> {
    for name in names {
        chroms.intern(name).map_err(|e| IoError::Detect {
            path: f.name().to_string(),
            msg: format!(
                "{} (index {}): {}",
                f.name(),
                f.info().path.display(),
                gloci_io::naming::unknown_chrom_err(chroms.order_source(), &e.0)
            ),
        })?;
    }
    Ok(())
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
