//! Per-chromosome parallel relations on indexed inputs (D41).

use super::files::{RelationFiles, RelationStep, RunOptions, relation_chroms};
use super::{RelateInputs, RelateOptions, StrandWatch, has_threshold};
use gloci_core::{ChromTable, StrandRule};
use gloci_engine::{SweepSummary, run_chunks};
use gloci_io::{IndexedFile, IoError, RecordReader, RecordWriter, UsageError};
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

/// Load each input's index once (the headers come from the readers already
/// open), then run `step` on chunks of whole chromosomes in parallel,
/// writing the chunks in chromosome-table order.
pub(super) fn run_indexed<R: Send>(
    files: RelationFiles,
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
        ..
    } = files;
    let mut chroms = match chroms {
        Some(c) => c,
        None => relation_chroms(&inputs, relate)?,
    };
    let open = |p: &str, r: &dyn RecordReader, ro| {
        IndexedFile::open(Path::new(p), r.format(), r.header().clone(), ro)
    };
    let file_a = open(&paths.a, inputs.a.as_ref(), &paths.ro_a)?;
    let files_b = paths
        .b
        .iter()
        .zip(&inputs.b)
        .map(|(p, r)| open(p, r.as_ref(), &paths.ro_b))
        .collect::<Result<Vec<_>, _>>()?;
    drop(inputs);

    for f in std::iter::once(&file_a).chain(&files_b) {
        if f.info().stale {
            log::warn!(
                "the index {} is older than {name}; it may be out of date (rebuild it with: gloci index {name})",
                f.info().path.display(),
                name = f.name()
            );
        }
    }
    let names_a = file_a.reference_names();
    intern_index_names(&mut chroms, &file_a, &names_a)?;
    if chroms.is_closed() {
        for f in &files_b {
            intern_index_names(&mut chroms, f, &f.reference_names())?;
        }
    }
    let mut ordered = names_a;
    ordered.sort_by(|x, y| {
        let (x, y) = (chroms.get(x), chroms.get(y));
        chroms.cmp(x.expect("interned"), y.expect("interned"))
    });
    let chunks: Vec<Vec<Vec<u8>>> = ordered.into_iter().map(|c| vec![c]).collect();

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
        let a = watch(file_a.chromosomes_reader(chunk)?, &flags[0]);
        let b = files_b
            .iter()
            .zip(&flags[1..])
            .map(|(f, flag)| Ok(watch(f.chromosomes_reader(chunk)?, flag)))
            .collect::<Result<Vec<_>, IoError>>()?;
        step(RelateInputs { a, b }, None, sink)
    };
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
