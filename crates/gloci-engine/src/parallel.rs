//! Per-chromosome parallel driver for indexed inputs (D41).

use crate::chunk::ChunkSink;
use gloci_io::{Compression, IoError, RecordWriter, detect_compression, find_index};
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::Read;
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

/// Returns true when every path in `paths` is a regular BGZF file with a `.tbi`
/// or `.csi` index, and no path is `"-"`. Returns false if `paths` is empty.
pub fn all_indexed(paths: &[&str]) -> bool {
    if paths.is_empty() {
        return false;
    }
    for &path in paths {
        if path == "-" {
            return false;
        }
        let p = Path::new(path);
        let Ok(meta) = fs::metadata(p) else {
            return false;
        };
        if !meta.is_file() {
            return false;
        }
        let Ok(mut file) = File::open(p) else {
            return false;
        };
        let mut head = [0u8; 18];
        let Ok(n) = file.read(&mut head) else {
            return false;
        };
        if detect_compression(&head[..n]) != Compression::Bgzf {
            return false;
        }
        if find_index(p).is_none() {
            return false;
        }
    }
    true
}

/// Run `work` for each chromosome in `chroms` across up to `threads` workers,
/// with at most `threads` chunks in flight. Drains finished chunks to `out` in
/// strict chromosome order. Returns results in chromosome order.
#[allow(clippy::type_complexity)]
pub fn run_per_chrom<R: Send>(
    chroms: Vec<Vec<u8>>,
    threads: NonZeroUsize,
    tmp_dir: &Path,
    work: &(dyn Fn(&[u8], &mut ChunkSink) -> Result<R, IoError> + Sync),
    out: &mut dyn RecordWriter,
) -> Result<Vec<R>, IoError> {
    if chroms.is_empty() {
        return Ok(Vec::new());
    }

    if threads.get() == 1 {
        let mut results = Vec::with_capacity(chroms.len());
        for chrom in &chroms {
            let mut sink = ChunkSink::new(tmp_dir);
            let r = work(chrom, &mut sink)?;
            sink.drain_to(out)?;
            results.push(r);
        }
        return Ok(results);
    }

    std::thread::scope(|s| {
        let num_workers = threads.get().min(chroms.len());
        let (job_tx, job_rx) = crossbeam_channel::bounded::<(usize, Vec<u8>)>(threads.get());
        let (result_tx, result_rx) =
            crossbeam_channel::bounded::<(usize, Result<(R, ChunkSink), IoError>)>(threads.get());
        let cancelled = Arc::new(AtomicBool::new(false));

        for _ in 0..num_workers {
            let job_rx = job_rx.clone();
            let result_tx = result_tx.clone();
            let cancelled = cancelled.clone();
            s.spawn(move || {
                while let Ok((idx, chrom)) = job_rx.recv() {
                    if cancelled.load(Ordering::Relaxed) {
                        break;
                    }
                    let mut sink = ChunkSink::new(tmp_dir);
                    let res = work(&chrom, &mut sink);
                    if cancelled.load(Ordering::Relaxed) {
                        break;
                    }
                    let to_send = match res {
                        Ok(r) => Ok((r, sink)),
                        Err(e) => Err(e),
                    };
                    let _ = result_tx.send((idx, to_send));
                }
            });
        }
        drop(result_tx);

        let mut next_to_launch = 0;
        let mut next_to_drain = 0;
        let mut pending: BTreeMap<usize, (R, ChunkSink)> = BTreeMap::new();
        let mut results = Vec::with_capacity(chroms.len());
        let mut first_error = None;

        while next_to_drain < chroms.len() {
            while next_to_launch < chroms.len()
                && (next_to_launch - next_to_drain) < threads.get()
                && first_error.is_none()
            {
                let _ = job_tx.send((next_to_launch, chroms[next_to_launch].clone()));
                next_to_launch += 1;
            }

            let Ok((idx, res)) = result_rx.recv() else {
                break;
            };

            match res {
                Ok((r, sink)) => {
                    if first_error.is_none() {
                        pending.insert(idx, (r, sink));
                        while let Some((r, mut sink)) = pending.remove(&next_to_drain) {
                            if let Err(e) = sink.drain_to(out) {
                                first_error = Some(e);
                                cancelled.store(true, Ordering::Relaxed);
                                break;
                            }
                            results.push(r);
                            next_to_drain += 1;
                        }
                    }
                }
                Err(e) => {
                    if first_error.is_none() {
                        first_error = Some(e);
                        cancelled.store(true, Ordering::Relaxed);
                    }
                }
            }

            if first_error.is_some() {
                break;
            }
        }

        if let Some(e) = first_error {
            cancelled.store(true, Ordering::Relaxed);
            drop(job_tx);
            drop(pending);
            return Err(e);
        }

        drop(job_tx);
        Ok(results)
    })
}
