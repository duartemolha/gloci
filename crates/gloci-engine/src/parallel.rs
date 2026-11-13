//! Per-chromosome parallel driver for indexed inputs (D41).

use crate::chunk::ChunkSink;
use gloci_io::{Compression, IndexKind, IoError, RecordWriter, detect_compression, find_index};
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::Read;
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Returns true when every path in `paths` is a regular BGZF file with a `.tbi`
/// or `.csi` index (or a `.bai`, for a BAM), or a CRAM with a `.crai`, and no path is `"-"`. Returns false if `paths` is empty.
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
        // A CRAM is indexed by a `.crai`.
        if head[..n].starts_with(b"CRAM") {
            if find_index(p).is_none_or(|i| i.kind != IndexKind::Crai) {
                return false;
            }
            continue;
        }
        if detect_compression(&head[..n]) != Compression::Bgzf {
            return false;
        }
        if find_index(p).is_none() {
            return false;
        }
    }
    true
}

/// The memory every in-flight chunk together may hold before spilling to
/// `--tmp-dir` (D41 as amended): each chunk gets this divided by the chunks
/// in flight, at least [`MIN_CHUNK_MEM`].
pub const CHUNK_MEM_BUDGET: usize = 64 << 20;
/// The most chunks in flight, whatever the thread count: each one is a
/// sweep with its own read buffers. Measured, indexed and plain-file runs
/// peak near 64 MiB plus about 20 MiB per chunk in flight (more with many
/// `-b` files, each of which has its own reader per chunk). Threads beyond
/// it still decompress and compress BGZF.
pub const MAX_CHUNKS_IN_FLIGHT: usize = 8;
/// The least a chunk keeps in memory.
pub const MIN_CHUNK_MEM: usize = 1 << 20;

/// TEST HOOK, not a user setting: overrides each chunk's in-memory limit in
/// bytes, so tests can force chunks to spill. Read by [`run_chunks`].
pub const TEST_CHUNK_MEM_ENV: &str = "GLOCI_TEST_CHUNK_MEM_LIMIT";

/// The in-memory limit of one chunk with `in_flight` chunks in flight.
pub fn chunk_mem_limit(in_flight: usize) -> usize {
    if let Some(n) = std::env::var(TEST_CHUNK_MEM_ENV)
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
    {
        return n;
    }
    (CHUNK_MEM_BUDGET / in_flight.max(1)).max(MIN_CHUNK_MEM)
}

/// A worker's panic as an error (a panic is always a bug).
fn panic_error(payload: Box<dyn std::any::Any + Send>) -> IoError {
    let msg = payload
        .downcast_ref::<&str>()
        .map(|s| s.to_string())
        .or_else(|| payload.downcast_ref::<String>().cloned())
        .unwrap_or_else(|| "unknown panic".to_string());
    IoError::Internal(format!("a worker thread panicked: {msg}"))
}

/// Run `work` for each chunk job in `jobs` (a run of consecutive
/// chromosomes, say) across `min(threads, MAX_CHUNKS_IN_FLIGHT)` workers,
/// with that many chunks in flight, which share [`CHUNK_MEM_BUDGET`] of memory before they
/// spill to `tmp_dir`. Finished chunks are drained to `out` in strict job
/// order, and results come back in job order. When chunks fail, the error
/// of the first failing job in job order is returned (whatever finished
/// first), so the message is the same for every thread count; a panicking
/// worker becomes an error.
#[allow(clippy::type_complexity)]
pub fn run_chunks<J: Clone + Send + Sync, R: Send>(
    jobs: Vec<J>,
    threads: NonZeroUsize,
    tmp_dir: &Path,
    work: &(dyn Fn(&J, &mut ChunkSink) -> Result<R, IoError> + Sync),
    out: &mut dyn RecordWriter,
) -> Result<Vec<R>, IoError> {
    if jobs.is_empty() {
        return Ok(Vec::new());
    }
    let in_flight = threads.get().min(MAX_CHUNKS_IN_FLIGHT).min(jobs.len());
    let mem_limit = chunk_mem_limit(in_flight);
    // One buffer per chunk in flight, reused from chunk to chunk: allocating
    // and freeing a buffer per chunk leaves freed buffers resident in the
    // allocator for a while, beyond the budget.
    let mut buffers: Vec<Vec<u8>> = (0..in_flight).map(|_| Vec::new()).collect();
    let cancel = gloci_io::cancel::Watch::current();
    let run_one = |job: &J, buf: Vec<u8>| -> Result<(R, ChunkSink), IoError> {
        let mut sink = ChunkSink::with_buffer(tmp_dir, mem_limit, buf);
        let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| work(job, &mut sink)))
            .map_err(panic_error)??;
        Ok((r, sink))
    };

    if in_flight == 1 {
        let mut results = Vec::with_capacity(jobs.len());
        let mut buf = buffers.pop().unwrap_or_default();
        for job in &jobs {
            cancel.check()?;
            let (r, mut sink) = run_one(job, buf)?;
            sink.drain_to(out)?;
            buf = sink.into_buffer();
            results.push(r);
        }
        return Ok(results);
    }

    // The smallest index of a failed job so far; later jobs are skipped.
    let fail_at = AtomicUsize::new(usize::MAX);
    std::thread::scope(|s| {
        let (job_tx, job_rx) = crossbeam_channel::bounded::<(usize, J, Vec<u8>)>(in_flight);
        // `None`: the job was skipped because an earlier one failed.
        type Done<R> = Option<Result<(R, ChunkSink), IoError>>;
        let (result_tx, result_rx) = crossbeam_channel::bounded::<(usize, Done<R>)>(in_flight);

        for _ in 0..in_flight {
            let job_rx = job_rx.clone();
            let result_tx = result_tx.clone();
            let (run_one, fail_at) = (&run_one, &fail_at);
            let cancel = cancel.clone();
            // CRAM 3.1's codecs (order-1 rANS, the name tokeniser) need more stack than a
            // spawned thread's default 2 MiB, and chunk workers decode CRAM.
            let spawned = std::thread::Builder::new()
                .stack_size(gloci_io::WORKER_STACK)
                .spawn_scoped(s, move || {
                    let _run = cancel.enter();
                    while let Ok((idx, job, buf)) = job_rx.recv() {
                        let res = if idx > fail_at.load(Ordering::Relaxed) {
                            None
                        } else {
                            let r = match cancel.check() {
                                Err(e) => Err(e),
                                Ok(()) => run_one(&job, buf),
                            };
                            if r.is_err() {
                                fail_at.fetch_min(idx, Ordering::Relaxed);
                            }
                            Some(r)
                        };
                        if result_tx.send((idx, res)).is_err() {
                            break;
                        }
                    }
                });
            if let Err(e) = spawned {
                return Err(IoError::Internal(format!(
                    "cannot start a worker thread: {e}"
                )));
            }
        }
        drop(result_tx);
        drop(job_rx);

        let mut next_to_launch = 0;
        let mut next_to_drain = 0;
        let mut pending: BTreeMap<usize, (R, ChunkSink)> = BTreeMap::new();
        let mut results = Vec::with_capacity(jobs.len());
        // The failing job with the smallest index so far, and its error.
        let mut first_error: Option<(usize, IoError)> = None;

        loop {
            // Done when every job before the first failing one has reported:
            // only those can hold an earlier error.
            let need = first_error.as_ref().map_or(jobs.len(), |(i, _)| *i);
            if next_to_drain >= need {
                break;
            }
            while first_error.is_none()
                && next_to_launch < jobs.len()
                && next_to_launch - next_to_drain < in_flight
            {
                // At most `in_flight` chunks are launched and not yet
                // drained, so a buffer is always free here.
                let buf = buffers.pop().unwrap_or_default();
                if job_tx
                    .send((next_to_launch, jobs[next_to_launch].clone(), buf))
                    .is_err()
                {
                    break;
                }
                next_to_launch += 1;
            }
            let Ok((idx, res)) = result_rx.recv() else {
                break;
            };
            match res {
                Some(Ok(done)) => {
                    pending.insert(idx, done);
                }
                Some(Err(e)) if first_error.as_ref().is_none_or(|(i, _)| idx < *i) => {
                    first_error = Some((idx, e));
                }
                // No result, or an error later than the one already kept.
                None | Some(Err(_)) => {}
            }
            while let Some((r, mut sink)) = pending.remove(&next_to_drain) {
                if first_error.is_none()
                    && let Err(e) = sink.drain_to(out)
                {
                    fail_at.fetch_min(next_to_drain, Ordering::Relaxed);
                    first_error = Some((next_to_drain, e));
                }
                buffers.push(sink.into_buffer());
                results.push(r);
                next_to_drain += 1;
            }
        }
        drop(job_tx);
        drop(pending);
        match first_error {
            Some((_, e)) => Err(e),
            None if results.len() == jobs.len() => Ok(results),
            // Unreachable while workers run under `catch_unwind`: a lost
            // result must never pass as a complete, shorter output.
            None => Err(IoError::Internal(format!(
                "parallel run returned {} results for {} chunks; this is a bug, \
please report it",
                results.len(),
                jobs.len()
            ))),
        }
    })
}
