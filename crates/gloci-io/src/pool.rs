//! A pool of workers whose results are released in submission order (D110): the CRAM decode
//! pool and the `.crai` builder use it.
//!
//! - Workers have 16 MiB stacks: noodles-cram's order-1 rANS decoders and its name tokeniser
//!   keep tables of 256 KB and more on the stack, past the 2 MiB a spawned thread gets.
//! - Dropping the pool cancels the jobs not yet started, so a reader that stops early (`| head`)
//!   does not wait for the whole queue to be decoded.
//! - A job that panics (a bug; release builds abort) comes back as an error, so the reader
//!   waiting for it does not hang.

use crate::error::IoError;
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, mpsc};
use std::thread::JoinHandle;

/// The stack of a pool worker.
pub const WORKER_STACK: usize = 16 << 20;

type Job<J> = (u64, J);

/// Work done on up to `threads` workers, results released in submission order.
pub(crate) struct Ordered<J: Send + 'static, R: Send + 'static> {
    jobs: Option<mpsc::Sender<Job<J>>>,
    results: mpsc::Receiver<(u64, Result<R, IoError>)>,
    workers: Vec<JoinHandle<()>>,
    /// Set on drop: workers skip the jobs left.
    cancel: Arc<AtomicBool>,
    pending: BTreeMap<u64, Result<R, IoError>>,
    submitted: u64,
    released: u64,
}

impl<J: Send + 'static, R: Send + 'static> Ordered<J, R> {
    /// One worker thread per state in `states`, each running `work` on the jobs it takes.
    pub fn new<W: Send + 'static>(
        states: Vec<W>,
        work: fn(&mut W, J) -> Result<R, IoError>,
    ) -> Result<Self, IoError> {
        let (job_tx, job_rx) = mpsc::channel::<Job<J>>();
        let job_rx = Arc::new(Mutex::new(job_rx));
        let (res_tx, results) = mpsc::channel();
        let mut workers = Vec::with_capacity(states.len());
        let cancel = Arc::new(AtomicBool::new(false));
        for (i, mut state) in states.into_iter().enumerate() {
            let rx = job_rx.clone();
            let tx = res_tx.clone();
            let cancel = cancel.clone();
            let h = std::thread::Builder::new()
                .name(format!("gloci-pool-{i}"))
                .stack_size(WORKER_STACK)
                .spawn(move || {
                    loop {
                        let job = match rx.lock() {
                            Ok(g) => g.recv(),
                            Err(_) => break,
                        };
                        let Ok((seq, j)) = job else { break };
                        if cancel.load(Ordering::Relaxed) {
                            break;
                        }
                        let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                            work(&mut state, j)
                        }))
                        .unwrap_or_else(|p| {
                            let msg = p
                                .downcast_ref::<&str>()
                                .map(|s| s.to_string())
                                .or_else(|| p.downcast_ref::<String>().cloned())
                                .unwrap_or_default();
                            Err(IoError::Internal(format!(
                                "a worker thread panicked: {msg}"
                            )))
                        });
                        if tx.send((seq, r)).is_err() {
                            break;
                        }
                    }
                })
                .map_err(|e| IoError::Internal(format!("cannot start a worker thread: {e}")))?;
            workers.push(h);
        }
        Ok(Ordered {
            jobs: Some(job_tx),
            results,
            workers,
            cancel,
            pending: BTreeMap::new(),
            submitted: 0,
            released: 0,
        })
    }

    /// Jobs submitted and not yet released.
    pub fn in_flight(&self) -> u64 {
        self.submitted - self.released
    }

    pub fn submit(&mut self, j: J) {
        let seq = self.submitted;
        self.submitted += 1;
        let sent = self
            .jobs
            .as_ref()
            .is_some_and(|tx| tx.send((seq, j)).is_ok());
        if !sent {
            self.pending.insert(
                seq,
                Err(IoError::Internal("the decoder threads stopped".into())),
            );
        }
    }

    /// Take a slot in the order for a result known already (a read error).
    pub fn submit_result(&mut self, r: Result<R, IoError>) {
        let seq = self.submitted;
        self.submitted += 1;
        self.pending.insert(seq, r);
    }

    /// The next result in order, waiting for it; `None` when nothing is in flight.
    pub fn next(&mut self) -> Option<Result<R, IoError>> {
        if self.released == self.submitted {
            return None;
        }
        loop {
            if let Some(r) = self.pending.remove(&self.released) {
                self.released += 1;
                return Some(r);
            }
            match self.results.recv() {
                Ok((seq, r)) => {
                    self.pending.insert(seq, r);
                }
                Err(_) => {
                    self.released += 1;
                    return Some(Err(IoError::Internal(
                        "a decoder thread stopped unexpectedly".into(),
                    )));
                }
            }
        }
    }
}

impl<J: Send + 'static, R: Send + 'static> Drop for Ordered<J, R> {
    fn drop(&mut self) {
        self.cancel.store(true, Ordering::Relaxed);
        self.jobs = None;
        for h in self.workers.drain(..) {
            let _ = h.join();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::AtomicUsize;
    use std::time::Duration;

    fn slow(done: &mut Arc<AtomicUsize>, j: u64) -> Result<u64, IoError> {
        std::thread::sleep(Duration::from_millis((j % 3) * 5));
        if j == 7 {
            panic!("job 7");
        }
        done.fetch_add(1, Ordering::Relaxed);
        Ok(j * 2)
    }

    #[test]
    fn results_come_in_order_and_a_panic_is_an_error() {
        let done = Arc::new(AtomicUsize::new(0));
        let mut p = Ordered::new(vec![done.clone(), done.clone(), done.clone()], slow).unwrap();
        for j in 0..20 {
            p.submit(j);
        }
        for j in 0..20 {
            match p.next().unwrap() {
                Ok(v) => assert_eq!(v, j * 2),
                Err(e) => {
                    assert_eq!(j, 7);
                    assert!(e.to_string().contains("panicked: job 7"), "{e}");
                }
            }
        }
        assert!(p.next().is_none());
    }

    /// The workers' stacks hold what noodles' codecs put there (order-1 rANS and the name
    /// tokeniser need more than a spawned thread's default 2 MiB, in release builds with fat
    /// LTO); this bites in debug builds too: a job with 6 MiB on its stack must run.
    #[test]
    fn workers_have_large_stacks() {
        const { assert!(WORKER_STACK >= 8 << 20) };
        fn deep(_: &mut (), j: u64) -> Result<u64, IoError> {
            let buf = std::hint::black_box([j as u8; 6 << 20]);
            Ok(buf[(j as usize) % buf.len()] as u64 + buf[buf.len() - 1] as u64)
        }
        let mut p = Ordered::new(vec![(), ()], deep).unwrap();
        for j in 0..4 {
            p.submit(j);
        }
        for j in 0..4 {
            assert_eq!(p.next().unwrap().unwrap(), 2 * j);
        }
    }

    /// Dropping the pool skips the jobs not started (a reader stopped by `| head`).
    #[test]
    fn drop_cancels_the_jobs_left() {
        fn wait(done: &mut Arc<AtomicUsize>, j: u64) -> Result<u64, IoError> {
            std::thread::sleep(Duration::from_millis(20));
            done.fetch_add(1, Ordering::Relaxed);
            Ok(j)
        }
        let done = Arc::new(AtomicUsize::new(0));
        let mut p = Ordered::new(vec![done.clone(), done.clone()], wait).unwrap();
        for j in 0..200 {
            p.submit(j);
        }
        assert_eq!(p.next().unwrap().unwrap(), 0);
        drop(p);
        let n = done.load(Ordering::Relaxed);
        assert!(n < 20, "{n} of 200 jobs ran after the pool was dropped");
    }
}
