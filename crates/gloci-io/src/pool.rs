//! A pool of workers whose results are released in submission order (D110): the CRAM decode
//! pool and the `.crai` builder use it.

use crate::error::IoError;
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex, mpsc};
use std::thread::JoinHandle;

type Job<J> = (u64, J);

/// Work done on up to `threads` workers, results released in submission order.
pub(crate) struct Ordered<J: Send + 'static, R: Send + 'static> {
    jobs: Option<mpsc::Sender<Job<J>>>,
    results: mpsc::Receiver<(u64, Result<R, IoError>)>,
    workers: Vec<JoinHandle<()>>,
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
        for (i, mut state) in states.into_iter().enumerate() {
            let rx = job_rx.clone();
            let tx = res_tx.clone();
            let h = std::thread::Builder::new()
                .name(format!("gloci-cram-{i}"))
                .spawn(move || {
                    loop {
                        let job = match rx.lock() {
                            Ok(g) => g.recv(),
                            Err(_) => break,
                        };
                        let Ok((seq, j)) = job else { break };
                        if tx.send((seq, work(&mut state, j))).is_err() {
                            break;
                        }
                    }
                })
                .map_err(|e| IoError::Internal(format!("cannot start a decoder thread: {e}")))?;
            workers.push(h);
        }
        Ok(Ordered {
            jobs: Some(job_tx),
            results,
            workers,
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
        self.jobs = None;
        for h in self.workers.drain(..) {
            let _ = h.join();
        }
    }
}
