//! Work handed from one driving thread to workers and back in order, with every wait asleep.
//!
//! [`Handoff`] is where a driver and its workers meet: a queue of jobs and a table of the
//! results, under one mutex with two condition variables. The driver queues jobs and takes
//! each result by its number, so results leave in the order their jobs went in and the
//! earliest error is the one the driver meets first. How much is in flight is the driver's to
//! bound: it takes a result before it queues past its limit.
//!
//! Whoever waits here sleeps until it is woken. A crossbeam channel's blocking calls first
//! retry around `sched_yield`, and so does a rayon worker with no job; on Linux a yield on a
//! CPU that another thread wants gives up a whole scheduler slice (milliseconds). With more
//! runnable threads than cores every hand-off then costs a slice or more, which is far more
//! than a BGZF block or a batch of records takes to work on.
//!
//! - The selection stage of the `vcf` commands uses [`Handoff`] with scoped workers that
//!   borrow the stage, and the guards [`Driving`] and [`Present`].
//! - [`Pool`] owns its worker threads, for a driver that outlives any scope: the BGZF reader
//!   and writer (`bgzf.rs`). Dropping it cancels the jobs not yet started, and a job that
//!   panics comes back as an error.

use std::collections::{BTreeMap, VecDeque};
use std::io;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Condvar, Mutex, MutexGuard};
use std::thread::JoinHandle;

/// A result's place in the order: the number of the job it answers.
pub trait Sequenced {
    fn seq(&self) -> u64;
}

impl<T> Sequenced for (u64, T) {
    fn seq(&self) -> u64 {
        self.0
    }
}

/// Where the driver and the workers of a parallel run meet: the jobs to work on, and the
/// results by number. One thread drives (it alone calls [`Handoff::take`]).
pub struct Handoff<J, D> {
    state: Mutex<State<J, D>>,
    /// A job was queued for a worker that sleeps, or the hand-off closed.
    job_ready: Condvar,
    /// The result the driver sleeps for came back, or a worker left.
    done_ready: Condvar,
}

struct State<J, D> {
    jobs: VecDeque<J>,
    done: BTreeMap<u64, D>,
    /// Workers still running.
    workers: usize,
    /// Workers asleep until a job comes.
    idle: usize,
    /// The result the driver is asleep for.
    wanted: Option<u64>,
    /// The driver takes no more results.
    closed: bool,
}

impl<J, D> Handoff<J, D> {
    /// A hand-off that `workers` workers will join, each holding a [`Present`].
    pub fn new(workers: usize) -> Handoff<J, D> {
        Handoff {
            state: Mutex::new(State {
                jobs: VecDeque::new(),
                done: BTreeMap::new(),
                workers,
                idle: 0,
                wanted: None,
                closed: false,
            }),
            job_ready: Condvar::new(),
            done_ready: Condvar::new(),
        }
    }

    /// The state. Nothing panics while holding the lock (a worker's panic is caught around
    /// its job), so a poisoned lock still holds a sound state.
    fn lock(&self) -> MutexGuard<'_, State<J, D>> {
        self.state.lock().unwrap_or_else(|e| e.into_inner())
    }

    /// Queue `job` for a worker. False, and the job is dropped, when no worker is left to
    /// take it.
    pub fn send(&self, job: J) -> bool {
        let mut s = self.lock();
        if s.workers == 0 {
            return false;
        }
        s.jobs.push_back(job);
        // A busy worker finds the job when it asks for its next one.
        let wake = s.idle > 0;
        drop(s);
        if wake {
            self.job_ready.notify_one();
        }
        true
    }

    /// A worker's next job, waiting for one; `None` once the driver has closed the hand-off.
    pub fn next_job(&self) -> Option<J> {
        let mut s = self.lock();
        loop {
            if s.closed {
                return None;
            }
            if let Some(job) = s.jobs.pop_front() {
                return Some(job);
            }
            s.idle += 1;
            s = self.job_ready.wait(s).unwrap_or_else(|e| e.into_inner());
            s.idle -= 1;
        }
    }

    /// Result `seq`, waiting for it; `None` when it is not back and no worker is left to
    /// bring it.
    pub fn take(&self, seq: u64) -> Option<D> {
        let mut s = self.lock();
        loop {
            if let Some(done) = s.done.remove(&seq) {
                s.wanted = None;
                return Some(done);
            }
            if s.workers == 0 {
                s.wanted = None;
                return None;
            }
            s.wanted = Some(seq);
            s = self.done_ready.wait(s).unwrap_or_else(|e| e.into_inner());
        }
    }

    /// Result `seq` if it is back, without waiting.
    pub fn try_take(&self, seq: u64) -> Option<D> {
        self.lock().done.remove(&seq)
    }

    /// The driver is done: waiting workers leave, busy ones after their job.
    pub fn close(&self) {
        self.lock().closed = true;
        self.job_ready.notify_all();
    }
}

impl<J, D: Sequenced> Handoff<J, D> {
    /// A worker hands a result back.
    pub fn finish(&self, done: D) {
        let seq = done.seq();
        let mut s = self.lock();
        s.done.insert(seq, done);
        // A result the driver is not asleep for waits in the table until it is asked for.
        let wake = s.wanted == Some(seq);
        drop(s);
        if wake {
            self.done_ready.notify_one();
        }
    }
}

/// The driver's end of the hand-off: dropped when the driver leaves, with its result or
/// unwinding from a panic, so that the workers are told to stop and whoever waits for them
/// (a thread scope, a join) can end.
pub struct Driving<'a, J, D> {
    pub handoff: &'a Handoff<J, D>,
    /// Set on drop: a worker that has taken a job and not started it skips it.
    pub abort: &'a AtomicBool,
}

impl<J, D> Drop for Driving<'_, J, D> {
    fn drop(&mut self) {
        self.abort.store(true, Ordering::Relaxed);
        self.handoff.close();
    }
}

/// A worker's presence at the hand-off: dropped when the worker returns, however it returns,
/// so the driver never waits for a worker that is gone.
pub struct Present<'a, J, D>(pub &'a Handoff<J, D>);

impl<J, D> Drop for Present<'_, J, D> {
    fn drop(&mut self) {
        self.0.lock().workers -= 1;
        self.0.done_ready.notify_one();
    }
}

type PoolHandoff<J, R> = Handoff<(u64, J), (u64, io::Result<R>)>;

/// Worker threads that take jobs from one driver and give their results back in submission
/// order.
///
/// - The driver bounds what is in flight: it calls [`Pool::take`] before it submits past its
///   limit ([`Pool::in_flight`]).
/// - Dropping the pool cancels the jobs not yet started and waits for the workers.
/// - A job that panics (a bug; release builds abort) comes back as an error, and its worker
///   takes no more jobs.
pub struct Pool<J: Send + 'static, R: Send + 'static> {
    handoff: Arc<PoolHandoff<J, R>>,
    abort: Arc<AtomicBool>,
    workers: Vec<JoinHandle<()>>,
    submitted: u64,
    released: u64,
}

impl<J: Send + 'static, R: Send + 'static> Pool<J, R> {
    /// One worker thread per state in `states`, named `name-N`, each running `work` on the
    /// jobs it takes.
    pub fn new<W: Send + 'static>(
        name: &str,
        states: Vec<W>,
        work: fn(&mut W, J) -> io::Result<R>,
    ) -> io::Result<Self> {
        let mut pool = Pool {
            handoff: Arc::new(Handoff::new(states.len())),
            abort: Arc::new(AtomicBool::new(false)),
            workers: Vec::with_capacity(states.len()),
            submitted: 0,
            released: 0,
        };
        // The run that builds the pool owns its workers (D91).
        let run = crate::cancel::Watch::current();
        for (i, mut state) in states.into_iter().enumerate() {
            let handoff = pool.handoff.clone();
            let abort = pool.abort.clone();
            let run = run.clone();
            // A worker that cannot start ends this function, and the pool's drop lets the
            // workers already started go.
            let worker = std::thread::Builder::new()
                .name(format!("{name}-{i}"))
                .spawn(move || {
                    let _present = Present(&handoff);
                    let _run = run.enter();
                    while let Some((seq, job)) = handoff.next_job() {
                        if abort.load(Ordering::Relaxed) {
                            break;
                        }
                        let ran = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                            work(&mut state, job)
                        }));
                        let panicked = ran.is_err();
                        // A panic must reach the driver, which waits for this result.
                        let result = ran.unwrap_or_else(|p| {
                            let msg = p
                                .downcast_ref::<&str>()
                                .map(|s| s.to_string())
                                .or_else(|| p.downcast_ref::<String>().cloned())
                                .unwrap_or_default();
                            Err(io::Error::other(format!(
                                "a worker thread panicked: {msg}; this is a bug in gloci, \
please report it"
                            )))
                        });
                        handoff.finish((seq, result));
                        if panicked {
                            break;
                        }
                    }
                })?;
            pool.workers.push(worker);
        }
        Ok(pool)
    }

    /// Jobs submitted whose results are not yet taken.
    pub fn in_flight(&self) -> u64 {
        self.submitted - self.released
    }

    /// Queue `job`. With no worker left it is dropped, and its place in the order holds an
    /// error.
    pub fn submit(&mut self, job: J) {
        let seq = self.submitted;
        self.submitted += 1;
        self.handoff.send((seq, job));
    }

    /// The next result in order, waiting for it; `None` when nothing is in flight.
    pub fn take(&mut self) -> Option<io::Result<R>> {
        if self.released == self.submitted {
            return None;
        }
        let taken = self.handoff.take(self.released);
        self.released += 1;
        Some(match taken {
            Some((_, result)) => result,
            None => Err(io::Error::other(
                "the worker threads stopped before their work was done; this is a bug in \
gloci, please report it",
            )),
        })
    }

    /// The next result in order if it is back, without waiting.
    pub fn try_take(&mut self) -> Option<io::Result<R>> {
        if self.released == self.submitted {
            return None;
        }
        let (_, result) = self.handoff.try_take(self.released)?;
        self.released += 1;
        Some(result)
    }
}

impl<J: Send + 'static, R: Send + 'static> Drop for Pool<J, R> {
    fn drop(&mut self) {
        drop(Driving {
            handoff: &*self.handoff,
            abort: &self.abort,
        });
        for worker in self.workers.drain(..) {
            let _ = worker.join();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::AtomicUsize;
    use std::time::Duration;

    fn slow(done: &mut Arc<AtomicUsize>, j: u64) -> io::Result<u64> {
        std::thread::sleep(Duration::from_millis((j % 3) * 5));
        if j == 7 {
            panic!("job 7");
        }
        if j == 11 {
            return Err(io::Error::other("job 11"));
        }
        done.fetch_add(1, Ordering::Relaxed);
        Ok(j * 2)
    }

    #[test]
    fn results_come_in_order_and_a_panic_is_an_error() {
        let done = Arc::new(AtomicUsize::new(0));
        let mut p =
            Pool::new("test", vec![done.clone(), done.clone(), done.clone()], slow).unwrap();
        for j in 0..20 {
            p.submit(j);
        }
        assert_eq!(p.in_flight(), 20);
        for j in 0..20 {
            match p.take().unwrap() {
                Ok(v) => assert_eq!(v, j * 2),
                Err(e) if j == 7 => assert!(e.to_string().contains("panicked: job 7"), "{e}"),
                Err(e) => assert_eq!((j, e.to_string().as_str()), (11, "job 11")),
            }
        }
        assert_eq!(p.in_flight(), 0);
        assert!(p.take().is_none());
        assert!(p.try_take().is_none());
        assert_eq!(done.load(Ordering::Relaxed), 18);
    }

    /// The worker that panicked takes no more jobs: with one worker, the jobs after the panic
    /// are answered with an error instead of being waited for.
    #[test]
    fn the_jobs_after_a_lone_workers_panic_are_errors() {
        let done = Arc::new(AtomicUsize::new(0));
        let mut p = Pool::new("test", vec![done.clone()], slow).unwrap();
        for j in 5..10 {
            p.submit(j);
        }
        assert_eq!(p.take().unwrap().unwrap(), 10);
        assert_eq!(p.take().unwrap().unwrap(), 12);
        let e = p.take().unwrap().unwrap_err();
        assert!(e.to_string().contains("panicked: job 7"), "{e}");
        for _ in 8..10 {
            let e = p.take().unwrap().unwrap_err();
            assert!(e.to_string().contains("stopped"), "{e}");
        }
        assert!(p.take().is_none());
        // A job submitted now has nobody to take it.
        p.submit(1);
        assert!(p.take().unwrap().is_err());
    }

    /// `try_take` hands out a result only when it is the next in order and back.
    #[test]
    fn try_take_does_not_wait_and_keeps_the_order() {
        fn timed(_: &mut (), ms: u64) -> io::Result<u64> {
            std::thread::sleep(Duration::from_millis(ms));
            Ok(ms)
        }
        let mut p = Pool::new("test", vec![(), ()], timed).unwrap();
        p.submit(150);
        p.submit(0);
        // The second job is back long before the first, and is not handed out ahead of it.
        std::thread::sleep(Duration::from_millis(50));
        assert!(p.try_take().is_none());
        assert_eq!(p.in_flight(), 2);
        assert_eq!(p.take().unwrap().unwrap(), 150);
        assert_eq!(p.try_take().unwrap().unwrap(), 0);
        assert!(p.try_take().is_none());
    }

    /// Dropping the pool skips the jobs not started (a reader stopped by `| head`).
    #[test]
    fn drop_cancels_the_jobs_left() {
        fn wait(done: &mut Arc<AtomicUsize>, j: u64) -> io::Result<u64> {
            std::thread::sleep(Duration::from_millis(20));
            done.fetch_add(1, Ordering::Relaxed);
            Ok(j)
        }
        let done = Arc::new(AtomicUsize::new(0));
        let mut p = Pool::new("test", vec![done.clone(), done.clone()], wait).unwrap();
        for j in 0..200 {
            p.submit(j);
        }
        assert_eq!(p.take().unwrap().unwrap(), 0);
        drop(p);
        let n = done.load(Ordering::Relaxed);
        assert!(n < 20, "{n} of 200 jobs ran after the pool was dropped");
    }
}
