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
//! The selection stage of the `vcf` commands uses [`Handoff`] with scoped workers that borrow
//! the stage, and the guards [`Driving`] and [`Present`].

use std::collections::{BTreeMap, VecDeque};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Condvar, Mutex, MutexGuard};

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
