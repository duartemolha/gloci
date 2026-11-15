//! The waits of a worker pool sleep (`gloci_io::handoff`).
//!
//! This test has a test binary to itself. It keeps every CPU busy and counts how often the
//! kernel switches its two threads out while they could run; threads of tests running beside
//! it wake up on the same CPUs and add to that count.
#![cfg(target_os = "linux")]

use gloci_io::handoff::Pool;
use std::io;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::time::Duration;

/// How often the kernel has switched this thread out while it could run.
fn switched_out() -> u64 {
    let status = std::fs::read_to_string("/proc/thread-self/status").unwrap();
    let line = status
        .lines()
        .find_map(|l| l.strip_prefix("nonvoluntary_ctxt_switches:"))
        .expect("the kernel reports context switches");
    line.trim().parse().unwrap()
}

/// Sets its flag when it is dropped, however the test ends: the busy threads watch it.
struct StopOnDrop(Arc<AtomicBool>);

impl Drop for StopOnDrop {
    fn drop(&mut self) {
        self.0.store(true, Ordering::Relaxed);
    }
}

/// Leaves the worker's count of switches in its state.
fn work(switched: &mut Arc<AtomicU64>, _: ()) -> io::Result<()> {
    // Long enough for the driver to be waiting for this result.
    std::thread::sleep(Duration::from_millis(2));
    switched.store(switched_out(), Ordering::Relaxed);
    Ok(())
}

/// The driver that waits for a result, and the worker that waits for a job, sleep until they
/// are woken; neither retries around `sched_yield`. With every CPU busy, a thread that
/// yields is switched out at each yield while it could still run, and the kernel counts
/// that; one that sleeps is not. (With four yields before each wait, as a crossbeam channel
/// makes, the driver is switched out about 240 times in these 100 waits and the worker
/// about 350 times.)
#[test]
fn the_driver_and_the_workers_sleep_while_they_wait() {
    let rounds = 100u64;
    let stop = Arc::new(AtomicBool::new(false));
    // Dropped however the test ends: a failed assertion that left the busy threads running
    // would keep every CPU busy for the rest of the test run.
    let stopping = StopOnDrop(stop.clone());
    let cpus = std::thread::available_parallelism().map_or(4, |n| n.get());
    let busy: Vec<_> = (0..cpus)
        .map(|_| {
            let stop = stop.clone();
            std::thread::spawn(move || {
                while !stop.load(Ordering::Relaxed) {
                    std::hint::spin_loop();
                }
            })
        })
        .collect();
    let worker_switched = Arc::new(AtomicU64::new(0));
    let mut pool = Pool::new("test", vec![worker_switched.clone()], work).unwrap();
    // The first job leaves the worker's count before the waits that are measured.
    pool.submit(());
    pool.take().unwrap().unwrap();
    let worker_before = worker_switched.load(Ordering::Relaxed);
    let driver_before = switched_out();
    for _ in 0..rounds {
        // Long enough for the worker to be waiting for this job.
        std::thread::sleep(Duration::from_millis(2));
        pool.submit(());
        pool.take().unwrap().unwrap();
    }
    let driver = switched_out() - driver_before;
    let worker = worker_switched.load(Ordering::Relaxed) - worker_before;
    drop(stopping);
    for b in busy {
        b.join().unwrap();
    }
    assert!(
        driver < rounds / 2,
        "the waiting driver was switched out {driver} times while runnable in {rounds} waits"
    );
    assert!(
        worker < rounds / 2,
        "the waiting worker was switched out {worker} times while runnable in {rounds} waits"
    );
}
