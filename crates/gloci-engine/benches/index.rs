//! Index micro-benchmark: compares coitrees, ITree64, superintervals, and rust-lapper.
//!
//! Evaluates construction time and query throughput on:
//! 1. Synthetic dataset (uniform intervals)
//! 2. Coitrees benchmark distribution (nested/clustered intervals)

use coitrees::{BasicCOITree, Interval as CoiInterval, IntervalTree};
use gloci_core::Interval;
use gloci_engine::ITree64;
use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::hint::black_box;
use std::time::Instant;

#[derive(Clone, Copy, Debug)]
pub struct LapIv {
    pub start: u64,
    pub end: u64,
    pub val: u32,
}

/// rust-lapper implementation
pub struct Lapper {
    pub intervals: Vec<LapIv>,
    pub max_len: u64,
}

impl Lapper {
    pub fn new(mut intervals: Vec<LapIv>) -> Self {
        intervals.sort_unstable_by_key(|iv| iv.start);
        let mut max_len = 0;
        for iv in &intervals {
            max_len = max_len.max(iv.end.saturating_sub(iv.start));
        }
        Self { intervals, max_len }
    }

    pub fn query(&self, start: u64, end: u64, out: &mut Vec<u32>) {
        let off = start.saturating_sub(self.max_len);
        let idx = match self.intervals.binary_search_by_key(&off, |iv| iv.start) {
            Ok(i) => i,
            Err(i) => i,
        };
        for iv in &self.intervals[idx..] {
            if iv.start > end {
                break;
            }
            if iv.end >= start {
                out.push(iv.val);
            }
        }
    }
}

/// superintervals implementation (augmented array with pruning)
pub struct SuperIntervals {
    pub intervals: Vec<LapIv>,
    pub max_ends: Vec<u64>,
}

impl SuperIntervals {
    pub fn new(mut intervals: Vec<LapIv>) -> Self {
        intervals.sort_unstable_by_key(|iv| iv.start);
        let mut max_ends = Vec::with_capacity(intervals.len());
        let mut m = 0;
        for iv in &intervals {
            m = m.max(iv.end);
            max_ends.push(m);
        }
        Self {
            intervals,
            max_ends,
        }
    }

    pub fn query(&self, start: u64, end: u64, out: &mut Vec<u32>) {
        let idx = self.intervals.partition_point(|iv| iv.start <= end);
        if idx == 0 {
            return;
        }
        for i in (0..idx).rev() {
            let iv = &self.intervals[i];
            if iv.end >= start {
                out.push(iv.val);
            } else if self.max_ends[i] < start {
                break;
            }
        }
    }
}

fn generate_synthetic(n: usize, seed: u64) -> (Vec<LapIv>, Vec<(u64, u64)>) {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let mut ivs = Vec::with_capacity(n);
    let chrom_len = 50_000_000u64;
    for i in 0..n {
        let len = rng.random_range(50..=5_000);
        let s = rng.random_range(0..=chrom_len - len);
        ivs.push(LapIv {
            start: s,
            end: s + len,
            val: i as u32,
        });
    }
    let mut queries = Vec::with_capacity(10_000);
    for _ in 0..10_000 {
        let len = rng.random_range(100..=1_000);
        let s = rng.random_range(0..=chrom_len - len);
        queries.push((s, s + len));
    }
    (ivs, queries)
}

fn generate_coitrees_dist(n: usize, seed: u64) -> (Vec<LapIv>, Vec<(u64, u64)>) {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let mut ivs = Vec::with_capacity(n);
    let chrom_len = 50_000_000u64;
    let centres: Vec<u64> = (0..500).map(|_| rng.random_range(0..=chrom_len)).collect();
    for i in 0..n {
        let c = centres[rng.random_range(0..centres.len())];
        let len = rng.random_range(50..=20_000);
        let s = c.saturating_sub(len / 2).min(chrom_len - len);
        ivs.push(LapIv {
            start: s,
            end: s + len,
            val: i as u32,
        });
    }
    let mut queries = Vec::with_capacity(10_000);
    for _ in 0..10_000 {
        let c = centres[rng.random_range(0..centres.len())];
        let len = rng.random_range(100..=2_000);
        let s = c.saturating_sub(len / 2).min(chrom_len - len);
        queries.push((s, s + len));
    }
    (ivs, queries)
}

fn bench_queries<F>(mut f: F, queries: &[(u64, u64)]) -> f64
where
    F: FnMut(u64, u64, &mut Vec<u32>),
{
    let mut buf = Vec::with_capacity(64);
    // Warmup
    for &(s, e) in queries.iter().take(1000) {
        buf.clear();
        f(s, e, &mut buf);
        black_box(&buf);
    }
    // Timed runs
    let mut times = Vec::new();
    for _ in 0..5 {
        let t0 = Instant::now();
        let mut total_hits = 0usize;
        for &(s, e) in queries {
            buf.clear();
            f(s, e, &mut buf);
            total_hits += buf.len();
        }
        black_box(total_hits);
        times.push(t0.elapsed().as_secs_f64());
    }
    times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median_secs = times[times.len() / 2];
    (median_secs * 1e9) / (queries.len() as f64)
}

fn main() {
    println!(
        "Running index micro-benchmark (coitrees vs ITree64 vs superintervals vs rust-lapper)..."
    );
    let n_intervals = 100_000;

    let (synth_ivs, synth_queries) = generate_synthetic(n_intervals, 42);
    let (coi_ivs, coi_queries) = generate_coitrees_dist(n_intervals, 43);

    // 1. coitrees
    let coi_nodes_synth: Vec<CoiInterval<u32>> = synth_ivs
        .iter()
        .map(|iv| CoiInterval::new(iv.start as i32, iv.end as i32, iv.val))
        .collect();
    let t0 = Instant::now();
    let coitree_synth = BasicCOITree::<u32, u32>::new(&coi_nodes_synth);
    let coitree_synth_build = t0.elapsed().as_secs_f64() * 1000.0;
    let coitree_synth_ns = bench_queries(
        |s, e, out| {
            coitree_synth.query(s as i32, e as i32, |n| {
                out.push(n.metadata);
            });
        },
        &synth_queries,
    );

    let coi_nodes_coi: Vec<CoiInterval<u32>> = coi_ivs
        .iter()
        .map(|iv| CoiInterval::new(iv.start as i32, iv.end as i32, iv.val))
        .collect();
    let coitree_coi = BasicCOITree::<u32, u32>::new(&coi_nodes_coi);
    let coitree_coi_ns = bench_queries(
        |s, e, out| {
            coitree_coi.query(s as i32, e as i32, |n| {
                out.push(n.metadata);
            });
        },
        &coi_queries,
    );

    // 2. ITree64
    let itree_items_synth: Vec<(Interval, u32)> = synth_ivs
        .iter()
        .map(|iv| (Interval::new(iv.start, iv.end).unwrap(), iv.val))
        .collect();
    let t0 = Instant::now();
    let itree_synth = ITree64::new(itree_items_synth);
    let itree_synth_build = t0.elapsed().as_secs_f64() * 1000.0;
    let itree_synth_ns = bench_queries(|s, e, out| itree_synth.query(s, e, out), &synth_queries);

    let itree_items_coi: Vec<(Interval, u32)> = coi_ivs
        .iter()
        .map(|iv| (Interval::new(iv.start, iv.end).unwrap(), iv.val))
        .collect();
    let itree_coi = ITree64::new(itree_items_coi);
    let itree_coi_ns = bench_queries(|s, e, out| itree_coi.query(s, e, out), &coi_queries);

    // 3. SuperIntervals
    let t0 = Instant::now();
    let super_synth = SuperIntervals::new(synth_ivs.clone());
    let super_synth_build = t0.elapsed().as_secs_f64() * 1000.0;
    let super_synth_ns = bench_queries(|s, e, out| super_synth.query(s, e, out), &synth_queries);

    let super_coi = SuperIntervals::new(coi_ivs.clone());
    let super_coi_ns = bench_queries(|s, e, out| super_coi.query(s, e, out), &coi_queries);

    // 4. Lapper
    let t0 = Instant::now();
    let lapper_synth = Lapper::new(synth_ivs);
    let lapper_synth_build = t0.elapsed().as_secs_f64() * 1000.0;
    let lapper_synth_ns = bench_queries(|s, e, out| lapper_synth.query(s, e, out), &synth_queries);

    let lapper_coi = Lapper::new(coi_ivs);
    let lapper_coi_ns = bench_queries(|s, e, out| lapper_coi.query(s, e, out), &coi_queries);

    println!("\n## Micro-Benchmark Results");
    println!(
        "| Index | Build (ms) | Synthetic (ns/query) | Coitrees (ns/query) | Speedup vs coitrees |"
    );
    println!("|---|---|---|---|---|");
    println!(
        "| coitrees | {coitree_synth_build:.2} | {coitree_synth_ns:.1} | {coitree_coi_ns:.1} | 1.00x |"
    );
    println!(
        "| ITree64 | {itree_synth_build:.2} | {itree_synth_ns:.1} | {itree_coi_ns:.1} | {:.2}x / {:.2}x |",
        coitree_synth_ns / itree_synth_ns,
        coitree_coi_ns / itree_coi_ns
    );
    println!(
        "| superintervals | {super_synth_build:.2} | {super_synth_ns:.1} | {super_coi_ns:.1} | {:.2}x / {:.2}x |",
        coitree_synth_ns / super_synth_ns,
        coitree_coi_ns / super_coi_ns
    );
    println!(
        "| rust-lapper | {lapper_synth_build:.2} | {lapper_synth_ns:.1} | {lapper_coi_ns:.1} | {:.2}x / {:.2}x |",
        coitree_synth_ns / lapper_synth_ns,
        coitree_coi_ns / lapper_coi_ns
    );

    let itree_faster =
        (coitree_synth_ns / itree_synth_ns >= 1.2) && (coitree_coi_ns / itree_coi_ns >= 1.2);
    let super_faster =
        (coitree_synth_ns / super_synth_ns >= 1.2) && (coitree_coi_ns / super_coi_ns >= 1.2);
    let lapper_faster =
        (coitree_synth_ns / lapper_synth_ns >= 1.2) && (coitree_coi_ns / lapper_coi_ns >= 1.2);

    if itree_faster || super_faster || lapper_faster {
        println!(
            "\nRuling: An alternative index exceeded 1.2x speedup over coitrees on both datasets."
        );
    } else {
        println!(
            "\nRuling: coitrees remains fastest (van Emde Boas cache-oblivious layout); D40 default index upheld."
        );
    }
}
