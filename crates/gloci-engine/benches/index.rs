//! Index micro-benchmark (criterion): `rust-lapper`, which unsorted mode
//! uses (D61), against the real `coitrees` (its predecessor, D40) and
//! `superintervals` crates, on two synthetic sets of 100,000 intervals over a
//! 50 Mb chromosome: uniform positions (lengths 50–5,000), and clustered
//! around 500 centres (lengths 50–20,000, deep overlaps). 10,000 queries
//! each, in random order (unsorted mode queries A in input order).
//!
//! Run: `cargo bench -p gloci-engine --bench index`. D40's rule: a crate
//! replaces the incumbent only if it is at least 1.2x faster on both sets.

use coitrees::{COITree, GenericInterval, Interval as CoiInterval, IntervalTree};
use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use rand::{RngExt, SeedableRng};
use rand_chacha::ChaCha8Rng;
use rust_lapper::{Interval as LapInterval, Lapper};
use std::hint::black_box;
use std::time::Duration;
use superintervals::IntervalMap;

const CHROM_LEN: u64 = 50_000_000;
const N_INTERVALS: usize = 100_000;
const N_QUERIES: usize = 10_000;

/// Half-open `[start, end)` intervals with their index, and the queries.
struct Set {
    name: &'static str,
    ivs: Vec<(u64, u64, u32)>,
    queries: Vec<(u64, u64)>,
}

fn uniform(seed: u64) -> Set {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let ivs = (0..N_INTERVALS)
        .map(|i| {
            let len = rng.random_range(50..=5_000);
            let s = rng.random_range(0..=CHROM_LEN - len);
            (s, s + len, i as u32)
        })
        .collect();
    let queries = (0..N_QUERIES)
        .map(|_| {
            let len = rng.random_range(100..=1_000);
            let s = rng.random_range(0..=CHROM_LEN - len);
            (s, s + len)
        })
        .collect();
    Set {
        name: "uniform",
        ivs,
        queries,
    }
}

fn clustered(seed: u64) -> Set {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let centres: Vec<u64> = (0..500).map(|_| rng.random_range(0..=CHROM_LEN)).collect();
    let around = |rng: &mut ChaCha8Rng, lo: u64, hi: u64| {
        let c = centres[rng.random_range(0..centres.len())];
        let len = rng.random_range(lo..=hi);
        let s = c.saturating_sub(len / 2).min(CHROM_LEN - len);
        (s, s + len)
    };
    let ivs = (0..N_INTERVALS)
        .map(|i| {
            let (s, e) = around(&mut rng, 50, 20_000);
            (s, e, i as u32)
        })
        .collect();
    let queries = (0..N_QUERIES)
        .map(|_| around(&mut rng, 100, 2_000))
        .collect();
    Set {
        name: "clustered",
        ivs,
        queries,
    }
}

/// The three indexes over one set. coitrees and superintervals take closed
/// `i32` coordinates, so `[s, e)` is stored as `[s, e - 1]` and a query
/// `[qs, qe)` asks `[qs, qe - 1]`: every index reports the same hits.
struct Indexes {
    coitrees: COITree<u32, u32>,
    superintervals: IntervalMap<u32>,
    lapper: Lapper<u64, u32>,
}

fn coi_nodes(set: &Set) -> Vec<CoiInterval<u32>> {
    set.ivs
        .iter()
        .map(|&(s, e, v)| CoiInterval::new(s as i32, e as i32 - 1, v))
        .collect()
}

fn superintervals(set: &Set) -> IntervalMap<u32> {
    let mut m = IntervalMap::new();
    for &(s, e, v) in &set.ivs {
        m.add(s as i32, e as i32 - 1, v);
    }
    m.build();
    m
}

fn lapper(set: &Set) -> Lapper<u64, u32> {
    Lapper::new(
        set.ivs
            .iter()
            .map(|&(start, stop, val)| LapInterval { start, stop, val })
            .collect(),
    )
}

fn build(set: &Set) -> Indexes {
    Indexes {
        coitrees: COITree::new(&coi_nodes(set)),
        superintervals: superintervals(set),
        lapper: lapper(set),
    }
}

/// Hits of every query, summed.
fn run_queries(name: &str, ix: &Indexes, queries: &[(u64, u64)], buf: &mut Vec<u32>) -> usize {
    let mut total = 0;
    for &(s, e) in queries {
        buf.clear();
        match name {
            "coitrees" => ix
                .coitrees
                .query(s as i32, e as i32 - 1, |n| buf.push(*n.metadata())),
            "superintervals" => ix.superintervals.search_values(s as i32, e as i32 - 1, buf),
            "rust-lapper" => buf.extend(ix.lapper.find(s, e).map(|iv| iv.val)),
            _ => unreachable!(),
        }
        total += buf.len();
    }
    total
}

const NAMES: [&str; 3] = ["coitrees", "superintervals", "rust-lapper"];

fn bench_build(c: &mut Criterion) {
    for set in [uniform(42), clustered(43)] {
        let mut g = c.benchmark_group(format!("build/{}", set.name));
        g.throughput(Throughput::Elements(set.ivs.len() as u64));
        let nodes = coi_nodes(&set);
        g.bench_function("coitrees", |b| {
            b.iter(|| COITree::<u32, u32>::new(black_box(&nodes)))
        });
        g.bench_function("superintervals", |b| {
            b.iter(|| superintervals(black_box(&set)))
        });
        g.bench_function("rust-lapper", |b| b.iter(|| lapper(black_box(&set))));
        g.finish();
    }
}

fn bench_query(c: &mut Criterion) {
    for set in [uniform(42), clustered(43)] {
        let ix = build(&set);
        let mut buf = Vec::with_capacity(256);
        // Every index finds the same hits.
        let hits: Vec<usize> = NAMES
            .iter()
            .map(|n| run_queries(n, &ix, &set.queries, &mut buf))
            .collect();
        assert!(
            hits.iter().all(|&h| h == hits[0]),
            "{}: hits differ: {hits:?}",
            set.name
        );
        let mut g = c.benchmark_group(format!("query/{}", set.name));
        g.throughput(Throughput::Elements(set.queries.len() as u64));
        for name in NAMES {
            g.bench_with_input(BenchmarkId::from_parameter(name), &name, |b, name| {
                b.iter(|| run_queries(name, &ix, black_box(&set.queries), &mut buf))
            });
        }
        g.finish();
    }
}

fn config() -> Criterion {
    Criterion::default()
        .sample_size(30)
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(3))
}

criterion_group! {
    name = benches;
    config = config();
    targets = bench_build, bench_query
}
criterion_main!(benches);
