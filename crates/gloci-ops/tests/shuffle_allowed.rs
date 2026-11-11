//! `shuffle` placement (D152): exact counts, uniform draws and the constraints every
//! output must satisfy, whatever the include and exclude sets.
use gloci_core::{ChromId, ChromTable, Genome, Interval, OrderSource, VcfSpan};
use gloci_io::{
    Format, InvalidPolicy, IoError, ReaderOptions, RecordWriter, UsageError, open_reader,
};
use gloci_ops::rng::Stream;
use gloci_ops::shuffle::{Allowed, ShuffleOptions, shuffle};
use proptest::prelude::*;
use std::num::NonZeroUsize;

fn genome(sizes: &[(&str, u64)]) -> Genome {
    Genome::new(
        sizes
            .iter()
            .map(|(n, l)| (n.as_bytes().to_vec(), Some(*l)))
            .collect(),
    )
    .expect("genome")
}

fn table(g: &Genome) -> ChromTable {
    ChromTable::closed(g, OrderSource::Lexicographic)
}

fn iv(s: u64, e: u64) -> Interval {
    Interval::new(s, e).expect("interval")
}

#[test]
fn allowed_counts() {
    let g = genome(&[("a", 10), ("b", 30)]);
    let t = table(&g);
    let a = Allowed::build(&g, &t, None, &[]);
    assert_eq!(a.placements(5, None), 32);
    assert_eq!(a.placements(0, None), 42);
    assert_eq!(a.placements(31, None), 0);
    assert_eq!(a.placements(10, None), 1 + 21);
    assert_eq!(a.placements(5, Some(ChromId(0))), 6);
    assert_eq!(a.placements(5, Some(ChromId(1))), 26);
    assert_eq!(a.placements(11, Some(ChromId(0))), 0);
    // Regions on a chromosome that is not in the table (or past it) count for nothing.
    assert_eq!(a.placements(1, Some(ChromId(9))), 0);
    // A genome cut by include and exclude: b [0,30) minus [10,20) is [0,10) and [20,30).
    let inc = [(ChromId(1), iv(0, 30))];
    let exc = [(ChromId(1), iv(10, 20)), (ChromId(0), iv(0, 4))];
    let a = Allowed::build(&g, &t, Some(&inc), &exc);
    assert_eq!(a.placements(5, None), 12);
    assert_eq!(a.placements(5, Some(ChromId(0))), 0);
    // Bookended include regions merge: a record may straddle the seam.
    let inc = [(ChromId(1), iv(0, 10)), (ChromId(1), iv(10, 30))];
    let a = Allowed::build(&g, &t, Some(&inc), &[]);
    assert_eq!(a.placements(25, None), 6);
    // A region past the chromosome end is cut to it.
    let inc = [(ChromId(0), iv(5, 500))];
    let a = Allowed::build(&g, &t, Some(&inc), &[]);
    assert_eq!(a.placements(5, None), 1);
}

fn chi2(obs: &[f64], exp: &[f64]) -> f64 {
    obs.iter()
        .zip(exp)
        .map(|(o, e)| (o - e) * (o - e) / e)
        .sum()
}

#[test]
fn shuffle_uniform() {
    let n = 100_000u64;
    // Length-1 records over 100/200/700 bp: chromosomes in proportion to length, and every
    // base of the long one equally often.
    let g = genome(&[("c1", 100), ("c2", 200), ("c3", 700)]);
    let t = table(&g);
    let a = Allowed::build(&g, &t, None, &[]);
    let mut per = [0f64; 3];
    let mut bins = [0f64; 10];
    for j in 0..n {
        let mut s = Stream::new(7, j);
        let (c, p) = a.place(1, None, &mut s).expect("place");
        per[c.0 as usize] += 1.0;
        if c.0 == 2 {
            bins[(p / 70) as usize] += 1.0;
        }
    }
    let nf = n as f64;
    let x = chi2(&per, &[nf * 0.1, nf * 0.2, nf * 0.7]);
    assert!(x < 18.42, "chromosome counts chi2 {x}");
    let x = chi2(&bins, &[per[2] / 10.0; 10]);
    assert!(x < 33.72, "chr3 start bins chi2 {x}");

    // Allowed intervals of 10 and 30 bases, L = 5: 6 and 26 placements.
    let g = genome(&[("a", 10), ("b", 30)]);
    let t = table(&g);
    let a = Allowed::build(&g, &t, None, &[]);
    let mut per = [0f64; 2];
    for j in 0..n {
        let mut s = Stream::new(8, j);
        let (c, _) = a.place(5, None, &mut s).expect("place");
        per[c.0 as usize] += 1.0;
    }
    let x = chi2(&per, &[nf * 6.0 / 32.0, nf * 26.0 / 32.0]);
    assert!(x < 10.83, "interval counts chi2 {x}");

    // --chrom-first: each chromosome with room equally often, whatever its length.
    let g = genome(&[("c1", 100), ("c2", 200), ("c3", 700)]);
    let t = table(&g);
    let a = Allowed::build(&g, &t, None, &[]);
    let mut per = [0f64; 3];
    for j in 0..n {
        let mut s = Stream::new(9, j);
        let (c, _) = a.place_chrom_first(1, &mut s).expect("place");
        per[c.0 as usize] += 1.0;
    }
    let x = chi2(&per, &[nf / 3.0; 3]);
    assert!(x < 18.42, "chrom-first chi2 {x}");
    // A chromosome too short for the record is never picked.
    let mut s = Stream::new(1, 0);
    for _ in 0..200 {
        let (c, _) = a.place_chrom_first(150, &mut s).expect("place");
        assert_ne!(c.0, 0);
    }
    assert!(a.place_chrom_first(701, &mut s).is_none());
}

#[test]
fn placements_are_exact_and_cover_every_place() {
    // Every counted placement is drawable and none lies outside: enumerate by r.
    let g = genome(&[("a", 12), ("b", 7), ("c", 20)]);
    let t = table(&g);
    let exc = [(ChromId(2), iv(5, 8)), (ChromId(0), iv(3, 4))];
    let a = Allowed::build(&g, &t, None, &exc);
    for len in 0..8u64 {
        let total = a.placements(len, None);
        let mut seen = std::collections::BTreeSet::new();
        for j in 0..20_000u64 {
            let mut s = Stream::new(3, j);
            if let Some((c, p)) = a.place(len, None, &mut s) {
                seen.insert((c.0, p));
            }
        }
        assert_eq!(seen.len() as u128, total, "len {len}");
    }
}

struct Sink(Vec<String>);

impl RecordWriter for Sink {
    fn write_record(&mut self, raw: &[u8]) -> std::io::Result<()> {
        self.0.push(String::from_utf8(raw.to_vec()).unwrap());
        Ok(())
    }
    fn finish(self: Box<Self>) -> Result<(), IoError> {
        Ok(())
    }
}

fn opts(seed: u64) -> ShuffleOptions {
    ShuffleOptions {
        seed,
        include: None,
        exclude: None,
        same_chrom: false,
        chrom_first: false,
        no_overlapping: false,
        max_tries: 1000,
        skip_unplaceable: true,
    }
}

fn run(g: &Genome, bed: &std::path::Path, o: &ShuffleOptions) -> Result<Vec<String>, IoError> {
    let ro = ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: None,
        threads: NonZeroUsize::MIN,
        format: Some(Format::Bed),
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    };
    let r = open_reader(bed.to_str().unwrap(), &ro)?;
    let mut t = ChromTable::declared(g, OrderSource::Lexicographic);
    let mut w = Sink(Vec::new());
    shuffle(r, &mut t, g, o, NonZeroUsize::MIN, &mut w)?;
    Ok(w.0)
}

const LENS: [u64; 3] = [120, 200, 90];

fn regions() -> impl Strategy<Value = Vec<(usize, u64, u64)>> {
    prop::collection::vec(
        (0usize..3).prop_flat_map(|c| {
            (Just(c), 0..LENS[c]).prop_flat_map(move |(c, s)| (Just(c), Just(s), s + 1..=LENS[c]))
        }),
        0..5,
    )
}

fn covered(inc: &[(usize, u64, u64)], c: usize, s: u64, e: u64) -> bool {
    if s == e {
        return inc.iter().any(|&(ic, a, b)| ic == c && a <= s && s <= b);
    }
    (s..e).all(|x| inc.iter().any(|&(ic, a, b)| ic == c && a <= x && x < b))
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(96))]

    #[test]
    fn shuffle_constraints(
        include in regions(),
        exclude in regions(),
        feats in prop::collection::vec((0usize..3, 0u64..40, 0u64..30), 1..30),
        same_chrom in any::<bool>(),
        no_overlapping in any::<bool>(),
        seed in any::<u64>(),
    ) {
        let g = genome(&[("c0", LENS[0]), ("c1", LENS[1]), ("c2", LENS[2])]);
        let dir = tempfile::tempdir().unwrap();
        let mut bed = String::new();
        let mut inputs = Vec::new();
        for (i, &(c, s, l)) in feats.iter().enumerate() {
            let s = s.min(LENS[c]);
            let e = (s + l).min(LENS[c]);
            bed.push_str(&format!("c{c}\t{s}\t{e}\tr{i}\t{i}\t+\textra\n"));
            inputs.push((c, e - s));
        }
        let bed_path = dir.path().join("in.bed");
        std::fs::write(&bed_path, &bed).unwrap();
        let write = |name: &str, r: &[(usize, u64, u64)]| {
            let p = dir.path().join(name);
            let text: String = r.iter().map(|(c, s, e)| format!("c{c}\t{s}\t{e}\n")).collect();
            std::fs::write(&p, text).unwrap();
            p.to_str().unwrap().to_string()
        };
        let mut o = opts(seed);
        o.same_chrom = same_chrom;
        o.no_overlapping = no_overlapping;
        o.max_tries = 50;
        let all: Vec<(usize, u64, u64)> = (0..3).map(|c| (c, 0, LENS[c])).collect();
        let inc: &[(usize, u64, u64)] = if include.is_empty() { &all } else { &include };
        if !include.is_empty() {
            o.include = Some(write("inc.bed", &include));
        }
        if !exclude.is_empty() {
            o.exclude = Some(write("exc.bed", &exclude));
        }
        let out = run(&g, &bed_path, &o).unwrap();
        prop_assert!(out.len() <= feats.len());
        let mut placed: Vec<(String, u64, u64)> = Vec::new();
        for line in &out {
            let f: Vec<&str> = line.split('\t').collect();
            prop_assert_eq!(f.len(), 7);
            let i: usize = f[3][1..].parse().unwrap();
            let (c0, len) = inputs[i];
            let (s, e): (u64, u64) = (f[1].parse().unwrap(), f[2].parse().unwrap());
            prop_assert_eq!(e - s, len, "length kept: {}", line);
            let rest = format!("r{i}\t{i}\t+\textra");
            prop_assert_eq!(f[3..].join("\t"), rest);
            let c: usize = f[0][1..].parse().unwrap();
            if same_chrom {
                prop_assert_eq!(c, c0, "chromosome kept: {}", line);
            }
            prop_assert!(covered(inc, c, s, e), "outside the include set: {}", line);
            for &(xc, xs, xe) in &exclude {
                prop_assert!(!(xc == c && s < xe && xs < e), "overlaps exclude: {}", line);
            }
            for (pc, ps, pe) in &placed {
                prop_assert!(
                    !(no_overlapping && *pc == f[0] && s < *pe && *ps < e),
                    "overlaps another output: {}", line
                );
            }
            placed.push((f[0].to_string(), s, e));
        }
    }
}

#[test]
fn same_output_whatever_the_threads() {
    let g = genome(&[("c0", 1000), ("c1", 2000)]);
    let dir = tempfile::tempdir().unwrap();
    let text: String = (0..3000)
        .map(|i| format!("c{}\t{}\t{}\tr{i}\n", i % 2, i % 50, i % 50 + 1 + i % 17))
        .collect();
    let p = dir.path().join("in.bed");
    std::fs::write(&p, text).unwrap();
    let base = run(&g, &p, &opts(5)).unwrap();
    let ro = ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: None,
        threads: NonZeroUsize::MIN,
        format: Some(Format::Bed),
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    };
    for threads in [2, 8] {
        let r = open_reader(p.to_str().unwrap(), &ro).unwrap();
        let mut t = ChromTable::declared(&g, OrderSource::Lexicographic);
        let mut w = Sink(Vec::new());
        shuffle(
            r,
            &mut t,
            &g,
            &opts(5),
            NonZeroUsize::new(threads).unwrap(),
            &mut w,
        )
        .unwrap();
        assert_eq!(w.0, base, "-t {threads}");
    }
}

#[test]
fn missing_lengths_are_refused() {
    let g = Genome::new(vec![(b"c0".to_vec(), None)]).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("in.bed");
    std::fs::write(&p, "c0\t0\t5\n").unwrap();
    let e = run(&g, &p, &opts(1)).expect_err("error");
    assert!(matches!(e, IoError::Usage(UsageError(_))), "{e}");
}
