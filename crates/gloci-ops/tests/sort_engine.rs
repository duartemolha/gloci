//! The sort engine: in-memory and external merge sort with bounded memory
//! and an order that does not depend on the thread count (spec §7.4, D13,
//! D24, Review Focus 5).

use gloci_core::{ChromTable, VcfSpan};
use gloci_io::{
    Format, InvalidPolicy, IoError, ReaderOptions, RecordBatch, RecordReader, open_reader_from,
};
use gloci_ops::sort::{Dedup, SortConfig, SortStats, Sorter, Tie, merge_sorted};
use gloci_testkit::naive;
use gloci_testkit::strategies::{self, TestRec};
use proptest::prelude::*;
use std::io::{self, Cursor};
use std::num::NonZeroUsize;
use std::path::Path;

/// Records read per batch: small and odd, so batches straddle chunk
/// boundaries.
const BATCH: usize = 7;

/// Chromosomes `gen_records` draws from (a superset of `strategies::CHROMS`),
/// in lexicographic byte order.
const GEN_CHROMS: &[&str] = &[
    "chr1",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr1_KI270706v1_random",
    "chr2",
    "chr20",
    "chr21",
    "chr22",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chrM",
    "chrX",
    "chrY",
];

/// Rank of a chromosome in lexicographic byte order.
fn lexicographic_rank(chrom: &str) -> u32 {
    GEN_CHROMS
        .binary_search(&chrom)
        .unwrap_or_else(|_| panic!("unknown chromosome {chrom}")) as u32
}

#[test]
fn gen_chroms_are_lexicographic_and_cover_the_strategies() {
    assert!(
        GEN_CHROMS
            .windows(2)
            .all(|w| w[0].as_bytes() < w[1].as_bytes())
    );
    for c in strategies::CHROMS {
        lexicographic_rank(c);
    }
}

fn cfg(max_mem: u64, tie: Tie, dedup: Dedup) -> SortConfig {
    SortConfig {
        max_mem,
        tmp_dir: std::env::temp_dir(),
        threads: NonZeroUsize::new(4).unwrap(),
        tie,
        dedup,
    }
}

/// A spilling configuration with `n` threads.
fn threads(n: usize) -> SortConfig {
    SortConfig {
        threads: NonZeroUsize::new(n).unwrap(),
        ..cfg(256 << 10, Tie::InputOrder, Dedup::All)
    }
}

fn reader(name: &str, text: impl Into<Vec<u8>>) -> Box<dyn RecordReader> {
    let opts = ReaderOptions {
        threads: NonZeroUsize::MIN,
        format: Some(Format::Bed),
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    };
    open_reader_from(name, Box::new(Cursor::new(text.into())), &opts).unwrap()
}

fn bed_text(recs: &[TestRec]) -> String {
    recs.iter().map(|r| r.bed_line() + "\n").collect()
}

fn lines_of(out: &[u8]) -> Vec<String> {
    String::from_utf8(out.to_vec())
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect()
}

/// Read every batch of `r` into a new `Sorter`.
fn fill(r: &mut dyn RecordReader, chroms: &mut ChromTable, cfg: SortConfig) -> Sorter {
    let mut s = Sorter::new(cfg).unwrap();
    let mut batch = RecordBatch::default();
    loop {
        batch.clear();
        let more = r.read_batch(chroms, &mut batch, BATCH).unwrap();
        s.push(chroms, &batch).unwrap();
        if !more {
            return s;
        }
    }
}

/// Sort BED text; output lines joined with `\n`.
fn sort_text(text: impl Into<Vec<u8>>, cfg: SortConfig) -> (Vec<u8>, SortStats) {
    let mut r = reader("in.bed", text);
    let mut chroms = ChromTable::lexicographic();
    let s = fill(r.as_mut(), &mut chroms, cfg);
    let mut out = Vec::new();
    let stats = s
        .finish(&chroms, &mut |raw: &[u8]| -> io::Result<()> {
            out.extend_from_slice(raw);
            out.push(b'\n');
            Ok(())
        })
        .unwrap();
    assert_eq!(
        stats.records,
        out.iter().filter(|&&b| b == b'\n').count() as u64
    );
    (out, stats)
}

fn sort_lines(recs: &[TestRec], cfg: SortConfig) -> Vec<String> {
    lines_of(&sort_text(bed_text(recs), cfg).0)
}

fn sort_stats(recs: &[TestRec], cfg: SortConfig) -> SortStats {
    sort_text(bed_text(recs), cfg).1
}

fn sort_bytes(recs: &[TestRec], cfg: SortConfig) -> Vec<u8> {
    sort_text(bed_text(recs), cfg).0
}

fn sort_lines_raw(lines: &[&str], cfg: SortConfig) -> Vec<String> {
    let text: String = lines.iter().map(|l| format!("{l}\n")).collect();
    lines_of(&sort_text(text, cfg).0)
}

/// `n` pseudo-random records (deterministic). With `ties`, coordinates come
/// from a small range so many records share `(chrom, start, end)`.
fn generate(n: usize, ties: bool) -> Vec<TestRec> {
    let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut next = move || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    };
    (0..n)
        .map(|i| {
            let chrom = GEN_CHROMS[(next() % GEN_CHROMS.len() as u64) as usize];
            let (start, len) = if ties {
                (next() % 20, next() % 3)
            } else {
                (next() % 1_000_000, next() % 1_000)
            };
            TestRec {
                chrom: chrom.to_string(),
                start,
                end: start + len,
                strand: ['+', '-', '.'][(next() % 3) as usize],
                name: format!("g{i}"),
            }
        })
        .collect()
}

fn gen_records(n: usize) -> Vec<TestRec> {
    generate(n, false)
}

fn naive_lines(recs: &[TestRec]) -> Vec<String> {
    naive::sort(recs.to_vec(), lexicographic_rank)
        .iter()
        .map(TestRec::bed_line)
        .collect()
}

fn merge_with(texts: &[&str], tie: Tie) -> Result<Vec<String>, IoError> {
    let mut readers: Vec<Box<dyn RecordReader>> = texts
        .iter()
        .enumerate()
        .map(|(i, t)| reader(&format!("in{i}.bed"), t.as_bytes()))
        .collect();
    let mut chroms = ChromTable::lexicographic();
    let mut out = Vec::new();
    let n = merge_sorted(
        &mut readers,
        &mut chroms,
        tie,
        &std::env::temp_dir(),
        &mut |raw: &[u8]| -> io::Result<()> {
            out.extend_from_slice(raw);
            out.push(b'\n');
            Ok(())
        },
    )?;
    let lines = lines_of(&out);
    assert_eq!(n, lines.len() as u64);
    Ok(lines)
}

fn merge(texts: &[&str]) -> Result<Vec<String>, IoError> {
    merge_with(texts, Tie::InputOrder)
}

fn entries(dir: &Path) -> usize {
    std::fs::read_dir(dir).unwrap().count()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn matches_naive_in_memory_and_spilled(recs in strategies::bed_records()) {
        let want = naive::sort(recs.clone(), lexicographic_rank).iter().map(TestRec::bed_line).collect::<Vec<_>>();
        prop_assert_eq!(sort_lines(&recs, cfg(64 << 20, Tie::InputOrder, Dedup::All)), want.clone());
        prop_assert_eq!(sort_lines(&recs, cfg(4 << 10, Tie::InputOrder, Dedup::All)), want.clone());   // forces many spills
        // A few records per chunk (the fan-in is 256 / longest line, 5 to 17,
        // so merges cascade), and one record per chunk (fan-in 2).
        let (out, stats) = sort_text(bed_text(&recs), cfg(256, Tie::InputOrder, Dedup::All));
        prop_assert_eq!(lines_of(&out), want.clone());
        // The shortest line is 15 bytes, so at most 4 records fit in 256 bytes.
        prop_assert!(recs.len() <= 4 || stats.spills >= 2);
        prop_assert_eq!(sort_lines(&recs, cfg(0, Tie::InputOrder, Dedup::All)), want);
    }

    #[test]
    fn unique_and_duplicates_like_uniq(recs in strategies::bed_records()) {
        let full = naive::sort_full_line(recs.clone(), lexicographic_rank);
        prop_assert_eq!(sort_lines(&recs, cfg(4 << 10, Tie::FullLine, Dedup::Unique)), naive::uniq(&full));
        prop_assert_eq!(sort_lines(&recs, cfg(4 << 10, Tie::FullLine, Dedup::Duplicates)), naive::uniq_d(&full));
        // Across chunk and run boundaries; --unique/--duplicates imply the
        // full-line order (D13).
        prop_assert_eq!(sort_lines(&recs, cfg(256, Tie::InputOrder, Dedup::Unique)), naive::uniq(&full));
        prop_assert_eq!(sort_lines(&recs, cfg(256, Tie::InputOrder, Dedup::Duplicates)), naive::uniq_d(&full));
        prop_assert_eq!(sort_lines(&recs, cfg(64 << 20, Tie::InputOrder, Dedup::Duplicates)), naive::uniq_d(&full));
    }

    #[test]
    fn full_line_ties(recs in strategies::bed_records()) {
        let full = naive::sort_full_line(recs.clone(), lexicographic_rank);
        prop_assert_eq!(sort_lines(&recs, cfg(64 << 20, Tie::FullLine, Dedup::All)), full.clone());
        prop_assert_eq!(sort_lines(&recs, cfg(256, Tie::FullLine, Dedup::All)), full);
    }

    #[test]
    fn merge_sorted_matches_naive(recs in strategies::bed_records(), cut1 in 0usize..40, cut2 in 0usize..40) {
        // Three sorted inputs; ties go to the earlier input.
        let (a, b) = (cut1.min(recs.len()), cut2.min(recs.len()));
        let (lo, hi) = (a.min(b), a.max(b));
        let parts = [&recs[..lo], &recs[lo..hi], &recs[hi..]];
        let texts: Vec<String> = parts
            .iter()
            .map(|p| bed_text(&naive::sort(p.to_vec(), lexicographic_rank)))
            .collect();
        let refs: Vec<&str> = texts.iter().map(String::as_str).collect();
        let want: Vec<String> = naive::sort(
            parts.iter().flat_map(|p| naive::sort(p.to_vec(), lexicographic_rank)).collect(),
            lexicographic_rank,
        )
        .iter()
        .map(TestRec::bed_line)
        .collect();
        prop_assert_eq!(merge(&refs).unwrap(), want);

        let full: Vec<String> = parts.iter().map(|p| {
            naive::sort_full_line(p.to_vec(), lexicographic_rank).iter().map(|l| format!("{l}\n")).collect()
        }).collect();
        let refs: Vec<&str> = full.iter().map(String::as_str).collect();
        prop_assert_eq!(merge_with(&refs, Tie::FullLine).unwrap(), naive::sort_full_line(recs.clone(), lexicographic_rank));
    }
}

#[test]
fn many_runs_cascade() {
    let recs = gen_records(200_000);
    let (out, s) = sort_text(bed_text(&recs), cfg(64 << 10, Tie::InputOrder, Dedup::All));
    assert!(s.spills > 128, "only {} spills", s.spills);
    assert_eq!(s.records, recs.len() as u64);
    assert_eq!(lines_of(&out), naive_lines(&recs));
    assert_eq!(
        sort_stats(&recs[..1000], cfg(64 << 20, Tie::InputOrder, Dedup::All)),
        SortStats {
            records: 1000,
            spills: 0
        }
    );
}

#[test]
fn thread_count_does_not_change_output() {
    let input = generate(50_000, true);
    let want = sort_bytes(&input, threads(1));
    assert_eq!(want, sort_bytes(&input, threads(8)));
    assert_eq!(lines_of(&want), naive_lines(&input));
    // In memory too.
    let in_memory = |n| SortConfig {
        max_mem: 64 << 20,
        ..threads(n)
    };
    assert_eq!(sort_bytes(&input, in_memory(1)), want);
    assert_eq!(sort_bytes(&input, in_memory(8)), want);
    let full = |n| SortConfig {
        tie: Tie::FullLine,
        ..threads(n)
    };
    assert_eq!(sort_bytes(&input, full(1)), sort_bytes(&input, full(8)));
}

#[test]
fn record_larger_than_budget() {
    // Review Focus 5
    let big = format!("chr1\t5\t6\t{}", "x".repeat(3 << 20));
    assert_eq!(
        sort_lines_raw(
            &[&big, "chr1\t1\t2"],
            cfg(1 << 20, Tie::InputOrder, Dedup::All)
        ),
        vec!["chr1\t1\t2".to_string(), big]
    );
    // Several, each in a chunk of its own, between small records.
    let bigs: Vec<String> = (0..4)
        .map(|i| format!("chr1\t{}\t9\t{}", 4 - i, "y".repeat((2 << 20) + i)))
        .collect();
    let mut input: Vec<&str> = vec!["chr1\t7\t8", "chr1\t0\t1"];
    input.extend(bigs.iter().map(String::as_str));
    input.push("chr1\t3\t4");
    let mut want: Vec<String> = vec!["chr1\t0\t1".into()];
    want.extend(bigs.iter().rev().cloned());
    want.insert(3, "chr1\t3\t4".into());
    want.push("chr1\t7\t8".into());
    assert_eq!(
        sort_lines_raw(&input, cfg(1 << 20, Tie::InputOrder, Dedup::All)),
        want
    );
}

#[test]
fn spill_dir_removed_after_finish_and_on_drop() {
    let tmp = tempfile::tempdir().unwrap();
    let small = SortConfig {
        tmp_dir: tmp.path().to_path_buf(),
        ..cfg(4 << 10, Tie::InputOrder, Dedup::All)
    };
    let recs = gen_records(5_000);

    // finish()
    let mut r = reader("in.bed", bed_text(&recs));
    let mut chroms = ChromTable::lexicographic();
    let s = fill(r.as_mut(), &mut chroms, small.clone());
    assert_eq!(entries(tmp.path()), 1, "one spill directory while sorting");
    let mut n = 0u64;
    let stats = s
        .finish(&chroms, &mut |_: &[u8]| -> io::Result<()> {
            n += 1;
            Ok(())
        })
        .unwrap();
    assert!(stats.spills > 1);
    assert_eq!(n, recs.len() as u64);
    assert_eq!(
        entries(tmp.path()),
        0,
        "spill directory removed after finish"
    );

    // dropped mid-push
    let mut r = reader("in.bed", bed_text(&recs));
    let mut chroms = ChromTable::lexicographic();
    let mut s = Sorter::new(small.clone()).unwrap();
    let mut batch = RecordBatch::default();
    r.read_batch(&mut chroms, &mut batch, 2_000).unwrap();
    s.push(&chroms, &batch).unwrap();
    let spill_dir = std::fs::read_dir(tmp.path())
        .unwrap()
        .next()
        .unwrap()
        .unwrap()
        .path();
    assert!(entries(&spill_dir) > 0, "runs were spilled");
    drop(s);
    assert_eq!(entries(tmp.path()), 0, "spill directory removed on drop");

    // finish() failing on the output
    let mut r = reader("in.bed", bed_text(&recs));
    let mut chroms = ChromTable::lexicographic();
    let s = fill(r.as_mut(), &mut chroms, small);
    let err = s
        .finish(&chroms, &mut |_: &[u8]| -> io::Result<()> {
            Err(io::ErrorKind::BrokenPipe.into())
        })
        .unwrap_err();
    assert!(
        matches!(&err, IoError::Write { source, .. } if source.kind() == io::ErrorKind::BrokenPipe)
    );
    assert_eq!(entries(tmp.path()), 0, "spill directory removed on error");
}

#[test]
fn unwritable_tmp_dir_fails_up_front() {
    let tmp = tempfile::tempdir().unwrap();
    let missing = SortConfig {
        tmp_dir: tmp.path().join("missing"),
        ..cfg(4 << 10, Tie::InputOrder, Dedup::All)
    };
    assert!(matches!(Sorter::new(missing), Err(IoError::Write { .. })));
}

#[test]
fn output_errors_keep_their_kind() {
    // In memory: the error reaches the caller with its io::ErrorKind intact
    // (so `| head` can be recognised as a broken pipe).
    let mut r = reader("in.bed", "chr1\t1\t2\nchr1\t0\t1\n");
    let mut chroms = ChromTable::lexicographic();
    let s = fill(
        r.as_mut(),
        &mut chroms,
        cfg(1 << 20, Tie::InputOrder, Dedup::All),
    );
    let err = s
        .finish(&chroms, &mut |_: &[u8]| -> io::Result<()> {
            Err(io::ErrorKind::BrokenPipe.into())
        })
        .unwrap_err();
    let IoError::Write { source, .. } = &err else {
        panic!("{err:?}")
    };
    assert_eq!(source.kind(), io::ErrorKind::BrokenPipe);

    let mut readers = vec![reader("a.bed", "chr1\t1\t2\n")];
    let err = merge_sorted(
        &mut readers,
        &mut chroms,
        Tie::InputOrder,
        &std::env::temp_dir(),
        &mut |_: &[u8]| -> io::Result<()> { Err(io::ErrorKind::BrokenPipe.into()) },
    )
    .unwrap_err();
    assert!(
        matches!(&err, IoError::Write { source, .. } if source.kind() == io::ErrorKind::BrokenPipe)
    );
}

#[test]
fn chromosomes_first_seen_after_a_spill() {
    // Ranks are recomputed at each chunk sort and at the merge: chr10 and
    // chr1 only appear after chunks holding chr2 were spilled.
    let mut lines = Vec::new();
    for c in ["chr2", "chr10", "chr1"] {
        for i in (0..300).rev() {
            lines.push(format!("{c}\t{i}\t{}", i + 1));
        }
    }
    let refs: Vec<&str> = lines.iter().map(String::as_str).collect();
    let mut want = lines.clone();
    want.sort_by_key(|l| {
        let mut f = l.split('\t');
        let c = f.next().unwrap().to_string();
        (c, f.next().unwrap().parse::<u64>().unwrap())
    });
    assert_eq!(
        sort_lines_raw(&refs, cfg(1 << 10, Tie::InputOrder, Dedup::All)),
        want
    );
}

#[test]
fn empty_input() {
    let (out, stats) = sort_text("", cfg(1 << 20, Tie::InputOrder, Dedup::All));
    assert!(out.is_empty());
    assert_eq!(stats, SortStats::default());
    assert_eq!(merge(&["", ""]).unwrap(), Vec::<String>::new());
}

#[test]
fn merge_sorted_rejects_unsorted_input() {
    let sorted_a = "chr1\t1\t2\nchr1\t5\t6\n";
    let unsorted_b = "chr1\t9\t10\nchr1\t3\t4\n";
    assert!(matches!(
        merge(&[sorted_a, unsorted_b]),
        Err(IoError::Order(_))
    ));
    let Err(IoError::Order(e)) = merge(&[sorted_a, unsorted_b]) else {
        unreachable!()
    };
    assert_eq!(e.path, "in1.bed");
    assert_eq!(
        merge(&[sorted_a, "chr1\t3\t4\nchr2\t0\t1\n"]).unwrap(),
        ["chr1\t1\t2", "chr1\t3\t4", "chr1\t5\t6", "chr2\t0\t1"]
    );
}

/// M2: more than 128 inputs spill to the caller's tmp dir, not the system
/// temp dir, and leave nothing there.
#[test]
fn merge_sorted_rounds_use_the_given_tmp_dir() {
    let make = || -> Vec<Box<dyn RecordReader>> {
        (0..130)
            .map(|i| {
                reader(
                    &format!("in{i}.bed"),
                    format!("chr1\t{i}\t{}\n", i + 1).as_bytes(),
                )
            })
            .collect()
    };
    let missing = std::path::Path::new("/nonexistent-gloci-tmp/dir");
    let err = merge_sorted(
        &mut make(),
        &mut ChromTable::lexicographic(),
        Tie::InputOrder,
        missing,
        &mut |_: &[u8]| -> io::Result<()> { Ok(()) },
    )
    .unwrap_err();
    assert!(
        err.to_string().contains("/nonexistent-gloci-tmp/dir"),
        "{err}"
    );
    let tmp = tempfile::tempdir().unwrap();
    let mut n = 0;
    merge_sorted(
        &mut make(),
        &mut ChromTable::lexicographic(),
        Tie::InputOrder,
        tmp.path(),
        &mut |_: &[u8]| -> io::Result<()> {
            n += 1;
            Ok(())
        },
    )
    .unwrap();
    assert_eq!(n, 130);
    assert_eq!(std::fs::read_dir(tmp.path()).unwrap().count(), 0);
}
