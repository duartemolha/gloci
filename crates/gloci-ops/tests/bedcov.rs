//! `bedcov` (D140): the totals of each region equal a naive per-base count, the output follows
//! the BED at every thread count, and an index is required.

use gloci_core::VcfSpan;
use gloci_io::{
    Accept, AlignOptions, InvalidPolicy, ReaderOptions, open_reader, write_alignment_index,
};
use gloci_ops::bedcov::{BedcovOptions, bedcov};
use gloci_ops::depth::CountRules;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

const LN1: u64 = 3000;
const LN2: u64 = 1500;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

#[derive(Clone)]
struct Read {
    rid: usize,
    /// 0-based.
    pos: u64,
    flag: u16,
    /// (op, length) over the reference and the query.
    ops: Vec<(char, u64)>,
}

impl Read {
    fn ref_len(&self) -> u64 {
        self.ops
            .iter()
            .filter(|(o, _)| matches!(o, 'M' | 'D' | 'N'))
            .map(|(_, n)| n)
            .sum()
    }
    fn query_len(&self) -> u64 {
        self.ops
            .iter()
            .filter(|(o, _)| matches!(o, 'M' | 'I'))
            .map(|(_, n)| n)
            .sum()
    }
    fn sam(&self, i: usize) -> String {
        let (rname, cigar, seq, qual) = if self.flag & 4 != 0 {
            (
                ["chr1", "chr2"][self.rid],
                "*".to_string(),
                "ACGT".to_string(),
                "IIII".to_string(),
            )
        } else {
            let cigar: String = self.ops.iter().map(|(o, n)| format!("{n}{o}")).collect();
            let q = self.query_len() as usize;
            (
                ["chr1", "chr2"][self.rid],
                cigar,
                "A".repeat(q),
                "I".repeat(q),
            )
        };
        format!(
            "r{i}\t{}\t{rname}\t{}\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}\n",
            self.flag,
            self.pos + 1
        )
    }
}

fn random_reads(seed: u64, n: usize) -> Vec<Read> {
    let mut rng = Rng(seed | 1);
    let mut reads = Vec::new();
    while reads.len() < n {
        let rid = usize::from(rng.below(3) == 0);
        let ln = [LN1, LN2][rid];
        let pos = rng.below(ln - 100);
        let mut ops = vec![('M', 1 + rng.below(30))];
        if rng.below(2) == 0 {
            ops.push((if rng.below(2) == 0 { 'D' } else { 'N' }, 1 + rng.below(40)));
            ops.push(('M', 1 + rng.below(30)));
        }
        if rng.below(4) == 0 {
            ops.push(('I', 1 + rng.below(5)));
            ops.push(('M', 1 + rng.below(10)));
        }
        let flag = match rng.below(12) {
            0 => 4,
            1 => 1024,
            2 => 16,
            _ => 0,
        };
        let r = Read {
            rid,
            pos,
            flag,
            ops,
        };
        if pos + r.ref_len() <= ln {
            reads.push(r);
        }
    }
    reads.sort_by_key(|r| (r.rid, r.pos));
    reads
}

fn write_bam(dir: &Path, name: &str, reads: &[Read]) -> PathBuf {
    let mut sam =
        format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{LN1}\n@SQ\tSN:chr2\tLN:{LN2}\n");
    for (i, r) in reads.iter().enumerate() {
        sam.push_str(&r.sam(i));
    }
    let bam = gloci_testkit::fixtures::sam_to_bam(&sam);
    let path = gloci_testkit::fixtures::write(dir, name, &bam);
    write_alignment_index(&path, &path, false, NonZeroUsize::MIN).unwrap();
    path
}

fn ro(exclude: u16, skip: bool) -> ReaderOptions {
    ReaderOptions {
        accept: Accept {
            alignments: Some(AlignOptions {
                exclude_flags: exclude,
                ..AlignOptions::default()
            }),
            ..Accept::default()
        },
        chrom_map: None,
        threads: NonZeroUsize::MIN,
        format: None,
        vcf_span: VcfSpan::Ref,
        invalid: if skip {
            InvalidPolicy::Skip
        } else {
            InvalidPolicy::Fail
        },
    }
}

fn run(
    bed: &Path,
    files: &[&Path],
    o: &BedcovOptions,
    t: usize,
) -> Result<String, gloci_io::IoError> {
    let r = ro(0x704, false);
    let reader = open_reader(
        bed.to_str().unwrap(),
        &ReaderOptions {
            accept: Accept::default(),
            ..r.clone()
        },
    )?;
    let files: Vec<String> = files.iter().map(|p| p.display().to_string()).collect();
    let mut out = Vec::new();
    bedcov(
        reader,
        &files,
        o,
        &r,
        NonZeroUsize::new(t).unwrap(),
        &mut out,
    )?;
    Ok(String::from_utf8(out).unwrap())
}

fn opts(split: bool, split_deletions: bool) -> BedcovOptions {
    BedcovOptions {
        rules: CountRules {
            deletions: !split_deletions,
            ref_skips: !split,
            min_baseq: 0,
        },
        bases_at_depth: Some(2),
        read_count: true,
        column_names: false,
    }
}

/// The sum, the bases at depth 2 and the reads of `iv` on `rid`, base by base.
fn naive(reads: &[Read], rid: usize, s: u64, e: u64, o: &BedcovOptions) -> (u64, u64, u64) {
    if s >= e {
        return (0, 0, 0);
    }
    let ln = [LN1, LN2][rid];
    let mut depth = vec![0u32; ln as usize];
    let mut n = 0;
    for r in reads.iter().filter(|r| r.rid == rid && r.flag & 0x704 == 0) {
        if r.pos < e && s < r.pos + r.ref_len() {
            n += 1;
        }
        let mut p = r.pos;
        for &(op, len) in &r.ops {
            match op {
                'M' => (p..p + len).for_each(|q| depth[q as usize] += 1),
                'D' if o.rules.deletions => (p..p + len).for_each(|q| depth[q as usize] += 1),
                'N' if o.rules.ref_skips => (p..p + len).for_each(|q| depth[q as usize] += 1),
                _ => {}
            }
            if matches!(op, 'M' | 'D' | 'N') {
                p += len;
            }
        }
    }
    let hi = e.min(ln);
    let (mut sum, mut at) = (0, 0);
    for q in s.min(hi)..hi {
        sum += u64::from(depth[q as usize]);
        at += u64::from(depth[q as usize] >= 2);
    }
    (sum, at, n)
}

fn bed_text(regions: &[(usize, u64, u64)]) -> String {
    regions
        .iter()
        .enumerate()
        .map(|(i, (rid, s, e))| format!("{}\t{s}\t{e}\tr{i}\n", ["chr1", "chr2"][*rid]))
        .collect()
}

#[test]
fn bedcov_equals_naive() {
    let dir = tempfile::tempdir().unwrap();
    let reads = random_reads(7, 500);
    let bam = write_bam(dir.path(), "a.bam", &reads);
    let mut rng = Rng(99);
    let mut regions = Vec::new();
    for _ in 0..300 {
        let rid = usize::from(rng.below(3) == 0);
        let ln = [LN1, LN2][rid];
        let s = rng.below(ln + 20);
        let e = match rng.below(5) {
            0 => s,
            1 => s + 1,
            2 => s + rng.below(400),
            3 => s + rng.below(40),
            _ => ln + rng.below(30),
        }
        .max(s);
        regions.push((rid, s, e));
    }
    // Bookends of reads: regions that end where a read starts and start where it ends.
    for r in reads.iter().take(20) {
        regions.push((r.rid, r.pos.saturating_sub(5), r.pos));
        regions.push((r.rid, r.pos + r.ref_len(), r.pos + r.ref_len() + 5));
    }
    let bed = gloci_testkit::fixtures::write(dir.path(), "r.bed", bed_text(&regions).as_bytes());
    for (split, split_deletions) in [(false, false), (true, false), (false, true), (true, true)] {
        let o = opts(split, split_deletions);
        let got = run(&bed, &[&bam], &o, 2).unwrap();
        let lines: Vec<&str> = got.lines().collect();
        assert_eq!(lines.len(), regions.len());
        for (i, (line, (rid, s, e))) in lines.iter().zip(&regions).enumerate() {
            let (sum, at, n) = naive(&reads, *rid, *s, *e, &o);
            let want = format!(
                "{}\t{s}\t{e}\tr{i}\t{sum}\t{at}\t{n}",
                ["chr1", "chr2"][*rid]
            );
            assert_eq!(
                *line, want,
                "split {split} split-deletions {split_deletions}"
            );
        }
    }
}

#[test]
fn bedcov_output_order_at_every_thread_count() {
    let dir = tempfile::tempdir().unwrap();
    let reads = random_reads(11, 400);
    let bam = write_bam(dir.path(), "a.bam", &reads);
    let mut rng = Rng(5);
    let regions: Vec<(usize, u64, u64)> = (0..10_000)
        .map(|_| {
            let rid = usize::from(rng.below(3) == 0);
            let ln = [LN1, LN2][rid];
            let s = rng.below(ln);
            (rid, s, s + rng.below(60))
        })
        .collect();
    let bed = gloci_testkit::fixtures::write(dir.path(), "r.bed", bed_text(&regions).as_bytes());
    let o = opts(false, false);
    let one = run(&bed, &[&bam], &o, 1).unwrap();
    for (i, line) in one.lines().enumerate() {
        let (rid, s, e) = regions[i];
        assert!(
            line.starts_with(&format!("{}\t{s}\t{e}\tr{i}\t", ["chr1", "chr2"][rid])),
            "line {i}: {line}"
        );
    }
    assert_eq!(one.lines().count(), 10_000);
    for t in [2, 8] {
        assert_eq!(run(&bed, &[&bam], &o, t).unwrap(), one, "-t {t}");
    }
}

#[test]
fn bedcov_needs_an_index() {
    let dir = tempfile::tempdir().unwrap();
    let reads = random_reads(3, 20);
    let bam = write_bam(dir.path(), "depth2.bam", &reads);
    std::fs::remove_file(format!("{}.bai", bam.display())).unwrap();
    let bed = gloci_testkit::fixtures::write(dir.path(), "r.bed", b"chr1\t0\t10\n");
    let e = run(&bed, &[&bam], &opts(false, false), 1).unwrap_err();
    assert!(matches!(e, gloci_io::IoError::Usage(_)), "{e:?}");
    let msg = e.to_string();
    assert!(
        msg.contains(&format!("gloci index {}", bam.display())),
        "{msg}"
    );
}
