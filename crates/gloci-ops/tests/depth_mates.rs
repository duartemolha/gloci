//! `depth --skip-mate-overlap` (D139): what the mate table clips, how large it grows, and that
//! per-chromosome runs equal one stream.

use gloci_core::{Interval, VcfSpan};
use gloci_io::{Accept, AlignOptions, InvalidPolicy, ReaderOptions};
use gloci_ops::depth::{AllPositions, CountRules, DepthOptions, DepthSummary, depth};
use std::collections::BTreeMap;
use std::num::NonZeroUsize;
use std::path::PathBuf;

fn mates() -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/samtools/mates.bam")
        .display()
        .to_string()
}

/// samtools depth's default exclusion: UNMAP, SECONDARY, QCFAIL, DUP.
const DEFAULT_EXCLUDE: u16 = 0x704;

fn ro() -> ReaderOptions {
    ReaderOptions {
        accept: Accept {
            alignments: Some(AlignOptions {
                exclude_flags: DEFAULT_EXCLUDE,
                ..AlignOptions::default()
            }),
            ..Accept::default()
        },
        chrom_map: None,
        threads: NonZeroUsize::MIN,
        format: None,
        vcf_span: VcfSpan::Ref,
        invalid: InvalidPolicy::Fail,
    }
}

fn opts(skip: bool, deletions: bool, min_baseq: u8) -> DepthOptions {
    DepthOptions {
        rules: CountRules {
            deletions,
            ref_skips: false,
            min_baseq,
        },
        min_read_length: 0,
        all: AllPositions::Covered,
        regions: None,
        skip_mate_overlap: skip,
        column_names: false,
    }
}

fn run(paths: &[&str], o: &DepthOptions, t: usize) -> (String, DepthSummary) {
    let paths: Vec<String> = paths.iter().map(|p| p.to_string()).collect();
    let mut out = Vec::new();
    let s = depth(&paths, o, &ro(), NonZeroUsize::new(t).unwrap(), &mut out).unwrap();
    (String::from_utf8(out).unwrap(), s)
}

/// `(chrom, 1-based position)` to the depth of the first input.
fn table(out: &str) -> BTreeMap<(String, u64), u32> {
    out.lines()
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (
                (f[0].to_string(), f[1].parse().unwrap()),
                f[2].parse().unwrap(),
            )
        })
        .collect()
}

/// Every position of `chrom` from `a` to `b` (1-based, inclusive) has depth `want`.
#[track_caller]
fn expect(t: &BTreeMap<(String, u64), u32>, chrom: &str, a: u64, b: u64, want: u32) {
    for p in a..=b {
        assert_eq!(
            t.get(&(chrom.to_string(), p)).copied(),
            Some(want),
            "{chrom}:{p} wants {want}"
        );
    }
}

#[test]
fn mate_overlap_table() {
    let m = mates();
    let plain = table(&run(&[&m], &opts(true, false, 0), 1).0);
    // a: partial overlap, once over [100, 250).
    expect(&plain, "chr1", 101, 250, 1);
    // b: both mates at one position, the longer one second: once over the longer read.
    expect(&plain, "chr1", 301, 380, 1);
    // c: the second mate inside the first: once over the first.
    expect(&plain, "chr1", 401, 600, 1);
    // d: no overlap, both counted, and nothing is held afterwards (see mate_table_holds_*).
    expect(&plain, "chr1", 701, 750, 1);
    expect(&plain, "chr1", 901, 950, 1);
    // f, with the mate position unknown, counts once.
    expect(&plain, "chr1", 1201, 1250, 1);
    // A mate on another reference is not paired with anything.
    expect(&plain, "chr1", 1801, 1840, 1);
    expect(&plain, "chr2", 5, 44, 1);
    expect(&plain, "chr2", 101, 170, 1);

    // e under -J: the deletion before the clip is not counted twice.
    let del = table(&run(&[&m], &opts(true, true, 0), 1).0);
    expect(&del, "chr1", 1001, 1040, 1);
    // e2: a deletion that starts before the clip and runs past it counts from the clip on.
    expect(&del, "chr1", 1101, 1145, 1);
    let nodel = table(&run(&[&m], &opts(true, false, 0), 1).0);
    expect(&nodel, "chr1", 1101, 1130, 1);
    expect(&nodel, "chr1", 1131, 1140, 0);
    expect(&nodel, "chr1", 1141, 1145, 1);

    // h under -q 10: the second mate's bases that the first does not cover use their own
    // qualities (10 bases of quality 5 then 10 of 40).
    let q = table(&run(&[&m], &opts(true, false, 10), 1).0);
    expect(&q, "chr1", 1601, 1650, 1);
    expect(&q, "chr1", 1651, 1660, 0);
    expect(&q, "chr1", 1661, 1670, 1);
    expect(&plain, "chr1", 1651, 1660, 1);

    // Without the option every base counts: a's overlap is 2.
    let all = table(&run(&[&m], &opts(false, false, 0), 1).0);
    expect(&all, "chr1", 151, 200, 2);
}

#[test]
fn stale_entry_template_is_not_clipped() {
    // g is a pair followed by two supplementary alignments of one name inside the second mate.
    // samtools keeps the first supplementary as a pending mate (its mate position is behind it)
    // and clips the second by its end; gloci holds neither (#depth-mate-overlap).
    let m = mates();
    let t = table(&run(&[&m], &opts(true, false, 0), 1).0);
    expect(&t, "chr1", 1401, 1500, 1);
    expect(&t, "chr1", 1501, 1519, 2);
    expect(&t, "chr1", 1520, 1530, 3);
    expect(&t, "chr1", 1531, 1539, 2);
    expect(&t, "chr1", 1540, 1550, 1);
    // m: a read waiting for a mate at the end of chr1 does not clip a read of that name on chr2.
    expect(&t, "chr2", 301, 350, 1);
}

fn sam_text(pairs: u64, spacing: u64, len: u64, gap: u64) -> String {
    let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:100000000\n");
    // Pair i: first mate at i * spacing, second `gap` bases later; sorted by merging the two.
    let mut next_second = 0u64;
    let mut line = |name: u64, flag: u16, pos: u64, mate: u64| {
        text.push_str(&format!(
            "p{name}\t{flag}\tc\t{}\t60\t{len}M\t=\t{}\t0\t*\t*\n",
            pos + 1,
            mate + 1
        ));
    };
    for i in 0..pairs {
        let pos = i * spacing;
        while next_second < i && next_second * spacing + gap <= pos {
            let s = next_second * spacing;
            line(next_second, 147, s + gap, s);
            next_second += 1;
        }
        line(i, 99, pos, pos + gap);
    }
    while next_second < pairs {
        let s = next_second * spacing;
        line(next_second, 147, s + gap, s);
        next_second += 1;
    }
    text
}

#[test]
fn mate_table_bounded() {
    // 1,000,000 pairs of two 100-base reads 50 apart, a pair every 5 bases: depth 30 (the
    // mates overlap), and at most 20 or so first mates wait for theirs at once.
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("pairs.sam");
    std::fs::write(&p, sam_text(1_000_000, 5, 100, 50)).unwrap();
    let paths = vec![p.display().to_string()];
    let mut sink = std::io::sink();
    let s = depth(
        &paths,
        &opts(true, false, 0),
        &ro(),
        NonZeroUsize::MIN,
        &mut sink,
    )
    .unwrap();
    assert!(s.rows > 5_000_000, "{}", s.rows);
    assert!(s.max_pending_mates > 0);
    assert!(s.max_pending_mates <= 60, "{}", s.max_pending_mates);
}

#[test]
fn mate_table_holds_only_reads_with_a_mate_ahead() {
    // d (a pair that does not overlap) and f (no mate position) are never held; the table of
    // a file of only those stays empty.
    let dir = tempfile::tempdir().unwrap();
    let text = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:2000\n\
d\t99\tc\t701\t60\t50M\t=\t901\t0\t*\t*\n\
f\t65\tc\t801\t60\t50M\t*\t0\t0\t*\t*\n\
d\t147\tc\t901\t60\t50M\t=\t701\t0\t*\t*\n";
    let p = dir.path().join("d.sam");
    std::fs::write(&p, text).unwrap();
    let paths = vec![p.display().to_string()];
    let (out, s) = {
        let mut out = Vec::new();
        let s = depth(
            &paths,
            &opts(true, false, 0),
            &ro(),
            NonZeroUsize::MIN,
            &mut out,
        )
        .unwrap();
        (String::from_utf8(out).unwrap(), s)
    };
    assert_eq!(s.max_pending_mates, 0);
    assert_eq!(out.lines().count(), 150);
}

#[test]
fn mates_per_chromosome_equals_streaming() {
    let m = mates();
    let both = [m.as_str(), m.as_str()];
    for paths in [&both[..1], &both[..]] {
        for mode in [
            AllPositions::Covered,
            AllPositions::RefsWithReads,
            AllPositions::AllRefs,
        ] {
            for regions in [
                None,
                Some(vec![
                    (
                        b"chr1".to_vec(),
                        Interval {
                            start: 90,
                            end: 700,
                        },
                    ),
                    (
                        b"chr1".to_vec(),
                        Interval {
                            start: 1000,
                            end: 1700,
                        },
                    ),
                    (b"chr2".to_vec(), Interval { start: 0, end: 400 }),
                ]),
            ] {
                for (deletions, q) in [(false, 0), (true, 0), (true, 10)] {
                    let mut o = opts(true, deletions, q);
                    o.all = mode;
                    o.regions = regions.clone();
                    let (want, ws) = run(paths, &o, 1);
                    for t in [2, 8] {
                        let (got, gs) = run(paths, &o, t);
                        assert_eq!(got, want, "-t {t} {mode:?} {:?}", o.regions);
                        assert_eq!(gs.rows, ws.rows);
                        assert_eq!(gs.max_pending_mates, ws.max_pending_mates);
                    }
                }
            }
        }
    }
}
