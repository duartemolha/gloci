//! The depth core (D137, D138): what counts per position, which rows print, the modes of
//! `--all-*`, several inputs, the order and reference checks, the memory bound, and per-chromosome
//! runs equal to one stream.

use gloci_core::{Interval, VcfSpan};
use gloci_io::{
    Accept, AlignOptions, Format, InvalidPolicy, ReaderOptions, aln_ref, open_aln_source_from,
};
use gloci_ops::depth::{
    AllPositions, CountRules, DepthOptions, DepthRing, DepthSummary, depth, depth_regions,
};
use std::io::Cursor;
use std::num::NonZeroUsize;
use std::path::PathBuf;

fn fixture(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/samtools")
        .join(name)
        .display()
        .to_string()
}

fn ro(exclude: u16) -> ReaderOptions {
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
        invalid: InvalidPolicy::Fail,
    }
}

/// samtools depth's default exclusion: UNMAP, SECONDARY, QCFAIL, DUP.
const DEFAULT_EXCLUDE: u16 = 0x704;

fn opts() -> DepthOptions {
    DepthOptions {
        rules: CountRules::default(),
        min_read_length: 0,
        all: AllPositions::Covered,
        regions: None,
        skip_mate_overlap: false,
        column_names: false,
    }
}

fn threads(n: usize) -> NonZeroUsize {
    NonZeroUsize::new(n).unwrap()
}

fn run(paths: &[&str], o: &DepthOptions, t: usize) -> (String, DepthSummary) {
    let paths: Vec<String> = paths.iter().map(|p| p.to_string()).collect();
    let mut out = Vec::new();
    let s = depth(&paths, o, &ro(DEFAULT_EXCLUDE), threads(t), &mut out).unwrap();
    (String::from_utf8(out).unwrap(), s)
}

fn run_err(paths: &[&str], o: &DepthOptions, t: usize) -> String {
    let paths: Vec<String> = paths.iter().map(|p| p.to_string()).collect();
    let mut out = Vec::new();
    depth(&paths, o, &ro(DEFAULT_EXCLUDE), threads(t), &mut out)
        .unwrap_err()
        .to_string()
}

/// The depths a lone read puts on the positions from its start to the end of its span.
fn ring_depths(cigar: &str, qual: &str, seq: &str, rules: CountRules) -> Vec<u32> {
    let sam = format!(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:100000\nr\t0\tc\t11\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}\n"
    );
    let mut src =
        open_aln_source_from("t.sam", Box::new(Cursor::new(sam.into_bytes())), &ro(0)).unwrap();
    let mut buf = Vec::new();
    src.next_record(&mut buf).unwrap().unwrap();
    let hdr = src.header().clone();
    let a = aln_ref(Format::Sam, &buf, &hdr).unwrap();
    let (span, _) = a.span(100_000).unwrap();
    let mut ring = DepthRing::new();
    ring.add(&a, span, &rules, 0);
    (span.start..ring.covered_end())
        .map(|p| ring.take(p))
        .collect()
}

fn rules(deletions: bool, min_baseq: u8) -> CountRules {
    CountRules {
        deletions,
        ref_skips: false,
        min_baseq,
    }
}

fn bases(n: usize) -> String {
    "ACGT".repeat(n / 4 + 1)[..n].to_string()
}

/// The query bases a CIGAR uses (M, I, S, =, X).
fn query_len(cigar: &str) -> usize {
    let mut n = 0;
    let mut total = 0;
    for c in cigar.chars() {
        if let Some(d) = c.to_digit(10) {
            n = n * 10 + d as usize;
        } else {
            if "MIS=X".contains(c) {
                total += n;
            }
            n = 0;
        }
    }
    total
}

struct Row {
    cigar: &'static str,
    /// `*` for a missing QUAL (the read keeps its sequence); `SEQ*` for a read with no SEQ.
    qual: String,
    rules: CountRules,
    want: Vec<u32>,
}

fn row(cigar: &'static str, qual: &str, rules: CountRules, want: Vec<u32>) -> Row {
    Row {
        cigar,
        qual: qual.to_string(),
        rules,
        want,
    }
}

#[test]
fn depth_counting_table() {
    let q30 = |n: usize| "?".repeat(n);
    let table = vec![
        row("5M", &q30(5), rules(false, 20), vec![1; 5]),
        row("5M", "??+??", rules(false, 20), vec![1, 1, 0, 1, 1]),
        row(
            "2M3D2M",
            &q30(4),
            rules(false, 0),
            vec![1, 1, 0, 0, 0, 1, 1],
        ),
        row("2M3D2M", &q30(4), rules(true, 0), vec![1; 7]),
        // The base after the deletion has quality 10: the deleted positions fail -q 20, as
        // does that base (the base before the deletion is fine).
        row("2M3D2M", "??+?", rules(true, 20), vec![1, 1, 0, 0, 0, 0, 1]),
        // A good base after it passes the deleted positions, though the one before fails.
        row("2M3D2M", "?+??", rules(true, 20), vec![1, 0, 1, 1, 1, 1, 1]),
        row(
            "3M100N3M",
            &q30(6),
            rules(false, 0),
            [vec![1; 3], vec![0; 100], vec![1; 3]].concat(),
        ),
        // N never counts, even with -J.
        row(
            "3M100N3M",
            &q30(6),
            rules(true, 0),
            [vec![1; 3], vec![0; 100], vec![1; 3]].concat(),
        ),
        row("2S3M2I3M", &q30(10), rules(false, 0), vec![1; 6]),
        // Soft clips and insertions advance the query: index 2 is the first M base.
        row(
            "2S3M2I3M",
            "??+???????",
            rules(false, 20),
            vec![0, 1, 1, 1, 1, 1],
        ),
        row(
            "2S3M2I3M",
            "???????+??",
            rules(false, 20),
            vec![1, 1, 1, 0, 1, 1],
        ),
        // An inserted base is no position: its quality gates nothing.
        row(
            "2S3M2I3M",
            "?????+????",
            rules(false, 20),
            vec![1, 1, 1, 1, 1, 1],
        ),
        row("2=1X2=", &q30(5), rules(false, 20), vec![1; 5]),
        // QUAL * is 0xFF for every base: it passes -q 50.
        row("4M", "*", rules(false, 50), vec![1; 4]),
        row("3M1D3M", "*", rules(true, 50), vec![1; 7]),
        // SEQ *: every base passes, samtools reads past the record there.
        row("4M", "SEQ*", rules(false, 20), vec![1; 4]),
        row("3M2D3M", "SEQ*", rules(true, 20), vec![1; 8]),
        // H and P touch neither the query nor the depth.
        row("2H3M1P3M2H", &q30(6), rules(false, 0), vec![1; 6]),
    ];
    for r in table {
        let (seq, qual) = match r.qual.as_str() {
            "SEQ*" => ("*".to_string(), "*".to_string()),
            "*" => (bases(query_len(r.cigar)), "*".to_string()),
            q => (bases(q.len()), q.to_string()),
        };
        let got = ring_depths(r.cigar, &qual, &seq, r.rules);
        assert_eq!(got, r.want, "{} {:?} {:?}", r.cigar, r.qual, r.rules);
    }
}

fn sam(sq: &[(&str, u64)], reads: &[(&str, u16, &str, u64, &str)]) -> String {
    let mut t = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for (n, l) in sq {
        t.push_str(&format!("@SQ\tSN:{n}\tLN:{l}\n"));
    }
    for (name, flag, chr, pos, cigar) in reads {
        let l = query_len(cigar);
        let (seq, qual) = if l == 0 {
            ("*".to_string(), "*".to_string())
        } else {
            (bases(l), "I".repeat(l))
        };
        t.push_str(&format!(
            "{name}\t{flag}\t{chr}\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}\n"
        ));
    }
    t
}

fn write(dir: &tempfile::TempDir, name: &str, text: &str) -> String {
    let p = dir.path().join(name);
    std::fs::write(&p, text).unwrap();
    p.display().to_string()
}

#[test]
fn depth_spliced_rows() {
    // A spliced read prints every intron position, with 0.
    let dir = tempfile::tempdir().unwrap();
    let f = write(
        &dir,
        "s.sam",
        &sam(&[("chr1", 1000)], &[("r", 0, "chr1", 10, "3M100N3M")]),
    );
    let (out, s) = run(&[&f], &opts(), 1);
    let rows: Vec<&str> = out.lines().collect();
    assert_eq!(rows.len(), 106);
    assert_eq!(rows[0], "chr1\t10\t1");
    assert_eq!(rows[3], "chr1\t13\t0");
    assert_eq!(rows[102], "chr1\t112\t0");
    assert_eq!(rows[103], "chr1\t113\t1");
    assert_eq!(s.rows, 106);

    // In depth.bam the spliced read starts at 60 and ends at 165: no position of its span is
    // left out, in order, each once.
    let (out, _) = run(&[&fixture("depth.bam")], &opts(), 1);
    let pos: Vec<u64> = out
        .lines()
        .filter(|l| l.starts_with("chr1\t"))
        .map(|l| l.split('\t').nth(1).unwrap().parse().unwrap())
        .collect();
    assert!(pos.windows(2).all(|w| w[0] < w[1]));
    for p in 60..=165 {
        assert!(pos.contains(&p), "position {p}");
    }
    // Position 150 has only the intron under it.
    assert!(out.contains("chr1\t150\t0\n"));
}

fn count_rows(out: &str, chrom: &str) -> usize {
    out.lines()
        .filter(|l| l.split('\t').next() == Some(chrom))
        .count()
}

#[test]
fn depth_all_modes() {
    let f = fixture("depth.bam");
    let mut o = opts();
    let (out, _) = run(&[&f], &o, 1);
    assert_eq!(count_rows(&out, "chr3"), 0);
    assert!(count_rows(&out, "chr1") < 500);
    o.all = AllPositions::RefsWithReads;
    let (out, _) = run(&[&f], &o, 1);
    assert_eq!(count_rows(&out, "chr1"), 500);
    assert_eq!(count_rows(&out, "chr2"), 300);
    assert_eq!(count_rows(&out, "chr3"), 0);
    o.all = AllPositions::AllRefs;
    let (out, _) = run(&[&f], &o, 1);
    assert_eq!(count_rows(&out, "chr3"), 200);
    assert!(out.ends_with("chr3\t200\t0\n"));

    // A region on the empty chr3 prints its positions, with -a.
    o.all = AllPositions::RefsWithReads;
    o.regions = Some(vec![(b"chr3".to_vec(), Interval { start: 0, end: 10 })]);
    let (out, _) = run(&[&f], &o, 1);
    let want: String = (1..=10).map(|p| format!("chr3\t{p}\t0\n")).collect();
    assert_eq!(out, want);
    // Without -a, nothing: no read is there.
    o.all = AllPositions::Covered;
    assert_eq!(run(&[&f], &o, 1).0, "");
}

#[test]
fn depth_inputs() {
    let (out, _) = run(&[&fixture("depth.bam"), &fixture("depth2.bam")], &opts(), 1);
    // chr3 has reads in the second file only.
    assert!(out.contains("chr3\t20\t0\t1\n"));
    assert!(out.contains("chr1\t12\t1\t1\n"));
    assert!(out.lines().all(|l| l.split('\t').count() == 4));
    // The names of the inputs are the column names.
    let mut o = opts();
    o.column_names = true;
    let (a, b) = (fixture("depth.bam"), fixture("depth2.bam"));
    let (out, _) = run(&[&a, &b], &o, 1);
    assert!(out.starts_with("#CHROM\tPOS\t"));
    assert!(out.lines().next().unwrap().ends_with("depth2.bam"));

    let e = run_err(
        &[&fixture("depth-reordered.bam"), &fixture("depth.bam")],
        &opts(),
        1,
    );
    assert!(e.contains("@SQ 1"), "{e}");
    assert!(e.contains("chr2") && e.contains("chr1"), "{e}");
    let e = run_err(&[&fixture("depth-unsorted.bam")], &opts(), 1);
    assert!(e.contains("gloci sort"), "{e}");
    assert!(e.contains("chr1:17") && e.contains("chr1:20"), "{e}");
}

#[test]
fn depth_past_ln() {
    let (out, _) = run(&[&fixture("depth-pastend.bam")], &opts(), 1);
    // The reads at chr2:290 for 20 bases and chr2:300 for 10 are cut at LN 300.
    assert!(out.contains("chr2\t300\t2\n"));
    assert!(out.contains("chr2\t299\t1\n"));
    assert!(!out.contains("chr2\t301\t"));
    // A read that starts at LN or beyond is dropped.
    let dir = tempfile::tempdir().unwrap();
    let f = write(
        &dir,
        "p.sam",
        &sam(
            &[("c", 10)],
            &[("r", 0, "c", 10, "5M"), ("s", 0, "c", 11, "5M")],
        ),
    );
    let (out, _) = run(&[&f], &opts(), 1);
    assert_eq!(out, "c\t10\t1\n");
}

#[test]
fn depth_memory_bounded() {
    // One read with a 500 kb intron, then 100,000 reads of 150 bases: the ring holds what the
    // longest read spans, not the reference.
    let dir = tempfile::tempdir().unwrap();
    let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:100000000\n");
    let seq = bases(150);
    let qual = "I".repeat(150);
    text.push_str("big\t0\tc\t1\t60\t5M500000N5M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n");
    for i in 0..100_000u64 {
        text.push_str(&format!(
            "r{i}\t0\tc\t{}\t60\t150M\t*\t0\t0\t{seq}\t{qual}\n",
            2 + i * 5
        ));
    }
    let f = write(&dir, "m.sam", &text);
    let paths = vec![f];
    let mut sink = std::io::sink();
    let s = depth(&paths, &opts(), &ro(DEFAULT_EXCLUDE), threads(1), &mut sink).unwrap();
    assert!(s.rows > 500_000);
    assert!(s.max_ring <= 1 << 20, "{}", s.max_ring);
    assert!(s.max_ring >= 500_000, "{}", s.max_ring);
}

/// `reads` (a long read first, then 20,000 short ones) in `inputs` files: no ring is larger
/// than `next_pow2(longest span + 1)` counters (at least 64), whatever reaches it.
fn assert_ring_bound(cigar: &str, span: u64, inputs: usize, deletions: bool) {
    let dir = tempfile::tempdir().unwrap();
    let l = query_len(cigar);
    let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:100000000\n");
    text.push_str(&format!(
        "big\t0\tc\t1\t60\t{cigar}\t*\t0\t0\t{}\t{}\n",
        bases(l),
        "I".repeat(l)
    ));
    for i in 0..20_000u64 {
        text.push_str(&format!(
            "r{i}\t0\tc\t{}\t60\t50M\t*\t0\t0\t{}\t{}\n",
            2 + i * 7,
            bases(50),
            "I".repeat(50)
        ));
    }
    let paths: Vec<String> = (0..inputs)
        .map(|k| write(&dir, &format!("b{k}.sam"), &text))
        .collect();
    let mut o = opts();
    o.rules.deletions = deletions;
    let mut sink = std::io::sink();
    let s = depth(&paths, &o, &ro(DEFAULT_EXCLUDE), threads(1), &mut sink).unwrap();
    let bound = usize::try_from(span + 1)
        .unwrap()
        .next_power_of_two()
        .max(64);
    assert!(
        s.max_ring <= bound,
        "{cigar} x{inputs}: {} > {bound}",
        s.max_ring
    );
    // And it really holds the span (the bound is not met by dropping counts).
    assert!(s.max_ring >= span as usize, "{cigar}: {}", s.max_ring);
}

#[test]
fn depth_ring_bound_holds_for_long_m_d_and_two_inputs() {
    // A long M read, a long D read (counted with -J), and the same with two inputs.
    assert_ring_bound("400000M", 400_000, 1, false);
    assert_ring_bound("5M300000D5M", 300_010, 1, true);
    assert_ring_bound("400000M", 400_000, 2, false);
    assert_ring_bound("5M300000D5M", 300_010, 2, true);
}

fn plans() -> Vec<DepthOptions> {
    let mut all = Vec::new();
    for mode in [
        AllPositions::Covered,
        AllPositions::RefsWithReads,
        AllPositions::AllRefs,
    ] {
        for regions in [
            None,
            Some(vec![
                (b"chr1".to_vec(), Interval { start: 5, end: 130 }),
                (
                    b"chr1".to_vec(),
                    Interval {
                        start: 120,
                        end: 180,
                    },
                ),
                (b"chr2".to_vec(), Interval { start: 0, end: 40 }),
                (b"chr3".to_vec(), Interval { start: 10, end: 20 }),
            ]),
        ] {
            for (rules, len) in [(CountRules::default(), 0), (rules(true, 20), 31)] {
                let mut o = opts();
                o.all = mode;
                o.regions = regions.clone();
                o.rules = rules;
                o.min_read_length = len;
                all.push(o);
            }
        }
    }
    all
}

#[test]
fn depth_per_chromosome_equals_streaming() {
    let one = [fixture("depth.bam")];
    let two = [fixture("depth.bam"), fixture("depth2.bam")];
    for files in [&one[..], &two[..]] {
        let paths: Vec<&str> = files.iter().map(String::as_str).collect();
        for o in plans() {
            // One thread streams; more read each run of references through the indexes.
            let (want, ws) = run(&paths, &o, 1);
            for t in [2, 8] {
                let (got, gs) = run(&paths, &o, t);
                assert_eq!(got, want, "-t {t} {:?} {:?}", o.all, o.regions);
                assert_eq!(gs.rows, ws.rows);
            }
        }
    }
}

#[test]
fn depth_regions_resolve_and_drop_unknown() {
    let f = fixture("depth.bam");
    let r = depth_regions(
        &f,
        &["chr1:100-200".into(), "nope:1-5".into(), "chr2".into()],
        None,
        &ro(0),
    )
    .unwrap()
    .unwrap();
    assert_eq!(
        r,
        vec![
            (
                b"chr1".to_vec(),
                Interval {
                    start: 99,
                    end: 200
                }
            ),
            (
                b"chr2".to_vec(),
                Interval {
                    start: 0,
                    end: u64::MAX
                }
            ),
        ]
    );
    assert!(depth_regions(&f, &[], None, &ro(0)).unwrap().is_none());
    assert!(depth_regions("-", &["chr1".into()], None, &ro(0)).is_err());
}

#[test]
fn depth_all_refs_huge_header() {
    // 200,000 references of 10 bases with --all-references: 2,000,000 rows, in time linear in
    // the references.
    let dir = tempfile::tempdir().unwrap();
    let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n");
    for i in 0..200_000 {
        text.push_str(&format!("@SQ\tSN:ctg{i}\tLN:10\n"));
    }
    text.push_str("a\t0\tctg5\t3\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
    text.push_str("b\t0\tctg150000\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
    let f = write(&dir, "h.sam", &text);
    let mut o = opts();
    o.all = AllPositions::AllRefs;
    let paths = vec![f];
    let started = std::time::Instant::now();
    let mut sink = std::io::sink();
    let s = depth(&paths, &o, &ro(DEFAULT_EXCLUDE), threads(1), &mut sink).unwrap();
    assert_eq!(s.rows, 2_000_000);
    assert!(
        started.elapsed() < std::time::Duration::from_secs(20),
        "{:?}",
        started.elapsed()
    );
}

#[test]
fn depth_ring_equals_naive_counts() {
    // Random reads (M, D, N, I, S, quality gates) added in order to one ring, which is read at
    // every position as the rows do, give the counts of a plain array.
    let mut x = 0x9E37_79B9_7F4A_7C15u64;
    let mut next = move |m: u64| {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x % m
    };
    for round in 0..40 {
        let min_baseq = [0u8, 20, 35][round % 3];
        let deletions = round % 2 == 0;
        let mut text = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c\tLN:100000\n");
        let mut pos = 1u64;
        let mut reads = Vec::new();
        for i in 0..60 {
            pos += next(40);
            let mut cigar = String::new();
            if next(3) == 0 {
                cigar.push_str(&format!("{}S", 1 + next(3)));
            }
            let blocks = 1 + next(3);
            for b in 0..blocks {
                cigar.push_str(&format!("{}M", 1 + next(30)));
                if b + 1 < blocks {
                    match next(3) {
                        0 => cigar.push_str(&format!("{}D", 1 + next(5))),
                        1 => cigar.push_str(&format!("{}N", 1 + next(200))),
                        _ => cigar.push_str(&format!("{}I", 1 + next(4))),
                    }
                }
            }
            let l = query_len(&cigar);
            let qual: String = (0..l)
                .map(|_| char::from(33 + 5 + next(40) as u8))
                .collect();
            text.push_str(&format!(
                "r{i}\t0\tc\t{pos}\t60\t{cigar}\t*\t0\t0\t{}\t{qual}\n",
                bases(l)
            ));
            reads.push((pos - 1, cigar, qual));
        }
        // The plain array.
        let mut want = vec![0u32; 20_000];
        let gate = |qual: &str, q: usize| min_baseq == 0 || (qual.as_bytes()[q] - 33) >= min_baseq;
        let mut span_end = Vec::new();
        for (start, cigar, qual) in &reads {
            let (mut p, mut q, mut n) = (*start as usize, 0usize, 0usize);
            for c in cigar.chars() {
                if let Some(d) = c.to_digit(10) {
                    n = n * 10 + d as usize;
                    continue;
                }
                match c {
                    'M' => {
                        for k in 0..n {
                            if gate(qual, q + k) {
                                want[p + k] += 1;
                            }
                        }
                        p += n;
                        q += n;
                    }
                    'D' => {
                        if deletions && gate(qual, q) {
                            for k in 0..n {
                                want[p + k] += 1;
                            }
                        }
                        p += n;
                    }
                    'N' => p += n,
                    _ => q += n,
                }
                n = 0;
            }
            span_end.push(p);
        }
        // The ring, with every position taken in order between the reads.
        let mut src =
            open_aln_source_from("r.sam", Box::new(Cursor::new(text.into_bytes())), &ro(0))
                .unwrap();
        let hdr = src.header().clone();
        let mut ring = DepthRing::new();
        let rules = rules(deletions, min_baseq);
        let mut buf = Vec::new();
        let mut taken = 0u64;
        let mut got = vec![0u32; 20_000];
        let mut k = 0;
        loop {
            buf.clear();
            if src.next_record(&mut buf).unwrap().is_none() {
                break;
            }
            let a = aln_ref(Format::Sam, &buf, &hdr).unwrap();
            let start = a.pos().unwrap();
            while taken < start {
                got[taken as usize] = ring.take(taken);
                taken += 1;
            }
            let span = Interval {
                start,
                end: span_end[k] as u64,
            };
            ring.add(&a, span, &rules, 0);
            k += 1;
        }
        while taken < ring.covered_end() {
            got[taken as usize] = ring.take(taken);
            taken += 1;
        }
        assert_eq!(got, want, "round {round}");
        assert!(ring.capacity() <= 1 << 12);
    }
}
