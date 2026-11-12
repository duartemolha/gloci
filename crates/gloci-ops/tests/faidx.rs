//! `faidx`: wrapping (also across fetch chunks), strand marks and the
//! complement table, against naive references.
use gloci_core::Interval;
use gloci_io::{FETCH_CHUNK, IndexedFasta, write_fasta_index};
use gloci_ops::faidx::{FaidxOptions, FaidxRegion, MarkStrand, faidx};
use gloci_ops::seq::complement;
use std::num::NonZeroUsize;
use std::path::PathBuf;

/// A FASTA of `(name, sequence, line width)` records, indexed.
fn fasta(dir: &tempfile::TempDir, recs: &[(&str, &[u8], usize)]) -> IndexedFasta {
    let p: PathBuf = dir.path().join("t.fa");
    let mut text = Vec::new();
    for (name, seq, w) in recs {
        text.extend_from_slice(format!(">{name}\n").as_bytes());
        for line in seq.chunks(*w) {
            text.extend_from_slice(line);
            text.push(b'\n');
        }
    }
    std::fs::write(&p, text).unwrap();
    write_fasta_index(&p, NonZeroUsize::MIN).unwrap();
    IndexedFasta::open(&p).unwrap()
}

fn pattern(n: usize) -> Vec<u8> {
    (0..n)
        .map(|i| b"ACGTNacgtRYKM"[(i * 7 + i / 13) % 13])
        .collect()
}

fn region(label: &str, chrom: &str, start: u64, end: u64) -> FaidxRegion {
    FaidxRegion {
        label: label.into(),
        chrom: chrom.as_bytes().to_vec(),
        iv: Interval { start, end },
    }
}

fn run(f: &mut IndexedFasta, regions: &[FaidxRegion], o: &FaidxOptions) -> Vec<u8> {
    let mut out = Vec::new();
    faidx(f, regions, o, &mut out).unwrap();
    out
}

fn wrapped(seq: &[u8], w: usize) -> Vec<u8> {
    let mut out = Vec::new();
    if w == 0 {
        out.extend_from_slice(seq);
        out.push(b'\n');
    } else {
        for l in seq.chunks(w) {
            out.extend_from_slice(l);
            out.push(b'\n');
        }
    }
    out
}

#[test]
fn faidx_wrapping() {
    let d = tempfile::tempdir().unwrap();
    let a = pattern(500);
    let b: Vec<u8> = pattern(333).into_iter().rev().collect();
    let mut f = fasta(&d, &[("a", &a, 60), ("b", &b, 70)]);
    let rs = [region("a", "a", 0, u64::MAX), region("b", "b", 0, u64::MAX)];
    // default: each sequence's own .fai line width
    let got = run(&mut f, &rs, &FaidxOptions::default());
    let mut want = b">a\n".to_vec();
    want.extend(wrapped(&a, 60));
    want.extend(b">b\n");
    want.extend(wrapped(&b, 70));
    assert_eq!(got, want);
    // an explicit width, and 0 for one line
    for w in [7u64, 500, 1000, 0] {
        let o = FaidxOptions {
            line_width: Some(w),
            ..Default::default()
        };
        let got = run(&mut f, &rs[..1], &o);
        let mut want = b">a\n".to_vec();
        want.extend(wrapped(&a, w as usize));
        assert_eq!(got, want, "width {w}");
    }
}

#[test]
fn faidx_wrapping_across_chunks() {
    // 3 MiB and a bit: every output line is full except the last, whatever
    // the fetch chunk and the line width (also one that divides neither).
    let d = tempfile::tempdir().unwrap();
    let len = 3 * FETCH_CHUNK as usize + 4321;
    let seq = pattern(len);
    let mut f = fasta(&d, &[("big", &seq, 60)]);
    let rs = [region("big", "big", 0, u64::MAX)];
    for w in [60u64, 61, 1000, 65_537] {
        let o = FaidxOptions {
            line_width: Some(w),
            ..Default::default()
        };
        let got = run(&mut f, &rs, &o);
        let body = &got[5..];
        let mut want = Vec::new();
        want.extend(wrapped(&seq, w as usize));
        assert_eq!(body, &want[..], "width {w}");
        let lines: Vec<&[u8]> = body.split(|&c| c == b'\n').collect();
        let (last, full) = lines[..lines.len() - 1].split_last().unwrap();
        assert!(full.iter().all(|l| l.len() == w as usize), "width {w}");
        assert!(last.len() <= w as usize && !last.is_empty());
    }
    // the reverse complement carries the column the same way
    let o = FaidxOptions {
        line_width: Some(77),
        reverse_complement: true,
        ..Default::default()
    };
    let got = run(&mut f, &rs, &o);
    let rc: Vec<u8> = seq.iter().rev().map(|&b| complement(b, false)).collect();
    let mut want = b">big/rc\n".to_vec();
    want.extend(wrapped(&rc, 77));
    assert_eq!(got, want);
}

#[test]
fn faidx_strand_marks() {
    let d = tempfile::tempdir().unwrap();
    let mut f = fasta(&d, &[("a", b"AACCGGTT", 8)]);
    let rs = [region("a:1-4", "a", 0, 4)];
    for (mark, fwd, rev) in [
        ("rc", ">a:1-4\nAACC\n", ">a:1-4/rc\nGGTT\n"),
        ("no", ">a:1-4\nAACC\n", ">a:1-4\nGGTT\n"),
        ("sign", ">a:1-4(+)\nAACC\n", ">a:1-4(-)\nGGTT\n"),
        ("custom,P,N", ">a:1-4P\nAACC\n", ">a:1-4N\nGGTT\n"),
    ] {
        for (reverse_complement, want) in [(false, fwd), (true, rev)] {
            let o = FaidxOptions {
                reverse_complement,
                mark: mark.parse().unwrap(),
                ..Default::default()
            };
            assert_eq!(
                run(&mut f, &rs, &o),
                want.as_bytes(),
                "{mark} {reverse_complement}"
            );
        }
    }
    for bad in ["", "RC", "custom", "custom,+", "plus"] {
        assert!(bad.parse::<MarkStrand>().is_err(), "{bad:?}");
    }
    assert_eq!(
        "custom,,x,y".parse::<MarkStrand>(),
        Ok(MarkStrand::Custom("".into(), "x,y".into()))
    );
}

#[test]
fn faidx_edges() {
    let d = tempfile::tempdir().unwrap();
    let mut f = fasta(&d, &[("a", b"ACGTACGTAC", 4)]);
    let o = FaidxOptions::default();
    // an end past the sequence is cut; a start at or past it gives the header
    let rs = [
        region("x", "a", 6, 99),
        region("y", "a", 10, 99),
        region("z", "a", 50, 60),
        region("w", "a", 7, 7),
    ];
    assert_eq!(run(&mut f, &rs, &o), b">x\nGTAC\n>y\n>z\n>w\n");
    // a sequence the index lacks is an error, not an empty record
    let mut out = Vec::new();
    let e = faidx(&mut f, &[region("q:1-2", "nope", 0, 2)], &o, &mut out).unwrap_err();
    assert!(e.to_string().contains("nope"), "{e}");
    assert!(out.is_empty());
}

/// samtools 1.21's reverse complement of the bytes `!`..`~` without `>`, in
/// order, read back from `samtools faidx -i` (so reversed here): the table
/// `comp_base` of faidx.c, restricted to the bytes a FASTA line can hold.
const SAMTOOLS_RC: &[u8] = b"~}|{zrxwbaasyqponklmjidcfehgvt`_^]\\[ZRXWBAASYQPONKLMJIDCFEHGVT@?=<;:9876543210/.-,+*)('&%$#\"!";

#[test]
fn faidx_complement_table() {
    let input: Vec<u8> = (33u8..127).filter(|&b| b != b'>').collect();
    assert_eq!(input.len(), SAMTOOLS_RC.len());
    for (i, &b) in input.iter().enumerate() {
        let want = SAMTOOLS_RC[SAMTOOLS_RC.len() - 1 - i];
        assert_eq!(complement(b, false), want, "byte {b:#x} ({})", b as char);
    }
    // outside the graphic range samtools leaves a byte alone, as gloci does
    for b in (0u8..33).chain(127..=255) {
        assert_eq!(complement(b, false), b, "byte {b:#x}");
    }
}
