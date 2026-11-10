//! `nuc`: the percentage formatting, the base classes, pattern counting across
//! chunk boundaries and large records streamed against the naive answer.
use gloci_core::ChromTable;
use gloci_io::{FETCH_CHUNK, IndexedFasta, InvalidPolicy, ReaderOptions, write_fasta_index};
use gloci_ops::nuc::{BaseCounts, NucOptions, nuc, pct_f32};
use gloci_ops::seq::SeqSummary;
use gloci_testkit::mem::mem_reader;
use std::num::NonZeroUsize;

fn ro() -> ReaderOptions {
    ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: None,
        threads: NonZeroUsize::MIN,
        format: None,
        vcf_span: gloci_core::VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    }
}

fn opts() -> NucOptions {
    NucOptions {
        strand: false,
        seq: false,
        pattern: None,
        ignore_case: false,
        skip_outside: false,
    }
}

#[test]
fn pct_f32_table() {
    for (c, l, want) in [
        (27, 29, "0.931035"),
        (39, 41, "0.951219"),
        (16, 51, "0.313726"),
        (1, 3, "0.333333"),
        (0, 10, "0.000000"),
        (10, 10, "1.000000"),
    ] {
        assert_eq!(pct_f32(c, l), want, "{c}/{l}");
    }
    // f64 would print 0.313725 and 0.931034 here.
    assert_ne!(format!("{:.6}", 16f64 / 51f64), pct_f32(16, 51));
}

#[test]
fn base_counts_table() {
    let mut c = BaseCounts::default();
    c.add(b"ACGTUNacgtunRYKM-*");
    assert_eq!(
        c,
        BaseCounts {
            a: 2,
            c: 2,
            g: 2,
            t: 4,
            n: 2,
            other: 6
        }
    );
    assert_eq!(c.total(), 18);
    // Counting is additive over chunks.
    let mut d = BaseCounts::default();
    d.add(b"ACGTUNac");
    d.add(b"gtunRYKM-*");
    assert_eq!(c, d);
}

fn make_fasta(dir: &std::path::Path, seq: &[u8]) -> IndexedFasta {
    let mut fa = b">big\n".to_vec();
    for line in seq.chunks(60) {
        fa.extend_from_slice(line);
        fa.push(b'\n');
    }
    let path = dir.join("big.fa");
    std::fs::write(&path, &fa).unwrap();
    write_fasta_index(&path, NonZeroUsize::MIN).unwrap();
    IndexedFasta::open(&path).unwrap()
}

fn run(fasta: &mut IndexedFasta, bed: &str, o: &NucOptions) -> (String, SeqSummary) {
    let mut out = Vec::new();
    let sum = nuc(
        mem_reader("in.bed", bed, &ro()),
        fasta,
        &mut ChromTable::first_appearance(),
        o,
        &mut out,
    )
    .unwrap();
    (String::from_utf8(out).unwrap(), sum)
}

fn naive_overlapping(s: &[u8], p: &[u8], ic: bool) -> u64 {
    s.windows(p.len())
        .filter(|w| {
            if ic {
                w.eq_ignore_ascii_case(p)
            } else {
                *w == p
            }
        })
        .count() as u64
}

#[test]
fn pattern_across_chunks() {
    // A pattern straddling every multiple of the fetch chunk, and a run that
    // matches itself at every offset.
    let chunk = FETCH_CHUNK as usize;
    let len = 3 * chunk + 1234;
    let mut seq = vec![b'A'; len];
    for k in 1..=3 {
        seq[k * chunk - 2..k * chunk + 2].copy_from_slice(b"GATT");
    }
    seq[100..110].copy_from_slice(b"ccccccgatt");
    let d = tempfile::tempdir().unwrap();
    let mut fa = make_fasta(d.path(), &seq);
    for (pat, ic) in [
        ("GATT", false),
        ("gatt", true),
        ("gatt", false),
        ("TTA", false),
        ("AA", false),
        ("G", false),
    ] {
        let o = NucOptions {
            pattern: Some(pat.as_bytes().to_vec()),
            ignore_case: ic,
            ..opts()
        };
        let (text, _) = run(&mut fa, &format!("big\t0\t{len}\n"), &o);
        let last = text.lines().nth(1).unwrap().rsplit('\t').next().unwrap();
        assert_eq!(
            last.parse::<u64>().unwrap(),
            naive_overlapping(&seq, pat.as_bytes(), ic),
            "{pat} ic={ic}"
        );
    }
    // The GATT pattern is found once at each boundary and once in the lowercase run.
    let o = NucOptions {
        pattern: Some(b"GATT".to_vec()),
        ..opts()
    };
    let (text, _) = run(&mut fa, &format!("big\t0\t{len}\n"), &o);
    assert!(text.trim_end().ends_with("\t3"), "{text}");
}

#[test]
fn large_records_match_the_naive_answer() {
    const LEN: usize = 3 << 20;
    let mut x = 88172645463325252u64;
    let seq: Vec<u8> = (0..LEN)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            b"ACGTacgtRYKMNnUBDHV"[(x % 19) as usize]
        })
        .collect();
    let d = tempfile::tempdir().unwrap();
    let mut fa = make_fasta(d.path(), &seq);
    let o = NucOptions {
        strand: true,
        seq: true,
        pattern: Some(b"CG".to_vec()),
        ..opts()
    };
    let bed = format!(
        "big\t5\t{}\tp\t0\t+\nbig\t5\t{}\tm\t0\t-\n",
        LEN - 3,
        LEN - 3
    );
    let (text, sum) = run(&mut fa, &bed, &o);
    assert_eq!(sum.written, 2);
    let body = &seq[5..LEN - 3];
    let rc: Vec<u8> = body
        .iter()
        .rev()
        .map(|&b| gloci_ops::seq::complement(b, false))
        .collect();
    for (line, s) in text.lines().skip(1).zip([body.to_vec(), rc]) {
        let f: Vec<&str> = line.split('\t').collect();
        let mut want = BaseCounts::default();
        want.add(&s);
        assert_eq!(
            f[8..14],
            [
                want.a.to_string(),
                want.c.to_string(),
                want.g.to_string(),
                want.t.to_string(),
                want.n.to_string(),
                want.other.to_string()
            ]
        );
        assert_eq!(f[14], s.len().to_string());
        assert_eq!(f[15].as_bytes(), s.as_slice(), "sequence");
        assert_eq!(f[16], naive_overlapping(&s, b"CG", false).to_string());
    }
    // Nothing the size of a record is ever written in one piece.
    assert!(FETCH_CHUNK < LEN as u64);
}

#[test]
fn points_are_counted_and_skipped_and_the_header_stays() {
    let d = tempfile::tempdir().unwrap();
    let mut fa = make_fasta(d.path(), b"ACGTACGTAC");
    let (text, sum) = run(&mut fa, "big\t3\t3\nbig\t0\t4\nbig\t5\t5\n", &opts());
    assert_eq!(
        sum,
        SeqSummary {
            written: 1,
            points: 2,
            skipped: 0
        }
    );
    assert_eq!(text.lines().count(), 2);
    assert!(text.starts_with("#1_usercol\t2_usercol\t3_usercol\t4_pct_at\t"));
    assert!(text.ends_with("big\t0\t4\t0.500000\t0.500000\t1\t1\t1\t1\t0\t0\t4\n"));
    // A file with only a skipped first record still has the header.
    let (text, _) = run(&mut fa, "big\t3\t3\n", &opts());
    assert_eq!(text.lines().count(), 1);
    // No records, no header.
    let (text, _) = run(&mut fa, "", &opts());
    assert_eq!(text, "");
}
