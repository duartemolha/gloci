//! The complement table, and `getfasta` streaming a large record against the
//! naive answer.
use gloci_core::ChromTable;
use gloci_io::{FETCH_CHUNK, IndexedFasta, InvalidPolicy, ReaderOptions, write_fasta_index};
use gloci_ops::getfasta::{FastaName, GetfastaOptions, SeqLayout, getfasta};
use gloci_ops::seq::{SeqSummary, complement, revcomp_in_place};
use gloci_testkit::mem::mem_reader;
use std::io::Write;
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

#[test]
fn revcomp_table() {
    // (byte, DNA complement); lower case follows the upper-case rule.
    let pairs = [
        (b'A', b'T'),
        (b'T', b'A'),
        (b'C', b'G'),
        (b'G', b'C'),
        (b'U', b'A'),
        (b'R', b'Y'),
        (b'Y', b'R'),
        (b'K', b'M'),
        (b'M', b'K'),
        (b'B', b'V'),
        (b'V', b'B'),
        (b'D', b'H'),
        (b'H', b'D'),
    ];
    let mut mapped = [false; 256];
    for (b, c) in pairs {
        assert_eq!(complement(b, false), c, "{}", b as char);
        assert_eq!(
            complement(b.to_ascii_lowercase(), false),
            c.to_ascii_lowercase()
        );
        mapped[b as usize] = true;
        mapped[b.to_ascii_lowercase() as usize] = true;
    }
    for b in 0..=255u8 {
        if !mapped[b as usize] {
            assert_eq!(complement(b, false), b, "byte {b} is unchanged");
        }
        // RNA differs from DNA only for A.
        let want = match b {
            b'A' => b'U',
            b'a' => b'u',
            _ => complement(b, false),
        };
        assert_eq!(complement(b, true), want, "rna byte {b}");
    }
    for b in *b"SWNXswnx.-*" {
        assert_eq!(complement(b, false), b);
    }
    let mut s = *b"ACGTRYkmNuX";
    revcomp_in_place(&mut s, false);
    assert_eq!(&s, b"XaNkmRYACGT");
    let mut r = *b"AAa";
    revcomp_in_place(&mut r, true);
    assert_eq!(&r, b"uUU");
    revcomp_in_place(&mut [], false);
}

/// Counts what is written and the largest single write.
#[derive(Default)]
struct Probe {
    data: Vec<u8>,
    biggest: usize,
}
impl Write for Probe {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.biggest = self.biggest.max(b.len());
        self.data.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

fn naive_revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| complement(b, false)).collect()
}

#[test]
fn getfasta_streams_large_records() {
    const LEN: usize = 5 << 20;
    let mut x = 88172645463325252u64;
    let seq: Vec<u8> = (0..LEN)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            b"ACGTacgtRYKMNnBDHV"[(x % 18) as usize]
        })
        .collect();
    let mut fa = b">big\n".to_vec();
    for line in seq.chunks(60) {
        fa.extend_from_slice(line);
        fa.push(b'\n');
    }
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("big.fa");
    std::fs::write(&path, &fa).unwrap();
    write_fasta_index(&path, NonZeroUsize::MIN).unwrap();
    let bed = format!(
        "big\t0\t{LEN}\tplus\t0\t+\nbig\t0\t{LEN}\tminus\t0\t-\nbig\t7\t{}\tpart\t0\t-\n",
        LEN - 11
    );
    let mut fasta = IndexedFasta::open(&path).unwrap();
    let o = GetfastaOptions {
        strand: true,
        name: FastaName::NameOnly,
        split: false,
        layout: SeqLayout::Fasta,
        rna: false,
        skip_outside: false,
    };
    let mut out = Probe::default();
    let sum = getfasta(
        mem_reader("big.bed", &bed, &ro()),
        &mut fasta,
        &mut ChromTable::first_appearance(),
        &o,
        &mut out,
    )
    .unwrap();
    assert_eq!(
        sum,
        SeqSummary {
            written: 3,
            points: 0,
            skipped: 0
        }
    );
    let mut want = Vec::new();
    for (name, body) in [
        ("plus(+)", seq.clone()),
        ("minus(-)", naive_revcomp(&seq)),
        ("part(-)", naive_revcomp(&seq[7..LEN - 11])),
    ] {
        want.extend_from_slice(format!(">{name}\n").as_bytes());
        want.extend_from_slice(&body);
        want.push(b'\n');
    }
    assert!(out.data == want, "sequence differs from the naive answer");
    // No write (so no buffer held) is larger than one fetch chunk.
    assert!(
        out.biggest as u64 <= FETCH_CHUNK,
        "largest write {} bytes",
        out.biggest
    );
}
