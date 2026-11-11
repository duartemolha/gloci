//! Alignment sort (D121): spilled sorts equal in-memory ones for every order and thread
//! count, a record larger than the memory budget, the header rules and `--check`.

use gloci_core::VcfSpan;
use gloci_io::{
    Accept, AlignOptions, AlnSortOrder, AlnSource, Format, InvalidPolicy, OutFormat, OutputOptions,
    OutputSpec, ReaderOptions, WriterOptions, open_aln_source_from, open_output,
};
use gloci_ops::{AlnSortOptions, check_alignments, sort_alignments};
use std::io::Cursor;
use std::num::NonZeroUsize;
use std::path::Path;

fn source(name: &str, text: Vec<u8>) -> Box<dyn AlnSource> {
    let opts = ReaderOptions {
        accept: Accept {
            alignments: Some(AlignOptions {
                full_records: true,
                ..AlignOptions::default()
            }),
            ..Accept::default()
        },
        chrom_map: None,
        threads: NonZeroUsize::MIN,
        format: None,
        vcf_span: VcfSpan::Ref,
        invalid: InvalidPolicy::Fail,
    };
    open_aln_source_from(name, Box::new(Cursor::new(text)), &opts).unwrap()
}

const HEADER: &str = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:c1\tLN:1000000\n@SQ\tSN:c2\tLN:1000000\n\
@SQ\tSN:c3\tLN:1000000\n";

/// `n` pseudo-random records: few positions (so many ties), both strands, unplaced reads, names
/// with numbers of different lengths, every flag rank, and sequence `len` long.
fn sam_text(n: usize, len: usize) -> Vec<u8> {
    let mut out = HEADER.as_bytes().to_vec();
    let mut x = 0x2545_F491_4F6C_DD1Du64;
    let mut next = move || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    };
    let seq = "ACGT".repeat(len / 4 + 1);
    let seq = &seq[..len];
    let qual = "I".repeat(len);
    let flags = [
        0u16, 16, 0x41, 0x81, 0x101, 0x141, 0x801, 0x841, 4, 20, 0x61, 0x91,
    ];
    for _ in 0..n {
        let (r, f) = (next(), flags[(next() % 12) as usize]);
        let unplaced = f & 4 != 0 && r % 3 == 0;
        let (ref_name, pos) = if unplaced {
            ("*".to_string(), 0)
        } else {
            (format!("c{}", 1 + r % 3), 1 + (r >> 8) % 300)
        };
        let cigar = if f & 4 != 0 {
            "*".to_string()
        } else {
            format!("{len}M")
        };
        let name = format!(
            "r{}{}",
            (r >> 20) % 5000,
            ["", "a", "b7"][(r >> 40) as usize % 3]
        );
        out.extend_from_slice(
            format!(
                "{name}\t{f}\t{ref_name}\t{pos}\t30\t{cigar}\t*\t0\t0\t{seq}\t{qual}\tNM:i:{}\n",
                r % 7
            )
            .as_bytes(),
        );
    }
    out
}

/// A budget that spills several times for 60,000 small records (about 9 MiB of records).
const MEM: u64 = 2 << 20;

fn sort_to(
    dir: &Path,
    text: &[u8],
    order: AlnSortOrder,
    max_mem: u64,
    threads: usize,
    out_name: &str,
) -> (Vec<u8>, u32) {
    let threads = NonZeroUsize::new(threads).unwrap();
    let path = dir.join(out_name);
    let spec = OutputSpec {
        path: Some(path.clone()),
        format: Some(OutFormat::Bam),
        bgzf: false,
        direct: false,
    };
    let out = open_output(
        &spec,
        OutFormat::Bam,
        &OutputOptions {
            threads,
            level: 1,
            write_index: false,
        },
    )
    .unwrap();
    let o = AlnSortOptions {
        order,
        max_mem,
        tmp_dir: dir.to_path_buf(),
        threads,
    };
    let w = WriterOptions {
        write_header: true,
        command_line: None,
        version: "test",
    };
    let stats =
        sort_alignments(source("in.sam", text.to_vec()), &o, out, OutFormat::Bam, &w).unwrap();
    // The spill directory is gone.
    let left: Vec<_> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|n| n.starts_with("gloci-sort-"))
        .collect();
    assert!(left.is_empty(), "{left:?}");
    (std::fs::read(path).unwrap(), stats.spills)
}

#[test]
fn aln_sort_spills_equal_in_memory() {
    let dir = tempfile::tempdir().unwrap();
    let text = sam_text(60_000, 40);
    for order in [
        AlnSortOrder::Coordinate,
        AlnSortOrder::Name,
        AlnSortOrder::NameBytes,
    ] {
        let (want, spills) = sort_to(dir.path(), &text, order, 2 << 30, 1, "mem.bam");
        assert_eq!(spills, 0);
        for t in [1, 2, 8] {
            let (got, spills) = sort_to(dir.path(), &text, order, MEM, t, "spill.bam");
            assert!(spills >= 4, "{order:?} -t {t}: {spills} spills");
            assert!(
                got == want,
                "{order:?} -t {t}: output differs from the in-memory sort"
            );
        }
    }
}

#[test]
fn aln_sort_huge_record() {
    let dir = tempfile::tempdir().unwrap();
    let mut text = sam_text(2_000, 40);
    // A record of more than 40 MiB (a Z tag), placed in the middle of the sort order.
    let huge = format!(
        "huge\t0\tc2\t150\t30\t40M\t*\t0\t0\t{}\t{}\tXX:Z:{}\n",
        "ACGT".repeat(10),
        "I".repeat(40),
        "A".repeat(41 << 20)
    );
    text.extend_from_slice(huge.as_bytes());
    text.extend_from_slice(&sam_text(2_000, 40)[HEADER.len()..]);
    for order in [AlnSortOrder::Coordinate, AlnSortOrder::Name] {
        let (want, _) = sort_to(dir.path(), &text, order, 2 << 30, 1, "mem.bam");
        let (got, spills) = sort_to(dir.path(), &text, order, 16 << 20, 2, "spill.bam");
        assert!(spills >= 1);
        assert!(got == want, "{order:?}: the huge record changed the order");
    }
}

fn header_of(bam: &[u8]) -> String {
    // The decompressed header text, through gloci's own reader.
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("x.bam");
    std::fs::write(&p, bam).unwrap();
    let src = gloci_io::open_aln_source(
        p.to_str().unwrap(),
        &ReaderOptions {
            accept: Accept {
                alignments: Some(AlignOptions::default()),
                ..Accept::default()
            },
            chrom_map: None,
            threads: NonZeroUsize::MIN,
            format: Some(Format::Bam),
            vcf_span: VcfSpan::Ref,
            invalid: InvalidPolicy::Fail,
        },
    )
    .unwrap();
    String::from_utf8(src.header().text.clone()).unwrap()
}

#[test]
fn aln_sort_header() {
    let dir = tempfile::tempdir().unwrap();
    let sq = "@SQ\tSN:c1\tLN:10\n";
    let rec = "r\t0\tc1\t1\t30\t1M\t*\t0\t0\tA\tI\n";
    let cases: [(&str, AlnSortOrder, &str); 7] = [
        // SO replaced where it stands, GO and SS removed
        (
            "@HD\tVN:1.6\tSO:unsorted\tGO:query\tSS:x\n",
            AlnSortOrder::Coordinate,
            "@HD\tVN:1.6\tSO:coordinate\n",
        ),
        (
            "@HD\tVN:1.6\tSO:queryname\tSS:queryname:natural\n",
            AlnSortOrder::Coordinate,
            "@HD\tVN:1.6\tSO:coordinate\n",
        ),
        // SO appended when absent, GO removed
        (
            "@HD\tVN:1.6\tGO:reference\n",
            AlnSortOrder::Coordinate,
            "@HD\tVN:1.6\tSO:coordinate\n",
        ),
        (
            "@HD\tVN:1.6\tSO:coordinate\n",
            AlnSortOrder::Name,
            "@HD\tVN:1.6\tSO:queryname\tSS:queryname:natural\n",
        ),
        (
            "@HD\tVN:1.6\tSS:coordinate:minhash\tSO:coordinate\n",
            AlnSortOrder::NameBytes,
            "@HD\tVN:1.6\tSS:queryname:lexicographical\tSO:queryname\n",
        ),
        // no @HD: one first
        ("", AlnSortOrder::Coordinate, "@HD\tVN:1.6\tSO:coordinate\n"),
        (
            "",
            AlnSortOrder::NameBytes,
            "@HD\tVN:1.6\tSO:queryname\tSS:queryname:lexicographical\n",
        ),
    ];
    for (hd, order, want) in cases {
        let text = format!("{hd}{sq}{rec}").into_bytes();
        let (bam, _) = sort_to(dir.path(), &text, order, 2 << 30, 1, "h.bam");
        assert_eq!(header_of(&bam), format!("{want}{sq}"), "{hd:?} {order:?}");
    }
}

#[test]
fn aln_sort_check() {
    let sorted =
        format!("{HEADER}a\t0\tc1\t5\t30\t1M\t*\t0\t0\tA\tI\nb\t0\tc1\t5\t30\t1M\t*\t0\t0\tA\tI\n");
    assert!(
        check_alignments(
            source("s.sam", sorted.into_bytes()),
            AlnSortOrder::Coordinate
        )
        .unwrap()
    );
    let bad =
        format!("{HEADER}a\t0\tc1\t5\t30\t1M\t*\t0\t0\tA\tI\nb\t0\tc1\t4\t30\t1M\t*\t0\t0\tA\tI\n");
    assert!(
        !check_alignments(
            source("s.sam", bad.clone().into_bytes()),
            AlnSortOrder::Coordinate
        )
        .unwrap()
    );
    // by name the same file is in order
    assert!(check_alignments(source("s.sam", bad.into_bytes()), AlnSortOrder::Name).unwrap());
}
