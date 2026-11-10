//! Indexed queries and index builds refuse non-interval formats without panicking.
use gloci_core::{Interval, VcfSpan};
use gloci_io::*;
use gloci_testkit::fixtures;
use std::num::NonZeroUsize;

fn opts(accept: Accept) -> ReaderOptions {
    ReaderOptions {
        accept,
        chrom_map: None,
        threads: NonZeroUsize::MIN,
        format: None,
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    }
}

/// A BGZF `x.sam.gz` with a real tabix index next to it (built from BED text).
fn sam_gz_with_tbi(dir: &std::path::Path) -> std::path::PathBuf {
    let bed = fixtures::write(dir, "a.bed.gz", &fixtures::bgzip(b"chr1\t1\t2\n"));
    write_index(&bed, &bed, Format::Bed, NonZeroUsize::MIN).unwrap();
    let sam = fixtures::write(
        dir,
        "x.sam.gz",
        &fixtures::bgzip(b"@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\t*\n"),
    );
    std::fs::copy(dir.join("a.bed.gz.tbi"), dir.join("x.sam.gz.tbi")).unwrap();
    sam
}

#[test]
fn indexed_query_refuses_alignment_and_convert_only_formats() {
    let dir = tempfile::tempdir().unwrap();
    let sam = sam_gz_with_tbi(dir.path());
    let regions = [(b"chr1".to_vec(), Interval { start: 0, end: 5 })];
    for f in [
        Format::Sam,
        Format::Cram,
        Format::Psl,
        Format::Rmsk,
        Format::Fasta,
    ] {
        let o = opts(Accept {
            alignments: Some(AlignOptions::default()),
            convert_only: true,
            alignments_planned: false,
        });
        let e = IndexedFile::open(&sam, f, Header::None, &o).err().unwrap();
        assert!(e.to_string().contains("tabix or CSI index"), "{f:?}: {e}");
    }
    // Through the path commands use, with a command that accepts alignments.
    let accept = Accept {
        alignments: Some(AlignOptions::default()),
        ..Accept::default()
    };
    assert!(open_region_reader(&sam, &regions, &opts(accept), false).is_err());
}

#[test]
fn index_build_refuses_alignments_and_psl() {
    let dir = tempfile::tempdir().unwrap();
    let p = fixtures::write(dir.path(), "a.bed.gz", &fixtures::bgzip(b"chr1\t1\t2\n"));
    let e = write_index(&p, &p, Format::Cram, NonZeroUsize::MIN).unwrap_err();
    assert!(
        e.to_string().contains("not available in this build yet"),
        "{e}"
    );
    let e = write_index(&p, &p, Format::Psl, NonZeroUsize::MIN).unwrap_err();
    assert!(e.to_string().contains("gloci convert"), "{e}");
}
