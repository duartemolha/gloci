//! P3c's shared foundation: the record accessors of every layout, FLAGS parsing and names, the
//! regions source, and the CRAM prelude.

use gloci_core::{Interval, RegionSet, VcfSpan};
use gloci_io::*;
use gloci_testkit::fixtures::sam_to_bam;
use std::io::Write;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

fn aln(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/aln")
        .join(name)
}

fn ropts() -> ReaderOptions {
    ReaderOptions {
        accept: Accept {
            alignments: Some(AlignOptions::default()),
            ..Accept::default()
        },
        chrom_map: None,
        threads: NonZeroUsize::MIN,
        format: None,
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    }
}

static WARNINGS: Mutex<Vec<String>> = Mutex::new(Vec::new());

struct Capture;
impl log::Log for Capture {
    fn enabled(&self, _: &log::Metadata) -> bool {
        true
    }
    fn log(&self, r: &log::Record) {
        WARNINGS.lock().unwrap().push(r.args().to_string());
    }
    fn flush(&self) {}
}

fn scan_warnings() -> usize {
    WARNINGS
        .lock()
        .unwrap()
        .iter()
        .filter(|w| w.contains("scanning the whole file"))
        .count()
}

#[test]
fn parse_flags_samtools_rules() {
    for (s, v) in [
        ("0x704", 0x704),
        ("1796", 1796),
        ("secondary,dup", 0x500),
        ("Supplementary", 0x800),
        ("0", 0),
        ("4096", 0x1000),
        ("65535", 0xFFFF),
    ] {
        assert_eq!(parse_flags(s), Ok(v), "{s}");
    }
    for bad in ["65536", "12abc", "-1", "0x", ""] {
        assert!(parse_flags(bad).is_err(), "{bad:?}");
    }
    let e = parse_flags("010").unwrap_err();
    assert!(e.contains("leading zero") && e.contains("0x8"), "{e}");
    assert!(parse_flags("bogus").unwrap_err().contains("PAIRED"));
}

#[test]
fn flag_names_table() {
    assert_eq!(flag_names(0x704), "UNMAP,SECONDARY,QCFAIL,DUP");
    assert_eq!(flag_names(0), "");
    assert_eq!(flag_names(0xFFF).split(',').count(), 12);
    assert_eq!(flag_names(0x1000), "");
    assert_eq!(SAM_FLAGS.len(), 12);
}

type Regions = Vec<(Vec<u8>, Interval)>;
type Fields = (Option<usize>, Option<u64>, usize, Vec<u8>);

fn fields(a: &AlnRef<'_>) -> Fields {
    (
        a.mate_ref_id(),
        a.mate_pos(),
        a.seq_len(),
        a.quals().to_vec(),
    )
}

fn all(path: &Path) -> Vec<Fields> {
    let mut src = open_aln_source(path.to_str().unwrap(), &ropts()).unwrap();
    let (layout, h) = (src.layout(), src.header().clone());
    let mut buf = Vec::new();
    let mut v = Vec::new();
    while src.next_record(&mut buf).unwrap().is_some() {
        v.push(fields(&aln_ref(layout, &buf, &h).unwrap()));
        buf.clear();
    }
    v
}

#[test]
fn aln_accessors_bam_equals_sam() {
    let sam = all(&aln("flags.sam"));
    assert_eq!(sam, all(&aln("flags.bam")));
    assert_eq!(sam.len(), 14);
    // r_pair 99: the mate is on chr1 at 400.
    assert_eq!(sam[4].0, Some(0));
    assert_eq!(sam[4].1, Some(400));
    // A read with no mate has none.
    assert_eq!((sam[0].0, sam[0].1), (None, None));
    assert_eq!(sam[0].2, 10);
    assert_eq!(sam[0].3, vec![b'I' - 33; 10]);
}

#[test]
fn aln_accessors_missing_seq_and_qual() {
    let sam = "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:1000\n@SQ\tSN:c2\tLN:1000\n\
        a\t0\tc1\t10\t60\t4M\t*\t0\t0\tACGT\t*\n\
        b\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n\
        c\t65\tc1\t20\t60\t4M\tc2\t30\t0\tACGT\t!#5I\n\
        d\t65\tc1\t20\t60\t4M\t=\t30\t0\tACGT\t!#5I\n\
        e\t0\tc1\t10\t60\t*\t*\t0\t0\t*\t*\n";
    let dir = tempfile::tempdir().unwrap();
    let sam_path = dir.path().join("t.sam");
    std::fs::write(&sam_path, sam).unwrap();
    let bam_path = dir.path().join("t.bam");
    std::fs::File::create(&bam_path)
        .unwrap()
        .write_all(&sam_to_bam(sam))
        .unwrap();
    let from_sam = all(&sam_path);
    assert_eq!(from_sam, all(&bam_path));
    assert_eq!(from_sam[0].3, vec![0xFF; 4], "QUAL * is 0xFF per base");
    assert_eq!(from_sam[1], (None, None, 0, vec![]), "SEQ * has no quals");
    assert_eq!(from_sam[2], (Some(1), Some(29), 4, vec![0, 2, 20, 40]));
    assert_eq!(from_sam[3].0, Some(0), "= is the read's reference");
    assert_eq!(from_sam[4].2, 0);
}

fn collect(mut src: Box<dyn AlnSource>) -> Vec<Vec<u8>> {
    let mut buf = Vec::new();
    let mut v = Vec::new();
    while src.next_record(&mut buf).unwrap().is_some() {
        v.push(std::mem::take(&mut buf));
    }
    v
}

fn iv(start: u64, end: u64) -> Interval {
    Interval { start, end }
}

/// What `open_aln_source` gives for `path`, filtered by §4.3 with the D116 span.
fn filtered(path: &str, regions: &[(Vec<u8>, Interval)]) -> Vec<Vec<u8>> {
    let src = open_aln_source(path, &ropts()).unwrap();
    let (layout, h) = (src.layout(), src.header().clone());
    let sets: Vec<Option<RegionSet>> = h
        .refs
        .iter()
        .map(|(n, _)| {
            let ivs: Vec<Interval> = regions
                .iter()
                .filter(|(c, _)| c == n)
                .map(|(_, i)| *i)
                .collect();
            (!ivs.is_empty()).then(|| RegionSet::new(ivs))
        })
        .collect();
    collect(src)
        .into_iter()
        .filter(|raw| {
            aln_ref(layout, raw, &h)
                .ok()
                .and_then(|a| alignment_query_span(&a))
                .is_some_and(|(rid, span)| {
                    sets[rid].as_ref().is_some_and(|s| s.matches(span, false))
                })
        })
        .collect()
}

#[test]
fn regions_source_equals_filtered_scan() {
    let _ = log::set_logger(&Capture);
    log::set_max_level(log::LevelFilter::Warn);
    let dir = tempfile::tempdir().unwrap();
    let c = |s: &str| s.as_bytes().to_vec();
    let cases: Vec<(&str, Regions)> = vec![
        // overlapping, and regions that only share an edge with reads.
        ("flags", vec![(c("chr1"), iv(105, 205))]),
        ("flags", vec![(c("chr1"), iv(110, 200))]),
        (
            "flags",
            vec![(c("chr1"), iv(0, 100)), (c("chr1"), iv(110, 111))],
        ),
        // a point inside a read, and on its edge.
        ("flags", vec![(c("chr1"), iv(105, 105))]),
        (
            "flags",
            vec![(c("chr1"), iv(100, 100)), (c("chr1"), iv(110, 110))],
        ),
        // a placed unmapped read, which spans [700, 701).
        ("flags", vec![(c("chr1"), iv(700, 701))]),
        (
            "flags",
            vec![(c("chr1"), iv(699, 700)), (c("chr1"), iv(701, 702))],
        ),
        // regions on both references, repeated, and one the header lacks.
        (
            "flags",
            vec![
                (c("chr2"), iv(40, 70)),
                (c("chr1"), iv(95, 700)),
                (c("chr1"), iv(100, 400)),
                (c("chrZ"), iv(0, 10)),
            ],
        ),
        (
            "order",
            vec![(c("chr1"), iv(0, 1_000_000)), (c("chr2"), iv(0, 1_000_000))],
        ),
        ("order", vec![(c("chr2"), iv(100, 200))]),
        ("order", vec![]),
    ];
    for set in ["flags", "order"] {
        let bam = aln(&format!("{set}.bam"));
        let indexed = dir.path().join(format!("{set}-ix.bam"));
        std::fs::copy(&bam, &indexed).unwrap();
        std::fs::copy(
            aln(&format!("{set}.bam.bai")),
            dir.path().join(format!("{set}-ix.bam.bai")),
        )
        .unwrap();
        let plain = dir.path().join(format!("{set}-plain.bam"));
        std::fs::copy(&bam, &plain).unwrap();
        for (s, regions) in cases.iter().filter(|(s, _)| *s == set) {
            let want = filtered(bam.to_str().unwrap(), regions);
            let via_index = collect(
                open_aln_source_regions(indexed.to_str().unwrap(), regions, &ropts()).unwrap(),
            );
            assert_eq!(via_index, want, "{s} {regions:?}: through the index");
            let before = scan_warnings();
            let scanned = collect(
                open_aln_source_regions(plain.to_str().unwrap(), regions, &ropts()).unwrap(),
            );
            assert_eq!(scanned, want, "{s} {regions:?}: scanned");
            assert_eq!(scan_warnings(), before + 1, "one warning for the scan");
        }
    }
    // The cases are not all empty: an overlap, and the placed unmapped read, give records.
    assert_eq!(
        filtered(aln("flags.bam").to_str().unwrap(), &cases[0].1).len(),
        3
    );
    assert_eq!(
        filtered(aln("flags.bam").to_str().unwrap(), &cases[5].1).len(),
        1
    );
    // SAM is scanned too.
    let sam = aln("flags.sam");
    let regions = &cases[0].1;
    assert_eq!(
        collect(open_aln_source_regions(sam.to_str().unwrap(), regions, &ropts()).unwrap()).len(),
        3
    );
}

#[test]
fn cram_prelude_reads_header() {
    for name in ["flags.cram", "flags.v31.cram"] {
        let path = aln(name);
        let p = read_cram_prelude(&path).unwrap();
        let (major, minor) = if name.contains("v31") { (3, 1) } else { (3, 0) };
        assert_eq!((p.major, p.minor), (major, minor), "{name}");
        // The header text is what the file's own header reader finds.
        let f = std::fs::File::open(&path).unwrap();
        let h = read_cram_header(&mut std::io::BufReader::new(f), name).unwrap();
        assert_eq!(p.header_text, h.text, "{name}");
        assert!(p.header_text.starts_with(b"@HD"));
        // data_start is where noodles stands after the header container, and what follows is
        // data containers, then the end-of-file container.
        let mut r = noodles_cram::io::Reader::new(std::fs::File::open(&path).unwrap());
        r.read_header().unwrap();
        assert_eq!(p.data_start, r.position().unwrap(), "{name}");
        let bytes = std::fs::read(&path).unwrap();
        assert!(p.data_start + 38 < bytes.len() as u64, "{name}");
        assert!(cram_eof_present(&path).unwrap(), "{name}");
        let dir = tempfile::tempdir().unwrap();
        let cut = dir.path().join("cut.cram");
        std::fs::write(&cut, &bytes[..bytes.len() - 38]).unwrap();
        assert!(!cram_eof_present(&cut).unwrap(), "{name}: cut");
        let cut1 = dir.path().join("cut1.cram");
        std::fs::write(&cut1, &bytes[..bytes.len() - 1]).unwrap();
        assert!(!cram_eof_present(&cut1).unwrap());
    }
    // Not CRAM: an error, not a panic.
    assert!(read_cram_prelude(&aln("flags.sam")).is_err());
    assert!(read_cram_prelude(&aln("nonexistent.cram")).is_err());
}

#[test]
fn sam_aln_refuses_a_qual_that_does_not_fit_its_seq() {
    let h = AlnHeader::parse_sam(b"@HD\tVN:1.6\n@SQ\tSN:c1\tLN:1000\n".to_vec()).unwrap();
    let line = |seq: &str, qual: &str| format!("a\t0\tc1\t10\t60\t*\t*\t0\t0\t{seq}\t{qual}");
    let e = SamAln::parse(line("ACGT", "II").as_bytes(), &h).unwrap_err();
    assert!(e.contains("2 characters") && e.contains("sequence has 4"), "{e}");
    let e = SamAln::parse(line("*", "II").as_bytes(), &h).unwrap_err();
    assert!(e.contains("sequence has 0"), "{e}");
    let e = SamAln::parse(line("ACGT", "II I").as_bytes(), &h).unwrap_err();
    assert!(e.contains("outside ! to ~"), "{e}");
    // The accessors of a line that parses always agree with each other.
    for (seq, qual) in [("ACGT", "!#5~"), ("ACGT", "*"), ("*", "*")] {
        let l = line(seq, qual);
        let a = SamAln::parse(l.as_bytes(), &h).unwrap();
        assert_eq!(a.quals().len(), a.seq_len(), "{seq} {qual}");
    }
}
