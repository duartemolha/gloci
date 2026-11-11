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
    assert!(
        e.contains("2 characters") && e.contains("sequence has 4"),
        "{e}"
    );
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

/// A record's name, flag and reference span (chromosome, start, end), if it has one.
type HandSpan = (String, u16, Option<(String, u64, u64)>);

/// `(name, flag, span)` of each record of flags.sam, worked out from the text by hand: the span
/// of a mapped read is its POS and the reference bases of its CIGAR (M, D, N, =, X; at least
/// one), that of a placed unmapped read is the one base at POS (D101, D116); a read with no
/// reference or position has none.
fn hand_spans() -> Vec<HandSpan> {
    let text = std::fs::read_to_string(aln("flags.sam")).unwrap();
    text.lines()
        .filter(|l| !l.starts_with('@'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let flag: u16 = f[1].parse().unwrap();
            let pos: u64 = f[3].parse().unwrap();
            let span = if f[2] == "*" || pos == 0 {
                None
            } else if flag & 4 != 0 {
                Some((f[2].to_string(), pos - 1, pos))
            } else {
                let mut covered = 0u64;
                let mut n = 0u64;
                for c in f[5].chars() {
                    match c {
                        '0'..='9' => n = n * 10 + u64::from(c.to_digit(10).unwrap()),
                        'M' | 'D' | 'N' | '=' | 'X' => {
                            covered += n;
                            n = 0;
                        }
                        _ => n = 0,
                    }
                }
                Some((f[2].to_string(), pos - 1, pos - 1 + covered.max(1)))
            };
            (f[0].to_string(), flag, span)
        })
        .collect()
}

/// §4.3 by hand: a point region overlaps a span that holds it strictly inside; two intervals
/// overlap when they share a base, so bookended ones do not.
fn hand_overlap(span: (u64, u64), region: (u64, u64)) -> bool {
    if region.0 == region.1 {
        span.0 < region.0 && region.0 < span.1
    } else {
        span.0 < region.1 && region.0 < span.1
    }
}

fn named(src: Box<dyn AlnSource>) -> Vec<(String, u16)> {
    let (layout, h) = (src.layout(), src.header().clone());
    collect(src)
        .iter()
        .map(|raw| {
            let a = aln_ref(layout, raw, &h).unwrap();
            (String::from_utf8(a.name().to_vec()).unwrap(), a.flag())
        })
        .collect()
}

#[test]
fn regions_source_matches_hand_oracle_on_every_index() {
    let dir = tempfile::tempdir().unwrap();
    let copy = |from: &str, to: &str| {
        std::fs::copy(aln(from), dir.path().join(to)).unwrap();
    };
    // BAM: .bai, .csi (made by gloci), none. CRAM: .crai, none. SAM: scanned.
    copy("flags.bam", "bai.bam");
    copy("flags.bam.bai", "bai.bam.bai");
    copy("flags.bam", "csi.bam");
    write_alignment_index(
        &dir.path().join("csi.bam"),
        &dir.path().join("csi.bam"),
        true,
        NonZeroUsize::MIN,
    )
    .unwrap();
    assert!(dir.path().join("csi.bam.csi").exists());
    copy("flags.bam", "none.bam");
    copy("flags.cram", "crai.cram");
    copy("flags.cram.crai", "crai.cram.crai");
    copy("flags.cram", "none.cram");
    copy("flags.sam", "none.sam");
    let fa = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs/fa/ref.fa");
    let o = || {
        let mut o = ropts();
        o.accept.alignments.as_mut().unwrap().reference = Some(fa.clone());
        o
    };
    let regions: Vec<(&str, u64, u64)> = vec![
        ("chr1", 105, 205), // overlaps three reads
        ("chr1", 110, 200), // bookended by r_plain (ends 110) and r_mapq0 (starts 200)
        ("chr1", 105, 105), // a point strictly inside two reads
        ("chr1", 100, 100), // a point on the start edge of two reads
        ("chr1", 110, 110), // a point on their end edge
        ("chr1", 700, 701), // the placed unmapped read
        ("chr1", 700, 700), // a point on its start
        ("chr1", 549, 551), // the soft-clipped read r_supp starts at 550
        ("chr2", 50, 51),   // two reads at once
        ("chr2", 0, 50),    // bookended by them
        ("chrZ", 0, 100),   // not in the header
    ];
    let spans = hand_spans();
    let mut non_empty = 0;
    for (chrom, s, e) in &regions {
        let want: Vec<(String, u16)> = spans
            .iter()
            .filter(|(_, _, sp)| {
                sp.as_ref()
                    .is_some_and(|(c, a, b)| c == chrom && hand_overlap((*a, *b), (*s, *e)))
            })
            .map(|(n, f, _)| (n.clone(), *f))
            .collect();
        non_empty += usize::from(!want.is_empty());
        let region = vec![(chrom.as_bytes().to_vec(), iv(*s, *e))];
        for file in [
            "bai.bam",
            "csi.bam",
            "none.bam",
            "crai.cram",
            "none.cram",
            "none.sam",
        ] {
            let p = dir.path().join(file);
            let got = named(open_aln_source_regions(p.to_str().unwrap(), &region, &o()).unwrap());
            assert_eq!(got, want, "{file} {chrom}:{s}-{e}");
        }
    }
    assert!(non_empty >= 5, "the regions must select something");
    // Known answers, so the oracle itself is checked.
    let one = |c: &str, s, e| -> Vec<String> {
        spans
            .iter()
            .filter(|(_, _, sp)| {
                sp.as_ref()
                    .is_some_and(|(cc, a, b)| cc == c && hand_overlap((*a, *b), (s, e)))
            })
            .map(|(n, _, _)| n.clone())
            .collect()
    };
    assert_eq!(one("chr1", 105, 205), ["r_plain", "r_rev", "r_mapq0"]);
    assert!(one("chr1", 110, 200).is_empty());
    assert_eq!(one("chr1", 700, 701), ["r_unmapped_placed"]);
    assert!(one("chr1", 700, 700).is_empty());
}

#[test]
fn cram_prelude_survives_damaged_files() {
    let dir = tempfile::tempdir().unwrap();
    let good = std::fs::read(aln("flags.cram")).unwrap();
    let start = read_cram_prelude(&aln("flags.cram")).unwrap().data_start as usize;
    let write = |name: &str, bytes: &[u8]| {
        let p = dir.path().join(name);
        std::fs::write(&p, bytes).unwrap();
        p
    };
    // Cut inside the file definition, and at several places inside the header container.
    for cut in [0, 3, 10, 26, 27, 30, 40, start / 2, start - 1] {
        let p = write("cut.cram", &good[..cut]);
        assert!(read_cram_prelude(&p).is_err(), "cut at {cut}");
    }
    // A header container that claims 2 GiB - 1 (or more) of data, in the container's length field
    // and in the header block's sizes: an error, quickly, without allocating that much.
    for at in [26usize, 26 + 4 + 5, 40, 44, 48, 52] {
        let mut bad = good[..start.min(good.len())].to_vec();
        for (i, b) in [0xff, 0xff, 0xff, 0x7f].iter().enumerate() {
            if at + i < bad.len() {
                bad[at + i] = *b;
            }
        }
        let p = write("huge.cram", &bad);
        let t = std::time::Instant::now();
        let r = read_cram_prelude(&p);
        assert!(r.is_err(), "length bytes at {at}: {r:?}");
        assert!(
            t.elapsed().as_secs() < 5,
            "length bytes at {at} took too long"
        );
    }
    // Garbage, an empty file, and a directory.
    let junk: Vec<u8> = (0..4096u32)
        .map(|i| (i.wrapping_mul(2654435761) >> 13) as u8)
        .collect();
    assert!(read_cram_prelude(&write("junk.cram", &junk)).is_err());
    let mut magic = b"CRAM\x03\x00".to_vec();
    magic.extend_from_slice(&junk);
    assert!(read_cram_prelude(&write("magic.cram", &magic)).is_err());
    assert!(read_cram_prelude(&write("empty.cram", b"")).is_err());
    assert!(read_cram_prelude(dir.path()).is_err());
    // The EOF check needs only the tail: junk and short files have none.
    assert!(!cram_eof_present(&write("junk.cram", &junk)).unwrap());
    assert!(!cram_eof_present(&write("short.cram", b"CRAM")).unwrap());
}
