//! The BAM reader for interval commands: CIGAR spans, unmapped reads and filters,
//! header order, malformed records, slim records and large headers (D101-D103, D115, D119).

use gloci_core::{ChromTable, Interval, OrderSource, Strand, VcfSpan};
use gloci_io::sorted::{OrderGuard, sort_hint};
use gloci_io::*;
use gloci_testkit::fixtures::{bgzip, sam_to_bam};
use std::io::{Cursor, Read};
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::sync::Arc;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/aln")
        .join(name)
}

fn read_fixture(name: &str) -> Vec<u8> {
    std::fs::read(fixture(name)).unwrap()
}

fn ropts(ao: AlignOptions, invalid: InvalidPolicy, map: Option<Arc<ChromMap>>) -> ReaderOptions {
    ReaderOptions {
        accept: Accept {
            alignments: Some(ao),
            ..Accept::default()
        },
        chrom_map: map,
        threads: NonZeroUsize::MIN,
        format: None,
        vcf_span: VcfSpan::Allele,
        invalid,
    }
}

fn open(bam: &[u8], o: &ReaderOptions) -> Result<Box<dyn RecordReader>, IoError> {
    open_reader_from("t.bam", Box::new(Cursor::new(bam.to_vec())), o)
}

/// One record as the reader returned it.
struct Got {
    chrom: String,
    span: Interval,
    strand: Strand,
    loc: u64,
    raw: Vec<u8>,
}

struct Run {
    recs: Vec<Got>,
    diag: Diagnostics,
    hdr: Arc<AlnHeader>,
    table: ChromTable,
}

fn run_with(bam: &[u8], o: &ReaderOptions) -> Result<Run, IoError> {
    let mut r = open(bam, o)?;
    let mut table = chrom_table_for(None, &[&*r], gloci_core::ChromOrder::Karyotype)?;
    let Header::Alignment(hdr) = r.header().clone() else {
        panic!("a BAM reader has an alignment header");
    };
    let mut recs = Vec::new();
    let mut batch = RecordBatch::default();
    loop {
        batch.clear();
        let more = r.read_batch(&mut table, &mut batch, 3)?;
        for i in 0..batch.len() {
            let m = batch.recs[i];
            recs.push(Got {
                chrom: String::from_utf8_lossy(table.name(m.chrom)).into_owned(),
                span: m.span,
                strand: m.strand,
                loc: m.loc,
                raw: batch.raw(i).to_vec(),
            });
        }
        if !more {
            break;
        }
    }
    Ok(Run {
        recs,
        diag: r.diagnostics().clone(),
        hdr,
        table,
    })
}

fn run(bam: &[u8], ao: AlignOptions) -> Run {
    run_with(bam, &ropts(ao, InvalidPolicy::Fail, None)).unwrap()
}

fn name_of(g: &Got, h: &AlnHeader) -> String {
    let a = aln_ref(Format::Bam, &g.raw, h).unwrap();
    String::from_utf8_lossy(a.name()).into_owned()
}

fn ops_of(g: &Got, h: &AlnHeader) -> Vec<(u8, u32)> {
    let a = aln_ref(Format::Bam, &g.raw, h).unwrap();
    let mut v = Vec::new();
    a.for_each_op(&mut |op, n| v.push((op, n)));
    v
}

fn sam_header(refs: &[(&str, u64)], so: &str) -> String {
    let mut s = format!("@HD\tVN:1.6\tSO:{so}\n");
    for (n, l) in refs {
        s.push_str(&format!("@SQ\tSN:{n}\tLN:{l}\n"));
    }
    s
}

fn sam_read(name: &str, flag: u16, chrom: &str, pos1: u64, cigar: &str, seq: usize) -> String {
    let (s, q) = if seq == 0 {
        ("*".to_string(), "*".to_string())
    } else {
        ("A".repeat(seq), "I".repeat(seq))
    };
    format!("{name}\t{flag}\t{chrom}\t{pos1}\t60\t{cigar}\t*\t0\t0\t{s}\t{q}\n")
}

#[test]
fn bam_spans_table() {
    // (CIGAR, flag, query length, expected span, strand)
    let cases: &[(&str, u16, usize, (u64, u64))] = &[
        ("10M", 0, 10, (100, 110)),
        ("5S10M5S", 0, 20, (100, 110)),
        ("3H10M2H", 0, 10, (100, 110)),
        ("5M3I5M", 0, 13, (100, 110)),
        ("4=1X5=", 0, 10, (100, 110)),
        ("5M100N5M", 0, 10, (100, 210)),
        ("5M2D5M", 0, 10, (100, 112)),
        ("*", 0, 10, (100, 101)),
        ("10S5I", 0, 15, (100, 101)),
        ("10M", 16, 10, (100, 110)),
    ];
    let mut sam = sam_header(&[("chr1", 1_000_000), ("chrM", 16_569)], "coordinate");
    for (i, (cigar, flag, q, _)) in cases.iter().enumerate() {
        sam.push_str(&sam_read(&format!("r{i}"), *flag, "chr1", 101, cigar, *q));
    }
    sam.push_str(&sam_read("wrap", 0, "chrM", 16_561, "20M", 20));
    let bam = sam_to_bam(&sam);
    let r = run(&bam, AlignOptions::default());
    assert_eq!(r.recs.len(), cases.len() + 1);
    for (g, (cigar, flag, _, (s, e))) in r.recs.iter().zip(cases) {
        assert_eq!(
            g.span,
            Interval { start: *s, end: *e },
            "{cigar} flag {flag}"
        );
        let want = if flag & 16 != 0 {
            Strand::Minus
        } else {
            Strand::Plus
        };
        assert_eq!(g.strand, want, "{cigar} flag {flag}");
        let a = aln_ref(Format::Bam, &g.raw, &r.hdr).unwrap();
        let (iv, clipped) = a.span(1_000_000).unwrap();
        assert_eq!((iv, clipped), (g.span, false), "{cigar}");
    }
    let wrap = r.recs.last().unwrap();
    assert_eq!(
        wrap.span,
        Interval {
            start: 16_560,
            end: 16_569
        }
    );
    let a = aln_ref(Format::Bam, &wrap.raw, &r.hdr).unwrap();
    assert_eq!(a.span(16_569), Some((wrap.span, true)));
    // An unclipped read of the same chromosome reports no clipping.
    assert_eq!(a.span(1_000_000).map(|s| s.1), Some(false));
}

/// A raw BAM record (without `block_size`).
#[allow(clippy::too_many_arguments)]
fn raw_record(
    ref_id: i32,
    pos: i32,
    name: &[u8],
    flag: u16,
    cigar: &[u32],
    l_seq: u32,
    tags: &[u8],
) -> Vec<u8> {
    let mut v = Vec::new();
    v.extend_from_slice(&ref_id.to_le_bytes());
    v.extend_from_slice(&pos.to_le_bytes());
    v.push(name.len() as u8 + 1);
    v.push(60);
    v.extend_from_slice(&0u16.to_le_bytes());
    v.extend_from_slice(&(cigar.len() as u16).to_le_bytes());
    v.extend_from_slice(&flag.to_le_bytes());
    v.extend_from_slice(&l_seq.to_le_bytes());
    v.extend_from_slice(&(-1i32).to_le_bytes());
    v.extend_from_slice(&(-1i32).to_le_bytes());
    v.extend_from_slice(&0i32.to_le_bytes());
    v.extend_from_slice(name);
    v.push(0);
    for c in cigar {
        v.extend_from_slice(&c.to_le_bytes());
    }
    v.extend(std::iter::repeat_n(0x11u8, (l_seq as usize).div_ceil(2)));
    v.extend(std::iter::repeat_n(30u8, l_seq as usize));
    v.extend_from_slice(tags);
    v
}

fn header(refs: &[(&str, u64)], so: &str) -> AlnHeader {
    AlnHeader {
        text: sam_header(refs, so).into_bytes(),
        refs: refs
            .iter()
            .map(|(n, l)| (n.as_bytes().to_vec(), *l))
            .collect(),
        sort_order: SortOrder::Coordinate,
    }
}

/// A BGZF BAM of `h` with `records`, each framed with its `block_size`.
fn bam_file(h: &AlnHeader, records: &[Vec<u8>]) -> Vec<u8> {
    let mut v = Vec::new();
    h.write_bam(&mut v).unwrap();
    for r in records {
        v.extend_from_slice(&(r.len() as u32).to_le_bytes());
        v.extend_from_slice(r);
    }
    bgzip(&v)
}

const M: u32 = 0;
const I: u32 = 1;
const N: u32 = 3;
const S: u32 = 4;

#[test]
fn bam_long_cigar_cg_tag() {
    // 70,000 operations (1M1I, 35,000 times): too many for the 16-bit count, so the
    // record holds a placeholder (70000S 35000N) and the CG:B:I tag holds the CIGAR.
    let h = header(&[("chr1", 1_000_000)], "coordinate");
    let mut tag = b"CGBI".to_vec();
    tag.extend_from_slice(&70_000u32.to_le_bytes());
    for _ in 0..35_000 {
        tag.extend_from_slice(&(1u32 << 4 | M).to_le_bytes());
        tag.extend_from_slice(&(1u32 << 4 | I).to_le_bytes());
    }
    let rec = raw_record(
        0,
        99,
        b"long",
        0,
        &[70_000 << 4 | S, 35_000 << 4 | N],
        70_000,
        &tag,
    );
    let bam = bam_file(&h, &[rec]);
    let r = run(&bam, AlignOptions::default());
    assert_eq!(r.recs.len(), 1);
    assert_eq!(
        r.recs[0].span,
        Interval {
            start: 99,
            end: 35_099
        }
    );
    assert_eq!(ops_of(&r.recs[0], &r.hdr).len(), 70_000);
    // Slim keeps the CG tag, so the span and operations survive.
    let full = AlignOptions {
        full_records: true,
        ..AlignOptions::default()
    };
    let f = run(&bam, full);
    assert!(f.recs[0].raw.len() > r.recs[0].raw.len());
    assert_eq!(ops_of(&f.recs[0], &f.hdr), ops_of(&r.recs[0], &r.hdr));
}

#[test]
fn bam_unmapped_and_filters() {
    let bam = read_fixture("flags.bam");
    let names = |r: &Run| -> Vec<String> { r.recs.iter().map(|g| name_of(g, &r.hdr)).collect() };

    let r = run(&bam, AlignOptions::default());
    // Every mapped read is kept, secondary, supplementary, duplicate and QC-fail included.
    assert_eq!(
        names(&r),
        [
            "r_plain", "r_rev", "r_mapq0", "r_mapq10", "r_pair", "r_pair", "r_sec", "r_supp",
            "r_dup", "r_qcfail", "read9", "read10"
        ]
    );
    // The two unmapped reads (one placed, one not) are skipped and counted.
    assert_eq!((r.diag.unmapped, r.diag.filtered), (2, 0));

    let q = run(
        &bam,
        AlignOptions {
            min_mapq: 20,
            ..AlignOptions::default()
        },
    );
    let kept = names(&q);
    assert_eq!(kept.len(), 10);
    assert!(!kept.contains(&"r_mapq0".to_string()) && !kept.contains(&"r_mapq10".to_string()));
    assert_eq!((q.diag.unmapped, q.diag.filtered), (2, 2));

    let x = run(
        &bam,
        AlignOptions {
            exclude_flags: 0x900,
            ..AlignOptions::default()
        },
    );
    let kept = names(&x);
    assert_eq!(kept.len(), 10);
    assert!(!kept.contains(&"r_sec".to_string()) && !kept.contains(&"r_supp".to_string()));
    assert_eq!(x.diag.filtered, 2);

    let p = run(
        &bam,
        AlignOptions {
            require_flags: 0x2,
            ..AlignOptions::default()
        },
    );
    assert_eq!(names(&p), ["r_pair", "r_pair"]);
    assert_eq!(p.diag.filtered, 10);
    // The two mates are one chromosome's records 5 and 6 of the file.
    assert_eq!(p.recs.iter().map(|g| g.loc).collect::<Vec<_>>(), [5, 6]);
}

#[test]
fn bam_header_is_an_order_source() {
    let bam = read_fixture("order.bam");
    let o = ropts(AlignOptions::default(), InvalidPolicy::Fail, None);
    let mut r = open_reader_from("order.bam", Box::new(Cursor::new(bam.clone())), &o).unwrap();
    let mut table = chrom_table_for(None, &[&*r], gloci_core::ChromOrder::Karyotype).unwrap();
    assert_eq!(
        table.order_source(),
        &OrderSource::AlignmentHeader("order.bam".into())
    );
    let mut batch = RecordBatch::default();
    let mut guard = OrderGuard::new("order.bam", Format::Bam);
    let mut n = 0;
    loop {
        batch.clear();
        let more = r.read_batch(&mut table, &mut batch, 4).unwrap();
        guard.check(&table, &batch).unwrap();
        n += batch.len();
        if !more {
            break;
        }
    }
    assert_eq!(n, 15);
    let order: Vec<_> = (0..5)
        .map(|i| String::from_utf8_lossy(table.name(gloci_core::ChromId(i))).into_owned())
        .collect();
    assert_eq!(order, ["chrM", "chr1", "chr2", "chr10", "chrX"]);
    // A BED file sorted lexicographically is told to sort by this header.
    assert_eq!(
        sort_hint("b.bed", table.order_source()),
        "run: gloci sort b.bed -g order.bam"
    );
}

#[test]
fn bam_out_of_order_hints_sort() {
    for (so, named) in [("coordinate", false), ("queryname", true)] {
        let mut sam = sam_header(&[("chr1", 1000), ("chr2", 1000)], so);
        sam.push_str(&sam_read("a", 0, "chr2", 11, "10M", 10));
        sam.push_str(&sam_read("b", 0, "chr1", 11, "10M", 10));
        let bam = sam_to_bam(&sam);
        let o = ropts(AlignOptions::default(), InvalidPolicy::Fail, None);
        let mut r = open_reader_from("x.bam", Box::new(Cursor::new(bam)), &o).unwrap();
        let mut table = chrom_table_for(None, &[&*r], gloci_core::ChromOrder::Karyotype).unwrap();
        let Header::Alignment(h) = r.header().clone() else {
            panic!()
        };
        let mut batch = RecordBatch::default();
        r.read_batch(&mut table, &mut batch, 10).unwrap();
        let e = OrderGuard::new("x.bam", Format::Bam)
            .with_sort_order(h.sort_order)
            .check(&table, &batch)
            .unwrap_err()
            .to_string();
        assert!(e.contains("run: gloci sort x.bam -o x.sorted.bam"), "{e}");
        assert!(e.contains("--sort"), "{e}");
        assert_eq!(e.contains("SO:queryname"), named, "{e}");
        assert!(e.contains("record 2") && e.contains("chr1 11"), "{e}");
    }
}

#[test]
fn bam_malformed_records() {
    let h = header(&[("chr1", 100_000), ("chr2", 1000)], "coordinate");
    let good = || raw_record(0, 10, b"ok", 0, &[10 << 4 | M], 10, b"");
    let frame = |r: &[u8]| {
        let mut v = (r.len() as u32).to_le_bytes().to_vec();
        v.extend_from_slice(r);
        v
    };
    let file = |h: &AlnHeader, tail: Vec<u8>| {
        let mut v = Vec::new();
        h.write_bam(&mut v).unwrap();
        v.extend(frame(&good()));
        v.extend(tail);
        bgzip(&v)
    };
    let mut cases: Vec<(&str, Vec<u8>, &str, &str)> = Vec::new();
    // (shape, bytes after the good record, place, text of the message)
    cases.push((
        "block_size below 32",
        frame(&good()[..20]),
        "record 2 of chr1",
        "block_size 20",
    ));
    cases.push((
        "block_size past the data",
        {
            let mut v = 500u32.to_le_bytes().to_vec();
            v.extend_from_slice(&good()[..40]);
            v
        },
        "record 2 of chr1",
        "truncated",
    ));
    let mut empty = raw_record(0, 10, b"", 0, &[10 << 4 | M], 10, b"");
    cases.push((
        "empty name",
        frame(&{
            empty.shrink_to_fit();
            empty
        }),
        "record 2 of chr1",
        "read name is empty",
    ));
    let mut unterminated = raw_record(0, 10, b"abc", 0, &[10 << 4 | M], 10, b"");
    unterminated[32 + 3] = b'x';
    cases.push((
        "unterminated name",
        frame(&unterminated),
        "record 2 of chr1",
        "NUL",
    ));
    cases.push((
        "refID out of range",
        frame(&raw_record(5, 10, b"r", 0, &[10 << 4 | M], 10, b"")),
        "record 1 of refID 5",
        "outside the 2 references",
    ));
    cases.push((
        "pos below -1",
        frame(&raw_record(0, -2, b"r", 0, &[10 << 4 | M], 10, b"")),
        "record 2 of chr1",
        "pos -2",
    ));
    cases.push((
        "mapped read with refID -1",
        frame(&raw_record(-1, 5, b"r", 0, &[10 << 4 | M], 10, b"")),
        "record 1 of unplaced reads",
        "refID -1",
    ));
    cases.push((
        "mapped read with pos -1",
        frame(&raw_record(0, -1, b"r", 0, &[10 << 4 | M], 10, b"")),
        "record 2 of chr1",
        "pos -1",
    ));
    cases.push((
        "CIGAR op above 8",
        frame(&raw_record(0, 10, b"r", 0, &[10 << 4 | 9], 10, b"")),
        "record 2 of chr1",
        "operation code 9",
    ));
    cases.push((
        "l_seq differs from the CIGAR",
        frame(&raw_record(0, 10, b"r", 0, &[10 << 4 | M], 12, b"")),
        "record 2 of chr1",
        "covers 10 bases",
    ));
    for (shape, tail, place, text) in &cases {
        let bam = file(&h, tail.clone());
        let e = run_with(
            &bam,
            &ropts(AlignOptions::default(), InvalidPolicy::Fail, None),
        )
        .err()
        .unwrap_or_else(|| panic!("{shape}: no error"));
        let msg = e.to_string();
        assert!(msg.contains(place), "{shape}: {msg}");
        assert!(msg.contains(text), "{shape}: {msg}");
        assert!(msg.contains("t.bam"), "{shape}: {msg}");
    }
    // Under --skip-invalid a record whose frame is intact is skipped and counted, with
    // the same place; a truncated file is still an error.
    for (shape, tail, place, _) in &cases {
        if *shape == "block_size past the data" {
            continue;
        }
        let bam = file(&h, tail.clone());
        let r = run_with(
            &bam,
            &ropts(AlignOptions::default(), InvalidPolicy::Skip, None),
        )
        .unwrap_or_else(|e| panic!("{shape}: {e}"));
        assert_eq!(r.diag.skipped, 1, "{shape}");
        assert!(
            r.diag.examples[0].contains(place),
            "{shape}: {:?}",
            r.diag.examples
        );
        assert_eq!(r.recs.len(), 1, "{shape}");
    }
    let bam = file(&h, cases[1].1.clone());
    assert!(
        run_with(
            &bam,
            &ropts(AlignOptions::default(), InvalidPolicy::Skip, None)
        )
        .is_err()
    );
}

#[test]
fn bam_huge_header_linear() {
    let n = 200_000usize;
    let names: Vec<String> = (0..n).map(|i| format!("contig{i}")).collect();
    let refs: Vec<(&str, u64)> = names.iter().map(|s| (s.as_str(), 10_000)).collect();
    let h = header(&refs, "coordinate");
    let records: Vec<Vec<u8>> = (0..n)
        .step_by(1000)
        .map(|i| raw_record(i as i32, 100, b"r", 0, &[10 << 4 | M], 10, b""))
        .collect();
    let bam = bam_file(&h, &records);
    let t = std::time::Instant::now();
    let r = run(&bam, AlignOptions::default());
    let took = t.elapsed();
    assert_eq!(r.recs.len(), 200);
    assert_eq!(r.table.len(), n);
    assert_eq!(r.recs[199].chrom, "contig199000");
    assert!(took.as_secs() < 10, "200,000 @SQ took {took:?}");
}

#[test]
fn bam_chrom_map_renames_header() {
    let bam = read_fixture("order.bam");
    let map = Arc::new(ChromMap::parse("m", b"chr1\t1\n").unwrap());
    let r = run_with(
        &bam,
        &ropts(AlignOptions::default(), InvalidPolicy::Fail, Some(map)),
    )
    .unwrap();
    assert_eq!(r.hdr.refs[1].0, b"1");
    assert_eq!(r.hdr.refs[0].0, b"chrM");
    let text = String::from_utf8_lossy(&r.hdr.text).into_owned();
    assert!(text.contains("@SQ\tSN:1\tLN:248956422\n"), "{text}");
    assert!(!text.contains("SN:chr1\t"), "{text}");
    assert!(text.contains("@SQ\tSN:chr2\t"), "{text}");
    let chroms: Vec<_> = r.recs.iter().map(|g| g.chrom.as_str()).collect();
    assert_eq!(chroms.iter().filter(|c| **c == "1").count(), 3);
    assert!(!chroms.contains(&"chr1"));
    // A map that would merge two contigs is a usage error naming both.
    let clash = Arc::new(ChromMap::parse("m", b"chr1\tchr2\n").unwrap());
    let e = run_with(
        &bam,
        &ropts(AlignOptions::default(), InvalidPolicy::Fail, Some(clash)),
    )
    .err()
    .unwrap();
    assert!(
        e.to_string().contains("chr1") && e.to_string().contains("chr2"),
        "{e}"
    );
}

#[test]
fn bam_slim_records() {
    let bam = read_fixture("cigars.bam");
    let slim = run(
        &bam,
        AlignOptions {
            score_tag: Some(*b"NM"),
            ..AlignOptions::default()
        },
    );
    let full = run(
        &bam,
        AlignOptions {
            full_records: true,
            ..AlignOptions::default()
        },
    );
    assert_eq!(slim.recs.len(), 9);
    assert_eq!(full.recs.len(), 9);
    for (s, f) in slim.recs.iter().zip(&full.recs) {
        assert!(s.raw.len() < f.raw.len(), "{}", name_of(f, &full.hdr));
        let (a, b) = (
            aln_ref(Format::Bam, &s.raw, &slim.hdr).unwrap(),
            aln_ref(Format::Bam, &f.raw, &full.hdr).unwrap(),
        );
        assert_eq!(a.name(), b.name());
        assert_eq!((a.flag(), a.mapq()), (b.flag(), b.mapq()));
        assert_eq!(ops_of(s, &slim.hdr), ops_of(f, &full.hdr));
        assert_eq!((s.span, s.strand), (f.span, f.strand));
        assert_eq!(a.span(100_000), b.span(100_000));
        assert_eq!(a.int_tag(*b"NM"), b.int_tag(*b"NM"));
        assert!(a.int_tag(*b"NM").is_some());
    }
    // Without a score tag, the tag is dropped.
    let bare = run(&bam, AlignOptions::default());
    let a = aln_ref(Format::Bam, &bare.recs[0].raw, &bare.hdr).unwrap();
    assert_eq!(a.int_tag(*b"NM"), None);
}

#[test]
fn blocks_follow_d105() {
    let bam = read_fixture("cigars.bam");
    let r = run(&bam, AlignOptions::default());
    let blocks_of = |name: &str, split_d: bool| -> Vec<(u64, u64)> {
        let g = r.recs.iter().find(|g| name_of(g, &r.hdr) == name).unwrap();
        let a = aln_ref(Format::Bam, &g.raw, &r.hdr).unwrap();
        let mut out = Vec::new();
        blocks(&a, split_d, &mut out);
        out.iter().map(|i| (i.start, i.end)).collect()
    };
    assert_eq!(blocks_of("c_10M", false), [(100, 110)]);
    assert_eq!(blocks_of("c_skip", false), [(600, 605), (705, 710)]);
    assert_eq!(blocks_of("c_del", false), [(700, 712)]);
    assert_eq!(blocks_of("c_del", true), [(700, 705), (707, 712)]);
    assert_eq!(blocks_of("c_ins", false), [(400, 410)]);
    // No aligned bases: the one-base span is the block.
    assert_eq!(blocks_of("c_allclip", false), [(800, 801)]);
}

#[test]
fn fixtures_match_sources() {
    use noodles_sam::alignment::io::Write as _;
    // samtools and noodles pick different integer widths for a tag, so the check is
    // on the decoded records: each BAM, read back as SAM, is its source.
    for set in ["cigars", "flags", "order"] {
        let sam = std::fs::read_to_string(fixture(&format!("{set}.sam"))).unwrap();
        let theirs = read_fixture(&format!("{set}.bam"));
        let mut rd = noodles_bam::io::Reader::new(&theirs[..]);
        let header = rd.read_header().unwrap();
        let mut w = noodles_sam::io::Writer::new(Vec::new());
        for rec in rd.records() {
            w.write_alignment_record(&header, &rec.unwrap()).unwrap();
        }
        let decoded = String::from_utf8(w.into_inner()).unwrap();
        let source: String = sam
            .lines()
            .filter(|l| !l.starts_with('@'))
            .map(|l| format!("{l}\n"))
            .collect();
        assert_eq!(decoded, source, "{set}.bam does not decode to {set}.sam");
        // The references match the @SQ lines.
        let want: Vec<(Vec<u8>, u64)> = sam
            .lines()
            .filter_map(|l| l.strip_prefix("@SQ\tSN:"))
            .map(|l| {
                let (n, len) = l.split_once("\tLN:").unwrap();
                (n.as_bytes().to_vec(), len.parse().unwrap())
            })
            .collect();
        let mut v = Vec::new();
        noodles_bgzf::io::Reader::new(&theirs[..])
            .read_to_end(&mut v)
            .unwrap();
        assert_eq!(
            AlnHeader::read_bam(&mut &v[..]).unwrap().refs,
            want,
            "{set}"
        );
        // And the encoder the tests use agrees with the fixture on the references.
        let ours = sam_to_bam(&sam);
        let mut v = Vec::new();
        noodles_bgzf::io::Reader::new(&ours[..])
            .read_to_end(&mut v)
            .unwrap();
        assert_eq!(
            AlnHeader::read_bam(&mut &v[..]).unwrap().refs,
            want,
            "{set}"
        );
    }
}

#[test]
fn bam_reads_identically_at_every_thread_count() {
    let bam = read_fixture("order.bam");
    let mut seen = Vec::new();
    for t in [1, 2, 8] {
        let mut o = ropts(AlignOptions::default(), InvalidPolicy::Fail, None);
        o.threads = NonZeroUsize::new(t).unwrap();
        let r = run_with(&bam, &o).unwrap();
        seen.push(
            r.recs
                .iter()
                .map(|g| (g.chrom.clone(), g.span, g.loc, g.raw.clone()))
                .collect::<Vec<_>>(),
        );
    }
    assert_eq!(seen[0], seen[1]);
    assert_eq!(seen[0], seen[2]);
}

#[test]
fn fuzz_bam() {
    // 2,000 seeded flips and truncations of the decoded BAM, re-wrapped in BGZF.
    let plain = {
        let mut v = Vec::new();
        noodles_bgzf::io::Reader::new(&read_fixture("cigars.bam")[..])
            .read_to_end(&mut v)
            .unwrap();
        v
    };
    let mut x = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = move || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    };
    for i in 0..2000 {
        let mut d = plain.clone();
        for _ in 0..1 + next() % 4 {
            let at = (next() % d.len() as u64) as usize;
            d[at] = next() as u8;
        }
        if i % 3 == 0 {
            d.truncate((next() % d.len() as u64) as usize);
        }
        let mut input = vec![(next() & 0xff) as u8];
        input.extend(bgzip(&d));
        gloci_io::fuzz_bodies::read_bam(&input);
    }
}

#[test]
fn bam_read_past_the_reference_end_is_skipped_and_counted() {
    // A mapped read at or past LN covers no reference base: skipped and counted, with one
    // warning per input; a read that only overhangs is still cut (D101).
    let h = header(&[("chr1", 100)], "coordinate");
    let rec = |pos: i32, name: &[u8]| raw_record(0, pos, name, 0, &[10 << 4 | M], 10, b"");
    let bam = bam_file(
        &h,
        &[
            rec(50, b"in"),
            rec(95, b"over"),
            rec(100, b"at"),
            rec(150, b"past"),
        ],
    );
    let r = run(&bam, AlignOptions::default());
    let names: Vec<_> = r.recs.iter().map(|g| name_of(g, &r.hdr)).collect();
    assert_eq!(names, ["in", "over"]);
    assert_eq!(
        r.recs[1].span,
        Interval {
            start: 95,
            end: 100
        }
    );
    assert_eq!(r.diag.off_reference, 2);
    assert!(r.recs.iter().all(|g| g.span.start < g.span.end));
}

#[test]
fn bam_error_place_names_the_file_record() {
    let h = header(&[("chr1", 100_000), ("chr2", 1000)], "coordinate");
    let good = |ref_id: i32| raw_record(ref_id, 10, b"ok", 0, &[10 << 4 | M], 10, b"");
    let bad = raw_record(1, -2, b"r", 0, &[10 << 4 | M], 10, b"");
    let bam = bam_file(&h, &[good(0), good(0), good(1), bad]);
    let e = run_with(
        &bam,
        &ropts(AlignOptions::default(), InvalidPolicy::Fail, None),
    )
    .err()
    .unwrap()
    .to_string();
    assert!(e.contains("record 2 of chr2 (record 4 of the file)"), "{e}");
    // Only refID -1 is "unplaced reads"; another negative id is printed as it is.
    let neg = bam_file(&h, &[raw_record(-5, 10, b"r", 4, &[], 0, b"")]);
    let e = run_with(
        &neg,
        &ropts(AlignOptions::default(), InvalidPolicy::Fail, None),
    )
    .err()
    .unwrap()
    .to_string();
    assert!(
        e.contains("record 1 of refID -5 (record 1 of the file)"),
        "{e}"
    );
}
