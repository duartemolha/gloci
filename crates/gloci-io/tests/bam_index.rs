//! BAM indexes (D107): `.bai` and `.csi` built by gloci against samtools' `.bai`, queries through
//! them against a scan of the file, metadata and long references.

use gloci_core::{Interval, RegionSet, VcfSpan};
use gloci_io::*;
use gloci_testkit::fixtures::sam_to_bam;
use noodles_csi::BinningIndex as _;
use noodles_csi::binning_index::ReferenceSequence as _;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/aln")
        .join(name)
}

fn threads() -> NonZeroUsize {
    NonZeroUsize::new(2).unwrap()
}

fn ropts() -> ReaderOptions {
    ReaderOptions {
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
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    }
}

type Regions = Vec<(Vec<u8>, Interval)>;

/// Every record of `bam` in file order, as raw bytes.
fn all_records(bam: &Path) -> (std::sync::Arc<AlnHeader>, Vec<Vec<u8>>) {
    let mut src = open_aln_source(bam.to_str().unwrap(), &ropts()).unwrap();
    let hdr = src.header().clone();
    let mut out = Vec::new();
    let mut buf = Vec::new();
    while src.next_record(&mut buf).unwrap().is_some() {
        out.push(std::mem::take(&mut buf));
    }
    (hdr, out)
}

/// The records a region query must return (§4.3 on the D116 span), by scanning.
fn scan(hdr: &AlnHeader, recs: &[Vec<u8>], regions: &Regions, inclusive: bool) -> Vec<Vec<u8>> {
    let mut by_name: std::collections::HashMap<&[u8], Vec<Interval>> = Default::default();
    for (n, iv) in regions {
        by_name.entry(n).or_default().push(*iv);
    }
    let sets: std::collections::HashMap<&[u8], RegionSet> = by_name
        .into_iter()
        .map(|(n, v)| (n, RegionSet::new(v)))
        .collect();
    recs.iter()
        .filter(|r| {
            let a = aln_ref(Format::Bam, r, hdr).unwrap();
            alignment_query_span(&a).is_some_and(|(rid, iv)| {
                sets.get(&hdr.refs[rid].0[..])
                    .is_some_and(|s| s.matches(iv, inclusive))
            })
        })
        .cloned()
        .collect()
}

fn through_index(bam: &Path, regions: &Regions, inclusive: bool) -> Vec<Vec<u8>> {
    let r = open_reader(bam.to_str().unwrap(), &ropts()).unwrap();
    let ix = IndexedFile::open(bam, Format::Bam, r.header().clone(), &ropts()).unwrap();
    let mut src = ix.aln_regions_source(regions, inclusive).unwrap();
    let mut out = Vec::new();
    let mut buf = Vec::new();
    while src.next_record(&mut buf).unwrap().is_some() {
        out.push(std::mem::take(&mut buf));
    }
    out
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n.max(1)
    }
}

/// 300 seeded regions around the reads' ends, bin edges and at random.
fn random_regions(hdr: &AlnHeader, recs: &[Vec<u8>], seed: u64) -> Vec<Regions> {
    let mut points: Vec<Vec<u64>> =
        vec![vec![0, 1, 16_383, 16_384, 16_385, 131_071, 131_072]; hdr.refs.len()];
    for r in recs {
        let a = aln_ref(Format::Bam, r, hdr).unwrap();
        if let Some((rid, iv)) = alignment_query_span(&a) {
            for p in [iv.start, iv.end] {
                points[rid].extend([p.saturating_sub(1), p, p + 1]);
            }
        }
    }
    let mut rng = Rng(seed);
    (0..300)
        .map(|_| {
            let mut regions = Regions::new();
            for _ in 0..1 + rng.below(3) {
                let rid = rng.below(hdr.refs.len() as u64) as usize;
                let pts = &points[rid];
                let pick = |rng: &mut Rng| {
                    if rng.below(4) == 0 {
                        rng.below(hdr.refs[rid].1.max(1) + 50)
                    } else {
                        pts[rng.below(pts.len() as u64) as usize]
                    }
                };
                let (a, b) = (pick(&mut rng), pick(&mut rng));
                let (start, end) = (a.min(b), a.max(b));
                regions.push((hdr.refs[rid].0.clone(), Interval { start, end }));
            }
            regions
        })
        .collect()
}

fn setup(dir: &Path, bam: &str, samtools_bai: bool) -> PathBuf {
    let p = dir.join("x.bam");
    std::fs::copy(fixture(bam), &p).unwrap();
    if samtools_bai {
        std::fs::copy(fixture(&format!("{bam}.bai")), dir.join("x.bam.bai")).unwrap();
    } else {
        write_alignment_index(&p, &p, false, threads()).unwrap();
    }
    p
}

#[test]
fn bai_queries_match_samtools_bai() {
    for bam in [
        "cigars.bam",
        "flags.bam",
        "order.bam",
        "pastend.bam",
        "tags.bam",
        "negtag.bam",
        "splice-end.bam",
    ] {
        let (hdr, recs) = all_records(&fixture(bam));
        let sets = random_regions(&hdr, &recs, 0x9e37_79b9 ^ bam.len() as u64);
        let (d1, d2) = (tempfile::tempdir().unwrap(), tempfile::tempdir().unwrap());
        let samtools = setup(d1.path(), bam, true);
        let gloci = setup(d2.path(), bam, false);
        assert_eq!(find_index(&samtools).unwrap().kind, IndexKind::Bai);
        assert_eq!(find_index(&gloci).unwrap().kind, IndexKind::Bai);
        let mut nonempty = 0;
        for (i, regions) in sets.iter().enumerate() {
            for inclusive in [false, true] {
                let want = scan(&hdr, &recs, regions, inclusive);
                nonempty += usize::from(!want.is_empty());
                for (who, p) in [("samtools", &samtools), ("gloci", &gloci)] {
                    assert_eq!(
                        through_index(p, regions, inclusive),
                        want,
                        "{bam} region set {i} {regions:?} through the {who} .bai, inclusive {inclusive}"
                    );
                }
            }
        }
        assert!(
            nonempty > 100,
            "{bam}: only {nonempty} of 600 queries found reads"
        );
    }
}

/// Reads around the bin edges 16,384 and 131,072, a read that spans several windows, zero
/// reference-length reads (an all-clip read and placed unmapped reads), through gloci's own
/// `.bai` and `.csi`.
#[test]
fn bai_and_csi_queries_match_a_scan_on_generated_reads() {
    let mut sam = String::from(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000000\n@SQ\tSN:chr2\tLN:500\n",
    );
    let mut reads: Vec<(u64, String)> = Vec::new();
    let mut rng = Rng(77);
    for edge in [16_384u64, 131_072, 262_144] {
        for d in -3i64..=3 {
            for len in [1u64, 10, 100] {
                let pos = (edge as i64 + d - 5) as u64;
                reads.push((pos, format!("{len}M")));
            }
        }
    }
    reads.push((1000, "5000M".into()));
    reads.push((20_000, "3S5I".into()));
    reads.push((140_000, "70000M".into()));
    for _ in 0..200 {
        reads.push((1 + rng.below(999_000), format!("{}M", 1 + rng.below(300))));
    }
    reads.sort_by_key(|r| r.0);
    for (i, (pos, cigar)) in reads.iter().enumerate() {
        let (seq, qual) = match cigar.as_str() {
            "3S5I" => ("ACGTACGT".to_string(), "IIIIIIII".to_string()),
            c => {
                let n: usize = c.trim_end_matches('M').parse().unwrap();
                ("A".repeat(n), "I".repeat(n))
            }
        };
        sam.push_str(&format!(
            "r{i}\t0\tchr1\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}\n"
        ));
    }
    sam.push_str("u1\t4\tchr1\t16385\t0\t*\t*\t0\t0\tACGT\tIIII\n");
    sam.push_str("u2\t4\tchr1\t131073\t0\t*\t*\t0\t0\tACGT\tIIII\n");
    // placed unmapped reads must follow reads at lower positions; sorting by pos put them wrong
    let mut lines: Vec<&str> = sam.lines().collect();
    let (header, body) = lines.split_at_mut(3);
    body.sort_by_key(|l| l.split('\t').nth(3).unwrap().parse::<u64>().unwrap());
    let mut sorted = header.join("\n");
    sorted.push('\n');
    sorted.push_str(&body.join("\n"));
    sorted.push_str("\nc1\t0\tchr2\t10\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\nz\t4\t*\t0\t0\t*\t*\t0\t0\tAC\tII\n");
    let bam = sam_to_bam(&sorted);
    for csi in [false, true] {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("g.bam");
        std::fs::write(&p, &bam).unwrap();
        let ip = write_alignment_index(&p, &p, csi, threads()).unwrap();
        assert_eq!(ip.extension().unwrap(), if csi { "csi" } else { "bai" });
        let (hdr, recs) = all_records(&p);
        for (i, regions) in random_regions(&hdr, &recs, 5 + csi as u64)
            .iter()
            .enumerate()
        {
            for inclusive in [false, true] {
                assert_eq!(
                    through_index(&p, regions, inclusive),
                    scan(&hdr, &recs, regions, inclusive),
                    "csi {csi}, region set {i} {regions:?}, inclusive {inclusive}"
                );
            }
        }
    }
}

fn bai(p: &Path) -> noodles_bam::bai::Index {
    noodles_bam::bai::fs::read(p).unwrap()
}

fn meta<R: noodles_csi::binning_index::ReferenceSequence>(
    r: &R,
) -> Option<(
    noodles_bgzf::VirtualPosition,
    noodles_bgzf::VirtualPosition,
    u64,
    u64,
)> {
    r.metadata().map(|m| {
        (
            m.start_position(),
            m.end_position(),
            m.mapped_record_count(),
            m.unmapped_record_count(),
        )
    })
}

#[test]
fn bai_metadata() {
    for bam in [
        "cigars.bam",
        "flags.bam",
        "order.bam",
        "pastend.bam",
        "tags.bam",
    ] {
        let dir = tempfile::tempdir().unwrap();
        let p = setup(dir.path(), bam, false);
        let (ours, theirs) = (
            bai(&dir.path().join("x.bam.bai")),
            bai(&fixture(&format!("{bam}.bai"))),
        );
        assert_eq!(
            ours.unplaced_unmapped_record_count(),
            theirs.unplaced_unmapped_record_count(),
            "{bam}: n_no_coor"
        );
        assert_eq!(
            ours.reference_sequences().len(),
            theirs.reference_sequences().len()
        );
        for (i, (a, b)) in ours
            .reference_sequences()
            .iter()
            .zip(theirs.reference_sequences())
            .enumerate()
        {
            assert_eq!(meta(a), meta(b), "{bam}: metadata of reference {i}");
        }
        // The fixtures are small enough that even the bin order is samtools' (on larger files
        // htslib's hash order of bins within a reference can differ from gloci's sorted one).
        assert_eq!(
            std::fs::read(dir.path().join("x.bam.bai")).unwrap(),
            std::fs::read(fixture(&format!("{bam}.bai"))).unwrap(),
            "{bam}: .bai bytes"
        );
        drop(p);
    }
    // flags.bam: chr1 has one placed unmapped read, and one read has no reference.
    let t = bai(&fixture("flags.bam.bai"));
    assert_eq!(t.unplaced_unmapped_record_count(), Some(1));
    let m = t.reference_sequences()[0].metadata().unwrap();
    assert_eq!(
        (m.mapped_record_count(), m.unmapped_record_count()),
        (10, 1)
    );
}

fn long_reference_bam() -> Vec<u8> {
    let mut sam = String::from(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:big\tLN:600000000\n@SQ\tSN:small\tLN:1000\n",
    );
    for p in [1u64, 100, 20_000, 300_000_000, 540_000_000, 599_999_990] {
        sam.push_str(&format!(
            "r{p}\t0\tbig\t{p}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
        ));
    }
    sam.push_str(
        "u1\t4\tsmall\t5\t0\t*\t*\t0\t0\tACGT\tIIII\nu2\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n",
    );
    sam_to_bam(&sam)
}

#[test]
fn csi_for_long_references() {
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("l.bam");
    std::fs::write(&p, long_reference_bam()).unwrap();
    // A stale .bai next to it goes away when the .csi is written.
    std::fs::write(dir.path().join("l.bam.bai"), b"stale").unwrap();
    let ip = write_alignment_index(&p, &p, false, threads()).unwrap();
    assert_eq!(ip, dir.path().join("l.bam.csi"));
    assert!(!dir.path().join("l.bam.bai").exists());
    let info = find_index(&p).unwrap();
    assert_eq!(info.kind, IndexKind::Csi);
    let csi = noodles_csi::fs::read(&ip).unwrap();
    assert_eq!(csi.min_shift(), 14);
    assert_eq!(csi.depth(), 6);
    assert_eq!(csi.unplaced_unmapped_record_count(), Some(1));
    let (hdr, recs) = all_records(&p);
    let name = |i: usize| hdr.refs[i].0.clone();
    for regions in [
        vec![(
            name(0),
            Interval {
                start: 299_999_990,
                end: 300_000_010,
            },
        )],
        vec![(
            name(0),
            Interval {
                start: 539_999_999,
                end: 540_000_005,
            },
        )],
        vec![(
            name(0),
            Interval {
                start: 599_999_980,
                end: 600_000_000,
            },
        )],
        vec![
            (
                name(0),
                Interval {
                    start: 0,
                    end: 1 << 40,
                },
            ),
            (name(1), Interval { start: 0, end: 10 }),
        ],
        vec![(name(1), Interval { start: 4, end: 5 })],
    ] {
        let want = scan(&hdr, &recs, &regions, false);
        assert_eq!(through_index(&p, &regions, false), want, "{regions:?}");
    }

    // `--csi` on a small file; a .bai is replaced.
    let dir = tempfile::tempdir().unwrap();
    let small = setup(dir.path(), "flags.bam", true);
    assert_eq!(find_index(&small).unwrap().kind, IndexKind::Bai);
    let ip = write_alignment_index(&small, &small, true, threads()).unwrap();
    assert_eq!(ip, dir.path().join("x.bam.csi"));
    assert!(!dir.path().join("x.bam.bai").exists());
    let (hdr, recs) = all_records(&small);
    let regions: Regions = vec![(
        hdr.refs[0].0.clone(),
        Interval {
            start: 100,
            end: 710,
        },
    )];
    assert_eq!(
        through_index(&small, &regions, false),
        scan(&hdr, &recs, &regions, false)
    );
    // and writes a .bai again without it
    let ip = write_alignment_index(&small, &small, false, threads()).unwrap();
    assert_eq!(ip, dir.path().join("x.bam.bai"));
    assert!(!dir.path().join("x.bam.csi").exists());
}

#[test]
fn find_index_names_for_alignments() {
    let dir = tempfile::tempdir().unwrap();
    let bam = dir.path().join("a.bam");
    let cram = dir.path().join("a.cram");
    std::fs::write(&bam, b"x").unwrap();
    std::fs::write(&cram, b"x").unwrap();
    assert!(find_index(&bam).is_none());
    for (data, idx, kind) in [
        (&bam, "a.bai", IndexKind::Bai),
        (&bam, "a.bam.bai", IndexKind::Bai),
        (&bam, "a.bam.csi", IndexKind::Csi),
        (&cram, "a.crai", IndexKind::Crai),
        (&cram, "a.cram.crai", IndexKind::Crai),
        (&cram, "a.cram.csi", IndexKind::Csi),
    ] {
        let p = dir.path().join(idx);
        std::fs::write(&p, b"i").unwrap();
        let mtime = std::fs::metadata(&p).unwrap().modified().unwrap();
        // Make it newer than the data, and newer than any earlier name.
        let f = std::fs::File::options().write(true).open(&p).unwrap();
        f.set_modified(mtime + std::time::Duration::from_secs(10))
            .unwrap();
        let info = find_index(data).unwrap();
        assert_eq!((info.kind, &info.path), (kind, &p), "{idx}");
        assert!(!info.stale);
        std::fs::remove_file(&p).unwrap();
    }
    // A BED file never takes a .bai.
    let bed = dir.path().join("a.bed.gz");
    std::fs::write(&bed, b"x").unwrap();
    std::fs::write(dir.path().join("a.bed.gz.bai"), b"i").unwrap();
    assert!(find_index(&bed).is_none());
}

/// Open `bam` with the index bytes `bad` next to it and query it: an error is fine, a panic
/// or an abort is not.
fn try_index(bam: &Path, name: &str, bad: &[u8]) {
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("x.bam");
    std::fs::copy(bam, &p).unwrap();
    std::fs::write(dir.path().join(name), bad).unwrap();
    let r = open_reader(p.to_str().unwrap(), &ropts()).unwrap();
    let Ok(ix) = IndexedFile::open(&p, Format::Bam, r.header().clone(), &ropts()) else {
        return;
    };
    let hdr_names: Vec<Vec<u8>> = ix.reference_names();
    for name in hdr_names.iter().take(3) {
        let regions: Regions = vec![(
            name.clone(),
            Interval {
                start: 0,
                end: 1 << 40,
            },
        )];
        let _ = ix.weight(name);
        let _ = ix.aln_regions_source(&regions, false).and_then(|mut s| {
            let mut buf = Vec::new();
            while s.next_record(&mut buf)?.is_some() {
                buf.clear();
            }
            Ok(())
        });
        let _ = ix
            .chromosomes_reader(std::slice::from_ref(name))
            .and_then(|mut r| {
                let mut batch = RecordBatch::default();
                let mut chroms = gloci_core::ChromTable::first_appearance();
                while r.read_batch(&mut chroms, &mut batch, 8)? {
                    batch.clear();
                }
                Ok(())
            });
    }
}

#[test]
fn corrupt_index_is_an_error_not_a_panic() {
    let bai = std::fs::read(fixture("flags.bam.bai")).unwrap();
    let mut cases: Vec<Vec<u8>> = vec![
        Vec::new(),
        b"BAI\x01".to_vec(),
        b"nope nope nope nope".to_vec(),
        bai[..bai.len() / 2].to_vec(),
    ];
    // a huge reference count, and a huge bin count in the first reference
    for at in [4, 8] {
        let mut huge = bai.clone();
        huge[at..at + 4].copy_from_slice(&i32::MAX.to_le_bytes());
        cases.push(huge);
    }
    let mut rng = Rng(0xbad_1dea);
    for _ in 0..150 {
        let mut m = bai.clone();
        for _ in 0..1 + rng.below(4) {
            let at = rng.below(m.len() as u64) as usize;
            m[at] = rng.next() as u8;
        }
        cases.push(m);
    }
    for bad in &cases {
        try_index(&fixture("flags.bam"), "x.bam.bai", bad);
    }
    // the same for a CSI
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("c.bam");
    std::fs::copy(fixture("flags.bam"), &p).unwrap();
    write_alignment_index(&p, &p, true, threads()).unwrap();
    let csi = std::fs::read(dir.path().join("c.bam.csi")).unwrap();
    let mut plain = Vec::new();
    {
        use std::io::Read as _;
        noodles_bgzf::io::Reader::new(&csi[..])
            .read_to_end(&mut plain)
            .unwrap();
    }
    let recompress = |b: &[u8]| {
        use std::io::Write as _;
        let mut w = noodles_bgzf::io::Writer::new(Vec::new());
        w.write_all(b).unwrap();
        w.finish().unwrap()
    };
    let mut csis = vec![
        Vec::new(),
        csi[..csi.len() / 2].to_vec(),
        recompress(b"CSI\x01"),
    ];
    for _ in 0..150 {
        let mut m = plain.clone();
        for _ in 0..1 + rng.below(4) {
            let at = rng.below(m.len() as u64) as usize;
            m[at] = rng.next() as u8;
        }
        csis.push(recompress(&m));
    }
    for bad in &csis {
        try_index(&fixture("flags.bam"), "x.bam.csi", bad);
    }
}

/// What the whole-chromosome reader (per-chromosome parallelism) yields equals what the
/// sequential reader yields for that chromosome: the same reads, spans and strands, unmapped
/// reads skipped, with slim or whole records.
#[test]
fn chromosome_reader_equals_the_sequential_reader() {
    for bam in [
        "cigars.bam",
        "flags.bam",
        "order.bam",
        "pastend.bam",
        "tags.bam",
    ] {
        for full in [false, true] {
            let mut o = ropts();
            o.accept.alignments.as_mut().unwrap().full_records = full;
            let path = fixture(bam);
            let mut seq = open_reader(path.to_str().unwrap(), &o).unwrap();
            let header = seq.header().clone();
            let mut chroms = gloci_core::ChromTable::first_appearance();
            let mut batch = RecordBatch::default();
            let mut want: std::collections::BTreeMap<String, Vec<(Interval, Vec<u8>)>> =
                Default::default();
            loop {
                batch.clear();
                let more = seq.read_batch(&mut chroms, &mut batch, 4).unwrap();
                for i in 0..batch.len() {
                    let m = batch.recs[i];
                    want.entry(String::from_utf8_lossy(chroms.name(m.chrom)).into_owned())
                        .or_default()
                        .push((m.span, batch.raw(i).to_vec()));
                }
                if !more {
                    break;
                }
            }
            let dir = tempfile::tempdir().unwrap();
            let p = setup(dir.path(), bam, true);
            let ix = IndexedFile::open(&p, Format::Bam, header.clone(), &o).unwrap();
            let Header::Alignment(aln) = &header else {
                panic!("alignment header")
            };
            // every chromosome with records is listed, none without
            let with_records: Vec<String> = want.keys().cloned().collect();
            let mut listed: Vec<String> = ix
                .non_empty_reference_names()
                .iter()
                .map(|n| String::from_utf8_lossy(n).into_owned())
                .collect();
            listed.retain(|n| with_records.contains(n) || ix.weight(n.as_bytes()) > 0);
            for n in &with_records {
                assert!(listed.contains(n), "{bam}: {n} has reads but is not listed");
            }
            for (name, _) in &aln.refs {
                let n = String::from_utf8_lossy(name).into_owned();
                let mut r = ix.chromosomes_reader(std::slice::from_ref(name)).unwrap();
                let mut tab = gloci_core::ChromTable::first_appearance();
                let mut got = Vec::new();
                let mut b = RecordBatch::default();
                loop {
                    b.clear();
                    let more = r.read_batch(&mut tab, &mut b, 3).unwrap();
                    for i in 0..b.len() {
                        got.push((b.recs[i].span, b.raw(i).to_vec()));
                    }
                    if !more {
                        break;
                    }
                }
                assert_eq!(
                    got,
                    want.get(&n).cloned().unwrap_or_default(),
                    "{bam} full {full}: reads of {n}"
                );
            }
        }
    }
}

/// A read with a reference and POS 0 (BAM pos -1) is an unmapped read of that reference
/// (htslib counts it in the reference's unmapped total, not in `n_no_coor`), and the index is
/// samtools' byte for byte. Where htslib refuses to index, gloci refuses too.
#[test]
fn reads_with_a_reference_and_no_position() {
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("x.bam");
    std::fs::copy(fixture("pos-minus-one.bam"), &p).unwrap();
    write_alignment_index(&p, &p, false, threads()).unwrap();
    assert_eq!(
        std::fs::read(dir.path().join("x.bam.bai")).unwrap(),
        std::fs::read(fixture("pos-minus-one.bam.bai")).unwrap()
    );
    let t = bai(&dir.path().join("x.bam.bai"));
    assert_eq!(t.unplaced_unmapped_record_count(), Some(1));
    let m = meta(&t.reference_sequences()[0]).unwrap();
    assert_eq!((m.2, m.3), (1, 2));

    let head = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:1000\n@SQ\tSN:c2\tLN:1000\n";
    let r = |n: &str, f: u16, c: &str, pos: u64| {
        let cig = if f & 4 != 0 { "*" } else { "4M" };
        format!("{n}\t{f}\t{c}\t{pos}\t0\t{cig}\t*\t0\t0\tACGT\tIIII\n")
    };
    let u4 = "u4\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n".to_string();
    // the three files htslib refuses: pos -1 after a read at a position (A, D), and a
    // reference read after unplaced reads (B)
    let a = [
        r("r1", 0, "c1", 10),
        r("r3", 0, "c2", 10),
        r("u3", 4, "c2", 0),
        u4.clone(),
    ]
    .concat();
    let b = [
        r("r1", 0, "c1", 10),
        r("r3", 0, "c2", 10),
        u4.clone(),
        r("u3", 4, "c2", 0),
    ]
    .concat();
    let d = [
        r("r1", 0, "c1", 10),
        r("r3", 0, "c2", 10),
        r("u3", 4, "c2", 0),
        r("u1", 4, "c2", 0),
        u4.clone(),
    ]
    .concat();
    for (name, body) in [("A", a), ("B", b), ("D", d)] {
        let p = dir.path().join(format!("{name}.bam"));
        std::fs::write(&p, sam_to_bam(&format!("{head}{body}"))).unwrap();
        let e = write_alignment_index(&p, &p, false, threads()).unwrap_err();
        let msg = e.to_string();
        assert!(matches!(e, IoError::Order(_)), "{name}: {msg}");
        assert!(
            msg.contains("record") && msg.contains("gloci sort"),
            "{name}: {msg}"
        );
        assert!(!dir.path().join(format!("{name}.bam.bai")).exists());
    }
    // pos -1 first in its reference, then reads at positions, in two references: fine
    let ok = [
        r("u0", 4, "c1", 0),
        r("u1", 4, "c1", 0),
        r("r1", 0, "c1", 10),
        r("u2", 4, "c2", 0),
        u4,
    ]
    .concat();
    let p = dir.path().join("ok.bam");
    std::fs::write(&p, sam_to_bam(&format!("{head}{ok}"))).unwrap();
    write_alignment_index(&p, &p, false, threads()).unwrap();
    let t = bai(&dir.path().join("ok.bam.bai"));
    let m0 = meta(&t.reference_sequences()[0]).unwrap();
    assert_eq!((m0.2, m0.3), (1, 2));
    assert_eq!(t.unplaced_unmapped_record_count(), Some(1));
}

fn recompress(b: &[u8]) -> Vec<u8> {
    use std::io::Write as _;
    let mut w = noodles_bgzf::io::Writer::new(Vec::new());
    w.write_all(b).unwrap();
    w.finish().unwrap()
}

fn decompress(b: &[u8]) -> Vec<u8> {
    use std::io::Read as _;
    let mut out = Vec::new();
    noodles_bgzf::io::Reader::new(b)
        .read_to_end(&mut out)
        .unwrap();
    out
}

/// A corrupt `.tbi` (a huge count, a bin or offset out of range, random bytes) is a prompt
/// error, never a long run, a panic or an abort.
#[test]
fn corrupt_tbi_is_an_error_quickly() {
    let dir = tempfile::tempdir().unwrap();
    let bed = dir.path().join("a.bed.gz");
    let text: String = (0..300)
        .map(|i| {
            format!(
                "c{}\t{}\t{}\n",
                i / 100 + 1,
                (i % 100) * 50,
                (i % 100) * 50 + 40
            )
        })
        .collect();
    std::fs::write(&bed, gloci_testkit::fixtures::bgzip(text.as_bytes())).unwrap();
    write_index(&bed, &bed, Format::Bed, threads()).unwrap();
    let tbi = dir.path().join("a.bed.gz.tbi");
    let plain = decompress(&std::fs::read(&tbi).unwrap());
    let o = ropts();
    let started = std::time::Instant::now();
    let mut rng = Rng(0x7b1);
    let mut cases: Vec<Vec<u8>> = vec![
        Vec::new(),
        b"TBI\x01".to_vec(),
        plain[..plain.len() / 2].to_vec(),
    ];
    // every 32-bit field set huge, one at a time, over the first 200 bytes
    for at in (4..plain.len().min(200)).step_by(4) {
        let mut m = plain.clone();
        m[at..at + 4].copy_from_slice(&i32::MAX.to_le_bytes());
        cases.push(m);
    }
    for _ in 0..150 {
        let mut m = plain.clone();
        for _ in 0..1 + rng.below(4) {
            let at = rng.below(m.len() as u64) as usize;
            m[at] = rng.next() as u8;
        }
        cases.push(m);
    }
    for bad in &cases {
        std::fs::write(&tbi, recompress(bad)).unwrap();
        let r = open_reader(bed.to_str().unwrap(), &o).unwrap();
        let Ok(ix) = IndexedFile::open(&bed, Format::Bed, r.header().clone(), &o) else {
            continue;
        };
        let regions: Regions = vec![(b"c1".to_vec(), Interval { start: 0, end: 100 })];
        if let Ok(mut rd) = ix.regions_reader(&regions, false) {
            let mut chroms = gloci_core::ChromTable::first_appearance();
            let mut batch = RecordBatch::default();
            while let Ok(true) = rd.read_batch(&mut chroms, &mut batch, 64) {
                batch.clear();
            }
        }
    }
    assert!(
        started.elapsed() < std::time::Duration::from_secs(60),
        "{} corrupt indexes took {:?}",
        cases.len(),
        started.elapsed()
    );
}

/// An index of another BAM (other reference count) or one whose offsets lie beyond the file is
/// refused with a message that says to rebuild it.
#[test]
fn index_of_another_file_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    // flags.bam has two references; give it the index of a one-reference BAM
    let p = dir.path().join("x.bam");
    std::fs::copy(fixture("flags.bam"), &p).unwrap();
    std::fs::copy(fixture("cigars.bam.bai"), dir.path().join("x.bam.bai")).unwrap();
    let r = open_reader(p.to_str().unwrap(), &ropts()).unwrap();
    let e = IndexedFile::open(&p, Format::Bam, r.header().clone(), &ropts())
        .err()
        .expect("an index of another BAM is refused");
    assert!(e.to_string().contains("belongs to another file"), "{e}");
    // the offsets of a long BAM's index point past the end of a shorter file
    let big = tempfile::tempdir().unwrap();
    let bq = big.path().join("y.bam");
    std::fs::copy(fixture("cat-big.bam"), &bq).unwrap();
    write_alignment_index(&bq, &bq, false, threads()).unwrap();
    let short = dir.path().join("z.bam");
    let bytes = std::fs::read(&bq).unwrap();
    std::fs::copy(big.path().join("y.bam.bai"), dir.path().join("z.bam.bai")).unwrap();
    // half of the file plus an EOF block: a complete BGZF file, shorter than the index says
    let mut complete = bytes[..bytes.len() / 2].to_vec();
    complete.extend_from_slice(&BGZF_EOF);
    std::fs::write(&short, complete).unwrap();
    let r = open_reader(short.to_str().unwrap(), &ropts()).unwrap();
    let e = IndexedFile::open(&short, Format::Bam, r.header().clone(), &ropts())
        .err()
        .expect("offsets past the end are refused");
    assert!(e.to_string().contains("beyond the end"), "{e}");
}
