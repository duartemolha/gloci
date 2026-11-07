//! Index build (TBI/CSI, D20), index discovery (D26) and the index-backed
//! region reader (§4.3, D9, D10).

use gloci_core::{ChromTable, Interval, VcfSpan};
use gloci_io::*;
use gloci_testkit::strategies::{self, TestRec};
use gloci_testkit::{fixtures, naive};
use noodles_csi::BinningIndex;
use proptest::prelude::*;
use std::io::Write;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::sync::atomic::{AtomicU64, Ordering};

const BED_SMALL: &[u8] = b"chr1\t10\t20\ta\nchr1\t15\t30\tb\nchr2\t5\t6\tc\n";

const HEADER_ONLY: &[u8] = b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=1000000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// SNV, insertion, <DEL> with END, <DEL> with SVLEN only, a deletion; on
/// two contigs, sorted in header order.
const VCF_FIXTURE: &[u8] = b"##fileformat=VCFv4.3\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n\
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"SV length\">\n\
##FILTER=<ID=q10,Description=\"Quality below 10\">\n\
##ALT=<ID=DEL,Description=\"Deletion\">\n\
##contig=<ID=chr1,length=1000000>\n\
##contig=<ID=chr2,length=1000000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t.\tA\tG\t.\tPASS\t.\n\
chr1\t200\t.\tA\tACGT\t.\tPASS\t.\n\
chr1\t300\t.\tN\t<DEL>\t.\tPASS\tEND=40000\n\
chr1\t1000\t.\tN\t<DEL>\t.\tPASS\tSVLEN=-50\n\
chr2\t50\t.\tACGT\tA\t.\tq10\t.\n";

fn nz(n: usize) -> NonZeroUsize {
    NonZeroUsize::new(n).unwrap()
}

fn iv(s: u64, e: u64) -> Interval {
    Interval::new(s, e).unwrap()
}

/// A fresh path in Cargo's scratch directory for integration tests. Unlike
/// a `TempDir`, it stays valid after the helper that made it returns; the
/// directory is emptied once per run.
fn tmp(name: &str) -> PathBuf {
    static DIR: OnceLock<PathBuf> = OnceLock::new();
    static N: AtomicU64 = AtomicU64::new(0);
    let dir = DIR.get_or_init(|| {
        let d = Path::new(env!("CARGO_TARGET_TMPDIR")).join("gloci-io-index");
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    });
    dir.join(format!("{}-{name}", N.fetch_add(1, Ordering::Relaxed)))
}

fn bgzf_file(name: &str, data: &[u8]) -> PathBuf {
    let p = tmp(name);
    std::fs::write(&p, fixtures::bgzip(data)).unwrap();
    p
}

fn bgzf_bed(data: &[u8]) -> PathBuf {
    bgzf_file("x.bed.gz", data)
}

fn bgzf_vcf(data: &[u8]) -> PathBuf {
    bgzf_file("x.vcf.gz", data)
}

fn plain_bed() -> PathBuf {
    let p = tmp("x.bed");
    std::fs::write(&p, BED_SMALL).unwrap();
    p
}

/// BGZF with many small blocks: a block boundary after every line, and in
/// the middle of every third line, so records straddle blocks.
fn write_bgzf_bed(recs: &[TestRec]) -> PathBuf {
    let mut w = noodles_bgzf::io::Writer::new(Vec::new());
    for (i, r) in recs.iter().enumerate() {
        let line = r.bed_line() + "\n";
        let cut = if i % 3 == 1 { line.len() / 2 } else { 0 };
        w.write_all(&line.as_bytes()[..cut]).unwrap();
        w.flush().unwrap();
        w.write_all(&line.as_bytes()[cut..]).unwrap();
        w.flush().unwrap();
    }
    let p = tmp("p.bed.gz");
    std::fs::write(&p, w.finish().unwrap()).unwrap();
    p
}

/// Every record of `r`, read in small batches.
fn lines(mut r: Box<dyn RecordReader>) -> Vec<String> {
    let mut chroms = ChromTable::first_appearance();
    let mut batch = RecordBatch::default();
    let mut out = Vec::new();
    loop {
        batch.clear();
        let more = r.read_batch(&mut chroms, &mut batch, 3).unwrap();
        assert!(batch.len() <= 3);
        for i in 0..batch.len() {
            out.push(String::from_utf8(batch.raw(i).to_vec()).unwrap());
        }
        if !more {
            break;
        }
    }
    out
}

fn opts(span: VcfSpan) -> ReaderOptions {
    ReaderOptions {
        threads: nz(1),
        format: None,
        vcf_span: span,
        invalid: InvalidPolicy::Fail,
    }
}

fn bed_opts() -> ReaderOptions {
    opts(VcfSpan::Allele)
}

/// Rank of a strategy chromosome in lexicographic order.
fn lexicographic_rank(chrom: &str) -> u32 {
    let mut names = strategies::CHROMS.to_vec();
    names.sort_unstable();
    names
        .iter()
        .position(|c| *c == chrom)
        .expect("a chromosome from the strategy pool") as u32
}

fn write_and_find_at(p: &Path, f: Format) -> IndexInfo {
    let written = write_index(p, p, f, nz(2)).unwrap();
    let info = find_index(p).unwrap();
    assert_eq!(info.path, written);
    assert!(!info.stale);
    info
}

fn write_and_find(data: &[u8], f: Format) -> IndexInfo {
    write_and_find_at(&bgzf_file("w.gz", data), f)
}

fn write_and_find_bcf(vcf: &[u8]) -> IndexInfo {
    let p = tmp("w.bcf");
    std::fs::write(&p, fixtures::vcf_to_bcf(vcf)).unwrap();
    write_and_find_at(&p, Format::Bcf)
}

fn indexed(bed: &[u8]) -> PathBuf {
    let p = bgzf_bed(bed);
    write_index(&p, &p, Format::Bed, nz(1)).unwrap();
    p
}

fn indexed_vcf(records: &str) -> PathBuf {
    let mut text = HEADER_ONLY.to_vec();
    text.extend_from_slice(records.as_bytes());
    let p = bgzf_vcf(&text);
    write_index(&p, &p, Format::Vcf, nz(1)).unwrap();
    p
}

fn query(p: &Path, regions: &[(&str, Interval)], o: &ReaderOptions, incl: bool) -> Vec<String> {
    let regions: Vec<(Vec<u8>, Interval)> = regions
        .iter()
        .map(|(c, iv)| (c.as_bytes().to_vec(), *iv))
        .collect();
    lines(open_region_reader(p, &regions, o, incl).unwrap())
}

fn count(p: &Path, chrom: &str, iv: Interval, incl: bool) -> usize {
    query(p, &[(chrom, iv)], &bed_opts(), incl).len()
}

fn count_vcf(p: &Path, iv: Interval, span: VcfSpan, incl: bool) -> usize {
    query(p, &[("chr1", iv)], &opts(span), incl).len()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]
    #[test]
    fn region_reader_matches_naive(recs in strategies::bed_records(), regions in strategies::regions(4)) {
        let sorted = naive::sort(recs, lexicographic_rank);
        let p = write_bgzf_bed(&sorted); write_index(&p, &p, Format::Bed, nz(2)).unwrap();
        for incl in [false, true] {
            let got = lines(open_region_reader(&p, &regions, &bed_opts(), incl).unwrap());
            // filtering the sorted list keeps each record once, in file order
            let want: Vec<_> = sorted.iter().filter(|r| regions.iter().any(|(c, iv)| naive::matches(r, c, *iv, incl))).map(TestRec::bed_line).collect();
            prop_assert_eq!(got, want);
        }
    }
}

#[test]
fn kinds() {
    assert_eq!(write_and_find(BED_SMALL, Format::Bed).kind, IndexKind::Tbi);
    assert_eq!(
        write_and_find(b"chr1\t600000000\t600000010\n", Format::Bed).kind,
        IndexKind::Csi
    ); // > 2^29
    assert_eq!(write_and_find_bcf(VCF_FIXTURE).kind, IndexKind::Csi);
    noodles_tabix::fs::read(write_and_find(BED_SMALL, Format::Bed).path).unwrap(); // interop
}

#[test]
fn zero_length_bed_record() {
    // the user's insertion-site case through an index
    let p = indexed(b"chr1\t10\t10\tins\n");
    assert_eq!(count(&p, "chr1", iv(5, 15), false), 1);
    assert_eq!(count(&p, "chr1", iv(10, 15), false), 0);
    assert_eq!(count(&p, "chr1", iv(10, 15), true), 1);
}

#[test]
fn vcf_insertion_at_region_edge() {
    // D9
    let p = indexed_vcf("chr1\t200\t.\tA\tACGT\t.\t.\t.\n"); // point [200,200)
    assert_eq!(count_vcf(&p, iv(99, 200), VcfSpan::Allele, false), 0);
    assert_eq!(count_vcf(&p, iv(99, 200), VcfSpan::Allele, true), 1);
    assert_eq!(count_vcf(&p, iv(99, 200), VcfSpan::Ref, false), 1);
}

#[test]
fn unsorted_and_non_bgzf_are_refused() {
    let p = bgzf_bed(b"chr1\t5\t6\nchr1\t1\t2\n");
    assert!(matches!(
        write_index(&p, &p, Format::Bed, nz(1)),
        Err(IoError::Order(_))
    ));
    let p = plain_bed();
    assert!(
        write_index(&p, &p, Format::Bed, nz(1))
            .unwrap_err()
            .to_string()
            .contains("not BGZF-compressed")
    );
}

#[test]
fn empty_inputs_index_fine() {
    // Review Focus 1
    let p = bgzf_bed(b"");
    write_index(&p, &p, Format::Bed, nz(1)).unwrap();
    assert_eq!(count(&p, "chr1", iv(0, 100), false), 0);
    let v = bgzf_vcf(HEADER_ONLY);
    write_index(&v, &v, Format::Vcf, nz(1)).unwrap();
}

#[test]
fn stale_index_flagged() {
    let p = indexed(BED_SMALL);
    assert!(!find_index(&p).unwrap().stale);
    let later = std::time::SystemTime::now() + std::time::Duration::from_secs(60);
    std::fs::File::options()
        .write(true)
        .open(&p)
        .unwrap()
        .set_modified(later)
        .unwrap();
    assert!(find_index(&p).unwrap().stale);
}

// ---- beyond the brief ----

#[test]
fn error_messages_name_the_file_and_fix() {
    let p = plain_bed();
    let name = p.display().to_string();
    let e = write_index(&p, &p, Format::Bed, nz(1)).unwrap_err();
    assert!(matches!(e, IoError::Usage(_)));
    assert_eq!(
        e.to_string(),
        format!("{name} is not BGZF-compressed; recompress with: gloci view {name} -o {name}.gz")
    );
    // the user-visible name is the final path, not the temp file read
    let shown = Path::new("final.bed.gz");
    let e = write_index(&p, shown, Format::Bed, nz(1)).unwrap_err();
    assert!(e.to_string().starts_with("final.bed.gz is not BGZF"), "{e}");
    // a chromosome that comes back after another one is unsorted too
    let p = bgzf_bed(b"chr1\t1\t2\nchr2\t1\t2\nchr1\t5\t6\n");
    let e = write_index(&p, &p, Format::Bed, nz(1)).unwrap_err();
    assert!(e.to_string().contains("line 3"), "{e}");
    assert!(!Path::new(&format!("{}.tbi", p.display())).exists());
}

#[test]
fn track_lines_are_skipped_by_the_tabix_header() {
    let p = indexed(b"#c\nbrowser position chr1\ntrack name=x\nchr1\t1\t5\tr\n");
    let idx = noodles_tabix::fs::read(find_index(&p).unwrap().path).unwrap();
    let h = idx.header().unwrap();
    assert_eq!(h.line_skip_count(), 3);
    assert_eq!(h.line_comment_prefix(), b'#');
    assert_eq!(count(&p, "chr1", iv(0, 10), false), 1);
}

#[test]
fn csi_for_large_coordinates_round_trips() {
    let big = 1u64 << 40;
    let bed = format!(
        "chr1\t5\t10\ta\nchr1\t{}\t{big}\tb\nchr2\t{big}\t{big}\tc\n",
        big - 100
    );
    let p = indexed(bed.as_bytes());
    let info = find_index(&p).unwrap();
    assert_eq!(info.kind, IndexKind::Csi);
    let idx = noodles_csi::fs::read(&info.path).unwrap();
    assert_eq!((idx.min_shift(), idx.depth()), (14, 9)); // 2^(14+27) >= 2^40
    assert_eq!(count(&p, "chr1", iv(big - 1, big), false), 1);
    assert_eq!(count(&p, "chr2", iv(big - 1, big + 1), false), 1);
    assert_eq!(count(&p, "chr2", iv(big, big + 1), false), 0);
    assert_eq!(count(&p, "chr1", iv(0, u64::MAX), false), 2);
}

#[test]
fn writing_one_kind_removes_the_other() {
    let p = bgzf_bed(BED_SMALL);
    let csi = PathBuf::from(format!("{}.csi", p.display()));
    std::fs::write(&csi, b"stale").unwrap();
    let tbi = write_index(&p, &p, Format::Bed, nz(1)).unwrap();
    assert_eq!(tbi, PathBuf::from(format!("{}.tbi", p.display())));
    assert!(tbi.exists() && !csi.exists());
    assert!(find_index(&bgzf_bed(BED_SMALL)).is_none());
}

#[test]
fn regions_are_merged_and_unknown_chromosomes_yield_nothing() {
    let p = indexed(BED_SMALL);
    let got = query(
        &p,
        &[
            ("chr2", iv(0, 100)),
            ("chr1", iv(16, 17)),
            ("chrUn", iv(0, 100)),
            ("chr1", iv(0, 12)),
            ("chr1", iv(0, 12)),
        ],
        &bed_opts(),
        false,
    );
    assert_eq!(got, ["chr1\t10\t20\ta", "chr1\t15\t30\tb", "chr2\t5\t6\tc"]);
    // a point between two bookended regions is in neither (no union of regions)
    let p = indexed(b"chr1\t10\t10\tins\n");
    assert_eq!(
        query(
            &p,
            &[("chr1", iv(5, 10)), ("chr1", iv(10, 15))],
            &bed_opts(),
            false
        )
        .len(),
        0
    );
}

#[test]
fn bcf_and_vcf_queries_agree() {
    let bcf = tmp("q.bcf");
    std::fs::write(&bcf, fixtures::vcf_to_bcf(VCF_FIXTURE)).unwrap();
    write_index(&bcf, &bcf, Format::Bcf, nz(2)).unwrap();
    let vcf = bgzf_vcf(VCF_FIXTURE);
    write_index(&vcf, &vcf, Format::Vcf, nz(2)).unwrap();
    let keys = |p: &Path, region: (&str, Interval), span: VcfSpan| -> Vec<(Vec<u8>, u64)> {
        let regions = [(region.0.as_bytes().to_vec(), region.1)];
        let mut r = open_region_reader(p, &regions, &opts(span), false).unwrap();
        assert!(matches!(r.header(), Header::Vcf(_)));
        let mut chroms = ChromTable::first_appearance();
        let mut batch = RecordBatch::default();
        let mut out = Vec::new();
        loop {
            batch.clear();
            let more = r.read_batch(&mut chroms, &mut batch, 2).unwrap();
            out.extend(
                batch
                    .recs
                    .iter()
                    .map(|m| (chroms.name(m.chrom).to_vec(), m.key)),
            );
            if !more {
                break;
            }
        }
        out
    };
    for span in [VcfSpan::Allele, VcfSpan::Ref] {
        for region in [
            ("chr1", iv(0, 1_000_000)),
            ("chr1", iv(199, 201)),
            ("chr1", iv(30_000, 30_001)), // inside <DEL> END=40000
            ("chr1", iv(1040, 1041)),     // inside <DEL> SVLEN=-50
            ("chr2", iv(50, 51)),
            ("chr3", iv(0, 10)),
        ] {
            assert_eq!(
                keys(&bcf, region, span),
                keys(&vcf, region, span),
                "{region:?} {span:?}"
            );
        }
    }
    assert_eq!(
        keys(&bcf, ("chr1", iv(30_000, 30_001)), VcfSpan::Ref),
        [(b"chr1".to_vec(), 299)]
    );
    // the SVLEN-only <DEL> [999, 1050), and the END=40000 one [299, 40000)
    assert_eq!(
        keys(&bcf, ("chr1", iv(1040, 1041)), VcfSpan::Ref),
        [(b"chr1".to_vec(), 299), (b"chr1".to_vec(), 999)]
    );
    assert_eq!(
        keys(&bcf, ("chr1", iv(1050, 1051)), VcfSpan::Ref),
        [(b"chr1".to_vec(), 299)]
    );
}

#[test]
fn bcf_out_of_dictionary_order_is_unsorted() {
    let swapped = std::str::from_utf8(VCF_FIXTURE).unwrap().replace(
        "##contig=<ID=chr1,length=1000000>\n##contig=<ID=chr2,length=1000000>\n",
        "##contig=<ID=chr2,length=1000000>\n##contig=<ID=chr1,length=1000000>\n",
    );
    let p = tmp("u.bcf");
    std::fs::write(&p, fixtures::vcf_to_bcf(swapped.as_bytes())).unwrap();
    let e = write_index(&p, &p, Format::Bcf, nz(1)).unwrap_err();
    assert!(matches!(e, IoError::Order(_)), "{e}");
    assert!(e.to_string().contains("record 5"), "{e}");
}

/// `(chromosome, POS)` of every record of `r` (BCF records are binary).
fn positions(mut r: Box<dyn RecordReader>) -> Vec<(String, u64)> {
    let mut chroms = ChromTable::first_appearance();
    let mut batch = RecordBatch::default();
    let mut out = Vec::new();
    loop {
        batch.clear();
        let more = r.read_batch(&mut chroms, &mut batch, 2).unwrap();
        for m in &batch.recs {
            let c = String::from_utf8(chroms.name(m.chrom).to_vec()).unwrap();
            out.push((c, m.key + 1));
        }
        if !more {
            return out;
        }
    }
}

#[test]
fn htslib_bcf_index_with_contigs_out_of_dictionary_order() {
    // bcftools indexes a BCF whose contigs are contiguous but not in
    // dictionary order (e.g. made from a sort -k1,1 VCF): here chr2 (id 1)
    // comes before chr1 (id 0). A query on both returns both, in file order.
    let vcf = b"##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=1000>\n\
##contig=<ID=chr2,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr2\t5\t.\tA\tC\t.\t.\t.\n\
chr2\t9\t.\tG\tT\t.\t.\t.\n\
chr1\t7\t.\tA\tG\t.\t.\t.\n";
    let bcf = fixtures::vcf_to_bcf(vcf);
    let p = tmp("any.bcf");
    std::fs::write(&p, &bcf).unwrap();
    std::fs::write(
        format!("{}.csi", p.display()),
        fixtures::bcf_csi_any_order(&bcf),
    )
    .unwrap();
    assert_eq!(find_index(&p).unwrap().kind, IndexKind::Csi);
    let o = opts(VcfSpan::Allele);
    let q = |regions: &[(&str, Interval)]| {
        let regions: Vec<(Vec<u8>, Interval)> = regions
            .iter()
            .map(|(c, iv)| (c.as_bytes().to_vec(), *iv))
            .collect();
        positions(open_region_reader(&p, &regions, &o, false).unwrap())
    };
    let all = iv(0, u64::MAX);
    let both = vec![
        ("chr2".to_string(), 5),
        ("chr2".to_string(), 9),
        ("chr1".to_string(), 7),
    ];
    assert_eq!(q(&[("chr1", all), ("chr2", all)]), both);
    assert_eq!(q(&[("chr2", all), ("chr1", all)]), both);
    assert_eq!(q(&[("chr1", all)]), [("chr1".to_string(), 7)]);
    assert_eq!(
        q(&[("chr1", iv(0, 100)), ("chr2", iv(8, 9))]),
        [("chr2".to_string(), 9), ("chr1".to_string(), 7)]
    );
}

#[test]
fn region_edges_on_index_window_boundaries() {
    // a ends exactly where the second 16 kb window starts, and b, the first
    // record to reach that window, sets its linear-index offset
    let p = indexed(b"chr1\t100\t16384\ta\nchr1\t16384\t16384\tb\nchr1\t20000\t20010\tc\n");
    let names = |regions: &[(&str, Interval)], incl| -> Vec<String> {
        query(&p, regions, &bed_opts(), incl)
            .iter()
            .map(|l| l.split('\t').nth(3).unwrap().to_string())
            .collect()
    };
    assert_eq!(names(&[("chr1", iv(16384, 16384))], false), ["b"]);
    assert_eq!(names(&[("chr1", iv(16384, 16384))], true), ["a", "b"]);
    assert_eq!(names(&[("chr1", iv(16000, 16384))], false), ["a"]);
    assert_eq!(names(&[("chr1", iv(16000, 16384))], true), ["a", "b"]);
    // a is bookended with [16384, 16385): never an overlap (§4.3)
    assert_eq!(names(&[("chr1", iv(16384, 16385))], true), ["b"]);
    assert!(names(&[("chr1", iv(16384, 16385))], false).is_empty());
    // an insertion whose index span [16383, 16384) ends at the boundary;
    // its allele point 16384 is on the left edge of [16384, 16390)
    let v = indexed_vcf("chr1\t16384\t.\tA\tAT\t.\t.\t.\nchr1\t16391\t.\tC\tG\t.\t.\t.\n");
    assert_eq!(count_vcf(&v, iv(16384, 16390), VcfSpan::Allele, true), 1);
    assert_eq!(count_vcf(&v, iv(16384, 16390), VcfSpan::Allele, false), 0);
    assert_eq!(count_vcf(&v, iv(16384, 16390), VcfSpan::Ref, false), 0);
    assert_eq!(count_vcf(&v, iv(16384, 16384), VcfSpan::Allele, false), 1);
}

/// Stored (level 0) BGZF, so compressed sizes follow the text and leaf bins
/// with over 64 kB of records are not folded into their parents.
fn bgzf_stored(data: &[u8]) -> PathBuf {
    use noodles_bgzf::io::writer::{Builder, CompressionLevel};
    let level = CompressionLevel::try_from(0u8).unwrap();
    let mut w = Builder::default()
        .set_compression_level(level)
        .build_from_writer(Vec::new());
    w.write_all(data).unwrap();
    let p = tmp("stored.bed.gz");
    std::fs::write(&p, w.finish().unwrap()).unwrap();
    p
}

#[test]
fn records_in_parent_bins_and_between_chunks_tbi_and_csi() {
    let pad = "x".repeat(90);
    let mut bed = String::from("chr1\t150\t151\tA\n");
    bed += "chr1\t16000\t17000\tR\n"; // crosses the first window edge: a parent bin
    // ~160 kB in each of the second and third 16 kb windows
    for i in 0..1600 {
        bed += &format!("chr1\t{}\t{}\tw1{pad}\n", 20_000 + i * 5, 20_001 + i * 5);
    }
    for i in 0..1600 {
        bed += &format!("chr1\t{}\t{}\tw2{pad}\n", 33_000 + i * 4, 33_001 + i * 4);
    }
    bed += "chr1\t40000\t40001\tQ\n"; // between two chunks of the first region's window
    bed += "chr1\t50000\t70000\tS\n"; // a parent bin again
    let regions = [
        ("chr1", iv(100, 200)),
        ("chr1", iv(16500, 16501)),
        ("chr1", iv(40000, 60000)),
    ];
    for (extra, kind) in [
        ("", IndexKind::Tbi),
        ("chr2\t1073741824\t1073741825\tbig\n", IndexKind::Csi),
    ] {
        let p = bgzf_stored(format!("{bed}{extra}").as_bytes());
        write_index(&p, &p, Format::Bed, nz(2)).unwrap();
        assert_eq!(find_index(&p).unwrap().kind, kind);
        let got: Vec<String> = query(&p, &regions, &bed_opts(), false)
            .iter()
            .map(|l| l.split('\t').nth(3).unwrap().to_string())
            .collect();
        assert_eq!(got, ["A", "R", "Q", "S"], "{kind:?}");
        // R alone: its window's leaf bin starts after R
        assert_eq!(count(&p, "chr1", iv(16500, 16501), false), 1, "{kind:?}");
    }
}

#[test]
fn index_bytes_do_not_depend_on_threads() {
    // write_bgzf_bed puts every line in its own BGZF block (and splits
    // some), so the 4-thread reader really decodes blocks in parallel
    for (big, kind) in [(false, IndexKind::Tbi), (true, IndexKind::Csi)] {
        let mut recs: Vec<TestRec> = (0..3000u64)
            .map(|i| TestRec {
                chrom: if i < 2000 { "chr1" } else { "chr2" }.to_string(),
                start: i * 37,
                end: i * 37 + 50 + (i % 7) * 1000,
                strand: '+',
                name: format!("r{i}"),
            })
            .collect();
        if big {
            recs.push(TestRec {
                chrom: "chr2".into(),
                start: 1 << 35,
                end: (1 << 35) + 10,
                strand: '-',
                name: "big".into(),
            });
        }
        let p = write_bgzf_bed(&recs);
        let bytes =
            |t: usize| std::fs::read(write_index(&p, &p, Format::Bed, nz(t)).unwrap()).unwrap();
        let one = bytes(1);
        assert_eq!(one, bytes(4), "{kind:?}");
        assert_eq!(find_index(&p).unwrap().kind, kind);
    }
    let bcf = tmp("t.bcf");
    std::fs::write(&bcf, fixtures::vcf_to_bcf(VCF_FIXTURE)).unwrap();
    let bytes =
        |t: usize| std::fs::read(write_index(&bcf, &bcf, Format::Bcf, nz(t)).unwrap()).unwrap();
    let one = bytes(1);
    assert_eq!(one, bytes(4), "BCF");
}
