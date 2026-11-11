use gloci_core::{ChromOrder, ChromTable, Genome, OrderSource};
use gloci_io::*;
use gloci_testkit::fixtures;
use std::io::Cursor;
use std::num::NonZeroUsize;
use std::path::PathBuf;

const DICT: &[u8] =
    b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:248956422\n@SQ\tSN:chr2\tLN:242193529\n";
const VCF_WITH_CONTIGS: &[u8] = b"##fileformat=VCFv4.3\n##contig=<ID=chr2,length=1000,IDX=0>\n##contig=<ID=chr1,length=2000,IDX=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

fn names(g: Genome) -> Vec<String> {
    (0..g.len())
        .map(|i| String::from_utf8(g.name(i).to_vec()).unwrap())
        .collect()
}

/// Writes into a leaked temp dir (tests are short-lived).
fn write(name: &str, bytes: &[u8]) -> PathBuf {
    let dir = tempfile::tempdir().unwrap().keep();
    fixtures::write(&dir, name, bytes)
}

fn opts(format: Format) -> ReaderOptions {
    ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: None,
        threads: NonZeroUsize::new(1).unwrap(),
        format: Some(format),
        vcf_span: gloci_core::VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    }
}

fn reader(name: &str, bytes: &[u8], f: Format) -> Box<dyn RecordReader> {
    open_reader_from(name, Box::new(Cursor::new(bytes.to_vec())), &opts(f)).unwrap()
}

fn vcf_reader() -> Box<dyn RecordReader> {
    reader("h.vcf", VCF_WITH_CONTIGS, Format::Vcf)
}
fn bed_reader() -> Box<dyn RecordReader> {
    reader("x.bed", b"chr1\t1\t2\n", Format::Bed)
}
fn vcf_without_contigs() -> Box<dyn RecordReader> {
    reader(
        "n.vcf",
        b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        Format::Vcf,
    )
}

fn guard_err_karyotype(bed: &[u8]) -> IoError {
    guard_err_with(bed, ChromTable::karyotype())
}

fn guard_err_lex(bed: &[u8]) -> IoError {
    guard_err_with(bed, ChromTable::lexicographic())
}

fn guard_err_with(bed: &[u8], mut chroms: ChromTable) -> IoError {
    let mut r = reader("x.bed", bed, Format::Bed);
    let mut guard = OrderGuard::new("x.bed", Format::Bed);
    let mut batch = RecordBatch::default();
    loop {
        batch.clear();
        let more = r.read_batch(&mut chroms, &mut batch, 1000).unwrap();
        if let Err(e) = guard.check(&chroms, &batch) {
            return e;
        }
        if !more {
            panic!("no violation");
        }
    }
}

#[test]
fn load_genome_by_content() {
    assert_eq!(
        names(
            load_genome(&write(
                "g.fai",
                b"chr1\t248956422\t52\t60\t61\nchr2\t242193529\t253404903\t60\t61\n"
            ))
            .unwrap()
        ),
        ["chr1", "chr2"]
    );
    assert_eq!(
        names(load_genome(&write("g.genome", b"chrA\t10\n")).unwrap()),
        ["chrA"]
    );
    assert_eq!(
        load_genome(&write("g.dict", DICT)).unwrap().length(0),
        Some(248956422)
    );
    assert_eq!(
        names(load_genome(&write("h.vcf.gz", &fixtures::bgzip(VCF_WITH_CONTIGS))).unwrap()),
        ["chr2", "chr1"]
    );
    let g = load_genome(&write("h.bcf", &fixtures::vcf_to_bcf(VCF_WITH_CONTIGS))).unwrap();
    assert_eq!(g.length(1), Some(2000));
    assert_eq!(names(g), ["chr2", "chr1"]);
    let g = load_genome(&write("h.bam", &fixtures::bam_header(&[("chrM", 16569)]))).unwrap();
    assert_eq!(g.length(0), Some(16569));
    assert_eq!(names(g), ["chrM"]);
    let e = load_genome(&write("bad.genome", b"chr1\tlong\n"))
        .unwrap_err()
        .to_string();
    assert!(e.contains("bad.genome line 1"), "{e}");
}

#[test]
fn bcf_header_order_is_dictionary_order() {
    let vcf: &[u8] = b"##fileformat=VCFv4.3\n##contig=<ID=chr1,length=2000,IDX=1>\n##contig=<ID=chr2,length=1000,IDX=0>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let bcf = fixtures::vcf_to_bcf(vcf);
    assert_eq!(
        names(load_genome(&write("o.bcf", &bcf)).unwrap()),
        ["chr2", "chr1"]
    );
    let r = reader("o.bcf", &bcf, Format::Bcf);
    let t = chrom_table_for(None, &[&*r], ChromOrder::Karyotype).unwrap();
    assert_eq!(t.get(b"chr2"), Some(gloci_core::ChromId(0)));
}

#[test]
fn precedence() {
    let g = write("g.genome", b"chrA\t10\n");
    let (v, b, n) = (vcf_reader(), bed_reader(), vcf_without_contigs());
    assert_eq!(
        *chrom_table_for(Some(g.to_str().unwrap()), &[&*v], ChromOrder::Karyotype)
            .unwrap()
            .order_source(),
        OrderSource::GenomeFile(g.clone())
    );
    assert!(matches!(
        chrom_table_for(None, &[&*v], ChromOrder::Lex)
            .unwrap()
            .order_source(),
        OrderSource::VcfHeader(_)
    ));
    assert_eq!(
        *chrom_table_for(None, &[&*b], ChromOrder::Karyotype)
            .unwrap()
            .order_source(),
        OrderSource::Karyotype
    );
    assert_eq!(
        *chrom_table_for(None, &[&*b], ChromOrder::Lex)
            .unwrap()
            .order_source(),
        OrderSource::Lexicographic
    );
    assert_eq!(
        *chrom_table_for(None, &[&*n], ChromOrder::Karyotype)
            .unwrap()
            .order_source(),
        OrderSource::Karyotype
    );
}

#[test]
fn guard_message() {
    let e = guard_err_lex(b"chr2\t5\t9\nchr10\t1\t2\n").to_string();
    assert_eq!(
        e,
        "x.bed is not sorted: line 2 (chr10\t1\t2) comes after line 1 (chr2\t5\t9)\n\
                   hint: run: gloci sort x.bed --chrom-order lex"
    );
}

#[test]
fn guard_karyotype_explains_old_sort_order() {
    let e = guard_err_karyotype(b"chr1\t5\t9\nchr10\t1\t2\nchr2\t0\t3\n").to_string();
    assert!(
        e.contains("line 3 (chr2\t0\t3) comes after line 2 (chr10\t1\t2)"),
        "{e}"
    );
    assert!(
        e.contains("looks sorted by sort -k1,1, sort-bed or bedtools sort")
            && e.contains("--chrom-order lex")
            && e.contains("docs/migration.md"),
        "{e}"
    );
    // Not a lexicographic pair: plain hint, no note.
    let e = guard_err_karyotype(b"chr2\t5\t9\nchr1\t1\t2\n").to_string();
    assert!(e.ends_with("hint: run: gloci sort x.bed"), "{e}");
    // Start-backwards never gets the note.
    let e = guard_err_karyotype(b"chr1\t5\t9\nchr1\t1\t2\n").to_string();
    assert!(e.ends_with("hint: run: gloci sort x.bed"), "{e}");
}

#[test]
fn guard_start_backwards_across_batches_and_hint() {
    let mut r = reader("x.bed", b"chr1\t50\t60\nchr1\t5\t9\n", Format::Bed);
    let mut chroms = ChromTable::karyotype();
    let mut guard = OrderGuard::new("x.bed", Format::Bed);
    let mut batch = RecordBatch::default();
    batch.clear();
    r.read_batch(&mut chroms, &mut batch, 1).unwrap();
    guard.check(&chroms, &batch).unwrap();
    batch.clear();
    r.read_batch(&mut chroms, &mut batch, 1).unwrap();
    let e = guard.check(&chroms, &batch).unwrap_err().to_string();
    assert_eq!(
        e,
        "x.bed is not sorted: line 2 (chr1\t5\t9) comes after line 1 (chr1\t50\t60)\nhint: run: gloci sort x.bed"
    );
}

/// Run a contiguity guard over `bed` read with `chroms`, in batches of `max`.
fn contiguous(bed: &[u8], mut chroms: ChromTable, max: usize) -> Result<(), IoError> {
    let mut r = reader("x.bed", bed, Format::Bed);
    let mut guard = OrderGuard::contiguous("x.bed", Format::Bed);
    let mut batch = RecordBatch::default();
    loop {
        batch.clear();
        let more = r.read_batch(&mut chroms, &mut batch, max).unwrap();
        guard.check(&chroms, &batch)?;
        if !more {
            return Ok(());
        }
    }
}

#[test]
fn contiguity_guard() {
    // any chromosome order, whatever the table's order
    let lex = b"chr1\t5\t9\nchr10\t1\t2\nchr10\t1\t5\nchr2\t0\t3\n";
    for max in [1, 1000] {
        contiguous(lex, ChromTable::karyotype(), max).unwrap();
        contiguous(
            b"chrX\t1\t2\nchr1\t1\t2\n",
            ChromTable::lexicographic(),
            max,
        )
        .unwrap();
    }
    // a chromosome that comes back, across batches too, and no lex note
    for max in [1, 1000] {
        let e = contiguous(
            b"chr1\t5\t9\nchr2\t1\t2\nchr1\t7\t8\n",
            ChromTable::karyotype(),
            max,
        )
        .unwrap_err()
        .to_string();
        assert_eq!(
            e,
            "x.bed is not sorted: line 3 (chr1\t7\t8) comes after line 2 (chr2\t1\t2)\nhint: run: gloci sort x.bed"
        );
    }
    let e = contiguous(
        b"chr1\t50\t60\nchr1\t5\t9\n",
        ChromTable::lexicographic(),
        1,
    )
    .unwrap_err()
    .to_string();
    assert!(
        e.contains("line 2 (chr1\t5\t9) comes after line 1 (chr1\t50\t60)"),
        "{e}"
    );
    assert!(
        e.ends_with("hint: run: gloci sort x.bed --chrom-order lex"),
        "{e}"
    );
}

#[test]
fn sort_hints() {
    let p = PathBuf::from("g.fai");
    assert_eq!(
        sort_hint("a.bed", &OrderSource::GenomeFile(p)),
        "run: gloci sort a.bed -g g.fai"
    );
    assert_eq!(
        sort_hint("a.bed", &OrderSource::VcfHeader("v".into())),
        "run: gloci sort a.bed"
    );
    assert_eq!(
        sort_hint("a.bed", &OrderSource::Lexicographic),
        "run: gloci sort a.bed --chrom-order lex"
    );
    assert_eq!(
        sort_hint("a.bed", &OrderSource::Karyotype),
        "run: gloci sort a.bed"
    );
}

#[test]
fn bcf_guard_renders_chrom_pos() {
    let vcf: Vec<u8> = [
        VCF_WITH_CONTIGS,
        b"chr1\t500\t.\tA\tG\t.\t.\t.\nchr2\t100\t.\tA\tG\t.\t.\t.\n",
    ]
    .concat();
    let bcf = fixtures::vcf_to_bcf(&vcf);
    let mut r = reader("c.bcf", &bcf, Format::Bcf);
    let mut chroms = chrom_table_for(None, &[&*r], ChromOrder::Karyotype).unwrap();
    let mut guard = OrderGuard::new("c.bcf", Format::Bcf);
    let mut batch = RecordBatch::default();
    r.read_batch(&mut chroms, &mut batch, 10).unwrap();
    let e = guard.check(&chroms, &batch).unwrap_err().to_string();
    assert_eq!(
        e,
        "c.bcf is not sorted: record 2 (chr2 100) comes after record 1 (chr1 500)\nhint: run: gloci sort c.bcf"
    );
}
