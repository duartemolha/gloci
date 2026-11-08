use gloci_core::{ChromTable, Strand, VcfSpan};
use gloci_io::*;
use gloci_testkit::fixtures;
use std::io::Cursor;
use std::num::NonZeroUsize;

fn strands(bytes: &[u8], format: Format) -> Vec<Strand> {
    let o = ReaderOptions {
        threads: NonZeroUsize::new(1).unwrap(),
        format: Some(format),
        vcf_span: VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    };
    let mut r = open_reader_from("x", Box::new(Cursor::new(bytes.to_vec())), &o).unwrap();
    let mut chroms = ChromTable::lexicographic();
    let mut batch = RecordBatch::default();
    let mut out = Vec::new();
    loop {
        batch.clear();
        let more = r.read_batch(&mut chroms, &mut batch, 16).unwrap();
        out.extend(batch.recs.iter().map(|m| m.strand));
        if !more {
            return out;
        }
    }
}

fn bed(line: &str) -> Strand {
    strands(format!("{line}\n").as_bytes(), Format::Bed)[0]
}

#[test]
fn bed_column_six() {
    assert_eq!(bed("chr1\t0\t5\tn\t0\t+"), Strand::Plus);
    assert_eq!(bed("chr1\t0\t5\tn\t0\t-"), Strand::Minus);
    assert_eq!(bed("chr1\t0\t5\tn\t0\t."), Strand::Unknown);
    assert_eq!(bed("chr1\t0\t5\tn"), Strand::Unknown);
    assert_eq!(bed("chr1\t0\t5\tn\t0\tx"), Strand::Unknown);
    assert_eq!(bed("chr1\t0\t5\tn\t0\t?"), Strand::Unknown);
    assert_eq!(bed("chr1\t0\t5\tn\t0\t+\t1\t2"), Strand::Plus);
}

#[test]
fn bedgraph_is_unknown() {
    assert_eq!(bed("chr1\t0\t5\t1.5"), Strand::Unknown);
}

#[test]
fn gff_column_seven() {
    let g = |s: &str| {
        strands(
            format!("chr1\tsrc\tgene\t1\t5\t.\t{s}\t.\tID=a\n").as_bytes(),
            Format::Gff3,
        )[0]
    };
    assert_eq!(g("+"), Strand::Plus);
    assert_eq!(g("-"), Strand::Minus);
    assert_eq!(g("?"), Strand::Unknown);
    assert_eq!(g("."), Strand::Unknown);
}

#[test]
fn vcf_and_bcf_are_unknown() {
    let vcf = b"##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t5\t.\tA\tC\t.\t.\t.\n";
    assert_eq!(strands(vcf, Format::Vcf), vec![Strand::Unknown]);
    assert_eq!(
        strands(&fixtures::vcf_to_bcf(vcf), Format::Bcf),
        vec![Strand::Unknown]
    );
}

#[test]
fn field_helpers() {
    assert_eq!(field(b"a\tb\tc", 1), Some(&b"b"[..]));
    assert_eq!(field(b"a\tb\tc", 3), None);
    assert_eq!(field(b"a\t\tc", 1), Some(&b""[..]));
    assert_eq!(column_count(b"a\tb\tc"), 3);
    assert_eq!(column_count(b""), 1);
    assert_eq!(
        split_fields(b"a\tb\t").collect::<Vec<_>>(),
        vec![&b"a"[..], b"b", b""]
    );
}
