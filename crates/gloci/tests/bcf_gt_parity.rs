//! Byte-level parity of native BCF GT encoding (D78) with bcftools 1.21: the
//! record bytes gloci writes for `bcf-gt-parity/input.vcf` equal the committed
//! `expected.bcf` (see the README beside them). The header differs by the
//! `##glociCommand` line, so only the records are compared.
mod common;
use common::*;
use std::io::Read;
use std::path::Path;

/// The record bytes of a BGZF BCF: everything after the magic, `l_text` and
/// the header text.
fn record_bytes(bcf: &[u8]) -> Vec<u8> {
    let mut all = Vec::new();
    noodles_bgzf::io::Reader::new(bcf)
        .read_to_end(&mut all)
        .unwrap();
    assert_eq!(&all[..5], b"BCF\x02\x02");
    let l_text = u32::from_le_bytes(all[5..9].try_into().unwrap()) as usize;
    all[9 + l_text..].to_vec()
}

#[test]
fn bcf_gt_records_equal_bcftools() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/bcf-gt-parity");
    let d = tempdir();
    let out = d.path().join("g.bcf");
    gloci()
        .args(["view", "-O", "bcf"])
        .arg(dir.join("input.vcf"))
        .arg("-o")
        .arg(&out)
        .assert()
        .success();
    let mine = record_bytes(&std::fs::read(&out).unwrap());
    let theirs = record_bytes(&std::fs::read(dir.join("expected.bcf")).unwrap());
    assert!(!theirs.is_empty());
    assert_eq!(
        mine, theirs,
        "gloci's BCF records differ from bcftools 1.21's"
    );
}

/// A BCF whose GT is 16-bit (an allele index of 63 or more) is written as text like any
/// other, as bcftools writes it. (It was an error while noodles rendered the text: M8.)
#[test]
fn a_16_bit_gt_is_written_as_text() {
    let dir = tempdir();
    let alts: Vec<String> = (1..=70).map(|i| format!("A{}", "C".repeat(i))).collect();
    let vcf = plain_file(
        dir.path(),
        "x.vcf",
        &format!(
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\nchr1\t10\t.\tA\t{}\t.\t.\t.\tGT\t1/65\n",
            alts.join(",")
        ),
    );
    let bcf = dir.path().join("x.bcf");
    gloci()
        .args(["view"])
        .arg(&vcf)
        .args(["-O", "bcf", "-o"])
        .arg(&bcf)
        .assert()
        .success();
    for (args, last) in [
        (vec!["view", "-O", "vcf"], "\tGT\t1/65"),
        (vec!["vcf", "filter", "-O", "vcf"], "\tGT\t1/65"),
        (vec!["vcf", "extract", "--format", "GT"], "\t1/65"),
        (vec!["convert"], ""),
    ] {
        let o = gloci().args(&args).arg(&bcf).output().unwrap();
        let err = String::from_utf8_lossy(&o.stderr);
        assert_eq!(o.status.code(), Some(0), "{args:?}: {err}");
        let out = String::from_utf8_lossy(&o.stdout);
        let line = out.lines().last().unwrap_or_default();
        assert!(
            line.starts_with("chr1\t") && line.ends_with(last),
            "{args:?}: {line}"
        );
    }
}
