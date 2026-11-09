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
