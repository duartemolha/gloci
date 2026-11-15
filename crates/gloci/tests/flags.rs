//! `gloci flags` (D136) and the D100 gate lists of the samtools-style commands (D125).
mod common;
use common::*;
use gloci_testkit::fixtures;
use predicates::prelude::PredicateBooleanExt;
use predicates::str::contains;

fn out(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success().get_output().clone();
    String::from_utf8(o.stdout).unwrap()
}

#[test]
fn flags_prints_number_and_names() {
    assert_eq!(
        out(&[
            "flags",
            "0x704",
            "1796",
            "secondary,dup",
            "UNMAP,MUNMAP",
            "0",
            "4096"
        ]),
        "0x704\t1796\tUNMAP,SECONDARY,QCFAIL,DUP\n\
         0x704\t1796\tUNMAP,SECONDARY,QCFAIL,DUP\n\
         0x500\t1280\tSECONDARY,DUP\n\
         0xc\t12\tUNMAP,MUNMAP\n\
         0x0\t0\t\n\
         0x1000\t4096\t\n"
    );
    assert_eq!(
        out(&["flags", "Supplementary"]),
        "0x800\t2048\tSUPPLEMENTARY\n"
    );
}

#[test]
fn flags_table_is_samtools_layout() {
    let t = out(&["flags"]);
    let lines: Vec<&str> = t.lines().collect();
    assert_eq!(lines.len(), 12);
    assert_eq!(
        lines[0],
        "   0x1     1  PAIRED         paired-end / multiple-segment sequencing technology"
    );
    assert_eq!(
        lines[11],
        " 0x800  2048  SUPPLEMENTARY  supplementary alignment"
    );
    assert_eq!(
        lines[4],
        "  0x10    16  REVERSE        SEQ is reverse complemented"
    );
}

#[test]
fn flags_usage_errors_exit_2_and_print_nothing() {
    gloci()
        .args(["flags", "1", "010"])
        .assert()
        .code(2)
        .stdout("")
        .stderr(contains("leading zero").and(contains("0x8")));
    gloci()
        .args(["flags", "bogus"])
        .assert()
        .code(2)
        .stderr(contains("PAIRED"));
    for bad in ["65536", "12abc", "-1", "0x", ""] {
        gloci().args(["flags", bad]).assert().code(2);
    }
    gloci().args(["flags", "0xFFFF"]).assert().success();
}

#[test]
fn gate_lists_p3c_commands() {
    let dir = tempdir();
    let bam = fixtures::write(
        dir.path(),
        "x.bam",
        &fixtures::bgzip(b"BAM\x01\0\0\0\0\0\0\0\0"),
    );
    let fa = plain_file(dir.path(), "ref.fa", ">chr1 d\nACGT\n");
    gloci()
        .args(["merge", bam.to_str().unwrap()])
        .assert()
        .code(2)
        .stderr(contains("depth").and(contains("flagstat")));
    gloci()
        .args(["sort", fa.to_str().unwrap()])
        .assert()
        .code(2)
        .stderr(contains("faidx").and(contains("dict")));
}
