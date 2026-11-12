//! `gloci index` on FASTA (D108) and `-g ref.fa`.
mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;
use std::path::{Path, PathBuf};

fn fa(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/fa")
        .join(name)
}

fn copy(dir: &Path, name: &str) -> PathBuf {
    let p = dir.join(name);
    std::fs::copy(fa(name), &p).unwrap();
    p
}

#[test]
fn index_fasta_is_samtools_identical_at_any_thread_count() {
    for name in ["ref.fa", "crlf.fa", "ref.fa.gz"] {
        for t in ["1", "2", "8"] {
            let d = tempdir();
            let p = copy(d.path(), name);
            gloci()
                .args(["index", "-t", t])
                .arg(&p)
                .assert()
                .success()
                .stdout("");
            let fai = std::fs::read(format!("{}.fai", p.display())).unwrap();
            assert_eq!(
                fai,
                std::fs::read(fa(&format!("{name}.fai"))).unwrap(),
                "{name} -t {t}"
            );
            if name.ends_with(".gz") {
                let gzi = std::fs::read(format!("{}.gzi", p.display())).unwrap();
                assert_eq!(
                    gzi,
                    std::fs::read(fa(&format!("{name}.gzi"))).unwrap(),
                    "{name} -t {t}"
                );
            }
        }
    }
}

#[test]
fn index_fasta_errors_and_exit_codes() {
    let d = tempdir();
    let irr = copy(d.path(), "irregular.fa");
    gloci()
        .arg("index")
        .arg(&irr)
        .assert()
        .code(1)
        .stderr(contains("irregular.fa line 4").and(contains("chr1")));
    // Plain gzip: exit 2 with the bgzip recipe.
    let plain = std::fs::read(fa("ref.fa")).unwrap();
    let gz =
        gloci_testkit::fixtures::write(d.path(), "p.fa.gz", &gloci_testkit::fixtures::gzip(&plain));
    gloci()
        .arg("index")
        .arg(&gz)
        .assert()
        .code(2)
        .stderr(contains("bgzip"));
    // Duplicate names: exit 1.
    let dup = plain_file(d.path(), "dup.fa", ">a\nAC\n>a\nGG\n");
    gloci()
        .arg("index")
        .arg(&dup)
        .assert()
        .code(1)
        .stderr(contains("more than once"));
}

#[test]
fn genome_from_fasta() {
    let d = tempdir();
    let p = copy(d.path(), "ref.fa");
    let bed = plain_file(
        d.path(),
        "a.bed",
        "chrX\t5\t9\nchrM\t1\t2\nchr10\t3\t4\nchr1\t7\t8\nchr2\t1\t2\n",
    );
    // No .fai: exit 2, naming the fix. Nothing is built implicitly.
    gloci()
        .args(["sort", "-g"])
        .arg(&p)
        .arg(&bed)
        .assert()
        .code(2)
        .stderr(contains("gloci index").and(contains("ref.fa")));
    assert!(!d.path().join("ref.fa.fai").exists());
    gloci().arg("index").arg(&p).assert().success();
    let by_fa = gloci()
        .args(["sort", "-g"])
        .arg(&p)
        .arg(&bed)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let fai = d.path().join("ref.fa.fai");
    let by_fai = gloci()
        .args(["sort", "-g"])
        .arg(&fai)
        .arg(&bed)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    assert_eq!(by_fa, by_fai);
    // File order of the FASTA: chrM, chr1, chr2, chr10, chrX (not karyotype, not lexicographic).
    assert_eq!(
        String::from_utf8(by_fa).unwrap(),
        "chrM\t1\t2\nchr1\t7\t8\nchr2\t1\t2\nchr10\t3\t4\nchrX\t5\t9\n"
    );
}

#[test]
fn gate_message_names_the_fasta_uses() {
    let d = tempdir();
    let p = copy(d.path(), "ref.fa");
    gloci()
        .arg("sort")
        .arg(&p)
        .assert()
        .code(2)
        .stderr(contains("getfasta, maskfasta, nuc, faidx or dict").and(contains("-g ref.fa")));
}
