//! `gloci flagstat` (D128): layouts, format inference, CRAM, thread counts; the byte-level
//! comparison with samtools is in the goldens (`flagstat-*-vs-samtools`).
mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;
use std::path::{Path, PathBuf};

fn inputs() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs")
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn out(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success().get_output().clone();
    String::from_utf8(o.stdout).unwrap()
}

fn fixture() -> String {
    s(&inputs().join("samtools/flagstat.bam"))
}

#[test]
fn flagstat_text_has_samtools_percentages() {
    let t = out(&["flagstat", &fixture()]);
    assert!(t.contains("306 + 10 mapped (95.63% : 83.33%)\n"), "{t}");
    assert!(t.contains("7 + 0 singletons (4.37% : N/A)\n"), "{t}");
    assert!(t.contains("40 + 1 secondary\n"), "{t}");
    assert!(
        t.contains("22 + 0 with mate mapped to a different chr (mapQ>=5)\n"),
        "{t}"
    );
}

#[test]
fn flagstat_format_inference() {
    let d = tempdir();
    let bam = fixture();
    let json = d.path().join("s.json");
    gloci()
        .args(["flagstat", &bam, "-o", &s(&json)])
        .assert()
        .success()
        .stdout("");
    let j = std::fs::read_to_string(&json).unwrap();
    assert!(j.starts_with("{\n \"QC-passed reads\": { \n"), "{j}");
    // -O wins over the extension.
    let tsv = d.path().join("t.json");
    gloci()
        .args(["flagstat", &bam, "-O", "tsv", "-o", &s(&tsv)])
        .assert()
        .success();
    assert!(
        std::fs::read_to_string(&tsv)
            .unwrap()
            .starts_with("320\t12\ttotal (QC-passed reads + QC-failed reads)\n")
    );
    // .tsv infers TSV; anything else is text.
    let t = d.path().join("u.tsv");
    gloci()
        .args(["flagstat", &bam, "-o", &s(&t)])
        .assert()
        .success();
    assert!(
        std::fs::read_to_string(&t)
            .unwrap()
            .contains("\tmapped %\n")
    );
    let txt = d.path().join("v.txt");
    gloci()
        .args(["flagstat", &bam, "-o", &s(&txt)])
        .assert()
        .success();
    assert!(
        std::fs::read_to_string(&txt)
            .unwrap()
            .starts_with("320 + 12 in total")
    );
    // An unknown -O is a usage error that lists the names.
    gloci()
        .args(["flagstat", &bam, "-O", "csv"])
        .assert()
        .code(2)
        .stdout("")
        .stderr(contains("text, tsv, json"));
}

#[test]
fn flagstat_cram_equals_bam() {
    let cram = s(&inputs().join("aln/flags.cram"));
    let reference = s(&inputs().join("fa/ref.fa"));
    let from_cram = out(&["flagstat", &cram, "--reference", &reference]);
    let from_sam = out(&["flagstat", &s(&inputs().join("aln/flags.sam"))]);
    assert_eq!(from_cram, from_sam);
    assert!(from_cram.starts_with("13 + 1 in total"), "{from_cram}");
    gloci()
        .args(["flagstat", &cram])
        .assert()
        .code(1)
        .stdout("")
        .stderr(contains("--reference"));
}

/// The mates fixtures, attached and detached mates in every layout samtools writes (CRAM 3.0
/// and 3.1, slices of several references, name order, no reference, raw blocks), give one
/// report; flagstat decodes only the fields it counts from each (`flagstat-cram-mates-*`
/// compare it with samtools).
#[test]
fn flagstat_cram_mates_layouts_agree() {
    let reference = s(&inputs().join("fa/ref.fa"));
    let want = out(&[
        "flagstat",
        &s(&inputs().join("aln/mates.cram")),
        "--reference",
        &reference,
    ]);
    assert!(want.starts_with("182 + 16 in total"), "{want}");
    assert!(
        want.contains("18 + 2 with mate mapped to a different chr\n"),
        "{want}"
    );
    for v in ["v31", "multi", "namesorted", "noref", "raw"] {
        let f = s(&inputs().join(format!("aln/mates.{v}.cram")));
        assert_eq!(
            out(&["flagstat", &f, "--reference", &reference]),
            want,
            "{v}"
        );
    }
}

/// From a CRAM, flagstat decodes no bases (#cram-limits): a FASTA with the right names and
/// lengths but other bases gives the same counts, where `view`, which decodes them, refuses it.
/// overhang.cram's read 1.5 Mb past chr1, past gloci's limit for decoding bases, is counted.
#[test]
fn flagstat_cram_reads_no_bases() {
    let d = tempdir();
    let fa = std::fs::read_to_string(inputs().join("fa/ref.fa")).unwrap();
    let other: String = fa
        .lines()
        .map(|l| {
            if l.starts_with('>') {
                format!("{l}\n")
            } else {
                format!("{}\n", "A".repeat(l.len()))
            }
        })
        .collect();
    let wrong = d.path().join("other.fa");
    std::fs::write(&wrong, other).unwrap();
    std::fs::copy(
        inputs().join("fa/ref.fa.fai"),
        d.path().join("other.fa.fai"),
    )
    .unwrap();
    let reference = s(&inputs().join("fa/ref.fa"));
    let mates = s(&inputs().join("aln/mates.cram"));
    assert_eq!(
        out(&["flagstat", &mates, "--reference", &s(&wrong)]),
        out(&["flagstat", &mates, "--reference", &reference])
    );
    gloci()
        .args(["view", &mates, "--reference", &s(&wrong)])
        .assert()
        .code(1)
        .stderr(contains("does not match"));
    let overhang = s(&inputs().join("aln/overhang.cram"));
    let t = out(&["flagstat", &overhang, "--reference", &reference]);
    assert!(t.starts_with("2 + 0 in total"), "{t}");
}

#[test]
fn flagstat_sam_equals_bam() {
    let sam = out(&["flagstat", &s(&inputs().join("samtools/flagstat.sam"))]);
    assert_eq!(sam, out(&["flagstat", &fixture()]));
}

#[test]
fn flagstat_is_identical_at_any_thread_count() {
    let big = s(&inputs().join("aln/cat-big.bam"));
    let cram = s(&inputs().join("aln/many.cram"));
    let mates = s(&inputs().join("aln/mates.cram"));
    let mates_multi = s(&inputs().join("aln/mates.multi.cram"));
    let reference = s(&inputs().join("fa/ref.fa"));
    let one_bam = fixture();
    for args in [
        vec!["flagstat", big.as_str()],
        vec!["flagstat", cram.as_str(), "--reference", reference.as_str()],
        vec![
            "flagstat",
            mates.as_str(),
            "--reference",
            reference.as_str(),
        ],
        vec![
            "flagstat",
            mates_multi.as_str(),
            "--reference",
            reference.as_str(),
        ],
        vec!["flagstat", one_bam.as_str()],
    ] {
        let mut one = args.clone();
        one.extend(["-t", "1"]);
        let want = out(&one);
        for t in ["2", "8"] {
            let mut a = args.clone();
            a.extend(["-t", t]);
            assert_eq!(out(&a), want, "{args:?} -t {t}");
        }
    }
}

#[test]
fn flagstat_stdin_and_missing_file() {
    let sam = std::fs::read(inputs().join("samtools/flagstat.sam")).unwrap();
    gloci()
        .args(["flagstat", "-"])
        .write_stdin(sam)
        .assert()
        .success()
        .stdout(contains("320 + 12 in total"));
    gloci()
        .args(["flagstat", "/nonexistent/x.bam"])
        .assert()
        .code(1)
        .stdout("")
        .stderr(contains("cannot open /nonexistent/x.bam"));
}

#[test]
fn flagstat_other_format_extension_has_its_own_message() {
    let d = tempfile::tempdir().unwrap();
    gloci()
        .args(["flagstat", &fixture(), "-o", &s(&d.path().join("x.bed"))])
        .assert()
        .failure()
        .code(2)
        .stderr(contains("flagstat writes text, tsv or json").and(contains("x.json")));
}

/// `-o x.json.gz` is JSON and `-o x.tsv.gz` TSV: the compression suffix is stripped before the
/// layout is read from the name, and a `.gz` name is written BGZF-compressed, as every TSV output
/// is (M17). depth, samples and bedcov compress a `.gz` name the same way.
#[test]
fn flagstat_gz_output_infers_the_layout_and_compresses() {
    use std::io::Read;
    let d = tempdir();
    let unzip = |p: &Path| {
        let bytes = std::fs::read(p).unwrap();
        assert!(
            bytes.starts_with(&[0x1f, 0x8b]),
            "{} is not gzip",
            p.display()
        );
        let mut text = String::new();
        noodles_bgzf::io::Reader::new(&bytes[..])
            .read_to_string(&mut text)
            .unwrap();
        text
    };
    let bam = fixture();
    for (name, layout) in [
        ("f.json.gz", "json"),
        ("f.tsv.gz", "tsv"),
        ("f.txt.gz", "text"),
        ("F.JSON.BGZ", "json"),
    ] {
        let p = d.path().join(name);
        gloci()
            .args(["flagstat", &bam, "-o", &s(&p)])
            .assert()
            .success();
        assert_eq!(unzip(&p), out(&["flagstat", &bam, "-O", layout]), "{name}");
    }
    // A compressed name of another format is still refused.
    gloci()
        .args(["flagstat", &bam, "-o", &s(&d.path().join("f.vcf.gz"))])
        .assert()
        .code(2)
        .stderr(contains("flagstat writes text, tsv or json"));
    let depth_bam = s(&inputs().join("samtools/depth.bam"));
    let bed = s(&inputs().join("samtools/targets.bed"));
    for (args, name) in [
        (vec!["depth", depth_bam.as_str()], "d.tsv.gz"),
        (vec!["samples", depth_bam.as_str()], "s.tsv.gz"),
        (vec!["bedcov", bed.as_str(), depth_bam.as_str()], "b.tsv.gz"),
    ] {
        let p = d.path().join(name);
        gloci().args(&args).arg("-o").arg(&p).assert().success();
        assert_eq!(unzip(&p), out(&args), "{args:?}");
    }
}
