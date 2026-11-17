//! `gloci liftover` (D174): the two outputs, warnings, usage errors and threads. The mapped and
//! unmapped records themselves are checked against CrossMap by the `liftover-*-vs-crossmap`
//! golden cases, whose fixtures are used here.
mod common;
use common::*;
use predicates::prelude::PredicateBooleanExt;
use predicates::str::contains;
use std::path::{Path, PathBuf};

const TOY_BED_OUT: &str = "chr1\t1120\t1180\tr1\t0\t+
chr1\t1150\t1200\tr2\t0\t+
chr1\t1250\t1300\tr2\t0\t+
chr1\t1330\t1350\tr3\t0\t+
chr1\t1350\t1380\tr3\t0\t+
chr1\t1100\t1150\tr5\t0\t+
chr2\t3940\t3980\tr6\t0\t-
chr2\t3920\t3940\tr7\t0\t+
chr2\t3880\t3900\tr7\t0\t+
chr1\t1420\t1440\tr8\t0\t+
chr3\t20\t40\tr8\t0\t+
chr4\t30\t50\tr8\t0\t+
chr1\t1150\t1150\tr9\t0\t+
chr5\t510\t520\tr11\t0\t+
5\t530\t540\tr12\t0\t-
1\t1120\t1180\tr13
chr1\t1120\t1130
chr2\t3940\t3980\tr16\t0\t.
chr1\t1550\t1600\tr17\t0\t+
chr2\t3950\t4000\tr17\t0\t-
chr1\t1100\t1101\tr19\t0\t+
chr1\t1100\t1200\tr21\t0\t+
chr1\t1250\t1350\tr21\t0\t+
chr1\t1350\t1600\tr21\t0\t+
chr3\t0\t100\tr21\t0\t+
chr4\t10\t60\tr21\t0\t+
chr2\t3920\t4000\tr21\t0\t-
chr2\t3800\t3900\tr21\t0\t-
chrU\t50\t100\tr22\t0\t+
chrU\t850\t900\tr22\t0\t-
chr1\t1160\t1170\tr23\t0\t+\tx\t1
chr1\t1420\t1420\tr24\t0\t+
chr3\t20\t20\tr24\t0\t+
chr4\t30\t30\tr24\t0\t+
chr1\t1180\t1200\tr25\t0\t+
";

const TOY_BED_UNMAPPED: &str = "chr1\t310\t340\tr4\t0\t+
chr1\t200\t200\tr10\t0\t+
chrZ\t1\t2\tr15\t0\t+
chr1\t100\t100\tr18\t0\t+
chr1\t600\t600\tr20\t0\t+
chr1\t600\t610\tr26\t0\t+
chrM\t10\t20\tr27\t0\t+
";

fn fixture(name: &str) -> String {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/liftover")
        .join(name)
        .display()
        .to_string()
}

/// A temp directory holding the fixtures a test names.
fn with(names: &[&str]) -> tempfile::TempDir {
    let d = tempdir();
    for n in names {
        std::fs::copy(fixture(n), d.path().join(n)).unwrap();
    }
    d
}

fn read(p: &Path) -> String {
    std::fs::read_to_string(p).unwrap()
}

#[test]
fn liftover_writes_mapped_and_unmapped() {
    let d = with(&["toy.bed", "toy.chain"]);
    gloci()
        .current_dir(d.path())
        .args([
            "liftover",
            "toy.bed",
            "--chain",
            "toy.chain",
            "-o",
            "m.bed",
            "--unmapped",
            "u.bed",
        ])
        .assert()
        .success()
        .stderr(contains("8 records were split").and(contains("6 lost bases")))
        .stderr(contains("did not lift").not());
    assert_eq!(read(&d.path().join("m.bed")), TOY_BED_OUT);
    assert_eq!(read(&d.path().join("u.bed")), TOY_BED_UNMAPPED);
}

#[test]
fn unmapped_records_are_counted_without_unmapped() {
    let d = with(&["toy.bed", "toy.chain"]);
    gloci()
        .current_dir(d.path())
        .args(["liftover", "toy.bed", "--chain", "toy.chain"])
        .assert()
        .success()
        .stdout(TOY_BED_OUT)
        .stderr(
            contains("7 of 27 records did not lift")
                .and(contains("on no chain"))
                .and(contains("--unmapped")),
        );
}

#[test]
fn both_outputs_on_stdout_is_a_usage_error() {
    let d = with(&["toy.bed", "toy.chain"]);
    gloci()
        .current_dir(d.path())
        .args([
            "liftover",
            "toy.bed",
            "--chain",
            "toy.chain",
            "--unmapped",
            "-",
        ])
        .assert()
        .code(2)
        .stderr(contains("--unmapped").and(contains("-o")));
}

/// `--unmapped` is an output: with the input on stdin, `--unmapped -` is stdout (`-o` names a
/// file), not a second use of standard input.
#[test]
fn stdin_input_with_unmapped_on_stdout() {
    let d = with(&["toy.bed", "toy.chain"]);
    gloci()
        .current_dir(d.path())
        .args([
            "liftover",
            "-",
            "--chain",
            "toy.chain",
            "-o",
            "m.bed",
            "--unmapped",
            "-",
        ])
        .write_stdin(read(&d.path().join("toy.bed")))
        .assert()
        .success()
        .stdout(TOY_BED_UNMAPPED);
    assert_eq!(read(&d.path().join("m.bed")), TOY_BED_OUT);
}

/// A chain on standard input, gzip-compressed, lifts as the file does.
#[test]
fn chain_on_stdin() {
    let d = with(&["toy.bed", "toy.chain"]);
    let gz = gloci_testkit::fixtures::gzip(&std::fs::read(d.path().join("toy.chain")).unwrap());
    gloci()
        .current_dir(d.path())
        .args(["liftover", "toy.bed", "--chain", "-"])
        .write_stdin(gz)
        .assert()
        .success()
        .stdout(TOY_BED_OUT);
}

#[test]
fn same_path_for_both_outputs_is_a_usage_error() {
    let d = with(&["toy.bed", "toy.chain"]);
    gloci()
        .current_dir(d.path())
        .args([
            "liftover",
            "toy.bed",
            "--chain",
            "toy.chain",
            "-o",
            "x.bed",
            "--unmapped",
            "x.bed",
        ])
        .assert()
        .code(2);
    assert!(!d.path().join("x.bed").exists());
}

#[test]
fn vcf_input_is_refused() {
    let d = with(&["toy.chain"]);
    plain_file(
        d.path(),
        "a.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t130\t.\tA\tC\t.\t.\t.\n",
    );
    gloci()
        .current_dir(d.path())
        .args(["liftover", "a.vcf", "--chain", "toy.chain"])
        .assert()
        .code(2)
        .stderr(contains("VCF liftover"));
}

#[test]
fn bad_chain_exits_1_naming_the_line() {
    let d = with(&["toy.bed", "bad-end.chain"]);
    gloci()
        .current_dir(d.path())
        .args(["liftover", "toy.bed", "--chain", "bad-end.chain"])
        .assert()
        .code(1)
        .stderr(contains("bad-end.chain line 4").and(contains("chain at line 1")));
}

#[test]
fn failed_run_leaves_no_outputs() {
    let d = with(&["toy-invalid.bed", "toy.chain"]);
    gloci()
        .current_dir(d.path())
        .args([
            "liftover",
            "toy-invalid.bed",
            "--chain",
            "toy.chain",
            "-o",
            "m.bed",
            "--unmapped",
            "u.bed",
        ])
        .assert()
        .code(1)
        .stderr(contains("line 2"));
    assert!(!d.path().join("m.bed").exists());
    assert!(!d.path().join("u.bed").exists());
}

#[test]
fn skip_invalid_lifts_the_rest() {
    let d = with(&["toy-invalid.bed", "toy.chain"]);
    gloci()
        .current_dir(d.path())
        .args([
            "liftover",
            "toy-invalid.bed",
            "--chain",
            "toy.chain",
            "--skip-invalid",
        ])
        .assert()
        .success()
        .stdout("chr1\t1120\t1180\tv1\nchr1\t1140\t1150\tv3\n");
}

#[test]
fn output_is_the_same_at_every_thread_count() {
    let d = with(&["toy.chain"]);
    let toy = read(Path::new(&fixture("toy.bed")));
    let mut big = String::new();
    for i in 0..2000 {
        for line in toy.lines() {
            let mut f: Vec<String> = line.split('\t').map(str::to_string).collect();
            if let Some(name) = f.get_mut(3) {
                *name = format!("{name}_{i}");
            }
            big.push_str(&f.join("\t"));
            big.push('\n');
        }
    }
    bgzf_file(d.path(), "big.bed.gz", &big);
    let mut runs = Vec::new();
    for t in ["1", "2", "8"] {
        let (m, u) = (format!("m{t}.bed.gz"), format!("u{t}.bed.gz"));
        gloci()
            .current_dir(d.path())
            .args([
                "liftover",
                "big.bed.gz",
                "--chain",
                "toy.chain",
                "-o",
                &m,
                "--unmapped",
                &u,
                "-t",
                t,
            ])
            .assert()
            .success();
        runs.push((
            std::fs::read(d.path().join(&m)).unwrap(),
            std::fs::read(d.path().join(&u)).unwrap(),
        ));
    }
    assert!(runs[0].0.len() > 1000 && runs[0].1.len() > 100);
    assert!(runs.iter().all(|r| *r == runs[0]));
}

#[test]
fn header_lines_only_with_header() {
    let d = with(&["toy.chain"]);
    let toy = read(Path::new(&fixture("toy.bed")));
    plain_file(d.path(), "h.bed", &format!("#chrom\tstart\tend\n{toy}"));
    let run = |extra: &[&str]| {
        gloci()
            .current_dir(d.path())
            .args([
                "liftover",
                "h.bed",
                "--chain",
                "toy.chain",
                "-o",
                "m.bed",
                "--unmapped",
                "u.bed",
            ])
            .args(extra)
            .assert()
            .success();
        (read(&d.path().join("m.bed")), read(&d.path().join("u.bed")))
    };
    let (m, u) = run(&[]);
    assert!(!m.contains("#chrom") && !u.contains("#chrom"));
    let (m, u) = run(&["--header"]);
    assert!(m.starts_with("#chrom\tstart\tend\n") && u.starts_with("#chrom\tstart\tend\n"));
}

#[test]
fn help_names_crossmap_and_the_ucsc_licence() {
    gloci()
        .args(["liftover", "--help"])
        .assert()
        .success()
        .stdout(
            contains("CrossMap")
                .and(contains("liftOver"))
                .and(contains("non-commercial"))
                .and(contains("gloci sort")),
        );
}

#[test]
fn unmapped_takes_any_name() {
    let d = with(&["toy.bed", "toy.chain"]);
    for name in ["m.bed.unmap", "u.txt", "noext"] {
        gloci()
            .current_dir(d.path())
            .args([
                "liftover",
                "toy.bed",
                "--chain",
                "toy.chain",
                "-o",
                "m.bed",
                "--unmapped",
                name,
            ])
            .assert()
            .success();
        assert_eq!(read(&d.path().join(name)), TOY_BED_UNMAPPED, "{name}");
    }
    assert_eq!(read(&d.path().join("m.bed")), TOY_BED_OUT);
}

#[test]
fn unmapped_gz_is_bgzf() {
    let d = with(&["toy.bed", "toy.chain"]);
    gloci()
        .current_dir(d.path())
        .args([
            "liftover",
            "toy.bed",
            "--chain",
            "toy.chain",
            "-o",
            "m.bed",
            "--unmapped",
            "u.unmap.gz",
        ])
        .assert()
        .success();
    let raw = std::fs::read(d.path().join("u.unmap.gz")).unwrap();
    assert_eq!(raw[..4], [0x1f, 0x8b, 0x08, 0x04]);
    let text = gloci()
        .args(["view", "-"])
        .write_stdin(raw)
        .output()
        .unwrap();
    // Read back through gloci itself: a BGZF BED of the unmapped records.
    assert_eq!(String::from_utf8_lossy(&text.stdout), TOY_BED_UNMAPPED);
}

#[cfg(unix)]
#[test]
fn a_link_to_the_output_is_the_same_file() {
    let d = with(&["toy.bed", "toy.chain"]);
    std::fs::write(d.path().join("x.bed"), "").unwrap();
    std::os::unix::fs::symlink("x.bed", d.path().join("y.bed")).unwrap();
    std::fs::hard_link(d.path().join("x.bed"), d.path().join("z.bed")).unwrap();
    for other in ["y.bed", "z.bed"] {
        gloci()
            .current_dir(d.path())
            .args([
                "liftover",
                "toy.bed",
                "--chain",
                "toy.chain",
                "-o",
                "x.bed",
                "--unmapped",
                other,
            ])
            .assert()
            .code(2);
    }
}

#[test]
fn invalid_blocks_name_the_line() {
    let d = with(&["toy.chain"]);
    std::fs::write(
        d.path().join("bad.bed"),
        "chr1\t120\t580\tx\t0\t+\t130\t560\t0\t2\t60,60,60,\t0,100,400,\n",
    )
    .unwrap();
    gloci()
        .current_dir(d.path())
        .args(["liftover", "bad.bed", "--chain", "toy.chain"])
        .assert()
        .code(1)
        .stderr(contains("line 1").and(contains("--ignore-blocks")));
    gloci()
        .current_dir(d.path())
        .args([
            "liftover",
            "bad.bed",
            "--chain",
            "toy.chain",
            "--skip-invalid",
        ])
        .assert()
        .success()
        .stderr(contains("skipped 1 invalid record in bad.bed").and(contains("line 1")));
}

/// The help's example shape: a GTF in, a GTF on stdout.
#[test]
fn gtf_lifts_to_stdout_as_the_help_shows() {
    let d = with(&["toy.gtf", "toy.chain"]);
    gloci()
        .current_dir(d.path())
        .args(["liftover", "toy.gtf", "--chain", "toy.chain"])
        .assert()
        .success()
        .stdout(
            "chr1\tgloci\texon\t1121\t1180\t.\t+\t.\tgene_id \"a\"; transcript_id \"a.1\";\n\
             chr2\tgloci\tCDS\t3941\t3980\t.\t+\t0\tgene_id \"b\"; transcript_id \"b.1\";\n",
        )
        .stderr(contains("1 of 3 records did not lift"));
}

#[test]
fn ignore_blocks_on_gff_exits_2() {
    let d = with(&["toy.gff", "toy.chain"]);
    gloci()
        .current_dir(d.path())
        .args([
            "liftover",
            "toy.gff",
            "--chain",
            "toy.chain",
            "--ignore-blocks",
        ])
        .assert()
        .code(2)
        .stderr(contains("--ignore-blocks"));
}

#[test]
fn min_ratio_needs_merge_segments() {
    let d = with(&["toy.bed", "toy.chain"]);
    gloci()
        .current_dir(d.path())
        .args([
            "liftover",
            "toy.bed",
            "--chain",
            "toy.chain",
            "--min-ratio",
            "0.5",
        ])
        .assert()
        .code(2)
        .stderr(contains("--merge-segments"));
}

#[test]
fn min_ratio_range() {
    let d = with(&["toy.bed", "toy.chain"]);
    for r in ["1.5", "-0.1", "nan"] {
        gloci()
            .current_dir(d.path())
            .args([
                "liftover",
                "toy.bed",
                "--chain",
                "toy.chain",
                "--merge-segments",
                "--min-ratio",
                r,
            ])
            .assert()
            .code(2)
            .stderr(contains("--min-ratio"));
    }
}

// Built-in chains (D174, R13).

#[test]
fn builtin_chain_equals_the_file() {
    // The golden was made from a slice of Ensembl's chain; the built-in one is the full file.
    let golden = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden/liftover-ensembl-bed-vs-crossmap/expected");
    gloci()
        .args([
            "liftover",
            &fixture("grch37-21-22.bed"),
            "--chain",
            "grch37-to-grch38",
        ])
        .assert()
        .success()
        .stdout(read(&golden));
}

#[test]
fn unknown_builtin_suggests() {
    gloci()
        .args([
            "liftover",
            &fixture("grch37-21-22.bed"),
            "--chain",
            "GRCh37-to-GRCh39",
        ])
        .assert()
        .code(2)
        .stderr(contains("did you mean GRCh37-to-GRCh38").and(contains("CHM13v2-to-GRCh37")));
}

#[test]
fn hg19_chrm_is_refused() {
    let d = tempdir();
    let bed = plain_file(
        d.path(),
        "m.bed",
        "chrM\t100\t200\tx\nchr21\t9411300\t9411400\ty\n",
    );
    gloci()
        .args([
            "liftover",
            bed.to_str().unwrap(),
            "--chain",
            "GRCh37-to-GRCh38",
        ])
        .args([
            "--unmapped",
            "-",
            "-o",
            d.path().join("o.bed").to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout("chrM\t100\t200\tx\n")
        .stderr(contains("hg19's chrM is not GRCh37's MT"));
    assert_eq!(read(&d.path().join("o.bed")).lines().count(), 1);
}

#[test]
fn hg19_chrm_warning_comes_once() {
    let d = tempdir();
    let bed = plain_file(
        d.path(),
        "m.bed",
        "chrM\t100\t200\nchrM\t300\t400\nchrM\t500\t600\n",
    );
    let out = gloci()
        .args([
            "liftover",
            bed.to_str().unwrap(),
            "--chain",
            "GRCh37-to-GRCh38",
        ])
        .output()
        .unwrap();
    let err = String::from_utf8_lossy(&out.stderr);
    assert_eq!(err.matches("hg19's chrM").count(), 1, "{err}");
}

#[test]
fn grch38_chrm_keeps_its_name() {
    let d = tempdir();
    let bed = plain_file(d.path(), "m.bed", "chrM\t100\t200\n");
    gloci()
        .args([
            "liftover",
            bed.to_str().unwrap(),
            "--chain",
            "GRCh38-to-GRCh37",
        ])
        .assert()
        .success()
        .stdout(predicates::str::starts_with("chrM\t"))
        .stderr(contains("hg19's chrM").not());
}

#[test]
fn a_file_named_like_a_builtin_is_the_file() {
    let d = with(&["toy.bed", "toy.chain"]);
    std::fs::copy(
        d.path().join("toy.chain"),
        d.path().join("GRCh37-to-GRCh38"),
    )
    .unwrap();
    gloci()
        .current_dir(d.path())
        .args(["liftover", "toy.bed", "--chain", "GRCh37-to-GRCh38"])
        .assert()
        .success()
        .stdout(TOY_BED_OUT);
}
