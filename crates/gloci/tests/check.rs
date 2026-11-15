mod common;
use common::*;
use predicates::str::contains;
use std::path::Path;

fn p(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

fn stdout_of(args: &[&str]) -> String {
    let o = gloci().args(args).assert().get_output().stdout.clone();
    String::from_utf8(o).unwrap()
}

const CLEAN: &str = "chr1\t100\t200\tgeneA\nchr1\t300\t400\tgeneB\nchr2\t10\t50\tgeneC\n";

fn clean_bed(d: &Path) -> String {
    p(&plain_file(d, "clean.bed", CLEAN))
}

/// Exactly 14 errors: an order violation (line 2), start > end (line 3) and
/// 12 non-integer starts (lines 4-15).
fn messy_bed(d: &Path) -> String {
    let mut t = String::from("chr1\t100\t200\nchr1\t50\t60\nchr1\t10\t5\n");
    for _ in 0..12 {
        t.push_str("chr1\tabc\t20\n");
    }
    p(&plain_file(d, "messy.bed", &t))
}

#[test]
fn clean_file() {
    let d = tempdir();
    let f = clean_bed(d.path());
    insta::assert_snapshot!(stdout_of(&["check", &f]).replace(&f, "clean.bed"));
    gloci().args(["check", &f]).assert().success();
}

#[test]
fn problems_listed_and_counted() {
    let d = tempdir();
    let f = messy_bed(d.path());
    let o = gloci()
        .args(["check", &f])
        .assert()
        .code(1)
        .get_output()
        .stdout
        .clone();
    let s = String::from_utf8(o).unwrap();
    assert!(
        s.contains("errors:       14")
            && s.contains("line 3: start 10 is greater than end 5")
            && s.contains("order:        NOT sorted"),
        "{s}"
    );
    assert_eq!(s.lines().filter(|l| l.starts_with("    line ")).count(), 10);
    // The first order violation carries the sort hint.
    assert!(s.contains("line 2: chr1 start 50 is before the previous start 100"));
    assert!(s.contains("run: gloci sort"));
}

#[test]
fn messy_snapshot() {
    let d = tempdir();
    let f = messy_bed(d.path());
    insta::assert_snapshot!(stdout_of(&["check", &f]).replace(&f, "messy.bed"));
}

#[test]
fn genome_bounds() {
    let d = tempdir();
    let g = p(&plain_file(d.path(), "g.genome", "chr1\t200\nchr2\t500\n"));
    let f = p(&plain_file(
        d.path(),
        "past.bed",
        "chr1\t100\t200\nchr1\t150\t300\nchr3\t1\t2\n",
    ));
    gloci()
        .args(["check", &f, "-g", &g])
        .assert()
        .code(1)
        .stdout(contains(
            "line 2: end 300 is past the end of chr1 (length 200)",
        ))
        .stdout(contains(
            "line 3: chromosome \"chr3\" is not in the genome file",
        ))
        .stdout(contains("errors:       2"));
}

#[test]
fn duplicates_warn_only() {
    let d = tempdir();
    let f = p(&plain_file(
        d.path(),
        "dup.bed",
        "chr1\t100\t200\tx\nchr1\t100\t300\ty\nchr1\t100\t200\tx\nchr1\t500\t600\tx\nchr1\t700\t800\tx\n",
    ));
    gloci()
        .args(["check", &f])
        .assert()
        .success()
        .stdout(contains("warnings:     1"))
        .stdout(contains("line 3: duplicate of an earlier line"))
        .stdout(contains("errors:       0"))
        .stdout(contains("result: OK"));
}

#[test]
fn missing_file_does_not_discard_other_reports() {
    let d = tempdir();
    let good = clean_bed(d.path());
    let good2 = p(&plain_file(d.path(), "good2.bed", "chrX\t1\t2\n"));
    let missing = p(&d.path().join("missing.bed"));
    let s = stdout_of(&["check", &good, &missing, &good2]);
    assert!(s.contains("records:      3"), "{s}");
    assert!(s.contains("chromosomes:  1 (chrX)"), "{s}");
    assert!(s.contains("missing.bed: cannot be checked"), "{s}");
    assert!(s.contains("cannot open"), "{s}");
    assert!(
        s.contains("result: FAILED (1 of 3 files have errors)"),
        "{s}"
    );
    gloci()
        .args(["check", &good, &missing, &good2])
        .assert()
        .code(1);
}

#[test]
fn hint_visible_after_many_parse_errors() {
    let d = tempdir();
    let mut t = String::from("chr1\t10\t20\n");
    for _ in 0..12 {
        t.push_str("chr1\tabc\t20\n");
    }
    t.push_str("chr1\t100\t200\nchr1\t50\t60\nchr2\t1\t2\nchr10\t1\t2\n");
    let f = p(&plain_file(d.path(), "h.bed", &t));
    let s = stdout_of(&["check", &f]);
    assert!(
        !s.contains("line 15:"),
        "order problems are past the listed 10: {s}"
    );
    assert!(
        s.contains("order:        NOT sorted\n                hint: run: gloci sort"),
        "{s}"
    );
    // The karyotype note comes from a later chromosome violation.
    assert!(s.contains("--chrom-order karyotype"), "{s}");
}

#[test]
fn blank_lines_start_middle_end_plain_and_stdin() {
    let d = tempdir();
    let text = "\nchr1\t1\t2\n\nchr1\t3\t4\n\n\n";
    let f = p(&plain_file(d.path(), "b.bed", text));
    let check = |s: String| {
        assert!(s.contains("warnings:     4"), "{s}");
        for n in [1, 3, 5, 6] {
            assert!(s.contains(&format!("line {n}: blank line")), "{s}");
        }
    };
    check(stdout_of(&["check", &f]));
    let o = gloci()
        .args(["check", "-", "--format", "bed"])
        .write_stdin(text)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    check(String::from_utf8(o).unwrap());
}

#[test]
fn comments_and_invalid_lines_are_not_blank_warnings() {
    let d = tempdir();
    let f = p(&plain_file(
        d.path(),
        "c.bed",
        "chr1\t1\t2\n# note\nchr1\tx\t4\nchr1\t5\t6\n",
    ));
    gloci()
        .args(["check", &f])
        .assert()
        .code(1)
        .stdout(contains("warnings:     0"));
}

#[test]
fn blank_lines_warn() {
    let d = tempdir();
    let f = p(&plain_file(
        d.path(),
        "blank.bed",
        "chr1\t1\t2\n\nchr1\t3\t4\n\nchr1\t5\t6\n",
    ));
    gloci()
        .args(["check", &f])
        .assert()
        .success()
        .stdout(contains("warnings:     2"))
        .stdout(contains("line 2: blank line"))
        .stdout(contains("line 4: blank line"));
}

#[test]
fn vcf_column_count() {
    let d = tempdir();
    let f = p(&plain_file(
        d.path(),
        "c.vcf",
        "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n\
         chr1\t5\t.\tA\tC\t.\t.\t.\tGT\t0/1\t1/1\n\
         chr1\t9\t.\tA\tC\t.\t.\t.\tGT\t0/1\n",
    ));
    gloci()
        .args(["check", &f])
        .assert()
        .code(1)
        .stdout(contains("line 5: expected 11 columns, found 10"))
        .stdout(contains("order:        sorted (##contig order of"));
}

#[test]
fn vcf_symbolic_alt_without_end_is_error() {
    let d = tempdir();
    let f = p(&plain_file(
        d.path(),
        "s.vcf",
        "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         chr1\t5\t.\tA\t<DEL>\t.\t.\t.\n",
    ));
    gloci()
        .args(["check", &f])
        .assert()
        .code(1)
        .stdout(contains("errors:       1"))
        .stdout(contains("line 4:"));
}

#[test]
fn bed12_blocks_checked() {
    let d = tempdir();
    let f = p(&plain_file(
        d.path(),
        "b12.bed",
        "chr1\t100\t200\tg\t0\t+\t100\t200\t0\t2\t10,10\t0,50\n\
         chr1\t300\t400\tg\t0\t+\t300\t400\t0\t2\t10,10\t0,50\n",
    ));
    gloci()
        .args(["check", &f])
        .assert()
        .code(1)
        .stdout(contains("errors:       2"))
        .stdout(contains("line 1:"));
}

#[test]
fn bedgraph_value_validated() {
    let d = tempdir();
    let f = p(&plain_file(
        d.path(),
        "v.bedgraph",
        "chr1\t0\t10\t1.5\nchr1\t10\t20\tabc\nchr1\t20\t30\n\
         chr1\t30\t40\t-2e3\nchr1\t40\t50\tnan\n",
    ));
    let s = stdout_of(&["check", &f, "--format", "bedgraph"]);
    assert!(
        s.contains("BEDGRAPH, plain") && s.contains("errors:       3"),
        "{s}"
    );
    assert!(s.contains("line 2: value \"abc\" is not a number"), "{s}");
    assert!(s.contains("line 3: value \"\" is not a number"), "{s}");
    assert!(s.contains("line 5: value \"nan\" is not a number"), "{s}");
    gloci()
        .args(["check", &f, "--format", "bedgraph"])
        .assert()
        .code(1);
}

#[test]
fn karyotype_sorted_file_hint_under_default_order() {
    let d = tempdir();
    let kar = p(&plain_file(
        d.path(),
        "kar.bed",
        "chr1\t1\t2\nchr2\t1\t2\nchr10\t1\t2\n",
    ));
    gloci()
        .args(["check", &kar])
        .assert()
        .code(1)
        .stdout(contains("order:        NOT sorted"))
        .stdout(contains("looks sorted in karyotype order"))
        .stdout(contains("--chrom-order karyotype"));
    gloci()
        .args(["check", &kar, "--chrom-order", "karyotype"])
        .assert()
        .success()
        .stdout(contains("sorted (karyotype order)"));
    let lex = p(&plain_file(
        d.path(),
        "lex.bed",
        "chr1\t1\t2\nchr10\t1\t2\nchr2\t1\t2\n",
    ));
    gloci()
        .args(["check", &lex])
        .assert()
        .success()
        .stdout(contains("sorted (lexicographic order)"));
    // The other direction keeps its hint.
    gloci()
        .args(["check", &lex, "--chrom-order", "karyotype"])
        .assert()
        .code(1)
        .stdout(contains(
            "looks sorted by sort -k1,1, sort-bed or bedtools sort",
        ))
        .stdout(contains("gloci sort --chrom-order karyotype"));
}

#[test]
fn many_chromosomes_summarised() {
    let d = tempdir();
    let f = p(&plain_file(
        d.path(),
        "many.bed",
        "chr1\t1\t2\nchr2\t1\t2\nchr3\t1\t2\nchr4\t1\t2\nchr5\t1\t2\n",
    ));
    gloci()
        .args(["check", &f])
        .assert()
        .success()
        .stdout(contains("chromosomes:  5 (chr1, chr2, chr3, … 2 more)"));
}

#[test]
fn multiple_files() {
    let d = tempdir();
    let (a, b) = (clean_bed(d.path()), messy_bed(d.path()));
    gloci()
        .args(["check", &a, &b])
        .assert()
        .code(1)
        .stdout(contains("result: FAILED (1 of 2 files have errors)"));
}

#[test]
fn empty_file() {
    let d = tempdir();
    let f = p(&plain_file(d.path(), "empty.bed", ""));
    gloci()
        .args(["check", &f])
        .assert()
        .success()
        .stdout(contains("records:      0"));
}

#[test]
fn bgzf_input() {
    let d = tempdir();
    let f = p(&bgzf_file(d.path(), "c.bed.gz", CLEAN));
    gloci()
        .args(["check", &f])
        .assert()
        .success()
        .stdout(contains("BED, BGZF"))
        .stdout(contains("records:      3"));
}
