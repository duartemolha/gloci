mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;

const IDS: [&str; 19] = [
    "hg19", "hg38", "hs1", "GRCh37", "GRCh38", "mm10", "mm39", "GRCm38", "GRCm39", "rn7",
    "danRer11", "GRCz11", "dm6", "BDGP6", "ce11", "WBcel235", "sacCer3", "R64", "TAIR10",
];

fn w(name: &str, text: &str) -> String {
    // Each call gets its own directory so equal names never collide.
    let d = tempfile::tempdir().unwrap().keep();
    plain_file(&d, name, text).to_string_lossy().into_owned()
}

fn stdout_of<const N: usize>(args: [&str; N]) -> String {
    let out = gloci()
        .args(args)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    String::from_utf8(out).unwrap()
}

/// hg38 has chrEBV (171,823 bp), as GRCh38 analysis-set files such as ENCODE's do.
#[test]
fn hg38_has_chrebv() {
    gloci()
        .args(["check", &w("ebv.bed", "chrEBV\t100\t200\n"), "-g", "hg38"])
        .assert()
        .success();
    gloci()
        .args([
            "check",
            &w("past.bed", "chrEBV\t171800\t171900\n"),
            "-g",
            "hg38",
        ])
        .assert()
        .code(1)
        .stdout(contains("past the end of chrEBV (length 171823)"))
        .stdout(contains("built-in genome hg38"));
    let d = tempdir();
    gloci()
        .args([
            "convert",
            &w("past.bed", "chrEBV\t171800\t171900\n"),
            "-g",
            "hg38",
        ])
        .arg("-o")
        .arg(d.path().join("x.bb"))
        .assert()
        .code(1)
        .stderr(contains("chromosome chrEBV (171823 bases, from -g hg38)"));
}

#[test]
fn sort_with_builtin_genome() {
    gloci()
        .args([
            "sort",
            &w("u.bed", "chr2\t5\t6\nchr1\t1\t2\n"),
            "-g",
            "hg38",
        ])
        .assert()
        .success()
        .stdout("chr1\t1\t2\nchr2\t5\t6\n");
    gloci()
        .args(["sort", &w("v.bed", "chrZ\t1\t2\n"), "-g", "hg38"])
        .assert()
        .code(1)
        .stderr(contains("is not in the built-in genome hg38"))
        .stderr(contains("gloci genomes show hg38"));
    gloci()
        .args(["sort", &w("x.bed", "1\t1\t2\n"), "-g", "hg37"])
        .assert()
        .code(2)
        .stderr(contains("did you mean"));
}

#[test]
fn genomes_command() {
    let out = stdout_of(["genomes"]);
    for id in IDS {
        assert!(out.contains(id), "{id}");
    }
    let hg38 = out.lines().find(|l| l.starts_with("hg38")).unwrap();
    let cols: Vec<_> = hg38.split_whitespace().collect();
    assert_eq!(cols[0], "hg38");
    assert!(hg38.contains("ucsc") && hg38.contains("456") && hg38.contains("3209457928"));
    assert!(
        stdout_of(["genomes", "show", "hg38"])
            .starts_with("chr1\t248956422\nchr2\t242193529\nchr3\t198295559\n")
    );
    assert!(stdout_of(["genomes", "show", "GRCH38"]).starts_with("1\t248956422\n2\t"));
    gloci()
        .args(["genomes", "show", "hg37"])
        .assert()
        .code(2)
        .stderr(contains("did you mean"));
}

/// `gloci genomes list` is the listing, shown under Commands in `gloci genomes --help`;
/// `gloci genomes` alone still prints it.
#[test]
fn genomes_list_is_the_listing() {
    let listed = stdout_of(["genomes", "list"]);
    assert_eq!(listed, stdout_of(["genomes"]));
    assert_eq!(listed.lines().count(), IDS.len());
    for (line, id) in listed.lines().zip(IDS) {
        assert_eq!(line.split('\t').next(), Some(id));
        assert_eq!(line.split('\t').count(), 6, "{line}");
    }
    let help = stdout_of(["genomes", "--help"]);
    assert!(
        help.lines().any(|l| l.trim_start().starts_with("list ")),
        "{help}"
    );
    gloci().args(["genomes", "list", "hg38"]).assert().code(2);
}

/// An unknown name says where the list is. `genomes show` takes only an id, so its message
/// offers no file; `-g` also takes files and says so.
#[test]
fn an_unknown_genome_points_at_the_list() {
    gloci()
        .args(["genomes", "show", "human"])
        .assert()
        .code(2)
        .stderr(contains("unknown genome \"human\""))
        .stderr(contains("built-in genomes: hg19, hg38,"))
        .stderr(contains(
            "gloci genomes list shows their species and assemblies",
        ))
        .stderr(contains(".fai").not());
    gloci()
        .args(["sort", &w("h.bed", "chr1\t1\t2\n"), "-g", "human"])
        .assert()
        .code(2)
        .stderr(contains("built-in genomes: hg19, hg38,"))
        .stderr(contains(
            "gloci genomes list shows their species and assemblies",
        ))
        .stderr(contains(
            "-g takes a built-in genome or an existing .fai, .genome, .dict, VCF/BCF or BAM file",
        ));
    gloci()
        .args(["genomes", "show", "hg37"])
        .assert()
        .code(2)
        .stderr(contains("unknown genome \"hg37\" (did you mean"));
}

#[test]
fn check_bounds_with_builtin() {
    gloci()
        .args(["check", &w("b.bed", "chrM\t0\t20000\n"), "-g", "hg38"])
        .assert()
        .code(1)
        .stdout(contains("past the end of chrM"));
}

#[test]
fn view_with_builtin_genome() {
    let f = w("r.bed", "chr1\t1\t100\nchr2\t5\t60\n");
    gloci()
        .args(["view", &f, "-g", "hg38", "chr2:1-100"])
        .assert()
        .success()
        .stdout("chr2\t5\t60\n");
}

#[test]
fn out_of_order_hint_names_the_builtin() {
    gloci()
        .args([
            "check",
            &w("o.bed", "chr2\t5\t6\nchr1\t1\t2\n"),
            "-g",
            "hg38",
        ])
        .assert()
        .code(1)
        .stdout(contains("-g hg38"));
}

#[test]
fn a_file_named_like_a_builtin_wins_with_a_warning() {
    let d = tempfile::tempdir().unwrap();
    plain_file(d.path(), "hg38", "chrA\t10\n");
    let bed = plain_file(d.path(), "a.bed", "chrA\t1\t2\n");
    gloci()
        .current_dir(d.path())
        .args(["sort", bed.to_str().unwrap(), "-g", "hg38"])
        .assert()
        .success()
        .stdout("chrA\t1\t2\n")
        .stderr(contains("rename or move that file"));
}

#[test]
fn a_directory_named_like_a_builtin_does_not_hijack_it() {
    let d = tempfile::tempdir().unwrap();
    std::fs::create_dir(d.path().join("hg38")).unwrap();
    let bed = plain_file(d.path(), "a.bed", "chr2\t5\t6\nchr1\t1\t2\n");
    gloci()
        .current_dir(d.path())
        .args(["sort", bed.to_str().unwrap(), "-g", "hg38"])
        .assert()
        .success()
        .stdout("chr1\t1\t2\nchr2\t5\t6\n");
}
