//! `gloci view FILE REGION...`: region queries through an index, by
//! scanning (checked), and by binary search on plain BED with
//! `--assume-sorted` (§4.3, §7.4, D8-D11, D26).
mod common;
use common::*;
use gloci_core::Interval;
use gloci_testkit::fixtures;
use gloci_testkit::naive;
use gloci_testkit::strategies::TestRec;
use predicates::prelude::*;
use predicates::str::contains;
use std::cell::Cell;
use std::path::PathBuf;

/// Test files in one temp dir that lives as long as the test.
struct Fx {
    dir: tempfile::TempDir,
    n: Cell<u32>,
}

impl Fx {
    fn new() -> Fx {
        Fx {
            dir: tempdir(),
            n: Cell::new(0),
        }
    }

    /// A fresh path ending in `suffix`.
    fn path(&self, suffix: &str) -> PathBuf {
        self.n.set(self.n.get() + 1);
        self.dir.path().join(format!("f{}{suffix}", self.n.get()))
    }

    fn write(&self, suffix: &str, bytes: &[u8]) -> String {
        let p = self.path(suffix);
        std::fs::write(&p, bytes).unwrap();
        p.to_string_lossy().into_owned()
    }

    fn plain_bed(&self, text: &str) -> String {
        self.write(".bed", text.as_bytes())
    }

    fn bgzf_bed_unindexed(&self, text: &str) -> String {
        self.write(".bed.gz", &fixtures::bgzip(text.as_bytes()))
    }

    fn bgzf_bed_indexed(&self, text: &str) -> String {
        let p = self.bgzf_bed_unindexed(text);
        index(&p);
        p
    }

    /// The same BED text through each access path: plain BED scanned (the
    /// default) and binary-searched (`--assume-sorted`), BGZF scanned, and
    /// BGZF through its index.
    fn every_path(&self, text: &str) -> [Access; 4] {
        let plain = self.plain_bed(text);
        [
            Access::new(&plain, &[]),
            Access::new(&plain, &["--assume-sorted"]),
            Access::new(&self.bgzf_bed_unindexed(text), &[]),
            Access::new(&self.bgzf_bed_indexed(text), &[]),
        ]
    }

    fn vcf_gz(&self, text: &str) -> String {
        self.write(".vcf.gz", &fixtures::bgzip(text.as_bytes()))
    }

    fn vcf_gz_indexed(&self, text: &str) -> String {
        let p = self.vcf_gz(text);
        index(&p);
        p
    }
}

/// A file and the flags that choose how view reads it.
struct Access {
    file: String,
    flags: Vec<&'static str>,
}

impl Access {
    fn new(file: &str, flags: &[&'static str]) -> Access {
        Access {
            file: file.to_string(),
            flags: flags.to_vec(),
        }
    }

    /// stdout of `gloci view FILE FLAGS... ARGS...`, which must succeed.
    fn view(&self, args: &[&str]) -> String {
        let mut all = vec!["view", self.file.as_str()];
        all.extend(&self.flags);
        all.extend(args);
        stdout_of(all)
    }
}

fn index(path: &str) {
    gloci().args(["index", path]).assert().success();
}

fn stdout_of<'a>(args: impl IntoIterator<Item = &'a str>) -> String {
    let out = gloci()
        .args(args)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    String::from_utf8(out).unwrap()
}

/// Number of record lines (not `#` lines) printed.
fn count<const N: usize>(args: [&str; N]) -> usize {
    stdout_of(args)
        .lines()
        .filter(|l| !l.starts_with('#'))
        .count()
}

fn records(text: &str) -> String {
    text.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| format!("{l}\n"))
        .collect()
}

fn rec(chrom: &str, start: u64, end: u64, name: &str) -> TestRec {
    TestRec {
        chrom: chrom.into(),
        start,
        end,
        strand: '.',
        name: name.into(),
    }
}

fn bed_text(recs: &[TestRec]) -> String {
    recs.iter().map(|r| r.bed_line() + "\n").collect()
}

/// Records around two overlapping regions, with points on their edges.
fn idx_recs() -> Vec<TestRec> {
    vec![
        rec("chr1", 50, 120, "a"),
        rec("chr1", 99, 99, "p0"),
        rec("chr1", 100, 100, "p1"),
        rec("chr1", 140, 160, "b"),
        rec("chr1", 149, 149, "p2"),
        rec("chr1", 199, 250, "c"),
        rec("chr1", 200, 200, "p3"),
        rec("chr1", 299, 300, "d"),
        rec("chr1", 300, 400, "e"),
        rec("chr2", 150, 160, "f"),
    ]
}

/// The records of `recs` that overlap one of the chr1 regions (1-based,
/// closed), each once, in file order.
fn naive_once_sorted(recs: &[TestRec], regions: &[(u64, u64)]) -> String {
    recs.iter()
        .filter(|r| {
            regions.iter().any(|&(s, e)| {
                naive::matches(
                    r,
                    b"chr1",
                    Interval {
                        start: s - 1,
                        end: e,
                    },
                    false,
                )
            })
        })
        .map(|r| r.bed_line() + "\n")
        .collect()
}

const VCF_HEAD: &str = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=100000>\n\
##contig=<ID=chr2,length=100000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

#[test]
fn insertion_site_inside_feature() {
    // the user's §4.3 example, every access path
    let fx = Fx::new();
    for f in fx.every_path("chr1\t10\t10\tins\n") {
        assert_eq!(f.view(&["chr1:6-15"]), "chr1\t10\t10\tins\n");
        // point on left edge of [10,15)
        assert_eq!(f.view(&["chr1:11-15"]), "");
        assert_eq!(
            f.view(&["chr1:11-15", "--inclusive-points"]),
            "chr1\t10\t10\tins\n"
        );
    }
}

#[test]
fn point_on_right_edge_and_equal_points() {
    // A point at a region's end is found with --inclusive-points, and two
    // equal points overlap: scanning must not stop at key == region end.
    let fx = Fx::new();
    let text = "chr1\t5\t9\ta\nchr1\t15\t15\tp\nchr1\t15\t20\tb\nchr1\t30\t40\tc\n";
    for f in fx.every_path(text) {
        assert_eq!(f.view(&["chr1:11-15"]), "");
        assert_eq!(
            f.view(&["chr1:11-15", "--inclusive-points"]),
            "chr1\t15\t15\tp\n"
        );
        // the point between bases 15 and 16 (1-based), as a region
        assert_eq!(f.view(&["chr1:16-15"]), "chr1\t15\t15\tp\n");
        let r = fx.plain_bed("chr1\t15\t15\n");
        assert_eq!(f.view(&["--regions-file", &r]), "chr1\t15\t15\tp\n");
    }
}

#[test]
fn vcf_span_controls_region_hits() {
    // insertion after base 200 vs chr1:100-200: its allele footprint is the
    // point 200, on the region's right edge, so it is not found; its ref
    // footprint [199,200) is.
    let fx = Fx::new();
    let text = format!(
        "{VCF_HEAD}chr1\t150\t.\tC\tT\t.\t.\t.\nchr1\t200\tins\tA\tACGT\t.\t.\t.\nchr1\t500\t.\tG\tA\t.\t.\t.\n"
    );
    let snv = "chr1\t150\t.\tC\tT\t.\t.\t.\n";
    let ins = "chr1\t200\tins\tA\tACGT\t.\t.\t.\n";
    let plain = fx.write(".vcf", text.as_bytes());
    for f in [fx.vcf_gz(&text), fx.vcf_gz_indexed(&text), plain] {
        assert_eq!(records(&stdout_of(["view", &f, "chr1:100-200"])), snv);
        assert_eq!(
            records(&stdout_of([
                "view",
                &f,
                "chr1:100-200",
                "--vcf-span",
                "ref"
            ])),
            format!("{snv}{ins}")
        );
        assert_eq!(
            records(&stdout_of([
                "view",
                &f,
                "chr1:100-200",
                "--inclusive-points"
            ])),
            format!("{snv}{ins}")
        );
        // the point itself is inside chr1:200-201
        assert_eq!(records(&stdout_of(["view", &f, "chr1:200-201"])), ins);
    }
}

#[test]
fn overlapping_regions_once_in_order() {
    let fx = Fx::new();
    let recs = idx_recs();
    let text = bed_text(&recs);
    let expected = naive_once_sorted(&recs, &[(150, 300), (100, 200)]);
    assert!(expected.contains("\tp1\t") && expected.contains("\tp3\t"));
    assert!(!expected.contains("\tp0\t") && !expected.contains("\te\t"));
    let idx = fx.bgzf_bed_indexed(&text);
    assert_eq!(
        stdout_of(["view", &idx, "chr1:150-300", "chr1:100-200"]),
        expected
    );
    for f in fx.every_path(&text) {
        assert_eq!(f.view(&["chr1:150-300", "chr1:100-200"]), expected);
        // repeated and nested regions change nothing
        assert_eq!(
            f.view(&[
                "chr1:100-200",
                "chr1:150-300",
                "chr1:160-170",
                "chr1:100-200"
            ]),
            expected
        );
    }
}

#[test]
fn regions_on_several_chromosomes_in_file_order() {
    let fx = Fx::new();
    let text = "chr1\t1\t5\ta\nchr1\t10\t20\tb\nchr2\t1\t5\tc\nchr2\t30\t40\td\nchr10\t5\t8\te\nchrX\t3\t4\tf\n";
    for f in fx.every_path(text) {
        assert_eq!(
            f.view(&["chrX", "chr2:31-35", "chr1:2-3", "chr10:1-100"]),
            "chr1\t1\t5\ta\nchr2\t30\t40\td\nchr10\t5\t8\te\nchrX\t3\t4\tf\n"
        );
        // a whole chromosome, and a region open to the chromosome end
        assert_eq!(f.view(&["chr2"]), "chr2\t1\t5\tc\nchr2\t30\t40\td\n");
        assert_eq!(f.view(&["chr2:6-"]), "chr2\t30\t40\td\n");
    }
}

#[test]
fn long_feature_starting_before_region_found_unindexed() {
    let fx = Fx::new();
    let text = "chr1\t0\t1000000\tbig\nchr1\t500\t600\tsmall\nchr1\t999990\t999999\tx\n";
    let f = fx.plain_bed(text);
    assert_eq!(
        stdout_of(["view", &f, "chr1:999991-999995", "--assume-sorted"]),
        "chr1\t0\t1000000\tbig\nchr1\t999990\t999999\tx\n"
    );
    for f in fx.every_path(text) {
        assert_eq!(
            f.view(&["chr1:999991-999995"]),
            "chr1\t0\t1000000\tbig\nchr1\t999990\t999999\tx\n"
        );
    }
}

/// A file in lexicographic chromosome order (chr1, chr10, chr2, ..., chr9,
/// chrM, chrX, chrY), as sort -k1,1 or sort-bed writes it. Karyotype
/// bisection for chr10 lands on chrM: chr2..chr9 sort before chr10.
fn lex_sorted_text() -> String {
    let mut names: Vec<String> = (1..=10).map(|i| format!("chr{i}")).collect();
    names.extend(["chrX".into(), "chrY".into(), "chrM".into()]);
    names.sort();
    let mut text = String::new();
    for n in &names {
        if n == "chr10" {
            text.push_str(CHR10_RECS);
            text.push_str("chr10\t200\t300\tfar\n");
        } else {
            for s in [5, 60, 200] {
                text.push_str(&format!("{n}\t{s}\t{}\tr\n", s + 10));
            }
        }
    }
    text
}

const CHR10_RECS: &str = "chr10\t5\t50\tx\nchr10\t60\t70\ty\n";

#[test]
fn other_chrom_order_without_flag_still_correct() {
    // file in lexicographic order (chr1, chr10, chr2, …) queried in
    // karyotype mode: bisection may miss → full scan, correct answer + warning
    let fx = Fx::new();
    let lex_sorted = fx.plain_bed(&lex_sorted_text());
    gloci()
        .args([
            "view",
            &lex_sorted,
            "chr10:1-100",
            "--assume-sorted",
            "--chrom-order",
            "karyotype",
        ])
        .assert()
        .success()
        .stdout(CHR10_RECS)
        .stderr(contains(
            "warning: could not locate chr10 by binary search; scanning the whole file (sort it with gloci sort, or index it)",
        ));
    // the default is lexicographic: bisection hits
    gloci()
        .args(["view", &lex_sorted, "chr10:1-100", "--assume-sorted"])
        .assert()
        .success()
        .stdout(CHR10_RECS)
        .stderr(contains("scanning").not());
    // without --assume-sorted the whole file is scanned and its order
    // checked: a lexicographically sorted file passes
    gloci()
        .args(["view", &lex_sorted, "chr10:1-100"])
        .assert()
        .success()
        .stdout(CHR10_RECS)
        .stderr(contains("no index for"));
}

#[test]
fn genome_order_for_bisection() {
    // A file in -g order (chr2 before chr1) is searched in that order.
    let fx = Fx::new();
    let g = fx.write(".genome", b"chr2\t1000\nchr1\t1000\n");
    let f = fx.plain_bed("chr2\t1\t5\ta\nchr2\t8\t9\tb\nchr1\t1\t5\tc\nchr1\t7\t9\td\n");
    gloci()
        .args([
            "view",
            &f,
            "chr1:6-9",
            "chr2:1-2",
            "-g",
            &g,
            "--assume-sorted",
        ])
        .assert()
        .success()
        .stdout("chr2\t1\t5\ta\nchr1\t7\t9\td\n")
        .stderr(contains("scanning").not());
    // a chromosome missing from -g is reported and yields nothing
    gloci()
        .args([
            "view",
            &f,
            "chr3:1-100",
            "chr1:1-2",
            "-g",
            &g,
            "--assume-sorted",
        ])
        .assert()
        .success()
        .stdout("chr1\t1\t5\tc\n")
        .stderr(contains(format!("chromosome \"chr3\" not found in {g}")));
}

#[test]
fn unsorted_plain_bed_errors_not_wrong() {
    let fx = Fx::new();
    let unsorted = fx.plain_bed("chr1\t50\t60\nchr1\t10\t20\n");
    gloci()
        .args(["view", &unsorted, "chr1:1-100"])
        .assert()
        .code(1)
        .stderr(contains("is not sorted"));
    gloci()
        .args(["view", &unsorted, "chr1:1-100", "--assume-sorted"])
        .assert()
        .code(1)
        .stderr(contains("is not sorted"));
    // a chromosome that comes back (as in a BED with several track
    // sections) is an error by default ...
    let back = "chr1\t10\t20\ta\nchr2\t5\t6\tb\nchr1\t30\t40\tc\n";
    let tracks =
        "track name=one\nchr1\t10\t20\ta\nchr2\t5\t6\tb\ntrack name=two\nchr1\t30\t40\tc\n";
    for text in [back, tracks] {
        let f = fx.plain_bed(text);
        gloci()
            .args(["view", &f, "chr1:1-100"])
            .assert()
            .code(1)
            .stderr(contains("is not sorted"))
            .stderr(contains("line"));
    }
    // ... while --assume-sorted reads and checks only the chr1 block that
    // binary search finds, as documented: the later chr1 record is missed
    let f = fx.plain_bed(back);
    gloci()
        .args(["view", &f, "chr1:1-100", "--assume-sorted"])
        .assert()
        .success()
        .stdout("chr1\t10\t20\ta\n");
    let f = fx.plain_bed(tracks);
    gloci()
        .args(["view", &f, "chr1:1-100", "--assume-sorted"])
        .assert()
        .success()
        .stdout("track name=one\nchr1\t10\t20\ta\n");
    // scanned BGZF files are checked too: out of order, and a chromosome
    // that comes back after another one
    for text in [
        "chr1\t50\t60\nchr1\t10\t20\n",
        "chr1\t1\t2\nchr2\t1\t2\nchr1\t5\t6\n",
    ] {
        let gz = fx.bgzf_bed_unindexed(text);
        gloci()
            .args(["view", &gz, "chr1:1-100"])
            .assert()
            .code(1)
            .stderr(contains("is not sorted"));
    }
}

#[test]
fn order_errors_name_the_real_line() {
    // Scanned, or with the chr2 block reached by binary search, errors name
    // the line's number in the whole file.
    let fx = Fx::new();
    let mut text = String::from("track name=t\n");
    for s in 0..50 {
        text.push_str(&format!("chr1\t{s}\t{}\n", s + 1));
    }
    text.push_str("chr2\t10\t20\nchr2\t30\t40\nchr2\t5\t6\nchr3\t1\t2\n");
    let f = fx.plain_bed(&text);
    let bad = fx.plain_bed(&text.replace("chr2\t30\t40", "chr2\t3x\t40"));
    let skip = fx.plain_bed(&text.replace("chr2\t5\t6", "chr2\t55\t5x"));
    for flags in [&[][..], &["--assume-sorted"][..]] {
        let run = |file: &str, extra: &[&str]| {
            let mut c = gloci();
            c.args(["view", file, "chr2:1-100"]).args(flags).args(extra);
            c.assert()
        };
        let a = run(&f, &[]).code(1).stderr(contains(
            "line 54 (chr2\t5\t6) comes after line 53 (chr2\t30\t40)",
        ));
        if !flags.is_empty() {
            a.stderr(contains("scanning").not());
        }
        // an invalid record
        run(&bad, &[])
            .code(1)
            .stderr(contains("line 53: start \"3x\" is not an integer"));
        run(&bad, &["--skip-invalid"]).code(1).stderr(contains(
            "line 54 (chr2\t5\t6) comes after line 52 (chr2\t10\t20)",
        ));
        run(&skip, &["--skip-invalid"])
            .success()
            .stdout("track name=t\nchr2\t10\t20\nchr2\t30\t40\n")
            .stderr(contains("skipped 1 invalid records"))
            .stderr(contains("line 54: end \"5x\" is not an integer"));
    }
}

#[test]
fn no_index_warns_and_scans() {
    let fx = Fx::new();
    let text = format!("{VCF_HEAD}chr1\t5\t.\tA\tC\t.\t.\t.\nchr1\t5000\t.\tA\tC\t.\t.\t.\n");
    let vcf_gz_no_index = fx.vcf_gz(&text);
    gloci()
        .args(["view", &vcf_gz_no_index, "chr1:1-1000"])
        .assert()
        .success()
        .stderr(contains("no index for"));
    gloci()
        .args(["view", &vcf_gz_no_index, "chr1:1-1000"])
        .assert()
        .stdout(contains("##glociCommand="))
        .stdout(contains("chr1\t5\t.\tA\tC").and(contains("chr1\t5000\t").not()))
        .stderr(contains(format!(
            "warning: no index for {vcf_gz_no_index}; scanning the whole file (create one with: gloci index {vcf_gz_no_index})"
        )));
    // plain BED is scanned unless --assume-sorted
    let bed = fx.plain_bed("chr1\t1\t5\n");
    gloci()
        .args(["view", &bed, "chr1:1-3"])
        .assert()
        .success()
        .stdout("chr1\t1\t5\n")
        .stderr(format!(
            "warning: no index for {bed}; scanning the whole file (sort and index it with: \
gloci sort {bed} -o {bed}.gz --write-index; or, if it is sorted, pass --assume-sorted to binary-search it)\n"
        ));
    gloci()
        .args(["view", &bed, "chr1:1-3", "--assume-sorted"])
        .assert()
        .success()
        .stdout("chr1\t1\t5\n")
        .stderr("");
    // an indexed file needs no warning
    let indexed = fx.vcf_gz_indexed(&text);
    gloci()
        .args(["view", &indexed, "chr1:1-1000"])
        .assert()
        .success()
        .stderr("");
}

/// A pipe or a descriptor path cannot have an index, as stdin cannot: it is
/// scanned as stdin is, without a fix command naming it (`/dev/stdin.gz`).
#[cfg(unix)]
#[test]
fn pipes_and_descriptors_are_scanned_without_an_index_hint() {
    for path in ["-", "/dev/stdin"] {
        gloci()
            .args(["view", path, "chr1:1-10"])
            .write_stdin("chr1\t1\t5\nchr2\t1\t5\n")
            .assert()
            .success()
            .stdout("chr1\t1\t5\n")
            .stderr("");
    }
}

#[test]
fn stale_index_warns() {
    // D26: an index older than its data file is used, with a warning.
    let fx = Fx::new();
    let f = fx.bgzf_bed_indexed("chr1\t1\t5\n");
    let tbi = format!("{f}.tbi");
    let old = std::time::SystemTime::now() - std::time::Duration::from_secs(3600);
    std::fs::File::options()
        .write(true)
        .open(&tbi)
        .unwrap()
        .set_modified(old)
        .unwrap();
    gloci()
        .args(["view", &f, "chr1:1-10"])
        .assert()
        .success()
        .stdout("chr1\t1\t5\n")
        .stderr(contains(format!("the index {tbi} is older than {f}")));
}

#[test]
fn chr_prefix_hint() {
    let fx = Fx::new();
    let idx_chr = fx.bgzf_bed_indexed("chr1\t1\t50\nchr2\t1\t50\n");
    gloci()
        .args(["view", &idx_chr, "1:1-100"])
        .assert()
        .success()
        .stderr(contains("did you mean \"chr1\""));
    gloci()
        .args(["view", &idx_chr, "1:1-100", "chr2:1-10"])
        .assert()
        .success()
        .stdout("chr2\t1\t50\n")
        .stderr(contains(format!(
            "warning: chromosome \"1\" not found in {idx_chr} (did you mean \"chr1\"?)"
        )));
    // the other way round, from VCF contigs, and from a scan of plain BED
    let vcf = fx.vcf_gz("##fileformat=VCFv4.3\n##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t5\t.\tA\tC\t.\t.\t.\n");
    gloci()
        .args(["view", &vcf, "chr1"])
        .assert()
        .success()
        .stderr(contains("chromosome \"chr1\" not found").and(contains("did you mean \"1\"")));
    let bed = fx.plain_bed("chr1\t1\t50\n");
    gloci()
        .args(["view", &bed, "1:1-100"])
        .assert()
        .success()
        .stdout("")
        .stderr(contains(format!(
            "chromosome \"1\" not found in {bed} (did you mean \"chr1\"?)"
        )));
}

#[test]
fn colon_contig_names() {
    // Review Focus 4
    let fx = Fx::new();
    let text = "##fileformat=VCFv4.3\n\
##contig=<ID=chr6,length=100000>\n\
##contig=<ID=HLA-A*01:01:01:01,length=3503>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr6\t50\t.\tA\tC\t.\t.\t.\n\
HLA-A*01:01:01:01\t50\t.\tA\tC\t.\t.\t.\n\
HLA-A*01:01:01:01\t500\t.\tA\tC\t.\t.\t.\n\
HLA-A*01:01:01:01\t3000\t.\tA\tC\t.\t.\t.\n";
    for hla_vcf_gz in [fx.vcf_gz(text), fx.vcf_gz_indexed(text)] {
        assert_eq!(count(["view", &hla_vcf_gz, "HLA-A*01:01:01:01"]), 3);
        assert_eq!(count(["view", &hla_vcf_gz, "{HLA-A*01:01:01:01}:1-100"]), 1);
        assert_eq!(
            count(["view", &hla_vcf_gz, "HLA-A*01:01:01:01:400-3,000"]),
            2
        );
    }
}

#[test]
fn regions_file() {
    // --regions-file r.bed (0-based, half-open) equals positional regions
    let fx = Fx::new();
    let recs = idx_recs();
    let text = bed_text(&recs);
    let expected = naive_once_sorted(&recs, &[(150, 300), (100, 200)]);
    let r = fx.plain_bed("track name=regions\nchr1\t149\t300\nchr1\t99\t200\tsecond\n");
    let r_gz = fx.bgzf_bed_unindexed("chr1\t149\t300\nchr1\t99\t200\n");
    for f in fx.every_path(&text) {
        assert_eq!(f.view(&["chr1:150-300", "chr1:100-200"]), expected);
        assert_eq!(f.view(&["--regions-file", &r]), expected);
        assert_eq!(f.view(&["--regions-file", &r_gz]), expected);
        // combined with positional regions
        assert_eq!(
            f.view(&["--regions-file", &r, "chr2:1-1000"]),
            format!("{expected}chr2\t150\t160\tf\t0\t.\n")
        );
        // an empty regions file selects nothing
        let empty = fx.plain_bed("");
        assert_eq!(f.view(&["--regions-file", &empty]), "");
    }
    // a bad regions file is an error naming its line
    let bad = fx.plain_bed("chr1\t10\t5\n");
    gloci()
        .args(["view", &fx.plain_bed(&text), "--regions-file", &bad])
        .assert()
        .code(1)
        .stderr(contains(format!("{bad} line 1")));
}

#[test]
fn headers_and_conversion_with_regions() {
    let fx = Fx::new();
    // BED header lines are kept, as in plain view
    for f in fx.every_path("track name=t\nchr1\t1\t5\nchr1\t10\t20\n") {
        assert_eq!(f.view(&["chr1:11-12"]), "track name=t\nchr1\t10\t20\n");
        assert_eq!(f.view(&["chr1:11-12", "--no-header"]), "chr1\t10\t20\n");
    }
    // VCF to BCF and back, restricted to a region
    let text = format!("{VCF_HEAD}chr1\t5\t.\tA\tC\t.\t.\t.\nchr2\t7\t.\tG\tT\t.\t.\t.\n");
    let bcf = fx.path(".bcf");
    gloci()
        .args(["view", &fx.vcf_gz_indexed(&text), "chr2", "-o"])
        .arg(&bcf)
        .arg("--write-index")
        .assert()
        .success();
    let bcf = bcf.to_string_lossy().into_owned();
    assert_eq!(
        records(&stdout_of(["view", &bcf, "-O", "vcf"])),
        "chr2\t7\t.\tG\tT\t.\t.\t.\n"
    );
    let full = fx.write(".bcf", &fixtures::vcf_to_bcf(text.as_bytes()));
    index(&full);
    for f in [
        full.clone(),
        fx.write(".bcf", &fixtures::vcf_to_bcf(text.as_bytes())),
    ] {
        let out = stdout_of(["view", &f, "chr1:1-10", "-O", "vcf"]);
        assert!(out.contains("##glociCommand="));
        assert_eq!(records(&out), "chr1\t5\t.\tA\tC\t.\t.\t.\n");
    }
}

#[test]
fn bad_regions_are_usage_errors() {
    let fx = Fx::new();
    let f = fx.plain_bed("chr1\t1\t5\n");
    for bad in ["chr1:0-5", "chr1:x-5", "chr1:9-5", ":5-6", "{chr1:1-5"] {
        gloci()
            .args(["view", &f, bad])
            .assert()
            .code(2)
            .stderr(contains("region"));
    }
}

#[test]
fn stdin_and_header_only_with_regions() {
    let fx = Fx::new();
    gloci()
        .args(["view", "-", "chr1:2-3"])
        .write_stdin("chr1\t1\t5\nchr1\t10\t20\n")
        .assert()
        .success()
        .stdout("chr1\t1\t5\n")
        .stderr("");
    let f = fx.plain_bed("track name=t\nchr1\t1\t5\n");
    assert_eq!(
        stdout_of(["view", &f, "chr1:2-3", "--header-only"]),
        "track name=t\n"
    );
}

#[test]
fn region_flag_equals_positional() {
    // spec §6.1: --region chr:s-e, repeatable, means the same as REGION
    let fx = Fx::new();
    let recs = idx_recs();
    let expected = naive_once_sorted(&recs, &[(150, 300), (100, 200)]);
    for f in fx.every_path(&bed_text(&recs)) {
        assert_eq!(
            f.view(&["--region", "chr1:150-300", "--region", "chr1:100-200"]),
            expected
        );
        assert_eq!(
            f.view(&["chr1:150-300", "--region", "chr1:100-200"]),
            expected
        );
        let r = fx.plain_bed("chr1\t149\t300\n");
        assert_eq!(
            f.view(&["--regions-file", &r, "--region", "chr1:100-200"]),
            expected
        );
    }
    let f = fx.plain_bed("chr1\t1\t5\n");
    gloci()
        .args(["view", &f, "--region", "chr1:0-5"])
        .assert()
        .code(2)
        .stderr(contains("invalid region"));
}

#[test]
fn htslib_bcf_index_with_contigs_out_of_dictionary_order() {
    // bcftools indexes a BCF whose contigs are contiguous but not in
    // dictionary order (chr2, id 1, before chr1, id 0): a query on both
    // returns both, in file order.
    let fx = Fx::new();
    let text = format!(
        "{VCF_HEAD}chr2\t5\t.\tA\tC\t.\t.\t.\nchr2\t9\t.\tG\tT\t.\t.\t.\nchr1\t7\t.\tA\tG\t.\t.\t.\n"
    );
    let bcf = fixtures::vcf_to_bcf(text.as_bytes());
    let f = fx.write(".bcf", &bcf);
    std::fs::write(format!("{f}.csi"), fixtures::bcf_csi_any_order(&bcf)).unwrap();
    let both = "chr2\t5\t.\tA\tC\t.\t.\t.\nchr2\t9\t.\tG\tT\t.\t.\t.\nchr1\t7\t.\tA\tG\t.\t.\t.\n";
    for args in [["chr1", "chr2"], ["chr2", "chr1"]] {
        gloci()
            .args(["view", &f, args[0], args[1], "-O", "vcf"])
            .assert()
            .success()
            .stdout(contains(both))
            .stderr("");
    }
    assert_eq!(
        records(&stdout_of(["view", &f, "chr1", "-O", "vcf"])),
        "chr1\t7\t.\tA\tG\t.\t.\t.\n"
    );
}
