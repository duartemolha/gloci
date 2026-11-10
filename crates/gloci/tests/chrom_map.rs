//! `--chrom-map` (D62): a rename table applied as records are read.
mod common;
use common::*;
use predicates::str::contains;
use std::path::Path;

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn stdout_of(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

const MAP: &str = "# the UCSC names\n1\tchr1\n2\tchr2\n";

const VCF_HEAD: &str = "##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n\
##contig=<ID=2,length=1000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n";
const VCF_ROWS: &str = "1\t10\t.\tA\tC\t50\tPASS\t.\tGT\t0/1\n\
1\t30\t.\tG\tT\t50\tPASS\t.\tGT\t1/1\n2\t5\t.\tC\tA\t50\tPASS\t.\tGT\t0/0\n";

struct Fx {
    dir: tempfile::TempDir,
}

impl Fx {
    fn new() -> Fx {
        Fx { dir: tempdir() }
    }
    fn plain(&self, name: &str, text: &str) -> String {
        s(&plain_file(self.dir.path(), name, text))
    }
    fn map(&self) -> String {
        self.plain("map.tsv", MAP)
    }
    fn indexed(&self, name: &str, text: &str) -> String {
        let p = s(&bgzf_file(self.dir.path(), name, text));
        gloci().args(["index", &p]).assert().success();
        p
    }
    fn vcf(&self) -> String {
        self.plain("in.vcf", &format!("{VCF_HEAD}{VCF_ROWS}"))
    }
    /// The VCF as a BCF with the file's own (unrenamed) names.
    fn bcf(&self) -> String {
        let p = s(&self.dir.path().join("in.bcf"));
        gloci()
            .args(["view", &self.vcf(), "-o", &p])
            .assert()
            .success();
        p
    }
}

#[test]
fn chrom_map_everywhere() {
    let fx = Fx::new();
    let map = fx.map();
    let bed = fx.plain("a.bed", "2\t5\t8\n1\t30\t40\n1\t10\t20\n");
    let b = fx.plain("b.bed", "1\t15\t35\n2\t1\t6\n");

    // sort: the renamed names are what -g sees and what is written
    assert_eq!(
        stdout_of(&["sort", &bed, "--chrom-map", &map, "-g", "hg38"]),
        "chr1\t10\t20\nchr1\t30\t40\nchr2\t5\t8\n"
    );
    // without the map the names are not in the genome
    gloci()
        .args(["sort", &bed, "-g", "hg38"])
        .assert()
        .failure();

    // intersect: both inputs are renamed, and pass the order checks
    let sorted = fx.plain("s.bed", "1\t10\t20\n1\t30\t40\n2\t5\t8\n");
    assert_eq!(
        stdout_of(&[
            "intersect",
            "-a",
            &sorted,
            "-b",
            &b,
            "--chrom-map",
            &map,
            "-g",
            "hg38"
        ]),
        "chr1\t10\t20\nchr1\t30\t40\nchr2\t5\t8\n"
    );
    // a mix of names: A already says chr1, B says 1
    let a_chr = fx.plain("a_chr.bed", "chr1\t10\t20\n");
    assert_eq!(
        stdout_of(&["intersect", "-a", &a_chr, "-b", &b, "--chrom-map", &map]),
        "chr1\t10\t20\n"
    );

    // merge
    assert_eq!(
        stdout_of(&["merge", &sorted, "--chrom-map", &map, "-g", "hg38"]),
        "chr1\t10\t20\nchr1\t30\t40\nchr2\t5\t8\n"
    );

    // view over an index built on the file's own names: the region is
    // written with the new name, and is found
    let idx = fx.indexed("idx.bed.gz", "1\t10\t20\n1\t30\t40\n2\t5\t8\n");
    assert_eq!(
        stdout_of(&["view", &idx, "chr1:25-100", "--chrom-map", &map]),
        "chr1\t30\t40\n"
    );
    assert_eq!(
        stdout_of(&["view", &idx, "chr2", "--chrom-map", &map]),
        "chr2\t5\t8\n"
    );
    // the file's own name no longer names the chromosome
    gloci()
        .args(["view", &idx, "1:1-100", "--chrom-map", &map])
        .assert()
        .success()
        .stdout("");
    // an unmapped name still works
    let idx2 = fx.indexed("idx2.bed.gz", "3\t1\t9\n");
    assert_eq!(
        stdout_of(&["view", &idx2, "3", "--chrom-map", &map]),
        "3\t1\t9\n"
    );

    // the same on a regions file, which is never renamed
    let regions = fx.plain("r.bed", "chr1\t25\t100\n");
    assert_eq!(
        stdout_of(&[
            "view",
            &idx,
            "--regions-file",
            &regions,
            "--chrom-map",
            &map
        ]),
        "chr1\t30\t40\n"
    );

    // VCF filter on VCF and BCF
    for input in [fx.vcf(), fx.bcf()] {
        let out = stdout_of(&["vcf", "filter", &input, "--chrom-map", &map, "-O", "vcf"]);
        let recs: Vec<&str> = out.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(recs.len(), 3, "{out}");
        assert!(recs[0].starts_with("chr1\t10\t"), "{out}");
        assert!(recs[2].starts_with("chr2\t5\t"), "{out}");
        assert!(out.contains("##contig=<ID=chr1,length=1000>"), "{out}");
        assert!(out.contains("##contig=<ID=chr2,length=1000>"), "{out}");
        assert!(!out.contains("##contig=<ID=1,"), "{out}");
    }

    // view region on an indexed VCF and BCF
    let vcf_gz = fx.indexed("v.vcf.gz", &format!("{VCF_HEAD}{VCF_ROWS}"));
    let out = stdout_of(&["view", &vcf_gz, "chr1:20-100", "--chrom-map", &map]);
    let recs: Vec<&str> = out.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(recs, ["chr1\t30\t.\tG\tT\t50\tPASS\t.\tGT\t1/1"], "{out}");
    let bcf = fx.bcf();
    gloci().args(["index", &bcf]).assert().success();
    let out = stdout_of(&[
        "view",
        &bcf,
        "chr1:20-100",
        "--chrom-map",
        &map,
        "-O",
        "vcf",
    ]);
    let recs: Vec<&str> = out.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(recs, ["chr1\t30\t.\tG\tT\t50\tPASS\t.\tGT\t1/1"], "{out}");
}

#[test]
fn indexed_relations_are_thread_deterministic_with_a_map() {
    let fx = Fx::new();
    let map = fx.map();
    let mut a = String::new();
    let mut b = String::new();
    for c in ["1", "2"] {
        for i in 0..200u64 {
            a.push_str(&format!("{c}\t{}\t{}\n", i * 10, i * 10 + 6));
            b.push_str(&format!("{c}\t{}\t{}\n", i * 10 + 3, i * 10 + 9));
        }
    }
    let (pa, pb) = (fx.plain("pa.bed", &a), fx.plain("pb.bed", &b));
    let (a, b) = (fx.indexed("a.bed.gz", &a), fx.indexed("b.bed.gz", &b));
    let run = |t: &str| {
        stdout_of(&[
            "intersect",
            "-a",
            &a,
            "-b",
            &b,
            "--chrom-map",
            &map,
            "-g",
            "hg38",
            "-t",
            t,
        ])
    };
    let one = run("1");
    assert_eq!(one.lines().count(), 400);
    assert!(one.starts_with("chr1\t0\t6\n") && one.ends_with("chr2\t1990\t1996\n"));
    assert_eq!(one, run("4"));
    // the streaming path agrees
    let streamed = stdout_of(&[
        "intersect",
        "-a",
        &pa,
        "-b",
        &pb,
        "--chrom-map",
        &map,
        "-g",
        "hg38",
        "-t",
        "1",
    ]);
    assert_eq!(one, streamed);
}

#[test]
fn chrom_map_errors() {
    let fx = Fx::new();
    let bed = fx.plain("a.bed", "1\t10\t20\n");
    let run = |text: &str| {
        let m = fx.plain("bad.tsv", text);
        let out = gloci()
            .args(["sort", &bed, "--chrom-map", &m])
            .assert()
            .code(2);
        (
            m,
            String::from_utf8(out.get_output().stderr.clone()).unwrap(),
        )
    };
    // a duplicate old name
    let (m, err) = run("1\tchr1\n2\tchr2\n1\tchrX\n");
    assert!(
        err.contains(&format!("{m}:3:")) && err.contains("line 1"),
        "{err}"
    );
    // a duplicate new name
    let (m, err) = run("1\tchr1\n2\tchr1\n");
    assert!(
        err.contains(&format!("{m}:2:"))
            && err.contains("\"chr1\" is already the new name on line 1"),
        "{err}"
    );
    // malformed lines name the line
    let (m, err) = run("1\tchr1\n2 chr2\n");
    assert!(
        err.contains(&format!("{m}:2: expected 2 tab-separated columns")),
        "{err}"
    );
    let (m, err) = run("1\tchr1\n\n3\t\n");
    assert!(
        err.contains(&format!("{m}:3: empty chromosome name")),
        "{err}"
    );
    let (m, err) = run("1\tchr1\tx\n");
    assert!(
        err.contains(&format!("{m}:1:")) && err.contains("found 3"),
        "{err}"
    );
    // the table is read even when the command fails for another reason
    gloci()
        .args(["sort", &bed, "--chrom-map", "/nonexistent/map.tsv"])
        .assert()
        .code(2)
        .stderr(contains(
            "cannot read the --chrom-map file /nonexistent/map.tsv",
        ));
}

#[test]
fn chrom_map_accepts_crlf_blank_lines_and_comments() {
    let fx = Fx::new();
    let bed = fx.plain("a.bed", "1\t10\t20\n");
    let m = fx.plain("m.tsv", "# c\r\n\r\n1\tchr1\r\n\n");
    assert_eq!(
        stdout_of(&["sort", &bed, "--chrom-map", &m]),
        "chr1\t10\t20\n"
    );
}

#[test]
fn genome_names_are_never_renamed() {
    let fx = Fx::new();
    let map = fx.map();
    let bed = fx.plain("a.bed", "1\t10\t20\n");
    // a genome that lists the file's old name does not match the new one
    let old = fx.plain("old.genome", "1\t1000\n");
    gloci()
        .args(["sort", &bed, "--chrom-map", &map, "-g", &old])
        .assert()
        .failure()
        .stderr(contains("chromosome \"chr1\" is not in the genome file"));
    // one that lists the new name does
    let new = fx.plain("new.genome", "chr1\t1000\n");
    assert_eq!(
        stdout_of(&["sort", &bed, "--chrom-map", &map, "-g", &new]),
        "chr1\t10\t20\n"
    );
}

#[test]
fn vcf_contig_header_is_renamed() {
    let fx = Fx::new();
    let map = fx.plain("m.tsv", "1\tchr1\n");
    let out = stdout_of(&["view", &fx.vcf(), "--chrom-map", &map]);
    assert!(out.contains("##contig=<ID=chr1,length=1000>\n"), "{out}");
    assert!(out.contains("##contig=<ID=2,length=1000>\n"), "{out}");
    assert!(out.contains("\nchr1\t10\t.\tA\tC"), "{out}");
    assert!(out.contains("\n2\t5\t.\tC\tA"), "{out}");
    // two contigs given one name by the map is an error naming the fix
    let map = fx.plain("m2.tsv", "1\t2\n");
    gloci()
        .args(["view", &fx.vcf(), "--chrom-map", &map])
        .assert()
        .failure()
        .stderr(contains("the same name \"2\""));
}

#[test]
fn bcf_to_bcf_renames_the_dictionary() {
    let fx = Fx::new();
    let map = fx.map();
    let out = s(&fx.dir.path().join("out.bcf"));
    gloci()
        .args(["view", &fx.bcf(), "--chrom-map", &map, "-o", &out])
        .assert()
        .success();
    // read back without a map: the dictionary, and so every record, has the
    // new names
    let text = stdout_of(&["view", &out, "-O", "vcf"]);
    assert!(text.contains("##contig=<ID=chr1,length=1000>\n"), "{text}");
    assert!(text.contains("##contig=<ID=chr2,length=1000>\n"), "{text}");
    let recs: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(
        recs,
        [
            "chr1\t10\t.\tA\tC\t50\tPASS\t.\tGT\t0/1",
            "chr1\t30\t.\tG\tT\t50\tPASS\t.\tGT\t1/1",
            "chr2\t5\t.\tC\tA\t50\tPASS\t.\tGT\t0/0",
        ]
    );
    // and it is a valid input for -g hg38
    gloci()
        .args(["vcf", "filter", &out, "-O", "vcf"])
        .assert()
        .success();
}

#[test]
fn vcf_to_bcf_conversion_uses_the_renamed_dictionary() {
    let fx = Fx::new();
    let map = fx.map();
    let out = s(&fx.dir.path().join("conv.bcf"));
    gloci()
        .args(["view", &fx.vcf(), "--chrom-map", &map, "-o", &out])
        .assert()
        .success();
    let text = stdout_of(&["view", &out, "-O", "vcf"]);
    assert!(text.contains("##contig=<ID=chr1,length=1000>\n"), "{text}");
    assert!(
        text.contains("\nchr1\t30\t.\tG\tT\t50\tPASS\t.\tGT\t1/1\n"),
        "{text}"
    );
}

#[test]
fn convert_renames() {
    let fx = Fx::new();
    let map = fx.map();
    assert_eq!(
        stdout_of(&["convert", &fx.vcf(), "--chrom-map", &map]),
        "chr1\t9\t10\t.\tA\tC\t50\tPASS\nchr1\t29\t30\t.\tG\tT\t50\tPASS\nchr2\t4\t5\t.\tC\tA\t50\tPASS\n"
    );
}

#[test]
fn gff_records_are_renamed() {
    let fx = Fx::new();
    let map = fx.map();
    let gff = fx.plain(
        "a.gff3",
        "##gff-version 3\n1\tsrc\tgene\t11\t20\t.\t+\t.\tID=g1\n",
    );
    assert_eq!(
        stdout_of(&["sort", &gff, "--chrom-map", &map]),
        "##gff-version 3\nchr1\tsrc\tgene\t11\t20\t.\t+\t.\tID=g1\n"
    );
}

#[test]
fn naming_warning_is_silenced_by_the_map() {
    let fx = Fx::new();
    let map = fx.map();
    let a = fx.plain("a.bed", "chr1\t10\t20\n");
    let b = fx.plain("b.bed", "1\t15\t25\n");
    gloci()
        .args(["intersect", "-a", &a, "-b", &b])
        .assert()
        .success()
        .stderr(contains("--chrom-map FILE"));
    let o = gloci()
        .args(["intersect", "-a", &a, "-b", &b, "--chrom-map", &map])
        .assert()
        .success()
        .stdout("chr1\t10\t20\n");
    assert!(!String::from_utf8_lossy(&o.get_output().stderr).contains("share no chromosome"));
}

#[test]
fn index_does_not_take_the_map() {
    let fx = Fx::new();
    let p = s(&bgzf_file(fx.dir.path(), "x.bed.gz", "1\t1\t2\n"));
    gloci()
        .args(["index", &p, "--chrom-map", &fx.map()])
        .assert()
        .code(2);
}

#[test]
fn index_with_both_old_and_new_name_is_refused() {
    let fx = Fx::new();
    let map = fx.map();
    let idx = fx.indexed("both.bed.gz", "1\t10\t20\nchr1\t30\t40\n");
    gloci()
        .args(["view", &idx, "chr1", "--chrom-map", &map])
        .assert()
        .code(2)
        .stderr(contains("both.bed.gz has both \"1\" and \"chr1\""))
        .stderr(contains("which would merge them"));
}

#[test]
fn assume_sorted_bisection_with_a_map() {
    let fx = Fx::new();
    let map = fx.map();
    let bed = fx.plain("s.bed", "1\t10\t20\n1\t30\t40\n2\t5\t8\n2\t50\t60\n");
    for region in ["chr1:25-100", "chr2", "chr2:1-6", "chr1"] {
        let scan = stdout_of(&["view", &bed, region, "--chrom-map", &map]);
        let bisect = stdout_of(&["view", &bed, region, "--chrom-map", &map, "--assume-sorted"]);
        assert!(!scan.is_empty(), "{region}");
        assert_eq!(scan, bisect, "{region}");
    }
    // a name the map renames away is missing, as on the indexed path
    gloci()
        .args([
            "view",
            &bed,
            "1:1-100",
            "--chrom-map",
            &map,
            "--assume-sorted",
        ])
        .assert()
        .success()
        .stdout("")
        .stderr(contains("could not locate 1 by binary search"));
}

/// A map that renames onto a name already in the file leaves a file that is
/// not sorted in new-name space. `--assume-sorted` trusts the claim that it
/// is, and silently answers from the renamed block only; the scan checks the
/// order and fails. Documented in docs/migration.md (a note, not a bug).
#[test]
fn assume_sorted_trusts_a_map_that_renames_onto_an_existing_name() {
    let fx = Fx::new();
    let bed = fx.plain("h.bed", "a\t10\t20\na\t30\t40\nb\t5\t8\nc\t50\t60\n");
    let map = fx.plain("h.map", "a\tc\n");
    let genome = fx.plain("h.genome", "c\t1000\nb\t1000\n");
    let base = ["view", &bed, "c:1-1000", "--chrom-map", &map, "-g", &genome];
    // the scan reports the unsorted file
    gloci()
        .args(base)
        .assert()
        .code(1)
        .stderr(contains("is not sorted"));
    // the bisection finds the renamed `a` block only; `c` (line 4) is missed
    gloci()
        .args(base)
        .arg("--assume-sorted")
        .assert()
        .success()
        .stdout("c\t10\t20\nc\t30\t40\n");
}

#[test]
fn order_reversing_map_falls_back_to_a_scan_under_assume_sorted() {
    let fx = Fx::new();
    let bed = fx.plain("r.bed", "a\t10\t20\na\t30\t40\nb\t5\t8\n");
    let map = fx.plain("r.map", "a\tz\nb\ty\n");
    for (region, want) in [("y", "y\t5\t8\n"), ("z", "z\t10\t20\nz\t30\t40\n")] {
        gloci()
            .args([
                "view",
                &bed,
                region,
                "--chrom-map",
                &map,
                "--chrom-order",
                "lex",
                "--assume-sorted",
            ])
            .assert()
            .success()
            .stdout(want)
            .stderr(contains("could not locate"))
            .stderr(contains("scanning the whole file"));
    }
}

#[test]
fn unreadable_map_is_a_usage_error() {
    let fx = Fx::new();
    let bed = fx.plain("a.bed", "1\t10\t20\n");
    gloci()
        .args(["sort", &bed, "--chrom-map", "/nonexistent/map.tsv"])
        .assert()
        .code(2)
        .stderr(contains("/nonexistent/map.tsv"))
        .stderr(contains("check the --chrom-map path"));
}

#[test]
fn unsorted_relation_with_a_map() {
    let fx = Fx::new();
    let map = fx.map();
    let a = fx.plain("a.bed", "2\t5\t8\n1\t10\t20\n");
    let b = fx.plain("b.bed", "1\t15\t25\n2\t1\t6\n");
    assert_eq!(
        stdout_of(&[
            "intersect",
            "-a",
            &a,
            "-b",
            &b,
            "--unsorted",
            "--chrom-map",
            &map
        ]),
        "chr2\t5\t8\nchr1\t10\t20\n"
    );
}

#[test]
fn gzip_and_stdin_inputs_with_a_map() {
    let fx = Fx::new();
    let map = fx.map();
    let gz = s(&bgzf_file(
        fx.dir.path(),
        "a.bed.gz",
        "1\t10\t20\n2\t5\t8\n",
    ));
    assert_eq!(
        stdout_of(&["sort", &gz, "--chrom-map", &map]),
        "chr1\t10\t20\nchr2\t5\t8\n"
    );
    let o = gloci()
        .args(["sort", "-", "--chrom-map", &map])
        .write_stdin("2\t5\t8\n1\t10\t20\n")
        .assert()
        .success();
    assert_eq!(
        String::from_utf8_lossy(&o.get_output().stdout),
        "chr1\t10\t20\nchr2\t5\t8\n"
    );
}

#[test]
fn spilling_sort_with_a_map() {
    let fx = Fx::new();
    let map = fx.map();
    let mut text = String::new();
    for i in (0..700_000u64).rev() {
        text.push_str(&format!("1\t{}\t{}\tname{i}\n", i * 10, i * 10 + 5));
    }
    let bed = fx.plain("big.bed", &text);
    let o = gloci()
        .args(["sort", &bed, "--chrom-map", &map, "--max-mem", "16M", "-v"])
        .assert()
        .success();
    let out = String::from_utf8_lossy(&o.get_output().stdout).into_owned();
    assert_eq!(out.lines().count(), 700_000);
    assert!(out.starts_with("chr1\t0\t5\tname0\n"));
    assert!(out.ends_with("chr1\t6999990\t6999995\tname699999\n"));
    assert!(out.lines().all(|l| l.starts_with("chr1\t")));
}

#[test]
fn bcf_to_vcf_with_a_map() {
    let fx = Fx::new();
    let map = fx.map();
    let text = stdout_of(&["view", &fx.bcf(), "--chrom-map", &map, "-O", "vcf"]);
    assert!(text.contains("##contig=<ID=chr1,length=1000>\n"), "{text}");
    assert!(text.contains("##contig=<ID=chr2,length=1000>\n"), "{text}");
    assert!(
        text.contains("\nchr1\t10\t.\tA\tC\t50\tPASS\t.\tGT\t0/1\n"),
        "{text}"
    );
    assert!(
        text.contains("\nchr2\t5\t.\tC\tA\t50\tPASS\t.\tGT\t0/0\n"),
        "{text}"
    );
}
