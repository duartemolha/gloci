mod common;
use common::*;
use gloci_testkit::fixtures;
use predicates::prelude::*;
use predicates::str::contains;
use std::io::Read;
use std::path::Path;
use std::time::Instant;

fn w(d: &Path, name: &str, text: &str) -> String {
    plain_file(d, name, text).to_string_lossy().into_owned()
}

fn stdout_of(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success().get_output().clone();
    String::from_utf8(o.stdout).unwrap()
}

fn bgzf_text(path: &Path) -> String {
    let mut s = String::new();
    noodles_bgzf::io::Reader::new(std::fs::File::open(path).unwrap())
        .read_to_string(&mut s)
        .unwrap();
    s
}

const VCF_HEAD_21: &str =
    "##fileformat=VCFv4.3\n##contig=<ID=chr2,length=1000>\n##contig=<ID=chr1,length=1000>\n";
const VCF_COLS: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

fn vcf(head: &str, cols: &str, recs: &[&str]) -> String {
    let mut s = format!("{head}{cols}");
    for r in recs {
        s.push_str(r);
        s.push('\n');
    }
    s
}

#[test]
fn bed_karyotype_default() {
    let d = tempdir();
    let u = w(
        d.path(),
        "u.bed",
        "chr2\t5\t10\nchrX\t1\t2\nchr10\t1\t2\nchr1\t7\t9\nchrM\t0\t5\nchr1\t3\t4\n",
    );
    gloci()
        .args(["sort", &u])
        .assert()
        .success()
        .stdout("chr1\t3\t4\nchr1\t7\t9\nchr2\t5\t10\nchr10\t1\t2\nchrX\t1\t2\nchrM\t0\t5\n");
    gloci()
        .args(["sort", &u, "--chrom-order", "lex"])
        .assert()
        .success()
        .stdout("chr1\t3\t4\nchr1\t7\t9\nchr10\t1\t2\nchr2\t5\t10\nchrM\t0\t5\nchrX\t1\t2\n");
}

#[test]
fn genome_order() {
    let d = tempdir();
    let g = w(d.path(), "g.genome", "chr1\t100\nchr2\t100\nchr10\t100\n");
    let u = w(
        d.path(),
        "u.bed",
        "chr10\t1\t2\nchr2\t5\t6\nchr1\t9\t10\nchr1\t3\t4\n",
    );
    gloci()
        .args(["sort", &u, "-g", &g])
        .assert()
        .success()
        .stdout("chr1\t3\t4\nchr1\t9\t10\nchr2\t5\t6\nchr10\t1\t2\n");
    let bad = w(d.path(), "bad.bed", "chrZ\t1\t2\n");
    gloci().args(["sort", &bad, "-g", &g]).assert().code(1);
}

#[test]
fn ties_keep_input_order_and_full_key() {
    let d = tempdir();
    let u = w(
        d.path(),
        "t.bed",
        "chr1\t1\t5\tb\nchr1\t1\t5\ta\nchr1\t0\t5\tz\n",
    );
    gloci()
        .args(["sort", &u])
        .assert()
        .success()
        .stdout("chr1\t0\t5\tz\nchr1\t1\t5\tb\nchr1\t1\t5\ta\n");
    gloci()
        .args(["sort", &u, "--full-key"])
        .assert()
        .success()
        .stdout("chr1\t0\t5\tz\nchr1\t1\t5\ta\nchr1\t1\t5\tb\n");
}

#[test]
fn vcf_contig_order_and_provenance() {
    let d = tempdir();
    let v = w(
        d.path(),
        "a.vcf",
        &vcf(
            VCF_HEAD_21,
            VCF_COLS,
            &[
                "chr1\t5\t.\tA\tC\t.\t.\t.",
                "chr2\t9\t.\tA\tC\t.\t.\t.",
                "chr2\t3\t.\tA\tC\t.\t.\t.",
            ],
        ),
    );
    let out = stdout_of(&["sort", &v]);
    assert!(out.contains("##glociCommand=") && out.contains("; Version=0.1.0\n#CHROM"));
    assert!(out.find("\nchr2\t").unwrap() < out.find("\nchr1\t").unwrap());
    assert!(out.find("\nchr2\t3").unwrap() < out.find("\nchr2\t9").unwrap());
}

#[test]
fn vcf_undeclared_contig_is_error() {
    let d = tempdir();
    let v = w(
        d.path(),
        "a.vcf",
        &vcf(VCF_HEAD_21, VCF_COLS, &["chrUn\t5\t.\tA\tC\t.\t.\t."]),
    );
    gloci()
        .args(["sort", &v])
        .assert()
        .code(1)
        .stderr(contains("chromosome \"chrUn\" is not in"));
}

#[test]
fn bcf_roundtrip() {
    let d = tempdir();
    let text = vcf(
        VCF_HEAD_21,
        VCF_COLS,
        &[
            "chr1\t5\t.\tA\tC\t.\t.\t.",
            "chr2\t9\t.\tA\tC\t.\t.\t.",
            "chr1\t2\t.\tG\tT\t.\t.\t.",
            "chr2\t3\t.\tA\tC\t.\t.\t.",
        ],
    );
    let bcf = fixtures::write(d.path(), "a.bcf", &fixtures::vcf_to_bcf(text.as_bytes()));
    let out = d.path().join("s.bcf");
    gloci()
        .args(["sort"])
        .arg(&bcf)
        .arg("-o")
        .arg(&out)
        .assert()
        .success();
    let mut r = noodles_bcf::io::Reader::new(std::fs::File::open(&out).unwrap());
    let header = r.read_header().unwrap();
    assert!(
        header
            .other_records()
            .keys()
            .any(|k| k.as_ref() == "glociCommand")
    );
    let got: Vec<(String, usize)> = r
        .records()
        .map(|rec| {
            let rec = rec.unwrap();
            (
                noodles_vcf::variant::Record::reference_sequence_name(&rec, &header)
                    .unwrap()
                    .to_string(),
                noodles_vcf::variant::Record::variant_start(&rec)
                    .unwrap()
                    .unwrap()
                    .get(),
            )
        })
        .collect();
    let want: Vec<(String, usize)> = [("chr2", 3), ("chr2", 9), ("chr1", 2), ("chr1", 5)]
        .iter()
        .map(|(c, p)| (c.to_string(), *p))
        .collect();
    assert_eq!(got, want);
}

#[test]
fn gff_and_header() {
    let d = tempdir();
    let g = w(
        d.path(),
        "a.gff3",
        "##gff-version 3\nchr2\tsrc\tgene\t5\t9\t.\t+\t.\tID=b\nchr1\tsrc\tgene\t20\t30\t.\t+\t.\tID=a\nchr1\tsrc\tgene\t3\t4\t.\t+\t.\tID=c\n",
    );
    gloci().args(["sort", &g]).assert().success().stdout(
        "##gff-version 3\nchr1\tsrc\tgene\t3\t4\t.\t+\t.\tID=c\nchr1\tsrc\tgene\t20\t30\t.\t+\t.\tID=a\nchr2\tsrc\tgene\t5\t9\t.\t+\t.\tID=b\n",
    );
    let b = w(d.path(), "h.bed", "#h\nchr2\t1\t2\nchr1\t1\t2\n");
    gloci()
        .args(["sort", &b, "--header"])
        .assert()
        .success()
        .stdout("#h\nchr1\t1\t2\nchr2\t1\t2\n");
    gloci()
        .args(["sort", &b])
        .assert()
        .success()
        .stdout("chr1\t1\t2\nchr2\t1\t2\n");
}

#[test]
fn bgzf_with_index() {
    let d = tempdir();
    let u = w(d.path(), "u.bed", "chr2\t5\t10\nchr1\t7\t9\nchr1\t3\t4\n");
    let out = d.path().join("s.bed.gz");
    gloci()
        .args(["sort", &u, "--write-index", "-o"])
        .arg(&out)
        .assert()
        .success()
        .stdout("");
    let bytes = std::fs::read(&out).unwrap();
    assert_eq!(&bytes[..4], &[0x1f, 0x8b, 0x08, 0x04]);
    assert_eq!(bgzf_text(&out), "chr1\t3\t4\nchr1\t7\t9\nchr2\t5\t10\n");
    noodles_tabix::fs::read(d.path().join("s.bed.gz.tbi")).unwrap();
}

#[test]
fn check_mode() {
    let d = tempdir();
    let sorted = w(d.path(), "s.bed", "chr1\t1\t2\nchr1\t5\t6\n");
    let unsorted = w(d.path(), "u.bed", "chr1\t5\t6\nchr1\t1\t2\n");
    gloci()
        .args(["sort", "--check", &sorted])
        .assert()
        .success()
        .stdout("");
    gloci()
        .args(["sort", "--check", &unsorted])
        .assert()
        .code(1)
        .stderr(contains("is not sorted: line 2").and(contains("hint: run: gloci sort")));
    // --check writes nothing, so an output file is a usage error.
    gloci()
        .args(["sort", "--check", &sorted, "-o"])
        .arg(d.path().join("x.bed"))
        .assert()
        .code(2);
}

#[test]
fn merge_sorted_inputs() {
    let d = tempdir();
    let a = w(d.path(), "a.bed", "chr1\t1\t2\nchr1\t9\t10\nchr2\t1\t2\n");
    let b = w(d.path(), "b.bed", "chr1\t5\t6\nchr2\t0\t1\nchr10\t1\t2\n");
    gloci()
        .args(["sort", "-m", &a, &b])
        .assert()
        .success()
        .stdout("chr1\t1\t2\nchr1\t5\t6\nchr1\t9\t10\nchr2\t0\t1\nchr2\t1\t2\nchr10\t1\t2\n");
    let u = w(d.path(), "u.bed", "chr1\t9\t10\nchr1\t1\t2\n");
    gloci()
        .args(["sort", "-m", &a, &u])
        .assert()
        .code(1)
        .stderr(contains("is not sorted"));
    // merge_sorted has no dedup.
    gloci()
        .args(["sort", "-m", "--unique", &a, &b])
        .assert()
        .code(2);
    gloci()
        .args(["sort", "-m", "--duplicates", &a, &b])
        .assert()
        .code(2);
}

#[test]
fn unique_duplicates() {
    let d = tempdir();
    let u = w(
        d.path(),
        "d.bed",
        "chr1\t1\t2\nchr1\t1\t2\nchr1\t1\t2\nchr2\t1\t2\nchr1\t5\t6\nchr1\t5\t6\n",
    );
    gloci()
        .args(["sort", &u, "--unique"])
        .assert()
        .success()
        .stdout("chr1\t1\t2\nchr1\t5\t6\nchr2\t1\t2\n");
    gloci()
        .args(["sort", &u, "--duplicates"])
        .assert()
        .success()
        .stdout("chr1\t1\t2\nchr1\t5\t6\n");
    gloci()
        .args(["sort", &u, "--unique", "--duplicates"])
        .assert()
        .code(2);
}

/// About 50 MB of deterministic pseudo-random BED lines.
fn big_bed() -> String {
    let mut s = String::with_capacity(52 << 20);
    let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut next = move || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    };
    while s.len() < 50 << 20 {
        let c = next() % 24 + 1;
        let start = next() % 200_000_000;
        let len = next() % 5000 + 1;
        s.push_str(&format!(
            "chr{c}\t{start}\t{}\tfeature_{}\t{}\t+\n",
            start + len,
            next() % 1_000_000,
            next() % 1000
        ));
    }
    s
}

#[test]
fn tiny_max_mem_matches_default() {
    let d = tempdir();
    let big = w(d.path(), "big.bed", &big_bed());
    let tmp = d.path().join("spill");
    std::fs::create_dir(&tmp).unwrap();
    let run = |args: &[&str]| {
        let o = gloci()
            .arg("sort")
            .arg(&big)
            .args(args)
            .assert()
            .success()
            .get_output()
            .clone();
        o.stdout
    };
    let small = {
        let o = gloci()
            .arg("sort")
            .arg(&big)
            .args(["--max-mem", "16M", "--tmp-dir"])
            .arg(&tmp)
            .assert()
            .success()
            .get_output()
            .clone();
        o.stdout
    };
    assert_eq!(small, run(&[]));
    assert_eq!(small, run(&["-t", "1", "--max-mem", "16M"]));
    // the spill directory is gone afterwards
    assert_eq!(std::fs::read_dir(&tmp).unwrap().count(), 0);
}

#[test]
fn usage_errors() {
    let d = tempdir();
    let bed = w(d.path(), "a.bed", "chr1\t1\t2\n");
    gloci()
        .args(["sort", &bed, "-o"])
        .arg(d.path().join("x.vcf.gz"))
        .assert()
        .code(2);
    gloci()
        .args(["sort", &bed, "--max-mem", "1M"])
        .assert()
        .code(2);
    gloci()
        .args(["sort", &bed, "--max-mem", "100"])
        .assert()
        .code(2);
    let a = w(
        d.path(),
        "a.vcf",
        &vcf(
            VCF_HEAD_21,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n",
            &[],
        ),
    );
    let c = w(
        d.path(),
        "c.vcf",
        &vcf(
            VCF_HEAD_21,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tC\n",
            &[],
        ),
    );
    gloci()
        .args(["sort", &a, &c])
        .assert()
        .code(1)
        .stderr(contains("sample"));
    // mixed formats
    gloci().args(["sort", &bed, &a]).assert().code(2);
    // an unwritable --tmp-dir fails up front
    gloci()
        .args(["sort", &bed, "--tmp-dir"])
        .arg(d.path().join("missing/dir"))
        .assert()
        .failure();
}

#[test]
fn tmp_dir_env_and_empty_env() {
    let d = tempdir();
    let bed = w(d.path(), "a.bed", "chr2\t1\t2\nchr1\t1\t2\n");
    // an empty GLOCI_TMPDIR is ignored
    gloci()
        .env("GLOCI_TMPDIR", "")
        .args(["sort", &bed])
        .assert()
        .success();
    // a bad GLOCI_TMPDIR is used (and fails)
    gloci()
        .env("GLOCI_TMPDIR", d.path().join("nope"))
        .args(["sort", &bed])
        .assert()
        .failure();
}

#[test]
fn in_place_rewrite() {
    let d = tempdir();
    let p = w(d.path(), "p.bed", "chr2\t1\t2\nchr1\t5\t6\nchr1\t1\t2\n");
    gloci().args(["sort", &p, "-o", &p]).assert().success();
    assert_eq!(
        std::fs::read_to_string(&p).unwrap(),
        "chr1\t1\t2\nchr1\t5\t6\nchr2\t1\t2\n"
    );
    // replacing an indexed .gz without --write-index leaves no stale index
    let gz = d.path().join("q.bed.gz");
    gloci()
        .args(["sort", &p, "--write-index", "-o"])
        .arg(&gz)
        .assert()
        .success();
    assert!(d.path().join("q.bed.gz.tbi").exists());
    gloci().args(["sort", &p, "-o"]).arg(&gz).assert().success();
    assert!(!d.path().join("q.bed.gz.tbi").exists());
}

#[test]
fn empty_inputs() {
    let d = tempdir();
    let empty_bed = w(d.path(), "e.bed", "");
    gloci()
        .args(["sort", &empty_bed])
        .assert()
        .success()
        .stdout("");
    let hv = w(d.path(), "h.vcf", &vcf(VCF_HEAD_21, VCF_COLS, &[]));
    assert!(stdout_of(&["sort", &hv]).ends_with("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"));
    let out = d.path().join("e.bed.gz");
    gloci()
        .args(["sort", &empty_bed, "--write-index", "-o"])
        .arg(&out)
        .assert()
        .success();
    assert!(d.path().join("e.bed.gz.tbi").exists());
    assert_eq!(bgzf_text(&out), "");
    gloci()
        .args(["sort", "-"])
        .write_stdin("")
        .assert()
        .success()
        .stdout("");
    // header-only BCF
    let bcf = fixtures::write(
        d.path(),
        "h.bcf",
        &fixtures::vcf_to_bcf(vcf(VCF_HEAD_21, VCF_COLS, &[]).as_bytes()),
    );
    let o = d.path().join("o.bcf");
    gloci()
        .arg("sort")
        .arg(&bcf)
        .arg("-o")
        .arg(&o)
        .assert()
        .success();
    noodles_bcf::io::Reader::new(std::fs::File::open(&o).unwrap())
        .read_header()
        .unwrap();
}

#[test]
fn many_contigs() {
    let d = tempdir();
    let n = 100_000u64;
    // a deterministic shuffle of 3 records for each of n scaffolds
    let mut ids: Vec<u64> = (0..n * 3).map(|i| i % n).collect();
    let mut x: u64 = 12345;
    for i in (1..ids.len()).rev() {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ids.swap(i, (x >> 33) as usize % (i + 1));
    }
    let mut text = String::new();
    for (i, id) in ids.iter().enumerate() {
        text.push_str(&format!("scaf{id}\t{}\t{}\n", i % 1000, i % 1000 + 5));
    }
    let p = w(d.path(), "many.bed", &text);
    let t = Instant::now();
    let out = stdout_of(&["sort", &p]);
    assert!(t.elapsed().as_secs() < 10, "took {:?}", t.elapsed());
    let names: Vec<u64> = out
        .lines()
        .map(|l| l[4..l.find('\t').unwrap()].parse().unwrap())
        .collect();
    assert_eq!(names.len() as u64, n * 3);
    assert!(names.windows(2).all(|p| p[0] <= p[1]), "natural order");
}

#[test]
fn check_on_old_lex_sorted_file_explains() {
    let d = tempdir();
    let l = w(d.path(), "l.bed", "chr1\t1\t2\nchr10\t1\t2\nchr2\t1\t2\n");
    gloci()
        .args(["sort", "--check", &l])
        .assert()
        .code(1)
        .stderr(contains(
            "looks sorted by sort -k1,1, sort-bed or bedtools sort",
        ));
    gloci()
        .args(["sort", "--check", "--chrom-order", "lex", &l])
        .assert()
        .success();
}

#[test]
fn skip_invalid_warns_at_the_end() {
    let d = tempdir();
    let p = w(d.path(), "s.bed", "chr2\t1\t2\nchr1\tx\t2\nchr1\t1\t2\n");
    gloci().args(["sort", &p]).assert().code(1);
    gloci()
        .args(["sort", &p, "--skip-invalid"])
        .assert()
        .success()
        .stdout("chr1\t1\t2\nchr2\t1\t2\n")
        .stderr(contains("warning: skipped 1 invalid records"));
}

#[test]
fn help_documents_tmp_dir() {
    let o = gloci().args(["sort", "--help"]).assert().success();
    let help = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    let flat = help.split_whitespace().collect::<Vec<_>>().join(" ");
    assert!(flat.contains("must exist and be writable"), "{help}");
}

fn bcf_in(d: &Path, name: &str, head: &str, recs: &[&str]) -> std::path::PathBuf {
    fixtures::write(
        d,
        name,
        &fixtures::vcf_to_bcf(vcf(head, VCF_COLS, recs).as_bytes()),
    )
}

#[test]
fn multi_bcf_needs_identical_dictionaries() {
    let d = tempdir();
    let swapped =
        "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n##contig=<ID=chr2,length=1000>\n";
    let a = bcf_in(
        d.path(),
        "a.bcf",
        VCF_HEAD_21,
        &["chr2\t9\t.\tA\tC\t.\t.\t."],
    );
    let b = bcf_in(d.path(), "b.bcf", swapped, &["chr1\t3\t.\tA\tC\t.\t.\t."]);
    let o = d.path().join("o.bcf");
    gloci()
        .arg("sort")
        .arg(&a)
        .arg(&b)
        .arg("-o")
        .arg(&o)
        .assert()
        .code(1)
        .stderr(
            contains("BCF headers of")
                .and(contains("a.bcf"))
                .and(contains("b.bcf"))
                .and(contains("differ"))
                .and(contains("gloci view FILE -O vcf.gz")),
        );
    assert!(!o.exists());

    // the same header in both: sorted and readable
    let c = bcf_in(
        d.path(),
        "c.bcf",
        VCF_HEAD_21,
        &["chr1\t3\t.\tA\tC\t.\t.\t."],
    );
    gloci()
        .arg("sort")
        .arg(&a)
        .arg(&c)
        .arg("-o")
        .arg(&o)
        .assert()
        .success();
    let mut r = noodles_bcf::io::Reader::new(std::fs::File::open(&o).unwrap());
    let header = r.read_header().unwrap();
    let got: Vec<(String, usize)> = r
        .records()
        .map(|rec| {
            let rec = rec.unwrap();
            (
                noodles_vcf::variant::Record::reference_sequence_name(&rec, &header)
                    .unwrap()
                    .to_string(),
                noodles_vcf::variant::Record::variant_start(&rec)
                    .unwrap()
                    .unwrap()
                    .get(),
            )
        })
        .collect();
    assert_eq!(got, [("chr2".to_string(), 9), ("chr1".to_string(), 3)]);
}

#[test]
fn two_vcfs_one_command_line_first_header() {
    let d = tempdir();
    let h1 = "##fileformat=VCFv4.3\n##source=first\n##contig=<ID=chr1,length=1000>\n";
    let h2 = "##fileformat=VCFv4.3\n##source=second\n##contig=<ID=chr1,length=1000>\n";
    let a = w(
        d.path(),
        "a.vcf",
        &vcf(h1, VCF_COLS, &["chr1\t5\t.\tA\tC\t.\t.\t."]),
    );
    let b = w(
        d.path(),
        "b.vcf",
        &vcf(h2, VCF_COLS, &["chr1\t2\t.\tA\tC\t.\t.\t."]),
    );
    let out = stdout_of(&["sort", &a, &b]);
    assert_eq!(out.matches("##glociCommand=").count(), 1);
    assert!(out.contains("##source=first\n") && !out.contains("##source=second"));
    assert!(out.find("\nchr1\t2").unwrap() < out.find("\nchr1\t5").unwrap());
}
