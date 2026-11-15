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
fn bed_default_is_lex() {
    let d = tempdir();
    let u = w(
        d.path(),
        "u.bed",
        "chr2\t5\t10\nchrX\t1\t2\nchr10\t1\t2\nchr1\t7\t9\nchrM\t0\t5\nchr1\t3\t4\n",
    );
    let lex = "chr1\t3\t4\nchr1\t7\t9\nchr10\t1\t2\nchr2\t5\t10\nchrM\t0\t5\nchrX\t1\t2\n";
    gloci().args(["sort", &u]).assert().success().stdout(lex);
    gloci()
        .args(["sort", &u, "--chrom-order", "lex"])
        .assert()
        .success()
        .stdout(lex);
    gloci()
        .args(["sort", &u, "--chrom-order", "karyotype"])
        .assert()
        .success()
        .stdout("chr1\t3\t4\nchr1\t7\t9\nchr2\t5\t10\nchr10\t1\t2\nchrX\t1\t2\nchrM\t0\t5\n");
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
    let b = w(d.path(), "b.bed", "chr1\t5\t6\nchr10\t1\t2\nchr2\t0\t1\n");
    gloci()
        .args(["sort", "-m", &a, &b])
        .assert()
        .success()
        .stdout("chr1\t1\t2\nchr1\t5\t6\nchr1\t9\t10\nchr10\t1\t2\nchr2\t0\t1\nchr2\t1\t2\n");
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

#[cfg(unix)]
#[test]
fn sort_300_inputs() {
    use std::os::unix::process::CommandExt;
    let d = tempdir();
    let mut files = Vec::new();
    let mut all_records = Vec::new();

    for i in 0..300 {
        let name = format!("in_{:03}.bed", i);
        let start1 = (300 - i) * 10;
        let start2 = start1 + 5;
        let content = format!(
            "chr1\t{}\t{}\nchr1\t{}\t{}\n",
            start1,
            start1 + 1,
            start2,
            start2 + 1
        );
        let path = w(d.path(), &name, &content);
        files.push(path);
        all_records.push((start1, start1 + 1));
        all_records.push((start2, start2 + 1));
    }

    all_records.sort();
    let mut want = String::new();
    for (s, e) in all_records {
        want.push_str(&format!("chr1\t{s}\t{e}\n"));
    }

    let run_cmd = |extra_args: &[&str]| -> String {
        let bin = assert_cmd::cargo::cargo_bin("gloci");
        let mut cmd = std::process::Command::new(bin);
        cmd.arg("sort");
        cmd.args(extra_args);
        cmd.args(&files);
        unsafe {
            cmd.pre_exec(|| {
                let rlim = libc::rlimit {
                    rlim_cur: 200,
                    rlim_max: 200,
                };
                if libc::setrlimit(libc::RLIMIT_NOFILE, &rlim) != 0 {
                    return Err(std::io::Error::last_os_error());
                }
                Ok(())
            });
        }
        let out = cmd.output().expect("failed to execute gloci");
        assert!(
            out.status.success(),
            "gloci exited with {}: {}",
            out.status,
            String::from_utf8_lossy(&out.stderr)
        );
        String::from_utf8(out.stdout).unwrap()
    };

    // Test 1: normal sort
    let got_sort = run_cmd(&[]);
    assert_eq!(
        got_sort, want,
        "gloci sort 300 inputs failed under ulimit -n 200"
    );

    // Test 2: sort -m (merge_sorted)
    let got_merge = run_cmd(&["-m"]);
    assert_eq!(
        got_merge, want,
        "gloci sort -m 300 inputs failed under ulimit -n 200"
    );
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
    let out = stdout_of(&["sort", &p, "--chrom-order", "karyotype"]);
    // 3 s alone, 12 s on a loaded 4-core box. A scan of the contig table per record would be
    // 3e10 comparisons, minutes; the bound only has to tell the two apart.
    assert!(t.elapsed().as_secs() < 60, "took {:?}", t.elapsed());
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
    gloci().args(["sort", "--check", &l]).assert().success();
    gloci()
        .args(["sort", "--check", "--chrom-order", "karyotype", &l])
        .assert()
        .code(1)
        .stderr(contains(
            "looks sorted by sort -k1,1, sort-bed or bedtools sort",
        ));
}

#[test]
fn check_on_karyotype_sorted_file_explains() {
    let d = tempdir();
    let k = w(d.path(), "k.bed", "chr1\t1\t2\nchr2\t1\t2\nchr10\t1\t2\n");
    gloci()
        .args(["sort", "--check", &k])
        .assert()
        .code(1)
        .stderr(contains("looks sorted in karyotype order"))
        .stderr(contains("-g GENOME or --chrom-order karyotype"));
    gloci()
        .args(["sort", "--check", "--chrom-order", "karyotype", &k])
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

// ---- C1: every input is read exactly once (stdin, pipes, FIFOs) ----

#[cfg(unix)]
#[test]
fn stdin_plus_file_keeps_every_record() {
    let d = tempdir();
    let b = w(d.path(), "b.bed", "chr1\t5\t6\nchr1\t1\t2\n");
    let out = bash_gloci(&format!("sort - '{b}'"), "chr2\t1\t2\nchr1\t3\t4\n");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(
        String::from_utf8(out.stdout).unwrap(),
        "chr1\t1\t2\nchr1\t3\t4\nchr1\t5\t6\nchr2\t1\t2\n"
    );
    // File first, stdin second.
    let out = bash_gloci(&format!("sort '{b}' -"), "chr2\t1\t2\nchr1\t3\t4\n");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(
        String::from_utf8(out.stdout).unwrap(),
        "chr1\t1\t2\nchr1\t3\t4\nchr1\t5\t6\nchr2\t1\t2\n"
    );
}

#[cfg(unix)]
#[test]
fn process_substitutions_keep_every_record() {
    let out = bash_gloci(
        "sort <(printf 'chr2\\t1\\t2\\nchr1\\t9\\t10\\n') <(printf 'chr1\\t1\\t2\\n')",
        "",
    );
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(
        String::from_utf8(out.stdout).unwrap(),
        "chr1\t1\t2\nchr1\t9\t10\nchr2\t1\t2\n"
    );
    let out = bash_gloci(
        "sort -m <(printf 'chr1\\t1\\t2\\nchr1\\t5\\t6\\n') <(printf 'chr1\\t3\\t4\\n')",
        "",
    );
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(
        String::from_utf8(out.stdout).unwrap(),
        "chr1\t1\t2\nchr1\t3\t4\nchr1\t5\t6\n"
    );
}

#[cfg(unix)]
#[test]
fn check_reads_piped_inputs() {
    let out = bash_gloci(
        "sort --check <(printf 'chr1\\t5\\t6\\nchr1\\t1\\t2\\n') <(printf 'chr1\\t3\\t4\\n')",
        "",
    );
    assert_eq!(
        out.status.code(),
        Some(1),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(
        String::from_utf8_lossy(&out.stderr).contains("is not sorted"),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    // Sorted pipes pass.
    let out = bash_gloci(
        "sort --check <(printf 'chr1\\t1\\t2\\n') <(printf 'chr1\\t3\\t4\\n')",
        "",
    );
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
}

/// A stdin of many batches (and more than the 64 KiB detection peek).
fn many_lines(n: usize, chrom: &str) -> String {
    (0..n)
        .map(|i| format!("{chrom}\t{}\t{}\n", n - i, n - i + 1))
        .collect()
}

#[cfg(unix)]
#[test]
fn large_stdin_plus_file() {
    let d = tempdir();
    let b = w(d.path(), "b.bed", "chr1\t0\t1\n");
    let input = many_lines(30_000, "chr1");
    let out = bash_gloci(&format!("sort - '{b}'"), &input);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let got = String::from_utf8(out.stdout).unwrap();
    assert_eq!(got.lines().count(), 30_001);
    assert!(got.starts_with("chr1\t0\t1\nchr1\t1\t2\nchr1\t2\t3\n"));
}

/// More inputs than one merge round holds open (`--max-mem 16M` gives a
/// fan-in of 42): stdin is spooled once and every record survives in
/// `sort`, `sort -m` and `--check`.
#[cfg(unix)]
#[test]
fn stdin_among_more_inputs_than_the_fan_in() {
    let d = tempdir();
    let mut files = Vec::new();
    let mut want = Vec::new();
    for i in 0..60u64 {
        let s = 1000 + i * 10;
        files.push(format!(
            "'{}'",
            w(
                d.path(),
                &format!("f{i:02}.bed"),
                &format!("chr1\t{s}\t{}\n", s + 1)
            )
        ));
        want.push((s, s + 1));
    }
    // stdin: sorted (for -m), many batches, interleaved with the files.
    let mut input = String::new();
    for i in 0..20_000u64 {
        input.push_str(&format!("chr1\t{}\t{}\n", i * 2, i * 2 + 1));
        want.push((i * 2, i * 2 + 1));
    }
    want.sort();
    let want: String = want
        .iter()
        .map(|(s, e)| format!("chr1\t{s}\t{e}\n"))
        .collect();
    let mut args = files.clone();
    args.insert(30, "-".into());
    let list = args.join(" ");
    let tmp = d.path().join("tmp");
    std::fs::create_dir(&tmp).unwrap();
    let tmp = tmp.display().to_string();
    for mode in ["", "-m"] {
        let out = bash_gloci(
            &format!("sort {mode} --max-mem 16M --tmp-dir '{tmp}' {list}"),
            &input,
        );
        assert!(
            out.status.success(),
            "{mode}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        assert!(
            String::from_utf8(out.stdout).unwrap() == want,
            "{mode}: output differs"
        );
    }
    let out = bash_gloci(&format!("sort --check --max-mem 16M {list}"), &input);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    // An unsorted stdin among them fails the check.
    let out = bash_gloci(
        &format!("sort --check --max-mem 16M {list}"),
        "chr1\t9\t10\nchr1\t1\t2\n",
    );
    assert_eq!(
        out.status.code(),
        Some(1),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    // No spool left behind.
    assert_eq!(std::fs::read_dir(&tmp).unwrap().count(), 0);
}

/// Sort `n` BGZF inputs of `records` records each at `-t 8`; returns the
/// most threads the process had at once (polled from `/proc/PID/status`
/// while it ran) and its peak RSS in bytes.
#[cfg(target_os = "linux")]
#[allow(clippy::zombie_processes)] // reaped by wait4, for its rusage
fn sort_bgzf_inputs_peaks(n: usize, records: usize) -> (u64, u64) {
    use std::sync::Arc;
    use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
    let d = tempdir();
    let files: Vec<String> = (0..n)
        .map(|i| {
            let mut text = String::new();
            for r in 0..records {
                // Interleaved across files so every input contributes to every part.
                let start = (r * n + (n - i)) * 10;
                text.push_str(&format!("chr1\t{start}\t{}\tf{i}_{r}\n", start + 5));
            }
            bgzf_file(d.path(), &format!("in_{i:03}.bed.gz"), &text)
                .to_string_lossy()
                .into_owned()
        })
        .collect();
    let mut child = std::process::Command::new(assert_cmd::cargo::cargo_bin("gloci"))
        .args(["sort", "-t", "8", "-o", "/dev/null"])
        .args(&files)
        .stderr(std::process::Stdio::piped())
        .spawn()
        .unwrap();
    let pid = child.id();
    let peak = Arc::new(AtomicU64::new(0));
    let stop = Arc::new(AtomicBool::new(false));
    let poller = {
        let (peak, stop) = (peak.clone(), stop.clone());
        std::thread::spawn(move || {
            while !stop.load(Ordering::Relaxed) {
                if let Ok(s) = std::fs::read_to_string(format!("/proc/{pid}/status"))
                    && let Some(l) = s.lines().find_map(|l| l.strip_prefix("Threads:"))
                    && let Ok(t) = l.trim().parse::<u64>()
                {
                    peak.fetch_max(t, Ordering::Relaxed);
                }
                std::thread::sleep(std::time::Duration::from_micros(200));
            }
        })
    };
    let mut status = 0;
    let mut ru: libc::rusage = unsafe { std::mem::zeroed() };
    let r = unsafe { libc::wait4(pid as libc::pid_t, &mut status, 0, &mut ru) };
    stop.store(true, Ordering::Relaxed);
    poller.join().unwrap();
    assert_eq!(r, pid as libc::pid_t, "wait4 failed");
    let mut err = String::new();
    child
        .stderr
        .take()
        .unwrap()
        .read_to_string(&mut err)
        .unwrap();
    assert!(
        libc::WIFEXITED(status) && libc::WEXITSTATUS(status) == 0,
        "sort failed: {err}"
    );
    (peak.load(Ordering::Relaxed), ru.ru_maxrss as u64 * 1024)
}

/// Regular files are reopened per pass, so 100 BGZF inputs
/// do not hold 100 readers (each with its own thread pool) open at once.
#[cfg(target_os = "linux")]
#[test]
fn sort_many_bgzf_inputs_bounded() {
    let (threads, _) = sort_bgzf_inputs_peaks(100, 50);
    assert!(threads < 120, "peak {threads} threads for 100 inputs");
}

/// The same with real data volume: peak RSS stays under 150 MB.
#[cfg(target_os = "linux")]
#[test]
#[ignore = "heavy: 100 inputs of 2k records; run with --ignored"]
fn sort_many_bgzf_inputs_bounded_rss() {
    let (threads, rss) = sort_bgzf_inputs_peaks(100, 2_000);
    assert!(threads < 120, "peak {threads} threads");
    assert!(rss < 150_000_000, "peak RSS {rss} bytes");
}
