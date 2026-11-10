//! Tests for `gloci vcf concat` (D74), checked mode and `--naive`.
mod common;
use common::*;
use gloci_testkit::fixtures;
use predicates::prelude::*;
use predicates::str::contains;
use std::io::Read;
use std::path::Path;

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

const COLS: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n";

fn vcf(extra_header: &str, recs: &[&str]) -> String {
    let mut t = format!(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n##contig=<ID=chr2,length=1000>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n{extra_header}{COLS}"
    );
    for r in recs {
        t.push_str(r);
        t.push('\n');
    }
    t
}

fn rec(chrom: &str, pos: u32) -> String {
    format!("{chrom}\t{pos}\t.\tA\tC\t50\tPASS\tDP={pos}\tGT\t0/1\t1/1")
}

fn text(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

/// `view -O vcf` output without the `##glociCommand` line `view` adds.
fn view_text(path: &str) -> String {
    text(&["view", path, "-O", "vcf"])
        .lines()
        .filter(|l| !l.starts_with("##glociCommand"))
        .map(|l| format!("{l}\n"))
        .collect()
}

fn gunzip(path: &Path) -> Vec<u8> {
    let mut out = Vec::new();
    noodles_bgzf::io::Reader::new(std::fs::File::open(path).unwrap())
        .read_to_end(&mut out)
        .unwrap();
    out
}

#[test]
fn concatenates_two_per_chromosome_files() {
    let d = tempdir();
    let a = s(&plain_file(
        d.path(),
        "a.vcf",
        &vcf("", &[&rec("chr1", 5), &rec("chr1", 9)]),
    ));
    let b = s(&plain_file(
        d.path(),
        "b.vcf",
        &vcf("", &[&rec("chr2", 3), &rec("chr2", 4)]),
    ));
    let out = text(&["vcf", "concat", &a, &b]);
    let want = vcf(
        "",
        &[
            &rec("chr1", 5),
            &rec("chr1", 9),
            &rec("chr2", 3),
            &rec("chr2", 4),
        ],
    );
    assert_eq!(out, want);
    // equal positions at the seam are fine
    let c = s(&plain_file(
        d.path(),
        "c.vcf",
        &vcf("", &[&rec("chr1", 9), &rec("chr2", 1)]),
    ));
    gloci().args(["vcf", "concat", &a, &c]).assert().success();
}

#[test]
fn threads_do_not_change_the_output() {
    let d = tempdir();
    let a = s(&bgzf_file(
        d.path(),
        "a.vcf.gz",
        &vcf("", &[&rec("chr1", 5)]),
    ));
    let b = s(&bgzf_file(
        d.path(),
        "b.vcf.gz",
        &vcf("", &[&rec("chr2", 3)]),
    ));
    for naive in [false, true] {
        let mut outs = Vec::new();
        for t in ["1", "4"] {
            let o = d.path().join(format!("o{naive}{t}.vcf.gz"));
            let mut args = vec!["vcf", "concat", &a, &b, "-t", t, "-o"];
            let os = s(&o);
            args.push(&os);
            if naive {
                args.push("--naive");
            }
            gloci().args(&args).assert().success();
            outs.push(std::fs::read(&o).unwrap());
        }
        assert_eq!(outs[0], outs[1]);
    }
}

#[test]
fn sample_mismatch_exits_1_naming_the_difference() {
    let d = tempdir();
    let a = s(&plain_file(d.path(), "a.vcf", &vcf("", &[&rec("chr1", 5)])));
    let b_text = vcf("", &[&rec("chr2", 3)]).replace("s2", "s3");
    let b = s(&plain_file(d.path(), "b.vcf", &b_text));
    gloci()
        .args(["vcf", "concat", &a, &b])
        .assert()
        .code(1)
        .stderr(contains(format!(
            "{b} has samples [s1,s3] but {a} has [s1,s2]"
        )));
}

#[test]
fn header_definition_mismatch_exits_1_with_the_multi_vcf_message() {
    let d = tempdir();
    let extra = "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Freq\">\n";
    let a = s(&plain_file(d.path(), "a.vcf", &vcf("", &[&rec("chr1", 5)])));
    let b = s(&plain_file(
        d.path(),
        "b.vcf",
        &vcf(extra, &[&rec("chr2", 3)]),
    ));
    gloci()
        .args(["vcf", "concat", &a, &b])
        .assert()
        .code(1)
        .stderr(contains(format!(
            "{b} defines ##INFO=<ID=AF> but {a} does not; give the files the same header definitions"
        )));
    let (ga, gb) = (
        s(&bgzf_file(
            d.path(),
            "a.vcf.gz",
            &vcf("", &[&rec("chr1", 5)]),
        )),
        s(&bgzf_file(
            d.path(),
            "b.vcf.gz",
            &vcf(extra, &[&rec("chr2", 3)]),
        )),
    );
    gloci()
        .args(["vcf", "concat", "--naive", &ga, &gb])
        .assert()
        .code(1)
        .stderr(contains("defines ##INFO=<ID=AF>"));
}

#[test]
fn order_across_files_exits_1_with_a_sort_hint() {
    let d = tempdir();
    let a = s(&plain_file(d.path(), "a.vcf", &vcf("", &[&rec("chr2", 5)])));
    let b = s(&plain_file(d.path(), "b.vcf", &vcf("", &[&rec("chr1", 3)])));
    gloci()
        .args(["vcf", "concat", &a, &b])
        .assert()
        .code(1)
        .stderr(
            contains(format!("{b} is not sorted after {a}: line 7 ("))
                .and(contains("comes before the last record of"))
                .and(contains(format!("hint: run: gloci sort {a} {b}"))),
        );
    // within the same chromosome too
    let c = s(&plain_file(d.path(), "c.vcf", &vcf("", &[&rec("chr2", 4)])));
    gloci().args(["vcf", "concat", &a, &c]).assert().code(1);
    // an unsorted file inside is still caught by checked mode
    let u = s(&plain_file(
        d.path(),
        "u.vcf",
        &vcf("", &[&rec("chr1", 9), &rec("chr1", 2)]),
    ));
    gloci()
        .args(["vcf", "concat", &u])
        .assert()
        .code(1)
        .stderr(contains("is not sorted"));
}

#[test]
fn order_across_files_is_checked_by_naive_too() {
    let d = tempdir();
    let a = s(&bgzf_file(
        d.path(),
        "a.vcf.gz",
        &vcf("", &[&rec("chr2", 5)]),
    ));
    let b = s(&bgzf_file(
        d.path(),
        "b.vcf.gz",
        &vcf("", &[&rec("chr1", 3)]),
    ));
    gloci()
        .args(["vcf", "concat", "--naive", &a, &b])
        .assert()
        .code(1)
        .stderr(
            contains(format!("{b} is not sorted after {a}: line 7 ("))
                .and(contains(format!("hint: run: gloci sort {a} {b}"))),
        );
}

#[test]
fn bcf_inputs_and_mixtures() {
    let d = tempdir();
    let at = vcf("", &[&rec("chr1", 5), &rec("chr1", 9)]);
    let bt = vcf("", &[&rec("chr2", 3)]);
    let a = s(&fixtures::write(
        d.path(),
        "a.bcf",
        &fixtures::vcf_to_bcf(at.as_bytes()),
    ));
    let b = s(&fixtures::write(
        d.path(),
        "b.bcf",
        &fixtures::vcf_to_bcf(bt.as_bytes()),
    ));
    let bv = s(&bgzf_file(d.path(), "b.vcf.gz", &bt));
    let av = s(&plain_file(d.path(), "a.vcf", &at));
    let want = vcf("", &[&rec("chr1", 5), &rec("chr1", 9), &rec("chr2", 3)]);

    // BCF in, BCF out; read it back as VCF
    let o = s(&d.path().join("o.bcf"));
    gloci()
        .args(["vcf", "concat", &a, &b, "-o", &o])
        .assert()
        .success();
    assert_eq!(view_text(&o), want);
    // BCF in, VCF out (the format follows -O)
    assert_eq!(text(&["vcf", "concat", &a, &b, "-O", "vcf"]), want);
    // VCF then BCF, to VCF and to BCF
    assert_eq!(text(&["vcf", "concat", &av, &b, "-O", "vcf"]), want);
    let o2 = s(&d.path().join("o2.bcf"));
    gloci()
        .args(["vcf", "concat", &av, &b, "-o", &o2])
        .assert()
        .success();
    assert_eq!(view_text(&o2), want);
    // BCF then VCF.gz
    let o3 = s(&d.path().join("o3.bcf"));
    gloci()
        .args(["vcf", "concat", &a, &bv, "-o", &o3])
        .assert()
        .success();
    assert_eq!(view_text(&o3), want);
    // the order across BCF files is checked
    gloci()
        .args(["vcf", "concat", &b, &a])
        .assert()
        .code(1)
        .stderr(
            contains(format!("{a} is not sorted after {b}: record 1 (chr1 5)"))
                .and(contains("comes before the last record of")),
        );
}

#[test]
fn naive_refuses_bcf_plain_and_stdin_with_exit_2() {
    let d = tempdir();
    let t = vcf("", &[&rec("chr1", 5)]);
    let b = s(&fixtures::write(
        d.path(),
        "a.bcf",
        &fixtures::vcf_to_bcf(t.as_bytes()),
    ));
    let p = s(&plain_file(d.path(), "p.vcf", &t));
    gloci()
        .args(["vcf", "concat", "--naive", &b])
        .assert()
        .code(2)
        .stderr(contains("BCF").and(contains("leave out --naive")));
    gloci()
        .args(["vcf", "concat", "--naive", &p])
        .assert()
        .code(2)
        .stderr(contains("must be BGZF-compressed VCF"));
    gloci()
        .args(["vcf", "concat", "--naive", "-"])
        .assert()
        .code(2);
    let g = s(&bgzf_file(d.path(), "g.vcf.gz", &t));
    gloci()
        .args(["vcf", "concat", "--naive", &g, "-O", "vcf"])
        .assert()
        .code(2)
        .stderr(contains("BGZF-compressed VCF"));
}

#[test]
fn naive_with_different_headers_exits_2_with_a_hint() {
    let d = tempdir();
    let a = s(&bgzf_file(
        d.path(),
        "a.vcf.gz",
        &vcf("", &[&rec("chr1", 5)]),
    ));
    // a meta line only: compatible for checked mode, not byte-identical
    let b = s(&bgzf_file(
        d.path(),
        "b.vcf.gz",
        &vcf("##source=other\n", &[&rec("chr2", 3)]),
    ));
    gloci().args(["vcf", "concat", &a, &b]).assert().success();
    gloci()
        .args(["vcf", "concat", "--naive", &a, &b])
        .assert()
        .code(2)
        .stderr(
            contains("byte-identical headers")
                .and(contains("##source=other"))
                .and(contains("leave out --naive")),
        );
}

/// A header of many INFO lines, so it spans several BGZF blocks and ends
/// inside the last one.
fn big_header_vcf(recs: &[String]) -> String {
    let mut extra = String::new();
    for i in 0..2500 {
        extra.push_str(&format!(
            "##INFO=<ID=X{i},Number=1,Type=Integer,Description=\"Filler line number {i} to make the header long\">\n"
        ));
    }
    let r: Vec<&str> = recs.iter().map(String::as_str).collect();
    vcf(&extra, &r)
}

#[test]
fn naive_concat_block_edges() {
    let d = tempdir();
    let many = |chrom: &str| -> Vec<String> { (1..=3000).map(|p| rec(chrom, p * 7)).collect() };
    let (at, bt, ct) = (
        big_header_vcf(&many("chr1")),
        big_header_vcf(&many("chr2")),
        big_header_vcf(&[]),
    );
    assert!(at.len() > 3 * 65280, "the header must span several blocks");
    let a = s(&bgzf_file(d.path(), "a.vcf.gz", &at));
    let b = s(&bgzf_file(d.path(), "b.vcf.gz", &bt));
    let c = s(&bgzf_file(d.path(), "c.vcf.gz", &ct));
    let n = s(&d.path().join("naive.vcf.gz"));
    let k = s(&d.path().join("checked.vcf.gz"));
    // with a header-only file in the middle
    gloci()
        .args(["vcf", "concat", "--naive", &a, &c, &b, "-o", &n])
        .assert()
        .success();
    gloci()
        .args(["vcf", "concat", &a, &c, &b, "-o", &k])
        .assert()
        .success();
    let (nd, kd) = (gunzip(Path::new(&n)), gunzip(Path::new(&k)));
    assert_eq!(
        nd, kd,
        "naive and checked output differ after decompression"
    );
    // the output ends with exactly one EOF marker and has none in the middle
    let raw = std::fs::read(&n).unwrap();
    let eof = &raw[raw.len() - 28..];
    assert_eq!(eof, gloci_io_eof());
    let markers = raw.windows(28).filter(|w| *w == gloci_io_eof()).count();
    assert_eq!(markers, 1);
    // gloci reads it back, with every record once
    let back = view_text(&n);
    assert_eq!(back.lines().filter(|l| !l.starts_with('#')).count(), 6000);
    assert_eq!(back.as_bytes(), nd.as_slice());
    // the first file's blocks are copied as they are
    let ra = std::fs::read(&a).unwrap();
    assert!(raw.starts_with(&ra[..ra.len() - 28]));
}

fn gloci_io_eof() -> &'static [u8] {
    &[
        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0, 0x1b, 0, 0x03,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
    ]
}

#[test]
fn naive_accepts_an_input_without_an_eof_marker_with_a_warning() {
    let d = tempdir();
    let at = vcf("", &[&rec("chr1", 5)]);
    let bt = vcf("", &[&rec("chr2", 3)]);
    let mut raw = fixtures::bgzip(at.as_bytes());
    raw.truncate(raw.len() - 28);
    let a = s(&fixtures::write(d.path(), "a.vcf.gz", &raw));
    let b = s(&bgzf_file(d.path(), "b.vcf.gz", &bt));
    let o = s(&d.path().join("o.vcf.gz"));
    gloci()
        .args(["vcf", "concat", "--naive", &a, &b, "-o", &o])
        .assert()
        .success()
        .stderr(contains(format!("{a}: BGZF end-of-file marker missing")));
    assert_eq!(
        String::from_utf8(gunzip(Path::new(&o))).unwrap(),
        vcf("", &[&rec("chr1", 5), &rec("chr2", 3)])
    );
    // checked mode still refuses a truncated file
    gloci()
        .args(["vcf", "concat", &a, &b])
        .assert()
        .failure()
        .stderr(contains("end-of-file marker missing"));
}

#[test]
fn naive_output_to_stdout_and_a_missing_final_newline() {
    let d = tempdir();
    let a = s(&bgzf_file(
        d.path(),
        "a.vcf.gz",
        &vcf("", &[&rec("chr1", 5)]),
    ));
    let b = s(&bgzf_file(
        d.path(),
        "b.vcf.gz",
        &vcf("", &[&rec("chr2", 3)]),
    ));
    let o = gloci()
        .args(["vcf", "concat", "--naive", &a, &b, "-O", "vcf.gz"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let mut got = Vec::new();
    noodles_bgzf::io::Reader::new(o.as_slice())
        .read_to_end(&mut got)
        .unwrap();
    assert_eq!(
        String::from_utf8(got).unwrap(),
        vcf("", &[&rec("chr1", 5), &rec("chr2", 3)])
    );
    // a file whose last record has no newline cannot be followed by another
    let nonl = vcf("", &[&rec("chr1", 5)]);
    let c = s(&bgzf_file(d.path(), "c.vcf.gz", nonl.trim_end()));
    gloci()
        .args(["vcf", "concat", "--naive", &c, &b])
        .assert()
        .code(1)
        .stderr(contains("no final newline"));
    // checked mode handles it
    gloci().args(["vcf", "concat", &c, &b]).assert().success();
}
