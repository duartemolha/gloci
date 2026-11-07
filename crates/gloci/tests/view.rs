mod common;
use common::*;
use gloci_testkit::fixtures;
use predicates::str::contains;
use std::io::Read;
use std::path::Path;

/// The `##glociCommand` line records argv as typed; under cargo test argv[0]
/// is the full path of the binary, so only the subcommand is checked.
fn has_command_line(text: &str) -> bool {
    text.lines().any(|l| {
        l.starts_with("##glociCommand=") && l.contains(" view ") && l.contains("; Version=")
    })
}

const VCF_HEAD: &str = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=1000>\n\
##contig=<ID=chr2,length=1000>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
##FILTER=<ID=q10,Description=\"Low quality\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n";

const VCF_RECS: &[&str] = &[
    "chr1\t5\trs1\tA\tC\t50\tPASS\tDP=10;AF=0.5\tGT\t0/1\t1/1",
    "chr1\t9\t.\tAG\tA,T\t.\tq10\tDP=3\tGT\t0/0\t0/2",
    "chr2\t3\t.\tG\tT\t12.5\t.\t.\tGT\t./.\t0|1",
];

fn fixture_vcf() -> String {
    format!("{VCF_HEAD}{}\n", VCF_RECS.join("\n"))
}

fn stdout_of(args: &[&str]) -> String {
    String::from_utf8(stdout_bytes(args)).unwrap()
}

fn stdout_bytes(args: &[&str]) -> Vec<u8> {
    gloci()
        .args(args)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone()
}

fn records_of(text: &str) -> Vec<&str> {
    text.lines().filter(|l| !l.starts_with('#')).collect()
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn bgzf_text(path: &Path) -> String {
    let mut out = String::new();
    noodles_bgzf::io::Reader::new(std::fs::File::open(path).unwrap())
        .read_to_string(&mut out)
        .unwrap();
    out
}

/// Records of a VCF text, re-rendered by noodles so formatting differences
/// (QUAL, float precision) cannot matter.
fn norm_vcf(text: &[u8]) -> Vec<String> {
    use noodles_vcf::variant::io::Write as _;
    let mut r = noodles_vcf::io::Reader::new(text);
    let h = r.read_header().unwrap();
    r.records()
        .map(|rec| {
            let mut w = noodles_vcf::io::Writer::new(Vec::new());
            w.write_variant_record(&h, &rec.unwrap()).unwrap();
            String::from_utf8(w.into_inner()).unwrap()
        })
        .collect()
}

fn norm_bcf(path: &Path) -> Vec<String> {
    use noodles_vcf::variant::io::Write as _;
    let mut r = noodles_bcf::io::Reader::new(std::fs::File::open(path).unwrap());
    let h = r.read_header().unwrap();
    r.records()
        .map(|rec| {
            let mut w = noodles_vcf::io::Writer::new(Vec::new());
            w.write_variant_record(&h, &rec.unwrap()).unwrap();
            String::from_utf8(w.into_inner()).unwrap()
        })
        .collect()
}

fn bcf_header_text(path: &Path) -> String {
    let mut buf = Vec::new();
    noodles_bgzf::io::Reader::new(std::fs::File::open(path).unwrap())
        .read_to_end(&mut buf)
        .unwrap();
    let l = u32::from_le_bytes(buf[5..9].try_into().unwrap()) as usize;
    String::from_utf8(buf[9..9 + l].to_vec())
        .unwrap()
        .trim_end_matches('\0')
        .to_string()
}

#[test]
fn vcf_passthrough_unchanged() {
    let d = tempdir();
    let text = fixture_vcf();
    let gz = s(&bgzf_file(d.path(), "a.vcf.gz", &text));
    let out = stdout_of(&["view", &gz]);
    assert_eq!(records_of(&out), records_of(&text)); // byte-identical record lines
    assert!(has_command_line(&out));
    assert!(out.starts_with(VCF_HEAD.split("#CHROM").next().unwrap()));
}

#[test]
fn bcf_to_vcf_and_back() {
    let d = tempdir();
    let text = fixture_vcf();
    let vcf = s(&plain_file(d.path(), "a.vcf", &text));
    let bcf = s(&fixtures::write(
        d.path(),
        "a.bcf",
        &fixtures::vcf_to_bcf(text.as_bytes()),
    ));
    let out = gloci()
        .args(["view", &bcf, "-O", "vcf"])
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    assert_eq!(norm_vcf(&out), norm_vcf(text.as_bytes()));
    let out_text = String::from_utf8(out).unwrap();
    assert!(has_command_line(&out_text));
    assert!(out_text.contains("##INFO=<ID=DP,"));

    let o = d.path().join("o.bcf");
    gloci()
        .args(["view", &vcf, "-o", &s(&o)])
        .assert()
        .success();
    assert_eq!(norm_bcf(&o), norm_vcf(text.as_bytes()));
    assert!(has_command_line(&bcf_header_text(&o)));

    // BCF to BCF keeps the records and adds the command line.
    let o2 = d.path().join("o2.bcf");
    gloci()
        .args(["view", &bcf, "-o", &s(&o2)])
        .assert()
        .success();
    assert_eq!(norm_bcf(&o2), norm_vcf(text.as_bytes()));
    assert!(has_command_line(&bcf_header_text(&o2)));

    // VCF to compressed VCF, and BCF to compressed VCF.
    let o3 = d.path().join("o3.vcf.gz");
    gloci()
        .args(["view", &bcf, "-o", &s(&o3)])
        .assert()
        .success();
    assert_eq!(
        norm_vcf(bgzf_text(&o3).as_bytes()),
        norm_vcf(text.as_bytes())
    );
}

#[test]
fn out_of_order_idx_survives_conversion() {
    // IDX= values that disagree with the line order: BCF records refer to
    // the IDX dictionaries, so the converted file must keep them.
    let head = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=1000,IDX=1>\n\
##contig=<ID=chr2,length=1000,IDX=0>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\",IDX=3>\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\",IDX=2>\n\
##FILTER=<ID=PASS,Description=\"All filters passed\",IDX=0>\n\
##FILTER=<ID=q10,Description=\"Low quality\",IDX=1>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let recs = [
        "chr1\t5\t.\tA\tC\t50\tPASS\tDP=10;AF=0.5",
        "chr1\t9\t.\tA\tT\t.\tq10\tDP=3",
        "chr2\t3\t.\tG\tT\t12.5\t.\t.",
    ];
    let d = tempdir();
    let text = format!("{head}{}\n", recs.join("\n"));
    let vcf = s(&plain_file(d.path(), "a.vcf", &text));
    let o = d.path().join("o.bcf");
    gloci()
        .args(["view", &vcf, "-o", &s(&o)])
        .assert()
        .success();
    assert_eq!(norm_bcf(&o), norm_vcf(text.as_bytes()));

    // BCF in dictionary order (chr2 first), converted back to VCF.
    let bcf_text = format!("{head}{}\n", [recs[2], recs[0], recs[1]].join("\n"));
    let bcf = s(&fixtures::write(
        d.path(),
        "b.bcf",
        &fixtures::vcf_to_bcf(bcf_text.as_bytes()),
    ));
    let out = stdout_bytes(&["view", &bcf, "-O", "vcf"]);
    assert_eq!(norm_vcf(&out), norm_vcf(bcf_text.as_bytes()));
    let b2 = d.path().join("b2.bcf");
    gloci()
        .args(["view", &bcf, "-o", &s(&b2)])
        .assert()
        .success();
    assert_eq!(norm_bcf(&b2), norm_vcf(bcf_text.as_bytes()));
}

#[test]
fn undeclared_info_is_a_clear_error() {
    let d = tempdir();
    let text = format!("{VCF_HEAD}chr1\t5\t.\tA\tC\t.\t.\tNOPE=1\tGT\t0/1\t1/1\n");
    let vcf = s(&plain_file(d.path(), "a.vcf", &text));
    let o = d.path().join("o.bcf");
    gloci()
        .args(["view", &vcf, "-o", &s(&o)])
        .assert()
        .code(1)
        .stderr(contains("line"));
    assert!(!o.exists());
}

#[test]
fn header_flags() {
    let d = tempdir();
    let vcf = s(&plain_file(d.path(), "a.vcf", &fixture_vcf()));
    let h = stdout_of(&["view", &vcf, "--header-only"]);
    assert!(!h.is_empty() && h.lines().all(|l| l.starts_with('#')));
    assert!(h.contains("#CHROM") && h.contains("##glociCommand"));
    let n = stdout_of(&["view", &vcf, "--no-header"]);
    assert!(!n.contains("##") && !n.contains("#CHROM"));
    assert_eq!(records_of(&n), VCF_RECS);
    gloci()
        .args(["view", &vcf, "--no-header", "-O", "bcf"])
        .assert()
        .code(2);
    gloci()
        .args(["view", &vcf, "--no-header", "--header-only"])
        .assert()
        .code(2);

    // BED header lines are kept by default and follow the flags.
    let bed = s(&plain_file(
        d.path(),
        "a.bed",
        "track name=x\nchr1\t1\t5\nchr1\t6\t9\n",
    ));
    assert_eq!(
        stdout_of(&["view", &bed]),
        "track name=x\nchr1\t1\t5\nchr1\t6\t9\n"
    );
    assert_eq!(
        stdout_of(&["view", &bed, "--no-header"]),
        "chr1\t1\t5\nchr1\t6\t9\n"
    );
    assert_eq!(
        stdout_of(&["view", &bed, "--header-only"]),
        "track name=x\n"
    );
}

#[test]
fn header_only_bcf_output() {
    let d = tempdir();
    let vcf = s(&plain_file(d.path(), "a.vcf", &fixture_vcf()));
    let o = d.path().join("h.bcf");
    gloci()
        .args(["view", &vcf, "--header-only", "-o", &s(&o)])
        .assert()
        .success();
    assert!(norm_bcf(&o).is_empty());
    assert!(bcf_header_text(&o).contains("#CHROM"));
}

#[test]
fn bed_recompress_with_index() {
    let d = tempdir();
    let bed = s(&plain_file(
        d.path(),
        "a.bed",
        "chr1\t1\t5\nchr1\t6\t9\nchr2\t0\t3\n",
    ));
    let o = d.path().join("o.bed.gz");
    gloci()
        .args(["view", &bed, "-o", &s(&o), "--write-index"])
        .assert()
        .success();
    assert_eq!(bgzf_text(&o), "chr1\t1\t5\nchr1\t6\t9\nchr2\t0\t3\n");
    assert!(d.path().join("o.bed.gz.tbi").exists());
}

#[test]
fn cross_family_refused() {
    let d = tempdir();
    let vcf = s(&plain_file(d.path(), "a.vcf", &fixture_vcf()));
    gloci()
        .args(["view", &vcf, "-O", "bed"])
        .assert()
        .code(2)
        .stderr(contains("gloci convert"));
    let bed = s(&plain_file(d.path(), "a.bed", "chr1\t1\t5\n"));
    gloci()
        .args(["view", &bed, "-O", "vcf"])
        .assert()
        .code(2)
        .stderr(contains("gloci convert"));
}

#[test]
fn bcf_bytes_to_pipe() {
    let d = tempdir();
    let bcf = s(&fixtures::write(
        d.path(),
        "a.bcf",
        &fixtures::vcf_to_bcf(fixture_vcf().as_bytes()),
    ));
    // not a TTY under assert_cmd
    assert!(stdout_bytes(&["view", &bcf]).starts_with(&[0x1f, 0x8b]));
    let vcf = s(&plain_file(d.path(), "a.vcf", &fixture_vcf()));
    assert!(stdout_bytes(&["view", &vcf, "-O", "bcf"]).starts_with(&[0x1f, 0x8b]));
}

#[test]
fn empty_inputs() {
    let d = tempdir();
    let empty_bed = s(&plain_file(d.path(), "e.bed", ""));
    gloci()
        .args(["view", &empty_bed])
        .assert()
        .success()
        .stdout("");
    let head_only = VCF_HEAD.to_string();
    let hb = s(&fixtures::write(
        d.path(),
        "h.bcf",
        &fixtures::vcf_to_bcf(head_only.as_bytes()),
    ));
    let out = stdout_of(&["view", &hb, "-O", "vcf"]);
    assert!(out.contains("#CHROM") && records_of(&out).is_empty());
    let hv = s(&plain_file(d.path(), "h.vcf", &head_only));
    let o = d.path().join("h2.bcf");
    gloci().args(["view", &hv, "-o", &s(&o)]).assert().success();
    assert!(norm_bcf(&o).is_empty());
    // empty stdin and an empty BGZF stream
    gloci()
        .args(["view", "-", "--format", "bed"])
        .write_stdin("")
        .assert()
        .success()
        .stdout("");
    let empty_gz = s(&fixtures::write(
        d.path(),
        "e.bed.gz",
        &fixtures::bgzip(b""),
    ));
    let o = d.path().join("e2.bed.gz");
    gloci()
        .args([
            "view",
            &empty_gz,
            "--format",
            "bed",
            "-o",
            &s(&o),
            "--write-index",
        ])
        .assert()
        .success();
    assert_eq!(bgzf_text(&o), "");
}

#[test]
fn stdin_input() {
    let out = gloci()
        .args(["view", "-"])
        .write_stdin(fixture_vcf())
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    assert_eq!(records_of(&String::from_utf8(out).unwrap()), VCF_RECS);
}

#[test]
fn unsorted_input_is_refused() {
    let d = tempdir();
    let bed = s(&plain_file(d.path(), "u.bed", "chr1\t10\t20\nchr1\t1\t5\n"));
    gloci()
        .args(["view", &bed])
        .assert()
        .code(1)
        .stderr(contains("gloci sort"));
}

#[test]
fn skip_invalid_warns() {
    let d = tempdir();
    let bed = s(&plain_file(
        d.path(),
        "b.bed",
        "chr1\t1\t5\nchr1\tx\t9\nchr1\t6\t9\n",
    ));
    gloci().args(["view", &bed]).assert().code(1);
    gloci()
        .args(["view", &bed, "--skip-invalid"])
        .assert()
        .success()
        .stdout("chr1\t1\t5\nchr1\t6\t9\n")
        .stderr(contains("skipped"));
}

#[test]
fn write_errors_name_the_output() {
    let d = tempdir();
    let bed = s(&plain_file(d.path(), "a.bed", "chr1\t1\t5\n"));
    gloci()
        .args(["view", &bed, "-o", "/nonexistent-dir/x.bed"])
        .assert()
        .code(1)
        .stderr(contains("/nonexistent-dir/x.bed"));
}

#[test]
fn interleaved_header_lines_without_idx() {
    // No IDX=, FILTER before INFO: the BCF dictionary follows line order.
    let head = "##fileformat=VCFv4.3\n\
##contig=<ID=chr1,length=1000>\n\
##FILTER=<ID=q10,Description=\"Low quality\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##FILTER=<ID=q20,Description=\"Lower quality\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let text =
        format!("{head}chr1\t5\t.\tA\tC\t50\tq20\tDP=10;AF=0.5\nchr1\t9\t.\tA\tT\t.\tq10\tDP=3\n");
    let d = tempdir();
    let vcf = s(&plain_file(d.path(), "a.vcf", &text));
    let o = d.path().join("o.bcf");
    gloci()
        .args(["view", &vcf, "-o", &s(&o)])
        .assert()
        .success();
    assert_eq!(norm_bcf(&o), norm_vcf(text.as_bytes()));
    let back = stdout_bytes(&["view", &s(&o), "-O", "vcf"]);
    assert_eq!(norm_vcf(&back), norm_vcf(text.as_bytes()));
}
