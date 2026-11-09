//! `gloci intersect`: modes, thresholds, strand, windows, VCF input,
//! `--match alleles`, empty inputs, order errors and determinism.
mod common;
use common::*;
use predicates::str::contains;
use std::path::Path;

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

/// Run gloci with `args` and return stdout, asserting success.
fn out(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

fn code(args: &[&str]) -> i32 {
    gloci().args(args).output().unwrap().status.code().unwrap()
}

const A: &str = "chr1\t10\t20\ta1\t0\t+\n\
chr1\t30\t40\ta2\t0\t-\n\
chr1\t50\t60\ta3\t0\t+\n\
chr2\t5\t15\ta4\t0\t+\n";
const B: &str = "chr1\t15\t35\tb1\t0\t+\n\
chr1\t18\t19\tb2\t0\t-\n\
chr2\t0\t8\tb3\t0\t+\n";

struct Fx {
    dir: tempfile::TempDir,
}

impl Fx {
    fn new() -> Fx {
        Fx { dir: tempdir() }
    }
    fn file(&self, name: &str, text: &str) -> String {
        s(&plain_file(self.dir.path(), name, text))
    }
}

#[test]
fn worked_examples_table() {
    // spec 4.3: (A, B, overlaps, overlaps with --inclusive-points)
    let rows: [(&str, &str, bool, bool); 7] = [
        ("chr1\t10\t10", "chr1\t5\t15", true, true),
        ("chr1\t5\t5", "chr1\t5\t15", false, true),
        ("chr1\t15\t15", "chr1\t5\t15", false, true),
        ("chr1\t10\t10", "chr1\t10\t10", true, true),
        ("chr1\t10\t10", "chr1\t11\t11", false, false),
        ("chr1\t5\t10", "chr1\t10\t15", false, false),
        ("chr1\t5\t11", "chr1\t10\t15", true, true),
    ];
    let fx = Fx::new();
    for (i, (a, b, plain, inclusive)) in rows.iter().enumerate() {
        let af = fx.file(&format!("a{i}.bed"), &format!("{a}\n"));
        let bf = fx.file(&format!("b{i}.bed"), &format!("{b}\n"));
        let want = |yes: bool| if yes { format!("{a}\n") } else { String::new() };
        let base = ["intersect", "-a", &af, "-b", &bf, "--mode", "any"];
        assert_eq!(out(&base), want(*plain), "row {i}: {a} vs {b}");
        let mut v = base.to_vec();
        v.push("--inclusive-points");
        assert_eq!(out(&v), want(*inclusive), "row {i} inclusive: {a} vs {b}");
    }
}

#[test]
fn modes() {
    let fx = Fx::new();
    let (a, b) = (fx.file("a.bed", A), fx.file("b.bed", B));
    let run = |extra: &[&str]| {
        let mut v = vec!["intersect", "-a", &a, "-b", &b];
        v.extend_from_slice(extra);
        out(&v)
    };
    let a_lines: Vec<&str> = A.lines().collect();
    let b_lines: Vec<&str> = B.lines().collect();
    // default is any
    let any = format!("{}\n{}\n{}\n", a_lines[0], a_lines[1], a_lines[3]);
    assert_eq!(run(&[]), any);
    assert_eq!(run(&["--mode", "any"]), any);
    assert_eq!(run(&["--mode", "none"]), format!("{}\n", a_lines[2]));
    let pairs = format!(
        "{a0}\t{b0}\n{a0}\t{b1}\n{a1}\t{b0}\n{a3}\t{b2}\n",
        a0 = a_lines[0],
        a1 = a_lines[1],
        a3 = a_lines[3],
        b0 = b_lines[0],
        b1 = b_lines[1],
        b2 = b_lines[2]
    );
    assert_eq!(run(&["--mode", "pairs"]), pairs);
    assert_eq!(
        run(&["--mode", "pairs", "--overlap-bases"]),
        "chr1\t10\t20\ta1\t0\t+\tchr1\t15\t35\tb1\t0\t+\t5\n\
chr1\t10\t20\ta1\t0\t+\tchr1\t18\t19\tb2\t0\t-\t1\n\
chr1\t30\t40\ta2\t0\t-\tchr1\t15\t35\tb1\t0\t+\t5\n\
chr2\t5\t15\ta4\t0\t+\tchr2\t0\t8\tb3\t0\t+\t3\n"
    );
    assert_eq!(
        run(&["--mode", "pairs", "--keep-unmatched"]),
        "chr1\t10\t20\ta1\t0\t+\tchr1\t15\t35\tb1\t0\t+\n\
chr1\t10\t20\ta1\t0\t+\tchr1\t18\t19\tb2\t0\t-\n\
chr1\t30\t40\ta2\t0\t-\tchr1\t15\t35\tb1\t0\t+\n\
chr1\t50\t60\ta3\t0\t+\t.\t-1\t-1\t.\t-1\t.\n\
chr2\t5\t15\ta4\t0\t+\tchr2\t0\t8\tb3\t0\t+\n"
    );
    assert_eq!(
        run(&["--mode", "pairs", "--keep-unmatched", "--overlap-bases"]),
        "chr1\t10\t20\ta1\t0\t+\tchr1\t15\t35\tb1\t0\t+\t5\n\
chr1\t10\t20\ta1\t0\t+\tchr1\t18\t19\tb2\t0\t-\t1\n\
chr1\t30\t40\ta2\t0\t-\tchr1\t15\t35\tb1\t0\t+\t5\n\
chr1\t50\t60\ta3\t0\t+\t.\t-1\t-1\t.\t-1\t.\t0\n\
chr2\t5\t15\ta4\t0\t+\tchr2\t0\t8\tb3\t0\t+\t3\n"
    );
    assert_eq!(
        run(&["--mode", "clip"]),
        "chr1\t15\t20\ta1\t0\t+\n\
chr1\t18\t19\ta1\t0\t+\n\
chr1\t30\t35\ta2\t0\t-\n\
chr2\t5\t8\ta4\t0\t+\n"
    );
    assert_eq!(
        run(&["--mode", "count"]),
        "chr1\t10\t20\ta1\t0\t+\t2\n\
chr1\t30\t40\ta2\t0\t-\t1\n\
chr1\t50\t60\ta3\t0\t+\t0\n\
chr2\t5\t15\ta4\t0\t+\t1\n"
    );
}

#[test]
fn flags_need_pairs_mode() {
    let fx = Fx::new();
    let (a, b) = (fx.file("a.bed", A), fx.file("b.bed", B));
    for extra in [
        &["--overlap-bases"][..],
        &["--keep-unmatched"],
        &["--mode", "count", "--overlap-bases"],
    ] {
        let mut v = vec!["intersect", "-a", &a, "-b", &b];
        v.extend_from_slice(extra);
        assert_eq!(code(&v), 2, "{extra:?}");
    }
}

#[test]
fn multiple_b_with_names() {
    let fx = Fx::new();
    let a = fx.file("a.bed", "chr1\t10\t20\ta1\nchr1\t50\t60\ta3\n");
    let x = fx.file("x.bed", "chr1\t15\t25\tx1\n");
    let y = fx.file("y.bed", "chr1\t5\t12\ty1\nchr1\t19\t30\ty2\n");
    let run = |extra: &[&str]| {
        let mut v = vec!["intersect", "-a", &a, "-b", &x, &y];
        v.extend_from_slice(extra);
        out(&v)
    };
    assert_eq!(
        run(&["--mode", "pairs", "--names", "X,Y"]),
        "chr1\t10\t20\ta1\tX\tchr1\t15\t25\tx1\n\
chr1\t10\t20\ta1\tY\tchr1\t5\t12\ty1\n\
chr1\t10\t20\ta1\tY\tchr1\t19\t30\ty2\n"
    );
    // labels default to basenames
    assert!(run(&["--mode", "pairs"]).starts_with("chr1\t10\t20\ta1\tx.bed\tchr1\t15\t25\tx1\n"));
    // count is the total over all B files
    assert_eq!(
        run(&["--mode", "count"]),
        "chr1\t10\t20\ta1\t3\nchr1\t50\t60\ta3\t0\n"
    );
    // R-P3: unmatched rows carry `.` as the label
    assert_eq!(
        run(&["--mode", "pairs", "--keep-unmatched", "--overlap-bases"]),
        "chr1\t10\t20\ta1\tx.bed\tchr1\t15\t25\tx1\t5\n\
chr1\t10\t20\ta1\ty.bed\tchr1\t5\t12\ty1\t2\n\
chr1\t10\t20\ta1\ty.bed\tchr1\t19\t30\ty2\t1\n\
chr1\t50\t60\ta3\t.\t.\t-1\t-1\t.\t0\n"
    );
    assert_eq!(
        code(&["intersect", "-a", &a, "-b", &x, &y, "--names", "only"]),
        2
    );
}

#[test]
fn thresholds() {
    let fx = Fx::new();
    let a = fx.file("a.bed", "chr1\t0\t10\ta\n");
    let b = fx.file("b.bed", "chr1\t5\t20\tb\n"); // 5 shared: 50% of A, 33% of B
    let any = |extra: &[&str]| {
        let mut v = vec!["intersect", "-a", &a, "-b", &b];
        v.extend_from_slice(extra);
        !out(&v).is_empty()
    };
    assert!(any(&["-f", "0.5"]));
    assert!(!any(&["-f", "0.6"]));
    assert!(!any(&["-F", "50%"]));
    assert!(any(&["-F", "30%"]));
    assert!(!any(&["-r", "-f", "0.5"]), "-r applies -f to B as well");
    assert!(any(&["-r", "-f", "0.3"]));
    assert!(any(&["-e", "-f", "0.9", "-F", "0.1"]));
    assert!(!any(&["-e", "-f", "0.9", "-F", "0.9"]));
    assert!(any(&["--min-bases", "5"]));
    assert!(!any(&["--min-bases", "6"]));
}

#[test]
fn strand() {
    let fx = Fx::new();
    let a = fx.file(
        "a.bed",
        "chr1\t0\t10\tap\t0\t+\nchr1\t20\t30\tam\t0\t-\nchr1\t40\t50\tau\t0\t.\n",
    );
    let b = fx.file(
        "b.bed",
        "chr1\t5\t8\tb1\t0\t+\nchr1\t22\t25\tb2\t0\t+\nchr1\t42\t45\tb3\t0\t+\n",
    );
    let names = |extra: &[&str]| -> (Vec<String>, String) {
        let mut v = vec!["intersect", "-a", &a, "-b", &b];
        v.extend_from_slice(extra);
        let o = gloci().args(&v).assert().success();
        let o = o.get_output();
        let text = String::from_utf8(o.stdout.clone()).unwrap();
        (
            text.lines()
                .map(|l| l.split('\t').nth(3).unwrap().to_string())
                .collect(),
            String::from_utf8(o.stderr.clone()).unwrap(),
        )
    };
    let (n, err) = names(&["-s"]);
    assert_eq!(n, ["ap"]);
    assert_eq!(err.matches("unknown strand").count(), 1, "{err}");
    let (n, err) = names(&["-S"]);
    assert_eq!(n, ["am"]);
    assert_eq!(err.matches("unknown strand").count(), 1, "{err}");
    let (n, err) = names(&[]);
    assert_eq!(n, ["ap", "am", "au"]);
    assert!(!err.contains("unknown strand"));
}

#[test]
fn windows() {
    let fx = Fx::new();
    let a = fx.file("a.bed", "chr1\t100\t110\ta\t0\t+\n");
    let am = fx.file("am.bed", "chr1\t100\t110\ta\t0\t-\n");
    let b = fx.file(
        "b.bed",
        "chr1\t85\t90\tfar_left\nchr1\t90\t95\tleft\nchr1\t115\t118\tright\nchr1\t150\t155\tfar\n",
    );
    let hit = |a: &str, extra: &[&str]| -> Vec<String> {
        let mut v = vec!["intersect", "-a", a, "-b", &b, "--mode", "pairs"];
        v.extend_from_slice(extra);
        out(&v)
            .lines()
            .map(|l| l.rsplit('\t').next().unwrap().to_string())
            .collect()
    };
    assert_eq!(hit(&a, &["-w", "10"]), ["left", "right"]);
    assert_eq!(hit(&a, &["-w", "20"]), ["far_left", "left", "right"]);
    // --left 5 --right 20: B [90,95) pads to [85,115) and [85,90) to
    // [80,110): both overlap A; [115,118) pads to [110,138): bookended.
    assert_eq!(
        hit(&a, &["--left", "5", "--right", "20"]),
        ["far_left", "left"]
    );
    // On a minus-strand A with --strand-window, left and right swap:
    // [115,118) pads to [95,123): overlaps; [90,95) pads to [70,100): no.
    assert_eq!(
        hit(&am, &["--left", "5", "--right", "20", "--strand-window"]),
        ["right"]
    );
    assert_eq!(
        hit(&am, &["--left", "5", "--right", "20"]),
        ["far_left", "left"]
    );
    assert_eq!(
        code(&["intersect", "-a", &a, "-b", &b, "-w", "10", "-f", "0.5"]),
        2
    );
}

const VHEAD: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

#[test]
fn vcf_a() {
    let fx = Fx::new();
    let a = fx.file(
        "a.vcf",
        &format!("{VHEAD}chr1\t12\t.\tA\tC\t.\t.\t.\nchr1\t100\t.\tG\tT\t.\t.\t.\n"),
    );
    let b = fx.file("b.bed", "chr1\t10\t20\tb\n");
    let any = out(&["intersect", "-a", &a, "-b", &b]);
    let lines: Vec<&str> = any.lines().collect();
    assert_eq!(lines[0], "##fileformat=VCFv4.2");
    assert!(lines.iter().any(|l| l.starts_with("##glociCommand=")));
    assert_eq!(lines.last().unwrap(), &"chr1\t12\t.\tA\tC\t.\t.\t.");
    assert_eq!(lines.len(), 4, "{any}");
    let pairs = out(&["intersect", "-a", &a, "-b", &b, "--mode", "pairs"]);
    assert_eq!(
        pairs,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tb_1\tb_2\tb_3\tb_4\n\
chr1\t12\t.\tA\tC\t.\t.\t.\tchr1\t10\t20\tb\n"
    );
    let count = out(&["intersect", "-a", &a, "-b", &b, "--mode", "count"]);
    assert!(count.starts_with("#CHROM\t"));
    assert!(count.ends_with("chr1\t100\t.\tG\tT\t.\t.\t.\t0\n"));
    gloci()
        .args(["intersect", "-a", &a, "-b", &b, "--mode", "clip"])
        .assert()
        .code(2)
        .stderr(contains("clip"));
}

fn vcf(recs: &[&str]) -> String {
    let mut t = VHEAD.to_string();
    for r in recs {
        t.push_str(r);
        t.push('\n');
    }
    t
}

#[test]
fn match_alleles() {
    let fx = Fx::new();
    let a = fx.file("a.vcf", &vcf(&["chr1\t100\t.\tA\tACGT\t.\t.\t."]));
    let hit = |b: &str| {
        let bf = fx.file("b.vcf", &vcf(&[b]));
        out(&["intersect", "-a", &a, "-b", &bf, "--match", "alleles"]).contains("chr1\t100\t")
    };
    assert!(hit("chr1\t100\t.\tA\tACGT,G\t.\t.\t."));
    assert!(!hit("chr1\t100\t.\tA\tAC\t.\t.\t."));
    assert!(!hit("chr1\t100\t.\tA\tG\t.\t.\t."));
    // the default match is by overlap, so the second record matches
    let bf = fx.file("b2.vcf", &vcf(&["chr1\t100\t.\tA\tAC\t.\t.\t."]));
    assert!(out(&["intersect", "-a", &a, "-b", &bf]).contains("chr1\t100\t"));

    // alleles are compared after trimming: ATG>AC at 100 is TG>C at 101
    let a2 = fx.file("a2.vcf", &vcf(&["chr1\t100\t.\tATG\tAC\t.\t.\t."]));
    let b2 = fx.file("b3.vcf", &vcf(&["chr1\t101\t.\tTG\tC\t.\t.\t."]));
    let o = out(&["intersect", "-a", &a2, "-b", &b2, "--match", "alleles"]);
    assert!(o.contains("chr1\t100\t.\tATG\tAC"), "{o}");

    let bed = fx.file("x.bed", "chr1\t0\t1000\n");
    assert_eq!(
        code(&["intersect", "-a", &a, "-b", &bed, "--match", "alleles"]),
        2
    );
    assert_eq!(
        code(&["intersect", "-a", &bed, "-b", &a, "--match", "alleles"]),
        2
    );
}

#[test]
fn empty_inputs() {
    let fx = Fx::new();
    let empty = fx.file("empty.bed", "");
    let b = fx.file("b.bed", B);
    let a = fx.file("a.bed", A);
    assert_eq!(
        out(&["intersect", "-a", &empty, "--a-format", "bed", "-b", &b]),
        ""
    );
    assert_eq!(
        out(&[
            "intersect",
            "-a",
            &a,
            "-b",
            &empty,
            "--b-format",
            "bed",
            "--mode",
            "none"
        ]),
        A
    );
    let av = fx.file("a.vcf", VHEAD);
    let o = out(&["intersect", "-a", &av, "-b", &b]);
    assert!(o.starts_with("##fileformat=VCFv4.2\n"));
    assert!(o.lines().all(|l| l.starts_with('#')), "{o}");
}

#[test]
fn order_errors() {
    let fx = Fx::new();
    let a = fx.file("a.bed", A);
    let bad_b = fx.file("badb.bed", "chr1\t50\t60\nchr1\t10\t20\n");
    gloci()
        .args(["intersect", "-a", &a, "-b", &bad_b])
        .assert()
        .code(1)
        .stderr(contains("badb.bed"));
    // --unsorted takes a shuffled A and keeps A's order
    let shuffled = fx.file(
        "shuf.bed",
        "chr2\t5\t15\ta4\t0\t+\nchr1\t30\t40\ta2\t0\t-\nchr1\t10\t20\ta1\t0\t+\n",
    );
    let b = fx.file("b.bed", B);
    assert_eq!(
        out(&["intersect", "-a", &shuffled, "-b", &b, "--unsorted"]),
        "chr2\t5\t15\ta4\t0\t+\nchr1\t30\t40\ta2\t0\t-\nchr1\t10\t20\ta1\t0\t+\n"
    );
    // an unsorted B is fine with --unsorted too
    let out_u = out(&["intersect", "-a", &a, "-b", &bad_b, "--unsorted"]);
    assert_eq!(out_u, "chr1\t10\t20\ta1\t0\t+\nchr1\t50\t60\ta3\t0\t+\n");
    gloci()
        .args(["intersect", "-a", &shuffled, "-b", &b])
        .assert()
        .code(1)
        .stderr(contains("shuf.bed"));
}

#[test]
fn determinism() {
    let fx = Fx::new();
    let a = fx.file("a.bed", A);
    let b = fx.file("b.bed", B);
    let b2 = fx.file("b2.bed", "chr1\t12\t14\tc1\t0\t+\n");
    for extra in [
        &["--mode", "any"][..],
        &["--mode", "none"],
        &["--mode", "pairs", "--keep-unmatched", "--overlap-bases"],
        &["--mode", "clip"],
        &["--mode", "count"],
        &["--mode", "count", "--unsorted"],
    ] {
        let run = |t: &str| {
            let mut v = vec!["intersect", "-a", &a, "-b", &b, &b2, "-t", t];
            v.extend_from_slice(extra);
            out(&v)
        };
        let one = run("1");
        assert!(!one.is_empty());
        assert_eq!(one, run("2"), "{extra:?}");
        assert_eq!(one, run("8"), "{extra:?}");
    }
}

#[test]
fn help_has_examples() {
    let h = out(&["intersect", "--help"]);
    assert!(h.contains("Examples:"));
    assert!(h.contains("default mode is any") || h.contains("The default mode is any"));
    assert!(h.contains("--unsorted"));
    assert!(h.contains("memory"));
}

#[test]
fn bcf_inputs_are_rendered_as_vcf_text() {
    let fx = Fx::new();
    let head = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let vcf_a = fx.file("a.vcf", &format!("{head}chr1\t100\t.\tA\tACGT\t.\t.\t.\n"));
    let vcf_b = fx.file(
        "b.vcf",
        &format!("{head}chr1\t100\t.\tA\tACGT,G\t.\t.\t.\n"),
    );
    let dir = fx.dir.path();
    let (a, b) = (s(&dir.join("a.bcf")), s(&dir.join("b.bcf")));
    out(&["view", &vcf_a, "-o", &a]);
    out(&["view", &vcf_b, "-o", &b]);
    // any with a BCF A writes BCF; pairs and count write TSV with VCF text
    let pairs = out(&[
        "intersect",
        "-a",
        &a,
        "-b",
        &b,
        "--mode",
        "pairs",
        "--match",
        "alleles",
    ]);
    assert_eq!(
        pairs,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tb_CHROM\tb_POS\tb_ID\tb_REF\tb_ALT\tb_QUAL\tb_FILTER\tb_INFO\n\
chr1\t100\t.\tA\tACGT\t.\t.\t.\tchr1\t100\t.\tA\tACGT,G\t.\t.\t.\n"
    );
    let count = out(&[
        "intersect",
        "-a",
        &a,
        "-b",
        &vcf_b,
        "--mode",
        "count",
        "--match",
        "alleles",
    ]);
    assert!(
        count.ends_with("chr1\t100\t.\tA\tACGT\t.\t.\t.\t1\n"),
        "{count}"
    );
    let o = dir.join("o.bcf");
    out(&[
        "intersect",
        "-a",
        &a,
        "-b",
        &b,
        "--match",
        "alleles",
        "-o",
        &s(&o),
    ]);
    let back = out(&["view", &s(&o), "-O", "vcf"]);
    assert!(back.contains("chr1\t100\t.\tA\tACGT"), "{back}");
}

/// The column-name line of a pairs report on a VCF A names every B column.
#[test]
fn pairs_header_names_every_b_column() {
    let fx = Fx::new();
    let a = fx.file("a.vcf", &vcf(&["chr1\t12\t.\tA\tC\t.\t.\t."]));
    let width = |text: &str| {
        let lines: Vec<&str> = text.lines().collect();
        let (h, r) = (lines[0], lines[1]);
        assert!(h.starts_with("#CHROM"), "{text}");
        assert_eq!(
            h.split('\t').count(),
            r.split('\t').count(),
            "names vs columns:\n{text}"
        );
        h.to_string()
    };
    let pairs = |b: &[&str], extra: &[&str]| {
        let mut v = vec!["intersect", "-a", &a, "-b"];
        v.extend_from_slice(b);
        v.extend_from_slice(&["--mode", "pairs"]);
        v.extend_from_slice(extra);
        out(&v)
    };
    // BED6 B
    let bed6 = fx.file("b6.bed", "chr1\t10\t20\tn\t0\t+\n");
    let h = width(&pairs(&[&bed6], &[]));
    assert!(h.ends_with("INFO\tb_1\tb_2\tb_3\tb_4\tb_5\tb_6"), "{h}");
    // VCF B with samples
    let vb = fx.file(
        "bs.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n\
chr1\t12\t.\tA\tC\t.\t.\t.\tGT\t0/1\t1/1\n",
    );
    let h = width(&pairs(&[&vb], &["--overlap-bases"]));
    assert!(
        h.ends_with("b_INFO\tb_FORMAT\tb_s1\tb_s2\toverlap_bases"),
        "{h}"
    );
    // several -b: label column first, widths from the first file
    let bed4 = fx.file("b4.bed", "chr1\t11\t13\tm\n");
    let h = width(&pairs(&[&bed6, &bed4], &[]));
    assert!(
        h.ends_with("INFO\tb_name\tb_1\tb_2\tb_3\tb_4\tb_5\tb_6"),
        "{h}"
    );
    // empty B: three names
    let empty = fx.file("e.bed", "");
    let o = out(&[
        "intersect",
        "-a",
        &a,
        "-b",
        &empty,
        "--b-format",
        "bed",
        "--mode",
        "pairs",
        "--keep-unmatched",
    ]);
    assert_eq!(
        o,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tb_1\tb_2\tb_3\n\
chr1\t12\t.\tA\tC\t.\t.\t.\t.\t-1\t-1\n"
    );
}

#[test]
fn match_alleles_bookended_and_multiallelic() {
    let fx = Fx::new();
    let run = |a: &str, b: &str, span: &str| {
        let (af, bf) = (fx.file("ba.vcf", &vcf(&[a])), fx.file("bb.vcf", &vcf(&[b])));
        out(&[
            "intersect",
            "-a",
            &af,
            "-b",
            &bf,
            "--match",
            "alleles",
            "--vcf-span",
            span,
        ])
        .lines()
        .any(|l| !l.starts_with('#'))
    };
    // the same insertion of G after position 100, written two ways whose
    // footprints are bookended under --vcf-span ref
    let (a, b) = (
        "chr1\t100\t.\tA\tAG\t.\t.\t.",
        "chr1\t101\t.\tC\tGC\t.\t.\t.",
    );
    assert!(run(a, b, "ref"));
    assert!(run(a, b, "allele"));
    // multi-allelic, hulls bookended under the allele span
    let (a, b) = (
        "chr1\t100\t.\tA\tAG,T\t.\t.\t.",
        "chr1\t101\t.\tC\tGC,T\t.\t.\t.",
    );
    assert!(run(a, b, "allele"));
    assert!(run(a, b, "ref"));
    // no shared allele: not a match, however close
    assert!(!run(a, "chr1\t101\t.\tC\tT,CC\t.\t.\t.", "allele"));
}

#[test]
fn match_alleles_refuses_overlap_options() {
    let fx = Fx::new();
    let a = fx.file("a.vcf", &vcf(&["chr1\t100\t.\tA\tC\t.\t.\t."]));
    for extra in [
        &["-f", "0.5"][..],
        &["-F", "0.5"],
        &["-r"],
        &["-e"],
        &["--min-bases", "1"],
        &["-s"],
        &["-S"],
        &["--split"],
        &["-w", "5"],
        &["--left", "5"],
        &["--right", "5"],
        &["--inclusive-points"],
    ] {
        let mut v = vec!["intersect", "-a", &a, "-b", &a, "--match", "alleles"];
        v.extend_from_slice(extra);
        gloci()
            .args(&v)
            .assert()
            .code(2)
            .stderr(contains("--match alleles"));
    }
}

/// A VCF A and a large BED6 B (more than the 64 KiB detection peek).
#[cfg(unix)]
fn vcf_a_and_big_bed(d: &Path) -> (String, String) {
    let mut vcf = String::from(
        "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=10000000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    let mut bed = String::new();
    for i in 0..5_000u64 {
        let p = 1_000 + i * 100;
        vcf.push_str(&format!("chr1\t{}\tv{i}\tA\tG\t.\tPASS\t.\n", p + 3));
        bed.push_str(&format!("chr1\t{p}\t{}\tfeature_{i}\t0\t+\n", p + 10));
    }
    (
        s(&plain_file(d, "a.vcf", &vcf)),
        s(&plain_file(d, "big.bed", &bed)),
    )
}

/// C4: a piped B with a VCF A loses no records (B's width is read ahead
/// from the reader already open, never by reopening the pipe).
#[cfg(unix)]
#[test]
fn piped_b_with_vcf_a_matches_regular_file() {
    let d = tempdir();
    let (a, b) = vcf_a_and_big_bed(d.path());
    for extra in ["--mode pairs", "--mode pairs --keep-unmatched"] {
        let want = bash_gloci(&format!("intersect -a '{a}' -b '{b}' {extra}"), "");
        assert!(
            want.status.success(),
            "{}",
            String::from_utf8_lossy(&want.stderr)
        );
        let got = bash_gloci(&format!("intersect -a '{a}' -b <(cat '{b}') {extra}"), "");
        assert!(
            got.status.success(),
            "{}",
            String::from_utf8_lossy(&got.stderr)
        );
        assert_eq!(
            String::from_utf8_lossy(&got.stdout),
            String::from_utf8_lossy(&want.stdout),
            "{extra}"
        );
        // stdin as B too.
        let bed = std::fs::read_to_string(&b).unwrap();
        let got = bash_gloci(&format!("intersect -a '{a}' -b - {extra}"), &bed);
        assert!(
            got.status.success(),
            "{}",
            String::from_utf8_lossy(&got.stderr)
        );
        assert_eq!(
            String::from_utf8_lossy(&got.stdout),
            String::from_utf8_lossy(&want.stdout),
            "stdin {extra}"
        );
    }
    let want = bash_gloci(&format!("intersect -a '{a}' -b '{b}' --mode pairs"), "");
    let text = String::from_utf8(want.stdout).unwrap();
    assert!(text.lines().next().unwrap().ends_with("\tb_6"), "{text}");
    assert_eq!(text.lines().count(), 5_001);
}
