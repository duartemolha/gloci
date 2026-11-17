//! bigWig and bigBed input (D171): read through the file's own index, in the command's
//! chromosome order.
mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;
use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};

use bigtools::bed::bedparser::BedValueError;
use bigtools::beddata::BedParserStreamingIterator;
use bigtools::{BedEntry, BigBedWrite, BigWigWrite, InputSortType, Value};

fn runtime() -> tokio::runtime::Runtime {
    tokio::runtime::Builder::new_current_thread()
        .build()
        .unwrap()
}

fn size_map(sizes: &[(&str, u32)]) -> HashMap<String, u32> {
    sizes.iter().map(|(c, n)| (c.to_string(), *n)).collect()
}

/// A bigWig written by bigtools directly (not gloci), records in the order given.
fn bigtools_bigwig(path: &Path, sizes: &[(&str, u32)], recs: &[(&str, u32, u32, f32)]) {
    let mut w = BigWigWrite::new(File::create(path).unwrap(), size_map(sizes));
    w.options.input_sort_type = InputSortType::START;
    let items: Vec<Result<(String, Value), BedValueError>> = recs
        .iter()
        .map(|&(c, start, end, value)| Ok((c.to_string(), Value { start, end, value })))
        .collect();
    w.write(
        BedParserStreamingIterator::wrap_iter(items.into_iter(), true),
        runtime(),
    )
    .unwrap();
}

fn bigtools_bigbed(path: &Path, sizes: &[(&str, u32)], recs: &[(&str, u32, u32, &str)]) {
    let mut w = BigBedWrite::new(File::create(path).unwrap(), size_map(sizes));
    w.options.input_sort_type = InputSortType::START;
    let items: Vec<Result<(String, BedEntry), BedValueError>> = recs
        .iter()
        .map(|&(c, start, end, rest)| {
            Ok((
                c.to_string(),
                BedEntry {
                    start,
                    end,
                    rest: rest.to_string(),
                },
            ))
        })
        .collect();
    w.write(
        BedParserStreamingIterator::wrap_iter(items.into_iter(), true),
        runtime(),
    )
    .unwrap();
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn k_file(dir: &Path) -> PathBuf {
    let p = dir.join("k.bw");
    bigtools_bigwig(
        &p,
        &[("chr1", 1000), ("chr2", 1000), ("chr10", 1000)],
        &[
            ("chr2", 0, 10, 0.1),
            ("chr2", 10, 20, 5.0),
            ("chr10", 5, 6, 1.0f32 / 3.0),
            ("chr1", 100, 200, -2.5),
            ("chr1", 300, 301, 1e-7),
        ],
    );
    p
}

fn b_file(dir: &Path) -> PathBuf {
    let p = dir.join("b.bb");
    bigtools_bigbed(
        &p,
        &[("chr1", 1000)],
        &[
            ("chr1", 10, 20, "a\t0\t+"),
            ("chr1", 15, 40, "b\t5\t-"),
            ("chr1", 500, 900, "c\t9\t."),
        ],
    );
    p
}

fn out(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

/// The chromosomes of `text`'s rows, each once, in order of appearance.
fn blocks(text: &str) -> Vec<&str> {
    let mut v: Vec<&str> = Vec::new();
    for l in text.lines() {
        let c = l.split('\t').next().unwrap();
        if v.last() != Some(&c) {
            v.push(c);
        }
    }
    v
}

const K_LEX: &str = "chr1\t100\t200\t-2.5\n\
chr1\t300\t301\t0.0000001\n\
chr10\t5\t6\t0.33333334\n\
chr2\t0\t10\t0.1\n\
chr2\t10\t20\t5\n";

#[test]
fn reads_bigwig_in_lex_order() {
    let d = tempdir();
    let k = k_file(d.path());
    assert_eq!(out(&["view", &s(&k)]), K_LEX);
}

#[test]
fn reads_in_the_commands_order() {
    let d = tempdir();
    let k = s(&k_file(d.path()));
    let karyo = out(&["view", "--chrom-order", "karyotype", &k]);
    assert_eq!(blocks(&karyo), ["chr1", "chr2", "chr10"]);
    let g = plain_file(
        d.path(),
        "order.genome",
        "chr10\t1000\nchr2\t1000\nchr1\t1000\n",
    );
    let by_g = out(&["view", "-g", &s(&g), &k]);
    assert_eq!(blocks(&by_g), ["chr10", "chr2", "chr1"]);
    let conv = out(&["convert", &k]);
    assert_eq!(blocks(&conv), ["chr2", "chr10", "chr1"]);
    assert_eq!(
        conv,
        "chr2\t0\t10\t0.1\nchr2\t10\t20\t5\nchr10\t5\t6\t0.33333334\n\
chr1\t100\t200\t-2.5\nchr1\t300\t301\t0.0000001\n"
    );
}

#[test]
fn map_reads_a_bigwig_b_in_lex_order() {
    let d = tempdir();
    let k = s(&k_file(d.path()));
    let a = plain_file(
        d.path(),
        "a.bed",
        "chr1\t0\t1000\nchr10\t0\t1000\nchr2\t0\t1000\n",
    );
    let text = out(&["map", "-a", &s(&a), "-b", &k, "-c", "4:count"]);
    let last: Vec<&str> = text
        .lines()
        .map(|l| l.rsplit('\t').next().unwrap())
        .collect();
    assert_eq!(last, ["2", "1", "2"], "{text}");
}

#[test]
fn empty_chromosomes_and_g_without_them() {
    let d = tempdir();
    let g3 = plain_file(
        d.path(),
        "g3.genome",
        "chr1\t5000\nchr2\t3000\nchr10\t2000\n",
    );
    let bw =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../xtask/diff/inputs/bigwig/pybigwig.bw");
    let text = out(&["view", "-g", &s(&g3), &s(&bw)]);
    assert_eq!(text.lines().count(), 9);
    assert_eq!(blocks(&text), ["chr1", "chr2", "chr10"]);
}

#[test]
fn bigbed_records_are_bed_lines() {
    let d = tempdir();
    let b = s(&b_file(d.path()));
    let want = "chr1\t10\t20\ta\t0\t+\nchr1\t15\t40\tb\t5\t-\nchr1\t500\t900\tc\t9\t.\n";
    assert_eq!(out(&["view", &b]), want);
    let q = plain_file(d.path(), "q.bed", "chr1\t0\t1000\tq\t0\t+\n");
    assert_eq!(
        out(&["intersect", "-a", &b, "-b", &s(&q), "-s"]),
        "chr1\t10\t20\ta\t0\t+\n"
    );
    assert_eq!(out(&["convert", &b]), want);
}

#[test]
fn chrom_map_renames_bigwig_chromosomes() {
    let d = tempdir();
    let k = s(&k_file(d.path()));
    let m = plain_file(d.path(), "map.tsv", "chr1\t1\nchr2\t2\nchr10\t10\n");
    let text = out(&["view", "--chrom-map", &s(&m), &k]);
    assert_eq!(blocks(&text), ["1", "10", "2"]);
}

#[test]
fn bigwig_must_come_from_a_file() {
    let d = tempdir();
    let k = k_file(d.path());
    gloci()
        .args(["view", "-"])
        .write_stdin(std::fs::read(&k).unwrap())
        .assert()
        .code(2)
        .stderr(contains("bigWig"))
        .stderr(contains("file"));
    gloci()
        .args(["view", "/dev/stdin"])
        .pipe_stdin(&k)
        .unwrap()
        .assert()
        .code(2)
        .stderr(contains("bigWig"))
        .stderr(contains("file"));
}

#[test]
fn format_override() {
    let d = tempdir();
    let k = s(&k_file(d.path()));
    let a = plain_file(d.path(), "a.bed", "chr1\t1\t2\n");
    let a = s(&a);
    gloci()
        .args(["view", "--format", "bigwig", &a])
        .assert()
        .failure()
        .stderr(contains(&a))
        .stderr(contains("not a bigWig"));
    gloci()
        .args(["view", "--format", "bigbed", &k])
        .assert()
        .failure()
        .stderr(contains("not a bigBed"));
    gloci()
        .args(["view", "--format", "bigwig", &k])
        .assert()
        .success();
}

#[test]
fn truncated_or_corrupt_bigwig_is_a_data_error() {
    let d = tempdir();
    let p = d.path().join("big.bw");
    let recs: Vec<(&str, u32, u32, f32)> = (0..3000u32)
        .map(|i| ("chr1", i * 10, i * 10 + 5, i as f32))
        .collect();
    bigtools_bigwig(&p, &[("chr1", 100_000)], &recs);
    let bytes = std::fs::read(&p).unwrap();
    for len in [bytes.len() / 2, 40] {
        let t = d.path().join(format!("cut{len}.bw"));
        std::fs::write(&t, &bytes[..len]).unwrap();
        let t = s(&t);
        gloci()
            .args(["view", &t])
            .assert()
            .code(1)
            .stderr(contains(&t))
            .stderr(contains("truncated"))
            .stderr(predicates::str::contains("internal error").not());
    }
    let mut bad = bytes.clone();
    bad[1024..1088].fill(0xFF);
    let c = d.path().join("bad.bw");
    std::fs::write(&c, &bad).unwrap();
    gloci()
        .args(["view", &s(&c)])
        .assert()
        .code(1)
        .stderr(predicates::str::contains("internal error").not());
}

#[test]
fn index_refuses_bigwig() {
    let d = tempdir();
    let k = s(&k_file(d.path()));
    gloci()
        .args(["index", &k])
        .assert()
        .code(2)
        .stderr(contains("own index"));
}

#[test]
fn reading_does_not_depend_on_threads() {
    let d = tempdir();
    let k = s(&k_file(d.path()));
    assert_eq!(out(&["view", "-t", "1", &k]), out(&["view", "-t", "8", &k]));
}

#[test]
fn golden_expected_is_the_hand_written_listing() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let listing =
        std::fs::read_to_string(root.join("xtask/diff/inputs/bigwig/pybigwig.lex.bedgraph"))
            .unwrap();
    let golden = root.join("crates/gloci/tests/golden/bigwig-view-pybigwig-file-vs-pybigwig");
    let expected = std::fs::read_to_string(golden.join("expected")).unwrap();
    assert_eq!(expected, listing);
}

#[test]
fn option_checks_and_alignment_commands_refuse_bbi_by_usage_error() {
    let d = tempdir();
    let k = s(&k_file(d.path()));
    for flag in ["--blocks", "--split-alleles", "--bed12"] {
        gloci()
            .args(["convert", &k, flag])
            .assert()
            .code(2)
            .stderr(predicates::str::contains("is empty").not());
    }
    gloci()
        .args(["flagstat", &k])
        .assert()
        .code(2)
        .stderr(contains("bigWig"))
        .stderr(predicates::str::contains("is empty").not());
}

#[test]
fn check_quick_reads_bbi() {
    let d = tempdir();
    let k = k_file(d.path());
    let b = b_file(d.path());
    for f in [&k, &b] {
        gloci().args(["check", "--quick", &s(f)]).assert().success();
    }
    let bytes = std::fs::read(&k).unwrap();
    let cut = d.path().join("cut.bw");
    std::fs::write(&cut, &bytes[..bytes.len() / 2]).unwrap();
    gloci()
        .args(["check", "--quick", &s(&cut)])
        .assert()
        .code(1)
        .stdout(contains("truncated").or(contains("not a valid")))
        .stderr(predicates::str::contains("internal error").not());
}

#[test]
fn errors_locate_bbi_records_by_record() {
    let d = tempdir();
    let k = s(&k_file(d.path()));
    let g = plain_file(
        d.path(),
        "small.genome",
        "chr1\t150\nchr2\t1000\nchr10\t1000\n",
    );
    gloci()
        .args(["genomecov", &k, "-g", &s(&g)])
        .assert()
        .failure()
        .stderr(contains("record "))
        .stderr(predicates::str::contains("line ").not());
}

#[test]
fn truncated_bigbed_and_format_override_and_missing_chromosome() {
    let d = tempdir();
    let b = b_file(d.path());
    gloci()
        .args(["view", "--format", "bigbed", &s(&b)])
        .assert()
        .success();
    let bytes = std::fs::read(&b).unwrap();
    let cut = d.path().join("cut.bb");
    std::fs::write(&cut, &bytes[..40]).unwrap();
    gloci()
        .args(["view", &s(&cut)])
        .assert()
        .code(1)
        .stderr(contains("not a valid bigBed"));
    let k = s(&k_file(d.path()));
    let g = plain_file(d.path(), "no10.genome", "chr1\t1000\nchr2\t1000\n");
    gloci()
        .args(["view", "-g", &s(&g), &k])
        .assert()
        .code(1)
        .stderr(contains("chr10"));
}

fn two_record_bigwig(dir: &Path) -> String {
    let p = dir.join("r.bw");
    bigtools_bigwig(
        &p,
        &[("chr1", 1000)],
        &[("chr1", 0, 100, 1.5), ("chr1", 100, 200, 0.25)],
    );
    s(&p)
}

const CLIPPED: &str = "chr1\t50\t100\t1.5\nchr1\t100\t150\t0.25\n";

/// stderr must never carry the advice for a text file.
fn no_text_advice(a: assert_cmd::assert::Assert) -> assert_cmd::assert::Assert {
    a.stderr(predicates::str::contains("scanning the whole file").not())
        .stderr(predicates::str::contains("gloci sort").not())
}

#[test]
fn region_read_clips_bigwig_records() {
    let d = tempdir();
    let r = two_record_bigwig(d.path());
    let a = gloci().args(["view", &r, "chr1:51-150"]).assert().success();
    let a = no_text_advice(a);
    assert_eq!(String::from_utf8_lossy(&a.get_output().stdout), CLIPPED);
}

#[test]
fn bigwig_record_spanning_two_regions_gives_two_pieces() {
    let d = tempdir();
    let p = d.path().join("s.bw");
    bigtools_bigwig(&p, &[("chr1", 1000)], &[("chr1", 0, 100, 1.5)]);
    assert_eq!(
        out(&["view", &s(&p), "chr1:11-20", "chr1:31-40"]),
        "chr1\t10\t20\t1.5\nchr1\t30\t40\t1.5\n"
    );
}

#[test]
fn bigbed_region_returns_whole_records_once() {
    let d = tempdir();
    let p = d.path().join("w.bb");
    bigtools_bigbed(&p, &[("chr1", 5000)], &[("chr1", 10, 2000, "big")]);
    let a = gloci()
        .args(["view", &s(&p), "chr1:21-30", "chr1:1001-1010"])
        .assert()
        .success();
    let a = no_text_advice(a);
    assert_eq!(
        String::from_utf8_lossy(&a.get_output().stdout),
        "chr1\t10\t2000\tbig\n"
    );
}

#[test]
fn bigbed_region_excludes_bookended() {
    let d = tempdir();
    let p = d.path().join("x.bb");
    bigtools_bigbed(&p, &[("chr1", 5000)], &[("chr1", 0, 10, "x")]);
    assert_eq!(out(&["view", &s(&p), "chr1:11-20"]), "");
}

#[test]
fn regions_file_reads_bigwig() {
    let d = tempdir();
    let r = two_record_bigwig(d.path());
    let reg = plain_file(d.path(), "reg.bed", "chr1\t50\t150\n");
    let a = gloci()
        .args(["view", "--regions-file", &s(&reg), &r])
        .assert()
        .success();
    let a = no_text_advice(a);
    assert_eq!(String::from_utf8_lossy(&a.get_output().stdout), CLIPPED);
}

#[test]
fn region_on_a_missing_chromosome_warns() {
    let d = tempdir();
    let r = two_record_bigwig(d.path());
    let a = gloci().args(["view", &r, "chrZ:1-10"]).assert().success();
    let a = no_text_advice(a)
        .stderr(contains("chrZ"))
        .stderr(contains("r.bw"));
    assert!(a.get_output().stdout.is_empty());
}

#[test]
fn region_names_go_through_the_chrom_map() {
    let d = tempdir();
    let r = two_record_bigwig(d.path());
    let m = plain_file(d.path(), "map.tsv", "chr1\t1\n");
    assert_eq!(
        out(&["view", "--chrom-map", &s(&m), &r, "1:51-150"]),
        CLIPPED.replace("chr1", "1")
    );
}

/// A bigBed whose only record is the point `chr1 50 50 "p"`. bigtools' writer does not
/// finish on a zero-length record, so it writes `[50, 51)` uncompressed and the stored end is
/// patched down to 50.
fn point_bigbed(path: &Path) {
    let mut w = BigBedWrite::new(File::create(path).unwrap(), size_map(&[("chr1", 5000)]));
    w.options.input_sort_type = InputSortType::START;
    w.options.compress = false;
    let items: Vec<Result<(String, BedEntry), BedValueError>> = vec![Ok((
        "chr1".to_string(),
        BedEntry {
            start: 50,
            end: 51,
            rest: "p".to_string(),
        },
    ))];
    w.write(
        BedParserStreamingIterator::wrap_iter(items.into_iter(), true),
        runtime(),
    )
    .unwrap();
    let mut bytes = std::fs::read(path).unwrap();
    let pat: Vec<u8> = [0u32, 50, 51]
        .iter()
        .flat_map(|v| v.to_le_bytes())
        .collect();
    let at = bytes
        .windows(pat.len())
        .position(|w| w == pat)
        .expect("the record is stored uncompressed");
    bytes[at + 8..at + 12].copy_from_slice(&50u32.to_le_bytes());
    std::fs::write(path, bytes).unwrap();
}

#[test]
fn bigbed_region_points() {
    let d = tempdir();
    let p = d.path().join("p.bb");
    point_bigbed(&p);
    assert_eq!(out(&["view", &s(&p)]), "chr1\t50\t50\tp\n");
    assert_eq!(out(&["view", &s(&p), "chr1:1-50"]), "");
    assert_eq!(
        out(&["view", "--inclusive-points", &s(&p), "chr1:1-50"]),
        "chr1\t50\t50\tp\n"
    );
}

#[test]
fn zero_length_region_on_a_bigwig_yields_nothing() {
    let d = tempdir();
    let r = two_record_bigwig(d.path());
    // `chr1:5-5` is the one base 5, not a zero-length region.
    assert_eq!(out(&["view", &r, "chr1:5-5"]), "chr1\t4\t5\t1.5\n");
    let reg = plain_file(d.path(), "zero.bed", "chr1\t50\t50\n");
    assert_eq!(out(&["view", "--regions-file", &s(&reg), &r]), "");
}

// ---- convert to and from bigWig and bigBed ----

fn inputs() -> PathBuf {
    diff_inputs().join("bigwig")
}

fn input(name: &str) -> String {
    s(&inputs().join(name))
}

fn text(name: &str) -> String {
    std::fs::read_to_string(inputs().join(name)).unwrap()
}

fn convert_ok(args: &[&str]) {
    gloci().arg("convert").args(args).assert().success();
}

/// The files `dir` holds, by name.
fn names(dir: &Path) -> Vec<String> {
    let mut v: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    v.sort();
    v
}

fn fail_with(code: i32, args: &[&str], needles: &[&str]) {
    let mut a = gloci().args(args).assert().code(code);
    for n in needles {
        a = a.stderr(contains(*n));
    }
}

#[test]
fn convert_bedgraph_to_bigwig_and_back() {
    let d = tempdir();
    let bw = s(&d.path().join("t.bw"));
    convert_ok(&[
        &input("signal.bedgraph"),
        "-g",
        &input("signal.genome"),
        "-o",
        &bw,
    ]);
    assert_eq!(out(&["view", &bw]), text("signal.bedgraph"));
    let back = d.path().join("t.bedgraph");
    convert_ok(&[&bw, "-o", &s(&back)]);
    assert_eq!(
        std::fs::read_to_string(&back).unwrap(),
        text("signal.bedgraph")
    );
}

#[test]
fn convert_wig_to_bigwig() {
    let d = tempdir();
    let wig = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/wig-golden");
    let genome = plain_file(d.path(), "g.genome", "chr1\t1000\nchr2\t1000\n");
    let bw = s(&d.path().join("w.bw"));
    convert_ok(&[&s(&wig.join("input.wig")), "-g", &s(&genome), "-o", &bw]);
    let expected = std::fs::read_to_string(wig.join("expected.bedgraph")).unwrap();
    // `view` lists chromosomes lexicographically.
    let mut rows: Vec<&str> = expected.lines().collect();
    rows.sort_by_key(|l| l.split('\t').next().unwrap().to_string());
    let want: String = rows.iter().map(|l| format!("{l}\n")).collect();
    assert_eq!(out(&["view", &bw]), want);
}

#[test]
fn convert_bedgraph_detected_as_bed() {
    let d = tempdir();
    let g = plain_file(d.path(), "g.genome", "chr1\t1000\n");
    let four = plain_file(d.path(), "four.bed", "chr1\t0\t10\t1.5\nchr1\t20\t30\t2\n");
    let bw = s(&d.path().join("a.bw"));
    convert_ok(&[&s(&four), "-g", &s(&g), "-o", &bw]);
    assert_eq!(out(&["view", &bw]), "chr1\t0\t10\t1.5\nchr1\t20\t30\t2\n");
    let name = plain_file(d.path(), "name.bed", "chr1\t0\t10\tgeneA\n");
    fail_with(
        1,
        &[
            "convert",
            &s(&name),
            "-g",
            &s(&g),
            "-o",
            &s(&d.path().join("b.bw")),
        ],
        &["line 1", "not a number"],
    );
    let six = plain_file(d.path(), "six.bed", "chr1\t0\t10\tn\t5\t+\n");
    fail_with(
        1,
        &[
            "convert",
            &s(&six),
            "-g",
            &s(&g),
            "-o",
            &s(&d.path().join("c.bw")),
        ],
        &["4 columns"],
    );
}

#[test]
fn convert_bigwig_to_bigwig_keeps_its_sizes() {
    let d = tempdir();
    let t = s(&d.path().join("t.bw"));
    let u = s(&d.path().join("u.bw"));
    convert_ok(&[
        &input("signal.bedgraph"),
        "-g",
        &input("signal.genome"),
        "-o",
        &t,
    ]);
    convert_ok(&[&t, "-o", &u]);
    let lens = |p: &str| {
        let mut v = gloci_io::bbi_chrom_sizes(p).unwrap();
        v.sort();
        v
    };
    assert_eq!(lens(&u), lens(&t));
    assert_eq!(
        lens(&u),
        vec![
            (b"chr1".to_vec(), 40000),
            (b"chr10".to_vec(), 6000),
            (b"chr2".to_vec(), 2000)
        ]
    );
    assert_eq!(out(&["view", &u]), text("signal.bedgraph"));
}

#[test]
fn convert_bed_to_bigbed_and_back() {
    let d = tempdir();
    for (bed, genome) in [
        ("peaks.bed", "peaks.genome"),
        ("genes12.bed", "peaks.genome"),
    ] {
        let bb = s(&d.path().join(format!("{bed}.bb")));
        convert_ok(&[&input(bed), "-g", &input(genome), "-o", &bb]);
        assert_eq!(out(&["view", &bb]), text(bed), "{bed} through view");
        let back = d.path().join(format!("{bed}.back"));
        convert_ok(&[&bb, "-o", &format!("{}.bed", s(&back))]);
        assert_eq!(
            std::fs::read_to_string(format!("{}.bed", s(&back))).unwrap(),
            text(bed),
            "{bed} through convert"
        );
    }
}

#[test]
fn convert_autosql_option() {
    let d = tempdir();
    let as8 = "table peak\n\"my peaks\"\n(\n    string chrom; \"c\"\n    uint chromStart; \"s\"\n    uint chromEnd; \"e\"\n    string name; \"n\"\n    uint score; \"sc\"\n    char[1] strand; \"st\"\n    float signal; \"sg\"\n    string f1; \"x\"\n)\n";
    let as_file = plain_file(d.path(), "my.as", as8);
    let bb = d.path().join("a.bb");
    convert_ok(&[
        &input("peaks.bed"),
        "-g",
        &input("peaks.genome"),
        "-o",
        &s(&bb),
        "--autosql",
        &s(&as_file),
    ]);
    let mut r = bigtools::BigBedRead::open_file(&bb).unwrap();
    assert_eq!(r.autosql().unwrap().unwrap(), as8);
    assert_eq!(out(&["view", &s(&bb)]), text("peaks.bed"));
    let five = plain_file(
        d.path(),
        "five.as",
        "table t\n\"t\"\n(\n    string chrom; \"c\"\n    uint chromStart; \"s\"\n    uint chromEnd; \"e\"\n    string name; \"n\"\n    uint score; \"sc\"\n)\n",
    );
    fail_with(
        1,
        &[
            "convert",
            &input("peaks.bed"),
            "-g",
            &input("peaks.genome"),
            "-o",
            &s(&d.path().join("b.bb")),
            "--autosql",
            &s(&five),
        ],
        &["declares 5 fields", "8 columns"],
    );
    // A bad autoSql is a usage error naming its file; nothing is left behind.
    let bad = plain_file(d.path(), "bad.as", "table bed\n(\n  nonsense\n)");
    let before = names(d.path());
    fail_with(
        2,
        &[
            "convert",
            &input("peaks.bed"),
            "-g",
            &input("peaks.genome"),
            "-o",
            &s(&d.path().join("c.bb")),
            "--autosql",
            &s(&bad),
        ],
        &["bad.as"],
    );
    assert_eq!(names(d.path()), before);
    fail_with(
        2,
        &[
            "convert",
            &input("signal.bedgraph"),
            "-g",
            &input("signal.genome"),
            "-o",
            &s(&d.path().join("x.bw")),
            "--autosql",
            &s(&as_file),
        ],
        &["--autosql"],
    );
}

#[test]
fn convert_to_bigwig_needs_sizes() {
    let d = tempdir();
    fail_with(
        2,
        &[
            "convert",
            &input("signal.bedgraph"),
            "-o",
            &s(&d.path().join("x.bw")),
        ],
        &["-g"],
    );
}

#[test]
fn convert_errors_name_the_input_line() {
    let d = tempdir();
    let g = plain_file(d.path(), "g.genome", "chr1\t1000\nchr2\t1000\n");
    let run = |name: &str, body: &str, needles: &[&str]| {
        let f = plain_file(d.path(), name, body);
        let o = s(&d.path().join(format!("{name}.bw")));
        fail_with(1, &["convert", &s(&f), "-g", &s(&g), "-o", &o], needles);
        assert!(!Path::new(&o).exists());
    };
    run(
        "overlap.bedgraph",
        "chr1\t0\t10\t1\nchr1\t20\t30\t2\nchr1\t25\t40\t3\n",
        &["line 3", "overlap"],
    );
    run(
        "split.bedgraph",
        "chr1\t0\t10\t1\nchr2\t0\t10\t1\nchr1\t20\t30\t3\n",
        &["gloci sort"],
    );
    run(
        "missing.bedgraph",
        "chr1\t0\t10\t1\nchrZ\t0\t10\t1\n",
        &["chrZ", "-g", "line 2"],
    );
}

#[test]
fn convert_track_option_conflicts() {
    let d = tempdir();
    let bg = input("signal.bedgraph");
    let gen_ = input("signal.genome");
    let x = |n: &str| s(&d.path().join(n));
    fail_with(
        2,
        &["convert", &input("peaks.bed"), "--blocks", "-o", &x("x.bb")],
        &["--blocks does not apply to a bigBed output"],
    );
    fail_with(
        2,
        &["convert", &bg, "--bam", "-g", "hg38", "-o", &x("x.bw")],
        &["--bam chooses another output format"],
    );
    fail_with(
        2,
        &[
            "convert",
            &bg,
            "-g",
            &gen_,
            "--write-index",
            "-o",
            &x("x.bw"),
        ],
        &["--write-index does not apply", "carries its own index"],
    );
    fail_with(
        2,
        &["convert", &bg, "-g", "hg38", "-o", &x("x.bed")],
        &["-g", "--bam", "bigWig"],
    );
    let vcf = plain_file(
        d.path(),
        "v.vcf",
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t5\t.\tA\tC\t.\t.\t.\n",
    );
    fail_with(
        2,
        &["convert", &s(&vcf), "-g", &gen_, "-o", &x("v.bw")],
        &[
            "a bigWig from bedGraph, a 4-column BED, WIG or bigWig",
            "not from VCF",
        ],
    );
    let typed = plain_file(
        d.path(),
        "typed.bedgraph",
        "track type=bedGraph\nchr1\t0\t10\t1\n",
    );
    fail_with(
        2,
        &["convert", &s(&typed), "-g", &gen_, "-o", &x("t.bb")],
        &["bedGraph"],
    );
    assert_eq!(names(d.path()), vec!["typed.bedgraph", "v.vcf"]);
}

#[test]
fn other_commands_refuse_bigwig_output() {
    let d = tempdir();
    fail_with(
        2,
        &[
            "sort",
            &input("signal.bedgraph"),
            "-o",
            &s(&d.path().join("y.bw")),
        ],
        &["gloci convert"],
    );
    fail_with(
        2,
        &[
            "view",
            &input("peaks.bed"),
            "-o",
            &s(&d.path().join("y.bb")),
        ],
        &["gloci convert"],
    );
}

#[test]
fn bigwig_output_must_be_a_file() {
    fail_with(
        2,
        &[
            "convert",
            &input("signal.bedgraph"),
            "-g",
            &input("signal.genome"),
            "-O",
            "bigwig",
        ],
        &["regular file"],
    );
}

#[test]
fn convert_threads_do_not_change_the_file() {
    let d = tempdir();
    for (src, genome, ext) in [
        ("signal.bedgraph", "signal.genome", "bw"),
        ("peaks.bed", "peaks.genome", "bb"),
    ] {
        let mut files = Vec::new();
        for t in ["1", "2", "8"] {
            let o = d.path().join(format!("t{t}.{ext}"));
            convert_ok(&[&input(src), "-g", &input(genome), "-t", t, "-o", &s(&o)]);
            files.push(std::fs::read(&o).unwrap());
        }
        assert_eq!(files[0], files[1], "{src}: -t 1 and -t 2 differ");
        assert_eq!(files[0], files[2], "{src}: -t 1 and -t 8 differ");
    }
}

#[test]
fn failed_convert_leaves_nothing() {
    let d = tempdir();
    let mut body = String::new();
    for i in 0..2999u64 {
        body.push_str(&format!("chr1\t{}\t{}\t1\n", i * 2, i * 2 + 1));
    }
    // Line 3,000 overlaps line 2,999.
    body.push_str(&format!("chr1\t{}\t{}\t1\n", 2998 * 2, 2998 * 2 + 5));
    let f = plain_file(d.path(), "long.bedgraph", &body);
    let g = plain_file(d.path(), "g.genome", "chr1\t100000\n");
    fail_with(
        1,
        &[
            "convert",
            &s(&f),
            "-g",
            &s(&g),
            "-o",
            &s(&d.path().join("x.bw")),
        ],
        &["line 3,000", "overlap"],
    );
    assert_eq!(names(d.path()), vec!["g.genome", "long.bedgraph"]);
}

#[test]
fn empty_input_to_bigwig_is_an_error() {
    let d = tempdir();
    let f = plain_file(d.path(), "empty.bedgraph", "");
    fail_with(
        1,
        &[
            "convert",
            &s(&f),
            "-g",
            &input("signal.genome"),
            "-o",
            &s(&d.path().join("e.bw")),
        ],
        &["no records"],
    );
    assert_eq!(names(d.path()), vec!["empty.bedgraph"]);
}

#[test]
fn golden_round_trips() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    for (case, source) in [
        ("bigwig-written-by-gloci-vs-pybigwig", "signal.bedgraph"),
        ("bigbed-written-by-gloci-vs-pybigwig", "peaks.bed"),
        ("bigbed12-written-by-gloci-vs-pybigwig", "genes12.bed"),
        ("bigbed-user-autosql-vs-pybigwig", "peaks.bed"),
    ] {
        let expected = std::fs::read_to_string(
            root.join("crates/gloci/tests/golden")
                .join(case)
                .join("expected"),
        )
        .unwrap();
        assert_eq!(expected, text(source), "{case}");
    }
}

#[test]
fn skip_invalid_warns_for_a_bigwig_output() {
    let d = tempdir();
    let f = plain_file(
        d.path(),
        "s.bedgraph",
        "chr1\t0\t10\t1\nchr1\tx\t30\t2\nchr1\t40\t50\t3\n",
    );
    let bw = s(&d.path().join("s.bw"));
    gloci()
        .args([
            "convert",
            "--skip-invalid",
            &s(&f),
            "-g",
            &input("signal.genome"),
            "-o",
            &bw,
        ])
        .assert()
        .success()
        .stderr(contains("skipped 1 invalid"));
    assert_eq!(out(&["view", &bw]), "chr1\t0\t10\t1\nchr1\t40\t50\t3\n");
}

#[test]
fn skip_invalid_never_hides_a_writer_refusal() {
    let d = tempdir();
    let f = plain_file(d.path(), "z.bedgraph", "chr1\t0\t10\t1\nchr1\t20\t20\t2\n");
    fail_with(
        1,
        &[
            "convert",
            "--skip-invalid",
            &s(&f),
            "-g",
            &input("signal.genome"),
            "-o",
            &s(&d.path().join("z.bw")),
        ],
        &["line 2", "zero-length"],
    );
    assert_eq!(names(d.path()), vec!["z.bedgraph"]);
}

#[test]
fn reference_with_a_track_output_warns() {
    let d = tempdir();
    let fa = plain_file(d.path(), "r.fa", ">chr1\nACGT\n");
    gloci()
        .args([
            "convert",
            &input("signal.bedgraph"),
            "-g",
            &input("signal.genome"),
            "--reference",
            &s(&fa),
            "-o",
            &s(&d.path().join("r.bw")),
        ])
        .assert()
        .success()
        .stderr(contains("is not used: no input is CRAM"));
}

#[test]
fn the_first_bad_line_is_the_one_reported() {
    let d = tempdir();
    let g = plain_file(d.path(), "g.genome", "chr1\t100\nchr2\t100\n");
    // Line 2 is past its chromosome; line 4 returns to chr1 after chr2.
    let f = plain_file(
        d.path(),
        "o.bedgraph",
        "chr1\t0\t10\t1\nchr1\t90\t200\t1\nchr2\t0\t10\t1\nchr1\t50\t60\t1\n",
    );
    fail_with(
        1,
        &[
            "convert",
            &s(&f),
            "-g",
            &s(&g),
            "-o",
            &s(&d.path().join("o.bw")),
        ],
        &["line 2", "past the end"],
    );
}

#[test]
fn chrom_map_renames_the_header_sizes() {
    let d = tempdir();
    let bw = s(&d.path().join("c.bw"));
    convert_ok(&[
        &input("signal.bedgraph"),
        "-g",
        &input("signal.genome"),
        "-o",
        &bw,
    ]);
    let map = plain_file(d.path(), "map.tsv", "chr1\tone\nchr10\tten\nchr2\ttwo\n");
    let out_bw = d.path().join("m.bw");
    convert_ok(&[&bw, "--chrom-map", &s(&map), "-o", &s(&out_bw)]);
    let mut sizes = gloci_io::bbi_chrom_sizes(&s(&out_bw)).unwrap();
    sizes.sort();
    assert_eq!(
        sizes,
        vec![
            (b"one".to_vec(), 40000),
            (b"ten".to_vec(), 6000),
            (b"two".to_vec(), 2000)
        ]
    );
    assert!(out(&["view", &s(&out_bw)]).starts_with("one\t0\t7\t-2.5\n"));
}

#[test]
fn bigbed_to_bigbed_keeps_the_input_autosql() {
    let d = tempdir();
    let sql = "table sig\n\"signals\"\n(\n    string chrom; \"c\"\n    uint chromStart; \"s\"\n    uint chromEnd; \"e\"\n    string name; \"n\"\n    uint score; \"sc\"\n    char[1] strand; \"st\"\n    float signal; \"sg\"\n    string tag; \"t\"\n)\n";
    let as_file = plain_file(d.path(), "sig.as", sql);
    let a = d.path().join("a.bb");
    convert_ok(&[
        &input("peaks.bed"),
        "-g",
        &input("peaks.genome"),
        "--autosql",
        &s(&as_file),
        "-o",
        &s(&a),
    ]);
    let b = d.path().join("b.bb");
    convert_ok(&[&s(&a), "-o", &s(&b)]);
    let mut r = bigtools::BigBedRead::open_file(&b).unwrap();
    assert_eq!(r.autosql().unwrap().unwrap(), sql);
    assert_eq!(out(&["view", &s(&b)]), text("peaks.bed"));
    // An explicit --autosql replaces it.
    let other = sql.replace("signal", "sigval");
    let o_file = plain_file(d.path(), "o.as", &other);
    let c = d.path().join("c.bb");
    convert_ok(&[&s(&a), "--autosql", &s(&o_file), "-o", &s(&c)]);
    let mut r = bigtools::BigBedRead::open_file(&c).unwrap();
    assert_eq!(r.autosql().unwrap().unwrap(), other);
}

#[test]
fn the_committed_user_autosql_bigbed_stores_its_autosql() {
    let mut r = bigtools::BigBedRead::open_file(inputs().join("peaks6.bb")).unwrap();
    assert_eq!(r.autosql().unwrap().unwrap(), text("peaks6.as"));
}

// ---- genomecov to bigWig ----

/// `gloci genomecov ARGS` stdout, sorted the way a bigWig lists its chromosomes.
fn sorted_genomecov(args: &[&str]) -> String {
    let text = out(&[&["genomecov"], args].concat());
    let mut lines: Vec<(&str, u64, &str)> = text
        .lines()
        .map(|l| {
            let mut f = l.splitn(3, '\t');
            let chrom = f.next().unwrap();
            let start = f.next().unwrap().parse().unwrap();
            (chrom, start, l)
        })
        .collect();
    lines.sort_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));
    lines.iter().map(|l| format!("{}\n", l.2)).collect()
}

fn genomecov_inputs() -> (String, String, String) {
    let i = diff_inputs();
    (
        s(&i.join("gc-a.bed")),
        s(&i.join("gc.genome")),
        s(&i.join("aln/cigars-bases.bam")),
    )
}

#[test]
fn genomecov_bga_to_bigwig_matches_bedgraph() {
    let d = tempdir();
    let (bed, genome, bam) = genomecov_inputs();
    let bw = s(&d.path().join("g.bw"));
    gloci()
        .args(["genomecov", &bed, "-g", &genome, "--bga", "-o", &bw])
        .assert()
        .success();
    assert_eq!(
        out(&["view", &bw]),
        sorted_genomecov(&[&bed, "-g", &genome, "--bga"])
    );
    let bw2 = s(&d.path().join("b.bw"));
    gloci()
        .args(["genomecov", &bam, "--bga", "-o", &bw2])
        .assert()
        .success();
    assert_eq!(out(&["view", &bw2]), sorted_genomecov(&[&bam, "--bga"]));
}

#[test]
fn genomecov_bedgraph_to_bigwig() {
    let d = tempdir();
    let (bed, genome, bam) = genomecov_inputs();
    for (n, extra) in [("a", vec![]), ("b", vec!["--scale", "0.5"])] {
        for (input, g) in [(&bed, vec!["-g", genome.as_str()]), (&bam, vec![])] {
            let bw = s(&d.path().join(format!("{n}-{}.bw", g.len())));
            let mut args = vec!["genomecov", input.as_str()];
            args.extend(&g);
            args.push("--bedgraph");
            args.extend(&extra);
            let mut with_out = args.clone();
            with_out.extend(["-o", &bw]);
            gloci().args(&with_out).assert().success();
            assert_eq!(
                out(&["view", &bw]),
                sorted_genomecov(&args[1..]),
                "{args:?}"
            );
        }
    }
}

#[test]
fn genomecov_bigwig_needs_a_bedgraph_mode() {
    let d = tempdir();
    let (bed, genome, _) = genomecov_inputs();
    let x = |n: &str| s(&d.path().join(n));
    let base = ["genomecov", bed.as_str(), "-g", genome.as_str()];
    let run = |extra: &[&str], needles: &[&str]| {
        let mut a: Vec<&str> = base.to_vec();
        a.extend(extra);
        fail_with(2, &a, needles);
    };
    run(&["-o", &x("x.bw")], &["--bedgraph", "--bga"]);
    run(&["--per-base", "-o", &x("x.bw")], &["--bedgraph", "--bga"]);
    run(&["--bga", "-o", &x("x.bb")], &["bigBed"]);
    run(
        &["--bga", "-o", &x("x.bw"), "--write-index"],
        &["--write-index"],
    );
    assert!(names(d.path()).is_empty());
}

#[test]
fn genomecov_bigwig_threads() {
    let d = tempdir();
    let bam = s(&diff_inputs().join("aln/cigars-bases.bam"));
    let at = |t: &str| {
        let bw = d.path().join(format!("t{t}.bw"));
        gloci()
            .args(["genomecov", &bam, "--bga", "-t", t, "-o", &s(&bw)])
            .assert()
            .success();
        std::fs::read(bw).unwrap()
    };
    assert_eq!(at("1"), at("8"));
}

#[test]
fn genomecov_empty_bedgraph_to_bigwig_is_an_error() {
    let d = tempdir();
    let (_, genome, _) = genomecov_inputs();
    let empty = plain_file(d.path(), "empty.bed", "");
    let bw = s(&d.path().join("x.bw"));
    fail_with(
        1,
        &[
            "genomecov",
            &s(&empty),
            "-g",
            &genome,
            "--bedgraph",
            "-o",
            &bw,
        ],
        &["no records"],
    );
    assert_eq!(names(d.path()), vec!["empty.bed"]);
}

#[test]
fn genomecov_refused_bigwig_record_is_a_data_error() {
    let d = tempdir();
    let (bed, genome, _) = genomecov_inputs();
    let bw = s(&d.path().join("x.bw"));
    let a = gloci()
        .args([
            "genomecov",
            &bed,
            "-g",
            &genome,
            "--bedgraph",
            "--scale",
            "1e39",
            "-o",
            &bw,
        ])
        .assert()
        .code(1)
        .stderr(contains("record 1 of the output"))
        .stderr(contains("32-bit float"));
    let err = String::from_utf8_lossy(&a.get_output().stderr).into_owned();
    assert!(!err.contains("free disk space"), "{err}");
    assert!(!Path::new(&bw).exists());
}

#[test]
fn map_over_bigwig_gives_bigwigaverageoverbed_columns() {
    let d = tempdir();
    let bw = d.path().join("r.bw");
    bigtools_bigwig(
        &bw,
        &[("chr1", 1000)],
        &[("chr1", 0, 100, 1.5), ("chr1", 100, 200, 0.25)],
    );
    let a = d.path().join("a.bed");
    std::fs::write(&a, "chr1\t50\t150\tx\nchr1\t500\t600\ty\n").unwrap();
    gloci()
        .args(["map", "-a", &s(&a), "-b", &s(&bw)])
        .args(["-c", "4:wsum", "-c", "0:bases_uniq", "-c", "4:wmean"])
        .args(["--empty", "0"])
        .assert()
        .success()
        .stdout("chr1\t50\t150\tx\t87.5\t100\t0.875\nchr1\t500\t600\ty\t0\t0\t0\n");
}

// ---- the writer's temporary files: the open-file limit and $TMPDIR (D171) ----

/// `text`'s rows sorted by chromosome then start, as `view` lists a bigWig.
fn lex_sorted(text: &str) -> String {
    let mut rows: Vec<(&str, u64, &str)> = text
        .lines()
        .map(|l| {
            let mut f = l.splitn(3, '\t');
            let chrom = f.next().unwrap();
            (chrom, f.next().unwrap().parse().unwrap(), l)
        })
        .collect();
    rows.sort_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));
    rows.iter().map(|r| format!("{}\n", r.2)).collect()
}

/// A genome of one 100 Mb chromosome, which gives every chromosome 10 zoom levels and so about
/// 11 temporary files in bigtools, and `n` scaffolds of 1 kb; a bedGraph and a BED with two
/// records on each chromosome. Returns the genome, the bedGraph and the BED.
fn scaffolds(dir: &Path, n: usize) -> (String, String, String) {
    let mut genome = String::from("big\t100000000\n");
    let mut bg = String::from("big\t0\t10\t1\nbig\t20\t30\t2\n");
    for i in 0..n {
        genome.push_str(&format!("s{i}\t1000\n"));
        bg.push_str(&format!("s{i}\t0\t10\t{i}\ns{i}\t500\t600\t0.5\n"));
    }
    let bed: String = bg
        .lines()
        .map(|l| format!("{}\n", l.rsplit_once('\t').unwrap().0))
        .collect();
    (
        s(&plain_file(dir, "g.genome", &genome)),
        s(&plain_file(dir, "in.bedgraph", &bg)),
        s(&plain_file(dir, "in.bed", &bed)),
    )
}

/// `gloci ARGS` run by `sh` after `limits` (`ulimit` commands); `None` when the shell cannot
/// set them here, which the caller reports as a skip.
#[cfg(unix)]
fn under_limits(limits: &str, args: &[&str]) -> Option<std::process::Output> {
    let out = std::process::Command::new("sh")
        .arg("-c")
        .arg(format!("{limits} || exit 77; exec \"$0\" \"$@\""))
        .arg(assert_cmd::cargo::cargo_bin("gloci"))
        .args(args)
        .output()
        .unwrap();
    if out.status.code() == Some(77) {
        eprintln!("skipped: sh cannot run `{limits}` here");
        return None;
    }
    Some(out)
}

/// The temporary files left in `dir`.
fn temp_files_in(dir: &Path) -> Vec<String> {
    names(dir)
        .into_iter()
        .filter(|n| n.contains(".tmp-"))
        .collect()
}

/// With the hard limit at 256 gloci cannot raise it, and bigtools needs about 4,400 files: each
/// run either writes a file that reads back whole or exits 1 naming the limit, never exit 0
/// with a bad file (as bigtools alone gives), and leaves no temporary file.
#[cfg(unix)]
#[test]
fn a_hard_open_file_limit_never_gives_a_corrupt_file() {
    let d = tempdir();
    let (g, bg, bed) = scaffolds(d.path(), 400);
    let want_bg = lex_sorted(&std::fs::read_to_string(&bg).unwrap());
    let want_bed = lex_sorted(&std::fs::read_to_string(&bed).unwrap());
    let want_cov = sorted_genomecov(&[&bed, "-g", &g, "--bga"]);
    let runs: [(&str, Vec<&str>, &String, &str); 4] = [
        (
            "convert -t 1",
            vec!["convert", &bg, "-g", &g, "-t", "1"],
            &want_bg,
            "bw",
        ),
        (
            "convert -t 4",
            vec!["convert", &bg, "-g", &g, "-t", "4"],
            &want_bg,
            "bw",
        ),
        (
            "convert to bigBed",
            vec!["convert", &bed, "-g", &g],
            &want_bed,
            "bb",
        ),
        (
            "genomecov --bga",
            vec!["genomecov", &bed, "-g", &g, "--bga"],
            &want_cov,
            "bw",
        ),
    ];
    for (i, (name, args, want, ext)) in runs.into_iter().enumerate() {
        let bw = d.path().join(format!("r{i}.{ext}"));
        let mut args = args;
        let bw_s = s(&bw);
        args.extend(["-o", &bw_s]);
        let Some(o) = under_limits("ulimit -n 256", &args) else {
            return;
        };
        let err = String::from_utf8_lossy(&o.stderr);
        match o.status.code() {
            Some(0) => assert_eq!(&out(&["view", &bw_s]), want, "{name}"),
            Some(1) => {
                assert!(
                    err.contains("ulimit -n: soft 256, hard 256"),
                    "{name}: {err}"
                );
                assert!(err.contains("$TMPDIR"), "{name}: {err}");
                assert!(!err.contains("internal error"), "{name}: {err}");
                assert!(!bw.exists(), "{name}");
            }
            other => panic!("{name}: exit {other:?}: {err}"),
        }
        assert!(temp_files_in(d.path()).is_empty(), "{name}");
    }
}

/// A soft limit of 256 under a hard one of 4,096: gloci raises the soft limit, and the file,
/// which needs about 1,100 open files, is written whole.
#[cfg(unix)]
#[test]
fn a_low_soft_open_file_limit_is_raised() {
    let d = tempdir();
    let (g, bg, _) = scaffolds(d.path(), 100);
    let want = lex_sorted(&std::fs::read_to_string(&bg).unwrap());
    for t in ["1", "4"] {
        let bw = s(&d.path().join(format!("t{t}.bw")));
        let Some(o) = under_limits(
            "ulimit -S -n 256 && ulimit -H -n 4096",
            &["-vv", "convert", &bg, "-g", &g, "-t", t, "-o", &bw],
        ) else {
            return;
        };
        let err = String::from_utf8_lossy(&o.stderr);
        assert!(o.status.success(), "-t {t}: {err}");
        assert!(
            err.contains("raised the open-file limit from 256 to 4096"),
            "-t {t}: {err}"
        );
        assert_eq!(out(&["view", &bw]), want, "-t {t}");
    }
    assert!(temp_files_in(d.path()).is_empty());
}

/// A `$TMPDIR` that does not exist is exit 1 naming it, before any output is created.
#[cfg(unix)]
#[test]
fn a_missing_tmpdir_is_an_error_naming_it() {
    let d = tempdir();
    let (g, bg, bed) = scaffolds(d.path(), 3);
    let missing = s(&d.path().join("no-such-dir"));
    let bw = s(&d.path().join("x.bw"));
    for args in [
        vec!["convert", &bg, "-g", &g, "-o", &bw],
        vec!["genomecov", &bed, "-g", &g, "--bga", "-o", &bw],
    ] {
        let a = gloci()
            .env("TMPDIR", &missing)
            .args(&args)
            .assert()
            .code(1)
            .stderr(contains(format!("$TMPDIR ({missing})")))
            .stderr(contains("set TMPDIR"));
        let err = String::from_utf8_lossy(&a.get_output().stderr).into_owned();
        assert!(!err.contains("internal error"), "{err}");
        assert_eq!(names(d.path()), vec!["g.genome", "in.bed", "in.bedgraph"]);
    }
}
