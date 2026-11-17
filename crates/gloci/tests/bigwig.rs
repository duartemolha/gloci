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
