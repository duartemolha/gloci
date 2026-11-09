//! `gloci merge`: signed gaps, points, strands, aggregation, k-way inputs,
//! VCF footprints, order errors and determinism.
mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;
use std::path::Path;

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn out(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

fn bed(rows: &[(u64, u64)]) -> String {
    rows.iter()
        .map(|(a, b)| format!("chr1\t{a}\t{b}\n"))
        .collect()
}

#[test]
fn gap_table() {
    let dir = tempdir();
    type Case = (&'static [(u64, u64)], &'static str, &'static str);
    let cases: [Case; 5] = [
        (&[(1, 5), (5, 8)], "0", "chr1\t1\t8\n"),
        (&[(1, 5), (6, 8)], "0", "chr1\t1\t5\nchr1\t6\t8\n"),
        (&[(1, 5), (6, 8)], "1", "chr1\t1\t8\n"),
        (&[(1, 5), (4, 8)], "-2", "chr1\t1\t5\nchr1\t4\t8\n"),
        (&[(1, 5), (3, 8)], "-2", "chr1\t1\t8\n"),
    ];
    for (i, (rows, d, want)) in cases.iter().enumerate() {
        let f = s(&plain_file(dir.path(), &format!("g{i}.bed"), &bed(rows)));
        assert_eq!(out(&["merge", &f, "-d", d]), *want, "case {i}");
        assert_eq!(out(&["merge", &f, &format!("--max-gap={d}")]), *want);
    }
}

#[test]
fn points() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "p.bed",
        "chr1\t1\t8\nchr1\t5\t5\nchr1\t9\t9\nchr1\t9\t9\n",
    ));
    assert_eq!(out(&["merge", &f, "-d", "0"]), "chr1\t1\t8\nchr1\t9\t9\n");
}

#[test]
fn per_strand() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "s.bed",
        "chr1\t1\t10\ta\t0\t+\nchr1\t5\t12\tb\t0\t-\nchr1\t8\t15\tc\t0\t+\nchr1\t9\t11\td\t0\t.\n",
    ));
    assert_eq!(out(&["merge", &f]), "chr1\t1\t15\n");
    assert_eq!(
        out(&["merge", &f, "-s"]),
        "chr1\t1\t15\nchr1\t5\t12\nchr1\t9\t11\n"
    );
    assert_eq!(
        out(&["merge", &f, "-s", "--count"]),
        "chr1\t1\t15\t2\nchr1\t5\t12\t1\nchr1\t9\t11\t1\n"
    );
}

#[test]
fn count_and_aggs() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "c.bed",
        "chr1\t1\t10\tx\t3\nchr1\t5\t12\ty\t7\nchr1\t5\t13\tx\t.\nchr1\t40\t50\tz\t2\nchr2\t1\t2\tq\t9\n",
    ));
    assert_eq!(
        out(&["merge", &f, "--count", "-c", "5:max", "-c", "4:distinct"]),
        "chr1\t1\t13\t3\t7\tx,y\nchr1\t40\t50\t1\t2\tz\nchr2\t1\t2\t1\t9\tq\n"
    );
    assert_eq!(
        out(&["merge", &f, "-c", "5:sum", "--delim", ";"]),
        "chr1\t1\t13\t10\nchr1\t40\t50\t2\nchr2\t1\t2\t9\n"
    );
    gloci()
        .args(["merge", &f, "-c", "5:bases"])
        .assert()
        .code(2)
        .stderr(contains("bases"));
    gloci()
        .args(["merge", &f, "-c", "0:bases"])
        .assert()
        .code(2);
}

#[test]
fn two_input_ops_are_refused() {
    let dir = tempdir();
    let f = s(&plain_file(dir.path(), "a.bed", "chr1\t1\t5\t3\n"));
    for op in [
        "wmean",
        "records",
        "min_record",
        "max_record",
        "frac_a",
        "frac_b",
        "bases_uniq",
    ] {
        gloci()
            .args(["merge", &f, "-c", &format!("4:{op}")])
            .assert()
            .code(2)
            .stderr(contains(format!("-c 4:{op}")))
            .stderr(contains("merge has none"));
    }
}

#[test]
fn single_input_ops_work_in_merge() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t1\t5\t1\nchr1\t3\t9\t2\nchr1\t4\t6\t2\nchr1\t5\t7\t10\n",
    ));
    gloci()
        .args(["merge", &f, "-c", "4:mode,4:mad,4:quantile:0.5,4:cv"])
        .assert()
        .success()
        .stdout("chr1\t1\t9\t2\t0.5\t2\t1.1181996111474777\n");
}

#[test]
fn multiple_inputs_kmerged() {
    let dir = tempdir();
    let a = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t1\t5\nchr1\t20\t30\nchr2\t1\t4\n",
    ));
    let b = s(&plain_file(
        dir.path(),
        "b.bed",
        "chr1\t4\t8\nchr1\t25\t26\nchr3\t1\t2\n",
    ));
    assert_eq!(
        out(&["merge", &a, &b, "--count"]),
        "chr1\t1\t8\t2\nchr1\t20\t30\t2\nchr2\t1\t4\t1\nchr3\t1\t2\t1\n"
    );
}

#[test]
fn vcf_spans_heap() {
    let dir = tempdir();
    let vcf = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t101\t.\tATTTT\tA\t.\t.\t.\n\
chr1\t102\t.\tA\tG\t.\t.\t.\n";
    let f = s(&plain_file(dir.path(), "v.vcf", vcf));
    assert_eq!(out(&["merge", &f]), "chr1\t101\t105\n");
}

#[test]
fn unsorted_input() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "u.bed",
        "chr1\t50\t60\nchr1\t1\t5\n",
    ));
    gloci()
        .args(["merge", &f])
        .assert()
        .code(1)
        .stderr(contains("gloci sort"));
    let g = s(&plain_file(dir.path(), "ok.bed", "chr1\t1\t5\n"));
    gloci().args(["merge", &g, &f]).assert().code(1);
    gloci().args(["merge", &f, "--unsorted"]).assert().code(2);
}

#[test]
fn tsv_output_and_empty_input() {
    let dir = tempdir();
    let f = s(&plain_file(dir.path(), "a.bed", "chr1\t1\t5\nchr1\t4\t9\n"));
    let o = dir.path().join("o.tsv");
    gloci().args(["merge", &f, "-o", &s(&o)]).assert().success();
    assert_eq!(std::fs::read_to_string(&o).unwrap(), "chr1\t1\t9\n");
    let e = s(&plain_file(dir.path(), "e.bed", ""));
    assert_eq!(out(&["merge", &e]), "");
}

#[test]
fn determinism() {
    let dir = tempdir();
    let a = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t1\t5\ta\t1\t+\nchr1\t4\t9\tb\t2\t-\nchr2\t1\t4\tc\t3\t+\n",
    ));
    let b = s(&plain_file(dir.path(), "b.bed", "chr1\t8\t12\td\t4\t+\n"));
    let run = |t: &str| out(&["merge", &a, &b, "-s", "--count", "-c", "5:mean", "-t", t]);
    let base = run("1");
    assert_eq!(run("2"), base);
    assert_eq!(run("8"), base);
}

#[test]
fn points_at_cluster_end() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "pe.bed",
        "chr1\t1\t8\nchr1\t8\t8\n",
    ));
    assert_eq!(out(&["merge", &f, "-d", "0"]), "chr1\t1\t8\n");
    assert_eq!(out(&["merge", &f, "-d", "-1"]), "chr1\t1\t8\nchr1\t8\t8\n");
}

const VCF_HEAD: &str = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

#[test]
fn vcf_heap_orders_by_span_start() {
    // The deletion has key 100 and span [101,105); the SNP after it has the
    // same key and span [100,101), so it must be taken first: the cluster
    // starts at 100, not 101.
    let dir = tempdir();
    let vcf = format!("{VCF_HEAD}chr1\t101\t.\tATTTT\tA\t.\t.\t.\nchr1\t101\t.\tA\tG\t.\t.\t.\n");
    let f = s(&plain_file(dir.path(), "h.vcf", &vcf));
    assert_eq!(out(&["merge", &f, "--count"]), "chr1\t100\t105\t2\n");
}

#[test]
fn vcf_and_bed_kmerged() {
    let dir = tempdir();
    let vcf = format!("{VCF_HEAD}chr1\t101\t.\tATTTT\tA\t.\t.\t.\nchr1\t300\t.\tA\tG\t.\t.\t.\n");
    let v = s(&plain_file(dir.path(), "k.vcf", &vcf));
    let b = s(&plain_file(
        dir.path(),
        "k.bed",
        "chr1\t104\t110\nchr1\t299\t300\nchr1\t500\t510\n",
    ));
    assert_eq!(
        out(&["merge", &v, &b, "--count"]),
        "chr1\t101\t110\t2\nchr1\t299\t300\t2\nchr1\t500\t510\t1\n"
    );
}

#[test]
fn c_column_past_the_end_of_a_line_is_an_error() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "short.bed",
        "chr1\t0\t100\tx\t5\nchr1\t50\t150\ty\n",
    ));
    for op in ["5:sum", "5:count"] {
        gloci()
            .args(["merge", &f, "-c", op])
            .assert()
            .code(1)
            .stderr(contains("line 2").and(contains("column 5")));
    }
    // column 0 is the whole record, an empty value is still just missing
    assert_eq!(out(&["merge", &f, "-c", "0:count"]), "chr1\t0\t150\t2\n");
}
