//! `gloci convert`: VCF/BCF, GFF3/GTF and BED12 to BED.
mod common;
use common::*;
use predicates::str::contains;
use std::path::Path;

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn out(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

const VCF_HEAD: &str = "##fileformat=VCFv4.2\n##contig=<ID=c1,length=100000>\n\
##contig=<ID=c2,length=100000>\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

fn vcf(dir: &Path, name: &str, rows: &str) -> String {
    s(&plain_file(dir, name, &format!("{VCF_HEAD}{rows}")))
}

/// The spec 4.2 table: one record per row, each on its own chromosome so the
/// file is sorted.
const TABLE: &str = "c1\t100\tr1\tA\tG\t10\tPASS\t.\n\
c1\t200\tr2\tAC\tGT\t.\t.\t.\n\
c1\t300\tr3\tACGT\tA\t.\t.\t.\n\
c1\t400\tr4\tA\tACGT\t.\t.\t.\n\
c1\t500\tr5\tATG\tAC\t.\t.\t.\n\
c1\t600\tr6\tACGT\tA,AC\t.\t.\t.\n\
c1\t700\tr7\tN\t<DEL>\t.\t.\tEND=1000\n\
c1\t800\tr8\tN\t<INS>\t.\t.\t.\n\
c1\t900\tr9\tN\t<DUP>\t.\t.\tEND=1200\n\
c1\t1000\tr10\tG\tG]17:198982]\t.\t.\t.\n";

#[test]
fn vcf_rows_allele_span() {
    let dir = tempdir();
    let f = vcf(dir.path(), "a.vcf", TABLE);
    assert_eq!(
        out(&["convert", &f]),
        "c1\t99\t100\tr1\tA\tG\t10\tPASS\n\
c1\t199\t201\tr2\tAC\tGT\t.\t.\n\
c1\t300\t303\tr3\tACGT\tA\t.\t.\n\
c1\t400\t400\tr4\tA\tACGT\t.\t.\n\
c1\t500\t502\tr5\tATG\tAC\t.\t.\n\
c1\t600\t603\tr6\tACGT\tA,AC\t.\t.\n\
c1\t700\t1000\tr7\tN\t<DEL>\t.\t.\n\
c1\t800\t800\tr8\tN\t<INS>\t.\t.\n\
c1\t900\t1200\tr9\tN\t<DUP>\t.\t.\n\
c1\t999\t1000\tr10\tG\tG]17:198982]\t.\t.\n"
    );
}

#[test]
fn vcf_rows_ref_span() {
    let dir = tempdir();
    let f = vcf(dir.path(), "a.vcf", TABLE);
    assert_eq!(
        out(&["convert", &f, "--vcf-span", "ref"]),
        "c1\t99\t100\tr1\tA\tG\t10\tPASS\n\
c1\t199\t201\tr2\tAC\tGT\t.\t.\n\
c1\t299\t303\tr3\tACGT\tA\t.\t.\n\
c1\t399\t400\tr4\tA\tACGT\t.\t.\n\
c1\t499\t502\tr5\tATG\tAC\t.\t.\n\
c1\t599\t603\tr6\tACGT\tA,AC\t.\t.\n\
c1\t699\t1000\tr7\tN\t<DEL>\t.\t.\n\
c1\t799\t800\tr8\tN\t<INS>\t.\t.\n\
c1\t899\t1200\tr9\tN\t<DUP>\t.\t.\n\
c1\t999\t1000\tr10\tG\tG]17:198982]\t.\t.\n"
    );
}

#[test]
fn split_alleles() {
    let dir = tempdir();
    let f = vcf(dir.path(), "a.vcf", "c1\t100\t.\tACGT\tA,AC\t.\t.\t.\n");
    assert_eq!(
        out(&["convert", &f, "--split-alleles"]),
        "c1\t100\t103\t.\tACGT\tA\t.\t.\nc1\t101\t103\t.\tACGT\tAC\t.\t.\n"
    );
}

#[test]
fn split_alleles_skips_star_and_gvcf_alleles() {
    let dir = tempdir();
    let f = vcf(
        dir.path(),
        "a.vcf",
        "c1\t100\t.\tA\tG,*\t.\t.\t.\nc1\t200\t.\tA\tT,<NON_REF>\t.\t.\t.\n\
c1\t300\t.\tA\t<NON_REF>\t.\t.\t.\n",
    );
    // A record with no ALT that has a footprint keeps one row, whole.
    assert_eq!(
        out(&["convert", &f, "--split-alleles"]),
        "c1\t99\t100\t.\tA\tG\t.\t.\nc1\t199\t200\t.\tA\tT\t.\t.\nc1\t299\t300\t.\tA\t<NON_REF>\t.\t.\n"
    );
}

#[test]
fn split_alleles_with_ref_span_gives_every_row_the_ref_footprint() {
    let dir = tempdir();
    let f = vcf(
        dir.path(),
        "a.vcf",
        "c1\t100\t.\tACGT\tA,AC\t.\t.\t.\nc1\t200\t.\tA\tG,*\t.\t.\t.\n\
c1\t300\t.\tN\t<DEL>,<INS>\t.\t.\tEND=500\n",
    );
    // `*` is still skipped; a symbolic ALT runs to END under ref.
    assert_eq!(
        out(&["convert", &f, "--split-alleles", "--vcf-span", "ref"]),
        "c1\t99\t103\t.\tACGT\tA\t.\t.\nc1\t99\t103\t.\tACGT\tAC\t.\t.\n\
c1\t199\t200\t.\tA\tG\t.\t.\nc1\t299\t500\t.\tN\t<DEL>\t.\t.\nc1\t299\t500\t.\tN\t<INS>\t.\t.\n"
    );
}

#[test]
fn types_del() {
    let dir = tempdir();
    let f = vcf(
        dir.path(),
        "a.vcf",
        "c1\t100\ta\tA\tG\t.\t.\t.\nc1\t200\tb\tAT\tA\t.\t.\t.\nc1\t300\tc\tA\tAT\t.\t.\t.\n\
c1\t400\td\tAT\tA,G\t.\t.\t.\nc1\t500\te\tA\tT,AG\t.\t.\t.\n",
    );
    // A record whose type set meets the list is kept whole (D50) ...
    assert_eq!(
        out(&["convert", &f, "--types", "del"]),
        "c1\t200\t201\tb\tAT\tA\t.\t.\nc1\t399\t401\td\tAT\tA,G\t.\t.\n"
    );
    // ... and with --split-alleles each ALT is judged by itself.
    assert_eq!(
        out(&["convert", &f, "--types", "del", "--split-alleles"]),
        "c1\t200\t201\tb\tAT\tA\t.\t.\nc1\t400\t401\td\tAT\tA\t.\t.\n"
    );
    assert_eq!(
        out(&["convert", &f, "--types", "ins,snv"]),
        "c1\t99\t100\ta\tA\tG\t.\t.\nc1\t300\t300\tc\tA\tAT\t.\t.\nc1\t499\t500\te\tA\tT,AG\t.\t.\n"
    );
}

#[test]
fn bcf_input_matches_vcf_input() {
    let dir = tempdir();
    let f = vcf(dir.path(), "a.vcf", TABLE);
    let b = s(&dir.path().join("a.bcf"));
    gloci().args(["view", &f, "-o", &b]).assert().success();
    for extra in [vec![], vec!["--split-alleles"], vec!["--vcf-span", "ref"]] {
        let mut v = vec!["convert", &f];
        v.extend(&extra);
        let mut bb = vec!["convert", &b];
        bb.extend(&extra);
        assert_eq!(out(&v), out(&bb), "{extra:?}");
    }
}

const GFF3: &str = "##gff-version 3\n\
c1\tsrc\tgene\t11\t20\t.\t+\t.\tID=g1;Name=Foo\n\
c1\tsrc\texon\t11\t15\t5.5\t-\t1\tParent=g1;Name=Bar\n\
c1\thav\tCDS\t16\t20\t.\t.\t2\tID=cds1;Parent=g1\n";

#[test]
fn gff3_id_and_name_attr() {
    let dir = tempdir();
    let f = s(&plain_file(dir.path(), "a.gff3", GFF3));
    assert_eq!(
        out(&["convert", &f]),
        "c1\t10\t20\tg1\t.\t+\tsrc\tgene\t.\tID=g1;Name=Foo\n\
c1\t10\t15\t.\t5.5\t-\tsrc\texon\t1\tParent=g1;Name=Bar\n\
c1\t15\t20\tcds1\t.\t.\thav\tCDS\t2\tID=cds1;Parent=g1\n"
    );
    assert_eq!(
        out(&["convert", &f, "--name-attr", "Name"]),
        "c1\t10\t20\tFoo\t.\t+\tsrc\tgene\t.\tID=g1;Name=Foo\n\
c1\t10\t15\tBar\t5.5\t-\tsrc\texon\t1\tParent=g1;Name=Bar\n\
c1\t15\t20\t.\t.\t.\thav\tCDS\t2\tID=cds1;Parent=g1\n"
    );
}

const GTF: &str = "c1\tsrc\texon\t100\t300\t5.5\t+\t1\tgene_id \"g1\"; transcript_id \"t1\"; gene_name \"Foo;Bar\";\n\
c1\tsrc\tgene\t100\t900\t.\t-\t.\tgene_id \"g2\";\n\
c2\tsrc\tgene\t10\t20\t.\t.\t.\tfoo \"bar\";\n";

#[test]
fn gtf_gene_id_default() {
    let dir = tempdir();
    let f = s(&plain_file(dir.path(), "a.gtf", GTF));
    let all = out(&["convert", &f]);
    assert_eq!(
        all,
        "c1\t99\t300\tg1\t5.5\t+\tsrc\texon\t1\tgene_id \"g1\"; transcript_id \"t1\"; gene_name \"Foo;Bar\";\n\
c1\t99\t900\tg2\t.\t-\tsrc\tgene\t.\tgene_id \"g2\";\n\
c2\t9\t20\t.\t.\t.\tsrc\tgene\t.\tfoo \"bar\";\n"
    );
    let names: Vec<_> = out(&["convert", &f, "--name-attr", "gene_name"])
        .lines()
        .map(|l| l.split('\t').nth(3).unwrap().to_string())
        .collect();
    // A quoted value keeps its semicolon.
    assert_eq!(names, ["Foo;Bar", ".", "."]);
    let names: Vec<_> = out(&["convert", &f, "--name-attr", "transcript_id"])
        .lines()
        .map(|l| l.split('\t').nth(3).unwrap().to_string())
        .collect();
    assert_eq!(names, ["t1", ".", "."]);
}

#[test]
fn feature_filter() {
    let dir = tempdir();
    let f = s(&plain_file(dir.path(), "a.gff3", GFF3));
    let types = |args: &[&str]| -> Vec<String> {
        let mut v = vec!["convert", &f];
        v.extend(args);
        out(&v)
            .lines()
            .map(|l| l.split('\t').nth(7).unwrap().to_string())
            .collect()
    };
    assert_eq!(types(&["--feature", "exon"]), ["exon"]);
    assert_eq!(
        types(&["--feature", "exon", "--feature", "CDS"]),
        ["exon", "CDS"]
    );
    assert_eq!(types(&["--feature", "gene,CDS"]), ["gene", "CDS"]);
    assert!(types(&["--feature", "nothing"]).is_empty());
}

#[test]
fn bed12_blocks() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "a.bed",
        "c1\t100\t200\tn1\t5\t+\t100\t200\t0\t2\t10,20,\t0,80,\n\
c1\t300\t340\tn2\t.\t-\t300\t340\t0\t1\t40\t0\n",
    ));
    assert_eq!(
        out(&["convert", &f, "--blocks"]),
        "c1\t100\t110\tn1\t5\t+\nc1\t180\t200\tn1\t5\t+\nc1\t300\t340\tn2\t.\t-\n"
    );
}

#[test]
fn usage_errors() {
    let dir = tempdir();
    let v = vcf(dir.path(), "a.vcf", "c1\t100\t.\tA\tG\t.\t.\t.\n");
    let g = s(&plain_file(dir.path(), "a.gff3", GFF3));
    let b3 = s(&plain_file(dir.path(), "a.bed", "c1\t1\t5\n"));
    let b12 = s(&plain_file(
        dir.path(),
        "b.bed",
        "c1\t100\t110\tn1\t5\t+\t100\t110\t0\t1\t10\t0\n",
    ));
    let fail = |args: &[&str]| gloci().args(args).assert().code(2);
    fail(&["convert", &b3]).stderr(contains("--blocks"));
    fail(&["convert", &b12]).stderr(contains("--blocks"));
    fail(&["convert", &v, "--blocks"]).stderr(contains("BED12"));
    fail(&["convert", &g, "--blocks"]).stderr(contains("BED12"));
    fail(&["convert", &v, "--feature", "exon"]).stderr(contains("GFF3/GTF"));
    fail(&["convert", &v, "--name-attr", "x"]).stderr(contains("GFF3/GTF"));
    fail(&["convert", &g, "--split-alleles"]).stderr(contains("VCF/BCF"));
    fail(&["convert", &g, "--types", "snv"]).stderr(contains("VCF/BCF"));
    fail(&["convert", &b12, "--blocks", "--types", "snv"]).stderr(contains("VCF/BCF"));
    fail(&["convert", &v, "--types", "bogus"]);
    // The format is not BED12.
    gloci()
        .args(["convert", &b3, "--blocks"])
        .assert()
        .code(1)
        .stderr(contains("12 columns"));
    // An empty attribute key is useless.
    fail(&["convert", &g, "--name-attr", ""]);
}

const BREAKS_ORDER: &str = "c1\t100\t.\tACGT\tA,AC\t.\t.\t.\nc1\t101\t.\tG\tT\t.\t.\t.\n";

#[test]
fn split_alleles_that_break_order_warn_once() {
    let dir = tempdir();
    let f = vcf(
        dir.path(),
        "a.vcf",
        &format!(
            "{BREAKS_ORDER}c1\t102\t.\tG\tT\t.\t.\t.\nc1\t103\t.\tACGT\tA,AC\t.\t.\t.\n\
c1\t104\t.\tG\tT\t.\t.\t.\n"
        ),
    );
    let o = gloci()
        .args(["convert", &f, "--split-alleles"])
        .assert()
        .success();
    let err = String::from_utf8(o.get_output().stderr.clone()).unwrap();
    assert_eq!(err.matches("out of order").count(), 1, "{err}");
    assert!(err.contains("output row 3"), "{err}");
    assert!(err.contains("gloci sort"), "{err}");
}

#[test]
fn sorted_output_and_unsorted_input_do_not_warn() {
    let dir = tempdir();
    let f = vcf(
        dir.path(),
        "a.vcf",
        "c1\t100\t.\tACGT\tA,AC\t.\t.\t.\nc1\t200\t.\tG\tT\t.\t.\t.\n",
    );
    gloci()
        .args(["convert", &f, "--split-alleles"])
        .assert()
        .success()
        .stderr("");
    // Input that is itself out of order is not the conversion's doing.
    let u = vcf(
        dir.path(),
        "u.vcf",
        "c1\t200\t.\tG\tT\t.\t.\t.\nc1\t100\t.\tG\tT\t.\t.\t.\n",
    );
    gloci().args(["convert", &u]).assert().success().stderr("");
}

#[test]
fn write_index_on_out_of_order_output_names_the_output() {
    let dir = tempdir();
    let f = vcf(dir.path(), "a.vcf", BREAKS_ORDER);
    let o = s(&dir.path().join("o.bed.gz"));
    gloci()
        .args(["convert", &f, "--split-alleles", "-o", &o, "--write-index"])
        .assert()
        .code(1)
        .stderr(contains("o.bed.gz"))
        .stderr(contains("not sorted"));
    assert!(!dir.path().join("o.bed.gz").exists());
}

#[test]
fn bgzf_output_with_index() {
    let dir = tempdir();
    let f = vcf(dir.path(), "a.vcf", TABLE);
    let o = s(&dir.path().join("o.bed.gz"));
    gloci()
        .args(["convert", &f, "-o", &o, "--write-index"])
        .assert()
        .success();
    assert!(dir.path().join("o.bed.gz.tbi").exists());
}

#[test]
fn reads_stdin() {
    let o = gloci()
        .args(["convert", "-", "--format", "vcf"])
        .write_stdin(format!("{VCF_HEAD}c1\t100\t.\tA\tG\t.\t.\t.\n"))
        .assert()
        .success();
    assert_eq!(
        String::from_utf8(o.get_output().stdout.clone()).unwrap(),
        "c1\t99\t100\t.\tA\tG\t.\t.\n"
    );
}

#[test]
fn determinism() {
    let dir = tempdir();
    let f = vcf(dir.path(), "a.vcf", TABLE);
    let g = s(&plain_file(dir.path(), "a.gff3", GFF3));
    for (file, extra) in [(&f, "--split-alleles"), (&g, "--name-attr")] {
        let run = |t: &str| {
            let mut a = vec!["convert", file.as_str()];
            if extra == "--name-attr" {
                a.extend(["--name-attr", "Name"]);
            } else {
                a.push(extra);
            }
            a.extend(["-t", t]);
            out(&a)
        };
        let one = run("1");
        assert_eq!(one, run("2"));
        assert_eq!(one, run("8"));
    }
}

#[test]
fn gff3_values_are_not_url_decoded() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "a.gff3",
        "##gff-version 3\nc1\tsrc\tgene\t1\t10\t.\t+\t.\tID=a%3Bb%2Cc\n",
    ));
    assert_eq!(
        out(&["convert", &f]),
        "c1\t0\t10\ta%3Bb%2Cc\t.\t+\tsrc\tgene\t.\tID=a%3Bb%2Cc\n"
    );
}

#[test]
fn invalid_bed12_is_a_data_error_even_with_skip_invalid() {
    let dir = tempdir();
    let bad = |name: &str, tail: &str| {
        s(&plain_file(
            dir.path(),
            name,
            &format!("c1\t100\t200\tn\t5\t+\t100\t200\t0\t{tail}\n"),
        ))
    };
    for (name, tail, msg) in [
        ("zero.bed", "0\t\t", "blockCount is 0"),
        ("count.bed", "2\t10,\t0,", "blockCount is 2"),
        ("first.bed", "1\t100\t5", "first block"),
        ("last.bed", "1\t50\t0", "last block ends at 150"),
        ("overlap.bed", "2\t60,60\t0,40", "overlap"),
    ] {
        let f = bad(name, tail);
        for extra in [vec![], vec!["--skip-invalid"]] {
            let mut a = vec!["convert", &f, "--blocks"];
            a.extend(&extra);
            gloci()
                .args(&a)
                .assert()
                .code(1)
                .stderr(contains(name))
                .stderr(contains("line 1"))
                .stderr(contains(msg));
        }
    }
}
