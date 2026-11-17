//! `gloci convert` on PSL and RepeatMasker `.out` (D114).
mod common;
use common::*;
use predicates::str::contains;

const HEADER: &str = "psLayout version 3\n\nmatch\tmis- \trep. \tN's\tQ gap\tQ gap\tT gap\tT gap\t\
strand\tQ        \tQ   \tQ    \tQ  \tT        \tT   \tT    \tT  \tblock\tblat\tblock\tblock\n\
match\tmatch\t   \tcount\tbases\tcount\tbases\t      \tname     \tsize\tstart\tend\tname     \t\
size\tstart\tend\tcount\tsizes \tqStarts\ttStarts\n\
---------------------------------------------------------------------------------------------\n";

/// A plus-strand single block, then a minus-strand three-block hit.
const ROWS: &str = "50\t1\t0\t0\t0\t0\t0\t0\t+\tq1\t100\t0\t51\tchr2\t1000\t100\t151\t1\t51,\t0,\t100,\n\
90\t2\t1\t0\t1\t3\t2\t10\t-\tq2\t120\t5\t100\tchr1\t5000\t2000\t2102\t3\t30,20,40,\t5,40,70,\t\
2000,2040,2070,\n";

const PSL_BED: &str = "chr2\t100\t151\tq1\t50\t+\t100\t1\t0\t0\t0\t0\t0\t0\t0\t51\t1000\t1\t51,\t0,\t100,\n\
chr1\t2000\t2102\tq2\t90\t-\t120\t2\t1\t0\t1\t3\t2\t10\t5\t100\t5000\t3\t30,20,40,\t5,40,70,\t\
2000,2040,2070,\n";

const RMSK_HEAD: &str = "   SW  perc perc perc  query      position in query           matching       \
repeat              position in repeat\nscore  div. del. ins.  sequence    begin     end    (left)    \
repeat         class/family         begin  end (left)   ID\n\n";

const RMSK_ROWS: &str = "  463  13.2  0.0  1.5  chr1       10001   10468 (248471868) +  (CCCTAA)n      \
Simple_repeat            1    468    (0)    1\n\
  239  28.4  1.2  3.0  chr1        9000    9300 (248482036) C  AluSx          SINE/Alu              \
(12)    290      1    2 *\n";

const RMSK_BED: &str = "chr1\t10000\t10468\t(CCCTAA)n\t463\t+\t13.2\t0.0\t1.5\t(248471868)\t\
Simple_repeat\t1\t468\t(0)\t1\n\
chr1\t8999\t9300\tAluSx\t239\t-\t28.4\t1.2\t3.0\t(248482036)\tSINE/Alu\t(12)\t290\t1\t2\t*\n";

fn convert(path: &std::path::Path) -> assert_cmd::assert::Assert {
    gloci().arg("convert").arg(path).assert()
}

#[test]
fn psl_spans_strands_and_blocks_with_and_without_header() {
    let dir = tempdir();
    for (name, text) in [
        ("bare.psl", ROWS.to_string()),
        ("head.psl", format!("{HEADER}{ROWS}")),
    ] {
        let f = plain_file(dir.path(), name, &text);
        convert(&f).success().stdout(PSL_BED).stderr("");
    }
}

#[test]
fn psl_without_final_newline_keeps_the_last_hit() {
    let dir = tempdir();
    let f = plain_file(dir.path(), "x.psl", ROWS.trim_end());
    convert(&f).success().stdout(PSL_BED);
}

#[test]
fn psl_translated_strand_is_written_as_is() {
    let dir = tempdir();
    let row = "20\t0\t0\t0\t0\t0\t0\t0\t+-\tq3\t30\t0\t20\tchrM\t16000\t10\t30\t1\t20,\t0,\t10,\n";
    let f = plain_file(dir.path(), "t.psl", row);
    convert(&f)
        .success()
        .stdout(contains("chrM\t10\t30\tq3\t20\t+-\t30\t"));
}

#[test]
fn psl_malformed_line_names_the_line() {
    let dir = tempdir();
    let bad = format!("{HEADER}{ROWS}50\t1\t0\t0\n");
    let f = plain_file(dir.path(), "bad.psl", &bad);
    let line = HEADER.lines().count() + ROWS.lines().count() + 1;
    convert(&f)
        .failure()
        .code(1)
        .stderr(contains(format!("line {line}")))
        .stderr(contains("expected 21 tab-separated PSL columns"));
    let swapped = format!(
        "{ROWS}{}",
        ROWS.lines()
            .next()
            .unwrap()
            .replace("\t100\t151\t", "\t151\t100\t")
    );
    let f = plain_file(dir.path(), "swap.psl", &swapped);
    convert(&f).failure().code(1).stderr(contains("line 3"));
    let nonnum = format!(
        "{ROWS}{}",
        ROWS.lines().next().unwrap().replace("\t151\t", "\tx\t")
    );
    let f = plain_file(dir.path(), "nonnum.psl", &nonnum);
    convert(&f).failure().code(1).stderr(contains("tEnd"));
    // --skip-invalid skips the bad row and keeps the others.
    let f = plain_file(dir.path(), "skip.psl", &format!("{ROWS}50\t1\n"));
    gloci()
        .arg("convert")
        .arg(&f)
        .arg("--skip-invalid")
        .assert()
        .success()
        .stdout(PSL_BED)
        .stderr(contains("skipped 1 invalid"));
}

#[test]
fn psl_out_of_order_output_warns_once() {
    let dir = tempdir();
    let first = ROWS.lines().nth(1).unwrap();
    let later = first.replace("\t2000\t2102\t", "\t100\t200\t");
    let latest = first.replace("\t2000\t2102\t", "\t50\t60\t");
    let f = plain_file(
        dir.path(),
        "o.psl",
        &format!("{first}\n{later}\n{latest}\n"),
    );
    let o = convert(&f).success();
    let err = String::from_utf8(o.get_output().stderr.clone()).unwrap();
    assert_eq!(err.matches("out of order").count(), 1, "{err}");
}

#[test]
fn psl_sorted_input_does_not_warn() {
    let dir = tempdir();
    let sorted = ROWS.lines().rev().collect::<Vec<_>>().join("\n") + "\n";
    let f = plain_file(dir.path(), "s.psl", &sorted);
    convert(&f).success().stderr("");
}

#[test]
fn rmsk_spans_strands_and_parentheses_with_and_without_header() {
    let dir = tempdir();
    for (name, text) in [
        ("bare.out", RMSK_ROWS.to_string()),
        ("head.out", format!("{RMSK_HEAD}{RMSK_ROWS}")),
    ] {
        let f = plain_file(dir.path(), name, &text);
        convert(&f).success().stdout(RMSK_BED);
    }
}

#[test]
fn rmsk_malformed_line_names_the_line() {
    let dir = tempdir();
    let bad = format!("{RMSK_HEAD}{RMSK_ROWS}  12  1.0  0.0  0.0  chr1  5  9\n");
    let f = plain_file(dir.path(), "bad.out", &bad);
    let line = RMSK_HEAD.lines().count() + 3;
    convert(&f)
        .failure()
        .code(1)
        .stderr(contains(format!("line {line}")))
        .stderr(contains("RepeatMasker columns"));
    // Begin 0 is not a 1-based position.
    let zero = RMSK_ROWS.replace("10001", "0");
    let f = plain_file(dir.path(), "zero.out", &zero);
    convert(&f).failure().code(1).stderr(contains("line 1"));
}

#[test]
fn rmsk_out_of_order_output_warns_once() {
    let dir = tempdir();
    let f = plain_file(dir.path(), "r.out", RMSK_ROWS);
    convert(&f).success().stderr(contains("out of order"));
}

#[test]
fn chrom_map_renames_psl_and_rmsk() {
    let dir = tempdir();
    let map = plain_file(dir.path(), "m.tsv", "chr1\t1\n");
    let p = plain_file(dir.path(), "a.psl", ROWS);
    let o = gloci()
        .arg("convert")
        .arg(&p)
        .arg("--chrom-map")
        .arg(&map)
        .assert()
        .success();
    let out = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    assert!(out.contains("\n1\t2000\t2102\tq2\t"), "{out}");
    let r = plain_file(dir.path(), "a.out", RMSK_ROWS);
    let o = gloci()
        .arg("convert")
        .arg(&r)
        .arg("--chrom-map")
        .arg(&map)
        .assert()
        .success();
    let out = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    assert!(out.starts_with("1\t10000\t10468\t"), "{out}");
}

#[test]
fn only_convert_reads_psl_and_rmsk() {
    let dir = tempdir();
    for (name, text) in [("x.psl", ROWS), ("x.out", RMSK_ROWS)] {
        let f = plain_file(dir.path(), name, text);
        let p = f.to_str().unwrap();
        gloci()
            .args(["merge", p])
            .assert()
            .code(2)
            .stderr(contains(format!("gloci convert {p}")));
    }
}

#[test]
fn compressed_psl_is_read() {
    let dir = tempdir();
    let f = bgzf_file(dir.path(), "x.psl.gz", &format!("{HEADER}{ROWS}"));
    convert(&f).success().stdout(PSL_BED);
}

#[test]
fn psl_with_a_22nd_column_is_refused_and_skippable() {
    let dir = tempdir();
    let row = ROWS.lines().next().unwrap();
    let f = plain_file(dir.path(), "w.psl", &format!("{ROWS}{row}\textra\n"));
    convert(&f)
        .failure()
        .code(1)
        .stderr(contains("line 3"))
        .stderr(contains("found 22"));
    gloci()
        .arg("convert")
        .arg(&f)
        .arg("--skip-invalid")
        .assert()
        .success()
        .stdout(PSL_BED);
}

#[test]
fn psl_block_count_must_match_the_block_lists() {
    let dir = tempdir();
    let row = ROWS
        .lines()
        .nth(1)
        .unwrap()
        .replace("\t3\t30,20,40,", "\t2\t30,20,40,");
    let f = plain_file(dir.path(), "b.psl", &format!("{row}\n"));
    convert(&f)
        .failure()
        .code(1)
        .stderr(contains("line 1"))
        .stderr(contains("blockCount is 2 but blockSizes lists 3"));
    let short = ROWS
        .lines()
        .nth(1)
        .unwrap()
        .replace("\t5,40,70,\t", "\t5,40,\t");
    let f = plain_file(dir.path(), "s.psl", &format!("{short}\n"));
    convert(&f)
        .failure()
        .code(1)
        .stderr(contains("qStarts lists 2"));
}

#[test]
fn rmsk_error_names_the_file_fields() {
    let dir = tempdir();
    let bad = RMSK_ROWS.replace("10001   10468", "10   5");
    let f = plain_file(dir.path(), "r.out", &bad);
    convert(&f)
        .failure()
        .code(1)
        .stderr(contains("begin 10 is greater than end 5"));
}

#[test]
fn bed_lookalikes_are_not_psl() {
    let dir = tempdir();
    // 21 columns on an integer chromosome, with a strand and a name.
    let look1 = "1\t100\t200\t3\t4\t5\t6\t7\t+\tname\t9\t10\t11\tx\t13\t14\t15\t16\t17\t18\t19\n";
    // The same, with PSL-looking sizes but no comma-terminated block lists.
    let look2 = "1\t2\t3\t4\t5\t6\t7\t8\t+\tq\t10\t1\t5\tt\t100\t1\t5\t1\t5\t1\t1\n";
    // A chromosome called psLayoutX.
    let look4 = "psLayoutX\t1\t2\n";
    for (name, text) in [
        ("look1.bed", look1),
        ("look2.tsv", look2),
        ("look4.bed", look4),
    ] {
        let f = plain_file(dir.path(), name, text);
        gloci().arg("merge").arg(&f).assert().success();
    }
}

#[test]
fn gate_message_mentions_format_override() {
    let dir = tempdir();
    let f = plain_file(dir.path(), "x.psl", ROWS);
    gloci()
        .arg("merge")
        .arg(&f)
        .assert()
        .code(2)
        .stderr(contains("--format bed"));
}
