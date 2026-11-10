//! `gloci view` on BAM and SAM: every record, formats, `@PG`, regions, filters, errors (D104,
//! D116, D120).
mod common;
use common::*;
use gloci_testkit::differential::compare;
use gloci_testkit::fixtures;
use predicates::str::contains;
use std::path::{Path, PathBuf};

fn aln(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/aln")
        .join(name)
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn run(args: &[&str]) -> Vec<u8> {
    gloci()
        .args(args)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone()
}

fn text(b: &[u8]) -> String {
    String::from_utf8(b.to_vec()).unwrap()
}

/// A SAM with `n` reads over chr1 and chr2, a few unmapped, long enough for many BGZF blocks.
fn big_sam(n: usize) -> String {
    let mut sam = String::from(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10000000\n@SQ\tSN:chr2\tLN:10000000\n",
    );
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let qual = "I".repeat(seq.len());
    for i in 0..n {
        let (chrom, pos) = if i < n / 2 {
            ("chr1", i * 50 + 1)
        } else {
            ("chr2", (i - n / 2) * 50 + 1)
        };
        let (flag, cigar) = if i % 97 == 0 {
            (4, "*")
        } else {
            (i % 2 * 16, "64M")
        };
        sam.push_str(&format!(
            "read{i}\t{flag}\t{chrom}\t{pos}\t{}\t{cigar}\t*\t0\t0\t{seq}\t{qual}\tNM:i:{}\tXF:f:0.1\n",
            i % 61,
            i % 7
        ));
    }
    sam
}

#[test]
fn view_bam_output_is_identical_at_every_thread_count() {
    let dir = tempdir();
    let sam = plain_file(dir.path(), "x.sam", &big_sam(30_000));
    // (noodles' encoder gives every integer tag the `i` type, samtools the smallest.)
    let bam = fixtures::write(dir.path(), "x.bam", &run(&["view", &s(&sam), "-O", "bam"]));
    let mut outs = Vec::new();
    for t in ["1", "2", "8"] {
        for (name, input) in [("bam", &bam), ("sam", &sam)] {
            let out = dir.path().join(format!("y{name}{t}.bam"));
            gloci()
                .args(["view", &s(input), "-o", &s(&out), "-t", t])
                .assert()
                .success();
            outs.push((name, std::fs::read(&out).unwrap()));
        }
    }
    for (name, o) in &outs {
        let (_, first) = outs.iter().find(|(n, _)| n == name).unwrap();
        // Only the @PG line differs: it records the command line, `-t` included.
        assert!(
            compare("bam", o, first),
            "{name}: differs with the thread count"
        );
    }
    // Every BGZF block after the first (which holds the @PG line, with `-t` in it) is the same
    // bytes at every thread count.
    let blocks = |b: &[u8]| -> Vec<Vec<u8>> {
        let (mut v, mut at) = (Vec::new(), 0);
        while at < b.len() {
            let size = usize::from(u16::from_le_bytes([b[at + 16], b[at + 17]])) + 1;
            v.push(b[at..at + size].to_vec());
            at += size;
        }
        v
    };
    for name in ["bam", "sam"] {
        let by_t: Vec<_> = outs
            .iter()
            .filter(|(n, _)| *n == name)
            .map(|(_, o)| blocks(o))
            .collect();
        assert!(by_t[0].len() > 10, "{name}: a big output");
        for t in &by_t[1..] {
            assert_eq!(t.len(), by_t[0].len(), "{name}");
            assert!(
                t[1..] == by_t[0][1..],
                "{name}: blocks differ with the thread count"
            );
        }
    }
    // SAM and BAM input give the same records.
    assert!(compare("bam", &outs[0].1, &outs[1].1));
    // SAM text out is the same at every thread count too.
    let texts: Vec<_> = ["1", "2", "8"]
        .iter()
        .map(|t| run(&["view", &s(&bam), "-O", "sam", "-t", t]))
        .collect();
    assert_eq!(texts[0], texts[1]);
    assert_eq!(texts[0], texts[2]);
}

#[test]
fn view_sam_is_a_passthrough_with_a_pg_line() {
    let sam = aln("order.sam");
    let out = text(&run(&["view", &s(&sam)]));
    let want = std::fs::read_to_string(&sam).unwrap();
    let pg: Vec<_> = out.lines().filter(|l| l.starts_with("@PG")).collect();
    assert_eq!(pg.len(), 1, "{out}");
    assert!(
        pg[0].starts_with("@PG\tID:gloci\tPN:gloci\tVN:"),
        "{}",
        pg[0]
    );
    assert!(
        pg[0].contains("\tCL:") && pg[0].contains(" view "),
        "{}",
        pg[0]
    );
    let rest: String = out
        .lines()
        .filter(|l| !l.starts_with("@PG"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(rest, want);
    // Through a pipe, a BAM goes back to SAM and gets a second @PG line.
    let bam = run(&["view", &s(&sam), "-O", "bam"]);
    let again = gloci()
        .args(["view", "-", "-O", "sam"])
        .write_stdin(bam)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let again = text(&again);
    assert_eq!(
        again.lines().filter(|l| l.starts_with("@PG")).count(),
        2,
        "{again}"
    );
    assert!(
        again.contains("@PG\tID:gloci.1\tPN:gloci\tPP:gloci\t"),
        "{again}"
    );
}

#[test]
fn view_keeps_unmapped_reads_and_the_filters_drop_them() {
    let bam = s(&aln("flags.bam"));
    let all = text(&run(&["view", &bam, "-O", "sam"]));
    assert_eq!(all.lines().filter(|l| !l.starts_with('@')).count(), 14);
    assert!(all.contains("r_unmapped_unplaced") && all.contains("r_unmapped_placed"));
    let mapped = text(&run(&[
        "view",
        &bam,
        "-O",
        "sam",
        "--exclude-flags",
        "unmap",
    ]));
    assert_eq!(mapped.lines().filter(|l| !l.starts_with('@')).count(), 12);
    let q = text(&run(&[
        "view",
        &bam,
        "-O",
        "sam",
        "--min-mapq",
        "20",
        "--require-flags",
        "0x1",
    ]));
    let names: Vec<_> = q
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(|l| &l[..6])
        .collect();
    assert_eq!(names, ["r_pair", "r_pair"]);
}

#[test]
fn view_regions_of_alignments() {
    let bam = s(&aln("flags.bam"));
    let names = |args: &[&str]| -> Vec<String> {
        let mut a = vec!["view", &bam, "-O", "sam", "--no-header"];
        a.extend_from_slice(args);
        text(&run(&a))
            .lines()
            .map(|l| l.split('\t').next().unwrap().to_string())
            .collect()
    };
    // 1-based closed regions; an unmapped read placed at chr1:701 counts as one base there.
    assert_eq!(names(&["chr1:701-701"]), ["r_unmapped_placed"]);
    assert_eq!(names(&["chr1:700-700"]), Vec::<String>::new());
    assert_eq!(names(&["chr1:101-101"]), ["r_plain", "r_rev"]);
    assert_eq!(names(&["chr2:51-60"]), ["read9", "read10"]);
    // Each read once, in file order, even when regions overlap or repeat.
    assert_eq!(
        names(&["chr2", "chr1:101-110", "chr1:105-120"]),
        ["r_plain", "r_rev", "read9", "read10"]
    );
    // A chromosome the header lacks is a warning, and the others still answer.
    gloci()
        .args([
            "view",
            &bam,
            "-O",
            "sam",
            "--no-header",
            "chrZ",
            "chr2:51-60",
        ])
        .assert()
        .success()
        .stderr(contains("chromosome \"chrZ\" not found in"))
        .stdout(contains("read9"));
    // Unmapped unplaced reads are in no region.
    assert!(!names(&["chr1", "chr2"]).contains(&"r_unmapped_unplaced".to_string()));
    // A regions file (BED, 0-based).
    let dir = tempdir();
    let bed = plain_file(dir.path(), "r.bed", "chr1\t700\t701\n");
    assert_eq!(
        text(&run(&[
            "view",
            &bam,
            "-O",
            "sam",
            "--no-header",
            "--regions-file",
            &s(&bed)
        ]))
        .lines()
        .count(),
        1
    );
    // SAM input is found the same way.
    let sam = s(&aln("flags.sam"));
    assert_eq!(
        text(&run(&["view", &sam, "--no-header", "chr1:701-701"])),
        text(&run(&[
            "view",
            &bam,
            "-O",
            "sam",
            "--no-header",
            "chr1:701-701"
        ]))
    );
}

#[test]
fn view_header_options_for_alignments() {
    let bam = s(&aln("flags.bam"));
    let h = text(&run(&["view", &bam, "-O", "sam", "--header-only"]));
    assert!(
        h.starts_with("@HD") && h.lines().all(|l| l.starts_with('@')),
        "{h}"
    );
    // --header-only to BAM is a header with no records.
    let b = run(&["view", &bam, "-O", "bam", "--header-only"]);
    let back = gloci()
        .args(["view", "-", "-O", "sam", "--no-header"])
        .write_stdin(b)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    assert!(back.is_empty());
    // BAM records need their header.
    gloci()
        .args(["view", &bam, "-O", "bam", "--no-header"])
        .assert()
        .code(2)
        .stderr(contains("--no-header cannot be used with BAM output"));
    gloci()
        .args(["view", &bam, "--no-header", "-o", "/dev/null"])
        .assert()
        .code(2)
        .stderr(contains("--no-header cannot be used with BAM output"));
}

#[test]
fn view_alignments_usage_errors() {
    let bam = s(&aln("flags.bam"));
    gloci()
        .args(["view", &bam, "-O", "bed"])
        .assert()
        .code(2)
        .stderr(contains("alignments are written as BAM or SAM"))
        .stderr(contains("gloci convert - -o x.bed"));
    gloci()
        .args(["view", &bam, "-O", "vcf"])
        .assert()
        .code(2)
        .stderr(contains("BAM alignments"));
    gloci()
        .args(["view", &bam, "--exclude-flags", "bogus", "-O", "sam"])
        .assert()
        .code(2)
        .stderr(contains("--exclude-flags").and(contains("unknown flag name 'bogus'")));
    gloci()
        .args(["view", &bam, "-o", "x.bam", "--write-index"])
        .assert()
        .code(2)
        .stderr(contains("--write-index is not available for BAM output"));
    // BED input cannot be written as BAM.
    let dir = tempdir();
    let bed = plain_file(dir.path(), "a.bed", "chr1\t1\t2\n");
    gloci()
        .args(["view", &s(&bed), "-O", "bam"])
        .assert()
        .code(2);
    // Binary output is refused on a terminal only; here stdout is a pipe.
    assert!(!run(&["view", &bam, "-O", "bam"]).is_empty());
    // A command that does not read alignments still refuses them, with its own hint.
    gloci().args(["merge", &bam]).assert().failure();
}

#[test]
fn view_reports_bad_records_where_they_are() {
    let dir = tempdir();
    let bad = plain_file(
        dir.path(),
        "bad.sam",
        "@SQ\tSN:chr1\tLN:100\nok\t0\tchr1\t5\t60\t4M\t*\t0\t0\tACGT\tIIII\nbad\t0\tchrZ\t5\t60\t4M\t*\t0\t0\tACGT\tIIII\nok2\t0\tchr1\t6\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
    );
    gloci()
        .args(["view", &s(&bad), "-O", "sam"])
        .assert()
        .failure()
        .stderr(contains("bad.sam line 3:").and(contains("RNAME 'chrZ'")));
    let out = gloci()
        .args(["view", &s(&bad), "-O", "sam", "--skip-invalid"])
        .assert()
        .success()
        .stderr(contains("skipped 1 invalid records").and(contains("line 3:")))
        .get_output()
        .stdout
        .clone();
    let out = text(&out);
    assert!(
        out.contains("\nok\t") && out.contains("\nok2\t") && !out.contains("\nbad\t"),
        "{out}"
    );
    // A BAM cut short is an error naming the record.
    let full = std::fs::read(aln("flags.bam")).unwrap();
    let cut = {
        let mut plain = Vec::new();
        std::io::Read::read_to_end(&mut noodles_bgzf::io::Reader::new(&full[..]), &mut plain)
            .unwrap();
        plain.truncate(plain.len() - 20);
        fixtures::bgzip(&plain)
    };
    let p = fixtures::write(dir.path(), "cut.bam", &cut);
    gloci()
        .args(["view", &s(&p), "-O", "sam"])
        .assert()
        .failure()
        .stderr(contains("cut.bam record").and(contains("of the file)")));
}

use predicates::prelude::*;

#[test]
fn view_sam_rejects_bad_lines_at_read_whatever_the_output() {
    let dir = tempdir();
    let bad = plain_file(
        dir.path(),
        "bad.sam",
        "@SQ\tSN:chr1\tLN:100\nok\t0\tchr1\t5\t60\t4M\t*\t0\t0\tACGT\tIIII\nbad\t0\tchr1\t6\t60\t4M\t*\t0\t0\tACGT\tII\nok2\t0\tchr1\t7\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:Q:1\n",
    );
    for out in ["sam", "bam"] {
        gloci()
            .args(["view", &s(&bad), "-O", out])
            .assert()
            .code(1)
            .stderr(contains("bad.sam line 3:").and(contains("quality string")))
            .stderr(contains("disk").not());
        let o = gloci()
            .args(["view", &s(&bad), "-O", out, "--skip-invalid"])
            .assert()
            .success()
            .stderr(contains("skipped 2 invalid records").and(contains("line 4:")))
            .get_output()
            .stdout
            .clone();
        assert!(!o.is_empty());
    }
}

#[test]
fn view_sam_with_crlf_matches_lf_and_chrom_map_is_a_usage_error() {
    let dir = tempdir();
    let lf = std::fs::read_to_string(aln("cigars.sam")).unwrap();
    let crlf = plain_file(dir.path(), "crlf.sam", &lf.replace('\n', "\r\n"));
    let from_crlf = run(&["view", &s(&crlf)]);
    let from_lf = run(&["view", &s(&aln("cigars.sam"))]);
    let strip = |b: &[u8]| -> String {
        text(b)
            .lines()
            .filter(|l| !l.starts_with("@PG"))
            .map(|l| format!("{l}\n"))
            .collect()
    };
    assert_eq!(strip(&from_crlf), strip(&from_lf));
    assert!(!from_crlf.contains(&b'\r'));
    let map = plain_file(dir.path(), "m.tsv", "chr1\t1\n");
    gloci()
        .args(["view", &s(&aln("cigars.sam")), "--chrom-map", &s(&map)])
        .assert()
        .code(2)
        .stderr(contains("convert it to BAM"));
    // BAM honours the map.
    let o = text(&run(&[
        "view",
        &s(&aln("cigars.bam")),
        "-O",
        "sam",
        "--chrom-map",
        &s(&map),
    ]));
    assert!(
        o.contains("@SQ\tSN:1\tLN:100000") && o.contains("\t1\t101\t"),
        "{o}"
    );
}

#[test]
fn view_reports_what_the_filters_dropped_at_v() {
    let bam = s(&aln("flags.bam"));
    let quiet = gloci()
        .args([
            "view",
            &bam,
            "-O",
            "sam",
            "--exclude-flags",
            "unmap",
            "--min-mapq",
            "20",
        ])
        .assert()
        .success();
    assert!(!String::from_utf8_lossy(&quiet.get_output().stderr).contains("dropped by"));
    gloci()
        .args([
            "-v",
            "view",
            &bam,
            "-O",
            "sam",
            "--exclude-flags",
            "unmap",
            "--min-mapq",
            "20",
        ])
        .assert()
        .success()
        .stderr(contains(
            "flags.bam: 4 reads dropped by --min-mapq, --require-flags or --exclude-flags",
        ));
    // A SAM with a blank line says so too.
    let dir = tempdir();
    let sam = plain_file(
        dir.path(),
        "b.sam",
        "@SQ\tSN:chr1\tLN:100\nok\t0\tchr1\t5\t60\t4M\t*\t0\t0\tACGT\tIIII\n\nok2\t0\tchr1\t6\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
    );
    gloci()
        .args(["-v", "view", &s(&sam)])
        .assert()
        .success()
        .stderr(contains("b.sam: 1 blank lines skipped"));
}

#[test]
fn view_write_index_is_refused_by_output_format() {
    let bam = s(&aln("flags.bam"));
    gloci()
        .args(["view", &bam, "-o", "x.sam", "--write-index"])
        .assert()
        .code(2)
        .stderr(contains("--write-index needs BAM output").and(contains("SAM text has no index")));
    gloci()
        .args(["view", &bam, "-o", "x.bam", "--write-index"])
        .assert()
        .code(2)
        .stderr(contains("not available for BAM output").and(contains("afterwards").not()));
}

#[test]
fn view_sam_to_sam_agrees_with_the_filters_on_normalised_reads() {
    let dir = tempdir();
    let sam = plain_file(
        dir.path(),
        "n.sam",
        "@SQ\tSN:chr1\tLN:1000\nlower\t0\tchr1\t5\t60\t4M\t*\t0\t0\tacgt\tIIII\tXP:i:+7\npos0\t0\tchr1\t0\t60\t4M\t=\t50\t0\tACGT\tIIII\n",
    );
    let out = text(&run(&[
        "view",
        &s(&sam),
        "--require-flags",
        "unmap",
        "--no-header",
    ]));
    assert_eq!(out, "pos0\t4\t*\t0\t60\t4M\t*\t50\t0\tACGT\tIIII\n");
    // A line that is not normalised is copied as written.
    let out = text(&run(&[
        "view",
        &s(&sam),
        "--exclude-flags",
        "unmap",
        "--no-header",
    ]));
    assert_eq!(
        out,
        "lower\t0\tchr1\t5\t60\t4M\t*\t0\t0\tacgt\tIIII\tXP:i:+7\n"
    );
}

#[test]
fn view_sam_holds_what_bam_cannot() {
    let dir = tempdir();
    let sam = plain_file(
        dir.path(),
        "t.sam",
        "@SQ\tSN:chr1\tLN:1000\nt\t1\tchr1\t5\t60\t4M\tchr1\t100\t3000000000\tACGT\tIIII\n",
    );
    let out = text(&run(&["view", &s(&sam), "--no-header"]));
    assert!(out.contains("\t3000000000\t"), "{out}");
    gloci()
        .args(["view", &s(&sam), "-O", "bam"])
        .assert()
        .code(1)
        .stderr(contains("this record cannot be written as BAM").and(contains("-O sam")));
}
