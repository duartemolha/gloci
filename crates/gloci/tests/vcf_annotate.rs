//! `gloci vcf annotate`: determinism, BCF and stdin, order errors, bounded
//! memory, unsupported options, headers and error messages. The record
//! results against bcftools are the `vcf-annotate-*` golden cases.
mod common;
use common::*;
use predicates::str::contains;
use std::path::{Path, PathBuf};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/vcf/annotate")
        .join(name)
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn out(args: &[&str]) -> String {
    let o = gloci().args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

fn records(text: &str) -> Vec<String> {
    text.lines()
        .filter(|l| !l.starts_with('#'))
        .map(String::from)
        .collect()
}

fn header(text: &str) -> Vec<String> {
    text.lines()
        .filter(|l| l.starts_with('#'))
        .map(String::from)
        .collect()
}

fn fx(name: &str) -> String {
    s(&fixture(name))
}

#[test]
fn output_is_identical_at_every_thread_count() {
    let (f, src, tab, bed) = (
        fx("file.vcf"),
        fx("src.vcf"),
        fx("ann.tab"),
        fx("regions.bed"),
    );
    let hdr = fx("tags.hdr");
    let runs: Vec<Vec<&str>> = vec![
        vec!["-a", &src, "-c", "INFO,ID,QUAL,FILTER", "-m", "+HIT"],
        vec![
            "-a",
            &tab,
            "--header-lines",
            &hdr,
            "-c",
            "CHROM,POS,REF,ALT,ID,-,-,-,-,INFO/TA,-,-,INFO/TR",
        ],
        vec![
            "-a",
            &bed,
            "--header-lines",
            &hdr,
            "-c",
            "CHROM,FROM,TO,INFO/TS",
        ],
        vec!["-x", "FORMAT,INFO/DP", "-I", "%CHROM:%POS"],
    ];
    for extra in runs {
        let one = out(&[&["vcf", "annotate", &f, "-t", "1"], &extra[..]].concat());
        assert!(!records(&one).is_empty());
        for t in ["2", "8"] {
            assert_eq!(
                out(&[&["vcf", "annotate", &f, "-t", t], &extra[..]].concat()),
                one,
                "{extra:?} at -t {t}"
            );
        }
    }
}

#[test]
fn bcf_bgzf_and_stdin_give_the_same_records() {
    let d = tempdir();
    let (f, src) = (fx("file.vcf"), fx("src.vcf"));
    let args = ["-a", src.as_str(), "-c", "INFO,ID,FILTER"];
    let want = records(&out(&[&["vcf", "annotate", &f], &args[..]].concat()));
    assert!(want.iter().any(|r| r.contains("rs1")));
    // BCF FILE and BCF source, BCF output.
    let fb = s(&d.path().join("f.bcf"));
    let sb = s(&d.path().join("s.bcf"));
    out(&["view", &f, "--vcf-span", "ref", "-o", &fb]);
    out(&["view", &src, "--vcf-span", "ref", "-o", &sb]);
    let got = out(&[
        "vcf",
        "annotate",
        &fb,
        "-a",
        &sb,
        "-c",
        "INFO,ID,FILTER",
        "-O",
        "vcf",
    ]);
    assert_eq!(records(&got), want);
    let ob = s(&d.path().join("out.bcf"));
    out(&[
        "vcf",
        "annotate",
        &fb,
        "-a",
        &sb,
        "-c",
        "INFO,ID,FILTER",
        "-o",
        &ob,
    ]);
    assert_eq!(
        records(&out(&["view", &ob, "--vcf-span", "ref", "-O", "vcf"])),
        want
    );
    // A BGZF FILE and source, FILE on stdin.
    let fz = bgzf_file(d.path(), "f.vcf.gz", &std::fs::read_to_string(&f).unwrap());
    let sz = bgzf_file(
        d.path(),
        "s.vcf.gz",
        &std::fs::read_to_string(&src).unwrap(),
    );
    let got = out(&[
        "vcf",
        "annotate",
        &s(&fz),
        "-a",
        &s(&sz),
        "-c",
        "INFO,ID,FILTER",
    ]);
    assert_eq!(records(&got), want);
    let o = gloci()
        .args(["vcf", "annotate", "-", "-a", &src, "-c", "INFO,ID,FILTER"])
        .write_stdin(std::fs::read_to_string(&f).unwrap())
        .assert()
        .success();
    assert_eq!(
        records(&String::from_utf8(o.get_output().stdout.clone()).unwrap()),
        want
    );
}

#[test]
fn unsupported_bcftools_options_are_usage_errors_naming_the_option() {
    let f = fx("file.vcf");
    for (args, named) in [
        (vec!["--rename-chrs", "map.txt"], "--rename-chrs"),
        (vec!["-s", "S1"], "-s/--samples"),
        (vec!["-S", "samples.txt"], "-s/--samples"),
        (vec!["-i", "QUAL>1"], "-i/--include"),
        (vec!["--exclude", "QUAL>1"], "-i/--include"),
        (vec!["-l", "AF:sum"], "-l/--merge-logic"),
        (vec!["--merge-logic", "AF:sum"], "-l/--merge-logic"),
        (vec!["--single-overlaps"], "--single-overlaps"),
        (vec!["--rename-annots", "m.txt"], "--rename-annots"),
        (vec!["--pair-logic", "all"], "--pair-logic"),
    ] {
        let assert = gloci()
            .args(["vcf", "annotate", &f])
            .args(&args)
            .assert()
            .failure()
            .stderr(contains(named));
        assert_eq!(assert.get_output().status.code(), Some(2), "{args:?}");
    }
    gloci()
        .args(["vcf", "annotate", &f, "--rename-chrs", "m"])
        .assert()
        .failure()
        .stderr(contains("--chrom-map"));
    gloci()
        .args(["vcf", "annotate", &f, "-i", "QUAL>1"])
        .assert()
        .failure()
        .stderr(contains("gloci vcf filter"));
}

#[test]
fn column_and_remove_mistakes_are_usage_errors_that_say_what_to_do() {
    let (f, src, tab) = (fx("file.vcf"), fx("src.vcf"), fx("ann.tab"));
    for (args, want) in [
        (
            vec!["-a", &src, "-c", "-INFO/AF"],
            "not implemented by bcftools",
        ),
        (vec!["-a", &src, "-c", "=INFO/AF"], "=TAG form"),
        (vec!["-a", &src, "-c", "FORMAT/DP"], "FORMAT"),
        (vec!["-a", &src, "-c", "ALT"], "REF or ALT"),
        (
            vec!["-a", &src, "-c", "INFO/NOPE"],
            "not defined in the annotation file",
        ),
        (vec!["-c", "ID"], "needs --annotations"),
        (vec!["-a", &tab], "needs --columns"),
        (
            vec!["-a", &tab, "-c", "CHROM,POS,INFO/TF"],
            "not defined in the input header",
        ),
        (vec!["-a", &tab, "-c", "ID"], "CHROM"),
        (vec!["-a", &tab, "-c", "CHROM,ID"], "POS"),
        (vec!["-a", &tab, "-c", "CHROM,POS,REF,ID"], "REF needs ALT"),
        (vec!["-x", "NOPE"], "no matching tag"),
        (vec!["-x", "DP"], "did you mean INFO/DP"),
        (vec!["-x", "FILTER/nope"], "not declared"),
        (
            vec!["-x", "INFO/DP", "-I", "%INFO/DP"],
            "runs before --set-id",
        ),
        (vec!["-I", "%POS_%REF"], "%POS_"),
        (
            vec!["-I", "%INFO/NOPE"],
            "no such tag defined in the VCF header",
        ),
        (vec!["-m", "-X"], "needs --annotations"),
        (
            vec!["-a", &src, "-c", "ID", "--min-overlap", "0.5"],
            "intervals",
        ),
        (vec!["-a", "-"], "needs a path"),
    ] {
        let assert = gloci()
            .args(["vcf", "annotate", &f])
            .args(&args)
            .assert()
            .failure()
            .stderr(contains(want));
        assert_eq!(assert.get_output().status.code(), Some(2), "{args:?}");
    }
}

#[test]
fn an_unsorted_file_or_source_is_an_error_with_the_fix() {
    let d = tempdir();
    let hdr = "##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n##contig=<ID=2,length=1000>\n##INFO=<ID=X,Number=1,Type=String,Description=\"x\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let good =
        format!("{hdr}1\t10\ta\tA\tG\t.\t.\t.\n1\t20\tb\tA\tG\t.\t.\t.\n2\t5\tc\tA\tG\t.\t.\t.\n");
    let bad = format!("{hdr}1\t20\tb\tA\tG\t.\t.\t.\n1\t10\ta\tA\tG\t.\t.\t.\n");
    let wrong_chrom = format!("{hdr}2\t5\tc\tA\tG\t.\t.\t.\n1\t10\ta\tA\tG\t.\t.\t.\n");
    let g = s(&plain_file(d.path(), "good.vcf", &good));
    let b = s(&plain_file(d.path(), "bad.vcf", &bad));
    let w = s(&plain_file(d.path(), "wrongchrom.vcf", &wrong_chrom));
    for (file, src) in [(&b, &g), (&g, &b), (&w, &g), (&g, &w)] {
        let assert = gloci()
            .args(["vcf", "annotate", file, "-a", src, "-c", "ID"])
            .assert()
            .failure()
            .stderr(contains("is not sorted"))
            .stderr(contains("gloci sort"));
        assert_eq!(assert.get_output().status.code(), Some(1));
    }
    // A tab source, positional and with intervals.
    let hdrl = d.path().join("x.hdr");
    std::fs::write(
        &hdrl,
        "##INFO=<ID=Y,Number=1,Type=String,Description=\"y\">\n",
    )
    .unwrap();
    let t1 = s(&plain_file(d.path(), "u1.tab", "1\t20\tx\n1\t10\ty\n"));
    let t2 = s(&plain_file(
        d.path(),
        "u2.tab",
        "1\t20\t30\tx\n1\t10\t12\ty\n",
    ));
    let t3 = s(&plain_file(
        d.path(),
        "u3.tab",
        "2\t1\t5\tx\n1\t10\t12\ty\n",
    ));
    for (t, cols) in [
        (&t1, "CHROM,POS,INFO/Y"),
        (&t2, "CHROM,FROM,TO,INFO/Y"),
        (&t3, "CHROM,FROM,TO,INFO/Y"),
    ] {
        let assert = gloci()
            .args(["vcf", "annotate", &g, "-a", t, "-c", cols, "--header-lines"])
            .arg(&hdrl)
            .assert()
            .failure()
            .stderr(contains("is not sorted"))
            .stderr(contains("gloci sort"));
        assert_eq!(assert.get_output().status.code(), Some(1), "{t}");
    }
    // FILE's own order violation is named as well.
    gloci()
        .args(["vcf", "annotate", &b])
        .args(["-m", "+X"])
        .assert()
        .success();
    gloci()
        .args(["vcf", "annotate", &b, "-a", &g, "-c", "ID"])
        .assert()
        .failure()
        .stderr(contains("bad.vcf"));
}

#[test]
fn a_source_chromosome_the_header_lacks_is_an_error_not_an_order_error() {
    let d = tempdir();
    let hdr = "##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n##INFO=<ID=Y,Number=1,Type=String,Description=\"y\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let file = s(&plain_file(
        d.path(),
        "f.vcf",
        &format!("{hdr}1\t10\ta\tA\tG\t.\t.\t.\n"),
    ));
    let svcf = s(&plain_file(
        d.path(),
        "s.vcf",
        &format!("{hdr}1\t10\tb\tA\tG\t.\t.\t.\ncA\t5\tc\tA\tG\t.\t.\t.\n"),
    ));
    let t1 = s(&plain_file(d.path(), "s1.tab", "1\t10\tx\ncA\t5\ty\n"));
    let t2 = s(&plain_file(
        d.path(),
        "s2.tab",
        "1\t10\t20\tx\ncA\t5\t6\ty\n",
    ));
    let h = d.path().join("y.hdr");
    std::fs::write(&h, "##INFO=<ID=Y,Number=1,Type=String,Description=\"y\">\n").unwrap();
    for (src, cols) in [
        (&svcf, "ID"),
        (&t1, "CHROM,POS,INFO/Y"),
        (&t2, "CHROM,FROM,TO,INFO/Y"),
    ] {
        let assert = gloci()
            .args([
                "vcf",
                "annotate",
                &file,
                "-a",
                src,
                "-c",
                cols,
                "--header-lines",
            ])
            .arg(&h)
            .assert()
            .failure()
            .stderr(contains(
                "chromosome \"cA\" is not in the ##contig lines of",
            ))
            .stderr(contains("add ##contig=<ID=cA>"));
        assert_eq!(assert.get_output().status.code(), Some(1), "{src}");
        let err = String::from_utf8_lossy(&assert.get_output().stderr).into_owned();
        assert!(!err.contains("not sorted"), "{err}");
    }
    // A source chromosome the header does declare is fine, even if FILE has no record on it.
    let ok = s(&plain_file(d.path(), "ok.tab", "1\t10\tx\n"));
    gloci()
        .args(["vcf", "annotate", &file, "-a", &ok, "-c", "CHROM,POS,ID"])
        .assert()
        .success();
}

#[test]
fn keep_sites_and_regions_are_refused_with_the_alternative() {
    let f = fx("file.vcf");
    for (args, named, hint) in [
        (vec!["-k"], "-k/--keep-sites", "gloci vcf filter"),
        (vec!["-r", "1:1-5"], "-r/--regions", "gloci view"),
        (vec!["-R", "r.bed"], "-R/--regions-file", "gloci view"),
        (vec!["-H", "##x=1"], "-H/--header-line", "--header-lines"),
        (vec!["-C", "cols.txt"], "-C/--columns-file", "--columns"),
    ] {
        let assert = gloci()
            .args(["vcf", "annotate", &f])
            .args(&args)
            .assert()
            .failure()
            .stderr(contains(named))
            .stderr(contains(hint));
        assert_eq!(assert.get_output().status.code(), Some(2), "{args:?}");
    }
}

#[test]
fn mark_sites_header_matches_bcftools_and_alt_without_ref_is_ignored() {
    let (f, tab, hdr) = (fx("file.vcf"), fx("ann.tab"), fx("tags.hdr"));
    let text = out(&["vcf", "annotate", &f, "-m", "+ALL"]);
    assert!(header(&text).contains(
        &"##INFO=<ID=ALL,Number=0,Type=Flag,Description=\"Sites marked with `bcftools annotate -m ALL`\">"
            .to_string()
    ));
    let text = out(&[
        "vcf",
        "annotate",
        &f,
        "-a",
        &tab,
        "-c",
        "CHROM,POS,ID",
        "-m",
        "-MISS",
    ]);
    assert!(
        header(&text)
            .iter()
            .any(|l| l.contains("Description=\"Sites not listed in MISS\""))
    );
    // ALT without REF changes nothing: rows match by position.
    let a = out(&[
        "vcf",
        "annotate",
        &f,
        "-a",
        &tab,
        "--header-lines",
        &hdr,
        "-c",
        "CHROM,POS,-,ALT,ID",
    ]);
    let b = out(&[
        "vcf",
        "annotate",
        &f,
        "-a",
        &tab,
        "--header-lines",
        &hdr,
        "-c",
        "CHROM,POS,-,-,ID",
    ]);
    assert_eq!(a, b);
    gloci()
        .args([
            "vcf",
            "annotate",
            &f,
            "-a",
            &tab,
            "--header-lines",
            &hdr,
            "-c",
            "CHROM,POS,REF,ID",
        ])
        .assert()
        .failure()
        .stderr(contains("REF needs ALT"));
    gloci()
        .args(["vcf", "annotate", &f, "-I", "%ALT"])
        .assert()
        .failure()
        .stderr(contains("%ALT"));
}

/// `n` records on two contigs, one per 10 bases.
fn big_vcf(n: u64, info: &str) -> String {
    let mut v = format!(
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=100000000>\n##contig=<ID=2,length=100000000>\n{info}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    );
    for c in ["1", "2"] {
        for i in 1..=n / 2 {
            v.push_str(&format!("{c}\t{}\t.\tA\tG\t.\t.\t.\n", i * 10));
        }
    }
    v
}

fn peak(stderr: &[u8]) -> u64 {
    let err = String::from_utf8_lossy(stderr);
    err.split("at most ")
        .nth(1)
        .and_then(|t| t.split(' ').next())
        .and_then(|t| t.parse().ok())
        .unwrap_or_else(|| panic!("no peak in: {err}"))
}

#[test]
fn neither_file_is_held_whole() {
    let d = tempdir();
    let n = 200_000u64;
    let f = s(&plain_file(d.path(), "f.vcf", &big_vcf(n, "")));
    let src_vcf = {
        let mut v = big_vcf(n, "##INFO=<ID=Z,Number=1,Type=Integer,Description=\"z\">\n");
        v = v.replace("\t.\t.\t.\n", "\t99\t.\tZ=1\n");
        v
    };
    let sv = s(&plain_file(d.path(), "s.vcf", &src_vcf));
    // A VCF source, a positional tab file and an interval file.
    let mut tab = String::new();
    let mut bed = String::new();
    for c in ["1", "2"] {
        for i in 1..=n / 2 {
            tab.push_str(&format!("{c}\t{}\t5\n", i * 10));
            bed.push_str(&format!("{c}\t{}\t{}\t5\n", i * 10 - 5, i * 10 + 5));
        }
    }
    let st = s(&plain_file(d.path(), "s.tab", &tab));
    let sb = s(&plain_file(d.path(), "s.bed", &bed));
    let hdr = d.path().join("h.hdr");
    std::fs::write(
        &hdr,
        "##INFO=<ID=Z,Number=1,Type=Integer,Description=\"z\">\n",
    )
    .unwrap();
    let hdr = s(&hdr);
    for args in [
        vec!["-a", &sv, "-c", "ID,INFO/Z"],
        vec!["-a", &st, "--header-lines", &hdr, "-c", "CHROM,POS,INFO/Z"],
        vec![
            "-a",
            &sb,
            "--header-lines",
            &hdr,
            "-c",
            "CHROM,FROM,TO,INFO/Z",
        ],
    ] {
        let o = gloci()
            .args(["-v", "vcf", "annotate", &f, "-o", "/dev/null"])
            .args(&args)
            .assert()
            .success();
        let held = peak(&o.get_output().stderr);
        assert!(held <= 16, "{args:?}: held {held} records of {n}");
        assert!(
            String::from_utf8_lossy(&o.get_output().stderr).contains(&format!("{n} records read"))
        );
    }
}

#[test]
fn header_lines_declared_and_removed() {
    let (f, src) = (fx("file.vcf"), fx("src.vcf"));
    let text = out(&[
        "vcf",
        "annotate",
        &f,
        "-a",
        &src,
        "-c",
        "INFO/AF,FILTER,NEW:=INFO/NAME",
    ]);
    let h = header(&text);
    assert!(
        h.iter()
            .any(|l| l.starts_with("##INFO=<ID=AF,Number=A,Type=Float"))
    );
    assert!(
        h.iter()
            .any(|l| l.starts_with("##INFO=<ID=NEW,Number=1,Type=String"))
    );
    assert!(h.iter().any(|l| l.starts_with("##FILTER=<ID=LowQ,")));
    assert!(!h.iter().any(|l| l.starts_with("##INFO=<ID=DB")));
    assert!(!h.iter().any(|l| l.contains("glociCommand")));
    // The header is FILE's, with no PASS line added.
    assert!(!h.iter().any(|l| l.starts_with("##FILTER=<ID=PASS")));
    // Removed tags lose their header lines; GT stays with -x FORMAT.
    let text = out(&["vcf", "annotate", &f, "-x", "INFO/DP,FORMAT,FILTER"]);
    let h = header(&text);
    assert!(!h.iter().any(|l| l.starts_with("##INFO=<ID=DP,")));
    assert!(h.iter().any(|l| l.starts_with("##INFO=<ID=OLD,")));
    assert!(h.iter().any(|l| l.starts_with("##FORMAT=<ID=GT,")));
    assert!(!h.iter().any(|l| l.starts_with("##FORMAT=<ID=DP,")));
    assert!(!h.iter().any(|l| l.starts_with("##FILTER=")));
    // --header-lines adds lines once, and the mark-sites tag is declared.
    let d = tempdir();
    let hl = d.path().join("extra.hdr");
    std::fs::write(
        &hl,
        "##INFO=<ID=DP,Number=1,Type=Float,Description=\"again\">\n##INFO=<ID=EXTRA,Number=0,Type=Flag,Description=\"e\">\n\n##source=test\n",
    )
    .unwrap();
    let text = out(&[
        "vcf",
        "annotate",
        &f,
        "--header-lines",
        &s(&hl),
        "-m",
        "+MARK",
    ]);
    let h = header(&text);
    assert_eq!(
        h.iter().filter(|l| l.starts_with("##INFO=<ID=DP,")).count(),
        1
    );
    assert!(h.iter().any(|l| l.starts_with("##INFO=<ID=EXTRA,")));
    assert!(h.iter().any(|l| l == "##source=test"));
    assert!(
        h.iter()
            .any(|l| l.starts_with("##INFO=<ID=MARK,Number=0,Type=Flag"))
    );
    // A line that is not a header line is an error.
    std::fs::write(&hl, "INFO=x\n").unwrap();
    gloci()
        .args(["vcf", "annotate", &f, "--header-lines", &s(&hl)])
        .assert()
        .failure()
        .stderr(contains("not a header line"));
}

#[test]
fn data_errors_say_where_and_how_to_fix() {
    let d = tempdir();
    let (f, bad) = (fx("file.vcf"), fx("src-badcount.vcf"));
    let assert = gloci()
        .args(["vcf", "annotate", &f, "-a", &bad, "-c", "INFO/AF"])
        .assert()
        .failure()
        .stderr(contains(
            "incorrect number of values (3) for the AF tag at 1:200",
        ));
    assert_eq!(assert.get_output().status.code(), Some(1));
    let hdr = fx("tags.hdr");
    let t = s(&plain_file(d.path(), "q.tab", "1\t100\tx\n"));
    let assert = gloci()
        .args(["vcf", "annotate", &f, "-a", &t, "-c", "CHROM,POS,QUAL"])
        .assert()
        .failure()
        .stderr(contains("could not parse QUAL at 1:100"))
        .stderr(contains("q.tab line 1"));
    assert_eq!(assert.get_output().status.code(), Some(1));
    let t = s(&plain_file(d.path(), "i.tab", "1\t100\tx\n"));
    gloci()
        .args(["vcf", "annotate", &f, "-a", &t, "--header-lines", &hdr])
        .args(["-c", "CHROM,POS,INFO/TI"])
        .assert()
        .failure()
        .stderr(contains("could not parse TI at 1:100"));
    let t = s(&plain_file(d.path(), "f.tab", "1\t100\tnotafilter\n"));
    gloci()
        .args(["vcf", "annotate", &f, "-a", &t, "-c", "CHROM,POS,FILTER"])
        .assert()
        .failure()
        .stderr(contains("FILTER `notafilter` is not defined"))
        .stderr(contains("--header-lines"));
    let t = s(&plain_file(d.path(), "short.tab", "1\t100\n"));
    gloci()
        .args(["vcf", "annotate", &f, "-a", &t, "-c", "CHROM,POS,ID"])
        .assert()
        .failure()
        .stderr(contains("short.tab line 1"))
        .stderr(contains("no column 3, which --columns names"));
}

#[test]
fn the_pairing_order_at_one_position_matches_the_synced_reader() {
    // The exact pair is written first, then the partial pair, then the
    // duplicate; the golden case pins the rest of the file against bcftools.
    let (f, src) = (fx("file.vcf"), fx("src.vcf"));
    let text = out(&["vcf", "annotate", &f, "-a", &src, "-c", "ID"]);
    let at200: Vec<(String, String)> = records(&text)
        .iter()
        .filter(|r| r.starts_with("1\t200\t"))
        .map(|r| {
            let c: Vec<&str> = r.split('\t').collect();
            (c[2].to_string(), c[4].to_string())
        })
        .collect();
    assert_eq!(
        at200,
        vec![
            ("rs2b".to_string(), "G".to_string()),
            ("rs2".to_string(), "T,G".to_string()),
            (".".to_string(), "G".to_string()),
        ]
    );
}

#[cfg(unix)]
#[test]
fn pipes_and_process_substitution_work_for_file_and_every_source_kind() {
    let d = tempdir();
    let (f, src, tab, reg, hdr) = (
        fx("file.vcf"),
        fx("src.vcf"),
        fx("ann.tab"),
        fx("regions.tab"),
        fx("tags.hdr"),
    );
    let sb = s(&d.path().join("s.bcf"));
    out(&["view", &src, "--vcf-span", "ref", "-o", &sb]);
    let run = |args: &str, stdin: &str| {
        let o = bash_gloci(args, stdin);
        assert!(
            o.status.success(),
            "{args}: {}",
            String::from_utf8_lossy(&o.stderr)
        );
        String::from_utf8(o.stdout).unwrap()
    };
    let file_text = std::fs::read_to_string(&f).unwrap();
    // (columns, plain source, ways to pipe it)
    for (cols, plain) in [
        ("-c INFO,ID,FILTER", src.clone()),
        ("-c INFO,ID,FILTER", sb.clone()),
        (
            "-c CHROM,POS,REF,ALT,ID,QUAL,-,-,-,INFO/TA --header-lines",
            tab.clone(),
        ),
        (
            "-c CHROM,FROM,TO,INFO/TS,INFO/TI --header-lines",
            reg.clone(),
        ),
    ] {
        let extra = if cols.ends_with("--header-lines") {
            format!(" {hdr}")
        } else {
            String::new()
        };
        let want = run(&format!("vcf annotate {f} -a {plain} {cols}{extra}"), "");
        assert!(want.lines().filter(|l| !l.starts_with('#')).count() > 5);
        // A BCF source is already BGZF; compressing it again is not a BCF.
        let mut pipes = vec![format!("<(cat {plain})")];
        if !plain.ends_with(".bcf") {
            pipes.push(format!("<(gzip -c {plain})"));
            pipes.push(format!("<(cat {plain} | gzip -c)"));
        }
        for pipe in pipes {
            for file in [format!("<(cat {f})"), "-".to_string()] {
                let stdin = if file == "-" { file_text.as_str() } else { "" };
                let got = run(
                    &format!("vcf annotate {file} -a {pipe} {cols}{extra}"),
                    stdin,
                );
                assert_eq!(got, want, "{file} {pipe} {cols}");
            }
        }
    }
}
