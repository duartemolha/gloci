//! `gloci vcf norm`: BCF, stdin, determinism, the reference index, bounded
//! reference reads, unsupported options and `--check-ref x`.
mod common;
use common::*;
use predicates::str::contains;
use std::path::{Path, PathBuf};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../xtask/diff/inputs/vcf/norm")
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

/// The reference, its index and the fixtures copied to a fresh directory.
fn work() -> (tempfile::TempDir, String, String) {
    let dir = tempdir();
    for f in ["ref.fa", "ref.fa.fai", "align.vcf", "mismatch.vcf"] {
        std::fs::copy(fixture(f), dir.path().join(f)).unwrap();
    }
    let r = s(&dir.path().join("ref.fa"));
    let v = s(&dir.path().join("align.vcf"));
    (dir, r, v)
}

#[test]
fn output_is_identical_at_every_thread_count() {
    let (_d, r, v) = work();
    let one = out(&["vcf", "norm", &v, "-f", &r, "-t", "1"]);
    assert_eq!(records(&one).len(), 33);
    for t in ["2", "8"] {
        assert_eq!(out(&["vcf", "norm", &v, "-f", &r, "-t", t]), one);
    }
    // Splitting first is the same at every thread count too.
    let sp = s(&fixture("split.vcf"));
    let one = out(&["vcf", "norm", &sp, "-f", &r, "--split", "-t", "1"]);
    for t in ["2", "8"] {
        assert_eq!(
            out(&["vcf", "norm", &sp, "-f", &r, "--split", "-t", t]),
            one
        );
    }
}

#[test]
fn split_here_is_split_then_normalise() {
    let (_d, r, _) = work();
    let sp = s(&fixture("split.vcf"));
    let joined = out(&["vcf", "norm", &sp, "-f", &r, "--split"]);
    let piped = bash_gloci(
        &format!(
            "vcf split {sp} | '{}' vcf norm - -f {r}",
            assert_cmd::cargo::cargo_bin("gloci").display()
        ),
        "",
    );
    assert!(
        piped.status.success(),
        "{}",
        String::from_utf8_lossy(&piped.stderr)
    );
    assert_eq!(
        records(&String::from_utf8(piped.stdout).unwrap()),
        records(&joined)
    );
}

#[test]
fn bcf_input_and_output_and_stdin() {
    let (d, r, v) = work();
    let want = records(&out(&["vcf", "norm", &v, "-f", &r]));
    let bcf = s(&d.path().join("in.bcf"));
    out(&["view", &v, "--vcf-span", "ref", "-o", &bcf]);
    assert_eq!(
        records(&out(&["vcf", "norm", &bcf, "-f", &r, "-O", "vcf"])),
        want
    );
    let obcf = s(&d.path().join("out.bcf"));
    out(&["vcf", "norm", &bcf, "-f", &r, "-o", &obcf]);
    assert_eq!(
        records(&out(&["view", &obcf, "--vcf-span", "ref", "-O", "vcf"])),
        want
    );
    let text = std::fs::read_to_string(&v).unwrap();
    let o = gloci()
        .args(["vcf", "norm", "-", "-f", &r])
        .write_stdin(text)
        .assert()
        .success();
    assert_eq!(
        records(&String::from_utf8(o.get_output().stdout.clone()).unwrap()),
        want
    );
}

#[test]
fn a_missing_fai_points_to_gloci_index() {
    let (d, _, v) = work();
    let bare = d.path().join("bare");
    std::fs::create_dir(&bare).unwrap();
    std::fs::copy(fixture("ref.fa"), bare.join("ref.fa")).unwrap();
    gloci()
        .args(["vcf", "norm", &v, "-f", &s(&bare.join("ref.fa"))])
        .assert()
        .failure()
        .stderr(contains("gloci index"));
    assert!(
        !bare.join("ref.fa.fai").exists(),
        "gloci must not build the index"
    );
}

#[test]
fn a_stale_index_warns_and_an_edited_fasta_fails() {
    let (d, r, v) = work();
    // The index is older than the FASTA: a warning, the run goes on.
    let fa = std::fs::File::options()
        .write(true)
        .open(d.path().join("ref.fa"))
        .unwrap();
    fa.set_modified(std::time::SystemTime::now() + std::time::Duration::from_secs(3600))
        .unwrap();
    gloci()
        .args(["vcf", "norm", &v, "-f", &r])
        .assert()
        .success()
        .stderr(contains("older than"));
    // An edited FASTA no longer matches its index: an error.
    let text = std::fs::read_to_string(d.path().join("ref.fa")).unwrap();
    std::fs::write(
        d.path().join("ref.fa"),
        text.replacen("GATCGTACGT", "GATCGTACGTTT", 1),
    )
    .unwrap();
    gloci()
        .args(["vcf", "norm", &v, "-f", &r])
        .assert()
        .failure()
        .stderr(contains("does not match"));
}

#[test]
fn unsupported_bcftools_options_are_usage_errors_naming_the_option() {
    let (_d, r, v) = work();
    for (args, named) in [
        (vec!["-a"], "-a/--atomize"),
        (vec!["--atomize"], "-a/--atomize"),
        (vec!["-m", "+any"], "-m+any"),
        (vec!["-m-snps"], "-m-snps"),
        (vec!["--old-rec-tag", "OLD"], "--old-rec-tag"),
        (vec!["-N"], "-N/--do-not-normalize"),
    ] {
        let assert = gloci()
            .args(["vcf", "norm", &v, "-f", &r])
            .args(&args)
            .assert()
            .failure()
            .stderr(contains(named));
        assert_eq!(assert.get_output().status.code(), Some(2), "{args:?}");
    }
}

#[test]
fn a_mismatch_error_says_where_and_how_to_fix_it() {
    let (d, r, _) = work();
    gloci()
        .args(["vcf", "norm", &s(&d.path().join("mismatch.vcf")), "-f", &r])
        .assert()
        .failure()
        .stderr(contains("chr1:50"))
        .stderr(contains("--check-ref"));
}

#[test]
fn check_ref_x_drops_only_the_failing_records() {
    let (d, r, _) = work();
    let m = s(&d.path().join("mismatch.vcf"));
    let all = out(&["vcf", "norm", &m, "-f", &r, "--check-ref", "w"]);
    let kept = out(&["vcf", "norm", &m, "-f", &r, "--check-ref", "x"]);
    let kept = records(&kept);
    assert_eq!(kept.len(), 4);
    for k in &kept {
        assert!(records(&all).contains(k));
    }
    let pos: Vec<_> = kept
        .iter()
        .map(|l| l.split('\t').take(5).collect::<Vec<_>>().join(" "))
        .collect();
    assert_eq!(
        pos,
        [
            "chr1 46 . GCA G",
            "chr1 49 . TG T",
            "chr1 60 . TC T",
            "chr2 4 . TN T"
        ]
    );
}

#[test]
fn an_unknown_contig_and_a_bad_window_are_errors() {
    let (d, r, _) = work();
    let v = d.path().join("other.vcf");
    std::fs::write(
        &v,
        "##fileformat=VCFv4.2\n##contig=<ID=chrZ,length=10>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrZ\t3\t.\tA\tC\t.\t.\t.\n",
    )
    .unwrap();
    gloci()
        .args(["vcf", "norm", &s(&v), "-f", &r])
        .assert()
        .failure()
        .stderr(contains("chrZ"))
        .stderr(contains("FASTA index"));
    gloci()
        .args(["vcf", "norm", &s(&v), "-f", &r, "--window", "0"])
        .assert()
        .failure();
}

#[test]
fn a_large_contig_is_read_in_bounded_windows() {
    // 24 MiB contig, a 200 kb homopolymer at its end and a deletion at the
    // end of that run, which must shift left across the whole run.
    let dir = tempdir();
    let total = 24usize << 20;
    let run = 200_000usize;
    let mut seq: Vec<u8> = Vec::with_capacity(total);
    let mut x: u64 = 0x9E3779B97F4A7C15;
    while seq.len() < total - run - 1 {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        seq.push(b"ACGT"[(x & 3) as usize]);
    }
    seq.push(b'C');
    seq.extend(std::iter::repeat_n(b'A', run));
    let mut fa = b">big\n".to_vec();
    for line in seq.chunks(60) {
        fa.extend_from_slice(line);
        fa.push(b'\n');
    }
    let fa_path = dir.path().join("big.fa");
    std::fs::write(&fa_path, &fa).unwrap();
    gloci().args(["index", &s(&fa_path)]).assert().success();
    let end = seq.len() as u64; // 1-based position of the last A
    let vcf = format!(
        "##fileformat=VCFv4.2\n##contig=<ID=big,length={}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nbig\t{}\t.\tAA\tA\t.\t.\t.\n",
        seq.len(),
        end - 1
    );
    let v = dir.path().join("big.vcf");
    std::fs::write(&v, vcf).unwrap();
    let o = gloci()
        .args(["-v", "vcf", "norm", &s(&v), "-f", &s(&fa_path)])
        .assert()
        .success();
    let stdout = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    let rec = records(&stdout);
    let anchor = (seq.len() - run) as u64; // the C before the run
    assert_eq!(
        rec[0].split('\t').take(5).collect::<Vec<_>>(),
        ["big", anchor.to_string().as_str(), ".", "CA", "C"]
    );
    let err = String::from_utf8(o.get_output().stderr.clone()).unwrap();
    let peak: u64 = err
        .split("peak reference window ")
        .nth(1)
        .and_then(|t| t.split(' ').next())
        .and_then(|t| t.parse().ok())
        .unwrap_or_else(|| panic!("no peak in: {err}"));
    assert!(peak <= 1 << 17, "held {peak} reference bases");
}

#[test]
fn the_largest_window_does_not_overflow() {
    let (_d, r, v) = work();
    let big = out(&[
        "vcf",
        "norm",
        &v,
        "-f",
        &r,
        "--window",
        "18446744073709551615",
    ]);
    assert_eq!(records(&big).len(), 33);
}

#[test]
fn skip_invalid_warns_once_with_a_count_and_the_first_line() {
    let (d, r, _v) = work();
    let mut t = std::fs::read_to_string(fixture("../split/bad-count.vcf")).unwrap();
    t.push_str("chr1\t20\t.\tA\tC,G\t.\t.\tRC=1,2\tGT\t0/1\t1/2\n");
    let m = s(&plain_file(d.path(), "m.vcf", &t));
    let run = |threads: &str| {
        let o = gloci()
            .args([
                "vcf",
                "norm",
                &m,
                "-f",
                &r,
                "--split",
                "--skip-invalid",
                "-t",
                threads,
            ])
            .assert()
            .success();
        String::from_utf8(o.get_output().stderr.clone()).unwrap()
    };
    let (e1, e4) = (run("1"), run("4"));
    assert_eq!(e1, e4);
    assert_eq!(e1.matches("skipped").count(), 1, "{e1}");
    assert!(
        e1.contains("m.vcf: 2 records skipped (first: line 19:"),
        "{e1}"
    );
}

/// Records that alternate between two chromosomes normalise exactly as each
/// chromosome does alone (the reference keeps the last two windows, M6).
#[test]
fn alternating_chromosomes_give_the_same_records_as_each_alone() {
    let (d, r, _v) = work();
    let hdr = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=120>\n##contig=<ID=chr2,length=52>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let fa = std::fs::read_to_string(fixture("ref.fa")).unwrap();
    let base = |chrom: &str, pos: usize| -> char {
        let mut it = fa.split('>').filter(|x| x.starts_with(chrom));
        let seq: String = it.next().unwrap().lines().skip(1).collect();
        (seq.as_bytes()[pos - 1] as char).to_ascii_uppercase()
    };
    let alt = |c: char| if c == 'A' { 'C' } else { 'A' };
    let (mut mixed, mut c1, mut c2) = (String::new(), String::new(), String::new());
    for p in 1..=40usize {
        for (c, acc) in [("chr1", &mut c1), ("chr2", &mut c2)] {
            let b = base(c, p);
            if !"ACGT".contains(b) {
                continue;
            }
            let line = format!("{c}\t{p}\t.\t{b}\t{}\t.\t.\t.\n", alt(b));
            mixed.push_str(&line);
            acc.push_str(&line);
        }
    }
    let run = |name: &str, body: &str| {
        let f = s(&plain_file(d.path(), name, &format!("{hdr}{body}")));
        records(&out(&["vcf", "norm", &f, "-f", &r, "-t", "1"]))
    };
    let mut want = run("c1.vcf", &c1);
    want.extend(run("c2.vcf", &c2));
    let mut got = run("mixed.vcf", &mixed);
    // Same records; the window sort keeps chromosomes together or in input order.
    want.sort();
    got.sort();
    assert_eq!(got, want);
}

/// Under `--split --skip-invalid` a multi-allelic record with failing pieces
/// is skipped whole and counted once, and none of its pieces is written.
#[test]
fn a_record_whose_pieces_fail_is_skipped_whole_and_counted_once() {
    let (d, r, _v) = work();
    // chr1:1 is G in the reference: REF A is a mismatch for both pieces of
    // the first record; the second record splits into two good pieces.
    let t = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=120>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tA\tC,T\t.\t.\t.\n\
chr1\t2\t.\tA\tC,T\t.\t.\t.\nchr1\t3\t.\tT\tC\t.\t.\t.\n";
    let m = s(&plain_file(d.path(), "w.vcf", t));
    let o = gloci()
        .args([
            "vcf",
            "norm",
            &m,
            "-f",
            &r,
            "--split",
            "--skip-invalid",
            "-v",
        ])
        .assert()
        .success();
    let err = String::from_utf8(o.get_output().stderr.clone()).unwrap();
    let got = records(&String::from_utf8(o.get_output().stdout.clone()).unwrap());
    assert!(err.contains("1 record skipped (first: line 4:"), "{err}");
    assert!(
        err.contains("3 records read") && err.contains(", 1 skipped;"),
        "{err}"
    );
    assert_eq!(got.len(), 3, "{got:?}");
    assert!(got.iter().all(|l| !l.starts_with("chr1\t1\t")), "{got:?}");
}
