//! Regression tests for the final-review fix wave (groups A-G).
mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;
use std::path::Path;

fn w(d: &Path, name: &str, text: &str) -> String {
    plain_file(d, name, text).to_string_lossy().into_owned()
}

fn p(d: &Path, name: &str) -> String {
    d.join(name).to_string_lossy().into_owned()
}

// ---- Group A: sort/view format handling ----

const BG: &str = "chr1\t10\t20\t0.5\nchr1\t30\t40\t1.5\n";

#[test]
fn trackless_bedgraph_to_bedgraph_output() {
    let d = tempdir();
    let i = w(d.path(), "sig.bedgraph", BG);
    for cmd in ["sort", "view"] {
        let o = p(d.path(), &format!("{cmd}.out.bedgraph"));
        gloci().args([cmd, &i, "-o", &o]).assert().success();
        assert_eq!(std::fs::read_to_string(&o).unwrap(), BG, "{cmd}");
    }
}

#[test]
fn bed_input_to_bedgraph_output_and_back() {
    let d = tempdir();
    let i = w(d.path(), "x.bed", "chr1\t1\t2\n");
    let o = p(d.path(), "x.bedgraph");
    gloci().args(["sort", &i, "-o", &o]).assert().success();
    let bg = w(d.path(), "y.bedgraph", "chr1\t1\t2\t3\n");
    let o2 = p(d.path(), "y.bed");
    gloci().args(["view", &bg, "-o", &o2]).assert().success();
}

#[test]
fn sort_accepts_format_flag() {
    let d = tempdir();
    let i = w(
        d.path(),
        "sig.txt",
        "track type=bedGraph\nchr1\t30\t40\t1\nchr1\t10\t20\t2\n",
    );
    gloci()
        .args(["sort", &i, "--format", "bedgraph", "--header"])
        .assert()
        .success()
        .stdout("track type=bedGraph\nchr1\t10\t20\t2\nchr1\t30\t40\t1\n");
    gloci()
        .args(["sort", &i, "--format", "nope"])
        .assert()
        .code(2);
}

#[test]
fn family_mismatch_message_is_identical_for_sort_and_view() {
    let d = tempdir();
    let v = w(
        d.path(),
        "a.vcf",
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t5\t.\tA\tC\t.\t.\t.\n",
    );
    let o = p(d.path(), "x.bed");
    let msg = |cmd: &str| {
        let out = gloci()
            .args([cmd, &v, "-o", &o])
            .assert()
            .code(2)
            .get_output()
            .stderr
            .clone();
        String::from_utf8(out).unwrap()
    };
    let (s, vw) = (msg("sort"), msg("view"));
    assert_eq!(s, vw);
    assert!(
        s.contains("the input is VCF but the output format is BED"),
        "{s}"
    );
    assert!(s.contains("-O"), "{s}");
    assert!(!Path::new(&o).exists());
}

// ---- Group B: undeclared VCF contigs ----

const VCF_COLS: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
/// The header declares chr1 only; the records use chr2 as well.
fn undeclared_vcf() -> String {
    format!(
        "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n{VCF_COLS}\
chr1\t5\t.\tA\tC\t.\t.\t.\nchr2\t7\t.\tA\tC\t.\t.\t.\nchr2\t9\t.\tG\tT\t.\t.\t.\n"
    )
}
const WARN2: &str =
    "chromosome \"chr2\" is not in the ##contig lines; add ##contig=<ID=chr2> to the header";
const SORT_HINT: &str =
    "add ##contig=<ID=chr2> to the header, or pass -g with a genome that lists it";

fn records(text: &str) -> Vec<&str> {
    text.lines().filter(|l| !l.starts_with('#')).collect()
}

#[test]
fn view_plain_and_recompress_accept_undeclared_contigs() {
    let d = tempdir();
    let x = w(d.path(), "x.vcf", &undeclared_vcf());
    let o = p(d.path(), "x.vcf.gz");
    gloci()
        .args(["view", &x, "-o", &o, "--write-index"])
        .assert()
        .success()
        .stderr(contains(WARN2).and(contains("x.vcf line 5")))
        // one warning per contig, not per record
        .stderr(predicate::function(|s: &str| {
            s.matches("warning:").count() == 1
        }));
    assert!(d.path().join("x.vcf.gz.tbi").exists());
    let out = gloci()
        .args(["view", &o])
        .assert()
        .success()
        .get_output()
        .clone();
    let text = String::from_utf8(out.stdout).unwrap();
    assert_eq!(records(&text).len(), 3);
}

#[test]
fn view_quiet_silences_the_undeclared_warning() {
    let d = tempdir();
    let x = w(d.path(), "x.vcf", &undeclared_vcf());
    gloci()
        .args(["-q", "view", &x])
        .assert()
        .success()
        .stderr("");
}

#[test]
fn view_region_on_undeclared_contig_via_index_and_scan() {
    let d = tempdir();
    let x = w(d.path(), "x.vcf", &undeclared_vcf());
    let o = p(d.path(), "x.vcf.gz");
    gloci().args(["view", &x, "-o", &o]).assert().success();
    gloci().args(["index", &o]).assert().success();
    for file in [&o, &x] {
        let out = gloci()
            .args(["view", file, "chr2:8-9"])
            .assert()
            .success()
            .stderr(contains(WARN2))
            .get_output()
            .clone();
        let text = String::from_utf8(out.stdout).unwrap();
        let recs = records(&text);
        assert_eq!(recs.len(), 1, "{file}: {text}");
        assert!(recs[0].starts_with("chr2\t9\t"), "{recs:?}");
    }
    // the whole contig, and a declared one alongside it
    let out = gloci()
        .args(["view", &o, "chr1", "chr2"])
        .assert()
        .success()
        .get_output()
        .clone();
    assert_eq!(records(&String::from_utf8(out.stdout).unwrap()).len(), 3);
}

#[test]
fn view_assume_sorted_accepts_undeclared_contigs() {
    let d = tempdir();
    let x = w(d.path(), "x.vcf", &undeclared_vcf());
    let out = gloci()
        .args(["view", &x, "--assume-sorted", "chr2"])
        .assert()
        .success()
        .get_output()
        .clone();
    assert_eq!(records(&String::from_utf8(out.stdout).unwrap()).len(), 2);
}

#[test]
fn index_accepts_undeclared_contigs() {
    let d = tempdir();
    let x = w(d.path(), "x.vcf", &undeclared_vcf());
    let o = p(d.path(), "x.vcf.gz");
    gloci().args(["view", &x, "-o", &o]).assert().success();
    gloci().args(["index", &o]).assert().success();
    assert!(d.path().join("x.vcf.gz.tbi").exists());
}

#[test]
fn sort_check_and_merge_refuse_undeclared_contigs_with_a_hint() {
    let d = tempdir();
    let x = w(d.path(), "x.vcf", &undeclared_vcf());
    for args in [
        vec!["sort", x.as_str()],
        vec!["sort", "--check", x.as_str()],
        vec!["sort", "-m", x.as_str()],
    ] {
        gloci()
            .args(&args)
            .assert()
            .code(1)
            .stderr(contains(SORT_HINT).and(contains("x.vcf line 5")));
    }
    let o = p(d.path(), "s.vcf");
    gloci().args(["sort", &x, "-o", &o]).assert().code(1);
    assert!(!Path::new(&o).exists());
}

#[test]
fn sort_accepts_an_undeclared_contig_that_g_lists() {
    let d = tempdir();
    let x = w(d.path(), "x.vcf", &undeclared_vcf());
    let g = w(d.path(), "g.genome", "chr2\t1000\nchr1\t1000\n");
    let out = gloci()
        .args(["sort", &x, "-g", &g])
        .assert()
        .success()
        .get_output()
        .clone();
    let text = String::from_utf8(out.stdout).unwrap();
    let recs = records(&text);
    assert!(
        recs[0].starts_with("chr2\t7") && recs[2].starts_with("chr1\t5"),
        "{recs:?}"
    );
}

#[test]
fn bcf_output_refuses_undeclared_contigs() {
    let d = tempdir();
    let x = w(d.path(), "x.vcf", &undeclared_vcf());
    let o = p(d.path(), "x.bcf");
    gloci()
        .args(["view", &x, "-o", &o])
        .assert()
        .code(1)
        .stderr(contains("BCF needs every contig").and(contains("add ##contig=<ID=chr2>")));
    assert!(!Path::new(&o).exists());
    // a declared file still converts
    let ok = w(
        d.path(),
        "ok.vcf",
        &format!("##fileformat=VCFv4.3\n##contig=<ID=chr1>\n{VCF_COLS}chr1\t5\t.\tA\tC\t.\t.\t.\n"),
    );
    gloci().args(["view", &ok, "-o", &o]).assert().success();
}

// ---- Group C: multi-VCF sort header definitions ----

fn vcf_with(defs: &str, recs: &[&str]) -> String {
    let mut s = format!("##fileformat=VCFv4.3\n{defs}##contig=<ID=chr1,length=1000>\n{VCF_COLS}");
    for r in recs {
        s.push_str(r);
        s.push('\n');
    }
    s
}

const DP: &str = "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n";
const REC1: &str = "chr1\t5\t.\tA\tC\t.\t.\tDP=3";
const REC2: &str = "chr1\t2\t.\tA\tC\t.\t.\t.";

#[test]
fn multi_vcf_sort_refuses_a_definition_only_one_file_has() {
    let d = tempdir();
    let a = w(d.path(), "a.vcf", &vcf_with("", &[REC2]));
    let b = w(d.path(), "b.vcf", &vcf_with(DP, &[REC1]));
    let (an, bn) = (a.as_str(), b.as_str());
    gloci()
        .args(["sort", an, bn])
        .assert()
        .code(1)
        .stderr(contains(format!(
            "{b} defines ##INFO=<ID=DP> but {a} does not; give the files the same header definitions, or sort them separately"
        )));
    // the other way round names the file that has it too
    gloci()
        .args(["sort", bn, an])
        .assert()
        .code(1)
        .stderr(contains(format!(
            "{b} defines ##INFO=<ID=DP> but {a} does not"
        )));
    // -m is a multi-file path too, and --check does not merge anything
    gloci().args(["sort", "-m", an, bn]).assert().code(1);
    gloci().args(["sort", "--check", an, bn]).assert().success();
}

#[test]
fn multi_vcf_sort_refuses_conflicting_definitions() {
    let d = tempdir();
    let dp_float = "##INFO=<ID=DP,Number=1,Type=Float,Description=\"Depth\">\n";
    let a = w(d.path(), "a.vcf", &vcf_with(DP, &[REC1]));
    let b = w(d.path(), "b.vcf", &vcf_with(dp_float, &[REC2]));
    gloci()
        .args(["sort", &a, &b])
        .assert()
        .code(1)
        .stderr(contains("define ##INFO=<ID=DP> differently").and(contains("Type=Float")));
}

#[test]
fn multi_vcf_sort_treats_pass_as_implicit_and_ignores_other_meta_lines() {
    let d = tempdir();
    let pass = "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    let a = w(
        d.path(),
        "a.vcf",
        &vcf_with(&format!("{pass}{DP}"), &[REC1]),
    );
    let b = w(
        d.path(),
        "b.vcf",
        &vcf_with(&format!("##source=other\n{DP}"), &[REC2]),
    );
    let out = gloci()
        .args(["sort", &a, &b])
        .assert()
        .success()
        .get_output()
        .clone();
    let text = String::from_utf8(out.stdout).unwrap();
    assert_eq!(records(&text), vec![REC2, REC1]);
    // the first file's header is kept
    assert!(text.contains("##FILTER=<ID=PASS") && !text.contains("##source=other"));
}

#[test]
fn multi_vcf_sort_accepts_identical_headers() {
    let d = tempdir();
    let a = w(d.path(), "a.vcf", &vcf_with(DP, &[REC1]));
    let b = w(d.path(), "b.vcf", &vcf_with(DP, &[REC2]));
    gloci().args(["sort", &a, &b]).assert().success();
    gloci().args(["sort", "-m", &a, &b]).assert().success();
}

// ---- Group E: device and special-file outputs ----

#[cfg(unix)]
mod special_outputs {
    use super::*;

    const BED: &str = "chr1\t30\t40\nchr1\t10\t20\n";

    #[test]
    fn dev_null_succeeds_for_sort_and_view() {
        let d = tempdir();
        let i = w(d.path(), "x.bed", BED);
        for cmd in ["sort", "view"] {
            gloci()
                .args([cmd, &i, "-o", "/dev/null"])
                .assert()
                .success()
                .stdout("")
                .stderr("");
        }
        // an explicit -O still picks the format
        gloci()
            .args(["sort", &i, "-o", "/dev/null", "-O", "bed.gz"])
            .assert()
            .success();
    }

    #[test]
    fn dev_stdout_writes_the_records() {
        let d = tempdir();
        let i = w(d.path(), "x.bed", BED);
        gloci()
            .args(["sort", &i, "-o", "/dev/stdout"])
            .assert()
            .success()
            .stdout("chr1\t10\t20\nchr1\t30\t40\n");
        gloci()
            .args(["view", &i, "-o", "/dev/stdout"])
            .assert()
            .success()
            .stdout(BED);
    }

    #[test]
    fn dev_fd_path_is_written_directly() {
        // The shell hands gloci fd 3 on a regular file; /dev/fd/3 must not
        // get a temp file beside it.
        let d = tempdir();
        let i = w(d.path(), "x.bed", BED);
        let target = d.path().join("captured");
        let bin = assert_cmd::cargo::cargo_bin("gloci");
        let st = std::process::Command::new("sh")
            .arg("-c")
            .arg(r#""$0" sort "$1" -o /dev/fd/3 3>"$2""#)
            .arg(&bin)
            .arg(&i)
            .arg(&target)
            .status()
            .unwrap();
        assert!(st.success());
        assert_eq!(
            std::fs::read_to_string(&target).unwrap(),
            "chr1\t10\t20\nchr1\t30\t40\n"
        );
    }

    #[test]
    fn fifo_destination_works() {
        let d = tempdir();
        let i = w(d.path(), "x.bed", BED);
        let fifo = d.path().join("pipe");
        assert!(
            std::process::Command::new("mkfifo")
                .arg(&fifo)
                .status()
                .unwrap()
                .success()
        );
        let reader = {
            let fifo = fifo.clone();
            std::thread::spawn(move || std::fs::read_to_string(fifo).unwrap())
        };
        gloci()
            .args(["sort", &i, "-o", &fifo.to_string_lossy()])
            .assert()
            .success();
        assert_eq!(reader.join().unwrap(), "chr1\t10\t20\nchr1\t30\t40\n");
        // no temp files beside the FIFO, and the FIFO is still a FIFO
        let names: Vec<_> = std::fs::read_dir(d.path())
            .unwrap()
            .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        assert!(names.iter().all(|n| !n.contains(".tmp-")), "{names:?}");
    }

    #[test]
    fn write_index_to_a_device_is_a_usage_error() {
        let d = tempdir();
        let i = w(d.path(), "x.bed", BED);
        for cmd in ["sort", "view"] {
            gloci()
                .args([cmd, &i, "-o", "/dev/null", "-O", "bed.gz", "--write-index"])
                .assert()
                .code(2)
                .stderr(contains("regular output file").and(contains("/dev/null")));
        }
    }

    #[test]
    fn missing_output_directory_says_cannot_create() {
        let d = tempdir();
        let i = w(d.path(), "x.bed", BED);
        let o = p(d.path(), "nope/out.bed");
        for cmd in ["sort", "view"] {
            gloci().args([cmd, &i, "-o", &o]).assert().code(1).stderr(
                contains(format!("cannot create {o}"))
                    .and(contains("check that the directory exists and is writable"))
                    .and(contains("readable").not()),
            );
        }
    }

    #[test]
    fn index_staging_failure_says_cannot_create() {
        use std::os::unix::fs::PermissionsExt;
        let d = tempdir();
        let i = w(d.path(), "x.bed", "chr1\t1\t2\n");
        let gz = p(d.path(), "x.bed.gz");
        gloci().args(["view", &i, "-o", &gz]).assert().success();
        let dir = d.path();
        std::fs::set_permissions(dir, std::fs::Permissions::from_mode(0o555)).unwrap();
        let denied = std::fs::File::create(dir.join("probe")).is_err();
        let out = gloci().args(["index", &gz]).assert().get_output().clone();
        std::fs::set_permissions(dir, std::fs::Permissions::from_mode(0o755)).unwrap();
        if denied {
            // (skipped when running as root, which can write anywhere)
            let err = String::from_utf8(out.stderr).unwrap();
            assert!(
                err.contains("cannot create") && !err.contains("readable"),
                "{err}"
            );
            assert_eq!(out.status.code(), Some(1));
        }
    }
}

// ---- Group G: minors ----

#[test]
fn view_write_index_on_unsorted_input_names_the_input() {
    let d = tempdir();
    let i = w(d.path(), "u.bed", "chr1\t2\t3\nchr2\t1\t2\nchr1\t5\t6\n");
    let o = p(d.path(), "out.bed.gz");
    gloci()
        .args(["view", &i, "-o", &o, "--write-index"])
        .assert()
        .code(1)
        .stderr(
            contains(format!("{i} is not sorted"))
                .and(contains(format!("gloci sort {i}")))
                .and(contains("out.bed.gz").not()),
        );
    assert!(!Path::new(&o).exists());
}

#[test]
fn genome_path_like_values_say_no_such_file() {
    let d = tempdir();
    let i = w(d.path(), "x.bed", "chr1\t1\t2\n");
    for g in [
        "./missing.fai",
        "missing.genome",
        "dir/none",
        "ref.dict",
        "x.BAM",
    ] {
        gloci()
            .args(["sort", &i, "-g", g])
            .assert()
            .code(2)
            .stderr(contains(format!("no such file: {g}")).and(contains("unknown genome").not()));
    }
    // a plain word is still a genome id
    gloci()
        .args(["sort", &i, "-g", "hg39"])
        .assert()
        .code(2)
        .stderr(contains("unknown genome"));
}

#[test]
fn check_reports_a_bad_genome_once() {
    let d = tempdir();
    let a = w(d.path(), "a.bed", "chr1\t1\t2\n");
    let b = w(d.path(), "b.bed", "chr1\t1\t2\n");
    let out = gloci()
        .args(["check", &a, &b, "-g", "./missing.fai"])
        .assert()
        .code(2)
        .stdout("")
        .get_output()
        .clone();
    let err = String::from_utf8(out.stderr).unwrap();
    assert_eq!(err.matches("no such file").count(), 1, "{err}");
    assert!(!err.contains(&a) && !err.contains(&b), "{err}");
}

#[test]
fn index_plain_gzip_hint_rewrites_in_place_and_stdin_is_a_usage_error() {
    let d = tempdir();
    let gz = p(d.path(), "p.bed.gz");
    gloci_testkit::fixtures::write(
        d.path(),
        "p.bed.gz",
        &gloci_testkit::fixtures::gzip(b"chr1\t1\t2\n"),
    );
    gloci()
        .args(["index", &gz])
        .assert()
        .code(2)
        .stderr(contains(format!("gloci view {gz} -o {gz}")).and(contains(".gz.gz").not()));
    gloci()
        .args(["index", "-"])
        .assert()
        .code(2)
        .stderr(contains("cannot read standard input"));
}

#[test]
fn view_help_describes_what_is_dropped() {
    let out = gloci()
        .args(["view", "--help"])
        .assert()
        .success()
        .get_output()
        .clone();
    let help = String::from_utf8(out.stdout).unwrap();
    assert!(!help.contains("exactly as read"), "{help}");
    assert!(
        help.contains("blank lines") && help.contains("CRLF"),
        "{help}"
    );
}

#[test]
fn bcf_index_refusal_for_contiguous_contigs_out_of_dictionary_order() {
    let d = tempdir();
    let x = w(
        d.path(),
        "x.vcf",
        &format!(
            "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n##contig=<ID=chr2,length=1000>\n{VCF_COLS}\
chr2\t7\t.\tA\tC\t.\t.\t.\nchr1\t5\t.\tA\tC\t.\t.\t.\n"
        ),
    );
    let bcf = p(d.path(), "x.bcf");
    gloci().args(["view", &x, "-o", &bcf]).assert().success();
    gloci().args(["index", &bcf]).assert().code(1).stderr(
        contains("not sorted")
            .and(contains(format!("gloci sort {bcf}")))
            .and(contains("contig dictionary")),
    );
    assert!(!Path::new(&format!("{bcf}.csi")).exists());
}

// ---- Group G10: output is byte-identical for every --threads value ----

/// `n` pseudo-random BED lines over chr1..chr4 (deterministic).
fn random_bed(n: usize, seed: u64) -> String {
    let mut x = seed;
    let mut next = move || {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        x >> 33
    };
    let mut s = String::new();
    for i in 0..n {
        let chrom = 1 + next() % 4;
        let start = next() % 5_000_000;
        let len = 1 + next() % 2_000;
        s.push_str(&format!(
            "chr{chrom}\t{start}\t{}\tf{i}\t{}\t+\n",
            start + len,
            next() % 1000
        ));
    }
    s
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

/// Run `args` with -t 1 and -t 4, to stdout and (when `file` names an output
/// path template) to a BGZF file; every result must be identical.
fn assert_thread_invariant(d: &Path, args: &[&str]) {
    let run = |t: &str| {
        let mut v: Vec<&str> = args.to_vec();
        v.extend(["-t", t]);
        stdout_bytes(&v)
    };
    let (one, four) = (run("1"), run("4"));
    assert!(!one.is_empty(), "{args:?} produced no output");
    assert_eq!(one, four, "{args:?}: stdout differs between -t 1 and -t 4");
    // and as a BGZF file
    let write = |t: &str, name: &str| {
        let o = p(d, name);
        let mut v: Vec<&str> = args.to_vec();
        v.extend(["-t", t, "-o", &o]);
        gloci().args(&v).assert().success();
        std::fs::read(&o).unwrap()
    };
    assert_eq!(
        write("1", "t1.bed.gz"),
        write("4", "t4.bed.gz"),
        "{args:?}: BGZF output differs between -t 1 and -t 4"
    );
}

#[test]
fn sort_merge_is_identical_across_threads() {
    let d = tempdir();
    let a = w(d.path(), "a.in.bed", &random_bed(60_000, 1));
    let b = w(d.path(), "b.in.bed", &random_bed(60_000, 2));
    let (sa, sb) = (p(d.path(), "a.bed"), p(d.path(), "b.bed"));
    gloci().args(["sort", &a, "-o", &sa]).assert().success();
    gloci().args(["sort", &b, "-o", &sb]).assert().success();
    assert_thread_invariant(d.path(), &["sort", "-m", &sa, &sb]);
}

#[test]
fn view_regions_are_identical_across_threads() {
    let d = tempdir();
    let raw = w(d.path(), "raw.bed", &random_bed(120_000, 3));
    let plain = p(d.path(), "x.bed");
    let gz = p(d.path(), "x.bed.gz");
    gloci()
        .args(["sort", &raw, "-o", &plain])
        .assert()
        .success();
    gloci()
        .args(["sort", &raw, "-o", &gz, "--write-index"])
        .assert()
        .success();
    let regions = ["chr1:100000-2000000", "chr3", "chr2:4000000-"];
    // indexed
    let mut args = vec!["view", gz.as_str()];
    args.extend(regions);
    assert_thread_invariant(d.path(), &args);
    // unindexed scan
    let mut args = vec!["view", plain.as_str()];
    args.extend(regions);
    assert_thread_invariant(d.path(), &args);
    // binary search
    let mut args = vec!["view", plain.as_str(), "--assume-sorted"];
    args.extend(regions);
    assert_thread_invariant(d.path(), &args);
    // all three agree with each other
    let by = |extra: &[&str], file: &str| {
        let mut v = vec!["view", file];
        v.extend(regions);
        v.extend(extra);
        stdout_bytes(&v)
    };
    let indexed = by(&[], &gz);
    assert_eq!(indexed, by(&[], &plain));
    assert_eq!(indexed, by(&["--assume-sorted"], &plain));
}
