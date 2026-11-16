//! `--sort` on `merge` and `complement` (D39): the result equals
//! `gloci sort INPUTS | gloci <cmd> -`, byte for byte, within `--max-mem`.
mod common;
use common::*;
use predicates::str::contains;
use std::path::Path;

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

/// Deterministic pseudo-random BED lines, `n` of them, over a few chromosomes.
fn shuffled_bed(seed: u64, n: usize) -> String {
    let mut x = seed | 1;
    let mut next = move || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    };
    let mut out = String::new();
    for i in 0..n {
        let c = ["chr1", "chr2", "chr10", "chrX", "chrM"][(next() % 5) as usize];
        let start = next() % 100_000;
        let len = next() % 300;
        out += &format!(
            "{c}\t{start}\t{}\tr{i}\t{}\t{}\n",
            start + len,
            next() % 100,
            ["+", "-", "."][(next() % 3) as usize]
        );
    }
    out
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

fn files(dir: &Path) -> Vec<String> {
    vec![
        s(&plain_file(dir, "a.bed", &shuffled_bed(1, 400))),
        s(&plain_file(dir, "b.bed", &shuffled_bed(2, 300))),
        s(&plain_file(dir, "c.bed", &shuffled_bed(3, 200))),
    ]
}

/// `gloci sort INPUTS > sorted.bed; gloci CMD ARGS sorted.bed`.
fn via_sort_pipe(dir: &Path, inputs: &[String], cmd: &[&str]) -> Vec<u8> {
    let mut sort_args = vec!["sort"];
    sort_args.extend(inputs.iter().map(String::as_str));
    // The same chromosome order as the command: -g goes to both.
    if let Some(i) = cmd.iter().position(|a| *a == "-g") {
        sort_args.extend(&cmd[i..=i + 1]);
    }
    let sorted = plain_file(dir, "sorted.bed", "");
    std::fs::write(&sorted, run(&sort_args)).unwrap();
    let mut args = cmd.to_vec();
    let sorted = s(&sorted);
    args.insert(1, &sorted);
    run(&args)
}

#[test]
fn merge_sort_equals_sort_pipe() {
    let dir = tempdir();
    let inputs = files(dir.path());
    for extra in [
        vec![],
        vec!["-d", "50"],
        vec!["-s", "--count", "-c", "5:max"],
    ] {
        let mut cmd = vec!["merge"];
        cmd.extend(&extra);
        let want = via_sort_pipe(dir.path(), &inputs, &cmd);
        assert!(!want.is_empty());
        let mut args = vec!["merge", "--sort"];
        args.extend(inputs.iter().map(String::as_str));
        args.extend(&extra);
        assert_eq!(run(&args), want, "merge {extra:?}");
    }
}

#[test]
fn complement_sort_equals_sort_pipe() {
    let dir = tempdir();
    let inputs = files(dir.path());
    let g = s(&plain_file(
        dir.path(),
        "g.genome",
        "chrX\t200000\nchr10\t200000\nchr2\t200000\nchr1\t200000\nchrM\t200000\n",
    ));
    for extra in [vec![], vec!["-g", g.as_str()]] {
        let mut cmd = vec!["complement"];
        cmd.extend(&extra);
        let want = via_sort_pipe(dir.path(), &inputs, &cmd);
        assert!(!want.is_empty());
        let mut args = vec!["complement", "--sort"];
        args.extend(inputs.iter().map(String::as_str));
        args.extend(&extra);
        assert_eq!(run(&args), want, "complement {extra:?}");
    }
}

#[test]
fn sort_flag_on_sorted_input_changes_nothing() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t1\t5\nchr1\t9\t12\n",
    ));
    assert_eq!(run(&["merge", &f, "--sort"]), run(&["merge", &f]));
    assert_eq!(run(&["complement", &f, "--sort"]), run(&["complement", &f]));
}

/// 330k shuffled BED lines (about 9 MB): with `--max-mem 16M` the in-memory
/// chunk (record bytes plus 32 bytes of index per record) fills up, so the
/// sorter spills at least once.
fn spill_bed() -> String {
    let mut s = String::with_capacity(10 << 20);
    let mut x: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut next = move || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    };
    for _ in 0..330_000 {
        let c = next() % 24 + 1;
        let start = next() % 200_000_000;
        s.push_str(&format!("chr{c}\t{start}\t{}\n", start + next() % 5000 + 1));
    }
    s
}

#[test]
fn sort_spills_within_max_mem() {
    let dir = tempdir();
    let big = s(&plain_file(dir.path(), "big.bed", &spill_bed()));
    let tmp = dir.path().join("spill");
    std::fs::create_dir(&tmp).unwrap();
    // The reference: sort once with `gloci sort`, then run the command on it.
    let sorted = s(&plain_file(dir.path(), "sorted.bed", ""));
    std::fs::write(&sorted, run(&["sort", &big])).unwrap();
    for cmd in ["merge", "complement"] {
        let small = gloci()
            .args([cmd, "--sort", "--max-mem", "16M", "-v", &big, "--tmp-dir"])
            .arg(&tmp)
            .assert()
            .success();
        let small = small.get_output().clone();
        let log = String::from_utf8_lossy(&small.stderr);
        assert!(
            log.contains("chunks spilled") && !log.contains(" 0 chunks spilled"),
            "{cmd} should have spilled: {log}"
        );
        let want = run(&[cmd, &sorted]);
        assert!(!want.is_empty());
        assert_eq!(small.stdout, want, "{cmd}");
        // Nothing is left in the temp directory.
        assert_eq!(std::fs::read_dir(&tmp).unwrap().count(), 0, "{cmd}");
    }
}

#[test]
fn past_end_under_sort_exits_1_and_cleans_up() {
    let dir = tempdir();
    let f = s(&plain_file(
        dir.path(),
        "a.bed",
        "chr1\t900\t1100\nchr1\t0\t10\n",
    ));
    let g = s(&plain_file(dir.path(), "g.genome", "chr1\t1000\n"));
    let tmp = dir.path().join("spill");
    std::fs::create_dir(&tmp).unwrap();
    gloci()
        .args(["complement", "--sort", &f, "-g", &g, "--tmp-dir"])
        .arg(&tmp)
        .assert()
        .code(1)
        .stderr(contains("extends past the end of chr1"))
        // Line numbers of the sorted stream are not the user's: records.
        .stderr(contains("the sorted inputs ("))
        .stderr(contains("record 2"));
    assert_eq!(std::fs::read_dir(&tmp).unwrap().count(), 0);
}

#[test]
fn merge_sort_determinism() {
    let dir = tempdir();
    let inputs = files(dir.path());
    let run_t = |t: &str| {
        let mut args = vec!["merge", "--sort", "-s", "--count", "-c", "5:mean", "-t", t];
        args.extend(inputs.iter().map(String::as_str));
        run(&args)
    };
    let base = run_t("1");
    assert!(!base.is_empty());
    assert_eq!(run_t("2"), base);
    assert_eq!(run_t("8"), base);
}

const VCF_HEAD: &str = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n\
##contig=<ID=chr2,length=100000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// Shuffled VCF records, deletions among them (their allele span starts
/// after POS-1, so allele and ref spans differ).
fn shuffled_vcf(seed: u64, n: usize) -> String {
    let mut x = seed | 1;
    let mut next = move || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        x
    };
    let mut out = VCF_HEAD.to_string();
    for _ in 0..n {
        let c = ["chr1", "chr2"][(next() % 2) as usize];
        let pos = next() % 2000 + 1;
        let (r, a) = match next() % 3 {
            0 => ("A", "G"),
            1 => ("ATTTT", "A"),
            _ => ("A", "ATT"),
        };
        out += &format!("{c}\t{pos}\t.\t{r}\t{a}\t.\t.\t.\n");
    }
    out
}

fn sort_flag_vs_pipe_vcf(ext: &str, convert: bool) {
    let dir = tempdir();
    let mut inputs = Vec::new();
    for (i, seed) in [11u64, 12].iter().enumerate() {
        let vcf = plain_file(dir.path(), &format!("in{i}.vcf"), &shuffled_vcf(*seed, 300));
        if convert {
            let bcf = dir.path().join(format!("in{i}.bcf"));
            gloci()
                .args(["view", &s(&vcf), "-o", &s(&bcf)])
                .assert()
                .success();
            inputs.push(s(&bcf));
        } else {
            inputs.push(s(&vcf));
        }
    }
    let sorted = s(&dir.path().join(format!("sorted.{ext}")));
    let mut sort_args = vec!["sort"];
    sort_args.extend(inputs.iter().map(String::as_str));
    sort_args.extend(["-o", &sorted]);
    gloci().args(&sort_args).assert().success();
    for cmd in ["merge", "complement"] {
        for span in ["allele", "ref"] {
            let want = run(&[cmd, &sorted, "--vcf-span", span]);
            assert!(!want.is_empty(), "{cmd} {span}");
            let mut args = vec![cmd, "--sort", "--vcf-span", span];
            args.extend(inputs.iter().map(String::as_str));
            assert_eq!(run(&args), want, "{ext} {cmd} --vcf-span {span}");
        }
    }
}

#[test]
fn vcf_sort_equals_sort_pipe_with_both_spans() {
    sort_flag_vs_pipe_vcf("vcf", false);
}

#[test]
fn bcf_sort_equals_sort_pipe() {
    sort_flag_vs_pipe_vcf("bcf", true);
}

/// `--max-mem` means something only with `--sort`; `--tmp-dir` is accepted
/// everywhere and ignored when nothing spills.
#[test]
fn max_mem_needs_sort_tmp_dir_does_not() {
    let dir = tempdir();
    let f = s(&plain_file(dir.path(), "a.bed", "chr1\t1\t5\nchr1\t3\t9\n"));
    for cmd in ["merge", "complement"] {
        gloci().args([cmd, &f, "--max-mem", "32M"]).assert().code(2);
        let want = gloci().args([cmd, &f]).assert().success();
        let want = want.get_output().stdout.clone();
        gloci()
            .args([cmd, &f, "--tmp-dir"])
            .arg(dir.path())
            .assert()
            .success()
            .stdout(want);
    }
}

#[test]
fn max_mem_too_small_is_usage_error() {
    let dir = tempdir();
    let f = s(&plain_file(dir.path(), "a.bed", "chr1\t1\t5\n"));
    gloci()
        .args(["merge", &f, "--sort", "--max-mem", "1M"])
        .assert()
        .code(2)
        .stderr(contains("at least 16M"));
}

#[test]
fn mixed_formats_with_sort_is_usage_error() {
    let dir = tempdir();
    let bed = s(&plain_file(dir.path(), "a.bed", "chr1\t1\t5\n"));
    let vcf = s(&plain_file(
        dir.path(),
        "a.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t5\t.\tA\tG\t.\t.\t.\n",
    ));
    gloci()
        .args(["merge", &bed, &vcf, "--sort"])
        .assert()
        .code(2)
        .stderr(contains("one format"));
}

/// The record number in a data error counts data records, not the header
/// lines of the sorted temp file or of the input.
#[test]
fn past_end_record_number_ignores_header_lines() {
    let dir = tempdir();
    let g = s(&plain_file(dir.path(), "g.genome", "chr1\t1000\n"));
    // A 3-line VCF header; the second record (in sorted order) is past the end.
    let vcf = s(&plain_file(
        dir.path(),
        "b.vcf",
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t1500\t.\tA\tG\t.\t.\t.\nchr1\t10\t.\tA\tG\t.\t.\t.\n",
    ));
    // A BED with track and comment lines before the data.
    let bed = s(&plain_file(
        dir.path(),
        "b.bed",
        "track name=x\n# comment\nchr1\t900\t1100\nchr1\t0\t10\n",
    ));
    for f in [&vcf, &bed] {
        gloci()
            .args(["complement", "--sort", f, "-g", &g])
            .assert()
            .code(1)
            .stderr(contains(format!("the sorted inputs ({f}) record 2:")));
    }
}
