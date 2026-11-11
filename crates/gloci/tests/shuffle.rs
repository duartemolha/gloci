//! `gloci shuffle` (D152, D153): placement rules, unplaceable records, BEDPE, thread
//! independence and errors. bedtools parity of the forced cases is in the golden cases
//! (`shuffle-forced-*-vs-bedtools`); its random placements differ by design
//! (docs/differences.md#shuffle-placement, #rng-stream).
mod common;
use common::*;
use predicates::str::contains;
use std::time::Instant;

const G: &str = "chr1\t100\nchr2\t200\nchr3\t700\n";

fn genome(dir: &std::path::Path) -> String {
    plain_file(dir, "g.genome", G)
        .to_string_lossy()
        .into_owned()
}

fn run(args: &[&str]) -> String {
    let o = gloci().arg("shuffle").args(args).assert().success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

fn fields(line: &str) -> Vec<&str> {
    line.split('\t').collect()
}

fn bed(n: usize) -> String {
    (0..n)
        .map(|i| format!("chr{}\t{}\t{}\tr{i}\t{i}\t-\n", i % 3 + 1, i, i + 5))
        .collect()
}

#[test]
fn shuffle_pinned() {
    let d = tempdir();
    let g = genome(d.path());
    let f = plain_file(
        d.path(),
        "a.bed",
        "chr1\t10\t20\ta\t5\t+\nchr2\t0\t50\tb\t6\t-\n",
    );
    let out = run(&[f.to_str().unwrap(), "-g", &g, "--seed", "42"]);
    assert_eq!(out, "chr3\t56\t66\ta\t5\t+\nchr3\t460\t510\tb\t6\t-\n");
}

#[test]
fn shuffle_keeps_columns_and_input_order() {
    let d = tempdir();
    let g = genome(d.path());
    let f = plain_file(d.path(), "a.bed", &bed(300));
    let out = run(&[f.to_str().unwrap(), "-g", &g, "--seed", "3"]);
    assert_eq!(out.lines().count(), 300);
    for (i, line) in out.lines().enumerate() {
        let f = fields(line);
        assert_eq!(&f[3..], [format!("r{i}"), format!("{i}"), "-".into()]);
        let (s, e): (u64, u64) = (f[1].parse().unwrap(), f[2].parse().unwrap());
        assert_eq!(e - s, 5);
        let len = match f[0] {
            "chr1" => 100,
            "chr2" => 200,
            "chr3" => 700,
            c => panic!("{c}"),
        };
        assert!(e <= len, "{line}");
    }
}

#[test]
fn shuffle_same_chrom_and_chrom_first() {
    let d = tempdir();
    let g = genome(d.path());
    let f = plain_file(d.path(), "a.bed", &bed(300));
    let f = f.to_str().unwrap();
    let out = run(&[f, "-g", &g, "--seed", "3", "--same-chrom"]);
    for (i, line) in out.lines().enumerate() {
        assert_eq!(fields(line)[0], format!("chr{}", i % 3 + 1));
    }
    let out = run(&[f, "-g", &g, "--seed", "3", "--chrom-first"]);
    let c1 = out.lines().filter(|l| l.starts_with("chr1\t")).count();
    // One third of 300 on the shortest chromosome (by length alone it would be ~12).
    assert!(c1 > 60, "chr1 got {c1}");
    gloci()
        .args(["shuffle", f, "-g", &g, "--same-chrom", "--chrom-first"])
        .assert()
        .code(2);
}

#[test]
fn shuffle_exclude_leaves_no_room() {
    let d = tempdir();
    let g = plain_file(d.path(), "g.genome", "chr1\t100\n");
    let bedf = plain_file(
        d.path(),
        "input.bed",
        "chr1\t0\t10\tx\t1\t+\nchr1\t20\t23\ty\t2\t+\n",
    );
    let ex = plain_file(d.path(), "ex.bed", "chr1\t0\t95\n");
    let out = d.path().join("out.bed");
    gloci()
        .args(["shuffle", bedf.to_str().unwrap(), "-g", g.to_str().unwrap()])
        .args(["--exclude", ex.to_str().unwrap(), "--seed", "1", "-o"])
        .arg(&out)
        .assert()
        .code(1)
        .stderr(contains("input.bed line 1"))
        .stderr(contains("length 10"))
        .stderr(contains("--exclude"))
        .stderr(contains("--skip-unplaceable"));
    assert!(!out.exists(), "no output file");
    gloci()
        .args(["shuffle", bedf.to_str().unwrap(), "-g", g.to_str().unwrap()])
        .args([
            "--exclude",
            ex.to_str().unwrap(),
            "--seed",
            "1",
            "--skip-unplaceable",
        ])
        .assert()
        .success()
        .stdout(contains("\ty\t2\t+"))
        .stderr(contains("1 record could not be placed"))
        .stderr(contains("input.bed line 1 length 10"));
}

#[test]
fn shuffle_whole_record_stays_in_include() {
    let d = tempdir();
    let g = genome(d.path());
    let f = plain_file(d.path(), "a.bed", &"chr1\t0\t8\tx\n".repeat(200));
    // Include intervals of 10 bp: an 8-bp record must fit inside one (bedtools would let it
    // start inside and run out).
    let inc = plain_file(
        d.path(),
        "inc.bed",
        "chr1\t10\t20\nchr2\t50\t60\nchr3\t0\t10\n",
    );
    let out = run(&[
        f.to_str().unwrap(),
        "-g",
        &g,
        "--include",
        inc.to_str().unwrap(),
        "--seed",
        "4",
    ]);
    let ok = [("chr1", 10, 20), ("chr2", 50, 60), ("chr3", 0, 10)];
    for line in out.lines() {
        let f = fields(line);
        let (s, e): (u64, u64) = (f[1].parse().unwrap(), f[2].parse().unwrap());
        assert!(
            ok.iter().any(|&(c, a, b)| c == f[0] && a <= s && e <= b),
            "{line}"
        );
    }
}

#[test]
fn shuffle_include_shorter_than_feature() {
    let d = tempdir();
    let g = genome(d.path());
    let f = plain_file(d.path(), "input.bed", "chr1\t0\t4\nchr1\t0\t9\n");
    let inc = plain_file(d.path(), "inc.bed", "chr1\t10\t15\nchr2\t0\t5\n");
    gloci()
        .args([
            "shuffle",
            f.to_str().unwrap(),
            "-g",
            &g,
            "--seed",
            "1",
            "--include",
        ])
        .arg(&inc)
        .assert()
        .code(1)
        .stderr(contains("input.bed line 2"))
        .stderr(contains("length 9"))
        .stderr(contains("--include"));
}

#[test]
fn shuffle_no_overlapping() {
    let d = tempdir();
    let g = plain_file(d.path(), "g.genome", "chr1\t100\nchr2\t100\n");
    let f = plain_file(d.path(), "a.bed", &"chr1\t0\t10\n".repeat(8));
    let out = run(&[
        f.to_str().unwrap(),
        "-g",
        g.to_str().unwrap(),
        "--no-overlapping",
        "--seed",
        "2",
    ]);
    let mut got: Vec<(String, u64, u64)> = out
        .lines()
        .map(|l| {
            let f = fields(l);
            (
                f[0].to_string(),
                f[1].parse().unwrap(),
                f[2].parse().unwrap(),
            )
        })
        .collect();
    assert_eq!(got.len(), 8);
    got.sort();
    for w in got.windows(2) {
        assert!(w[0].0 != w[1].0 || w[0].2 <= w[1].1, "{w:?}");
    }
}

#[test]
fn shuffle_no_overlapping_max_tries() {
    let d = tempdir();
    let g = plain_file(d.path(), "g.genome", "chr1\t100\n");
    let f = plain_file(d.path(), "input.bed", "chr1\t0\t60\nchr1\t0\t60\n");
    gloci()
        .args(["shuffle", f.to_str().unwrap(), "-g", g.to_str().unwrap()])
        .args(["--no-overlapping", "--max-tries", "5", "--seed", "1"])
        .assert()
        .code(1)
        .stderr(contains("input.bed line 2"))
        .stderr(contains("5 draws"))
        .stderr(contains("--max-tries"));
}

#[test]
fn shuffle_vcf_refused() {
    let d = tempdir();
    let g = genome(d.path());
    let v = plain_file(
        d.path(),
        "x.vcf",
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t5\t.\tA\tC\t.\t.\t.\n",
    );
    gloci()
        .args(["shuffle", v.to_str().unwrap(), "-g", &g, "--seed", "1"])
        .assert()
        .code(2)
        .stderr(contains("gloci convert"));
}

#[test]
fn shuffle_gff_is_rewritten_one_based() {
    let d = tempdir();
    let g = plain_file(d.path(), "g.genome", "chrA\t10\n");
    // chrB's feature has length 10 and chrA is 10 bases: the only place is the whole chromosome.
    let f = plain_file(
        d.path(),
        "a.gff",
        "chrB\tsrc\tgene\t3\t12\t.\t+\t.\tID=g1\n",
    );
    let out = run(&[
        f.to_str().unwrap(),
        "-g",
        g.to_str().unwrap(),
        "--seed",
        "1",
    ]);
    assert_eq!(out, "chrA\tsrc\tgene\t1\t10\t.\t+\t.\tID=g1\n");
}

#[test]
fn shuffle_header_and_gzip_output() {
    let d = tempdir();
    let g = genome(d.path());
    let f = plain_file(d.path(), "a.bed", "#note\ntrack name=x\nchr1\t0\t5\n");
    let f = f.to_str().unwrap();
    let out = run(&[f, "-g", &g, "--seed", "1", "--header"]);
    assert!(out.starts_with("#note\ntrack name=x\nchr"), "{out}");
    let out = run(&[f, "-g", &g, "--seed", "1"]);
    assert!(out.starts_with("chr"), "{out}");
    let gz = d.path().join("o.bed.gz");
    gloci()
        .args(["shuffle", f, "-g", &g, "--seed", "1", "-o"])
        .arg(&gz)
        .assert()
        .success();
    assert_eq!(&std::fs::read(&gz).unwrap()[..2], [0x1f, 0x8b]);
}

#[test]
fn shuffle_refuses_write_index_up_front() {
    let d = tempdir();
    let g = genome(d.path());
    let f = plain_file(d.path(), "a.bed", "chr1\t0\t5\n");
    let out = d.path().join("x.bed.gz");
    gloci()
        .args([
            "shuffle",
            f.to_str().unwrap(),
            "-g",
            &g,
            "--seed",
            "1",
            "--write-index",
            "-o",
        ])
        .arg(&out)
        .assert()
        .code(2)
        .stderr(contains("shuffle output is unsorted"))
        .stderr(contains("gloci sort"));
    assert!(!out.exists());
}

#[test]
fn shuffle_same_bytes_at_any_thread_count() {
    let d = tempdir();
    let g = genome(d.path());
    let f = plain_file(d.path(), "a.bed", &bed(20_000));
    let small = plain_file(d.path(), "small.bed", &bed(40));
    let f = f.to_str().unwrap();
    let small = small.to_str().unwrap();
    // Packing 20,000 records without overlap into 1,000 bases is impossible: that path
    // gets a small input.
    let (out, _) = (
        run(&[
            small,
            "-g",
            &g,
            "--seed",
            "9",
            "--no-overlapping",
            "-t",
            "1",
        ]),
        0,
    );
    for t in ["2", "8"] {
        assert!(out == run(&[small, "-g", &g, "--seed", "9", "--no-overlapping", "-t", t]));
    }
    for extra in [&[][..], &["--chrom-first"], &["--same-chrom"]] {
        let mut args = vec![f, "-g", &g, "--seed", "9", "-t", "1"];
        args.extend_from_slice(extra);
        let base = run(&args);
        for t in ["2", "8"] {
            let mut args = vec![f, "-g", &g, "--seed", "9", "-t", t];
            args.extend_from_slice(extra);
            assert!(base == run(&args), "-t {t} {extra:?} differs from -t 1");
        }
    }
}

#[test]
fn seed_warning_reproduces() {
    let d = tempdir();
    let g = genome(d.path());
    let f = plain_file(d.path(), "a.bed", &bed(50));
    let f = f.to_str().unwrap();
    let o = gloci().args(["shuffle", f, "-g", &g]).assert().success();
    let first = String::from_utf8(o.get_output().stdout.clone()).unwrap();
    let err = String::from_utf8(o.get_output().stderr.clone()).unwrap();
    let seed = err
        .split("using --seed ")
        .nth(1)
        .and_then(|s| s.split(' ').next())
        .unwrap_or_else(|| panic!("no seed in {err:?}"))
        .to_string();
    assert_eq!(first, run(&[f, "-g", &g, "--seed", &seed]));
}

#[test]
fn shuffle_huge_contig_count() {
    let d = tempdir();
    let mut text = String::new();
    for i in 0..200_000 {
        text.push_str(&format!("ctg{i}\t{}\n", 1000 + i % 7));
    }
    let g = plain_file(d.path(), "many.genome", &text);
    let f = plain_file(d.path(), "a.bed", &bed(10_000).replace("chr", "ctg"));
    let t = Instant::now();
    let o = gloci()
        .args(["shuffle", f.to_str().unwrap(), "-g"])
        .arg(&g)
        .args(["--seed", "5", "--same-chrom"])
        .assert()
        .success();
    assert_eq!(o.get_output().stdout.split(|&b| b == b'\n').count(), 10_001);
    let o = gloci()
        .args(["shuffle", f.to_str().unwrap(), "-g"])
        .arg(&g)
        .args(["--seed", "5", "--chrom-first"])
        .assert()
        .success();
    assert_eq!(o.get_output().stdout.split(|&b| b == b'\n').count(), 10_001);
    assert!(t.elapsed().as_secs_f64() < 20.0, "took {:?}", t.elapsed());
}

// --- BEDPE (D153) ---

const BP: &str = "chr1\t10\t20\tchr1\t50\t80\tp1\t7\t+\t-\n\
chr1\t10\t20\tchr2\t30\t45\tp2\t8\t-\t+\textra\n\
chr3\t100\t110\t.\t-1\t-1\tp3\t9\t+\t.\n";

fn bedpe(args: &[&str]) -> String {
    let o = gloci()
        .arg("shuffle")
        .args(args)
        .args(["--format", "bedpe"])
        .assert()
        .success();
    String::from_utf8(o.get_output().stdout.clone()).unwrap()
}

#[test]
fn shuffle_bedpe() {
    let d = tempdir();
    let g = genome(d.path());
    let f = plain_file(d.path(), "a.bedpe", BP);
    let f = f.to_str().unwrap();
    for seed in 0..60u64 {
        let out = bedpe(&[f, "-g", &g, "--seed", &seed.to_string()]);
        let rows: Vec<Vec<&str>> = out.lines().map(fields).collect();
        assert_eq!(rows.len(), 3);
        // A pair on one chromosome keeps its end lengths and its gap.
        let r = &rows[0];
        assert_eq!(r[0], r[3]);
        let n = |i: usize| r[i].parse::<u64>().unwrap();
        assert_eq!(n(2) - n(1), 10);
        assert_eq!(n(5) - n(4), 30);
        assert_eq!(n(4) - n(1), 40);
        assert_eq!(&r[6..], ["p1", "7", "+", "-"]);
        // Two chromosomes: each end placed on its own, lengths and columns kept.
        let r = &rows[1];
        let n = |i: usize| r[i].parse::<u64>().unwrap();
        assert_eq!(n(2) - n(1), 10);
        assert_eq!(n(5) - n(4), 15);
        assert_eq!(&r[6..], ["p2", "8", "-", "+", "extra"]);
        // An absent end stays absent.
        assert_eq!(rows[2][3..6], [".", "-1", "-1"]);
        assert_eq!(rows[2][6..], ["p3", "9", "+", "."]);
    }
}

#[test]
fn shuffle_bedpe_ends_are_independent_across_chromosomes() {
    let d = tempdir();
    let g = genome(d.path());
    let f = plain_file(d.path(), "a.bedpe", "chr1\t10\t20\tchr2\t30\t45\n");
    let f = f.to_str().unwrap();
    let mut same = 0;
    for seed in 0..200u64 {
        let out = bedpe(&[f, "-g", &g, "--seed", &seed.to_string()]);
        let r = fields(out.lines().next().unwrap())
            .into_iter()
            .map(String::from)
            .collect::<Vec<_>>();
        if r[0] == r[3] {
            same += 1;
        }
    }
    // Independent ends land on one chromosome by chance (about 54% for 100/200/700).
    assert!(same > 50 && same < 190, "{same} of 200 on one chromosome");
}

#[test]
fn shuffle_bedpe_exclude_keeps_both_ends_clear() {
    let d = tempdir();
    let g = genome(d.path());
    let f = plain_file(d.path(), "a.bedpe", BP);
    // Everything of chr3 above 300 and all of chr1 is excluded.
    let ex = plain_file(d.path(), "ex.bed", "chr1\t0\t100\nchr3\t300\t700\n");
    for seed in 0..40u64 {
        let out = bedpe(&[
            f.to_str().unwrap(),
            "-g",
            &g,
            "--exclude",
            ex.to_str().unwrap(),
            "--seed",
            &seed.to_string(),
        ]);
        for line in out.lines() {
            let f = fields(line);
            for (c, s, e) in [(0, 1, 2), (3, 4, 5)] {
                if f[c] == "." {
                    continue;
                }
                let (s, e): (u64, u64) = (f[s].parse().unwrap(), f[e].parse().unwrap());
                assert!(f[c] != "chr1", "{line}");
                assert!(f[c] != "chr3" || e <= 300, "{line}");
                assert!(s < e);
            }
        }
    }
}

#[test]
fn shuffle_bedpe_flags_are_refused() {
    let d = tempdir();
    let g = genome(d.path());
    let f = plain_file(d.path(), "a.bedpe", BP);
    let inc = plain_file(d.path(), "inc.bed", "chr1\t0\t100\n");
    for extra in [
        vec!["--include", inc.to_str().unwrap()],
        vec!["--same-chrom"],
        vec!["--chrom-first"],
        vec!["--no-overlapping"],
    ] {
        let mut args = vec![
            "shuffle",
            f.to_str().unwrap(),
            "-g",
            &g,
            "--format",
            "bedpe",
        ];
        args.extend(extra.iter());
        gloci()
            .args(&args)
            .assert()
            .code(2)
            .stderr(contains(extra[0]));
    }
}

#[test]
fn shuffle_bedpe_unplaceable_and_output_format() {
    let d = tempdir();
    let g = plain_file(d.path(), "g.genome", "chr1\t20\n");
    let f = plain_file(
        d.path(),
        "input.bedpe",
        "chr1\t0\t5\tchr1\t10\t18\nchr1\t0\t5\tchr1\t10\t40\n",
    );
    gloci()
        .args([
            "shuffle",
            f.to_str().unwrap(),
            "-g",
            g.to_str().unwrap(),
            "--format",
            "bedpe",
            "--seed",
            "1",
        ])
        .assert()
        .code(1)
        .stderr(contains("input.bedpe line 2"))
        .stderr(contains("--skip-unplaceable"));
    let out = d.path().join("o.bedpe.gz");
    gloci()
        .args([
            "shuffle",
            f.to_str().unwrap(),
            "-g",
            g.to_str().unwrap(),
            "--format",
            "bedpe",
            "--seed",
            "1",
            "--skip-unplaceable",
            "-o",
        ])
        .arg(&out)
        .assert()
        .success()
        .stderr(contains("1 record could not be placed"));
    assert_eq!(&std::fs::read(&out).unwrap()[..2], [0x1f, 0x8b]);
    gloci()
        .args([
            "shuffle",
            f.to_str().unwrap(),
            "-g",
            g.to_str().unwrap(),
            "--format",
            "bedpe",
            "-o",
        ])
        .arg(d.path().join("o.bed"))
        .assert()
        .code(2)
        .stderr(contains("BEDPE"));
}

#[test]
fn shuffle_bedpe_pinned() {
    let d = tempdir();
    let g = genome(d.path());
    let f = plain_file(d.path(), "a.bedpe", BP);
    let out = bedpe(&[f.to_str().unwrap(), "-g", &g, "--seed", "42"]);
    assert_eq!(
        out,
        "chr3\t114\t124\tchr3\t154\t184\tp1\t7\t+\t-\nchr3\t473\t483\tchr2\t107\t122\tp2\t8\t-\t+\textra\n\
chr2\t130\t140\t.\t-1\t-1\tp3\t9\t+\t.\n"
    );
}

#[test]
fn shuffle_chrom_map_applies_to_include_and_exclude() {
    let d = tempdir();
    let g = genome(d.path());
    let map = plain_file(d.path(), "map.tsv", "one\tchr1\ntwo\tchr2\n");
    let f = plain_file(d.path(), "a.bed", &"one\t0\t5\n".repeat(50));
    let inc = plain_file(d.path(), "inc.bed", "two\t100\t110\n");
    let out = run(&[
        f.to_str().unwrap(),
        "-g",
        &g,
        "--chrom-map",
        map.to_str().unwrap(),
        "--include",
        inc.to_str().unwrap(),
        "--seed",
        "1",
    ]);
    for line in out.lines() {
        let f = fields(line);
        let (s, e): (u64, u64) = (f[1].parse().unwrap(), f[2].parse().unwrap());
        assert!(f[0] == "chr2" && s >= 100 && e <= 110, "{line}");
    }
    assert_eq!(out.lines().count(), 50);
    // The exclude file uses the source names too: chr1 [0, 95) is off limits.
    let g1 = plain_file(d.path(), "g1.genome", "chr1\t100\n");
    let ex = plain_file(d.path(), "ex.bed", "one\t0\t95\n");
    let out = run(&[
        f.to_str().unwrap(),
        "-g",
        g1.to_str().unwrap(),
        "--chrom-map",
        map.to_str().unwrap(),
        "--exclude",
        ex.to_str().unwrap(),
        "--seed",
        "1",
    ]);
    assert_eq!(out.lines().next().unwrap(), "chr1\t95\t100");
    // BEDPE too.
    let p = plain_file(d.path(), "a.bedpe", "one\t0\t5\t.\t-1\t-1\n");
    let out = bedpe(&[
        p.to_str().unwrap(),
        "-g",
        g1.to_str().unwrap(),
        "--chrom-map",
        map.to_str().unwrap(),
        "--exclude",
        ex.to_str().unwrap(),
        "--seed",
        "1",
    ]);
    assert_eq!(out, "chr1\t95\t100\t.\t-1\t-1\n");
}

#[test]
fn shuffle_bedpe_exclude_warns_about_unknown_chromosomes() {
    let d = tempdir();
    let g = genome(d.path());
    let f = plain_file(d.path(), "a.bedpe", BP);
    let ex = plain_file(d.path(), "ex.bed", "chrZ\t0\t10\nchr1\t0\t5\n");
    gloci()
        .args([
            "shuffle",
            f.to_str().unwrap(),
            "-g",
            &g,
            "--format",
            "bedpe",
            "--seed",
            "1",
            "--exclude",
        ])
        .arg(&ex)
        .assert()
        .success()
        .stderr(contains("1 of the 2 regions"))
        .stderr(contains("not in -g"));
}

#[test]
fn shuffle_bedpe_refuses_flags_before_creating_output() {
    let d = tempdir();
    let g = genome(d.path());
    let f = plain_file(d.path(), "a.bedpe", BP);
    let out = d.path().join("o.bedpe");
    gloci()
        .args([
            "shuffle",
            f.to_str().unwrap(),
            "-g",
            &g,
            "--format",
            "bedpe",
            "--same-chrom",
            "-o",
        ])
        .arg(&out)
        .assert()
        .code(2)
        .stderr(contains("--same-chrom"));
    let mut left: Vec<String> = std::fs::read_dir(d.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    left.sort();
    assert_eq!(left, ["a.bedpe", "g.genome"], "no output or temp file");
}

#[test]
fn shuffle_reads_wig_as_intervals_and_writes_bedgraph() {
    let d = tempdir();
    let g = plain_file(d.path(), "g.genome", "chr1\t5\n");
    let f = plain_file(
        d.path(),
        "a.wig",
        "track type=wiggle_0\nfixedStep chrom=chr1 start=11 step=5 span=5\n1.5\n2.5\n",
    );
    // Every 5-base interval fills the 5-base chromosome: one place each.
    let out = run(&[
        f.to_str().unwrap(),
        "-g",
        g.to_str().unwrap(),
        "--seed",
        "1",
    ]);
    assert_eq!(out, "chr1\t0\t5\t1.5\nchr1\t0\t5\t2.5\n");
}
