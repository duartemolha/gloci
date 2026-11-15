//! Sorting, large fixtures and output formats of the BED and VCF commands.

use crate::common::*;
use crate::fixtures_bed_vcf::*;
use crate::helpers::*;
use std::fmt::Write as _;
use std::io::Read;

/// `--sort` over a budget smaller than the input spills to disk and gives the
/// bytes of the in-memory sort. (The input is asserted to be more than twice
/// the budget, so the spill is not optional; the spill directory is empty
/// again afterwards.)
#[test]
fn p2_sort_spill_matches_in_memory() {
    // --sort over a budget smaller than the input spills to disk and gives
    // the bytes of the in-memory sort.
    let fx = P2::new();
    let big = fx.plain("big.bed", &bed_text(1_600_000));
    assert!(std::fs::metadata(&big).unwrap().len() > 2 * (16 << 20));
    let tmp = tempdir();
    let spill = [
        "cluster",
        &big,
        "--sort",
        "--max-mem",
        "16M",
        "--tmp-dir",
        tmp.path().to_str().unwrap(),
    ];
    let spilled: Vec<_> = ["1", "8"].iter().map(|t| stdout_with(&spill, t)).collect();
    assert_eq!(spilled[0], spilled[1], "cluster --sort spilled");
    assert_eq!(
        spilled[0],
        stdout_with(&["cluster", &big, "--sort"], "2"),
        "spilled vs in-memory"
    );
    assert_eq!(
        std::fs::read_dir(tmp.path()).unwrap().count(),
        0,
        "spill files left"
    );
}

/// Fixtures large enough to be split: 60k records is more than one 16384-record
/// batch, the BGZF copies span dozens of blocks, and chr1's 40k records make a
/// groupby run that crosses batch boundaries. Output at `-t 1` must equal
/// `-t 8`.
#[test]
fn p2_determinism_large_fixture() {
    let fx = P2::new();
    let genome = fx.plain("g.genome", "chr1\t1000000\nchr2\t1000000\n");
    let mut a = String::new();
    let mut b = String::new();
    for (chrom, n) in [("chr1", 40_000usize), ("chr2", 20_000)] {
        for i in 0..n {
            let st = 100 + i * 23 + (i % 5) * 3;
            writeln!(
                a,
                "{chrom}\t{st}\t{}\ta{i}\t{}\t{}",
                st + 10 + (i % 7) * 9,
                i % 100,
                if i % 3 == 0 { '-' } else { '+' }
            )
            .unwrap();
            if i % 2 == 0 {
                let bs = 105 + i * 31;
                writeln!(b, "{chrom}\t{bs}\t{}\tb{i}\t1\t+", bs + 40).unwrap();
            }
        }
    }
    let a_plain = fx.plain("a.bed", &a);
    let a_gz = fx.bgzf("a.bed.gz", &a);
    let a_idx = fx.indexed("a_idx.bed.gz", &a);
    let b_plain = fx.plain("b.bed", &b);
    let b_gz = fx.bgzf("b.bed.gz", &b);
    assert!(
        std::fs::metadata(&a_gz).unwrap().len() > 3 * 65_536,
        "several BGZF blocks"
    );
    let run = |args: &[&str], t: &str| -> Vec<u8> {
        let mut c = gloci();
        c.args(args).args(["-t", t]);
        let o = c.output().unwrap();
        assert_eq!(
            o.status.code(),
            Some(0),
            "{args:?}: {}",
            String::from_utf8_lossy(&o.stderr)
        );
        o.stdout
    };
    for (name, tpl) in p2_bed_cmds(&genome) {
        // per-base output is one line per genome base: too long to repeat
        if name == "genomecov --per-base" {
            continue;
        }
        let sub = |ai: &str, bi: &str| -> Vec<String> {
            tpl.iter()
                .map(|x| match *x {
                    "{a}" => ai.to_string(),
                    "{b}" => bi.to_string(),
                    o => o.to_string(),
                })
                .collect()
        };
        let plain_args = sub(&a_plain, &b_plain);
        let refs: Vec<&str> = plain_args.iter().map(String::as_str).collect();
        let base = run(&refs, "1");
        // the default genomecov is a short histogram; the others are long
        let want = if name == "genomecov" { 50 } else { 100_000 };
        assert!(base.len() > want, "{name}: small output");
        assert_eq!(base, run(&refs, "8"), "{name}: -t 8");
        for (ai, bi) in [(&a_gz, &b_gz), (&a_idx, &b_gz)] {
            let args = sub(ai, bi);
            let refs: Vec<&str> = args.iter().map(String::as_str).collect();
            assert_eq!(base, run(&refs, "8"), "{name}: {ai} -t 8");
        }
    }
    // the groupby run over chr1 crosses the batch boundary and is one line
    let gb =
        String::from_utf8(run(&["groupby", &a_plain, "-g", "1", "-c", "4:count"], "8")).unwrap();
    assert_eq!(gb, "chr1\t40000\nchr2\t20000\n");

    // the same for a VCF of 40k records
    let mut v = String::from(P2_VCF_HEAD);
    for i in 0..40_000usize {
        writeln!(
            v,
            "chr1\t{}\t.\tA\tC\t{}\tPASS\tDP={}\tGT:DP\t0/1:{}\t1/1:3\t./.:4",
            100 + i * 7,
            20 + i % 50,
            i % 30,
            i % 9
        )
        .unwrap();
    }
    let v_plain = fx.plain("v.vcf", &v);
    let v_gz = fx.bgzf("v.vcf.gz", &v);
    let v_idx = fx.indexed("v_idx.vcf.gz", &v);
    for tpl in [
        vec!["cluster", "-d", "3"],
        vec!["spacing"],
        vec!["vcf", "extract", "--info", "DP"],
        vec!["vcf", "extract", "--format", "GT,DP"],
    ] {
        let with = |f: &str| -> Vec<String> {
            let mut a: Vec<String> = tpl.iter().map(|s| s.to_string()).collect();
            a.insert(if tpl[0] == "vcf" { 2 } else { 1 }, f.to_string());
            a
        };
        let args = with(&v_plain);
        let refs: Vec<&str> = args.iter().map(String::as_str).collect();
        let base = run(&refs, "1");
        assert!(base.len() > 100_000);
        assert_eq!(base, run(&refs, "8"), "{tpl:?}");
        for f in [&v_gz, &v_idx] {
            let args = with(f);
            let refs: Vec<&str> = args.iter().map(String::as_str).collect();
            assert_eq!(base, run(&refs, "8"), "{tpl:?} {f}");
        }
    }
    // vcf concat: a second file that starts after the first ends, so the
    // 40k-record first file spans several batches and BGZF blocks
    let mut tail = String::from(P2_VCF_HEAD);
    for i in 0..3_000usize {
        writeln!(
            tail,
            "chr1\t{}\t.\tA\tC\t{}\tPASS\tDP={}\tGT:DP\t0/1:{}\t1/1:3\t./.:4",
            1_000_000 + i * 7,
            20 + i % 50,
            i % 30,
            i % 9
        )
        .unwrap();
    }
    let tail_plain = fx.plain("tail.vcf", &tail);
    let tail_gz = fx.bgzf("tail.vcf.gz", &tail);
    for (first, second) in [
        (&v_plain, &tail_plain),
        (&v_gz, &tail_gz),
        (&v_idx, &tail_gz),
    ] {
        let args = ["vcf", "concat", first.as_str(), second.as_str()];
        let base = run(&args, "1");
        assert!(base.len() > 100_000, "concat output");
        assert_eq!(base, run(&args, "8"), "concat {first}");
    }
    // vcf stats files
    let stats = |f: &str, t: &str| -> Vec<Vec<u8>> {
        let d = tempdir();
        let prefix = s_of(&d.path().join("o"));
        gloci()
            .args([
                "vcf",
                "stats",
                f,
                "--out-prefix",
                &prefix,
                "--freq",
                "--het",
                "--site-quality",
                "-t",
                t,
            ])
            .assert()
            .success();
        ["frq", "het", "lqual"]
            .iter()
            .map(|e| std::fs::read(format!("{prefix}.{e}")).unwrap())
            .collect()
    };
    let base = stats(&v_plain, "1");
    assert_eq!(base[0].iter().filter(|&&c| c == b'\n').count(), 40_001);
    assert_eq!(base, stats(&v_plain, "8"));
    assert_eq!(base, stats(&v_gz, "8"));
    assert_eq!(base, stats(&v_idx, "8"));
}

/// `vcf concat` writes the format `-O` or the `-o` extension names; the first
/// input's format is only the default.
#[test]
fn p2_concat_output_format() {
    let fx = P2::new();
    let (t1, t2) = (p2_vcf(1..=1, 5), p2_vcf(2..=2, 5));
    let (v1, v2) = (fx.plain("1.vcf", &t1), fx.plain("2.vcf", &t2));
    let (b1, b2) = (fx.bcf("1.bcf", &t1), fx.bcf("2.bcf", &t2));
    let (z1, z2) = (fx.bgzf("1.vcf.gz", &t1), fx.bgzf("2.vcf.gz", &t2));
    let d = tempdir();
    let to = |name: &str, inputs: [&str; 2], extra: &[&str]| -> Vec<u8> {
        let out = s_of(&d.path().join(name));
        let mut full = vec!["vcf", "concat", inputs[0], inputs[1]];
        full.extend_from_slice(extra);
        full.extend(["-o", &out]);
        gloci().args(&full).assert().success();
        std::fs::read(&out).unwrap()
    };
    // -o extension decides, whatever the inputs are
    for inputs in [[&v1, &v2], [&z1, &z2], [&b1, &b2]] {
        let inputs = [inputs[0].as_str(), inputs[1].as_str()];
        assert!(is_bcf(&to("x.bcf", inputs, &[])), "{inputs:?} -o x.bcf");
        assert!(
            is_bcf(&to("y.out", inputs, &["-O", "bcf"])),
            "{inputs:?} -O bcf"
        );
        let gz = to("x.vcf.gz", inputs, &[]);
        assert!(!is_bcf(&gz), "{inputs:?} -o x.vcf.gz is VCF");
        let mut text = Vec::new();
        noodles_bgzf::io::Reader::new(&gz[..])
            .read_to_end(&mut text)
            .unwrap();
        assert!(
            text.starts_with(b"##fileformat=VCF"),
            "{inputs:?}: BGZF VCF text"
        );
        // -O vcf is plain text even for BCF input
        let plain = to("z.out", inputs, &["-O", "vcf"]);
        assert!(plain.starts_with(b"##fileformat=VCF"), "{inputs:?} -O vcf");
    }
    // the default (stdout) follows the first input
    let o = gloci().args(["vcf", "concat", &b1, &b2]).output().unwrap();
    assert!(o.status.success());
    assert!(is_bcf(&o.stdout), "BCF in, BCF out by default");
    let o = gloci().args(["vcf", "concat", &v1, &v2]).output().unwrap();
    assert!(o.stdout.starts_with(b"##fileformat=VCF"));
}
