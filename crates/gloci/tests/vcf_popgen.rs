//! Tests for `gloci vcf popgen` (decisions D158 to D161): every file matches
//! the vcftools 0.1.16 file committed in `gloci-vcf/tests/golden/popgen` (made
//! by its `regen.sh`), one pass equals separate runs, the thread count and the
//! input format change nothing, the windowed statistics need sorted input, the
//! population files of Fst follow their rules, and the sites LD holds stay
//! within --max-mem.
mod common;
use common::*;
use predicates::str::contains;
use std::path::{Path, PathBuf};

fn golden_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../gloci-vcf/tests/golden/popgen")
}

fn golden(name: &str) -> String {
    std::fs::read_to_string(golden_dir().join(name)).unwrap_or_else(|e| panic!("{name}: {e}"))
}

fn fixture() -> String {
    golden_dir()
        .join("popgen.vcf")
        .to_string_lossy()
        .into_owned()
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

/// A file of the golden directory, as a path argument.
fn gpath(name: &str) -> String {
    s(&golden_dir().join(name))
}

/// The statistics vcftools was run for (regen.sh): gloci's flags, the tag
/// of the golden's prefix and the file's extension.
const STATS: [(&[&str], &str, &str); 6] = [
    (&["--site-pi"], "", "sites.pi"),
    (&["--window-pi", "10"], ".w10", "windowed.pi"),
    (
        &["--window-pi", "10", "--window-pi-step", "5"],
        ".w10s5",
        "windowed.pi",
    ),
    (
        &["--window-pi", "10", "--window-pi-step", "3"],
        ".w10s3",
        "windowed.pi",
    ),
    (&["--tajima-d", "49"], ".t49", "Tajima.D"),
    (&["--tajima-d", "100"], ".t100", "Tajima.D"),
];

/// The masking and sample selection of each run (regen.sh's MASK).
const MASKS: [(&str, &[&str]); 3] = [
    ("", &[]),
    (".dp5", &["--min-dp", "5"]),
    (".indv", &["--samples", "s2,s3"]),
];

/// Run `vcf popgen FILE --out-prefix DIR/o ARGS`.
fn popgen(file: &str, args: &[&str]) -> (tempfile::TempDir, String) {
    let d = tempdir();
    let prefix = s(&d.path().join("o"));
    gloci()
        .args(["vcf", "popgen", file, "--out-prefix", &prefix])
        .args(args)
        .assert()
        .success();
    (d, prefix)
}

fn read(prefix: &str, ext: &str) -> Vec<u8> {
    std::fs::read(format!("{prefix}.{ext}")).unwrap_or_else(|e| panic!("{prefix}.{ext}: {e}"))
}

#[test]
fn every_file_matches_vcftools() {
    for (mask, mask_args) in MASKS {
        for (flags, tag, ext) in STATS {
            let mut args = flags.to_vec();
            args.extend_from_slice(mask_args);
            let (_d, prefix) = popgen(&fixture(), &args);
            let name = format!("popgen{mask}{tag}.{ext}");
            assert_eq!(
                String::from_utf8(read(&prefix, ext)).unwrap(),
                golden(&name),
                "{name} ({args:?})"
            );
        }
    }
}

/// Every statistic at once: the files vcftools wrote in separate runs.
#[test]
fn popgen_one_pass_equals_separate_runs() {
    for (mask, mask_args) in MASKS {
        let mut args = vec![
            "--site-pi",
            "--window-pi",
            "10",
            "--window-pi-step",
            "5",
            "--tajima-d",
            "49",
            "--threads",
            "8",
        ];
        args.extend_from_slice(mask_args);
        let (d, prefix) = popgen(&fixture(), &args);
        for (tag, ext) in [
            ("", "sites.pi"),
            (".w10s5", "windowed.pi"),
            (".t49", "Tajima.D"),
        ] {
            assert_eq!(
                String::from_utf8(read(&prefix, ext)).unwrap(),
                golden(&format!("popgen{mask}{tag}.{ext}")),
                "{mask} {ext}"
            );
        }
        assert_eq!(std::fs::read_dir(d.path()).unwrap().count(), 3, "{mask}");
    }
}

/// The same bytes at -t 1, 2 and 8, from VCF, BGZF VCF and BCF.
#[test]
fn popgen_output_is_the_same_at_every_thread_count() {
    let d = tempdir();
    let text = std::fs::read_to_string(fixture()).unwrap();
    let gz = s(&bgzf_file(d.path(), "p.vcf.gz", &text));
    let bcf = s(&d.path().join("p.bcf"));
    gloci()
        .args(["vcf", "filter", &fixture(), "--vcf-span", "ref", "-o", &bcf])
        .assert()
        .success();
    let all = [
        "--site-pi",
        "--window-pi",
        "10",
        "--window-pi-step",
        "3",
        "--tajima-d",
        "49",
        "--min-dp",
        "5",
    ];
    let files = |file: &str, t: &str| -> Vec<Vec<u8>> {
        let mut args = all.to_vec();
        args.extend(["-t", t]);
        let (_d, prefix) = popgen(file, &args);
        ["sites.pi", "windowed.pi", "Tajima.D"]
            .iter()
            .map(|e| read(&prefix, e))
            .collect()
    };
    let base = files(&fixture(), "1");
    assert_eq!(
        String::from_utf8(base[2].clone()).unwrap(),
        golden("popgen.dp5.t49.Tajima.D")
    );
    for (file, t) in [
        (fixture(), "2"),
        (fixture(), "8"),
        (gz.clone(), "1"),
        (gz, "8"),
        (bcf.clone(), "1"),
        (bcf, "8"),
    ] {
        assert_eq!(files(&file, t), base, "{file} -t {t}");
    }

    // Fst reads every member's genotype, with every other statistic beside it.
    let fst_text = std::fs::read_to_string(golden_dir().join("fst.vcf")).unwrap();
    let fst_gz = s(&bgzf_file(d.path(), "f.vcf.gz", &fst_text));
    let fst_bcf = s(&d.path().join("f.bcf"));
    gloci()
        .args([
            "vcf",
            "filter",
            &gpath("fst.vcf"),
            "--vcf-span",
            "ref",
            "-o",
            &fst_bcf,
        ])
        .assert()
        .success();
    let (p1, p2, p3) = (gpath("pop1.txt"), gpath("pop2.txt"), gpath("pop3.txt"));
    let fst_args = [
        "--site-pi",
        "--window-pi",
        "100",
        "--tajima-d",
        "100",
        "--fst-pop",
        &p1,
        "--fst-pop",
        &p2,
        "--fst-pop",
        &p3,
        "--fst-window-size",
        "100",
        "--fst-window-step",
        "30",
        "--min-dp",
        "5",
    ];
    let exts = [
        "sites.pi",
        "windowed.pi",
        "Tajima.D",
        "weir.fst",
        "weir.fst.summary",
        "windowed.weir.fst",
    ];
    let fst_files = |file: &str, t: &str| -> Vec<Vec<u8>> {
        let mut args = fst_args.to_vec();
        args.extend(["-t", t]);
        let (_d, prefix) = popgen(file, &args);
        exts.iter().map(|e| read(&prefix, e)).collect()
    };
    let base = fst_files(&gpath("fst.vcf"), "1");
    assert!(base[5].len() > 200, "the windows were written");
    for (file, t) in [
        (gpath("fst.vcf"), "2"),
        (gpath("fst.vcf"), "8"),
        (fst_gz.clone(), "1"),
        (fst_gz, "8"),
        (fst_bcf.clone(), "2"),
        (fst_bcf, "8"),
    ] {
        assert_eq!(fst_files(&file, t), base, "{file} -t {t}");
    }

    // LD holds and pairs the sites on the main thread, with pi beside it.
    for (fixture, args, exts) in [
        (
            "ld.vcf",
            &[
                "--geno-r2",
                "--interchrom-geno-r2",
                "--ld-window-bp",
                "150",
                "--min-dp",
                "5",
                "--site-pi",
            ][..],
            &["geno.ld", "interchrom.geno.ld", "sites.pi"][..],
        ),
        (
            "phased.vcf",
            &[
                "--hap-r2",
                "--geno-r2",
                "--interchrom-hap-r2",
                "--ld-window",
                "3",
            ],
            &["hap.ld", "geno.ld", "interchrom.hap.ld"],
        ),
    ] {
        let text = std::fs::read_to_string(golden_dir().join(fixture)).unwrap();
        let gz = s(&bgzf_file(d.path(), &format!("{fixture}.gz"), &text));
        let bcf = s(&d.path().join(format!("{fixture}.bcf")));
        gloci()
            .args([
                "vcf",
                "filter",
                &gpath(fixture),
                "--vcf-span",
                "ref",
                "-o",
                &bcf,
            ])
            .assert()
            .success();
        let files = |file: &str, t: &str| -> Vec<Vec<u8>> {
            let mut a = args.to_vec();
            a.extend(["-t", t]);
            let (_d, prefix) = popgen(file, &a);
            exts.iter().map(|e| read(&prefix, e)).collect()
        };
        let base = files(&gpath(fixture), "1");
        assert!(
            base.iter().all(|f| f.len() > 100),
            "{fixture}: rows written"
        );
        for (file, t) in [
            (gpath(fixture), "2"),
            (gpath(fixture), "8"),
            (gz.clone(), "1"),
            (gz, "8"),
            (bcf.clone(), "2"),
            (bcf, "8"),
        ] {
            assert_eq!(files(&file, t), base, "{file} -t {t}");
        }
    }
}

const UNSORTED: &str = "##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2
1\t500\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/0
2\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/0
1\t200\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/0
";

const BACKWARDS: &str = "##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2
1\t500\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/0
1\t499\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/0
";

/// Windowed statistics stream their rows, so unsorted input is exit 1 with
/// the fix and no file (decision D158, docs/differences.md#popgen-sorted-input);
/// --site-pi alone takes it. A step above the window and a missing statistic
/// are usage errors (#popgen-window-step).
#[test]
fn popgen_unsorted_window_stats() {
    let d = tempdir();
    let unsorted = s(&plain_file(d.path(), "u.vcf", UNSORTED));
    let backwards = s(&plain_file(d.path(), "b.vcf", BACKWARDS));
    for (file, why) in [
        (&unsorted, "the chromosome reappears at 1:200"),
        (&backwards, "a position goes backwards at 1:499"),
    ] {
        for stat in [&["--window-pi", "10"][..], &["--tajima-d", "10"]] {
            let out = tempdir();
            let prefix = s(&out.path().join("o"));
            gloci()
                .args(["vcf", "popgen", file, "--out-prefix", &prefix, "--site-pi"])
                .args(stat)
                .assert()
                .code(1)
                .stderr(contains(why))
                .stderr(contains("gloci sort"));
            assert_eq!(
                std::fs::read_dir(out.path()).unwrap().count(),
                0,
                "{stat:?}: a failed run writes no file"
            );
        }
    }
    let (_d, prefix) = popgen(&unsorted, &["--site-pi"]);
    assert_eq!(
        String::from_utf8(read(&prefix, "sites.pi")).unwrap(),
        "CHROM\tPOS\tPI\n1\t500\t0.5\n2\t100\t0.5\n1\t200\t0.5\n"
    );

    let usage = |args: &[&str], msg: &str| {
        let out = tempdir();
        let prefix = s(&out.path().join("o"));
        gloci()
            .args(["vcf", "popgen", &unsorted, "--out-prefix", &prefix])
            .args(args)
            .assert()
            .code(2)
            .stderr(contains(msg));
        assert_eq!(std::fs::read_dir(out.path()).unwrap().count(), 0);
    };
    usage(
        &["--window-pi", "10", "--window-pi-step", "20"],
        "--window-pi-step 20 is larger than the window",
    );
    usage(
        &["--window-pi", "0"],
        "--window-pi needs a window size above 0",
    );
    usage(
        &["--window-pi", "10", "--window-pi-step", "0"],
        "--window-pi-step needs a step above 0",
    );
    usage(&["--tajima-d", "0"], "--tajima-d needs a bin size above 0");
    usage(&["--window-pi-step", "5", "--site-pi"], "--window-pi");
    usage(
        &[],
        "add at least one of --site-pi, --window-pi SIZE, --tajima-d BIN, --fst-pop FILE",
    );
}

/// Every statistic is taken over the selected samples, so a run without one
/// is a usage error (vcftools stops too).
#[test]
fn popgen_needs_a_sample() {
    let d = tempdir();
    let sites_only = s(&plain_file(
        d.path(),
        "s.vcf",
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t5\t.\tA\tG\t.\t.\t.\n",
    ));
    let prefix = s(&d.path().join("o"));
    gloci()
        .args([
            "vcf",
            "popgen",
            &sites_only,
            "--out-prefix",
            &prefix,
            "--site-pi",
        ])
        .assert()
        .code(2)
        .stderr(contains("no sample"));
    assert!(!Path::new(&format!("{prefix}.sites.pi")).exists());
}

/// One Fst run of regen.sh.
struct FstRun {
    /// The prefix of its goldens.
    tag: &'static str,
    file: &'static str,
    pops: &'static [&'static str],
    /// gloci's other flags.
    args: &'static [&'static str],
    /// For a windowed run, the prefix of the run without the window, whose
    /// `.weir.fst` gloci writes too.
    base: Option<&'static str>,
}

const P2: &[&str] = &["pop1.txt", "pop2.txt"];
const P3: &[&str] = &["pop1.txt", "pop2.txt", "pop3.txt"];

const FST_RUNS: [FstRun; 9] = [
    FstRun {
        tag: "fst",
        file: "fst.vcf",
        pops: P2,
        args: &[],
        base: None,
    },
    FstRun {
        tag: "fst.p3",
        file: "fst.vcf",
        pops: P3,
        args: &[],
        base: None,
    },
    FstRun {
        tag: "fst.dp5",
        file: "fst.vcf",
        pops: P2,
        args: &["--min-dp", "5"],
        base: None,
    },
    FstRun {
        tag: "fst.w100",
        file: "fst.vcf",
        pops: P2,
        args: &["--fst-window-size", "100"],
        base: Some("fst"),
    },
    FstRun {
        tag: "fst.w100s50",
        file: "fst.vcf",
        pops: P2,
        args: &["--fst-window-size", "100", "--fst-window-step", "50"],
        base: Some("fst"),
    },
    FstRun {
        tag: "fst.p3.w100s50",
        file: "fst.vcf",
        pops: P3,
        args: &["--fst-window-size", "100", "--fst-window-step", "50"],
        base: Some("fst.p3"),
    },
    FstRun {
        tag: "fst.dp5.w100s50",
        file: "fst.vcf",
        pops: P2,
        args: &[
            "--min-dp",
            "5",
            "--fst-window-size",
            "100",
            "--fst-window-step",
            "50",
        ],
        base: Some("fst.dp5"),
    },
    FstRun {
        tag: "missingpop",
        file: "missingpop.vcf",
        pops: P2,
        args: &[],
        base: None,
    },
    FstRun {
        tag: "missingpop.w100",
        file: "missingpop.vcf",
        pops: P2,
        args: &["--fst-window-size", "100"],
        base: Some("missingpop"),
    },
];

/// `vcf popgen FILE --fst-pop ...` for one of the runs: the prefix.
fn fst_run(file: &str, pops: &[&str], args: &[&str]) -> (tempfile::TempDir, String) {
    let mut all: Vec<String> = Vec::new();
    for p in pops {
        all.push("--fst-pop".into());
        all.push(gpath(p));
    }
    all.extend(args.iter().map(|a| a.to_string()));
    let refs: Vec<&str> = all.iter().map(String::as_str).collect();
    popgen(&gpath(file), &refs)
}

/// C++'s `ostream << x` at precision 5, as vcftools logs the estimates.
fn g5(x: f64) -> String {
    let mut b = Vec::new();
    gloci_core::cfmt::write_g(&mut b, x, 5);
    String::from_utf8(b).unwrap()
}

/// The summary against what vcftools logged: N_SITES is the sites of its
/// `.weir.fst` that are not nan, and the two estimates, which gloci prints
/// to 6 digits, are vcftools' when written to the 5 digits of its log.
fn check_summary(summary: &str, log: &str, weir: &str, name: &str) {
    let rows: Vec<&str> = summary.lines().collect();
    assert_eq!(rows.len(), 2, "{name}: {summary}");
    assert_eq!(rows[0], "N_SITES\tMEAN_FST\tWEIGHTED_FST", "{name}");
    let cols: Vec<&str> = rows[1].split('\t').collect();
    let sites = weir
        .lines()
        .skip(1)
        .filter(|l| !l.ends_with("\tnan"))
        .count();
    assert_eq!(cols[0], sites.to_string(), "{name}: N_SITES");
    let logged = |what: &str| -> String {
        let line = log
            .lines()
            .find(|l| l.contains(&format!("{what} Fst estimate: ")))
            .unwrap_or_else(|| panic!("{name}: no {what} estimate in {log}"));
        line.rsplit(' ').next().unwrap().to_string()
    };
    for (col, what) in [(cols[1], "mean"), (cols[2], "weighted")] {
        let x: f64 = col.parse().unwrap();
        let want = logged(what);
        // Rounding the 6-digit value again could differ in the fifth digit
        // only when its sixth digit is a 5; compare within that unit then.
        if g5(x) != want {
            let w: f64 = want.parse().unwrap();
            let unit = 10f64.powi(w.abs().log10().floor() as i32 - 4);
            assert!(
                col.ends_with('5') && (x - w).abs() <= unit,
                "{name}: {what} {col} is not vcftools' {want}"
            );
        }
    }
}

#[test]
fn every_fst_file_matches_vcftools() {
    for FstRun {
        tag,
        file,
        pops,
        args,
        base,
    } in FST_RUNS
    {
        let (d, prefix) = fst_run(file, pops, args);
        let weir_golden = golden(&format!("{}.weir.fst", base.unwrap_or(tag)));
        let weir = String::from_utf8(read(&prefix, "weir.fst")).unwrap();
        assert_eq!(weir, weir_golden, "{tag}: .weir.fst");
        if base.is_some() {
            assert_eq!(
                String::from_utf8(read(&prefix, "windowed.weir.fst")).unwrap(),
                golden(&format!("{tag}.windowed.weir.fst")),
                "{tag}: .windowed.weir.fst"
            );
        }
        check_summary(
            &String::from_utf8(read(&prefix, "weir.fst.summary")).unwrap(),
            &golden(&format!("{tag}.weir.fst.log")),
            &weir_golden,
            tag,
        );
        let n = if base.is_some() { 3 } else { 2 };
        assert_eq!(std::fs::read_dir(d.path()).unwrap().count(), n, "{tag}");
    }
}

/// The populations select no samples: pi and Tajima's D are the same with
/// --fst-pop beside them, and Fst the same with them (decision D158,
/// docs/differences.md#popgen-independent-options).
#[test]
fn fst_populations_select_no_samples() {
    let pi = ["--site-pi", "--window-pi", "100", "--tajima-d", "100"];
    let (_a, alone) = popgen(&gpath("fst.vcf"), &pi);
    let mut args = pi.to_vec();
    args.extend(["--fst-pop", "POP1", "--fst-pop", "POP2"]);
    let (p1, p2) = (gpath("pop1.txt"), gpath("pop2.txt"));
    let args: Vec<&str> = args
        .into_iter()
        .map(|a| match a {
            "POP1" => p1.as_str(),
            "POP2" => p2.as_str(),
            a => a,
        })
        .collect();
    let (_b, both) = popgen(&gpath("fst.vcf"), &args);
    for ext in ["sites.pi", "windowed.pi", "Tajima.D"] {
        assert_eq!(read(&both, ext), read(&alone, ext), "{ext}");
    }
    assert_eq!(
        String::from_utf8(read(&both, "weir.fst")).unwrap(),
        golden("fst.weir.fst")
    );
    // pi of all 11 samples is not pi of the 7 in the populations
    let (_c, members) = popgen(
        &gpath("fst.vcf"),
        &["--site-pi", "--samples", "a1,a2,a3,a4,b1,b2,b3"],
    );
    assert_ne!(read(&members, "sites.pi"), read(&alone, "sites.pi"));
}

/// The population rules of decision D160 (docs/differences.md#fst-populations):
/// fewer than two files, a file with no selected sample and a sample in two
/// files are exit 2 naming them, with no file written; names that are not
/// selected samples are one warning with their count.
#[test]
fn fst_population_rules() {
    let d = tempdir();
    let (p1, p2, p3) = (gpath("pop1.txt"), gpath("pop2.txt"), gpath("pop3.txt"));
    let refuse = |args: &[&str], msgs: &[&str]| {
        let out = tempdir();
        let prefix = s(&out.path().join("o"));
        let mut a = gloci();
        a.args(["vcf", "popgen", &gpath("fst.vcf"), "--out-prefix", &prefix])
            .args(args);
        let mut assert = a.assert().code(2);
        for m in msgs {
            assert = assert.stderr(contains(*m));
        }
        assert_eq!(
            std::fs::read_dir(out.path()).unwrap().count(),
            0,
            "{args:?}: a refused run writes no file"
        );
    };
    refuse(
        &["--fst-pop", &p1],
        &["--fst-pop needs at least two population files"],
    );
    refuse(
        &["--fst-pop", &p1, "--fst-pop", &p2, "--samples", "a1,a2,c1"],
        &[&format!("population file {p2} names no selected sample")],
    );
    refuse(
        &[
            "--fst-pop",
            &p1,
            "--fst-pop",
            &p2,
            "--exclude-samples",
            "b1,b2,b3",
        ],
        &[&format!("population file {p2} names no selected sample")],
    );
    let both = s(&plain_file(d.path(), "both.txt", "c1\nb2\n"));
    refuse(
        &["--fst-pop", &p1, "--fst-pop", &p2, "--fst-pop", &both],
        &[&format!(
            "sample b2 is in two population files, {p2} and {both}"
        )],
    );
    refuse(
        &["--fst-pop", &p1, "--fst-pop", &p1],
        &[&format!(
            "sample a1 is in two population files, {p1} and {p1}"
        )],
    );
    refuse(
        &["--fst-pop", &p1, "--fst-pop", &p2, "--fst-window-size", "0"],
        &["--fst-window-size needs a window size above 0"],
    );
    refuse(
        &[
            "--fst-pop",
            &p1,
            "--fst-pop",
            &p2,
            "--fst-window-size",
            "100",
            "--fst-window-step",
            "101",
        ],
        &["--fst-window-step 101 is larger than the window (100)"],
    );
    refuse(
        &[
            "--fst-pop",
            &p1,
            "--fst-pop",
            &p2,
            "--fst-window-size",
            "100",
            "--fst-window-step",
            "0",
        ],
        &["--fst-window-step needs a step above 0"],
    );
    refuse(&["--fst-window-size", "100", "--site-pi"], &["--fst-pop"]);

    // pop2.txt names zz9, which is no sample: one warning, the run goes on.
    let out = tempdir();
    let prefix = s(&out.path().join("o"));
    let o = gloci()
        .args(["vcf", "popgen", &gpath("fst.vcf"), "--out-prefix", &prefix])
        .args(["--fst-pop", &p1, "--fst-pop", &p2, "--fst-pop", &p3])
        .output()
        .unwrap();
    assert!(o.status.success());
    let err = String::from_utf8(o.stderr).unwrap();
    assert_eq!(err.matches("warning").count(), 1, "{err}");
    assert!(
        err.contains("zz9 in the --fst-pop files is not a selected sample"),
        "{err}"
    );
    // Unselected samples count as unknown names, all in one warning.
    let o = gloci()
        .args(["vcf", "popgen", &gpath("fst.vcf"), "--out-prefix", &prefix])
        .args(["--fst-pop", &p1, "--fst-pop", &p2, "--fst-pop", &p3])
        .args(["--exclude-samples", "a4,c3"])
        .output()
        .unwrap();
    assert!(o.status.success());
    let err = String::from_utf8(o.stderr).unwrap();
    assert_eq!(err.matches("warning").count(), 1, "{err}");
    assert!(
        err.contains("3 names in the --fst-pop files are not selected samples")
            && err.contains("(the first is a4)"),
        "{err}"
    );
    // A population file that cannot be read is exit 1 naming it.
    let absent = s(&d.path().join("absent.txt"));
    gloci()
        .args(["vcf", "popgen", &gpath("fst.vcf"), "--out-prefix", &prefix])
        .args(["--fst-pop", &p1, "--fst-pop", &absent])
        .assert()
        .code(1)
        .stderr(contains(format!("cannot open {absent}")));
}

/// Fst windows stream their rows, so unsorted input is exit 1 with the fix
/// and no file (docs/differences.md#popgen-sorted-input); per-site Fst
/// alone takes any order.
#[test]
fn fst_windows_need_sorted_input() {
    let d = tempdir();
    let pops = [
        s(&plain_file(d.path(), "p1.txt", "s1\n")),
        s(&plain_file(d.path(), "p2.txt", "s2\n")),
    ];
    let unsorted = s(&plain_file(d.path(), "u.vcf", UNSORTED));
    let backwards = s(&plain_file(d.path(), "b.vcf", BACKWARDS));
    for (file, why) in [
        (&unsorted, "the chromosome reappears at 1:200"),
        (&backwards, "a position goes backwards at 1:499"),
    ] {
        let out = tempdir();
        let prefix = s(&out.path().join("o"));
        gloci()
            .args(["vcf", "popgen", file, "--out-prefix", &prefix])
            .args(["--fst-pop", &pops[0], "--fst-pop", &pops[1]])
            .args(["--fst-window-size", "1000"])
            .assert()
            .code(1)
            .stderr(contains("--fst-window-size needs input sorted"))
            .stderr(contains(why))
            .stderr(contains("gloci sort"));
        assert_eq!(std::fs::read_dir(out.path()).unwrap().count(), 0);
    }
    let (_o, prefix) = popgen(&unsorted, &["--fst-pop", &pops[0], "--fst-pop", &pops[1]]);
    let rows = String::from_utf8(read(&prefix, "weir.fst")).unwrap();
    assert_eq!(rows.lines().count(), 4, "{rows}");
}

/// A genotype of ploidy above two among the samples a statistic reads stops
/// the run with the record's position and writes no file, as vcftools and
/// `vcf stats` do: the selected samples for pi and Tajima's D, the population
/// members for Fst. Haploid genotypes still only drop their site.
#[test]
fn popgen_polyploid_stops_like_vcftools() {
    let d = tempdir();
    let ploidy =
        s(&Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../gloci-vcf/tests/golden/stats/ploidy.vcf"));
    // s1 is 0/0/0 at 1:10
    let p12 = s(&plain_file(d.path(), "p12.txt", "s1\ns2\n"));
    let p2 = s(&plain_file(d.path(), "p2.txt", "s2\n"));
    let p34 = s(&plain_file(d.path(), "p34.txt", "s3\ns4\n"));
    let stats: [&[&str]; 7] = [
        &["--site-pi"],
        &["--window-pi", "10"],
        &["--tajima-d", "100"],
        &["--fst-pop", &p12, "--fst-pop", &p34],
        &["--geno-r2"],
        &["--hap-r2"],
        &["--interchrom-geno-r2"],
    ];
    for stat in stats {
        let out = tempdir();
        let prefix = s(&out.path().join("o"));
        gloci()
            .args(["vcf", "popgen", &ploidy, "--out-prefix", &prefix])
            .args(stat)
            .assert()
            .code(1)
            .stderr(contains("1:10: polyploid genotypes are not supported"));
        assert_eq!(
            std::fs::read_dir(out.path()).unwrap().count(),
            0,
            "{stat:?}: a failed run writes no file"
        );
    }
    // Not among the samples read: the selection removes s1, masking removes
    // its genotype (DP 5), or s1 is in no population. vcftools stops in all
    // three (docs/differences.md#polyploid-unused-samples).
    popgen(&ploidy, &["--site-pi", "--exclude-samples", "s1"]);
    popgen(&ploidy, &["--site-pi", "--min-dp", "6"]);
    let (_o, prefix) = popgen(&ploidy, &["--fst-pop", &p2, "--fst-pop", &p34]);
    let rows = String::from_utf8(read(&prefix, "weir.fst")).unwrap();
    assert_eq!(rows.lines().count(), 4, "{rows}");
}

/// The LD runs of regen.sh: the golden's prefix, the fixture, gloci's flags
/// and the file's extension.
const LD_RUNS: [(&str, &str, &[&str], &str); 17] = [
    ("ld", "ld.vcf", &["--geno-r2"], "geno.ld"),
    (
        "ld.w2",
        "ld.vcf",
        &["--geno-r2", "--ld-window", "2"],
        "geno.ld",
    ),
    (
        "ld.bp100",
        "ld.vcf",
        &["--geno-r2", "--ld-window-bp", "100"],
        "geno.ld",
    ),
    (
        "ld.min2",
        "ld.vcf",
        &["--geno-r2", "--ld-window-min", "2"],
        "geno.ld",
    ),
    (
        "ld.bpmin10",
        "ld.vcf",
        &["--geno-r2", "--ld-window-bp-min", "10"],
        "geno.ld",
    ),
    (
        "ld.mix",
        "ld.vcf",
        &[
            "--geno-r2",
            "--ld-window",
            "3",
            "--ld-window-bp",
            "150",
            "--ld-window-min",
            "2",
            "--ld-window-bp-min",
            "5",
        ],
        "geno.ld",
    ),
    (
        "ld.r2",
        "ld.vcf",
        &["--geno-r2", "--min-r2", "0.2"],
        "geno.ld",
    ),
    (
        "ld.dp5",
        "ld.vcf",
        &["--geno-r2", "--min-dp", "5"],
        "geno.ld",
    ),
    (
        "ld.indv",
        "ld.vcf",
        &["--geno-r2", "--exclude-samples", "s1,s5"],
        "geno.ld",
    ),
    (
        "ld",
        "ld.vcf",
        &["--interchrom-geno-r2"],
        "interchrom.geno.ld",
    ),
    (
        "ld.r2",
        "ld.vcf",
        &["--interchrom-geno-r2", "--min-r2", "0.2"],
        "interchrom.geno.ld",
    ),
    ("phased", "phased.vcf", &["--hap-r2"], "hap.ld"),
    ("phased", "phased.vcf", &["--geno-r2"], "geno.ld"),
    (
        "phased.w1",
        "phased.vcf",
        &["--hap-r2", "--ld-window", "1"],
        "hap.ld",
    ),
    (
        "phased.dp5",
        "phased.vcf",
        &["--hap-r2", "--min-dp", "5"],
        "hap.ld",
    ),
    (
        "phased.r2",
        "phased.vcf",
        &["--hap-r2", "--min-r2", "0.2"],
        "hap.ld",
    ),
    (
        "phased",
        "phased.vcf",
        &["--interchrom-hap-r2"],
        "interchrom.hap.ld",
    ),
];

#[test]
fn every_ld_file_matches_vcftools() {
    for (tag, file, args, ext) in LD_RUNS {
        let (d, prefix) = popgen(&gpath(file), args);
        let name = format!("{tag}.{ext}");
        assert_eq!(
            String::from_utf8(read(&prefix, ext)).unwrap(),
            golden(&name),
            "{name} ({args:?})"
        );
        assert_eq!(std::fs::read_dir(d.path()).unwrap().count(), 1, "{name}");
    }
}

/// Every LD statistic in one run with Fst: each file is vcftools' alone, so
/// hap-r2's phase rule removes no site from the others; the unphased
/// chr1:40, chr1:60 and chr1:80 are in geno.ld and weir.fst and not in
/// hap.ld (decision D158, docs/differences.md#popgen-independent-options).
#[test]
fn ld_one_pass_keeps_its_rules_to_itself() {
    let d = tempdir();
    // p6 is haploid throughout, so pi would have no site: Fst over the
    // diploid samples stands in for the other statistics.
    let p1 = s(&plain_file(d.path(), "p1.txt", "p1\np2\n"));
    let p2 = s(&plain_file(d.path(), "p2.txt", "p3\np4\np5\n"));
    let (o, prefix) = popgen(
        &gpath("phased.vcf"),
        &[
            "--geno-r2",
            "--hap-r2",
            "--interchrom-hap-r2",
            "--fst-pop",
            &p1,
            "--fst-pop",
            &p2,
            "-t",
            "4",
        ],
    );
    let geno = String::from_utf8(read(&prefix, "geno.ld")).unwrap();
    let hap = String::from_utf8(read(&prefix, "hap.ld")).unwrap();
    assert_eq!(geno, golden("phased.geno.ld"));
    assert_eq!(hap, golden("phased.hap.ld"));
    assert_eq!(
        String::from_utf8(read(&prefix, "interchrom.hap.ld")).unwrap(),
        golden("phased.interchrom.hap.ld")
    );
    let fst = String::from_utf8(read(&prefix, "weir.fst")).unwrap();
    for pos in ["40", "60", "80"] {
        assert!(geno.contains(&format!("chr1\t10\t{pos}\t")), "{pos}");
        assert!(!hap.contains(&format!("\t{pos}\t")), "{pos}");
        assert!(fst.contains(&format!("chr1\t{pos}\t")), "{pos}: {fst}");
    }
    assert_eq!(std::fs::read_dir(o.path()).unwrap().count(), 5);
}

/// `LINES` phased sites of 8 samples on chromosome 1, 10 bases apart, and
/// one on chromosome 2.
fn ld_many(lines: u64) -> String {
    let mut t = String::from(
        "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
    );
    for i in 1..=8 {
        t.push_str(&format!("\ts{i}"));
    }
    t.push('\n');
    let gts = ["0|0", "0|1", "1|0", "1|1"];
    for i in 0..=lines {
        let chrom = if i == lines { "2" } else { "1" };
        t.push_str(&format!("{chrom}\t{}\t.\tA\tG\t.\t.\t.\tGT", 10 * (i + 1)));
        for k in 0..8u64 {
            t.push('\t');
            t.push_str(gts[((i * 7 + k * 3 + i * k) % 4) as usize]);
        }
        t.push('\n');
    }
    t
}

/// The sites LD holds stay within --max-mem: past it the run is exit 1,
/// naming the fix, and writes no file; a window holds fewer
/// (decision D161, docs/differences.md#ld-rules).
#[test]
fn ld_memory_guard() {
    let d = tempdir();
    // 8 samples: 24 bytes a site for genotypes, 32 for haplotypes, so 1K
    // holds 42 or 32 sites.
    let vcf = s(&plain_file(d.path(), "many.vcf", &ld_many(200)));
    for (stat, fixes) in [
        (
            &["--geno-r2"][..],
            &["--ld-window-bp", "--ld-window", "raise --max-mem"][..],
        ),
        (&["--hap-r2", "--min-r2", "0.5"], &["--ld-window-bp"]),
        (&["--interchrom-geno-r2"], &["raise --max-mem", "--region"]),
    ] {
        let out = tempdir();
        let prefix = s(&out.path().join("o"));
        let mut a = gloci()
            .args(["vcf", "popgen", &vcf, "--out-prefix", &prefix, "--site-pi"])
            .args(stat)
            .args(["--max-mem", "1K"])
            .assert()
            .code(1)
            .stderr(contains("pass --max-mem 1K"));
        for f in fixes {
            a = a.stderr(contains(*f));
        }
        assert_eq!(
            std::fs::read_dir(out.path()).unwrap().count(),
            0,
            "{stat:?}: a failed run writes no file"
        );
    }
    // Pairs at most 300 bases apart hold 31 sites each: 744 bytes for
    // genotypes and 992 for haplotypes, within 2K.
    let (_o, prefix) = popgen(
        &vcf,
        &[
            "--geno-r2",
            "--hap-r2",
            "--ld-window-bp",
            "300",
            "--max-mem",
            "2K",
        ],
    );
    let rows = String::from_utf8(read(&prefix, "geno.ld")).unwrap();
    assert_eq!(rows.lines().count(), 1 + 30 * 170 + 30 * 29 / 2);
    // The two windowed statistics share the budget: 1,736 bytes pass 1K.
    let out = tempdir();
    let prefix = s(&out.path().join("o"));
    gloci()
        .args(["vcf", "popgen", &vcf, "--out-prefix", &prefix])
        .args(["--geno-r2", "--hap-r2", "--ld-window-bp", "300"])
        .args(["--max-mem", "1K"])
        .assert()
        .code(1)
        .stderr(contains("holds the sites of a pair window"))
        .stderr(contains("pass --max-mem 1K"));
}

/// The windowed LD statistics stream their rows, so unsorted input is exit
/// 1 with the fix and no file (docs/differences.md#popgen-sorted-input); the
/// interchromosomal pairs take any order.
#[test]
fn ld_windows_need_sorted_input() {
    let d = tempdir();
    let unsorted = s(&plain_file(d.path(), "u.vcf", UNSORTED));
    let backwards = s(&plain_file(d.path(), "b.vcf", BACKWARDS));
    for (file, why) in [
        (&unsorted, "the chromosome reappears at 1:200"),
        (&backwards, "a position goes backwards at 1:499"),
    ] {
        let out = tempdir();
        let prefix = s(&out.path().join("o"));
        gloci()
            .args(["vcf", "popgen", file, "--out-prefix", &prefix])
            .args(["--geno-r2", "--ld-window-bp", "1000"])
            .assert()
            .code(1)
            .stderr(contains("--geno-r2 needs input sorted"))
            .stderr(contains(why))
            .stderr(contains("gloci sort"));
        assert_eq!(std::fs::read_dir(out.path()).unwrap().count(), 0);
    }
    let (_o, prefix) = popgen(&unsorted, &["--interchrom-geno-r2"]);
    assert_eq!(
        String::from_utf8(read(&prefix, "interchrom.geno.ld")).unwrap(),
        "CHR1\tPOS1\tCHR2\tPOS2\tN_INDV\tR^2\n1\t500\t2\t100\t2\t1\n\
         2\t100\t1\t200\t2\t1\n"
    );
}

/// Window options without a statistic that uses them, and windows that
/// admit no pair, are usage errors (docs/differences.md#popgen-window-step).
#[test]
fn ld_window_rules() {
    let usage = |args: &[&str], msg: &str| {
        let out = tempdir();
        let prefix = s(&out.path().join("o"));
        gloci()
            .args(["vcf", "popgen", &gpath("ld.vcf"), "--out-prefix", &prefix])
            .args(args)
            .assert()
            .code(2)
            .stderr(contains(msg));
        assert_eq!(std::fs::read_dir(out.path()).unwrap().count(), 0);
    };
    usage(
        &["--geno-r2", "--ld-window", "0"],
        "--ld-window 0 admits no pair",
    );
    usage(
        &["--hap-r2", "--ld-window", "2", "--ld-window-min", "3"],
        "--ld-window-min 3 is above --ld-window 2",
    );
    usage(
        &[
            "--geno-r2",
            "--ld-window-bp",
            "10",
            "--ld-window-bp-min",
            "11",
        ],
        "--ld-window-bp-min 11 is above --ld-window-bp 10",
    );
    usage(&["--interchrom-geno-r2", "--ld-window", "5"], "--ld-window");
    usage(&["--site-pi", "--min-r2", "0.2"], "--min-r2");
    usage(&["--geno-r2", "--max-mem", "2X"], "invalid size");
    usage(
        &[],
        "--geno-r2, --hap-r2, --interchrom-geno-r2 or --interchrom-hap-r2",
    );
}
