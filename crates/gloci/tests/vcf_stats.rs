//! Tests for `gloci vcf stats` (decision D72): the statistics match the
//! vcftools 0.1.16 files committed in `gloci-vcf/tests/golden/stats` (made by
//! its `regen.sh`), one pass equals separate runs, and a failed run publishes
//! nothing.
mod common;
use common::*;
use predicates::prelude::*;
use predicates::str::contains;
use std::path::{Path, PathBuf};

fn golden_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../gloci-vcf/tests/golden/stats")
}

fn golden(prefix: &str, ext: &str) -> Option<String> {
    std::fs::read_to_string(golden_dir().join(format!("{prefix}.{ext}"))).ok()
}

fn fixture(name: &str) -> String {
    golden_dir()
        .join(format!("{name}.vcf"))
        .to_string_lossy()
        .into_owned()
}

/// gloci flag and vcftools extension of every statistic.
const STATS: [(&str, &str); 15] = [
    ("--freq", "frq"),
    ("--counts", "frq.count"),
    ("--depth", "idepth"),
    ("--site-depth", "ldepth"),
    ("--site-mean-depth", "ldepth.mean"),
    ("--geno-depth", "gdepth"),
    ("--site-quality", "lqual"),
    ("--missing-site", "lmiss"),
    ("--missing-indv", "imiss"),
    ("--het", "het"),
    ("--hardy", "hwe"),
    ("--singletons", "singletons"),
    ("--filter-summary", "FILTER.summary"),
    ("--tstv-summary", "TsTv.summary"),
    ("--indel-hist", "indel.hist"),
];

const FIXTURES: [&str; 16] = [
    "main",
    "mixed",
    "quality",
    "ploidy",
    "filters",
    "filtersorted",
    "density",
    "hwegrid",
    "hwebig",
    "masked",
    "ftmask",
    "persample",
    "persmall",
    "permask",
    "nan",
    "mq",
];

fn stats_args<'a>(file: &'a str, prefix: &'a str, rest: &[&'a str]) -> Vec<&'a str> {
    let mut v = vec!["vcf", "stats", file, "--out-prefix", prefix];
    v.extend_from_slice(rest);
    v
}

/// Run one statistic and compare its file with vcftools' golden. Where
/// vcftools stopped without writing the file, gloci must fail too and write
/// nothing.
fn check_one(file: &str, golden_prefix: &str, flag: &str, ext: &str, mask: &[&str]) {
    let d = tempdir();
    let prefix = d.path().join("o");
    let prefix = prefix.to_string_lossy().into_owned();
    let mut rest = vec![flag];
    rest.extend_from_slice(mask);
    let args = stats_args(file, &prefix, &rest);
    match golden(golden_prefix, ext) {
        Some(want) => {
            gloci().args(&args).assert().success();
            let got = std::fs::read_to_string(format!("{prefix}.{ext}")).unwrap();
            assert_eq!(got, want, "{golden_prefix}.{ext} ({flag} {mask:?})");
        }
        None => {
            assert!(
                golden_prefix.starts_with("ploidy"),
                "no vcftools golden for {golden_prefix}.{ext}"
            );
            gloci()
                .args(&args)
                .assert()
                .code(1)
                .stderr(contains("polyploid"));
            assert_eq!(std::fs::read_dir(d.path()).unwrap().count(), 0);
        }
    }
}

#[test]
fn every_statistic_matches_vcftools_on_every_fixture() {
    for f in FIXTURES {
        let file = fixture(f);
        for (flag, ext) in STATS {
            // The per-site files of the extra fixtures are exercised too.
            check_one(&file, f, flag, ext, &[]);
        }
    }
}

#[test]
fn snp_density_matches_vcftools() {
    for f in ["main", "density"] {
        for bin in ["1000", "7", "2500"] {
            let d = tempdir();
            let prefix = d.path().join("o").to_string_lossy().into_owned();
            gloci()
                .args(stats_args(&fixture(f), &prefix, &["--snp-density", bin]))
                .assert()
                .success();
            assert_eq!(
                std::fs::read_to_string(format!("{prefix}.snpden")).unwrap(),
                golden(&format!("{f}.b{bin}"), "snpden").unwrap(),
                "{f} bin {bin}"
            );
        }
    }
}

/// A fixture, its golden prefix, the mask flags and the statistics (or carrier
/// substitutions) to check.
type MaskRun<'a> = (&'a str, &'a str, &'a [&'a str], &'a [(&'a str, &'a str)]);

/// The masking and sample-selection runs vcftools was given (regen.sh).
/// Per-sample statistics exist for each; the per-site ones only for the
/// `site` prefix.
struct Masked {
    fixture: &'static str,
    tag: &'static str,
    args: &'static [&'static str],
}

const FT: &str = "--mask-filtered-genotypes";

fn masked_runs() -> Vec<Masked> {
    let masks: [(&str, &[&str]); 6] = [
        ("gq10", &["--min-gq", "10"]),
        ("gq20", &["--min-gq", "20"]),
        ("dp5", &["--min-dp", "5"]),
        ("maxdp12", &["--max-dp", "12"]),
        ("ft", &[FT]),
        (
            "combo",
            &["--min-dp", "4", "--max-dp", "20", "--min-gq", "10", FT],
        ),
    ];
    let indv: (&str, &[&str]) = ("indv", &["--samples", "s2,s1"]);
    let mut v = Vec::new();
    // Fixtures whose DP and GQ values are all present: every mask is
    // comparable. A `.` value is masked by vcftools and not by gloci (D50),
    // so persample, persmall, nan and mq only run the sample selection (nan
    // also a threshold that no `.` value reaches), and
    // `missing_dp_gq_is_not_masked_unlike_vcftools` pins the difference.
    for fixture in ["masked", "ftmask", "permask"] {
        for (tag, args) in masks.iter().chain([&indv]) {
            v.push(Masked { fixture, tag, args });
        }
    }
    for fixture in ["persample", "persmall", "nan"] {
        v.push(Masked {
            fixture,
            tag: indv.0,
            args: indv.1,
        });
    }
    v
}

#[test]
fn per_sample_statistics_match_vcftools_under_masking_and_selection() {
    let runs = masked_runs();
    assert_eq!(runs.len(), 3 * 7 + 3);
    for m in runs {
        let file = fixture(m.fixture);
        for (flag, ext) in [
            ("--depth", "idepth"),
            ("--missing-indv", "imiss"),
            ("--het", "het"),
        ] {
            let prefix = format!("{}.{}", m.fixture, m.tag);
            check_one(&file, &prefix, flag, ext, m.args);
        }
    }
}

#[test]
fn site_statistics_match_vcftools_under_masking() {
    let gt = [
        ("--freq", "frq"),
        ("--counts", "frq.count"),
        ("--site-depth", "ldepth"),
        ("--site-mean-depth", "ldepth.mean"),
        ("--geno-depth", "gdepth"),
        ("--missing-site", "lmiss"),
        ("--hardy", "hwe"),
    ];
    let runs: [MaskRun; 5] = [
        ("masked", "gq20", &["--min-gq", "20"], &gt),
        (
            "masked",
            "dp5",
            &["--min-dp", "5"],
            &[
                ("--site-mean-depth", "ldepth.mean"),
                ("--geno-depth", "gdepth"),
                ("--missing-site", "lmiss"),
            ],
        ),
        ("ftmask", "ft", &[FT], &gt),
        ("ftmask", "maxdp6", &["--max-dp", "6"], &gt),
        (
            "ftmask",
            "combo",
            &["--min-dp", "4", "--max-dp", "7", "--min-gq", "10", FT],
            &gt,
        ),
    ];
    for (f, prefix, mask, stats) in runs {
        for (flag, ext) in stats {
            check_one(&fixture(f), prefix, flag, ext, mask);
        }
    }
}

/// Intended difference, anchored in docs/differences.md#vcf-stats-singletons:
/// vcftools names a genotype that masking removed as the carrier of a
/// singleton, and decides a doubleton from a raw homozygote. gloci looks only
/// at unmasked genotypes.
#[test]
fn singletons_name_an_unmasked_carrier_unlike_vcftools() {
    let cases: [MaskRun; 4] = [
        (
            "masked",
            "gq20",
            &["--min-gq", "20"],
            &[("1\t30\tS\tG\ts1\n", "1\t30\tS\tG\ts2\n")],
        ),
        (
            "ftmask",
            "ft",
            &[FT],
            &[
                ("1\t50\tD\tG\ts1\n", "1\t50\tD\tG\ts2\n"),
                ("1\t60\tD\tA\ts4\n", ""),
            ],
        ),
        (
            "ftmask",
            "maxdp6",
            &["--max-dp", "6"],
            &[("1\t60\tS\tG\ts1\n", "1\t60\tS\tG\ts2\n")],
        ),
        (
            "ftmask",
            "combo",
            &["--min-dp", "4", "--max-dp", "7", "--min-gq", "10", FT],
            &[],
        ),
    ];
    for (f, prefix, mask, fixes) in cases {
        let d = tempdir();
        let out = d.path().join("o").to_string_lossy().into_owned();
        let mut rest = vec!["--singletons"];
        rest.extend_from_slice(mask);
        gloci()
            .args(stats_args(&fixture(f), &out, &rest))
            .assert()
            .success();
        let got = std::fs::read_to_string(format!("{out}.singletons")).unwrap();
        let mut want = golden(prefix, "singletons").unwrap();
        for (from, to) in fixes {
            assert!(want.contains(from), "{prefix}: golden lacks {from:?}");
            want = want.replace(from, to);
        }
        assert_eq!(got, want, "{prefix}.singletons");
    }
}

/// Anchored in docs/differences.md#vcf-stats-mask-missing: a genotype whose
/// DP or GQ is `.` is masked by vcftools (when the FORMAT lists the field)
/// and kept by gloci (D50: a missing value never masks).
#[test]
fn missing_dp_gq_is_not_masked_unlike_vcftools() {
    let file = fixture("mq");
    for (tag, mask, vcftools_filtered, gloci_filtered) in [
        ("gq10", ["--min-gq", "10"], [1, 1, 1, 0], [0, 0, 1, 0]),
        ("dp4", ["--min-dp", "4"], [2, 0, 2, 0], [1, 0, 0, 0]),
    ] {
        let d = tempdir();
        let p = s(&d.path().join("o"));
        let mut rest = vec!["--missing-indv"];
        rest.extend(mask);
        gloci()
            .args(stats_args(&file, &p, &rest))
            .assert()
            .success();
        let got = std::fs::read_to_string(format!("{p}.imiss")).unwrap();
        let col = |text: &str| -> Vec<u64> {
            text.lines()
                .skip(1)
                .map(|l| l.split('\t').nth(2).unwrap().parse().unwrap())
                .collect()
        };
        let vcftools = golden(&format!("mq.{tag}"), "imiss").unwrap();
        assert_eq!(col(&vcftools), vcftools_filtered, "{tag} vcftools");
        assert_eq!(col(&got), gloci_filtered, "{tag} gloci");
    }
}

/// A sample cell cut short before DP or GQ (`0/1:9` under GT:DP:GQ, or a bare
/// `0/0`) is masked by vcftools and kept by gloci, like a `.` (D50).
#[test]
fn truncated_sample_cells_are_not_masked_unlike_vcftools() {
    let file = fixture("mqt");
    for (tag, mask, vcftools_filtered) in [
        ("gq10", ["--min-gq", "10"], [1, 1, 1, 1]),
        ("dp4", ["--min-dp", "4"], [0, 0, 1, 1]),
    ] {
        let d = tempdir();
        let p = s(&d.path().join("o"));
        let mut rest = vec!["--missing-indv"];
        rest.extend(mask);
        gloci()
            .args(stats_args(&file, &p, &rest))
            .assert()
            .success();
        let got = std::fs::read_to_string(format!("{p}.imiss")).unwrap();
        let col = |text: &str| -> Vec<u64> {
            text.lines()
                .skip(1)
                .map(|l| l.split('\t').nth(2).unwrap().parse().unwrap())
                .collect()
        };
        let vcftools = golden(&format!("mqt.{tag}"), "imiss").unwrap();
        assert_eq!(col(&vcftools), vcftools_filtered, "{tag} vcftools");
        assert_eq!(col(&got), [0, 0, 0, 0], "{tag} gloci");
    }
}

/// The type rules of the summary, as the `--summary` help and the docs state
/// them: mixed SNP and indel ALTs are other, other takes sv, an all-SNP
/// multi-allelic record is snp, and `<*>` is ignored beside a real ALT.
#[test]
fn summary_record_type_rules() {
    let d = tempdir();
    let rows = [
        "1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1",     // snp
        "1\t20\t.\tA\tG,C\t.\t.\t.\tGT\t1/2",   // multi-allelic, all SNPs: snp
        "1\t30\t.\tA\tG,AT\t.\t.\t.\tGT\t1/2",  // SNP and indel: other
        "1\t40\t.\tA\t<DEL>\t.\t.\t.\tGT\t0/1", // sv: other
        "1\t50\t.\tA\tG,<*>\t.\t.\t.\tGT\t0/1", // <*> ignored: snp
        "1\t60\t.\tA\t<*>\t.\t.\t.\tGT\t0/0",   // only <*>: ref
        "1\t70\t.\tA\tAT\t.\t.\t.\tGT\t0/1",    // indel
        "1\t80\t.\tAC\tGT\t.\t.\t.\tGT\t0/1",   // mnp
    ];
    let text = format!(
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n{}\n",
        rows.join("\n")
    );
    let file = s(&plain_file(d.path(), "types.vcf", &text));
    let p = s(&d.path().join("o"));
    gloci()
        .args(stats_args(&file, &p, &["--json"]))
        .assert()
        .success();
    let json = std::fs::read_to_string(format!("{p}.summary.json")).unwrap();
    let want = "\"counts_by_type\": {\n    \"indel\": 1,\n    \"mnp\": 1,\n    \"other\": 2,\n    \"ref\": 1,\n    \"snp\": 3\n  }";
    assert!(json.contains(want), "{json}");
}

fn read(p: &str) -> Vec<u8> {
    std::fs::read(p).unwrap()
}

#[test]
fn one_pass_equals_separate_runs() {
    let file = fixture("persample");
    let mask = [
        "--samples",
        "s1,s2,s3,s5,s6",
        "--min-dp",
        "4",
        "--min-gq",
        "10",
    ];
    let five = [
        "--freq",
        "--het",
        "--missing-indv",
        "--site-mean-depth",
        "--hardy",
    ];
    let d = tempdir();
    let one = d.path().join("one").to_string_lossy().into_owned();
    let mut rest: Vec<&str> = five.to_vec();
    rest.extend(mask);
    gloci()
        .args(stats_args(&file, &one, &rest))
        .assert()
        .success();
    for (flag, ext) in STATS.iter().filter(|(f, _)| five.contains(f)) {
        let sep = d
            .path()
            .join(format!("sep{ext}"))
            .to_string_lossy()
            .into_owned();
        let mut rest = vec![*flag];
        rest.extend(mask);
        gloci()
            .args(stats_args(&file, &sep, &rest))
            .assert()
            .success();
        assert_eq!(
            read(&format!("{one}.{ext}")),
            read(&format!("{sep}.{ext}")),
            "{flag}"
        );
    }
    // Only the five requested files exist.
    let n = std::fs::read_dir(d.path())
        .unwrap()
        .filter(|e| {
            e.as_ref()
                .unwrap()
                .file_name()
                .to_string_lossy()
                .starts_with("one.")
        })
        .count();
    assert_eq!(n, 5);
}

#[test]
fn bcf_input_gives_the_same_files() {
    let d = tempdir();
    let all: Vec<&str> = STATS
        .iter()
        .map(|(f, _)| *f)
        .chain(["--summary", "--json"])
        .collect();
    // (`ploidy` stops the genotype statistics; the persample fixtures have
    // an all-`.` FT column, which the BCF converter rejects.)
    for f in ["main", "mixed", "masked", "ftmask"] {
        let bcf = d.path().join(format!("{f}.bcf"));
        gloci()
            .args([
                "vcf",
                "filter",
                &fixture(f),
                "--vcf-span",
                "ref",
                "-o",
                &bcf.to_string_lossy(),
            ])
            .assert()
            .success();
        let (a, b) = (
            s(&d.path().join(format!("{f}.a"))),
            s(&d.path().join(format!("{f}.b"))),
        );
        let mut rest = all.clone();
        rest.extend(["--min-gq", "10"]);
        gloci()
            .args(stats_args(&fixture(f), &a, &rest))
            .assert()
            .success();
        gloci()
            .args(stats_args(&bcf.to_string_lossy(), &b, &rest))
            .assert()
            .success();
        for (_, ext) in STATS
            .iter()
            .chain(&[("", "summary.txt"), ("", "summary.json")])
        {
            assert_eq!(
                read(&format!("{a}.{ext}")),
                read(&format!("{b}.{ext}")),
                "{f}.{ext}"
            );
        }
    }
}

#[test]
fn output_is_identical_at_every_thread_count() {
    let d = tempdir();
    let all: Vec<&str> = STATS
        .iter()
        .map(|(f, _)| *f)
        .chain(["--summary", "--json"])
        .collect();
    let mut runs = Vec::new();
    for t in ["1", "2", "7"] {
        let p = d
            .path()
            .join(format!("t{t}"))
            .to_string_lossy()
            .into_owned();
        let mut rest = all.clone();
        rest.extend(["--threads", t]);
        gloci()
            .args(stats_args(&fixture("persample"), &p, &rest))
            .assert()
            .success();
        runs.push(p);
    }
    for (_, ext) in STATS
        .iter()
        .chain(&[("", "summary.txt"), ("", "summary.json")])
    {
        for p in &runs[1..] {
            assert_eq!(
                read(&format!("{}.{ext}", runs[0])),
                read(&format!("{p}.{ext}")),
                "{ext}"
            );
        }
    }
}

#[test]
fn selection_options_apply_before_the_statistics() {
    // `--include` and a region run through the shared selection stage; the
    // result equals filtering first and computing the statistics after.
    let d = tempdir();
    let file = fixture("persample");
    let filtered = d.path().join("f.vcf").to_string_lossy().into_owned();
    gloci()
        .args([
            "vcf",
            "filter",
            &file,
            "--include",
            "QUAL>40 && N_PASS(GT=\"het\")>1",
            "--chrom",
            "1",
            "--min-dp",
            "5",
            "-o",
            &filtered,
        ])
        .assert()
        .success();
    let (a, b) = (
        d.path().join("a").to_string_lossy().into_owned(),
        d.path().join("b").to_string_lossy().into_owned(),
    );
    let stats = ["--freq", "--het", "--hardy", "--site-quality"];
    let mut rest = stats.to_vec();
    rest.extend([
        "--include",
        "QUAL>40 && N_PASS(GT=\"het\")>1",
        "--chrom",
        "1",
        "--min-dp",
        "5",
    ]);
    gloci()
        .args(stats_args(&file, &a, &rest))
        .assert()
        .success();
    // The filter wrote masked genotypes as missing, so the masking is
    // already applied; the stats of the filtered file need no mask flag.
    gloci()
        .args(stats_args(&filtered, &b, &stats))
        .assert()
        .success();
    for ext in ["frq", "het", "hwe", "lqual"] {
        assert_eq!(
            read(&format!("{a}.{ext}")),
            read(&format!("{b}.{ext}")),
            "{ext}"
        );
    }
}

fn summary_fixture(d: &Path) -> String {
    let text = "\
##fileformat=VCFv4.2
##contig=<ID=1,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb
1\t10\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:10\t0/0:20
1\t20\t.\tC\tA\t50\tPASS\t.\tGT:DP\t0/1:30\t./.:.
1\t30\t.\tAT\tA\t50\tPASS\t.\tGT:DP\t0/0:10\t0/0:10
1\t40\t.\tAC\tCA\t50\tPASS\t.\tGT:DP\t0/1:10\t0/0:10
1\t50\t.\tA\tG,T\t50\tPASS\t.\tGT:DP\t1/2:10\t0/0:10
1\t60\t.\tA\t<DEL>\t50\tPASS\tEND=65\tGT:DP\t0/1:10\t0/0:10
1\t70\t.\tA\t.\t50\tPASS\t.\tGT:DP\t0/0:10\t0/0:10
";
    s(&plain_file(d, "sum.vcf", text))
}

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

#[test]
fn summary_text_and_json_follow_the_schema() {
    let d = tempdir();
    let file = summary_fixture(d.path());
    let p = s(&d.path().join("o"));
    gloci()
        .args(stats_args(&file, &p, &["--summary", "--json"]))
        .assert()
        .success();
    let txt = std::fs::read_to_string(format!("{p}.summary.txt")).unwrap();
    assert_eq!(
        txt,
        "Records:        7\n\
Samples:        2\n\
Multiallelic:   1 records\n\
Counts by type:\n\
\x20 snp   3\n\
\x20 mnp   1\n\
\x20 indel 1\n\
\x20 other 1\n\
\x20 ref   1\n\
Ts/Tv:          1 (1 transitions, 1 transversions, biallelic SNPs)\n\
Mean depth:     12.3077 (over 13 genotypes with a depth)\n\
Missingness:    0.0714286 (1 of 14 genotypes missing, 0 of them masked)\n"
    );
    let json = std::fs::read_to_string(format!("{p}.summary.json")).unwrap();
    assert_eq!(
        json,
        "{\n  \"counts_by_type\": {\n    \"indel\": 1,\n    \"mnp\": 1,\n    \"other\": 1,\n    \
\"ref\": 1,\n    \"snp\": 3\n  },\n  \"depth\": {\n    \"mean\": 12.3077,\n    \"values\": 13\n  },\n  \
\"genotypes\": {\n    \"filtered\": 0,\n    \"missing\": 1,\n    \"missing_fraction\": 0.0714286,\n    \
\"total\": 14\n  },\n  \"multiallelic_sites\": 1,\n  \"records\": 7,\n  \"samples\": 2,\n  \
\"schema_version\": 1,\n  \"ts_tv\": {\n    \"ratio\": 1,\n    \"transitions\": 1,\n    \
\"transversions\": 1\n  }\n}\n"
    );
}

#[test]
fn summary_alone_and_json_alone_write_one_file_each() {
    let d = tempdir();
    let file = summary_fixture(d.path());
    for (flag, ext, other) in [
        ("--summary", "summary.txt", "summary.json"),
        ("--json", "summary.json", "summary.txt"),
    ] {
        let p = s(&d.path().join(ext));
        gloci()
            .args(stats_args(&file, &p, &[flag]))
            .assert()
            .success();
        assert!(Path::new(&format!("{p}.{ext}")).exists());
        assert!(!Path::new(&format!("{p}.{other}")).exists());
    }
}

#[test]
fn summary_counts_masked_genotypes_and_nulls_undefined_ratios() {
    let d = tempdir();
    let file = summary_fixture(d.path());
    let p = s(&d.path().join("o"));
    gloci()
        .args(stats_args(
            &file,
            &p,
            &["--json", "--min-dp", "15", "--types", "indel,mnv"],
        ))
        .assert()
        .success();
    let json = std::fs::read_to_string(format!("{p}.summary.json")).unwrap();
    // 2 records, 4 genotypes: all but the DP-20/30 ones are masked; no SNPs
    // pass the type filter, so Ts/Tv has no denominator.
    assert!(json.contains("\"filtered\": 4"), "{json}");
    assert!(json.contains("\"ratio\": null"), "{json}");
    assert!(json.contains("\"mean\": null"), "{json}");
}

#[test]
fn usage_errors_exit_2() {
    let d = tempdir();
    let file = fixture("main");
    let p = s(&d.path().join("o"));
    // No statistic.
    gloci()
        .args(["vcf", "stats", &file, "--out-prefix", &p])
        .assert()
        .code(2)
        .stderr(
            contains("no statistic requested; add at least one of --freq")
                .and(contains("--summary or --json"))
                .and(contains("required arguments").not()),
        );
    // No prefix.
    gloci()
        .args(["vcf", "stats", &file, "--freq"])
        .assert()
        .code(2)
        .stderr(contains("--out-prefix"));
    // A zero bin.
    gloci()
        .args(stats_args(&file, &p, &["--snp-density", "0"]))
        .assert()
        .code(2)
        .stderr(contains("bin size"));
    // A bad threshold from the shared selection.
    gloci()
        .args(stats_args(&file, &p, &["--freq", "--min-maf", "2"]))
        .assert()
        .code(2);
    assert_eq!(std::fs::read_dir(d.path()).unwrap().count(), 0);
}

#[test]
fn unsorted_input_stops_snp_density_with_a_hint_and_publishes_nothing() {
    let d = tempdir();
    let text = std::fs::read_to_string(fixture("density")).unwrap();
    // Reverse the data lines: chromosome order and positions go backwards.
    let (head, rows): (Vec<&str>, Vec<&str>) = text.lines().partition(|l| l.starts_with('#'));
    let mut rows = rows;
    rows.reverse();
    let unsorted = format!("{}\n{}\n", head.join("\n"), rows.join("\n"));
    let file = s(&plain_file(d.path(), "unsorted.vcf", &unsorted));
    let p = s(&d.path().join("o"));
    // Earlier records already wrote rows to --freq's temp file.
    gloci()
        .args(stats_args(&file, &p, &["--freq", "--snp-density", "1000"]))
        .assert()
        .code(1)
        .stderr(contains("run: gloci sort FILE"));
    let left: Vec<_> = std::fs::read_dir(d.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    assert_eq!(left, vec!["unsorted.vcf".to_string()], "{left:?}");
}

#[test]
fn an_error_in_the_middle_publishes_nothing_and_keeps_old_files() {
    let d = tempdir();
    let text = std::fs::read_to_string(fixture("persmall")).unwrap();
    let mut lines: Vec<String> = text.lines().map(String::from).collect();
    let half = lines.len() / 2;
    lines[half] = "1\tnot-a-position\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/0\t0/0".into();
    let file = s(&plain_file(d.path(), "bad.vcf", &(lines.join("\n") + "\n")));
    let p = s(&d.path().join("o"));
    std::fs::write(format!("{p}.frq"), "old").unwrap();
    gloci()
        .args(stats_args(
            &file,
            &p,
            &["--freq", "--het", "--summary", "--site-quality"],
        ))
        .assert()
        .failure();
    let mut left: Vec<_> = std::fs::read_dir(d.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    left.sort();
    assert_eq!(left, vec!["bad.vcf", "o.frq"], "{left:?}");
    assert_eq!(std::fs::read_to_string(format!("{p}.frq")).unwrap(), "old");
}

#[test]
fn a_failed_publish_keeps_the_previous_outputs_and_names_the_squatter() {
    let d = tempdir();
    let p = s(&d.path().join("o"));
    // A previous run's file, and a directory on the destination of --hardy:
    // the publish stops before moving anything, and says what is there.
    std::fs::write(format!("{p}.frq"), "previous run").unwrap();
    std::fs::create_dir(format!("{p}.hwe")).unwrap();
    let a = gloci()
        .args(stats_args(&fixture("main"), &p, &["--freq", "--hardy"]))
        .assert()
        .code(1)
        .stderr(contains("o.hwe").and(contains("a directory with that name exists")));
    let err = String::from_utf8_lossy(&a.get_output().stderr).into_owned();
    assert!(!err.contains("free disk space"), "{err}");
    let mut left: Vec<_> = std::fs::read_dir(d.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    left.sort();
    assert_eq!(left, vec!["o.frq", "o.hwe"], "{left:?}");
    assert_eq!(
        std::fs::read(format!("{p}.frq")).unwrap(),
        b"previous run".to_vec()
    );
}

#[test]
fn a_successful_run_replaces_the_previous_outputs_and_leaves_no_backup() {
    let d = tempdir();
    let p = s(&d.path().join("o"));
    std::fs::write(format!("{p}.frq"), "previous run").unwrap();
    gloci()
        .args(stats_args(&fixture("main"), &p, &["--freq", "--hardy"]))
        .assert()
        .success();
    let mut left: Vec<_> = std::fs::read_dir(d.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    left.sort();
    assert_eq!(left, vec!["o.frq", "o.hwe"], "{left:?}");
    assert_ne!(std::fs::read(format!("{p}.frq")).unwrap(), b"previous run");
}

#[test]
fn a_missing_output_directory_is_reported_with_the_path() {
    let d = tempdir();
    let p = s(&d.path().join("nope/o"));
    gloci()
        .args(stats_args(&fixture("main"), &p, &["--freq"]))
        .assert()
        .code(1)
        .stderr(contains("nope/o.frq"));
}
