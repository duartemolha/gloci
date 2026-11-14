//! Tests for `gloci vcf popgen` (decisions D158, D159): every file matches
//! the vcftools 0.1.16 file committed in `gloci-vcf/tests/golden/popgen` (made
//! by its `regen.sh`), one pass equals separate runs, the thread count and the
//! input format change nothing, and the windowed statistics need sorted input.
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
        "add at least one of --site-pi, --window-pi SIZE or --tajima-d BIN",
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
