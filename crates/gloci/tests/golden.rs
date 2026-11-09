//! Golden reference outputs: for every case directory under `tests/golden`,
//! run gloci with the stored `args` on the stored input and compare with
//! `expected`, the output of the pinned reference tool (bedtools, bedops, GNU
//! sort, htslib, bcftools), written by `cargo xtask diff --update-golden`.
//!
//! The rules match the live harness: a case without `expected_difference`
//! must match, one with it must differ (a stale entry fails), and the anchor
//! must exist in docs/differences.md.
use gloci_testkit::differential::{Verdict, check_anchor, compare, observe, split_args, verdict};
use std::path::{Path, PathBuf};

struct GoldenCase {
    name: String,
    dir: PathBuf,
    args: Vec<String>,
    compare: String,
    expected_difference: Option<String>,
    expected: Vec<u8>,
    /// gloci's own output where it differs from the reference tool, pinned
    /// exactly (present exactly when `expected_difference` is).
    expected_gloci: Option<Vec<u8>>,
}

fn read_line(dir: &Path, file: &str) -> Option<String> {
    std::fs::read_to_string(dir.join(file))
        .ok()
        .map(|s| s.trim().to_string())
}

fn golden_cases() -> Vec<GoldenCase> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    let mut v: Vec<_> = std::fs::read_dir(&root)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.is_dir())
        .map(|dir| GoldenCase {
            name: dir.file_name().unwrap().to_string_lossy().into_owned(),
            args: split_args(
                &read_line(&dir, "args")
                    .unwrap_or_else(|| panic!("{}: no args file", dir.display())),
            )
            .unwrap_or_else(|| panic!("{}: unbalanced quote in args", dir.display())),
            compare: read_line(&dir, "compare").unwrap_or_else(|| "exact".into()),
            expected_difference: read_line(&dir, "expected_difference"),
            expected: std::fs::read(dir.join("expected")).unwrap(),
            expected_gloci: std::fs::read(dir.join("expected_gloci")).ok(),
            dir,
        })
        .collect();
    v.sort_by(|a, b| a.name.cmp(&b.name));
    v
}

impl GoldenCase {
    /// Run gloci in the case directory (inputs are referenced by name).
    fn run_gloci(&self) -> Vec<u8> {
        let out = assert_cmd::Command::cargo_bin("gloci")
            .unwrap()
            .current_dir(&self.dir)
            .args(&self.args)
            .output()
            .unwrap();
        assert!(
            self.compare == "exit-status" || out.status.success(),
            "{}: gloci {:?} failed: {}",
            self.name,
            self.args,
            String::from_utf8_lossy(&out.stderr)
        );
        observe(&self.compare, out.status.code(), &out.stdout)
    }

    fn verdict(&self, got: &[u8]) -> Verdict {
        verdict(
            self.expected_difference.as_deref(),
            compare(&self.compare, got, &self.expected),
        )
    }
}

#[test]
fn every_golden_case_agrees_with_its_reference_tool() {
    let doc = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/differences.md"),
    )
    .unwrap();
    let cases = golden_cases();
    for want in [
        "sort-bed-vs-bedtools-sort",
        "sort-bed-vs-sort-bed",
        "sort-bed-vs-gnu-sort",
        "sort-vcf-vs-bcftools-sort",
        "view-region-ref-inclusive-vs-tabix",
        "view-region-default-vs-tabix",
    ] {
        assert!(
            cases.iter().any(|c| c.name == want),
            "missing golden case {want}"
        );
    }
    for c in cases {
        if let Some(a) = &c.expected_difference {
            check_anchor(&doc, a).unwrap_or_else(|e| panic!("{}: {e}", c.name));
        }
        let got = c.run_gloci();
        assert_eq!(
            c.expected_gloci.is_some(),
            c.expected_difference.is_some(),
            "{}: expected_gloci must exist exactly for documented differences",
            c.name
        );
        if let Some(pinned) = &c.expected_gloci {
            assert_eq!(
                String::from_utf8_lossy(&got),
                String::from_utf8_lossy(pinned),
                "{}: gloci's output changed (it still differs from the reference tool, but differently)",
                c.name
            );
        }
        assert_eq!(
            c.verdict(&got),
            Verdict::Pass,
            "{}: {}",
            c.name,
            if c.expected_difference.is_some() {
                "stale expected_difference (gloci now equals the reference tool)"
            } else {
                "gloci output differs from the reference tool"
            }
        );
    }
}
