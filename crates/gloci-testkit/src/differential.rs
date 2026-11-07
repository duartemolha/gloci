//! Comparison rules shared by the live differential harness (`cargo xtask
//! diff`) and the committed golden tests (`crates/gloci/tests/golden.rs`).

/// Do two outputs agree under the compare mode?
///
/// - `exact`: byte for byte.
/// - `lines-unordered`: the same lines, in any order.
/// - `vcf-records`: header lines (starting with `#`) are ignored, the records
///   must be identical and in the same order.
/// - `exit-status`: the observed values are exit codes (see [`observe`]); they
///   must be equal. This pins error behaviour (gloci refuses, the reference
///   tool accepts, or the reverse).
pub fn compare(mode: &str, a: &[u8], b: &[u8]) -> bool {
    match mode {
        "exact" => a == b,
        "lines-unordered" => sorted_lines(a) == sorted_lines(b),
        "vcf-records" => records(a) == records(b),
        "exit-status" => a == b,
        other => panic!("unknown compare mode `{other}`"),
    }
}

/// What a command produced, as compared and stored: the exit code for the
/// `exit-status` mode, else standard output with the machine-specific
/// `##glociCommand` header line normalised.
pub fn observe(mode: &str, code: Option<i32>, stdout: &[u8]) -> Vec<u8> {
    if mode == "exit-status" {
        return format!("{}\n", code.unwrap_or(-1)).into_bytes();
    }
    let mut out = Vec::with_capacity(stdout.len());
    for line in stdout.split_inclusive(|&c| c == b'\n') {
        if line.starts_with(b"##glociCommand=") {
            out.extend_from_slice(b"##glociCommand=<normalised>\n");
        } else {
            out.extend_from_slice(line);
        }
    }
    out
}

fn sorted_lines(x: &[u8]) -> Vec<&[u8]> {
    let mut v: Vec<&[u8]> = x.split(|&c| c == b'\n').collect();
    v.sort_unstable();
    v
}

fn records(x: &[u8]) -> Vec<&[u8]> {
    x.split(|&c| c == b'\n')
        .filter(|l| !l.starts_with(b"#"))
        .collect()
}

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum Verdict {
    /// Equal where equality was expected, or different where a documented
    /// difference was expected.
    Pass,
    /// Different with no documented difference.
    Mismatch,
    /// A documented difference whose outputs are now equal.
    StaleDifference,
}

pub fn verdict(expected_difference: Option<&str>, equal: bool) -> Verdict {
    match (expected_difference.is_some(), equal) {
        (false, true) | (true, false) => Verdict::Pass,
        (false, false) => Verdict::Mismatch,
        (true, true) => Verdict::StaleDifference,
    }
}

/// The anchors a document defines: `{#name}` at the end of a heading line.
pub fn anchors(doc: &str) -> Vec<&str> {
    doc.lines()
        .filter(|l| l.starts_with('#'))
        .filter_map(|l| {
            let l = l.trim_end().strip_suffix('}')?;
            let (_, a) = l.rsplit_once("{#")?;
            Some(a)
        })
        .collect()
}

/// Check that `anchor` (with or without a leading `#`) is defined in `doc`.
pub fn check_anchor(doc: &str, anchor: &str) -> Result<(), String> {
    let a = anchor.strip_prefix('#').unwrap_or(anchor);
    if anchors(doc).contains(&a) {
        Ok(())
    } else {
        Err(format!(
            "expected_difference `{anchor}` has no `{{#{a}}}` heading in docs/differences.md"
        ))
    }
}
