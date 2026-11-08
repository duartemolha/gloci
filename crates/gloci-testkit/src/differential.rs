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
/// - `numeric`: line counts match; tab fields are compared pairwise. Two
///   numbers are equal within a relative tolerance of 1e-5 (absolute when
///   both are below 1e-9); any other field must match exactly.
pub fn compare(mode: &str, a: &[u8], b: &[u8]) -> bool {
    match mode {
        "exact" => a == b,
        "lines-unordered" => sorted_lines(a) == sorted_lines(b),
        "vcf-records" => records(a) == records(b),
        "exit-status" => a == b,
        "numeric" => numeric_equal(a, b),
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

fn numeric_equal(a: &[u8], b: &[u8]) -> bool {
    let (la, lb): (Vec<&[u8]>, Vec<&[u8]>) = (
        a.split(|&c| c == b'\n').collect(),
        b.split(|&c| c == b'\n').collect(),
    );
    la.len() == lb.len()
        && la.iter().zip(&lb).all(|(x, y)| {
            let (fx, fy): (Vec<&[u8]>, Vec<&[u8]>) = (
                x.split(|&c| c == b'\t').collect(),
                y.split(|&c| c == b'\t').collect(),
            );
            fx.len() == fy.len() && fx.iter().zip(&fy).all(|(p, q)| field_equal(p, q))
        })
}

fn parse_f64(x: &[u8]) -> Option<f64> {
    let f: f64 = std::str::from_utf8(x).ok()?.parse().ok()?;
    f.is_finite().then_some(f)
}

fn is_int_literal(x: &[u8]) -> bool {
    let d = match x.first() {
        Some(b'-' | b'+') => &x[1..],
        _ => x,
    };
    !d.is_empty() && d.iter().all(u8::is_ascii_digit)
}

fn field_equal(p: &[u8], q: &[u8]) -> bool {
    if p == q {
        return true;
    }
    // Integer literals (coordinates, counts) must match exactly.
    if is_int_literal(p) && is_int_literal(q) {
        return false;
    }
    match (parse_f64(p), parse_f64(q)) {
        (Some(x), Some(y)) => {
            let diff = (x - y).abs();
            if x.abs() < 1e-9 && y.abs() < 1e-9 {
                diff < 1e-9
            } else {
                diff <= 1e-5 * x.abs().max(y.abs())
            }
        }
        _ => false,
    }
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

#[cfg(test)]
mod tests {
    use super::compare;

    #[test]
    fn numeric_mode() {
        assert!(compare(
            "numeric",
            b"a\t1.0000001\t2\n",
            b"a\t1.0000002\t2\n"
        ));
        assert!(compare("numeric", b"x\t3\n", b"x\t3.0\n"));
        assert!(compare("numeric", b"1e-12\n", b"3e-12\n"));
        assert!(!compare("numeric", b"1\n", b"1.001\n"));
        assert!(!compare("numeric", b"a\t1\n", b"b\t1\n"));
        assert!(!compare("numeric", b"a\t1\n", b"a\t1\t2\n"));
        assert!(!compare("numeric", b"1\n2\n", b"1\n"));
        assert!(!compare("numeric", b"100000000\n", b"100000500\n"));
        assert!(!compare("numeric", b"100000\n", b"100001\n"));
        assert!(compare("numeric", b"100000000\n", b"100000000\n"));
        assert!(!compare("numeric", b".\n", b"NA\n"));
        assert!(!compare("numeric", b".\n", b"0\n"));
        assert!(!compare("numeric", b"nan\n", b"nan2\n"));
        assert!(!compare("numeric", b"inf\n", b"1e308\n"));
        assert!(!compare("numeric", b"inf\n", b"-inf\n"));
        assert!(compare("numeric", b"nan\n", b"nan\n"));
        assert!(!compare("exact", b"x\t3\n", b"x\t3.0\n"));
    }
}
