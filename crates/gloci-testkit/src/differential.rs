//! Comparison rules shared by the live differential harness (`cargo xtask
//! diff`) and the committed golden tests (`crates/gloci/tests/golden.rs`),
//! and the record-level match counts of the real-data concordance run (D97).

use std::collections::HashMap;
use std::io::{self, BufRead};

/// Do two outputs agree under the compare mode?
///
/// - `exact`: byte for byte.
/// - `lines-unordered`: the same lines, in any order.
/// - `vcf-records`: meta-header lines (starting with `##`) are ignored; the
///   `#CHROM` line and data lines must be identical and in the same order.
/// - `exit-status`: the observed values are exit codes (see [`observe`]); they
///   must be equal. This pins error behaviour (gloci refuses, the reference
///   tool accepts, or the reverse).
/// - `numeric`: line counts match; tab fields are compared pairwise. Two
///   numbers are equal within a relative tolerance of 1e-5 (absolute when
///   both are below 1e-9); any other field must match exactly.
/// - `sam`: SAM text; the `@PG` header lines are ignored.
/// - `vcf`: VCF text, header included; only the lines each tool stamps with its command
///   (`##glociCommand=`, `##bcftools_*Version=`, `##bcftools_*Command=`) are ignored.
/// - `bam` (D118): both are BAM files. They are decompressed, and then the
///   header text without its `@PG` lines, the reference list and the record
///   bytes must be identical; the BGZF block layout and level do not matter.
/// - `cram` (D141): both are CRAM files with embedded or no references. They are decoded through
///   `open_aln_source`; the header text without its `@PG` lines, the references and the
///   records (as BAM bytes) must be equal. The container layout and codecs do not matter.
pub fn compare(mode: &str, a: &[u8], b: &[u8]) -> bool {
    match mode {
        "exact" => a == b,
        "lines-unordered" => sorted_lines(a) == sorted_lines(b),
        "vcf-records" => vcf_records_equal(a, b),
        "exit-status" => a == b,
        "numeric" => numeric_equal(a, b),
        "sam" => {
            let keep = |x: &[u8]| -> Vec<Vec<u8>> {
                x.split(|&c| c == b'\n')
                    .filter(|l| !l.starts_with(b"@PG\t"))
                    .map(<[u8]>::to_vec)
                    .collect()
            };
            keep(a) == keep(b)
        }
        "vcf" => {
            let stamp = |l: &[u8]| {
                l.starts_with(b"##glociCommand=")
                    || (l.starts_with(b"##bcftools_")
                        && l.splitn(2, |&c| c == b'=')
                            .next()
                            .is_some_and(|k| k.ends_with(b"Version") || k.ends_with(b"Command")))
            };
            let keep = |x: &[u8]| -> Vec<Vec<u8>> {
                x.split(|&c| c == b'\n')
                    .filter(|l| !stamp(l))
                    .map(<[u8]>::to_vec)
                    .collect()
            };
            keep(a) == keep(b)
        }
        "bam" => match (bam_parts(a), bam_parts(b)) {
            (Some(x), Some(y)) => x == y,
            _ => false,
        },
        "cram" => match (cram_parts(a), cram_parts(b)) {
            (Some(x), Some(y)) => x == y,
            _ => false,
        },
        other => panic!("unknown compare mode `{other}`"),
    }
}

/// Split a golden `args` line into arguments: whitespace separates, and a
/// single-quoted run is one argument's text verbatim (so an expression such as
/// `'QUAL>30 && ID~"^rs"'` survives). `None` for an unbalanced quote.
pub fn split_args(line: &str) -> Option<Vec<String>> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let (mut in_word, mut quoted) = (false, false);
    for c in line.chars() {
        match (c, quoted) {
            ('\'', _) => {
                quoted = !quoted;
                in_word = true;
            }
            (c, false) if c.is_whitespace() => {
                if in_word {
                    out.push(std::mem::take(&mut cur));
                    in_word = false;
                }
            }
            (c, _) => {
                cur.push(c);
                in_word = true;
            }
        }
    }
    if quoted {
        return None;
    }
    if in_word {
        out.push(cur);
    }
    Some(out)
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
        if matches!(mode, "bam" | "sam") && line.starts_with(b"@PG\tID:gloci") {
            // gloci's own @PG line names the machine's command line (D118).
            continue;
        }
        if line.starts_with(b"##glociCommand=") {
            out.extend_from_slice(b"##glociCommand=<normalised>\n");
        } else {
            out.extend_from_slice(line);
        }
    }
    out
}

/// A BAM file with the `@PG` lines and the compression taken out:
/// (header text, references as name and length, record bytes).
type BamParts = (Vec<u8>, Vec<(Vec<u8>, u32)>, Vec<u8>);

fn bam_parts(bgzf: &[u8]) -> Option<BamParts> {
    use std::io::Read;
    let mut raw = Vec::new();
    noodles_bgzf::io::Reader::new(bgzf)
        .read_to_end(&mut raw)
        .ok()?;
    let mut at = 0usize;
    let mut take = |n: usize| -> Option<&[u8]> {
        let s = raw.get(at..at.checked_add(n)?)?;
        at += n;
        Some(s)
    };
    let u32_of = |s: &[u8]| u32::from_le_bytes([s[0], s[1], s[2], s[3]]);
    if take(4)? != b"BAM\x01" {
        return None;
    }
    let l_text = u32_of(take(4)?) as usize;
    let text = take(l_text)?;
    let text = text.split(|&c| c == 0).next().unwrap_or(&[]);
    let mut header = Vec::new();
    for line in text.split_inclusive(|&c| c == b'\n') {
        if !line.starts_with(b"@PG\t") {
            header.extend_from_slice(line);
        }
    }
    let n_ref = u32_of(take(4)?);
    let mut refs = Vec::new();
    for _ in 0..n_ref {
        let l_name = u32_of(take(4)?) as usize;
        let name = take(l_name)?.to_vec();
        let l_ref = u32_of(take(4)?);
        refs.push((name, l_ref));
    }
    // The records must be framed: each is `block_size` and that many bytes (at least the
    // 32 fixed bytes), with nothing left over, so a truncated file does not compare equal.
    let records = raw.get(at..)?.to_vec();
    let mut i = 0usize;
    while i < records.len() {
        let size = u32::from_le_bytes(records.get(i..i + 4)?.try_into().ok()?) as usize;
        if size < 32 || i.checked_add(4 + size)? > records.len() {
            return None;
        }
        i += 4 + size;
    }
    Some((header, refs, records))
}

/// A CRAM file with the `@PG` lines and the encoding taken out: (header text, references as
/// name and length, each record's BAM bytes with its length). `None` when it does not decode
/// (a reference that is needed but not embedded included).
type CramParts = (Vec<u8>, Vec<(Vec<u8>, u64)>, Vec<u8>);

fn cram_parts(cram: &[u8]) -> Option<CramParts> {
    use gloci_io::{AlignOptions, open_aln_source_from};
    let opts = gloci_io::ReaderOptions {
        threads: std::num::NonZeroUsize::MIN,
        format: Some(gloci_io::Format::Cram),
        vcf_span: gloci_core::VcfSpan::Allele,
        invalid: gloci_io::InvalidPolicy::Fail,
        chrom_map: None,
        accept: gloci_io::Accept {
            alignments: Some(AlignOptions {
                full_records: true,
                ..AlignOptions::default()
            }),
            ..gloci_io::Accept::default()
        },
    };
    let mut src =
        open_aln_source_from("cram", Box::new(std::io::Cursor::new(cram.to_vec())), &opts).ok()?;
    let header = src.header().clone();
    let mut text = Vec::new();
    for line in header.text.split_inclusive(|&c| c == b'\n') {
        if !line.starts_with(b"@PG\t") {
            text.extend_from_slice(line);
        }
    }
    let mut records = Vec::new();
    let mut buf = Vec::new();
    loop {
        buf.clear();
        match src.next_record(&mut buf) {
            Ok(Some(_)) => {
                records.extend_from_slice(&(buf.len() as u64).to_le_bytes());
                records.extend_from_slice(&buf);
            }
            Ok(None) => break,
            Err(_) => return None,
        }
    }
    Some((text, header.refs.clone(), records))
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
        .filter(|l| !l.is_empty() && !l.starts_with(b"##"))
        .collect()
}

fn vcf_records_equal(a: &[u8], b: &[u8]) -> bool {
    let mut ra = records(a);
    let mut rb = records(b);
    let a_has_header = ra.first().is_some_and(|l| l.starts_with(b"#"));
    let b_has_header = rb.first().is_some_and(|l| l.starts_with(b"#"));
    if a_has_header != b_has_header {
        if a_has_header {
            ra.remove(0);
        } else {
            rb.remove(0);
        }
    }
    ra == rb
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

/// The compare modes of the real-data concordance run (D97), as [`line_key`]
/// reads them.
pub const CONCORDANCE_MODES: &[&str] = &["exact", "numeric", "vcf-records"];

/// How far gloci's output agrees with a reference tool's, record by record
/// (D97): the lines of each, and how many pair up as multisets (a line of
/// either side is used at most once, so order does not matter).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MatchStats {
    pub n_gloci: u64,
    pub n_reference: u64,
    pub n_matching: u64,
    /// Both sides have the same lines in the same order.
    pub same_order: bool,
}

impl MatchStats {
    /// `n_matching` over the larger side; 1.0 when both are empty.
    pub fn rate(&self) -> f64 {
        match self.n_gloci.max(self.n_reference) {
            0 => 1.0,
            n => self.n_matching as f64 / n as f64,
        }
    }
}

/// A number of a `numeric` line that is compared with a tolerance: its
/// value and the decimal places its text gives (digits after the point, less
/// the exponent: `0.25` has 2, `1.5e-3` has 4).
type Approx = (f64, i32);

/// One output line under a concordance compare mode: the key that must be
/// equal, and for `numeric` the numbers compared with a tolerance.
struct Parsed {
    key: u64,
    approx: Vec<Approx>,
}

/// Relative difference two doubles may have from computing in another order.
const REL_TOL: f64 = 1e-9;

fn decimals(text: &[u8]) -> i32 {
    let s = std::str::from_utf8(text).unwrap_or("");
    let (mant, exp) = match s.find(['e', 'E']) {
        Some(i) => (&s[..i], s[i + 1..].parse::<i32>().unwrap_or(0)),
        None => (s, 0),
    };
    let frac = mant.find('.').map_or(0, |i| (mant.len() - i - 1) as i32);
    frac - exp
}

/// Two numbers agree when they differ by no more than half a unit of the last
/// decimal place of the less precise text (a tool printing fewer digits
/// rounded there), or by a relative [`REL_TOL`].
fn approx_equal(x: Approx, y: Approx) -> bool {
    let diff = (x.0 - y.0).abs();
    let mag = x.0.abs().max(y.0.abs());
    let printed = 0.5 * 10f64.powi(-x.1.min(y.1));
    diff <= printed + 4.0 * f64::EPSILON * mag || diff <= REL_TOL * mag
}

fn parse_line(mode: &str, line: &[u8]) -> Option<Parsed> {
    use std::hash::Hasher;
    let mut h = std::collections::hash_map::DefaultHasher::new();
    let mut approx = Vec::new();
    match mode {
        "exact" => h.write(line),
        "vcf-records" => {
            if line.starts_with(b"#") {
                return None;
            }
            h.write(line)
        }
        "numeric" => {
            for (i, f) in line.split(|&c| c == b'\t').enumerate() {
                if i > 0 {
                    h.write_u8(b'\t');
                }
                match numeric_field(f) {
                    Numeric::Int(v) => h.write(v.to_string().as_bytes()),
                    Numeric::Approx(x) => {
                        h.write_u8(0x01);
                        approx.push((x, decimals(f)));
                    }
                    Numeric::Text => h.write(f),
                }
            }
        }
        other => panic!("unknown concordance compare mode `{other}`"),
    }
    Some(Parsed {
        key: h.finish(),
        approx,
    })
}

enum Numeric {
    /// An integer literal, or a number within [`REL_TOL`] of an integer.
    Int(i128),
    Approx(f64),
    Text,
}

fn numeric_field(f: &[u8]) -> Numeric {
    if is_int_literal(f) {
        return match std::str::from_utf8(f)
            .ok()
            .and_then(|s| s.parse::<i128>().ok())
        {
            Some(v) => Numeric::Int(v),
            None => Numeric::Text,
        };
    }
    let Some(x) = parse_f64(f) else {
        return Numeric::Text;
    };
    let r = x.round();
    if (x - r).abs() <= REL_TOL * x.abs().max(1.0) && r.abs() < 1e15 {
        Numeric::Int(r as i128)
    } else {
        Numeric::Approx(x)
    }
}

/// The key of one output line (without its newline) under a concordance
/// compare mode, or `None` for a line the mode skips. Lines can match only
/// when their keys are equal:
///
/// - `exact`: the bytes of the line.
/// - `vcf-records`: the same, and every line that starts with `#` (meta lines
///   and the column-name line) is skipped, as tools write different headers.
/// - `numeric`: tab fields; an integer literal counts by its value (`007` is
///   7), and so does a number within a relative 1e-9 of an integer (`3.0`,
///   `2.9999999999999996`); any other number is compared apart, with the
///   tolerance of the less precise text: half a unit of its last decimal
///   place, so `0.1929825` (a tool printing 7 decimals) equals
///   `0.19298245614035087`, while `0.012345` does not equal `0.012346`. Other
///   fields count by their bytes.
///
/// Panics on another mode ([`match_stats`] checks it first).
pub fn line_key(mode: &str, line: &[u8]) -> Option<u64> {
    parse_line(mode, line).map(|p| p.key)
}

/// Call `f` with the key and the text of every line of `r` that `mode` keeps,
/// until `f` returns false. Returns the number of lines kept and a digest of
/// their keys in order.
pub fn for_each_line_key(
    mode: &str,
    r: &mut dyn BufRead,
    mut f: impl FnMut(u64, &[u8]) -> bool,
) -> io::Result<(u64, u64)> {
    for_each_parsed(mode, r, |p, line| f(p.key, line))
}

fn for_each_parsed(
    mode: &str,
    r: &mut dyn BufRead,
    mut f: impl FnMut(Parsed, &[u8]) -> bool,
) -> io::Result<(u64, u64)> {
    use std::hash::Hasher;
    let mut order = std::collections::hash_map::DefaultHasher::new();
    let (mut n, mut buf) = (0u64, Vec::new());
    loop {
        buf.clear();
        if r.read_until(b'\n', &mut buf)? == 0 {
            break;
        }
        let line = buf.strip_suffix(b"\n").unwrap_or(&buf);
        if let Some(p) = parse_line(mode, line) {
            n += 1;
            order.write_u64(p.key);
            if !f(p, line) {
                break;
            }
        }
    }
    Ok((n, order.finish()))
}

/// [`MatchStats`] of two outputs under a concordance compare mode
/// ([`CONCORDANCE_MODES`]), streamed: what is kept is a multiset of line keys
/// (with the tolerance-compared numbers of `numeric` lines).
pub fn match_stats(
    mode: &str,
    gloci: &mut dyn BufRead,
    reference: &mut dyn BufRead,
) -> io::Result<MatchStats> {
    match_stats_with_residual(mode, gloci, reference).map(|(s, _)| s)
}

/// The lines [`match_stats_with_residual`] left unmatched on each side, which a
/// second pass over the outputs can pick out with [`Residual::take`].
pub struct Residual {
    mode: String,
    /// Per key: lines only gloci wrote (positive) or only the reference
    /// wrote (negative); `exact` and `vcf-records`.
    counts: HashMap<u64, i64>,
    /// Per key: the numbers of unmatched `numeric` lines, gloci's and the
    /// reference's.
    gloci: HashMap<u64, Vec<Box<[Approx]>>>,
    reference: HashMap<u64, Vec<Box<[Approx]>>>,
}

impl Residual {
    /// Whether `line` (without its newline) is one of the unmatched lines of
    /// gloci's output (`gloci_side`) or of the reference's. Each unmatched
    /// line answers true once, so repeated lines are counted out.
    pub fn take(&mut self, gloci_side: bool, line: &[u8]) -> bool {
        let Some(p) = parse_line(&self.mode, line) else {
            return false;
        };
        if self.mode != "numeric" {
            let sign = if gloci_side { 1 } else { -1 };
            return match self.counts.get_mut(&p.key) {
                Some(c) if *c * sign > 0 => {
                    *c -= sign;
                    true
                }
                _ => false,
            };
        }
        let side = if gloci_side {
            &mut self.gloci
        } else {
            &mut self.reference
        };
        let Some(list) = side.get_mut(&p.key) else {
            return false;
        };
        match list.iter().position(|v| **v == *p.approx) {
            Some(i) => {
                list.swap_remove(i);
                true
            }
            None => false,
        }
    }
}

/// [`match_stats`], and the unmatched lines ([`Residual`]).
pub fn match_stats_with_residual(
    mode: &str,
    gloci: &mut dyn BufRead,
    reference: &mut dyn BufRead,
) -> io::Result<(MatchStats, Residual)> {
    use std::collections::hash_map::Entry;
    if !CONCORDANCE_MODES.contains(&mode) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "unknown concordance compare mode `{mode}` (use exact, numeric or vcf-records)"
            ),
        ));
    }
    let mut res = Residual {
        mode: mode.to_string(),
        counts: HashMap::new(),
        gloci: HashMap::new(),
        reference: HashMap::new(),
    };
    let numeric = mode == "numeric";
    let (n_gloci, order_g) = for_each_parsed(mode, gloci, |p, _| {
        if numeric {
            res.gloci
                .entry(p.key)
                .or_default()
                .push(p.approx.into_boxed_slice());
        } else {
            *res.counts.entry(p.key).or_insert(0) += 1;
        }
        true
    })?;
    let mut n_matching = 0u64;
    let (n_reference, order_r) = for_each_parsed(mode, reference, |p, _| {
        if numeric {
            let found = res.gloci.get_mut(&p.key).and_then(|list| {
                let i = list.iter().position(|g| {
                    g.len() == p.approx.len()
                        && g.iter().zip(&p.approx).all(|(&a, &b)| approx_equal(a, b))
                })?;
                list.swap_remove(i);
                Some(list.is_empty())
            });
            match found {
                Some(emptied) => {
                    n_matching += 1;
                    if emptied {
                        res.gloci.remove(&p.key);
                    }
                }
                None => res
                    .reference
                    .entry(p.key)
                    .or_default()
                    .push(p.approx.into_boxed_slice()),
            }
        } else {
            match res.counts.entry(p.key) {
                Entry::Occupied(mut e) => {
                    if *e.get() > 0 {
                        n_matching += 1;
                    }
                    *e.get_mut() -= 1;
                    if *e.get() == 0 {
                        e.remove();
                    }
                }
                Entry::Vacant(e) => {
                    e.insert(-1);
                }
            }
        }
        true
    })?;
    let stats = MatchStats {
        n_gloci,
        n_reference,
        n_matching,
        same_order: n_gloci == n_reference && n_matching == n_gloci && order_g == order_r,
    };
    Ok((stats, res))
}

#[cfg(test)]
mod tests {
    use super::{compare, observe, split_args};

    /// A one-record BAM (`@HD`, one `@SQ`, record `r1` 4M ACGT on chr1:1): `pg` is the `@PG`
    /// lines, `level` the BGZF level and `qual` the first quality byte.
    fn bam(pg: &str, level: u8, qual: u8) -> Vec<u8> {
        use noodles_bgzf::io::writer::{Builder, CompressionLevel};
        use std::io::Write;
        let text = format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n{pg}");
        let mut raw = b"BAM\x01".to_vec();
        raw.extend_from_slice(&(text.len() as u32).to_le_bytes());
        raw.extend_from_slice(text.as_bytes());
        raw.extend_from_slice(&1u32.to_le_bytes());
        raw.extend_from_slice(&5u32.to_le_bytes());
        raw.extend_from_slice(b"chr1\0");
        raw.extend_from_slice(&1000u32.to_le_bytes());
        // block: refID, pos, l_read_name, mapq, bin, n_cigar, flag, l_seq, next_refID,
        // next_pos, tlen, read_name, cigar, seq, qual.
        let mut rec = Vec::new();
        rec.extend_from_slice(&0i32.to_le_bytes());
        rec.extend_from_slice(&0i32.to_le_bytes());
        rec.extend_from_slice(&[3, 60]);
        rec.extend_from_slice(&4681u16.to_le_bytes());
        rec.extend_from_slice(&1u16.to_le_bytes());
        rec.extend_from_slice(&0u16.to_le_bytes());
        rec.extend_from_slice(&4u32.to_le_bytes());
        rec.extend_from_slice(&(-1i32).to_le_bytes());
        rec.extend_from_slice(&(-1i32).to_le_bytes());
        rec.extend_from_slice(&0i32.to_le_bytes());
        rec.extend_from_slice(b"r1\0");
        rec.extend_from_slice(&(4u32 << 4).to_le_bytes());
        rec.extend_from_slice(&[0x12, 0x48]);
        rec.extend_from_slice(&[qual, 30, 30, 30]);
        raw.extend_from_slice(&(rec.len() as u32).to_le_bytes());
        raw.extend_from_slice(&rec);
        let mut w = Builder::default()
            .set_compression_level(CompressionLevel::try_from(level).unwrap())
            .build_from_writer(Vec::new());
        w.write_all(&raw).unwrap();
        w.finish().unwrap()
    }

    /// `bam` cut short by `n` decompressed bytes, recompressed.
    fn truncated(bgzf: &[u8], n: usize) -> Vec<u8> {
        use std::io::{Read, Write};
        let mut raw = Vec::new();
        noodles_bgzf::io::Reader::new(bgzf)
            .read_to_end(&mut raw)
            .unwrap();
        raw.truncate(raw.len() - n);
        let mut w = noodles_bgzf::io::Writer::new(Vec::new());
        w.write_all(&raw).unwrap();
        w.finish().unwrap()
    }

    #[test]
    fn compare_bam_ignores_pg_and_compression() {
        let a = bam("@PG\tID:gloci\tPN:gloci\n", 1, 7);
        let b = bam("@PG\tID:samtools\tPN:samtools\tVN:1.21\n", 9, 7);
        assert_ne!(a, b);
        assert!(compare("bam", &a, &b));
        assert!(compare("bam", &a, &a));
        assert!(!compare("bam", &a, &bam("@PG\tID:gloci\n", 1, 8)));
        assert!(!compare("bam", &a, b"not a bam"));
        // A truncated record fails framing, even compared with itself.
        let cut = truncated(&a, 3);
        assert!(!compare("bam", &cut, &cut));
    }

    const FLAGS_CRAM: &[u8] = include_bytes!("../../../xtask/diff/inputs/aln/flags.embed.cram");
    const CIGARS_CRAM: &[u8] = include_bytes!("../../../xtask/diff/inputs/aln/cigars.embed.cram");

    /// `cram` with its header container replaced by one holding `text`.
    fn with_header(cram: &[u8], text: &[u8]) -> Vec<u8> {
        let path = std::env::temp_dir().join(format!(
            "gloci-cram-compare-{}-{}.cram",
            std::process::id(),
            cram.len()
        ));
        std::fs::write(&path, cram).unwrap();
        let p = gloci_io::read_cram_prelude(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        let mut out = cram[..26].to_vec();
        gloci_io::write_cram_header_container(text, &mut out).unwrap();
        out.extend_from_slice(&cram[p.data_start as usize..]);
        out
    }

    fn header_text(cram: &[u8]) -> Vec<u8> {
        let path = std::env::temp_dir().join(format!(
            "gloci-cram-header-{}-{}.cram",
            std::process::id(),
            cram.len()
        ));
        std::fs::write(&path, cram).unwrap();
        let p = gloci_io::read_cram_prelude(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        let mut t = p.header_text;
        if t.last() != Some(&b'\n') {
            t.push(b'\n');
        }
        t
    }

    #[test]
    fn compare_cram_ignores_pg_and_container_layout() {
        let text = header_text(FLAGS_CRAM);
        assert!(compare("cram", FLAGS_CRAM, FLAGS_CRAM));
        // The header text is the same but for @PG lines: equal, though the bytes differ.
        let mut pg = text.clone();
        pg.extend_from_slice(b"@PG\tID:other\tPN:other\n");
        let b = with_header(FLAGS_CRAM, &pg);
        assert_ne!(b, FLAGS_CRAM);
        assert!(compare("cram", FLAGS_CRAM, &b));
        // Any other header line counts.
        let mut co = text.clone();
        co.extend_from_slice(b"@CO\textra\n");
        assert!(!compare("cram", FLAGS_CRAM, &with_header(FLAGS_CRAM, &co)));
        // Same header, other records.
        let other = with_header(CIGARS_CRAM, &text);
        assert!(!compare("cram", FLAGS_CRAM, &other));
        // Not a CRAM file.
        assert!(!compare("cram", FLAGS_CRAM, b"not a cram"));
    }

    #[test]
    fn at_pg_filter_is_scoped_to_bam_and_sam_modes() {
        let text = b"@HD\tVN:1.6\n@PG\tID:gloci\tPN:gloci\nr1\t0\n";
        assert_eq!(observe("sam", Some(0), text), b"@HD\tVN:1.6\nr1\t0\n");
        assert_eq!(observe("exact", Some(0), text), text);
    }

    #[test]
    fn split_args_keeps_single_quoted_runs() {
        assert_eq!(
            split_args("vcf filter  input.vcf --include 'A>1 && B~\"x\"'").unwrap(),
            ["vcf", "filter", "input.vcf", "--include", "A>1 && B~\"x\""]
        );
        assert_eq!(split_args("a ''").unwrap(), ["a", ""]);
        assert!(split_args("a 'b").is_none());
    }

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

    fn stats(mode: &str, g: &str, r: &str) -> super::MatchStats {
        super::match_stats(mode, &mut g.as_bytes(), &mut r.as_bytes()).unwrap()
    }

    #[test]
    fn match_stats_modes() {
        let ten: String = (0..10).map(|i| format!("chr1\t{i}\t{}\n", i + 5)).collect();
        // Identical inputs.
        let s = stats("exact", &ten, &ten);
        assert_eq!((s.n_gloci, s.n_reference, s.n_matching), (10, 10, 10));
        assert_eq!(s.rate(), 1.0);
        assert!(s.same_order);
        // One changed line in 10.
        let changed = ten.replace("chr1\t3\t8\n", "chr1\t3\t9\n");
        let s = stats("exact", &changed, &ten);
        assert_eq!((s.n_gloci, s.n_reference, s.n_matching), (10, 10, 9));
        assert!((s.rate() - 0.9).abs() < 1e-12);
        // Reordered lines: all match, the order does not.
        let mut lines: Vec<&str> = ten.split_inclusive('\n').collect();
        lines.reverse();
        let s = stats("exact", &lines.concat(), &ten);
        assert_eq!(s.n_matching, 10);
        assert_eq!(s.rate(), 1.0);
        assert!(!s.same_order);
        // A line repeated on one side matches once.
        let s = stats("exact", "a\na\n", "a\nb\n");
        assert_eq!((s.n_gloci, s.n_reference, s.n_matching), (2, 2, 1));
        // A missing final newline does not change the line.
        assert_eq!(stats("exact", "a\nb", "a\nb\n").rate(), 1.0);
        // vcf-records ignores every `#` line, the column-name line included.
        let s = stats(
            "vcf-records",
            "##fileformat=VCFv4.2\n##glociCommand=x\n#CHROM\tPOS\nchr1\t5\n",
            "##fileformat=VCFv4.2\nchr1\t5\n",
        );
        assert_eq!((s.n_gloci, s.n_reference, s.n_matching), (1, 1, 1));
        assert!(s.same_order);
        assert!(stats("exact", "#h\nx\n", "x\n").rate() < 1.0);
        // numeric: 3 and 3.0 are equal, and a number printed with fewer digits
        // equals the one it rounds.
        assert_eq!(stats("numeric", "x\t3\n", "x\t3.0\n").rate(), 1.0);
        assert_eq!(stats("numeric", "x\t19\n", "x\t19.000000\n").rate(), 1.0);
        assert_eq!(
            stats("numeric", "x\t0.0123456789\n", "x\t0.0123457\n").rate(),
            1.0
        );
        // 7 decimals, as bedtools coverage prints, of a value at a tie of the 6th digit.
        assert_eq!(
            stats("numeric", "x\t0.19298245614035087\n", "x\t0.1929825\n").rate(),
            1.0
        );
        assert_eq!(
            stats("numeric", "x\t2.9044049999999997\n", "x\t2.904405\n").rate(),
            1.0
        );
        assert_eq!(
            stats("numeric", "x\t2.9999999999999996\n", "x\t3\n").rate(),
            1.0
        );
        assert_eq!(stats("numeric", "x\t1.5e-3\n", "x\t0.0015\n").rate(), 1.0);
        assert_eq!(stats("numeric", "x\t-0.0\n", "x\t0\n").rate(), 1.0);
        assert!(stats("numeric", "x\t0.012345\n", "x\t0.012346\n").rate() < 1.0);
        assert!(stats("numeric", "x\t0.0028571\n", "x\t0.0028572\n").rate() < 1.0);
        assert!(stats("numeric", "x\t2.9999999\n", "x\t3\n").rate() < 1.0);
        assert!(stats("numeric", "x\t3.4\n", "x\t3\n").rate() < 1.0);
        // Numbers are matched within their line, any order of lines.
        let s = stats("numeric", "a\t0.25\nb\t0.5\n", "b\t0.50000\na\t0.250\n");
        assert_eq!((s.n_matching, s.same_order), (2, false));
        let s = stats("numeric", "a\t0.25\na\t0.5\n", "a\t0.5\na\t0.6\n");
        assert_eq!(s.n_matching, 1);
        // Integer literals must match exactly, other text exactly.
        assert!(stats("numeric", "x\t100000000\n", "x\t100000500\n").rate() < 1.0);
        assert!(stats("numeric", "x\t.\n", "x\tNA\n").rate() < 1.0);
        assert!(stats("numeric", "x\t3\n", "y\t3\n").rate() < 1.0);
        assert!(stats("exact", "x\t3\n", "x\t3.0\n").rate() < 1.0);
        // Two empty inputs agree.
        let s = stats("exact", "", "");
        assert_eq!((s.n_gloci, s.n_reference, s.n_matching), (0, 0, 0));
        assert_eq!(s.rate(), 1.0);
        assert!(s.same_order);
        // One empty side matches nothing.
        assert_eq!(stats("exact", "", "a\n").rate(), 0.0);
        // Unknown modes are an error, not a guess.
        assert!(super::match_stats("bam", &mut &b""[..], &mut &b""[..]).is_err());
    }

    #[test]
    fn residual_names_the_unmatched_lines() {
        let (s, mut residual) =
            super::match_stats_with_residual("exact", &mut &b"a\nb\nb\n"[..], &mut &b"b\nc\n"[..])
                .unwrap();
        assert_eq!(s.n_matching, 1);
        assert!(residual.take(true, b"a"));
        assert!(!residual.take(true, b"a"));
        assert!(residual.take(true, b"b"));
        assert!(!residual.take(true, b"b"));
        assert!(!residual.take(false, b"b"));
        assert!(residual.take(false, b"c"));
        let (s, mut residual) = super::match_stats_with_residual(
            "numeric",
            &mut &b"x\t0.25\nx\t0.75\n"[..],
            &mut &b"x\t0.250\nx\t0.5\n"[..],
        )
        .unwrap();
        assert_eq!(s.n_matching, 1);
        assert!(!residual.take(true, b"x\t0.25"));
        assert!(residual.take(true, b"x\t0.75"));
        assert!(residual.take(false, b"x\t0.5"));
        assert_eq!(super::line_key("vcf-records", b"#CHROM\tPOS"), None);
        assert_eq!(
            super::line_key("numeric", b"x\t3"),
            super::line_key("numeric", b"x\t3.000")
        );
    }
}
