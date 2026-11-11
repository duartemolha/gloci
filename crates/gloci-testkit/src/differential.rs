//! Comparison rules shared by the live differential harness (`cargo xtask
//! diff`) and the committed golden tests (`crates/gloci/tests/golden.rs`).

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
        gloci_io::write_cram_header_container(p.major, text, &mut out).unwrap();
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
}
