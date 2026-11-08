//! VCF record and header rewriting for `vcf filter`.

use super::samples::SampleSelection;
use crate::stats::{MaskRules, SiteStats};
use gloci_io::VcfHeader;

/// A plan describing how to rewrite VCF records and headers.
#[derive(Clone, Debug, PartialEq)]
pub struct RewritePlan {
    pub selection: Option<SampleSelection>,
    pub mask: MaskRules,
    pub keep_info: Option<Vec<Vec<u8>>>,
    pub drop_info: Vec<Vec<u8>>,
    pub recalc_ac: bool,
}

impl RewritePlan {
    /// True when the plan makes no changes to records or headers:
    /// all samples kept, no masking, no INFO pruning, no recalc.
    pub fn is_passthrough(&self) -> bool {
        self.selection
            .as_ref()
            .is_none_or(|s| s.keep.iter().all(|&k| k))
            && self.mask == MaskRules::default()
            && self.keep_info.is_none()
            && self.drop_info.is_empty()
            && !self.recalc_ac
    }
}

/// Format a float like C `printf("%g", x)`: 6 significant digits,
/// trailing zeros trimmed, exponential notation with at least two
/// exponent digits when exponent < -4 or >= 6.
pub fn fmt_g6(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x.is_sign_negative() { "-inf" } else { "inf" }.to_string();
    }
    if x == 0.0 {
        return if x.is_sign_negative() { "-0" } else { "0" }.to_string();
    }

    let s = format!("{:.5e}", x);
    let (mantissa_str, exp_str) = s.split_once('e').expect("e in scientific format");
    let exp: i32 = exp_str.parse().expect("valid exp");

    let is_neg = mantissa_str.starts_with('-');
    let mantissa_clean = if is_neg {
        &mantissa_str[1..]
    } else {
        mantissa_str
    };

    let mut digits = String::with_capacity(6);
    for ch in mantissa_clean.chars() {
        if ch != '.' {
            digits.push(ch);
        }
    }

    let mut out = String::new();
    if is_neg {
        out.push('-');
    }

    if (-4..6).contains(&exp) {
        if exp >= 0 {
            let int_len = (exp as usize) + 1;
            out.push_str(&digits[..int_len]);
            let frac = &digits[int_len..];
            let trimmed = frac.trim_end_matches('0');
            if !trimmed.is_empty() {
                out.push('.');
                out.push_str(trimmed);
            }
        } else {
            out.push('0');
            let zeros = (-exp - 1) as usize;
            let trimmed = digits.trim_end_matches('0');
            if !trimmed.is_empty() {
                out.push('.');
                for _ in 0..zeros {
                    out.push('0');
                }
                out.push_str(trimmed);
            }
        }
    } else {
        out.push(digits.chars().next().unwrap());
        let frac = &digits[1..];
        let trimmed = frac.trim_end_matches('0');
        if !trimmed.is_empty() {
            out.push('.');
            out.push_str(trimmed);
        }
        out.push('e');
        if exp >= 0 {
            out.push('+');
            if exp < 10 {
                out.push('0');
            }
            out.push_str(&exp.to_string());
        } else {
            out.push('-');
            let abs_exp = -exp;
            if abs_exp < 10 {
                out.push('0');
            }
            out.push_str(&abs_exp.to_string());
        }
    }
    out
}

fn should_drop_info(key: &[u8], plan: &RewritePlan) -> bool {
    if plan.recalc_ac && (key == b"AC" || key == b"AN" || key == b"AF") {
        return false;
    }
    if let Some(keep) = &plan.keep_info
        && !keep.iter().any(|k| k.as_slice() == key)
    {
        return true;
    }
    plan.drop_info.iter().any(|d| d.as_slice() == key)
}

fn parse_info_id(line: &[u8]) -> Option<&[u8]> {
    let body = line.strip_prefix(b"##INFO=<")?;
    let body = body.strip_suffix(b"\n").unwrap_or(body);
    let body = body.strip_suffix(b"\r").unwrap_or(body);
    let body = body.strip_suffix(b">").unwrap_or(body);
    let mut i = 0;
    while i < body.len() {
        let start = i;
        let mut in_quote = false;
        while i < body.len() {
            if body[i] == b'"' {
                in_quote = !in_quote;
            } else if body[i] == b',' && !in_quote {
                break;
            }
            i += 1;
        }
        let field = &body[start..i];
        if i < body.len() {
            i += 1;
        }
        if let Some(id) = field.strip_prefix(b"ID=") {
            let id = match id {
                [b'"', inner @ .., b'"'] => inner,
                v => v,
            };
            return Some(id);
        }
    }
    None
}

/// Rewrite a VCF header according to `plan`.
///
/// Pruned INFO keys have their `##INFO` lines removed.
/// If `recalc_ac` is true, declarations for `AC`, `AN`, `AF` are inserted if missing.
/// If sample selection is present, the `#CHROM` column line lists only the kept samples.
pub fn rewrite_header(h: &VcfHeader, plan: &RewritePlan) -> Vec<u8> {
    let mut out = Vec::with_capacity(h.raw.len());
    let defs = h.definitions();
    let mut added_recalc_headers = false;

    for line in h.raw.split_inclusive(|&b| b == b'\n') {
        if let Some(id) = parse_info_id(line)
            && should_drop_info(id, plan)
        {
            continue;
        }
        if line.starts_with(b"#CHROM") {
            if plan.recalc_ac && !added_recalc_headers {
                if !defs.contains_key(&("INFO", b"AC".to_vec())) {
                    out.extend_from_slice(b"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">\n");
                }
                if !defs.contains_key(&("INFO", b"AN".to_vec())) {
                    out.extend_from_slice(b"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n");
                }
                if !defs.contains_key(&("INFO", b"AF".to_vec())) {
                    out.extend_from_slice(b"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency for each ALT allele in the same order as listed\">\n");
                }
                added_recalc_headers = true;
            }
            if let Some(sel) = &plan.selection {
                let any_kept = sel.keep.iter().any(|&k| k);
                if any_kept {
                    out.extend_from_slice(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
                    for (i, &keep) in sel.keep.iter().enumerate() {
                        if keep {
                            out.push(b'\t');
                            if let Some(s) = h.samples.get(i) {
                                out.extend_from_slice(s);
                            }
                        }
                    }
                } else {
                    out.extend_from_slice(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
                }
                out.push(b'\n');
                continue;
            }
        }
        out.extend_from_slice(line);
    }
    out
}

fn mask_gt(gt: &[u8], out: &mut Vec<u8>) {
    let mut in_allele = false;
    for &b in gt {
        if b == b'/' || b == b'|' {
            if in_allele {
                out.push(b'.');
                in_allele = false;
            } else {
                out.push(b'.');
            }
            out.push(b);
        } else {
            in_allele = true;
        }
    }
    if in_allele || out.is_empty() {
        out.push(b'.');
    }
}

fn mask_sample(sample: &[u8], gt_idx: usize, out: &mut Vec<u8>) {
    let mut i = 0;
    let mut sub_idx = 0;
    while i < sample.len() {
        let start = i;
        while i < sample.len() && sample[i] != b':' {
            i += 1;
        }
        let field = &sample[start..i];
        if sub_idx == gt_idx {
            mask_gt(field, out);
        } else {
            out.extend_from_slice(field);
        }
        if i < sample.len() {
            out.push(b':');
            i += 1;
        }
        sub_idx += 1;
    }
}

fn rewrite_info(
    info_col: &[u8],
    alt_col: &[u8],
    plan: &RewritePlan,
    stats: Option<&SiteStats>,
    out: &mut Vec<u8>,
) -> Result<(), String> {
    let mut new_ac = None;
    let mut new_an = None;
    let mut new_af = None;

    if plan.recalc_ac {
        let s =
            stats.ok_or_else(|| "recalc_ac requested but stats were not provided".to_string())?;
        let n_alts = if alt_col == b"." || alt_col.is_empty() {
            0
        } else {
            alt_col.split(|&b| b == b',').count()
        };
        if n_alts > 0 {
            let mut ac_vals = Vec::with_capacity(n_alts);
            let mut af_vals = Vec::with_capacity(n_alts);
            for k in 1..=n_alts {
                let c = s.allele_counts.get(k).copied().unwrap_or(0);
                ac_vals.push(c.to_string());
                if s.called_alleles == 0 {
                    af_vals.push("0".to_string());
                } else {
                    af_vals.push(fmt_g6(c as f64 / s.called_alleles as f64));
                }
            }
            new_ac = Some(format!("AC={}", ac_vals.join(",")));
            new_af = Some(format!("AF={}", af_vals.join(",")));
        }
        new_an = Some(format!("AN={}", s.called_alleles));
    }

    let mut kept_items: Vec<Vec<u8>> = Vec::new();
    if info_col != b"." && !info_col.is_empty() {
        for item in info_col.split(|&b| b == b';') {
            if item.is_empty() {
                continue;
            }
            let key = match memchr::memchr(b'=', item) {
                Some(eq) => &item[..eq],
                None => item,
            };
            if should_drop_info(key, plan) {
                continue;
            }
            if plan.recalc_ac {
                if key == b"AC" {
                    if let Some(ac) = new_ac.take() {
                        kept_items.push(ac.into_bytes());
                    }
                    continue;
                }
                if key == b"AN" {
                    if let Some(an) = new_an.take() {
                        kept_items.push(an.into_bytes());
                    }
                    continue;
                }
                if key == b"AF" {
                    if let Some(af) = new_af.take() {
                        kept_items.push(af.into_bytes());
                    }
                    continue;
                }
            }
            kept_items.push(item.to_vec());
        }
    }

    if plan.recalc_ac {
        if let Some(ac) = new_ac.take() {
            kept_items.push(ac.into_bytes());
        }
        if let Some(an) = new_an.take() {
            kept_items.push(an.into_bytes());
        }
        if let Some(af) = new_af.take() {
            kept_items.push(af.into_bytes());
        }
    }

    if kept_items.is_empty() {
        out.push(b'.');
    } else {
        for (i, item) in kept_items.iter().enumerate() {
            if i > 0 {
                out.push(b';');
            }
            out.extend_from_slice(item);
        }
    }
    Ok(())
}

/// Rewrite a VCF text record line according to `plan`.
///
/// Unselected sample columns are dropped. Masked samples have their GT alleles
/// replaced with `.`, keeping ploidy, phase and other FORMAT subfields.
/// Pruned INFO keys are removed. If `recalc_ac` is true, AC, AN and AF are
/// updated or appended.
pub fn rewrite_text(
    line: &[u8],
    plan: &RewritePlan,
    masked: &[bool],
    stats: Option<&SiteStats>,
    out: &mut Vec<u8>,
) -> Result<(), String> {
    let mut raw = line;
    while let [rest @ .., b'\n' | b'\r'] = raw {
        raw = rest;
    }
    if plan.is_passthrough() {
        out.clear();
        out.extend_from_slice(raw);
        return Ok(());
    }

    let cols: Vec<&[u8]> = raw.split(|&b| b == b'\t').collect();
    if cols.len() < 8 {
        return Err(format!(
            "malformed VCF line: expected at least 8 columns, found {}",
            cols.len()
        ));
    }

    out.clear();
    // CHROM, POS, ID, REF, ALT, QUAL, FILTER
    for col in &cols[..7] {
        out.extend_from_slice(col);
        out.push(b'\t');
    }

    // INFO
    rewrite_info(cols[7], cols[4], plan, stats, out)?;

    if cols.len() == 8 {
        return Ok(());
    }

    // FORMAT and samples
    let format_col = cols[8];
    let sample_cols = &cols[9..];
    let n_samples = sample_cols.len();

    let any_kept = if let Some(sel) = &plan.selection {
        (0..n_samples).any(|i| sel.keep.get(i).copied().unwrap_or(true))
    } else {
        n_samples > 0
    };

    if !any_kept {
        return Ok(());
    }

    out.push(b'\t');
    out.extend_from_slice(format_col);

    let gt_idx = format_col.split(|&b| b == b':').position(|k| k == b"GT");

    for (i, &sample) in sample_cols.iter().enumerate() {
        let keep = plan
            .selection
            .as_ref()
            .is_none_or(|s| s.keep.get(i).copied().unwrap_or(true));
        if !keep {
            continue;
        }
        out.push(b'\t');
        let is_masked = masked.get(i).copied().unwrap_or(false);
        if is_masked {
            if let Some(idx) = gt_idx {
                mask_sample(sample, idx, out);
            } else {
                out.extend_from_slice(sample);
            }
        } else {
            out.extend_from_slice(sample);
        }
    }

    Ok(())
}
