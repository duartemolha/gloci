//! VCF record and header rewriting for `vcf filter`.

use super::samples::SampleSelection;
use crate::stats::{MaskRules, SiteStats};
use gloci_io::VcfHeader;
use std::io::Write;

/// A plan describing how to rewrite VCF records and headers.
#[derive(Clone, Debug, Default, PartialEq)]
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
    let mut v = Vec::with_capacity(16);
    crate::stats::fmt::write_g6(&mut v, x);
    String::from_utf8(v).expect("ascii")
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

/// Append GT `gt` with every allele missing (D50): one `.` per allele, the
/// separators kept, and a VCF 4.4 leading phase indicator (`|0|1`) kept in
/// front. An empty GT becomes `.`.
fn mask_gt(gt: &[u8], out: &mut Vec<u8>) {
    let body = match gt.first() {
        Some(&p @ (b'|' | b'/')) => {
            out.push(p);
            &gt[1..]
        }
        _ => gt,
    };
    for &b in body {
        if b == b'/' || b == b'|' {
            out.push(b'.');
            out.push(b);
        }
    }
    out.push(b'.');
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

fn write_ac(out: &mut Vec<u8>, s: &SiteStats, n_alts: usize) {
    out.extend_from_slice(b"AC=");
    for k in 1..=n_alts {
        if k > 1 {
            out.push(b',');
        }
        let c = s.allele_counts.get(k).copied().unwrap_or(0);
        write!(out, "{c}").unwrap();
    }
}

fn write_an(out: &mut Vec<u8>, s: &SiteStats) {
    write!(out, "AN={}", s.called_alleles).unwrap();
}

fn write_af(out: &mut Vec<u8>, s: &SiteStats, n_alts: usize) {
    out.extend_from_slice(b"AF=");
    for k in 1..=n_alts {
        if k > 1 {
            out.push(b',');
        }
        if s.called_alleles == 0 {
            // No called allele: the frequency is unknown, not 0.
            out.push(b'.');
        } else {
            let c = s.allele_counts.get(k).copied().unwrap_or(0);
            let freq = c as f64 / s.called_alleles as f64;
            out.extend_from_slice(fmt_g6(freq).as_bytes());
        }
    }
}

fn rewrite_info(
    info_col: &[u8],
    alt_col: &[u8],
    plan: &RewritePlan,
    stats: Option<&SiteStats>,
    out: &mut Vec<u8>,
) -> Result<(), String> {
    let s = if plan.recalc_ac {
        Some(stats.ok_or_else(|| "recalc_ac requested but stats were not provided".to_string())?)
    } else {
        None
    };

    let n_alts = if alt_col == b"." || alt_col.is_empty() {
        0
    } else {
        alt_col.split(|&b| b == b',').count()
    };

    // If n_alts == 0, AC and AF are not applicable
    let mut wrote_ac = !plan.recalc_ac || n_alts == 0;
    let mut wrote_an = !plan.recalc_ac;
    let mut wrote_af = !plan.recalc_ac || n_alts == 0;
    let mut written_any = false;

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
                    if !wrote_ac {
                        if written_any {
                            out.push(b';');
                        }
                        write_ac(out, s.unwrap(), n_alts);
                        wrote_ac = true;
                        written_any = true;
                    }
                    continue;
                }
                if key == b"AN" {
                    if !wrote_an {
                        if written_any {
                            out.push(b';');
                        }
                        write_an(out, s.unwrap());
                        wrote_an = true;
                        written_any = true;
                    }
                    continue;
                }
                if key == b"AF" {
                    if !wrote_af {
                        if written_any {
                            out.push(b';');
                        }
                        write_af(out, s.unwrap(), n_alts);
                        wrote_af = true;
                        written_any = true;
                    }
                    continue;
                }
            }
            if written_any {
                out.push(b';');
            }
            out.extend_from_slice(item);
            written_any = true;
        }
    }

    if plan.recalc_ac {
        let stats_ref = s.unwrap();
        if !wrote_ac {
            if written_any {
                out.push(b';');
            }
            write_ac(out, stats_ref, n_alts);
            written_any = true;
        }
        if !wrote_an {
            if written_any {
                out.push(b';');
            }
            write_an(out, stats_ref);
            written_any = true;
        }
        if !wrote_af {
            if written_any {
                out.push(b';');
            }
            write_af(out, stats_ref, n_alts);
            written_any = true;
        }
    }

    if !written_any {
        out.push(b'.');
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

    let mut iter = raw.split(|&b| b == b'\t');
    let chrom = match iter.next() {
        Some(c) => c,
        None => return Err("malformed VCF line: expected at least 8 columns".to_string()),
    };
    let pos = match iter.next() {
        Some(c) => c,
        None => return Err("malformed VCF line: expected at least 8 columns".to_string()),
    };
    let id = match iter.next() {
        Some(c) => c,
        None => return Err("malformed VCF line: expected at least 8 columns".to_string()),
    };
    let ref_col = match iter.next() {
        Some(c) => c,
        None => return Err("malformed VCF line: expected at least 8 columns".to_string()),
    };
    let alt_col = match iter.next() {
        Some(c) => c,
        None => return Err("malformed VCF line: expected at least 8 columns".to_string()),
    };
    let qual = match iter.next() {
        Some(c) => c,
        None => return Err("malformed VCF line: expected at least 8 columns".to_string()),
    };
    let filter = match iter.next() {
        Some(c) => c,
        None => return Err("malformed VCF line: expected at least 8 columns".to_string()),
    };
    let info = match iter.next() {
        Some(c) => c,
        None => return Err("malformed VCF line: expected at least 8 columns".to_string()),
    };

    out.clear();
    // CHROM, POS, ID, REF, ALT, QUAL, FILTER
    for col in [chrom, pos, id, ref_col, alt_col, qual, filter] {
        out.extend_from_slice(col);
        out.push(b'\t');
    }

    // INFO
    rewrite_info(info, alt_col, plan, stats, out)?;

    let Some(format_col) = iter.next() else {
        return Ok(());
    };

    // If sample selection was provided and selects 0 samples, omit FORMAT.
    if plan.selection.as_ref().is_some_and(|s| s.count() == 0) {
        return Ok(());
    }

    let info_end = out.len();
    out.push(b'\t');
    out.extend_from_slice(format_col);

    let gt_idx = format_col.split(|&b| b == b':').position(|k| k == b"GT");

    let mut samples_kept = 0;
    for (i, sample) in iter.enumerate() {
        let keep = plan
            .selection
            .as_ref()
            .is_none_or(|s| s.keep.get(i).copied().unwrap_or(true));
        if !keep {
            continue;
        }
        samples_kept += 1;
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

    if samples_kept == 0 {
        out.truncate(info_end);
    }

    Ok(())
}
