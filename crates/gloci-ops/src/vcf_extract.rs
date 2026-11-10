//! `gloci vcf extract`: INFO and FORMAT tables and the 012 genotype matrix
//! (decision D73). The records pass the shared selection stage first.

use crate::vcf_select::{SelectStage, SelectedInput, VcfSelect, run_selection};
use gloci_io::{BcfText, IoError, Loc, UsageError};
use gloci_vcf::{
    InfoKey, InfoKind, RewritePlan, format_cell, format_declared, format_index, genotype_012,
    info_cell, info_kind, rewrite_text, split_columns,
};
use std::io::Write;

const COMMAND: &str = "vcf extract";

/// What `vcf extract` writes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ExtractMode {
    /// One row per site: `CHROM POS ID REF ALT` and one column per INFO key.
    Info(Vec<String>),
    /// FORMAT keys, one row per site (`wide`) or per site and sample (`long`).
    Format { keys: Vec<String>, long: bool },
    /// The 012 matrix, `.indv` and `.pos` files.
    Matrix012,
}

/// Options for `vcf extract`.
#[derive(Clone, Debug)]
pub struct ExtractOptions {
    pub select: VcfSelect,
    pub mode: ExtractMode,
}

/// Counts from an extract run.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct ExtractSummary {
    pub records_in: u64,
    pub sites_out: u64,
    /// Multi-allelic sites the 012 matrix leaves out.
    pub multiallelic_skipped: u64,
}

fn usage(msg: impl Into<String>) -> IoError {
    UsageError(msg.into()).into()
}

fn write_err(source: std::io::Error) -> IoError {
    IoError::Write {
        path: "output".to_string(),
        source,
    }
}

/// Reject an empty or duplicated key list before any input is read.
pub fn check_keys(flag: &str, keys: &[String]) -> Result<(), IoError> {
    if keys.is_empty() || keys.iter().any(String::is_empty) {
        return Err(usage(format!(
            "{flag} needs at least one key, such as {flag} DP; list keys separated by commas"
        )));
    }
    for (i, k) in keys.iter().enumerate() {
        if keys[..i].contains(k) {
            return Err(usage(format!("{flag} names {k} twice; list each key once")));
        }
    }
    Ok(())
}

/// The record `hit` as VCF text with the selection's sample choice and
/// masking applied: a BCF record is rendered by the converter `view -O vcf`
/// uses, so VCF and BCF input give identical lines.
struct Lines<'a> {
    bcf: Option<BcfText>,
    plan: &'a RewritePlan,
    rewritten: Vec<u8>,
}

impl<'a> Lines<'a> {
    fn line<'b>(
        &'b mut self,
        name: &str,
        hit: &crate::vcf_select::Hit<'b>,
    ) -> Result<&'b [u8], IoError> {
        let text: &[u8] = match self.bcf.as_mut() {
            Some(b) => b.render(hit.raw, hit.rec_no)?,
            None => hit.raw,
        };
        if self.plan.is_passthrough() {
            return Ok(text);
        }
        self.rewritten.clear();
        rewrite_text(text, self.plan, hit.masked, hit.stats, &mut self.rewritten).map_err(
            |msg| IoError::Parse {
                path: name.to_string(),
                loc: Loc::Record(hit.rec_no),
                msg,
            },
        )?;
        Ok(&self.rewritten)
    }
}

fn kept_samples(stage: &SelectStage<'_>) -> Vec<Vec<u8>> {
    let all = &stage.header().samples;
    match &stage.plan().selection {
        Some(sel) => all
            .iter()
            .zip(&sel.keep)
            .filter(|(_, k)| **k)
            .map(|(s, _)| s.clone())
            .collect(),
        None => all.clone(),
    }
}

fn warn_undeclared(flag: &str, section: &str, keys: &[String], declared: impl Fn(&str) -> bool) {
    for k in keys {
        if !declared(k) {
            log::warn!(
                "{flag} {k}: the header does not declare {section}/{k}, so every value is \
missing unless the records carry it; check the key's spelling"
            );
        }
    }
}

/// Write the first five columns of a site: `CHROM POS ID REF ALT`.
fn site_prefix(cols: &[&[u8]], out: &mut Vec<u8>) {
    for (i, c) in [cols[0], cols[1], cols[2], cols[3], cols[4]]
        .iter()
        .enumerate()
    {
        if i > 0 {
            out.push(b'\t');
        }
        out.extend_from_slice(c);
    }
}

/// Extract INFO or FORMAT values of the records of `input` as a TSV to `w`.
pub fn extract_table(
    input: SelectedInput,
    o: &ExtractOptions,
    w: &mut dyn Write,
) -> Result<ExtractSummary, IoError> {
    o.select.validate()?;
    let (info_keys, format_keys, long): (&[String], &[String], bool) = match &o.mode {
        ExtractMode::Info(k) => (k, &[], false),
        ExtractMode::Format { keys, long } => (&[], keys, *long),
        ExtractMode::Matrix012 => {
            return Err(usage("extract_table does not write the 012 matrix"));
        }
    };
    let base = RewritePlan::default();
    let (reader, mut chroms) = input.into_parts();
    let mut stage = SelectStage::new(COMMAND, reader.as_ref(), &o.select, &base, false)?;
    let name = stage.name().to_string();
    let header_raw = stage.header().raw.clone();
    let samples = kept_samples(&stage);
    let infos: Vec<InfoKey> = info_keys
        .iter()
        .map(|k| InfoKey {
            name: k.as_bytes().to_vec(),
            kind: info_kind(&header_raw, k.as_bytes()),
        })
        .collect();
    warn_undeclared("--info", "INFO", info_keys, |k| {
        info_kind(&header_raw, k.as_bytes()) != InfoKind::Undeclared
    });
    warn_undeclared("--format", "FORMAT", format_keys, |k| {
        format_declared(&header_raw, k.as_bytes())
    });

    let mut head = b"CHROM\tPOS\tID\tREF\tALT".to_vec();
    if !infos.is_empty() {
        for k in &infos {
            head.push(b'\t');
            head.extend_from_slice(&k.name);
        }
    } else if long {
        head.extend_from_slice(b"\tSAMPLE");
        for k in format_keys {
            head.push(b'\t');
            head.extend_from_slice(k.as_bytes());
        }
    } else {
        for s in &samples {
            for k in format_keys {
                head.push(b'\t');
                head.extend_from_slice(s);
                head.push(b':');
                head.extend_from_slice(k.as_bytes());
            }
        }
    }
    head.push(b'\n');
    w.write_all(&head).map_err(write_err)?;

    let mut lines = Lines {
        bcf: if stage.input_is_bcf() {
            Some(BcfText::for_header(&name, stage.header())?)
        } else {
            None
        },
        plan: &stage.plan().clone(),
        rewritten: Vec::new(),
    };
    let mut row = Vec::new();
    let mut sites_out = 0u64;
    let records_in = run_selection(reader, &mut chroms, &mut stage, |hit| {
        let line = lines.line(&name, &hit)?;
        let parse_err = |msg| IoError::Parse {
            path: name.clone(),
            loc: Loc::Record(hit.rec_no),
            msg,
        };
        let mut cols: Vec<&[u8]> = Vec::new();
        split_columns(line, &mut cols).map_err(parse_err)?;
        let indexes: Vec<Option<usize>> = if format_keys.is_empty() {
            Vec::new()
        } else {
            format_keys
                .iter()
                .map(|k| format_index(cols.get(8).copied().unwrap_or(b"."), k.as_bytes()))
                .collect()
        };
        let sample_col = |i: usize| cols.get(9 + i).copied().unwrap_or(b".");
        if long && !format_keys.is_empty() {
            for (i, sample) in samples.iter().enumerate() {
                row.clear();
                site_prefix(&cols, &mut row);
                row.push(b'\t');
                row.extend_from_slice(sample);
                for idx in &indexes {
                    row.push(b'\t');
                    row.extend_from_slice(format_cell(sample_col(i), *idx));
                }
                row.push(b'\n');
                w.write_all(&row).map_err(write_err)?;
            }
        } else {
            row.clear();
            site_prefix(&cols, &mut row);
            if format_keys.is_empty() {
                for k in &infos {
                    row.push(b'\t');
                    row.extend_from_slice(info_cell(cols[7], k));
                }
            } else {
                for i in 0..samples.len() {
                    for idx in &indexes {
                        row.push(b'\t');
                        row.extend_from_slice(format_cell(sample_col(i), *idx));
                    }
                }
            }
            row.push(b'\n');
            w.write_all(&row).map_err(write_err)?;
        }
        sites_out += 1;
        Ok(())
    })?;
    Ok(ExtractSummary {
        records_in,
        sites_out,
        multiallelic_skipped: 0,
    })
}

/// Write vcftools' `--012` files for the records of `input`: `matrix` (one
/// row per kept sample: its 0-based index, then one code per site), `indv`
/// (the kept sample names) and `pos` (`CHROM POS` of each site written).
/// Sites with more than one ALT allele are left out, as vcftools does.
pub fn extract_matrix(
    input: SelectedInput,
    o: &ExtractOptions,
    matrix: &mut dyn Write,
    indv: &mut dyn Write,
    pos: &mut dyn Write,
) -> Result<ExtractSummary, IoError> {
    o.select.validate()?;
    let base = RewritePlan::default();
    let (reader, mut chroms) = input.into_parts();
    let mut stage = SelectStage::new(COMMAND, reader.as_ref(), &o.select, &base, false)?;
    let name = stage.name().to_string();
    let samples = kept_samples(&stage);
    for s in &samples {
        indv.write_all(s).map_err(write_err)?;
        indv.write_all(b"\n").map_err(write_err)?;
    }
    // One column of codes per sample, filled site by site.
    let mut codes: Vec<Vec<i8>> = vec![Vec::new(); samples.len()];
    let mut lines = Lines {
        bcf: if stage.input_is_bcf() {
            Some(BcfText::for_header(&name, stage.header())?)
        } else {
            None
        },
        plan: &stage.plan().clone(),
        rewritten: Vec::new(),
    };
    let mut sites_out = 0u64;
    let mut multiallelic_skipped = 0u64;
    let records_in = run_selection(reader, &mut chroms, &mut stage, |hit| {
        let line = lines.line(&name, &hit)?;
        let parse_err = |msg| IoError::Parse {
            path: name.clone(),
            loc: Loc::Record(hit.rec_no),
            msg,
        };
        let mut cols: Vec<&[u8]> = Vec::new();
        split_columns(line, &mut cols).map_err(parse_err)?;
        let n_alt = if cols[4] == b"." || cols[4].is_empty() {
            0
        } else {
            cols[4].split(|&b| b == b',').count()
        };
        if n_alt > 1 {
            multiallelic_skipped += 1;
            return Ok(());
        }
        let gt_idx = format_index(cols.get(8).copied().unwrap_or(b"."), b"GT");
        for (i, column) in codes.iter_mut().enumerate() {
            let gt = format_cell(cols.get(9 + i).copied().unwrap_or(b"."), gt_idx);
            let code = genotype_012(gt, n_alt).map_err(|msg| IoError::Parse {
                path: name.clone(),
                loc: Loc::Record(hit.rec_no),
                msg: format!(
                    "sample {} at {}:{}: {msg}",
                    String::from_utf8_lossy(&samples[i]),
                    String::from_utf8_lossy(cols[0]),
                    String::from_utf8_lossy(cols[1])
                ),
            })?;
            column.push(code);
        }
        pos.write_all(cols[0]).map_err(write_err)?;
        pos.write_all(b"\t").map_err(write_err)?;
        pos.write_all(cols[1]).map_err(write_err)?;
        pos.write_all(b"\n").map_err(write_err)?;
        sites_out += 1;
        Ok(())
    })?;
    if multiallelic_skipped > 0 {
        log::warn!(
            "{name}: the 012 matrix leaves out {multiallelic_skipped} multi-allelic site(s), as \
vcftools does; split them first (bcftools norm -m-) to include their alleles"
        );
    }
    let mut row = Vec::new();
    for (i, column) in codes.iter().enumerate() {
        row.clear();
        row.extend_from_slice(i.to_string().as_bytes());
        for c in column {
            row.push(b'\t');
            row.extend_from_slice(c.to_string().as_bytes());
        }
        row.push(b'\n');
        matrix.write_all(&row).map_err(write_err)?;
    }
    Ok(ExtractSummary {
        records_in,
        sites_out,
        multiallelic_skipped,
    })
}
