//! `gloci check`: validate a file without stopping at the first problem and
//! summarise it on one screen (§7.4).
//!
//! The reader runs with `InvalidPolicy::Skip`, so unparseable records are
//! counted and the rest is still read. On top of that every record is
//! checked for order, for `end` past the chromosome length (with `-g`), for
//! BED12 block consistency, for a numeric bedGraph value and for the VCF
//! column count. Exact duplicate lines and blank lines are warnings (D25).

use gloci_core::{ChromId, ChromOrder, OrderChecker, OrderSource, OrderViolation, thousands};
use gloci_io::{
    Compression, Format, Header, IoError, LEX_SORTED_NOTE, Loc, ReaderOptions, RecordBatch,
    bed12_blocks, chrom_table_for, load_genome, next_batch_size, open_reader, sort_hint,
};
use std::collections::{HashMap, HashSet};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::path::Path;

/// Problems listed per category.
const MAX_LISTED: usize = 10;

/// A count of problems and the first few, each `"line N: msg"`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Problems {
    pub count: u64,
    /// At most 10.
    pub first: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CheckReport {
    pub name: String,
    pub format: Format,
    pub compression: Compression,
    pub records: u64,
    /// Chromosomes in order of first appearance, with their record counts.
    pub chroms: Vec<(Vec<u8>, u64)>,
    /// `sorted (karyotype order)` or `NOT sorted`.
    pub order: String,
    /// How to fix a file that is NOT sorted (one line).
    pub order_hint: Option<String>,
    /// The file could not be opened or read at all; `errors` says why.
    pub unreadable: bool,
    pub errors: Problems,
    pub warnings: Problems,
}

fn format_name(f: Format) -> &'static str {
    match f {
        Format::Bed => "BED",
        Format::BedGraph => "BEDGRAPH",
        Format::Vcf => "VCF",
        Format::Bcf => "BCF",
        Format::Gff3 => "GFF3",
        Format::Gtf => "GTF",
    }
}

fn compression_name(c: Compression) -> &'static str {
    match c {
        Compression::None => "plain",
        Compression::Gzip => "gzip",
        Compression::Bgzf => "BGZF",
    }
}

fn order_source_name(s: &OrderSource) -> String {
    match s {
        OrderSource::GenomeFile(p) => format!("genome file {}", p.display()),
        OrderSource::VcfHeader(n) => format!("##contig order of {n}"),
        OrderSource::Lexicographic => "lexicographic order".to_string(),
        OrderSource::Karyotype => "karyotype order".to_string(),
        OrderSource::FirstAppearance => "first appearance".to_string(),
    }
}

impl CheckReport {
    /// A report for a file that could not be checked.
    pub fn unreadable(name: &str, msg: String) -> CheckReport {
        CheckReport {
            name: name.to_string(),
            format: Format::Bed,
            compression: Compression::None,
            records: 0,
            chroms: Vec::new(),
            order: String::new(),
            order_hint: None,
            unreadable: true,
            errors: Problems {
                count: 1,
                first: vec![msg],
            },
            warnings: Problems::default(),
        }
    }

    pub fn has_errors(&self) -> bool {
        self.errors.count > 0
    }

    /// The summary of one file, ending with a newline.
    pub fn render(&self) -> String {
        if self.unreadable {
            let mut s = format!("{}: cannot be checked\n  errors:       1\n", self.name);
            for l in &self.errors.first {
                s.push_str(&format!("    {l}\n"));
            }
            return s;
        }
        let mut s = format!(
            "{}: {}, {}\n  records:      {}\n",
            self.name,
            format_name(self.format),
            compression_name(self.compression),
            thousands(self.records)
        );
        let n = self.chroms.len();
        let names: Vec<String> = self
            .chroms
            .iter()
            .take(3)
            .map(|(c, _)| String::from_utf8_lossy(c).into_owned())
            .collect();
        let more = if n > 3 {
            format!(", … {} more", thousands((n - 3) as u64))
        } else {
            String::new()
        };
        if n == 0 {
            s.push_str("  chromosomes:  0\n");
        } else {
            s.push_str(&format!(
                "  chromosomes:  {} ({}{more})\n",
                thousands(n as u64),
                names.join(", ")
            ));
        }
        s.push_str(&format!("  order:        {}\n", self.order));
        if let Some(h) = &self.order_hint {
            s.push_str(&format!("                hint: {h}\n"));
        }
        for (label, p) in [
            ("errors:      ", &self.errors),
            ("warnings:    ", &self.warnings),
        ] {
            s.push_str(&format!("  {label} {}\n", thousands(p.count)));
            for l in &p.first {
                s.push_str(&format!("    {l}\n"));
            }
        }
        s
    }
}

/// The record number a listed problem starts with, to merge lists in file order.
fn problem_line(p: &str) -> u64 {
    let rest = p
        .strip_prefix("line ")
        .or_else(|| p.strip_prefix("record "))
        .unwrap_or("");
    let digits: String = rest
        .chars()
        .take_while(|c| c.is_ascii_digit() || *c == ',')
        .filter(|c| *c != ',')
        .collect();
    digits.parse().unwrap_or(u64::MAX)
}

/// Problems found while reading, kept in discovery order.
#[derive(Default)]
struct Found {
    count: u64,
    first: Vec<(u64, String)>,
}

impl Found {
    fn add(&mut self, loc: u64, msg: String) {
        self.count += 1;
        if self.first.len() < MAX_LISTED {
            self.first.push((loc, msg));
        }
    }
}

/// Is `v` an integer or decimal number with an optional sign and exponent?
fn is_number(v: &[u8]) -> bool {
    let v = v.strip_prefix(b"-").or(v.strip_prefix(b"+")).unwrap_or(v);
    let (mant, exp) = match v.iter().position(|b| *b == b'e' || *b == b'E') {
        Some(i) => (&v[..i], Some(&v[i + 1..])),
        None => (v, None),
    };
    let (int, frac) = match mant.iter().position(|b| *b == b'.') {
        Some(i) => (&mant[..i], &mant[i + 1..]),
        None => (mant, &mant[..0]),
    };
    let digits = |d: &[u8]| d.iter().all(u8::is_ascii_digit);
    if !digits(int) || !digits(frac) || (int.is_empty() && frac.is_empty()) {
        return false;
    }
    match exp {
        None => true,
        Some(e) => {
            let e = e.strip_prefix(b"-").or(e.strip_prefix(b"+")).unwrap_or(e);
            !e.is_empty() && digits(e)
        }
    }
}

fn hash_line(l: &[u8]) -> u64 {
    let mut h = DefaultHasher::new();
    l.hash(&mut h);
    h.finish()
}

fn loc_of(binary: bool, n: u64) -> Loc {
    if binary { Loc::Record(n) } else { Loc::Line(n) }
}

/// Check `path`, reading all of it. `genome` (`-g`) supplies chromosome
/// order, names and lengths; without it `order` is the fallback order.
pub fn check_file(
    path: &str,
    genome: Option<&Path>,
    order: ChromOrder,
    ro: &ReaderOptions,
) -> Result<CheckReport, IoError> {
    let mut ro = ro.clone();
    ro.invalid = gloci_io::InvalidPolicy::Skip;
    let mut reader = open_reader(path, &ro)?;
    let genome_data = match genome {
        Some(p) => Some(load_genome(p)?),
        None => None,
    };
    let mut chroms = chrom_table_for(genome, &[reader.as_ref()], order)?;
    let format = reader.format();
    let binary = format == Format::Bcf;
    let expected_cols = match reader.header() {
        Header::Vcf(h) if format == Format::Vcf => Some(if h.samples.is_empty() {
            8
        } else {
            9 + h.samples.len()
        }),
        _ => None,
    };
    let source = chroms.order_source().clone();

    let mut checker = OrderChecker::new();
    let mut counts: Vec<u64> = Vec::new();
    let mut appearance: Vec<ChromId> = Vec::new();
    let mut lengths: Vec<Option<Option<u64>>> = Vec::new();
    let genome_lens: Option<HashMap<Vec<u8>, Option<u64>>> = genome_data.as_ref().map(|g| {
        (0..g.len())
            .map(|k| (g.name(k).to_vec(), g.length(k)))
            .collect()
    });
    let mut records = 0u64;
    let (mut order_errs, mut other_errs, mut dups) =
        (Found::default(), Found::default(), Found::default());
    let mut lex_looking = false;
    let mut group: Option<(ChromId, u64)> = None;
    let mut seen: HashSet<u64> = HashSet::new();
    let mut fatal: Option<String> = None;

    let mut batch = RecordBatch::default();
    let mut max = 1;
    loop {
        batch.clear();
        let more = match reader.read_batch(&mut chroms, &mut batch, max) {
            Ok(m) => m,
            Err(e) => {
                fatal = Some(match &e {
                    IoError::Parse { loc, msg, .. } => format!("{loc}: {msg}"),
                    other => other.to_string(),
                });
                false
            }
        };
        for (i, m) in batch.recs.iter().enumerate() {
            let raw = batch.raw(i);
            let at = |msg: String| format!("{}: {msg}", loc_of(binary, m.loc));
            records += 1;

            let idx = m.chrom.0 as usize;
            if counts.len() <= idx {
                counts.resize(idx + 1, 0);
            }
            if counts[idx] == 0 {
                appearance.push(m.chrom);
            }
            counts[idx] += 1;

            if let Err(v) = checker.check(&chroms, m.chrom, m.key) {
                let msg = match &v {
                    OrderViolation::ChromBackwards { prev, next } => format!(
                        "chromosome {} is out of order: it sorts before {}, which came earlier",
                        String::from_utf8_lossy(chroms.name(*next)),
                        String::from_utf8_lossy(chroms.name(*prev))
                    ),
                    OrderViolation::StartBackwards { chrom, prev, next } => format!(
                        "{} start {} is before the previous start {}",
                        String::from_utf8_lossy(chroms.name(*chrom)),
                        thousands(*next),
                        thousands(*prev)
                    ),
                };
                if let OrderViolation::ChromBackwards { prev, next } = v
                    && source == OrderSource::Karyotype
                    && chroms.name(prev) < chroms.name(next)
                {
                    lex_looking = true;
                }
                order_errs.add(m.loc, at(msg));
            }

            if let Some(g) = &genome_lens {
                if lengths.len() <= idx {
                    lengths.resize(idx + 1, None);
                }
                let len = *lengths[idx]
                    .get_or_insert_with(|| g.get(chroms.name(m.chrom)).copied().flatten());
                if let Some(len) = len
                    && m.span.end > len
                {
                    other_errs.add(
                        m.loc,
                        at(format!(
                            "end {} is past the end of {} (length {})",
                            m.span.end,
                            String::from_utf8_lossy(chroms.name(m.chrom)),
                            len
                        )),
                    );
                }
            }

            match format {
                Format::Bed => {
                    if raw.iter().filter(|&&b| b == b'\t').count() >= 11
                        && let Err(msg) = bed12_blocks(raw)
                    {
                        other_errs.add(m.loc, at(msg));
                    }
                }
                Format::BedGraph => {
                    let v = raw.split(|&b| b == b'\t').nth(3).unwrap_or(b"");
                    if !is_number(v) {
                        other_errs.add(
                            m.loc,
                            at(format!(
                                "value \"{}\" is not a number",
                                String::from_utf8_lossy(v)
                            )),
                        );
                    }
                }
                Format::Vcf => {
                    if let Some(want) = expected_cols {
                        let found = raw.iter().filter(|&&b| b == b'\t').count() + 1;
                        if found != want {
                            other_errs
                                .add(m.loc, at(format!("expected {want} columns, found {found}")));
                        }
                    }
                }
                _ => {}
            }

            if group != Some((m.chrom, m.key)) {
                group = Some((m.chrom, m.key));
                seen.clear();
            }
            if !seen.insert(hash_line(raw)) {
                dups.add(m.loc, at("duplicate of an earlier line".to_string()));
            }
        }
        max = next_batch_size(&batch, max);
        if !more {
            break;
        }
    }

    // Parse problems, merged with the others in file order.
    let diag = reader.diagnostics();
    let mut errs: Vec<(u64, String)> = diag
        .examples
        .iter()
        .map(|e| (problem_line(e), e.clone()))
        .collect();
    errs.extend(order_errs.first);
    errs.extend(other_errs.first);
    errs.sort_by_key(|(l, _)| *l);
    errs.truncate(MAX_LISTED);
    let mut errors = Problems {
        count: diag.skipped + order_errs.count + other_errs.count,
        first: errs.into_iter().map(|(_, m)| m).collect(),
    };
    if let Some(msg) = fatal {
        errors.count += 1;
        if errors.first.len() < MAX_LISTED {
            errors.first.push(msg);
        }
    }

    let mut warnings = Problems {
        count: dups.count,
        first: dups.first.into_iter().map(|(_, m)| m).collect(),
    };
    // Blank lines are dropped by the text readers, which count them.
    warnings.count += diag.blank_lines;
    warnings.first.extend(
        diag.blank_examples
            .iter()
            .map(|n| format!("{}: blank line", Loc::Line(*n))),
    );
    warnings.first.sort_by_key(|m| problem_line(m));
    warnings.first.truncate(MAX_LISTED);

    let order_text = if order_errs.count == 0 {
        format!("sorted ({})", order_source_name(&source))
    } else {
        "NOT sorted".to_string()
    };
    Ok(CheckReport {
        name: reader.name().to_string(),
        format,
        compression: reader.compression(),
        records,
        chroms: appearance
            .iter()
            .map(|&id| (chroms.name(id).to_vec(), counts[id.0 as usize]))
            .collect(),
        order: order_text,
        order_hint: (order_errs.count > 0).then(|| {
            let mut h = sort_hint(reader.name(), &source);
            if lex_looking {
                h.push_str(LEX_SORTED_NOTE);
            }
            h
        }),
        unreadable: false,
        errors,
        warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numbers() {
        for ok in ["1", "-2", "+3.5", ".5", "5.", "1e5", "1.5E-3", "0"] {
            assert!(is_number(ok.as_bytes()), "{ok}");
        }
        for bad in ["", "abc", "1.2.3", "e5", "1e", "-", ".", "1,5", "nan"] {
            assert!(!is_number(bad.as_bytes()), "{bad}");
        }
    }

    #[test]
    fn problem_lines() {
        assert_eq!(problem_line("line 1,234: x"), 1234);
        assert_eq!(problem_line("record 7: x"), 7);
        assert_eq!(problem_line("other"), u64::MAX);
    }
}
