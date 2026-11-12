//! `gloci check`: validate a file without stopping at the first problem and
//! summarise it on one screen (§7.4).
//!
//! The reader runs with `InvalidPolicy::Skip`, so unparseable records are
//! counted and the rest is still read. On top of that every record is
//! checked for order, for `end` past the chromosome length (with `-g`), for
//! BED12 block consistency, for a numeric bedGraph value and for the VCF
//! column count. Exact duplicate lines and blank lines are warnings (D25).

use gloci_core::{
    ChromId, ChromOrder, ChromTable, Genome, OrderChecker, OrderSource, OrderViolation, thousands,
};
use gloci_io::{
    AlignOptions, AlnSource, BGZF_EOF, Compression, Format, Header, InvalidPolicy, IoError,
    LEX_SORTED_NOTE, Loc, Opened, ReaderOptions, RecordBatch, SortOrder, UsageError, aln_ref,
    bed12_blocks, check_bcf_record, chrom_table_for, cram_eof_present, cram_first_container,
    detect_compression, next_batch_size, open_reader_or_alignments,
    open_reader_or_alignments_from_nonempty, read_cram_prelude, sort_hint,
};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

/// Problems listed per category.
const MAX_LISTED: usize = 10;

const CHROM_MAP_NOT_FOR_ALIGNMENTS: &str = "check describes an alignment file as it is; \
--chrom-map does not apply to BAM, CRAM or SAM";

fn chrom_map_usage() -> IoError {
    IoError::Usage(UsageError(CHROM_MAP_NOT_FOR_ALIGNMENTS.to_string()))
}

/// Is `e` the refusal of a SAM source to rename references (one message for all three formats)?
fn is_sam_chrom_map_refusal(e: &IoError) -> bool {
    matches!(e, IoError::Usage(u) if u.0.contains("--chrom-map cannot rename"))
}

/// What a full check of an alignment file adds to the summary (D129).
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct AlnSummary {
    /// Mapped reads (placed, flag 0x4 clear).
    pub mapped: u64,
    /// Unmapped reads that carry a reference and position.
    pub unmapped_placed: u64,
    /// Reads with no reference.
    pub unplaced: u64,
}

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
    /// Set for an alignment file read in full.
    pub aln: Option<AlnSummary>,
    /// A usage error found while checking this file (an option that does not apply to its
    /// format); the command stops with it, exit status 2.
    pub usage: Option<String>,
    /// The report of `check --quick`: only the header, the first record and the end of the file
    /// were looked at.
    pub quick: bool,
}

fn format_name(f: Format) -> &'static str {
    match f {
        Format::Bed => "BED",
        Format::BedGraph => "BEDGRAPH",
        Format::Wig => "WIG",
        Format::Vcf => "VCF",
        Format::Bcf => "BCF",
        Format::Gff3 => "GFF3",
        Format::Gtf => "GTF",
        Format::Bam => "BAM",
        Format::Cram => "CRAM",
        Format::Sam => "SAM",
        Format::Fasta => "FASTA",
        Format::Psl => "PSL",
        Format::Rmsk => "RMSK",
        Format::Bedpe => "BEDPE",
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
        OrderSource::BuiltinGenome(id) => format!("built-in genome {id}"),
        OrderSource::VcfHeader(n) => format!("##contig order of {n}"),
        OrderSource::AlignmentHeader(n) => format!("@SQ order of {n}"),
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
            aln: None,
            usage: None,
            quick: false,
        }
    }

    /// A file that cannot be checked because of an option that does not apply to it.
    pub fn usage_error(name: &str, msg: String) -> CheckReport {
        CheckReport {
            usage: Some(msg.clone()),
            ..CheckReport::quick_failed(name, msg)
        }
    }

    /// A file that failed `check --quick`.
    pub fn quick_failed(name: &str, msg: String) -> CheckReport {
        CheckReport {
            quick: true,
            ..CheckReport::unreadable(name, msg)
        }
    }

    pub fn has_errors(&self) -> bool {
        self.errors.count > 0
    }

    /// The summary of one file, ending with a newline.
    pub fn render(&self) -> String {
        if self.quick {
            if self.unreadable {
                let mut s = format!("{}: failed the quick check\n  errors:       1\n", self.name);
                for l in &self.errors.first {
                    s.push_str(&format!("    {l}\n"));
                }
                return s;
            }
            return format!(
                "{}: {}{}\n  quick check:  passed (header, first record, end of file)\n",
                self.name,
                format_name(self.format),
                if self.format == Format::Cram {
                    String::new()
                } else {
                    format!(", {}", compression_name(self.compression))
                }
            );
        }
        if self.unreadable {
            let mut s = format!("{}: cannot be checked\n  errors:       1\n", self.name);
            for l in &self.errors.first {
                s.push_str(&format!("    {l}\n"));
            }
            return s;
        }
        let mut s = format!(
            "{}: {}{}\n  records:      {}\n",
            self.name,
            format_name(self.format),
            if self.format == Format::Cram {
                String::new()
            } else {
                format!(", {}", compression_name(self.compression))
            },
            thousands(self.records)
        );
        if let Some(a) = &self.aln {
            s.push_str(&format!(
                "  mapped:       {}\n  unmapped placed: {}\n  unplaced:     {}\n",
                thousands(a.mapped),
                thousands(a.unmapped_placed),
                thousands(a.unplaced)
            ));
        }
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
        let label = if self.aln.is_some() {
            "sorted:"
        } else {
            "order:"
        };
        s.push_str(&format!("  {label:<13} {}\n", self.order));
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
    genome: Option<&(Genome, OrderSource)>,
    order: ChromOrder,
    ro: &ReaderOptions,
) -> Result<CheckReport, IoError> {
    let mut ro = ro.clone();
    ro.invalid = InvalidPolicy::Skip;
    // An alignment file is read as whole records (D129); everything else as intervals.
    let alignments = ro.accept.alignments.get_or_insert_with(|| AlignOptions {
        full_records: true,
        ..AlignOptions::default()
    });
    alignments.full_records = true;
    let opened = open_reader_or_alignments(path, &ro).map_err(|e| {
        if ro.chrom_map.is_some() && is_sam_chrom_map_refusal(&e) {
            chrom_map_usage()
        } else {
            e
        }
    })?;
    let mut reader = match opened {
        Opened::Alignments(src) => {
            if ro.chrom_map.is_some() {
                return Err(chrom_map_usage());
            }
            if genome.is_some() {
                static ONCE: std::sync::Once = std::sync::Once::new();
                ONCE.call_once(|| {
                    log::warn!(
                        "-g and --chrom-order do not apply to BAM, CRAM or SAM: the header's \
@SQ lines give the references and the sort order is the @HD SO"
                    )
                });
            }
            return check_alignments(path, src);
        }
        Opened::Records(r) => r,
    };
    let mut chroms = match genome {
        Some((g, source)) => ChromTable::closed(g, source.clone()),
        None => chrom_table_for(None, &[reader.as_ref()], order)?,
    };
    let genome_data = genome.map(|(g, _)| g);
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
                Format::Bcf => {
                    if let Err(msg) = check_bcf_record(raw) {
                        other_errs.add(m.loc, at(msg));
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
        aln: None,
        usage: None,
        quick: false,
    })
}

/// The 1-based position `pos` shows as; a read with none shows 0, as SAM does.
fn shown(pos: Option<u64>) -> u64 {
    pos.map_or(0, |p| p + 1)
}

/// Every record of the alignment file `src`, through the checks of D129.
fn check_alignments(path: &str, mut src: Box<dyn AlnSource>) -> Result<CheckReport, IoError> {
    let header = src.header().clone();
    let name = src.name().to_string();
    let layout = src.layout();
    let format = src.format();
    let ref_name = |i: usize| String::from_utf8_lossy(&header.refs[i].0).into_owned();
    let declared = match header.sort_order {
        SortOrder::Coordinate => Some("coordinate"),
        SortOrder::Queryname => Some("queryname"),
        SortOrder::Unsorted => Some("unsorted"),
        SortOrder::Unknown => None,
    };
    let coordinate = header.sort_order == SortOrder::Coordinate;

    let mut summary = AlnSummary::default();
    let mut counts = vec![0u64; header.refs.len()];
    let mut appearance: Vec<usize> = Vec::new();
    let mut records = 0u64;
    let (mut order_errs, mut other_errs, mut warns) =
        (Found::default(), Found::default(), Found::default());
    let mut prev: Option<(usize, Option<u64>)> = None;
    let mut seen_unplaced = false;
    let mut group: Option<(usize, Option<u64>)> = None;
    let mut seen: HashSet<u64> = HashSet::new();
    let mut fatal: Option<String> = None;
    let mut buf = Vec::new();
    loop {
        buf.clear();
        let n = match src.next_record(&mut buf) {
            Ok(Some(n)) => n,
            Ok(None) => break,
            Err(e) => {
                fatal = Some(match &e {
                    IoError::Parse { loc, msg, .. } => format!("{loc}: {msg}"),
                    other => other.to_string(),
                });
                break;
            }
        };
        records += 1;
        let r = match aln_ref(layout, &buf, &header) {
            Ok(r) => r,
            Err(msg) => {
                other_errs.add(n, format!("record {n}: {msg}"));
                continue;
            }
        };
        let (rid, pos, flag) = (r.ref_id(), r.pos(), r.flag());
        match rid {
            None => {
                summary.unplaced += 1;
                seen_unplaced = true;
                group = None;
            }
            Some(i) => {
                if flag & 0x4 != 0 {
                    summary.unmapped_placed += 1;
                } else {
                    summary.mapped += 1;
                }
                if counts[i] == 0 {
                    appearance.push(i);
                }
                counts[i] += 1;
                let at = |msg: String| format!("record {n} of {}: {msg}", ref_name(i));
                if coordinate {
                    if seen_unplaced {
                        order_errs.add(n, at("a placed read follows an unplaced one".to_string()));
                    } else if let Some((pr, pp)) = prev {
                        if i < pr {
                            order_errs.add(
                                n,
                                at(format!(
                                    "reference {} comes before {}, which came earlier in the file",
                                    ref_name(i),
                                    ref_name(pr)
                                )),
                            );
                        } else if i == pr && pos < pp {
                            order_errs.add(
                                n,
                                at(format!(
                                    "position {} is before the previous position {}",
                                    thousands(shown(pos)),
                                    thousands(shown(pp))
                                )),
                            );
                        }
                    }
                    prev = Some((i, pos));
                }
                let ref_len = header.refs[i].1;
                if let Some((_, true)) = r.span(ref_len) {
                    warns.add(
                        n,
                        at(format!(
                            "the alignment at position {} runs past the end of {} (length {})",
                            thousands(shown(pos)),
                            ref_name(i),
                            thousands(ref_len)
                        )),
                    );
                }
                if group != Some((i, pos)) {
                    group = Some((i, pos));
                    seen.clear();
                }
                if !seen.insert(hash_line(&buf)) {
                    warns.add(n, at("duplicate of an earlier record".to_string()));
                }
            }
        }
    }

    let diag = src.diagnostics();
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
    let mut add_error = |msg: String| {
        errors.count += 1;
        if errors.first.len() < MAX_LISTED {
            errors.first.push(msg);
        }
    };
    if let Some(msg) = fatal {
        add_error(msg);
    }
    // The end-of-file marker of a regular file. BGZF is checked on open; a CRAM's container here.
    let p = Path::new(path);
    let head = if path != "-" && p.is_file() {
        let mut h = Vec::new();
        File::open(p)
            .and_then(|f| f.take(18).read_to_end(&mut h))
            .ok();
        h
    } else {
        Vec::new()
    };
    if format == Format::Cram && !head.is_empty() && !cram_eof_present(p).unwrap_or(true) {
        add_error(format!(
            "{name} has no CRAM end-of-file container; the file is probably truncated"
        ));
    }
    let warnings = Problems {
        count: warns.count,
        first: warns.first.into_iter().map(|(_, m)| m).collect(),
    };

    let order_text = match (declared, order_errs.count) {
        (Some("coordinate"), 0) => "yes (@HD SO:coordinate)".to_string(),
        (Some("coordinate"), _) => "NOT sorted (@HD SO:coordinate)".to_string(),
        (Some(so), _) => format!("no (@HD SO:{so}; the order is not checked)"),
        (None, _) => "no (no @HD SO; the order is not checked)".to_string(),
    };
    let compression = match format {
        Format::Bam => Compression::Bgzf,
        Format::Cram => Compression::None,
        _ => detect_compression(&head),
    };
    Ok(CheckReport {
        name,
        format,
        compression,
        records,
        chroms: appearance
            .iter()
            .map(|&i| (header.refs[i].0.clone(), counts[i]))
            .collect(),
        order: order_text,
        order_hint: (order_errs.count > 0).then(|| {
            format!(
                "sort it with: gloci sort {path} -o sorted.{}",
                if format == Format::Sam { "sam" } else { "bam" }
            )
        }),
        unreadable: false,
        errors,
        warnings,
        aln: Some(summary),
        usage: None,
        quick: false,
    })
}

/// What `check --quick` is asked to accept.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct QuickOptions {
    /// Accept an alignment file whose header has no `@SQ` line (samtools quickcheck `-u`).
    pub unaligned: bool,
}

/// samtools quickcheck, for any format `check` reads (D129): the header parses, an alignment
/// file lists a reference (unless `q.unaligned`), the first record parses (for a CRAM, the
/// first container, so that no reference is needed), and a regular file ends with its BGZF or
/// CRAM end-of-file marker. Standard input is not checked for the marker, as htslib. Only the
/// header, the first record and the last 28 bytes are read of a BGZF file (the header, the first
/// container and the last 38 bytes of a CRAM). A problem is recorded in the report, which never
/// fails to be made; an option that does not apply to the format is a usage error in it.
pub fn quick_check(path: &str, o: &ReaderOptions, q: &QuickOptions) -> CheckReport {
    let r = if path == "-" {
        quick_open(path, Box::new(std::io::stdin()), o, q)
    } else if gloci_io::memio::is_mem_path(path) {
        // An in-memory input has no end-of-file marker to look at: it is read as a stream.
        gloci_io::open_raw(path)
            .map_err(|e| e.to_string())
            .and_then(|src| quick_open(&gloci_io::memio::display_name(path), src, o, q))
    } else {
        File::open(path)
            .map_err(|e| format!("cannot open {path}: {e}"))
            .and_then(|mut f| {
                let mut magic = [0u8; 4];
                let n = f
                    .read(&mut magic)
                    .map_err(|e| format!("cannot read {path}: {e}"))?;
                f.seek(SeekFrom::Start(0))
                    .map_err(|e| format!("cannot read {path}: {e}"))?;
                if &magic[..n] == b"CRAM" {
                    if o.chrom_map.is_some() {
                        return Err(CHROM_MAP_NOT_FOR_ALIGNMENTS.to_string());
                    }
                    return quick_cram(path, q).map(|()| (Format::Cram, Compression::None));
                }
                let regular = f.metadata().is_ok_and(|m| m.is_file());
                if regular {
                    quick_regular(path, f, o, q)
                } else {
                    quick_open(path, Box::new(f), o, q)
                }
            })
    };
    match r {
        Ok((format, compression)) => quick_ok(path, format, compression),
        Err(msg) if msg == CHROM_MAP_NOT_FOR_ALIGNMENTS => CheckReport::usage_error(path, msg),
        Err(msg) => CheckReport::quick_failed(path, msg),
    }
}

/// [`quick_check`] of a regular file already open as `r`, named `name`; reads nothing but the
/// header, the first record and the last 28 bytes. A CRAM needs a path ([`quick_check`]).
pub fn quick_check_reader<R: Read + Seek + Send + 'static>(
    name: &str,
    r: R,
    o: &ReaderOptions,
    q: &QuickOptions,
) -> CheckReport {
    match quick_regular(name, r, o, q) {
        Ok((format, compression)) => quick_ok(name, format, compression),
        Err(msg) if msg == CHROM_MAP_NOT_FOR_ALIGNMENTS => CheckReport::usage_error(name, msg),
        Err(msg) => CheckReport::quick_failed(name, msg),
    }
}

fn quick_ok(name: &str, format: Format, compression: Compression) -> CheckReport {
    CheckReport {
        name: name.to_string(),
        format,
        compression,
        records: 0,
        chroms: Vec::new(),
        order: String::new(),
        order_hint: None,
        unreadable: false,
        errors: Problems::default(),
        warnings: Problems::default(),
        aln: None,
        usage: None,
        quick: true,
    }
}

/// A BGZF file must end with the 28-byte end-of-file block, then it is opened.
fn quick_regular<R: Read + Seek + Send + 'static>(
    name: &str,
    mut r: R,
    o: &ReaderOptions,
    q: &QuickOptions,
) -> Result<(Format, Compression), String> {
    let read_err = |e: std::io::Error| format!("cannot read {name}: {e}");
    let len = r.seek(SeekFrom::End(0)).map_err(read_err)?;
    if len == 0 {
        return Err(IoError::empty_input(name).to_string());
    }
    r.seek(SeekFrom::Start(0)).map_err(read_err)?;
    let mut head = Vec::new();
    (&mut r).take(18).read_to_end(&mut head).map_err(read_err)?;
    if detect_compression(&head) == Compression::Bgzf {
        let mut tail = [0u8; 28];
        let ok = len >= 28 && {
            r.seek(SeekFrom::Start(len - 28)).map_err(read_err)?;
            r.read_exact(&mut tail).is_ok() && tail == BGZF_EOF
        };
        if !ok {
            return Err(IoError::Truncated {
                path: name.to_string(),
            }
            .to_string());
        }
    }
    r.seek(SeekFrom::Start(0)).map_err(read_err)?;
    let (format, _) = quick_open(name, Box::new(r), o, q)?;
    Ok((format, detect_compression(&head)))
}

/// The header and the first record of `src`. No bytes at all, before or after decompression,
/// fails: there is no header, though content detection would call it an empty BED (which the
/// full check accepts).
fn quick_open(
    name: &str,
    src: Box<dyn Read + Send>,
    o: &ReaderOptions,
    q: &QuickOptions,
) -> Result<(Format, Compression), String> {
    let mut src = BufReader::new(src);
    let empty = loop {
        match src.fill_buf() {
            Ok(b) => break b.is_empty(),
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
            Err(e) => return Err(format!("cannot read {name}: {e}")),
        }
    };
    let shown = if name == "-" { "stdin" } else { name };
    if empty {
        return Err(IoError::empty_input(shown).to_string());
    }
    let src: Box<dyn Read + Send> = Box::new(src);
    let mut o = o.clone();
    o.invalid = InvalidPolicy::Fail;
    // One thread reads one block at a time and no further.
    o.threads = std::num::NonZeroUsize::MIN;
    o.accept
        .alignments
        .get_or_insert_with(AlignOptions::default);
    let opened = open_reader_or_alignments_from_nonempty(name, src, &o).map_err(|e| {
        if o.chrom_map.is_some() && is_sam_chrom_map_refusal(&e) {
            CHROM_MAP_NOT_FOR_ALIGNMENTS.to_string()
        } else {
            e.to_string()
        }
    })?;
    // Nothing once decompressed (a BGZF end-of-file block alone, a gzip of nothing) is as
    // empty as no bytes at all.
    let Some(opened) = opened else {
        return Err(IoError::empty_input(shown).to_string());
    };
    match opened {
        Opened::Alignments(mut s) => {
            if o.chrom_map.is_some() {
                return Err(CHROM_MAP_NOT_FOR_ALIGNMENTS.to_string());
            }
            if s.header().refs.is_empty() && !q.unaligned {
                return Err(format!(
                    "{name} has no @SQ line in its header; if it holds unaligned reads, pass \
--unaligned"
                ));
            }
            let mut buf = Vec::new();
            if s.next_record(&mut buf)
                .map_err(|e| e.to_string())?
                .is_some()
            {
                aln_ref(s.layout(), &buf, s.header())
                    .map_err(|msg| format!("{name}: record 1: {msg}"))?;
            }
            let c = if s.format() == Format::Bam {
                Compression::Bgzf
            } else {
                Compression::None
            };
            Ok((s.format(), c))
        }
        Opened::Records(mut r) => {
            let mut chroms = chrom_table_for(None, &[r.as_ref()], ChromOrder::default())
                .map_err(|e| e.to_string())?;
            let mut batch = RecordBatch::default();
            r.read_batch(&mut chroms, &mut batch, 1)
                .map_err(|e| e.to_string())?;
            Ok((r.format(), r.compression()))
        }
    }
}

/// A CRAM: the file definition and header, a first container that reads, and the end-of-file
/// container. No reference is needed.
fn quick_cram(path: &str, q: &QuickOptions) -> Result<(), String> {
    let p = Path::new(path);
    let prelude = read_cram_prelude(p).map_err(|e| e.to_string())?;
    let has_sq = prelude
        .header_text
        .split(|&b| b == b'\n')
        .any(|l| l.starts_with(b"@SQ\t") || l == b"@SQ");
    if !has_sq && !q.unaligned {
        return Err(format!(
            "{path} has no @SQ line in its header; if it holds unaligned reads, pass --unaligned"
        ));
    }
    cram_first_container(p, prelude.data_start).map_err(|e| e.to_string())?;
    match cram_eof_present(p) {
        Ok(true) => Ok(()),
        Ok(false) => Err(format!(
            "{path} has no CRAM end-of-file container; the file is probably truncated"
        )),
        Err(e) => Err(format!("cannot read {path}: {e}")),
    }
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
