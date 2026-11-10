//! `gloci vcf concat`: write VCF/BCF files one after another (D74).
//!
//! The inputs must have identical samples and compatible header
//! definitions (the rule `sort` applies to several VCFs), and must follow
//! each other in the order of the chromosome table: the last record of one
//! file sorts at or before the first record of the next. Otherwise it is an
//! error with a `gloci sort` hint.
//!
//! Checked mode streams every record through the normal reader and writer,
//! so the output format follows `-o`/`-O` and inputs may mix VCF, VCF.gz and
//! BCF. `--naive` instead copies the BGZF blocks of identical-header
//! BGZF VCFs without recompressing them (see `gloci_vcf::concat`). It checks
//! the order across files from each file's first and last record but, unlike
//! checked mode, trusts the order inside each file.
//!
//! Neither mode adds a `##glociCommand` line: the output header is the
//! first file's, so both modes write the same bytes (before compression).

use crate::sort::{HeaderRuleCmd, check_definitions_for, check_samples_for};
use gloci_core::{ChromId, ChromOrder, ChromTable};
use gloci_io::{
    AtomicFile, BGZF_EOF, BcfDicts, BcfText, Compression, Format, Header, IoError, OrderGuard,
    OutFormat, OutputOptions, OutputSpec, ReaderOptions, RecordBatch, RecordReader, TextToBcf,
    UsageError, WriterOptions, chrom_table_lenient, next_batch_size, open_output, open_reader,
    open_reader_from, open_record_writer, refuse_binary_to_tty, sort_hint,
};
use gloci_vcf::concat::{NaiveError, NaiveFile, copy_blocks};
use std::fs::File;
use std::io::{self, IsTerminal, Write};

#[derive(Debug, Clone)]
pub struct ConcatOptions {
    /// The order of chromosomes when the headers declare none.
    pub chrom_order: ChromOrder,
    /// Copy BGZF blocks instead of reading records.
    pub naive: bool,
}

fn usage(msg: impl Into<String>) -> IoError {
    UsageError(msg.into()).into()
}

/// Longest record excerpt shown in an error, in characters.
const SHOW_CHARS: usize = 120;

fn show(bytes: &[u8]) -> String {
    let s = String::from_utf8_lossy(bytes);
    match s.char_indices().nth(SHOW_CHARS) {
        Some((cut, _)) => format!("{}...", &s[..cut]),
        None => s.into_owned(),
    }
}

/// The last record of the previous non-empty file, for the order check at
/// the start of the next one.
#[derive(Default)]
struct Boundary {
    prev: Option<Prev>,
}

struct Prev {
    file: String,
    chrom: ChromId,
    key: u64,
    shown: String,
}

impl Boundary {
    /// The first record of `file` (at `loc`, shown as `shown`) must not sort
    /// before the last record of the file before it.
    fn check_first(
        &self,
        chroms: &ChromTable,
        file: &str,
        chrom: ChromId,
        key: u64,
        loc: &str,
        shown: &str,
    ) -> Result<(), IoError> {
        let Some(p) = &self.prev else { return Ok(()) };
        let before = match chroms.cmp(chrom, p.chrom) {
            std::cmp::Ordering::Less => true,
            std::cmp::Ordering::Equal => key < p.key,
            std::cmp::Ordering::Greater => false,
        };
        if !before {
            return Ok(());
        }
        let hint = sort_hint(&format!("{} {file}", p.file), chroms.order_source());
        Err(IoError::Detect {
            path: file.to_string(),
            msg: format!(
                "{file} is not sorted after {}: {loc} ({shown}) comes before the last record of {} ({}); \
vcf concat writes the files one after another, so each must start where the one before ends\nhint: {hint}",
                p.file, p.file, p.shown
            ),
        })
    }

    fn set_last(&mut self, file: &str, chrom: ChromId, key: u64, shown: String) {
        self.prev = Some(Prev {
            file: file.to_string(),
            chrom,
            key,
            shown,
        });
    }
}

/// Open `path`. With `lenient`, a BGZF file without its EOF marker is
/// accepted (`--naive` warns about it instead).
fn open(path: &str, ro: &ReaderOptions, lenient: bool) -> Result<Box<dyn RecordReader>, IoError> {
    if !lenient || path == "-" {
        return open_reader(path, ro);
    }
    let f = File::open(path).map_err(|source| IoError::Open {
        path: path.to_string(),
        source,
    })?;
    open_reader_from(path, Box::new(f), ro)
}

pub fn vcf_concat(
    files: &[String],
    opts: &ConcatOptions,
    ro: &ReaderOptions,
    spec: &OutputSpec,
    oo: &OutputOptions,
    wo: &WriterOptions,
) -> Result<(), IoError> {
    if files.is_empty() {
        return Err(usage("vcf concat needs at least one input file"));
    }
    if opts.naive {
        if oo.write_index {
            return Err(usage(
                "--write-index cannot be combined with --naive; index the output afterwards with: gloci index OUT.vcf.gz",
            ));
        }
        if spec.format.is_some_and(|f| f != OutFormat::Vcf) || (spec.format.is_some() && !spec.bgzf)
        {
            return Err(usage(
                "--naive writes BGZF-compressed VCF; use an output name ending in .vcf.gz (or -O vcf.gz), or leave out --naive",
            ));
        }
        if files.iter().any(|f| f == "-") {
            return Err(usage(
                "--naive reads whole files and cannot read stdin; name the files, or leave out --naive",
            ));
        }
    }
    if opts.naive && ro.chrom_map.is_some() {
        return Err(usage(
            "--chrom-map cannot be combined with --naive, which copies BGZF blocks without changing them; leave out --naive (checked mode renames the contigs), or rename the files first with: gloci view FILE --chrom-map MAP -o OUT.vcf.gz",
        ));
    }
    let ro = ro.clone();
    let mut readers: Vec<Box<dyn RecordReader>> = Vec::new();
    for f in files {
        let r = open(f, &ro, opts.naive)?;
        if !matches!(r.format(), Format::Vcf | Format::Bcf) {
            return Err(usage(format!(
                "{f} is not a VCF or BCF file; vcf concat joins VCF and BCF files"
            )));
        }
        readers.push(r);
    }
    if opts.naive {
        for r in &readers {
            let (name, what) = (r.name(), "--naive copies BGZF blocks, so");
            if r.format() == Format::Bcf {
                return Err(usage(format!(
                    "{what} it cannot join the BCF file {name}: the order across BCF files cannot be checked without reading their records; leave out --naive (checked mode handles BCF)"
                )));
            }
            if r.compression() != Compression::Bgzf {
                return Err(usage(format!(
                    "{what} every input must be BGZF-compressed VCF, but {name} is not; compress it with: gloci view {name} -o {name}.gz, or leave out --naive"
                )));
            }
        }
    }
    check_samples_for(HeaderRuleCmd::CONCAT, &readers)?;
    check_definitions_for(HeaderRuleCmd::CONCAT, &readers)?;
    let refs: Vec<&dyn RecordReader> = readers.iter().map(|r| r.as_ref()).collect();
    let chroms = chrom_table_lenient(None, &refs, opts.chrom_order)?;
    if opts.naive {
        naive(files, readers, chroms, spec)
    } else {
        checked(readers, chroms, spec, oo, wo)
    }
}

/// How a file's records reach an output of the chosen format.
enum Conv {
    /// Written as read.
    Raw,
    /// A BCF record rendered as a VCF line.
    ToText(Box<BcfText>),
    /// A VCF line encoded as a BCF record.
    ToBcf(Box<TextToBcf>),
    /// A BCF record with other dictionaries: rendered, then encoded again.
    Recode(Box<BcfText>, Box<TextToBcf>),
}

fn checked(
    mut readers: Vec<Box<dyn RecordReader>>,
    mut chroms: ChromTable,
    spec: &OutputSpec,
    oo: &OutputOptions,
    wo: &WriterOptions,
) -> Result<(), IoError> {
    let first_format = readers[0].format();
    let out_format = spec.format.unwrap_or(OutFormat::family(first_format));
    if !matches!(out_format, OutFormat::Vcf | OutFormat::Bcf) {
        return Err(gloci_io::family_mismatch(first_format, out_format).into());
    }
    refuse_binary_to_tty(
        spec.path.is_none() && io::stdout().is_terminal(),
        out_format,
        spec.bgzf,
    )?;
    let Header::Vcf(first_header) = readers[0].header().clone() else {
        return Err(IoError::Internal(
            "a VCF/BCF input has no VCF header".into(),
        ));
    };
    let dest = spec
        .path
        .as_ref()
        .map_or_else(|| "<stdout>".to_string(), |p| p.display().to_string());
    let named = |e: IoError| e.with_output_name(&dest);

    // The dictionaries BCF output is encoded against, when it is BCF.
    let out_dicts = if out_format == OutFormat::Bcf {
        Some(match &first_header.bcf {
            Some(d) => d.clone(),
            None => BcfDicts::from_header(&first_header.raw)
                .map_err(|e| IoError::Internal(format!("{}: {e}", readers[0].name())))?,
        })
    } else {
        None
    };
    let mut convs = Vec::with_capacity(readers.len());
    for r in &readers {
        let Header::Vcf(h) = r.header() else {
            return Err(IoError::Internal(
                "a VCF/BCF input has no VCF header".into(),
            ));
        };
        let from_bcf = r.format() == Format::Bcf;
        convs.push(match (&out_dicts, from_bcf) {
            (None, false) => Conv::Raw,
            (None, true) => Conv::ToText(Box::new(
                BcfText::for_reader(r.as_ref())?.expect("a BCF reader"),
            )),
            (Some(_), false) => Conv::ToBcf(Box::new(TextToBcf::new(r.name(), &first_header)?)),
            (Some(d), true) if h.bcf.as_ref() == Some(d) => Conv::Raw,
            (Some(_), true) => Conv::Recode(
                Box::new(BcfText::for_reader(r.as_ref())?.expect("a BCF reader")),
                Box::new(TextToBcf::new(r.name(), &first_header)?),
            ),
        });
    }

    let out = open_output(spec, out_format, oo)?;
    let (records, header) = (
        if out_format == OutFormat::Bcf {
            Format::Bcf
        } else {
            Format::Vcf
        },
        Header::Vcf(first_header),
    );
    let wo = WriterOptions {
        write_header: true,
        command_line: None,
        ..wo.clone()
    };
    let mut writer = open_record_writer(out, records, out_format, &header, &wo).map_err(named)?;
    let mut boundary = Boundary::default();
    let mut batch = RecordBatch::default();
    let mut line = Vec::new();
    for (reader, conv) in readers.iter_mut().zip(convs.iter_mut()) {
        let name = reader.name().to_string();
        let binary = reader.format() == Format::Bcf;
        let mut guard = OrderGuard::new(&name, reader.format());
        let mut max = 1;
        let mut first = true;
        let mut last: Option<(ChromId, u64, String)> = None;
        loop {
            batch.clear();
            let more = reader.read_batch(&mut chroms, &mut batch, max)?;
            guard.check(&chroms, &batch)?;
            if let Some(m) = batch.recs.first().filter(|_| first) {
                first = false;
                let shown = shown_record(&chroms, &batch, 0, binary);
                let loc = if binary {
                    format!("record {}", m.loc)
                } else {
                    format!("line {}", m.loc)
                };
                boundary.check_first(&chroms, &name, m.chrom, m.key, &loc, &shown)?;
            }
            for i in 0..batch.len() {
                let raw = batch.raw(i);
                let loc = batch.recs[i].loc;
                let res = match conv {
                    Conv::Raw => writer.write_record(raw),
                    Conv::ToText(t) => {
                        let text = t.render(raw, loc)?;
                        writer.write_record(text)
                    }
                    Conv::ToBcf(e) => {
                        line.clear();
                        e.encode(raw, loc, &mut line)?;
                        writer.write_record(&line)
                    }
                    Conv::Recode(t, e) => {
                        let text = t.render(raw, loc)?.to_vec();
                        line.clear();
                        e.encode(&text, loc, &mut line)?;
                        writer.write_record(&line)
                    }
                };
                res.map_err(|source| IoError::Write {
                    path: dest.clone(),
                    source,
                })?;
            }
            if let Some(i) = batch.len().checked_sub(1) {
                let m = &batch.recs[i];
                last = Some((m.chrom, m.key, shown_record(&chroms, &batch, i, binary)));
            }
            max = next_batch_size(&batch, max);
            if !more {
                break;
            }
        }
        if let Some((chrom, key, shown)) = last {
            boundary.set_last(&name, chrom, key, shown);
        }
        if let Some(s) = reader.diagnostics().summary(&name) {
            log::warn!("{s}");
        }
    }
    writer.finish().map_err(named)
}

/// Record `i` of `batch` as an error shows it: the (cut) line, or
/// `CHROM POS` for BCF.
fn shown_record(chroms: &ChromTable, batch: &RecordBatch, i: usize, binary: bool) -> String {
    let m = &batch.recs[i];
    if binary {
        format!(
            "{} {}",
            String::from_utf8_lossy(chroms.name(m.chrom)),
            m.key + 1
        )
    } else {
        show(batch.raw(i))
    }
}

/// Where `--naive` writes: a temp file published atomically, stdout, or a
/// device or pipe written directly.
enum Sink {
    Atomic(AtomicFile),
    Stdout(io::BufWriter<io::Stdout>),
    Direct(io::BufWriter<File>),
}

impl Write for Sink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Sink::Atomic(w) => w.write(buf),
            Sink::Stdout(w) => w.write(buf),
            Sink::Direct(w) => w.write(buf),
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match self {
            Sink::Atomic(w) => w.flush(),
            Sink::Stdout(w) => w.flush(),
            Sink::Direct(w) => w.flush(),
        }
    }
}

/// The CHROM and POS-1 of a VCF line.
fn chrom_pos(line: &[u8]) -> Option<(&[u8], u64)> {
    let mut f = line.split(|&b| b == b'\t');
    let chrom = f.next().filter(|c| !c.is_empty())?;
    let pos: u64 = std::str::from_utf8(f.next()?).ok()?.parse().ok()?;
    Some((chrom, pos.checked_sub(1)?))
}

/// The first differing header line, for the message.
fn first_difference(a: &[u8], b: &[u8]) -> String {
    let (mut la, mut lb) = (a.split(|&c| c == b'\n'), b.split(|&c| c == b'\n'));
    loop {
        match (la.next(), lb.next()) {
            (Some(x), Some(y)) if x == y => {}
            (x, y) => {
                let text = |l: Option<&[u8]>| {
                    l.filter(|l| !l.is_empty())
                        .map_or_else(|| "(no line)".to_string(), show)
                };
                return format!("{} vs {}", text(x), text(y));
            }
        }
    }
}

fn naive(
    files: &[String],
    readers: Vec<Box<dyn RecordReader>>,
    mut chroms: ChromTable,
    spec: &OutputSpec,
) -> Result<(), IoError> {
    drop(readers);
    refuse_binary_to_tty(
        spec.path.is_none() && io::stdout().is_terminal(),
        OutFormat::Vcf,
        true,
    )?;
    let dest = spec
        .path
        .as_ref()
        .map_or_else(|| "<stdout>".to_string(), |p| p.display().to_string());
    let werr = |source| IoError::Write {
        path: dest.clone(),
        source,
    };
    let mut sink = match &spec.path {
        None => Sink::Stdout(io::BufWriter::new(io::stdout())),
        Some(p) if spec.direct => {
            let f = std::fs::OpenOptions::new()
                .write(true)
                .truncate(std::fs::metadata(p).is_ok_and(|m| m.is_file()))
                .open(p)
                .map_err(|source| IoError::Create {
                    path: p.display().to_string(),
                    source,
                })?;
            refuse_binary_to_tty(f.is_terminal(), OutFormat::Vcf, true)?;
            Sink::Direct(io::BufWriter::new(f))
        }
        Some(p) => Sink::Atomic(AtomicFile::create(p).map_err(|source| IoError::Create {
            path: p.display().to_string(),
            source,
        })?),
    };

    let mut header: Option<Vec<u8>> = None;
    let mut boundary = Boundary::default();
    let mut missing_newline: Option<&str> = None;
    for (idx, path) in files.iter().enumerate() {
        if let Some(prev) = missing_newline.take() {
            return Err(IoError::Detect {
                path: prev.to_string(),
                msg: format!(
                    "the last record of {prev} has no final newline, so --naive would join it to the next record of {path}; leave out --naive"
                ),
            });
        }
        let mut src = File::open(path).map_err(|source| IoError::Open {
            path: path.clone(),
            source,
        })?;
        let info = copy_blocks(
            &mut src,
            &mut sink,
            header.as_deref(),
            idx + 1 < files.len(),
        ).map_err(|e| match e {
            NaiveError::Write(source) => IoError::Write {
                path: dest.clone(),
                source,
            },
            NaiveError::Read(source) => IoError::Read {
                path: path.clone(),
                source,
            },
            NaiveError::HeaderDiffers(other) => usage(format!(
                "--naive needs byte-identical headers, but the header of {path} differs from that of {} ({}); make the headers identical, or leave out --naive to merge them",
                files[0],
                first_difference(header.as_deref().unwrap_or_default(), &other)
            )),
            e => IoError::Detect {
                path: path.clone(),
                msg: format!(
                    "{path}: {e}; --naive needs well-formed BGZF-compressed VCF files, so check the file or leave out --naive"
                ),
            },
        })?;
        if !info.eof_marker {
            log::warn!(
                "{path}: BGZF end-of-file marker missing; the file is probably truncated (its blocks are copied as they are)"
            );
        }
        order_check(&mut boundary, &mut chroms, path, &info)?;
        if !info.ends_with_newline && idx + 1 < files.len() {
            missing_newline = Some(path);
        }
        header.get_or_insert(info.header);
    }
    sink.write_all(&BGZF_EOF).map_err(werr)?;
    sink.flush().map_err(werr)?;
    match sink {
        Sink::Atomic(af) => af.commit().map_err(werr),
        _ => Ok(()),
    }
}

/// Check the first record of `info` against the last of the file before,
/// and remember its last record.
fn order_check(
    boundary: &mut Boundary,
    chroms: &mut ChromTable,
    path: &str,
    info: &NaiveFile,
) -> Result<(), IoError> {
    let Some(first) = &info.first_record else {
        return Ok(());
    };
    let bad = |what: &str, line: &[u8]| IoError::Detect {
        path: path.to_string(),
        msg: format!(
            "{path}: cannot read the chromosome and position of its {what} record ({}); check the file",
            show(line)
        ),
    };
    let mut resolve = |what: &str, line: &[u8]| -> Result<(ChromId, u64), IoError> {
        let (name, key) = chrom_pos(line).ok_or_else(|| bad(what, line))?;
        let id = chroms.intern(name).map_err(|e| IoError::Detect {
            path: path.to_string(),
            msg: format!("{path}: unknown chromosome {}", e.0),
        })?;
        Ok((id, key))
    };
    let (fc, fk) = resolve("first", first)?;
    let last = match &info.last_record {
        Some(l) => Some((resolve("last", l)?, l)),
        None => None,
    };
    // The first record is on the line after the header.
    let line = info.header.iter().filter(|&&b| b == b'\n').count() + 1;
    boundary.check_first(chroms, path, fc, fk, &format!("line {line}"), &show(first))?;
    if let Some(((lc, lk), l)) = last {
        boundary.set_last(path, lc, lk, show(l));
    }
    Ok(())
}
