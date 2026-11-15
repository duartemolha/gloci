//! Sorting a set of readers into one sorted reader (`--sort`, D39).
//!
//! [`sort_readers`] is the one place that turns unsorted inputs into a sorted
//! record stream for `merge`, `complement` and the relation
//! commands. It runs the external sorter (same chromosome table and tie
//! rules as `gloci sort`) and writes the sorted records to a file in a
//! registered temp directory under `tmp_dir`, then opens that file as a
//! [`RecordReader`]. Memory stays within `max_mem`; the temp directory is
//! removed when the returned reader is dropped, and on SIGINT/SIGTERM.

use super::{Order, SortConfig, Sorter, Tie};
use crate::sort::Dedup;
use gloci_core::{ChromTable, Strand, VcfSpan};
use gloci_io::{
    AlnHeader, AlnSortOrder, Compression, Diagnostics, Format, Header, InvalidPolicy, IoError,
    OutFormat, OutputOptions, OutputSpec, ReaderOptions, RecordBatch, RecordReader, SORTED_INPUTS,
    TempGuard, UsageError, WriterOptions, create_temp_dir_in, next_batch_size, open_output,
    open_reader, open_record_writer,
};
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::sync::Arc;
use tempfile::TempDir;

/// Settings of `--sort`: the memory budget, where spill files go, and how
/// the sorted stream is read back.
#[derive(Debug, Clone)]
pub struct SortSettings {
    pub max_mem: u64,
    pub tmp_dir: PathBuf,
    /// Threads for the sorter and the sorted file's writer and reader.
    pub threads: NonZeroUsize,
    /// The VCF/BCF span the sorted stream is read with: the caller's own
    /// `--vcf-span`. (Sorting itself always keys VCF records by REF span.)
    pub vcf_span: VcfSpan,
}

/// The options to open inputs with when they will be passed to
/// [`sort_readers`] (or to a relation run with `RelateOptions::sort`):
/// the same as `user`, except that VCF/BCF records are read by their REF
/// span, the span `gloci sort` keys them by. Opening them with another
/// span would sort them in an order the sorted file's reader then cannot
/// vouch for, silently. The user's span goes in [`SortSettings::vcf_span`].
pub fn sort_input_options(user: &ReaderOptions) -> ReaderOptions {
    ReaderOptions {
        vcf_span: VcfSpan::Ref,
        ..user.clone()
    }
}

/// Sort `inputs` together into one sorted reader, as a one-element list
/// (the shape `merge` and `complement` take). The inputs must have been
/// opened with [`sort_input_options`].
pub fn sort_together(
    inputs: Vec<Box<dyn RecordReader>>,
    chroms: &mut ChromTable,
    settings: &SortSettings,
) -> Result<Vec<Box<dyn RecordReader>>, IoError> {
    Ok(vec![sort_readers(inputs, chroms, settings)?])
}

/// `--max-mem` below this is a usage error (D23).
pub const MIN_MAX_MEM: u64 = 16 << 20;

/// A reader over a sorted temp file; removes the file's directory on drop.
struct SortedReader {
    // Fields drop in order: the open file first, then the directory, then
    // its registration.
    inner: Box<dyn RecordReader>,
    name: String,
    /// Data records returned so far.
    records: u64,
    _dir: TempDir,
    _guard: TempGuard,
}

impl RecordReader for SortedReader {
    fn name(&self) -> &str {
        &self.name
    }
    fn format(&self) -> Format {
        self.inner.format()
    }
    fn compression(&self) -> Compression {
        self.inner.compression()
    }
    fn header(&self) -> &Header {
        self.inner.header()
    }
    fn read_batch(
        &mut self,
        chroms: &mut ChromTable,
        batch: &mut RecordBatch,
        max: usize,
    ) -> Result<bool, IoError> {
        let from = batch.recs.len();
        let more = self.inner.read_batch(chroms, batch, max)?;
        // The file's own line numbers include the header lines written to it
        // and mean nothing to the user: number the data records instead.
        for m in &mut batch.recs[from..] {
            self.records += 1;
            m.loc = self.records;
        }
        Ok(more)
    }
    fn diagnostics(&self) -> &Diagnostics {
        self.inner.diagnostics()
    }
}

/// Sort all `inputs` together and return one reader of the sorted records.
///
/// The inputs must be of one format family, with identical VCF/BCF samples
/// and definitions (as for `gloci sort`). They must have been opened with
/// [`sort_input_options`] (REF span, the span `gloci sort` keys VCF records
/// by); the sorted file is read back with `settings.vcf_span`. `chroms` is the
/// table the caller will use for the sorted stream. Skipped-record warnings
/// of the inputs are logged here. The result is named `the sorted inputs (A, B)`
/// so errors do not show a temp path; `IoError::with_sorted_input_loc` turns
/// its line numbers into record numbers.
pub fn sort_readers(
    mut inputs: Vec<Box<dyn RecordReader>>,
    chroms: &mut ChromTable,
    settings: &SortSettings,
) -> Result<Box<dyn RecordReader>, IoError> {
    if settings.max_mem < MIN_MAX_MEM {
        return Err(UsageError(format!(
            "--max-mem {} bytes is too small; use at least 16M",
            settings.max_mem
        ))
        .into());
    }
    let first = inputs
        .first()
        .ok_or_else(|| IoError::Internal("sort_readers needs at least one input".into()))?;
    let format = first.format();
    for r in &inputs[1..] {
        if OutFormat::family(r.format()) != OutFormat::family(format) {
            return Err(UsageError(format!(
                "{} is {} but {} is {}; --sort needs inputs of one format: convert one of them, \
or run the command on each format separately",
                inputs[0].name(),
                format.name(),
                r.name(),
                r.format().name()
            ))
            .into());
        }
    }
    check_samples(&inputs)?;
    check_definitions(&inputs)?;

    let orig_names: Vec<String> = inputs.iter().map(|r| r.name().to_string()).collect();
    let mut chroms_per_input: Vec<Vec<Vec<u8>>> = vec![Vec::new(); inputs.len()];
    let mut marks = vec![gloci_io::naming::SeenChroms::default(); inputs.len()];
    let mut sorter = Sorter::new(SortConfig {
        max_mem: settings.max_mem,
        tmp_dir: settings.tmp_dir.clone(),
        threads: settings.threads,
        order: Order::Interval(Tie::InputOrder),
        dedup: Dedup::All,
    })?;
    let mut batch = RecordBatch::default();
    let aln = format.is_alignment();
    for (i, r) in inputs.iter_mut().enumerate() {
        let mut last_chrom: Option<gloci_core::ChromId> = None;
        let mut max = 1;
        loop {
            batch.clear();
            let more = r.read_batch(chroms, &mut batch, max)?;
            for m in &batch.recs {
                if last_chrom != Some(m.chrom) {
                    last_chrom = Some(m.chrom);
                    if marks[i].first_time(m.chrom) {
                        chroms_per_input[i].push(chroms.name(m.chrom).to_vec());
                    }
                }
            }
            if aln {
                // D121's coordinate order, the one `gloci sort` writes: (reference, position,
                // strand), then input order. The table is the header's (D124), so its rank is
                // the reference's.
                sorter.push_keyed(chroms, &batch, |m| {
                    (
                        (m.span.start << 1) | u64::from(m.strand == Strand::Minus),
                        0,
                    )
                })?;
            } else {
                sorter.push(chroms, &batch)?;
            }
            max = next_batch_size(&batch, max);
            if !more {
                break;
            }
        }
    }
    if let Some(w) = gloci_io::naming::nary_naming_warning(&orig_names, &chroms_per_input) {
        log::warn!("{w}");
    }
    for r in &inputs {
        if let Some(s) = r.diagnostics().summary(r.name()) {
            log::warn!("{s}");
        }
        if r.format().is_alignment() {
            r.diagnostics().log_alignment_counts(r.name());
        }
    }

    // Alignments: BAM (level 1) for BAM and CRAM, SAM for SAM, so the records come back in
    // the layout the caller already prepared to print (D124).
    let out_format = OutFormat::family(format);
    let whole = inputs[0].whole_records();
    let ext = match out_format {
        OutFormat::Vcf => "vcf",
        OutFormat::Bcf => "bcf",
        OutFormat::Bam => "bam",
        OutFormat::Sam => "sam",
        OutFormat::Gff3 | OutFormat::Gtf => "gff",
        _ => "bed",
    };
    let (dir, guard) =
        create_temp_dir_in(&settings.tmp_dir, "gloci-sorted-").map_err(|source| {
            IoError::Write {
                path: settings.tmp_dir.display().to_string(),
                source,
            }
        })?;
    let path = dir.path().join(format!("sorted.{ext}"));
    let spec = OutputSpec {
        path: Some(path.clone()),
        format: Some(out_format),
        bgzf: false,
        direct: false,
    };
    let output = open_output(
        &spec,
        out_format,
        &OutputOptions {
            threads: settings.threads,
            level: 1,
            write_index: false,
        },
    )?;
    let wopts = WriterOptions {
        write_header: format.is_alignment() || matches!(format, Format::Gff3 | Format::Gtf),
        command_line: None,
        version: env!("CARGO_PKG_VERSION"),
    };
    let header = match inputs[0].header() {
        Header::Alignment(h) => coordinate_header(h, chroms),
        h => h.clone(),
    };
    let names: Vec<&str> = inputs.iter().map(|r| r.name()).collect();
    let name = format!("{SORTED_INPUTS}{})", names.join(", "));
    let mut w = open_record_writer(output, format, out_format, &header, &wopts)?;
    let dest = path.display().to_string();
    let stats = sorter
        .finish(chroms, &mut |raw| w.write_record(raw))
        .map_err(|e| e.with_output_name(&dest))?;
    log::info!(
        "--sort: sorted {} records, {} chunks spilled to disk",
        stats.records,
        stats.spills
    );
    w.finish().map_err(|e| e.with_output_name(&dest))?;
    drop(inputs);

    let ro = ReaderOptions {
        accept: if format.is_alignment() {
            // The filters were applied when the records were first read. The records come
            // back as whole or as slim as they went in.
            gloci_io::Accept {
                alignments: Some(gloci_io::AlignOptions {
                    full_records: whole,
                    ..Default::default()
                }),
                ..Default::default()
            }
        } else {
            gloci_io::Accept::default()
        },
        chrom_map: None,
        threads: settings.threads,
        format: None,
        vcf_span: settings.vcf_span,
        invalid: InvalidPolicy::Fail,
    };
    let inner = open_reader(&path.to_string_lossy(), &ro)?;
    Ok(Box::new(SortedReader {
        inner,
        name,
        records: 0,
        _dir: dir,
        _guard: guard,
    }))
}

/// The header of alignments that `--sort` sorted with `chroms`: `h` with `SO:coordinate`, as
/// `gloci sort` leaves it, since the records come out in its order (D121). When the table
/// (from `-g`, say) orders the references unlike the header, the order is not one
/// `SO:coordinate` describes, and the header says `SO:unsorted`, whatever it said before.
pub(crate) fn coordinate_header(h: &AlnHeader, chroms: &ChromTable) -> Header {
    Header::Alignment(Arc::new(if follows_header(h, chroms) {
        h.with_sort_order(AlnSortOrder::Coordinate)
    } else {
        h.with_unsorted_order()
    }))
}

/// Whether `chroms` ranks the references of `h` that it knows in the header's order.
fn follows_header(h: &AlnHeader, chroms: &ChromTable) -> bool {
    let ranks = chroms.ranks();
    let mut seen = h
        .refs
        .iter()
        .filter_map(|(name, _)| chroms.get(name))
        .map(|id| ranks[id.0 as usize]);
    let Some(mut prev) = seen.next() else {
        return true;
    };
    seen.all(|r| {
        let ok = prev < r;
        prev = r;
        ok
    })
}

/// The command whose multi-file header rule is being checked, for messages.
#[derive(Debug, Clone, Copy)]
pub struct HeaderRuleCmd {
    /// `sort needs identical sample columns`
    pub verb: &'static str,
    /// `... cannot be sorted together as they are`
    pub together: &'static str,
    /// `give the files the same header definitions, or ...`
    pub fix: &'static str,
}

impl HeaderRuleCmd {
    pub const SORT: HeaderRuleCmd = HeaderRuleCmd {
        verb: "sort",
        together: "sorted together",
        fix: "sort them separately",
    };
    pub const CONCAT: HeaderRuleCmd = HeaderRuleCmd {
        verb: "concat",
        together: "concatenated",
        fix: "concat only files that share them",
    };
}

/// Inputs that are VCF/BCF need identical sample columns (D29).
pub fn check_samples_pair(
    first_name: &str,
    first: &gloci_io::VcfHeader,
    other_name: &str,
    other: &gloci_io::VcfHeader,
) -> Result<(), IoError> {
    check_samples_pair_for(HeaderRuleCmd::SORT, first_name, first, other_name, other)
}

/// [`check_samples_pair`] with the messages of `cmd`.
pub fn check_samples_pair_for(
    cmd: HeaderRuleCmd,
    first_name: &str,
    first: &gloci_io::VcfHeader,
    other_name: &str,
    other: &gloci_io::VcfHeader,
) -> Result<(), IoError> {
    if let (Some(a), Some(b)) = (&first.bcf, &other.bcf)
        && a != b
    {
        return Err(IoError::Detect {
            path: other_name.to_string(),
            msg: format!(
                "the BCF headers of {first_name} and {other_name} differ (contig or INFO/FILTER/FORMAT dictionaries); \
BCF records refer to them by index, so they cannot be {} as they are. \
Convert them to VCF first (gloci view FILE -O vcf.gz) or use inputs with identical headers",
                cmd.together
            ),
        });
    }
    if other.samples != first.samples {
        let list = |s: &[Vec<u8>]| {
            s.iter()
                .map(|n| String::from_utf8_lossy(n).into_owned())
                .collect::<Vec<_>>()
                .join(",")
        };
        return Err(IoError::Detect {
            path: other_name.to_string(),
            msg: format!(
                "{other_name} has samples [{}] but {first_name} has [{}]; {} needs identical sample columns in every input",
                list(&other.samples),
                list(&first.samples),
                cmd.verb
            ),
        });
    }
    Ok(())
}

/// Inputs that are VCF/BCF need identical sample columns (D29).
pub fn check_samples(readers: &[Box<dyn RecordReader>]) -> Result<(), IoError> {
    check_samples_for(HeaderRuleCmd::SORT, readers)
}

/// [`check_samples`] with the messages of `cmd`.
pub fn check_samples_for(
    cmd: HeaderRuleCmd,
    readers: &[Box<dyn RecordReader>],
) -> Result<(), IoError> {
    let Header::Vcf(first) = readers[0].header() else {
        return Ok(());
    };
    for r in &readers[1..] {
        if let Header::Vcf(h) = r.header() {
            check_samples_pair_for(cmd, readers[0].name(), first, r.name(), h)?;
        }
    }
    Ok(())
}

/// Inputs that are VCF/BCF must define the same INFO, FILTER, FORMAT, ALT and
/// contig lines: the output keeps the first file's header, so a definition
/// only another file has would be lost, leaving records that use an
/// undeclared key. Other meta lines may differ. `##FILTER=<ID=PASS>` is
/// implicitly present in every header.
pub fn check_definitions_pair(
    first_name: &str,
    first: &gloci_io::VcfHeader,
    other_name: &str,
    other: &gloci_io::VcfHeader,
) -> Result<(), IoError> {
    check_definitions_pair_for(HeaderRuleCmd::SORT, first_name, first, other_name, other)
}

/// [`check_definitions_pair`] with the messages of `cmd`.
pub fn check_definitions_pair_for(
    cmd: HeaderRuleCmd,
    first_name: &str,
    first: &gloci_io::VcfHeader,
    other_name: &str,
    other: &gloci_io::VcfHeader,
) -> Result<(), IoError> {
    let base = first.definitions();
    let implicit = |kind: &str, id: &[u8]| kind == "FILTER" && id == b"PASS";
    let show = |kind: &str, id: &[u8]| format!("##{kind}=<ID={}>", String::from_utf8_lossy(id));
    let fail = |path: &str, msg: String| IoError::Detect {
        path: path.to_string(),
        msg,
    };
    let fix = format!("give the files the same header definitions, or {}", cmd.fix);
    let theirs = other.definitions();
    for ((kind, id), line) in &theirs {
        match base.get(&(*kind, id.clone())) {
            None if implicit(kind, id) => {}
            None => {
                return Err(fail(
                    other_name,
                    format!(
                        "{other_name} defines {} but {first_name} does not; {fix}",
                        show(kind, id)
                    ),
                ));
            }
            Some(mine) if mine != line => {
                return Err(fail(
                    other_name,
                    format!(
                        "{first_name} and {other_name} define {} differently ({} vs {}); {fix}",
                        show(kind, id),
                        String::from_utf8_lossy(mine),
                        String::from_utf8_lossy(line)
                    ),
                ));
            }
            Some(_) => {}
        }
    }
    for (kind, id) in base.keys() {
        if !theirs.contains_key(&(*kind, id.clone())) && !implicit(kind, id) {
            return Err(fail(
                first_name,
                format!(
                    "{first_name} defines {} but {other_name} does not; {fix}",
                    show(kind, id)
                ),
            ));
        }
    }
    Ok(())
}

/// Inputs that are VCF/BCF must define the same INFO, FILTER, FORMAT, ALT and
/// contig lines: the output keeps the first file's header, so a definition
/// only another file has would be lost, leaving records that use an
/// undeclared key. Other meta lines may differ. `##FILTER=<ID=PASS>` is
/// implicitly present in every header.
pub fn check_definitions(readers: &[Box<dyn RecordReader>]) -> Result<(), IoError> {
    check_definitions_for(HeaderRuleCmd::SORT, readers)
}

/// [`check_definitions`] with the messages of `cmd`.
pub fn check_definitions_for(
    cmd: HeaderRuleCmd,
    readers: &[Box<dyn RecordReader>],
) -> Result<(), IoError> {
    let Header::Vcf(first) = readers[0].header() else {
        return Ok(());
    };
    for r in &readers[1..] {
        if let Header::Vcf(h) = r.header() {
            check_definitions_pair_for(cmd, readers[0].name(), first, r.name(), h)?;
        }
    }
    Ok(())
}
