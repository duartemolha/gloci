//! Output specification (D17), format/level validation, and the atomic/BGZF output writer.

use crate::atomic::AtomicFile;
use crate::detect::Format;
use crate::error::{IoError, UsageError};
use crate::index::TsvColumns;
use crate::staged::StagedSet;
use noodles_bgzf::io::multithreaded_writer::Builder as BgzfBuilder;
use noodles_bgzf::io::writer::CompressionLevel;
use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter, IsTerminal, Stdout, StdoutLock, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutFormat {
    Bed,
    BedGraph,
    Vcf,
    Bcf,
    Gff3,
    Gtf,
    Tsv,
    /// BAM output (D100): BGZF always.
    Bam,
    Sam,
    /// FASTA output (`.fa`, `.fasta`, `.fna`; `.gz` is BGZF).
    Fasta,
}

impl OutFormat {
    /// The name users know the format by, for messages.
    pub fn name(self) -> &'static str {
        match self {
            OutFormat::Bed => "BED",
            OutFormat::BedGraph => "bedGraph",
            OutFormat::Vcf => "VCF",
            OutFormat::Bcf => "BCF",
            OutFormat::Gff3 => "GFF3",
            OutFormat::Gtf => "GTF",
            OutFormat::Tsv => "TSV",
            OutFormat::Bam => "BAM",
            OutFormat::Sam => "SAM",
            OutFormat::Fasta => "FASTA",
        }
    }

    /// The record format an index of this output describes (`None` for TSV, SAM and FASTA).
    pub fn index_format(self) -> Option<Format> {
        match self {
            OutFormat::Bed | OutFormat::BedGraph => Some(Format::Bed),
            OutFormat::Vcf => Some(Format::Vcf),
            OutFormat::Bcf => Some(Format::Bcf),
            OutFormat::Gff3 => Some(Format::Gff3),
            OutFormat::Gtf => Some(Format::Gtf),
            OutFormat::Bam => Some(Format::Bam),
            OutFormat::Tsv | OutFormat::Sam | OutFormat::Fasta => None,
        }
    }

    /// The output format that corresponds to an input format.
    pub fn family(f: Format) -> OutFormat {
        match f {
            Format::Bed => OutFormat::Bed,
            // WIG is input only; its records are bedGraph rows.
            Format::BedGraph | Format::Wig => OutFormat::BedGraph,
            Format::Vcf => OutFormat::Vcf,
            Format::Bcf => OutFormat::Bcf,
            Format::Gff3 => OutFormat::Gff3,
            Format::Gtf => OutFormat::Gtf,
            Format::Bam | Format::Cram => OutFormat::Bam,
            Format::Sam => OutFormat::Sam,
            Format::Fasta => OutFormat::Fasta,
            // PSL and RepeatMasker are read by `convert` only, which writes BED.
            Format::Psl | Format::Rmsk => OutFormat::Bed,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OutputSpec {
    /// `None` means stdout.
    pub path: Option<PathBuf>,
    /// `None` means "same as the input".
    pub format: Option<OutFormat>,
    pub bgzf: bool,
    /// `path` names an existing device, pipe or other non-regular file
    /// (`/dev/null`, `/dev/fd/N`, a FIFO): it is written to directly, with
    /// no temp file, rename or index.
    pub direct: bool,
}

#[derive(Debug, Clone, Copy)]
pub struct OutputOptions {
    pub threads: NonZeroUsize,
    pub level: u8,
    pub write_index: bool,
}

const ACCEPTED: &str = ".bed, .bed.gz, .bedgraph, .bedgraph.gz, .vcf, .vcf.gz, .bcf, \
.tsv, .tsv.gz, .gff, .gff.gz, .gff3, .gff3.gz, .gtf, .gtf.gz, .bam, .sam, .fa, .fa.gz, \
.fasta, .fasta.gz, .fna, .fna.gz";

/// `(suffix without dot, format, bgzf)`, matched case-insensitively.
fn lookup(s: &str) -> Option<(OutFormat, bool)> {
    use OutFormat::*;
    Some(match s {
        "bed" => (Bed, false),
        "bed.gz" => (Bed, true),
        "bedgraph" => (BedGraph, false),
        "bedgraph.gz" => (BedGraph, true),
        "vcf" => (Vcf, false),
        "vcf.gz" => (Vcf, true),
        "bcf" => (Bcf, true),
        "tsv" => (Tsv, false),
        "tsv.gz" => (Tsv, true),
        "gff" | "gff3" => (Gff3, false),
        "gff.gz" | "gff3.gz" => (Gff3, true),
        "gtf" => (Gtf, false),
        "gtf.gz" => (Gtf, true),
        "bam" => (Bam, true),
        "sam" => (Sam, false),
        "fa" | "fasta" | "fna" => (Fasta, false),
        "fa.gz" | "fasta.gz" | "fna.gz" => (Fasta, true),
        _ => return None,
    })
}

fn from_extension(path: &Path) -> Option<(OutFormat, bool)> {
    let name = path.file_name()?.to_string_lossy().to_ascii_lowercase();
    // try two-part suffix first (x.bed.gz), then one-part
    let mut dots = name.match_indices('.').map(|(i, _)| i).rev();
    let last = dots.next()?;
    if let Some(prev) = dots.next()
        && let Some(r) = lookup(&name[prev + 1..])
    {
        return Some(r);
    }
    lookup(&name[last + 1..])
}

/// `output` as a destination: `None` for stdout (`-`, `/dev/stdout`), and whether it is written
/// directly. An existing path that is not a regular file once symlinks are followed cannot take
/// a temp file and a rename; its format is chosen as for stdout. A descriptor path (`/dev/fd/N`,
/// as in `-o >(cmd)` or `3> file`) counts even when it leads to a regular file, because the
/// directory it sits in cannot hold a temp file.
fn output_target(output: Option<&Path>) -> (Option<&Path>, bool) {
    let path = output.filter(|p| p.as_os_str() != "-" && *p != Path::new("/dev/stdout"));
    let direct = path.is_some_and(|p| {
        p.starts_with("/dev/fd")
            || p.starts_with("/proc/self/fd")
            || p == Path::new("/dev/stderr")
            || std::fs::metadata(p).is_ok_and(|m| !m.is_file() && !m.is_dir())
    });
    (path, direct)
}

/// The spec of an output that holds plain text whatever its name (`view --count`): no format
/// is read from the extension, and nothing is compressed.
pub fn plain_text_output_spec(output: Option<&Path>) -> OutputSpec {
    let (path, direct) = output_target(output);
    OutputSpec {
        path: path.map(Path::to_path_buf),
        format: Some(OutFormat::Tsv),
        bgzf: false,
        direct,
    }
}

pub fn parse_output_spec(
    output: Option<&Path>,
    output_format: Option<&str>,
) -> Result<OutputSpec, UsageError> {
    // "-" and /dev/stdout mean stdout. (Opening /dev/stdout again would
    // restart at offset 0 of a file the shell redirected to.)
    let (path, direct) = output_target(output);
    let cram_error = || {
        UsageError(
            "gloci does not write CRAM; write BAM with -o FILE.bam and convert it with \
`samtools view -C -T ref.fa`"
                .into(),
        )
    };
    if path
        .and_then(|p| p.extension())
        .is_some_and(|e| e.eq_ignore_ascii_case("cram"))
        || output_format.is_some_and(|s| s.trim_start_matches('.').eq_ignore_ascii_case("cram"))
    {
        return Err(cram_error());
    }
    let from_o = match output_format {
        None => None,
        Some(s) => {
            let key = s.trim_start_matches('.').to_ascii_lowercase();
            Some(lookup(&key).ok_or_else(|| {
                UsageError(format!(
                    "unknown output format '{s}' for -O; use one of: {ACCEPTED}"
                ))
            })?)
        }
    };
    let from_ext = match path {
        None => None,
        Some(p) => match from_extension(p) {
            Some(r) => Some(r),
            // An unrecognised extension is acceptable when -O names the format.
            None if from_o.is_some() || direct => None,
            None => {
                return Err(UsageError(format!(
                    "cannot tell the output format from '{}'; use an extension from: {ACCEPTED}, or pass -O",
                    p.display()
                )));
            }
        },
    };
    let (format, bgzf) = match (from_ext, from_o) {
        (Some(a), Some(b)) if a != b => {
            return Err(UsageError(format!(
                "-o {} and -O {} disagree about the output format; drop one of them",
                path.map(|p| p.display().to_string()).unwrap_or_default(),
                output_format.unwrap_or_default()
            )));
        }
        (Some(a), _) | (None, Some(a)) => (Some(a.0), a.1),
        (None, None) => (None, false),
    };
    Ok(OutputSpec {
        path: path.map(Path::to_path_buf),
        format,
        bgzf,
        direct,
    })
}

pub fn check_level(n: u8) -> Result<u8, UsageError> {
    if n <= 12 {
        Ok(n)
    } else {
        Err(UsageError(format!(
            "--compression-level {n} is out of range; use a value from 0 to 12"
        )))
    }
}

pub fn refuse_binary_to_tty(is_tty: bool, f: OutFormat, bgzf: bool) -> Result<(), UsageError> {
    if is_tty && (matches!(f, OutFormat::Bcf | OutFormat::Bam) || bgzf) {
        Err(UsageError(
            "refusing to write binary (BCF or BGZF) output to a terminal; \
use -o FILE, or choose a text format with -O (e.g. -O vcf or -O bed)"
                .into(),
        ))
    } else {
        Ok(())
    }
}

type FileBuf = BufWriter<AtomicFile>;

enum Sink {
    StdoutPlain(BufWriter<StdoutLock<'static>>),
    StdoutBgzf(noodles_bgzf::io::MultithreadedWriter<BufWriter<Stdout>>),
    FilePlain(FileBuf),
    FileBgzf(noodles_bgzf::io::MultithreadedWriter<FileBuf>),
    /// A device or pipe, written to directly.
    DirectPlain(BufWriter<File>),
    DirectBgzf(noodles_bgzf::io::MultithreadedWriter<BufWriter<File>>),
}

pub struct Output {
    sink: Sink,
    dest: Option<PathBuf>,
    index: IndexPlan,
    /// FASTA output without `--write-index`: an old `.fai`/`.gzi` beside the
    /// destination describes other data and is removed when the file is replaced.
    drop_fasta_index: bool,
}

/// What `--write-index` builds at `finish`.
enum IndexPlan {
    None,
    /// The record format and thread count to index with.
    Format(Format, NonZeroUsize),
    /// TSV output: the columns are not known until the TSV writer sets them
    /// (D43, from A's format); finishing without them is a bug.
    TsvPending(NonZeroUsize),
    Tsv(TsvColumns, NonZeroUsize),
    /// FASTA output: `.fai`, and `.gzi` when the output is BGZF (D108).
    Fasta(NonZeroUsize),
}

const BUF: usize = 1 << 20;

pub fn open_output(
    spec: &OutputSpec,
    f: OutFormat,
    opts: &OutputOptions,
) -> Result<Output, IoError> {
    open_output_as(spec, f, opts, false)
}

/// An output for `f` (BGZF, so a binary format) that does not compress what is written to
/// it: the caller writes finished BGZF blocks (`gloci cat`, D123), and ends the stream
/// with [`BGZF_EOF`](crate::BGZF_EOF) itself. `--write-index` indexes the finished file.
pub fn open_raw_output(
    spec: &OutputSpec,
    f: OutFormat,
    opts: &OutputOptions,
) -> Result<Output, IoError> {
    open_output_as(spec, f, opts, true)
}

fn open_output_as(
    spec: &OutputSpec,
    f: OutFormat,
    opts: &OutputOptions,
    raw: bool,
) -> Result<Output, IoError> {
    // A raw output holds BGZF blocks it did not make: it is binary, but not compressed here.
    let bgzf = !raw && (spec.bgzf || matches!(f, OutFormat::Bcf | OutFormat::Bam));
    let level = CompressionLevel::try_from(check_level(opts.level)?)
        .map_err(|e| UsageError(e.to_string()))?;
    let mut index = IndexPlan::None;
    if opts.write_index {
        if spec.path.is_none() {
            return Err(UsageError(
                "--write-index needs -o FILE; an index cannot be written for stdout".into(),
            )
            .into());
        }
        if spec.direct {
            return Err(UsageError(format!(
                "--write-index needs a regular output file, but {} is a device or pipe; \
write to a file such as -o out.bed.gz and index that",
                spec.path.as_deref().unwrap_or(Path::new("")).display()
            ))
            .into());
        }
        if f == OutFormat::Sam {
            return Err(UsageError(
                "--write-index needs BAM output, as SAM text has no index; write BAM with: \
gloci view FILE -o x.bam --write-index"
                    .into(),
            )
            .into());
        }
        // A raw output holds finished BGZF blocks, which is what an index needs.
        if !bgzf && !raw && f != OutFormat::Fasta {
            return Err(UsageError(
                "--write-index needs BGZF output; use a .gz output name (or .bcf, or .bam)".into(),
            )
            .into());
        }
        index = match f.index_format() {
            _ if f == OutFormat::Fasta => IndexPlan::Fasta(opts.threads),
            Some(format) => IndexPlan::Format(format, opts.threads),
            // The columns depend on A's format (D43); `open_tsv_writer` decides.
            None => IndexPlan::TsvPending(opts.threads),
        };
    }
    let builder = || {
        BgzfBuilder::default()
            .set_compression_level(level)
            .set_worker_count(opts.threads)
    };
    let sink = match &spec.path {
        None if bgzf => Sink::StdoutBgzf(
            builder().build_from_writer(BufWriter::with_capacity(BUF, io::stdout())),
        ),
        None => Sink::StdoutPlain(BufWriter::with_capacity(BUF, io::stdout().lock())),
        Some(p) if spec.direct => {
            let create = |source| IoError::Create {
                path: p.display().to_string(),
                source,
            };
            let truncate = std::fs::metadata(p).is_ok_and(|m| m.is_file());
            let file = OpenOptions::new()
                .write(true)
                .truncate(truncate)
                .open(p)
                .map_err(create)?;
            refuse_binary_to_tty(file.is_terminal(), f, bgzf || raw)?;
            let w = BufWriter::with_capacity(BUF, file);
            if bgzf {
                Sink::DirectBgzf(builder().build_from_writer(w))
            } else {
                Sink::DirectPlain(w)
            }
        }
        Some(p) => {
            let af = AtomicFile::create(p).map_err(|source| IoError::Create {
                path: p.display().to_string(),
                source,
            })?;
            let w = BufWriter::with_capacity(BUF, af);
            if bgzf {
                Sink::FileBgzf(builder().build_from_writer(w))
            } else {
                Sink::FilePlain(w)
            }
        }
    };
    Ok(Output {
        sink,
        dest: spec.path.clone(),
        drop_fasta_index: f == OutFormat::Fasta && !opts.write_index,
        index,
    })
}

impl Output {
    /// Index this TSV output by `cols` at `finish` (D43). Only meaningful
    /// when `--write-index` was given to `open_output`.
    pub fn set_tsv_index(&mut self, cols: TsvColumns) {
        if let IndexPlan::TsvPending(t) | IndexPlan::Tsv(_, t) = self.index {
            self.index = IndexPlan::Tsv(cols, t);
        }
    }

    /// Whether `--write-index` was asked for.
    pub(crate) fn wants_index(&self) -> bool {
        !matches!(self.index, IndexPlan::None)
    }

    fn dest_name(&self) -> String {
        self.dest
            .as_ref()
            .map_or_else(|| "<stdout>".into(), |p| p.display().to_string())
    }

    /// Finishes the stream (BGZF EOF block), builds the index from the temp
    /// file when `--write-index` is set, removes stale `PATH.tbi`/`PATH.csi`
    /// (D19), and publishes the file, then the index, atomically. If indexing
    /// fails, nothing is published.
    pub fn finish(self) -> Result<(), IoError> {
        let name = self.dest_name();
        let werr = |source| IoError::Write {
            path: name.clone(),
            source,
        };
        let atomic: Option<AtomicFile> = match self.sink {
            Sink::StdoutPlain(mut w) => {
                w.flush().map_err(werr)?;
                None
            }
            Sink::StdoutBgzf(mut w) => {
                let mut inner = w.finish().map_err(werr)?;
                inner.flush().map_err(werr)?;
                None
            }
            Sink::DirectPlain(mut w) => {
                w.flush().map_err(werr)?;
                None
            }
            Sink::DirectBgzf(mut w) => {
                let mut inner = w.finish().map_err(werr)?;
                inner.flush().map_err(werr)?;
                None
            }
            Sink::FilePlain(w) => Some(w.into_inner().map_err(|e| werr(e.into_error()))?),
            Sink::FileBgzf(mut w) => {
                let inner = w.finish().map_err(werr)?;
                Some(inner.into_inner().map_err(|e| werr(e.into_error()))?)
            }
        };
        if let (Some(af), Some(dest)) = (atomic, self.dest.as_ref()) {
            // Publish order: build (and so validate) the index from the temp
            // file, remove the old PATH.tbi/PATH.csi, rename the data, then
            // rename the index. Interrupted at any point, PATH is left with
            // no index, never beside an index of other data (D19).
            let built = match &self.index {
                IndexPlan::None => None,
                IndexPlan::Format(format, threads) => Some(crate::index::build_index(
                    af.temp_path(),
                    dest,
                    *format,
                    *threads,
                )?),
                IndexPlan::Tsv(cols, threads) => Some(crate::index::build_tsv_index(
                    af.temp_path(),
                    dest,
                    cols,
                    *threads,
                )?),
                IndexPlan::Fasta(_) => None,
                IndexPlan::TsvPending(_) => {
                    return Err(IoError::Detect {
                        path: name.clone(),
                        msg:
                            "internal error: TSV output asked for an index but no columns were set"
                                .into(),
                    });
                }
            };
            let mut staged = Vec::new();
            if let Some((kind, bytes)) = built {
                let path = crate::index::index_path(dest, kind);
                staged.push((crate::index::stage_index(&path, &bytes)?, path));
            }
            if let IndexPlan::Fasta(threads) = &self.index {
                let (fai, gzi) = crate::fasta::fasta_index_bytes(af.temp_path(), *threads)?;
                let fai_dest = crate::fasta::fai_path(dest);
                staged.push((crate::index::stage_index(&fai_dest, &fai)?, fai_dest));
                if let Some(gzi) = gzi {
                    let gzi_dest = crate::fasta::gzi_path(dest);
                    staged.push((crate::index::stage_index(&gzi_dest, &gzi)?, gzi_dest));
                }
            }
            crate::index::remove_indexes(dest)?;
            if self.drop_fasta_index {
                for stale in [crate::fasta::fai_path(dest), crate::fasta::gzi_path(dest)] {
                    match std::fs::remove_file(&stale) {
                        Ok(()) => {}
                        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
                        Err(source) => {
                            return Err(IoError::Write {
                                path: stale.display().to_string(),
                                source,
                            });
                        }
                    }
                }
            }
            af.commit().map_err(werr)?;
            for (idx, path) in staged {
                idx.commit().map_err(|source| IoError::Write {
                    path: path.display().to_string(),
                    source,
                })?;
            }
        }
        Ok(())
    }
}

impl Output {
    /// Finish the stream like [`Output::finish`], but hand a regular file to
    /// `set` instead of renaming it into place, so that several outputs can be
    /// published together by [`StagedSet::publish`]. Outputs that are not
    /// regular files, and outputs with an index, are finished as usual.
    pub fn finish_staged(self, set: &mut StagedSet) -> Result<(), IoError> {
        if self.wants_index() || !matches!(self.sink, Sink::FilePlain(_) | Sink::FileBgzf(_)) {
            return self.finish();
        }
        let name = self.dest_name();
        let werr = |source| IoError::Write {
            path: name.clone(),
            source,
        };
        let af = match self.sink {
            Sink::FilePlain(w) => w.into_inner().map_err(|e| werr(e.into_error()))?,
            Sink::FileBgzf(mut w) => {
                let inner = w.finish().map_err(werr)?;
                inner.into_inner().map_err(|e| werr(e.into_error()))?
            }
            _ => return Err(IoError::Internal("finish_staged on a stream".into())),
        };
        if let Some(dest) = self.dest {
            set.push(af, dest);
        }
        Ok(())
    }
}

impl Write for Output {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match &mut self.sink {
            Sink::StdoutPlain(w) => w.write(buf),
            Sink::StdoutBgzf(w) => w.write(buf),
            Sink::FilePlain(w) => w.write(buf),
            Sink::FileBgzf(w) => w.write(buf),
            Sink::DirectPlain(w) => w.write(buf),
            Sink::DirectBgzf(w) => w.write(buf),
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match &mut self.sink {
            Sink::StdoutPlain(w) => w.flush(),
            Sink::StdoutBgzf(w) => w.flush(),
            Sink::FilePlain(w) => w.flush(),
            Sink::FileBgzf(w) => w.flush(),
            Sink::DirectPlain(w) => w.flush(),
            Sink::DirectBgzf(w) => w.flush(),
        }
    }
}
