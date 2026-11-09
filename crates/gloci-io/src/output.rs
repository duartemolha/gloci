//! Output specification (D17), format/level validation, and the atomic/BGZF output writer.

use crate::atomic::AtomicFile;
use crate::detect::Format;
use crate::error::{IoError, UsageError};
use crate::index::TsvColumns;
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
}

impl OutFormat {
    /// The record format an index of this output describes (`None` for TSV).
    pub fn index_format(self) -> Option<Format> {
        match self {
            OutFormat::Bed | OutFormat::BedGraph => Some(Format::Bed),
            OutFormat::Vcf => Some(Format::Vcf),
            OutFormat::Bcf => Some(Format::Bcf),
            OutFormat::Gff3 => Some(Format::Gff3),
            OutFormat::Gtf => Some(Format::Gtf),
            OutFormat::Tsv => None,
        }
    }

    /// The output format that corresponds to an input format.
    pub fn family(f: Format) -> OutFormat {
        match f {
            Format::Bed => OutFormat::Bed,
            Format::BedGraph => OutFormat::BedGraph,
            Format::Vcf => OutFormat::Vcf,
            Format::Bcf => OutFormat::Bcf,
            Format::Gff3 => OutFormat::Gff3,
            Format::Gtf => OutFormat::Gtf,
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
.tsv, .tsv.gz, .gff, .gff.gz, .gff3, .gff3.gz, .gtf, .gtf.gz";

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

pub fn parse_output_spec(
    output: Option<&Path>,
    output_format: Option<&str>,
) -> Result<OutputSpec, UsageError> {
    // "-" and /dev/stdout mean stdout. (Opening /dev/stdout again would
    // restart at offset 0 of a file the shell redirected to.)
    let path = output.filter(|p| p.as_os_str() != "-" && *p != Path::new("/dev/stdout"));
    // An existing path that is not a regular file once symlinks are followed
    // cannot take a temp file and a rename; its format is chosen as for stdout.
    // A descriptor path (`/dev/fd/N`, as in `-o >(cmd)` or `3> file`) counts
    // even when it leads to a regular file, because the directory it sits in
    // cannot hold a temp file.
    let direct = path.is_some_and(|p| {
        p.starts_with("/dev/fd")
            || p.starts_with("/proc/self/fd")
            || p == Path::new("/dev/stderr")
            || std::fs::metadata(p).is_ok_and(|m| !m.is_file() && !m.is_dir())
    });
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
    if is_tty && (f == OutFormat::Bcf || bgzf) {
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
}

const BUF: usize = 1 << 20;

pub fn open_output(
    spec: &OutputSpec,
    f: OutFormat,
    opts: &OutputOptions,
) -> Result<Output, IoError> {
    let bgzf = spec.bgzf || f == OutFormat::Bcf;
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
        if !bgzf {
            return Err(UsageError(
                "--write-index needs BGZF output; use a .gz output name (or .bcf)".into(),
            )
            .into());
        }
        index = match f.index_format() {
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
            refuse_binary_to_tty(file.is_terminal(), f, bgzf)?;
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
                IndexPlan::TsvPending(_) => {
                    return Err(IoError::Detect {
                        path: name.clone(),
                        msg:
                            "internal error: TSV output asked for an index but no columns were set"
                                .into(),
                    });
                }
            };
            let staged = match built {
                Some((kind, bytes)) => {
                    let path = crate::index::index_path(dest, kind);
                    Some((crate::index::stage_index(&path, &bytes)?, path))
                }
                None => None,
            };
            crate::index::remove_indexes(dest)?;
            af.commit().map_err(werr)?;
            if let Some((idx, path)) = staged {
                idx.commit().map_err(|source| IoError::Write {
                    path: path.display().to_string(),
                    source,
                })?;
            }
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
