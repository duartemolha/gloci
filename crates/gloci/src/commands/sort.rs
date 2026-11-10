//! `gloci sort`: external sort, `-m` merge of sorted inputs, `--check`.
use crate::Ctx;
use crate::cli::SortArgs;
use gloci_core::VcfSpan;
use gloci_io::{
    Format, Header, InputSpool, IoError, OrderGuard, OutFormat, OutputOptions, ReaderOptions,
    RecordBatch, RecordReader, UsageError, WriterOptions, check_family, chrom_table_for,
    is_regular_file, next_batch_size, open_output, open_reader, open_reader_from,
    open_record_writer, parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::{
    Dedup, MIN_MAX_MEM, SortConfig, Sorter, Tie, check_definitions_pair, check_samples_pair,
    merge_fan_in, merge_sorted_batches,
};
use std::fs::File;
use std::io::IsTerminal;
use std::path::{Path, PathBuf};

fn usage(msg: impl Into<String>) -> anyhow::Error {
    UsageError(msg.into()).into()
}

pub fn run(ctx: &Ctx, a: &SortArgs) -> anyhow::Result<()> {
    if a.max_mem < MIN_MAX_MEM {
        return Err(usage(format!(
            "--max-mem {} bytes is too small; use at least 16M",
            a.max_mem
        )));
    }
    if a.merge_sorted && (a.unique || a.duplicates) {
        return Err(usage(
            "-m cannot be combined with --unique or --duplicates; \
sort the inputs without -m to remove or find duplicates",
        ));
    }
    if a.check && (a.out.output.is_some() || a.out.write_index) {
        return Err(usage(
            "--check writes nothing, so -o and --write-index do not apply",
        ));
    }
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;

    let files = if a.files.is_empty() {
        vec!["-".to_string()]
    } else {
        a.files.clone()
    };
    let format = match a.format.format.as_deref() {
        None => None,
        Some(s) => Some(s.parse::<Format>().map_err(UsageError)?),
    };
    let opts = ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: a.chrom_map.load()?,
        threads: ctx.threads,
        format,
        vcf_span: VcfSpan::Ref,
        invalid: if a.invalid.skip_invalid {
            gloci_io::InvalidPolicy::Skip
        } else {
            gloci_io::InvalidPolicy::Fail
        },
    };
    // Every input is opened exactly once. When the inputs fit in one merge
    // round, the header-pass readers are kept for the record pass. Otherwise
    // regular files are reopened per round and any input that cannot be
    // reopened (stdin, a pipe, a FIFO) is first spooled once to --tmp-dir.
    let mut sources = Sources::open(&files, &opts, merge_fan_in(a.max_mem), &ctx.tmp_dir)?;
    let (format, header, mut chroms) = if files.len() == 1 {
        let r0 = sources.peek(0)?;
        let format = r0.format();
        let header = r0.header().clone();
        let chroms = chrom_table_for(a.genome.genome.as_deref(), &[r0], a.chrom_order.chrom_order)?;
        (format, header, chroms)
    } else {
        let r0 = sources.peek(0)?;
        let format = r0.format();
        let header = r0.header().clone();
        let mut chroms = None;
        if a.genome.genome.is_some() {
            chroms = Some(chrom_table_for(
                a.genome.genome.as_deref(),
                &[],
                a.chrom_order.chrom_order,
            )?);
        } else if let Header::Vcf(vh) = r0.header()
            && !vh.contigs.is_empty()
        {
            chroms = Some(chrom_table_for(None, &[r0], a.chrom_order.chrom_order)?);
        }
        sources.release(0);

        for i in 1..files.len() {
            let r = sources.peek(i)?;
            if OutFormat::family(r.format()) != OutFormat::family(format) {
                return Err(usage(format!(
                    "{} is {} but {} is {}; sort needs inputs of one format: sort each format \
separately",
                    files[0],
                    format.name(),
                    r.name(),
                    r.format().name()
                )));
            }
            if let (Header::Vcf(first), Header::Vcf(other)) = (&header, r.header()) {
                check_samples_pair(&files[0], first, r.name(), other)?;
                if !a.check {
                    check_definitions_pair(&files[0], first, r.name(), other)?;
                }
                if chroms.is_none() && !other.contigs.is_empty() {
                    chroms = Some(chrom_table_for(None, &[r], a.chrom_order.chrom_order)?);
                }
            }
            sources.release(i);
        }
        let chroms = match chroms {
            Some(c) => c,
            None => chrom_table_for(None, &[], a.chrom_order.chrom_order)?,
        };
        (format, header, chroms)
    };

    if a.check {
        let mut batch = RecordBatch::default();
        for i in 0..files.len() {
            let mut r = sources.take(i)?;
            let mut guard = OrderGuard::new(r.name(), r.format());
            let mut max = 1;
            loop {
                batch.clear();
                let more = r.read_batch(&mut chroms, &mut batch, max)?;
                guard.check(&chroms, &batch)?;
                max = next_batch_size(&batch, max);
                if !more {
                    break;
                }
            }
            if let Some(s) = r.diagnostics().summary(r.name()) {
                log::warn!("{s}");
            }
        }
        return Ok(());
    }

    let out_format = spec.format.unwrap_or_else(|| OutFormat::family(format));
    check_family(format, out_format)?;
    refuse_binary_to_tty(
        spec.path.is_none() && std::io::stdout().is_terminal(),
        out_format,
        spec.bgzf,
    )?;
    let dest = spec
        .path
        .as_ref()
        .map_or_else(|| "<stdout>".to_string(), |p| p.display().to_string());
    let (tie, dedup) = (
        if a.full_key {
            Tie::FullLine
        } else {
            Tie::InputOrder
        },
        if a.unique {
            Dedup::Unique
        } else if a.duplicates {
            Dedup::Duplicates
        } else {
            Dedup::All
        },
    );
    // Creating the spill directory first makes a bad --tmp-dir fail before
    // any input is read.
    let sorter = if a.merge_sorted {
        None
    } else {
        Some(Sorter::new(SortConfig {
            max_mem: a.max_mem,
            tmp_dir: ctx.tmp_dir.clone(),
            threads: ctx.threads,
            tie,
            dedup,
        })?)
    };
    let output = open_output(
        &spec,
        out_format,
        &OutputOptions {
            threads: ctx.threads,
            level: a.out.compression_level,
            write_index: a.out.write_index,
        },
    )?;
    let wopts = WriterOptions {
        // GFF headers carry ##gff-version; BED headers are opt-in (D14).
        write_header: a.header.header || matches!(format, Format::Gff3 | Format::Gtf),
        command_line: Some(ctx.command_line.clone()),
        version: ctx.version,
    };

    match sorter {
        Some(mut sorter) => {
            let mut batch = RecordBatch::default();
            for i in 0..files.len() {
                let mut r = sources.take(i)?;
                let mut max = 1;
                loop {
                    batch.clear();
                    let more = r.read_batch(&mut chroms, &mut batch, max)?;
                    sorter.push(&chroms, &batch)?;
                    max = next_batch_size(&batch, max);
                    if !more {
                        break;
                    }
                }
                if let Some(s) = r.diagnostics().summary(r.name()) {
                    log::warn!("{s}");
                }
            }
            let mut w = open_record_writer(output, format, out_format, &header, &wopts)
                .map_err(|e| e.with_output_name(&dest))?;
            sorter
                .finish(&chroms, &mut |raw| w.write_record(raw))
                .map_err(|e| e.with_output_name(&dest))?;
            w.finish().map_err(|e| e.with_output_name(&dest))?;
        }
        None => {
            let mut w = open_record_writer(output, format, out_format, &header, &wopts)
                .map_err(|e| e.with_output_name(&dest))?;
            merge_sorted_batches(
                files.len(),
                &mut chroms,
                tie,
                a.max_mem,
                &ctx.tmp_dir,
                |start, count| (start..start + count).map(|i| sources.take(i)).collect(),
                &mut |raw| w.write_record(raw),
            )
            .map_err(|e| e.with_output_name(&dest))?;
            w.finish().map_err(|e| e.with_output_name(&dest))?;
        }
    }
    Ok(())
}

/// The inputs of one run. Each is read from its source exactly once: a
/// regular file is reopened for each pass, so the number of open readers and
/// their BGZF thread pools stays bounded however many files there are. An
/// input that cannot be reopened (stdin, a pipe, a FIFO, `/dev/fd/N`) keeps
/// its header-pass reader when every input fits in one merge round
/// (`fan_in`), and is otherwise spooled to --tmp-dir before anything else is
/// read.
struct Sources {
    paths: Vec<String>,
    opts: ReaderOptions,
    /// Per input: whether it is a regular file, which is reopened for each
    /// pass rather than held open.
    regular: Vec<bool>,
    /// Whether the inputs that cannot be reopened fit in one merge round, so
    /// their readers are kept for the record pass instead of being spooled.
    keep: bool,
    readers: Vec<Option<Box<dyn RecordReader>>>,
    spooled: Vec<Option<PathBuf>>,
    // Declared last: the spool files outlive every reader of them.
    _spool: Option<InputSpool>,
}

impl Sources {
    fn open(
        paths: &[String],
        opts: &ReaderOptions,
        fan_in: usize,
        tmp_dir: &Path,
    ) -> Result<Sources, IoError> {
        let keep = paths.len() <= fan_in;
        let regular: Vec<bool> = paths.iter().map(|p| is_regular_file(p)).collect();
        let mut spool = None;
        let mut spooled = vec![None; paths.len()];
        if !keep {
            for (i, p) in paths.iter().enumerate() {
                if !regular[i] {
                    let s = match &mut spool {
                        Some(s) => s,
                        None => spool.insert(InputSpool::create(tmp_dir)?),
                    };
                    spooled[i] = Some(s.spool(p)?);
                }
            }
        }
        Ok(Sources {
            paths: paths.to_vec(),
            opts: opts.clone(),
            regular,
            keep,
            readers: (0..paths.len()).map(|_| None).collect(),
            spooled,
            _spool: spool,
        })
    }

    fn reopen(&self, i: usize) -> Result<Box<dyn RecordReader>, IoError> {
        match &self.spooled[i] {
            None => open_reader(&self.paths[i], &self.opts),
            Some(copy) => {
                let name = if self.paths[i] == "-" {
                    "stdin"
                } else {
                    &self.paths[i]
                };
                let f = File::open(copy).map_err(|source| IoError::Read {
                    path: name.to_string(),
                    source,
                })?;
                open_reader_from(name, Box::new(f), &self.opts)
            }
        }
    }

    /// The reader of input `i` for the header pass.
    fn peek(&mut self, i: usize) -> Result<&dyn RecordReader, IoError> {
        if self.readers[i].is_none() {
            self.readers[i] = Some(self.reopen(i)?);
        }
        Ok(self.readers[i].as_deref().expect("just opened"))
    }

    /// End the header pass for input `i`: a regular file's reader is closed
    /// (the file is reopened for the record pass, so at most one reader of
    /// regular files is open at a time); any other reader is kept for the
    /// record pass, or closed when it was spooled.
    fn release(&mut self, i: usize) {
        if self.regular[i] || !self.keep {
            self.readers[i] = None;
        }
    }

    /// The reader of input `i` for the record pass.
    fn take(&mut self, i: usize) -> Result<Box<dyn RecordReader>, IoError> {
        match self.readers[i].take() {
            Some(r) => Ok(r),
            None => self.reopen(i),
        }
    }
}
