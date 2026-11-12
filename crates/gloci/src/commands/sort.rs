//! `gloci sort`: external sort, `-m` merge of sorted inputs, `--check`.
use crate::Ctx;
use crate::cli::SortArgs;
use gloci_core::VcfSpan;
use gloci_io::{
    AlignOptions, AlnSortOrder, AlnSource, Format, Header, InputSpool, IoError, Opened, OrderGuard,
    OutFormat, OutputOptions, ReaderOptions, RecordBatch, RecordReader, UsageError, WriterOptions,
    check_family, chrom_table_for, is_regular_file, next_batch_size, open_aln_source, open_output,
    open_reader, open_reader_from, open_reader_or_alignments, open_record_writer,
    parse_output_spec, refuse_binary_to_tty,
};
use gloci_io::{open_aln_source_from, open_reader_or_alignments_from};
use gloci_ops::threads_per_reader;
use gloci_ops::{
    AlnSortOptions, Dedup, MIN_MAX_MEM, Order, SortConfig, Sorter, Tie, check_alignments,
    check_definitions_pair, check_samples_pair, merge_alignments_lazy, merge_fan_in,
    merge_sorted_batches, sort_alignments,
};
use std::fs::File;
use std::io::IsTerminal;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

/// The reader options of alignment input: every record, and the CRAM reference.
fn align_options(reference: Option<&Path>) -> AlignOptions {
    AlignOptions {
        full_records: true,
        reference: reference.map(Path::to_path_buf),
        ..AlignOptions::default()
    }
}

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
    // Each input opened for a merge round gets its share of the threads (see run_alignments).
    let probe_threads = threads_per_reader(ctx.threads, files.len().min(merge_fan_in(a.max_mem)));
    let mut sources = Sources::open(&files, &opts, merge_fan_in(a.max_mem), &ctx.tmp_dir)?;
    if let Some(src) = sources.probe_alignments(a.reference.as_deref(), probe_threads)? {
        // (With -m the other inputs may be CRAM: only a single input is known to be none.)
        if !a.merge_sorted || files.len() == 1 {
            gloci_io::warn_unused_reference(a.reference.as_deref(), src.format() == Format::Cram);
        }
        return run_alignments(ctx, a, &sources, &files, src, &spec);
    }
    gloci_io::warn_unused_reference(a.reference.as_deref(), false);
    if a.by_name || a.by_name_bytes {
        return Err(usage(format!(
            "{} is not an alignment file, and -n and -N sort alignments (BAM, SAM, CRAM) by read \
name; drop {}, or give a BAM, SAM or CRAM file",
            gloci_io::memio::display_name(&files[0]),
            if a.by_name { "-n" } else { "-N" }
        )));
    }
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
    let dest = spec.path.as_ref().map_or_else(
        || "<stdout>".to_string(),
        |p| gloci_io::memio::display_name(&p.to_string_lossy()),
    );
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
            order: Order::Interval(tie),
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

/// Sort (or `--check`) one alignment input (D121), or with `-m` merge sorted ones (D122).
fn run_alignments(
    ctx: &Ctx,
    a: &SortArgs,
    sources: &Sources,
    files: &[String],
    src: Box<dyn AlnSource>,
    spec: &gloci_io::OutputSpec,
) -> anyhow::Result<()> {
    if files.len() > 1 && !a.merge_sorted {
        return Err(usage(format!(
            "sort takes one alignment file, but {} were given; merge sorted BAM files with \
gloci sort -m, or join files already in order with gloci cat",
            files.len()
        )));
    }
    for (given, flag) in [
        (a.genome.genome.is_some(), "-g"),
        (a.full_key, "--full-key"),
        (a.unique, "--unique"),
        (a.duplicates, "--duplicates"),
    ] {
        if given {
            return Err(usage(format!(
                "{flag} does not apply to alignments, which are ordered by the header's \
references; drop {flag}"
            )));
        }
    }
    let order = if a.by_name {
        AlnSortOrder::Name
    } else if a.by_name_bytes {
        AlnSortOrder::NameBytes
    } else {
        AlnSortOrder::Coordinate
    };
    // Header pass: every input's header, for the merged header and the checks. Input 0 is open
    // already. Inputs that fit in one merge round stay open for the record pass (a pipe cannot
    // be read twice); with more, each is closed and reopened (a regular file by its path, any
    // other input from its spool) when its round starts, so at most one round's inputs are
    // open, each with its share of the threads.
    let n = files.len();
    let fan = merge_fan_in(a.max_mem);
    let keep = n <= fan;
    let mut ropts = sources.opts.clone();
    ropts.accept.alignments = Some(align_options(a.reference.as_deref()));
    ropts.threads = threads_per_reader(ctx.threads, n.min(fan));
    let mut src0 = Some(src);
    let mut held: Vec<Option<Box<dyn AlnSource>>> = (0..n).map(|_| None).collect();
    let mut headers = Vec::with_capacity(n);
    let mut all_sam = true;
    for (i, slot) in held.iter_mut().enumerate() {
        let s = match src0.take() {
            Some(s) => s,
            None => sources.open_aln(i, &ropts)?,
        };
        headers.push((s.name().to_string(), s.header().clone()));
        all_sam &= s.layout() == Format::Sam;
        if keep {
            *slot = Some(s);
        }
    }
    let mut open = |i: usize| -> Result<Box<dyn AlnSource>, IoError> {
        match held[i].take() {
            Some(s) => Ok(s),
            None => sources.open_aln(i, &ropts),
        }
    };
    if a.check {
        for i in 0..n {
            if !check_alignments(open(i)?, order)? {
                // check_alignments already logged the error naming both records.
                return Err(crate::Reported.into());
            }
        }
        return Ok(());
    }
    if a.out.write_index && order != AlnSortOrder::Coordinate {
        return Err(usage(
            "--write-index needs coordinate order, as a BAM index finds reads by position; \
drop --write-index, or drop -n/-N",
        ));
    }
    // The inputs' format (D104): SAM when every input is SAM, else BAM (a CRAM, or a mix of
    // formats, writes BAM); -o and -O override.
    let out_format = spec.format.unwrap_or(if all_sam {
        OutFormat::Sam
    } else {
        OutFormat::Bam
    });
    if !matches!(out_format, OutFormat::Bam | OutFormat::Sam) {
        return Err(usage(format!(
            "{} output is not possible for sorted alignments; write BAM or SAM (-o x.bam, \
-o x.sam); CRAM input is written as BAM",
            out_format.name()
        )));
    }
    refuse_binary_to_tty(
        spec.path.is_none() && std::io::stdout().is_terminal(),
        out_format,
        spec.bgzf,
    )?;
    let dest = spec.path.as_ref().map_or_else(
        || "<stdout>".to_string(),
        |p| gloci_io::memio::display_name(&p.to_string_lossy()),
    );
    let output = open_output(
        spec,
        out_format,
        &OutputOptions {
            threads: ctx.threads,
            level: a.out.compression_level,
            write_index: a.out.write_index,
        },
    )?;
    let wopts = WriterOptions {
        write_header: true,
        command_line: Some(ctx.command_line.clone()),
        version: ctx.version,
    };
    let o = AlnSortOptions {
        order,
        max_mem: a.max_mem,
        tmp_dir: ctx.tmp_dir.clone(),
        threads: ctx.threads,
    };
    if a.merge_sorted {
        let n = merge_alignments_lazy(
            &headers,
            &mut open,
            order,
            fan,
            &ctx.tmp_dir,
            output,
            out_format,
            &wopts,
        )
        .map_err(|e| e.with_output_name(&dest))?;
        log::info!("merged {n} records from {} files", files.len());
        return Ok(());
    }
    let src = open(0)?;
    let stats = sort_alignments(src, &o, output, out_format, &wopts)
        .map_err(|e| e.with_output_name(&dest))?;
    log::info!(
        "sorted {} records, {} chunks spilled to disk",
        stats.records,
        stats.spills
    );
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
                    "stdin".to_string()
                } else {
                    gloci_io::memio::display_name(&self.paths[i])
                };
                let f = File::open(copy).map_err(|source| IoError::Read {
                    path: name.clone(),
                    source,
                })?;
                open_reader_from(&name, Box::new(f), &self.opts)
            }
        }
    }

    /// Open input 0 as every-record alignments when it is a BAM, SAM or CRAM file. Any other
    /// format is kept as the reader of the header pass, so the input is not opened twice.
    fn probe_alignments(
        &mut self,
        a_ref: Option<&Path>,
        threads: NonZeroUsize,
    ) -> Result<Option<Box<dyn AlnSource>>, IoError> {
        let mut o = self.opts.clone();
        o.accept.alignments = Some(align_options(a_ref));
        o.threads = threads;
        let opened = match &self.spooled[0] {
            None => open_reader_or_alignments(&self.paths[0], &o)?,
            Some(copy) => {
                let name = if self.paths[0] == "-" {
                    "stdin"
                } else {
                    &self.paths[0]
                };
                let f = File::open(copy).map_err(|source| IoError::Read {
                    path: name.to_string(),
                    source,
                })?;
                open_reader_or_alignments_from(name, Box::new(f), &o)?
            }
        };
        match opened {
            Opened::Alignments(src) => Ok(Some(src)),
            Opened::Records(r) => {
                self.readers[0] = Some(r);
                Ok(None)
            }
        }
    }

    /// Open input `i` as every-record alignments: a regular file by its path, any other input
    /// from its spooled copy (a pipe cannot be read twice).
    fn open_aln(&self, i: usize, o: &ReaderOptions) -> Result<Box<dyn AlnSource>, IoError> {
        match &self.spooled[i] {
            None => open_aln_source(&self.paths[i], o),
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
                open_aln_source_from(name, Box::new(f), o)
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
