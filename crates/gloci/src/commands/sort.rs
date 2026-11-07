//! `gloci sort`: external sort, `-m` merge of sorted inputs, `--check`.
use crate::Ctx;
use crate::cli::SortArgs;
use gloci_core::VcfSpan;
use gloci_io::{
    Format, Header, IoError, OrderGuard, OutFormat, OutputOptions, ReaderOptions, RecordBatch,
    RecordReader, UsageError, WriterOptions, check_family, chrom_table_for, open_output,
    open_reader, open_record_writer, parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::{Dedup, SortConfig, Sorter, Tie, merge_sorted};
use std::io::IsTerminal;

/// `--max-mem` below this is a usage error (D23).
const MIN_MAX_MEM: u64 = 16 << 20;
/// Bytes of records to read per batch, and the bounds on the record count.
const BATCH_BYTES: usize = 4 << 20;
const MAX_BATCH: usize = 16384;

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
    let opts = ReaderOptions {
        threads: ctx.threads,
        format: None,
        vcf_span: VcfSpan::Ref,
        invalid: if a.invalid.skip_invalid {
            gloci_io::InvalidPolicy::Skip
        } else {
            gloci_io::InvalidPolicy::Fail
        },
    };
    let mut readers = files
        .iter()
        .map(|f| open_reader(f, &opts))
        .collect::<Result<Vec<_>, _>>()?;
    let format = readers[0].format();
    for r in &readers[1..] {
        if OutFormat::family(r.format()) != OutFormat::family(format) {
            return Err(usage(format!(
                "{} is {:?} but {} is {:?}; sort needs inputs of one format",
                readers[0].name(),
                format,
                r.name(),
                r.format()
            )));
        }
    }
    check_samples(&readers)?;
    let refs: Vec<&dyn RecordReader> = readers.iter().map(|r| r.as_ref()).collect();
    let mut chroms = chrom_table_for(a.genome.genome.as_deref(), &refs, a.chrom_order.chrom_order)?;

    if a.check {
        let mut batch = RecordBatch::default();
        for r in &mut readers {
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
        }
        warn_skipped(&readers);
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
    let header = readers[0].header().clone();

    match sorter {
        Some(mut sorter) => {
            let mut batch = RecordBatch::default();
            for r in &mut readers {
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
            }
            let mut w = open_record_writer(output, format, out_format, &header, &wopts)
                .map_err(|e| super::name_output(e, &dest))?;
            sorter
                .finish(&chroms, &mut |raw| w.write_record(raw))
                .map_err(|e| super::name_output(e, &dest))?;
            w.finish()?;
        }
        None => {
            let mut w = open_record_writer(output, format, out_format, &header, &wopts)
                .map_err(|e| super::name_output(e, &dest))?;
            merge_sorted(&mut readers, &mut chroms, tie, &mut |raw| {
                w.write_record(raw)
            })
            .map_err(|e| super::name_output(e, &dest))?;
            w.finish()?;
        }
    }
    warn_skipped(&readers);
    Ok(())
}

/// Records to request next, so a batch is about `BATCH_BYTES` whatever the
/// record length (a VCF line can be megabytes).
fn next_batch_size(batch: &RecordBatch, current: usize) -> usize {
    if batch.is_empty() {
        return current;
    }
    let per_rec = batch.buf.len() / batch.len() + 1;
    (BATCH_BYTES / per_rec).clamp(1, MAX_BATCH)
}

/// Inputs that are VCF/BCF need identical sample columns (D29).
fn check_samples(readers: &[Box<dyn RecordReader>]) -> Result<(), IoError> {
    let Header::Vcf(first) = readers[0].header() else {
        return Ok(());
    };
    for r in &readers[1..] {
        if let Header::Vcf(h) = r.header()
            && let (Some(a), Some(b)) = (&first.bcf, &h.bcf)
            && a != b
        {
            let (name, other) = (readers[0].name(), r.name());
            return Err(IoError::Detect {
                path: other.to_string(),
                msg: format!(
                    "the BCF headers of {name} and {other} differ (contig or INFO/FILTER/FORMAT dictionaries); \
BCF records refer to them by index, so they cannot be sorted together as they are. \
Convert them to VCF first (gloci view FILE -O vcf.gz) or use inputs with identical headers"
                ),
            });
        }
        if let Header::Vcf(h) = r.header()
            && h.samples != first.samples
        {
            let list = |s: &[Vec<u8>]| {
                s.iter()
                    .map(|n| String::from_utf8_lossy(n).into_owned())
                    .collect::<Vec<_>>()
                    .join(",")
            };
            let (name, other) = (readers[0].name(), r.name());
            return Err(IoError::Detect {
                path: other.to_string(),
                msg: format!(
                    "{other} has samples [{}] but {name} has [{}]; sort needs identical sample columns in every input",
                    list(&h.samples),
                    list(&first.samples)
                ),
            });
        }
    }
    Ok(())
}

fn warn_skipped(readers: &[Box<dyn RecordReader>]) {
    for r in readers {
        if let Some(s) = r.diagnostics().summary(r.name()) {
            log::warn!("{s}");
        }
    }
}
