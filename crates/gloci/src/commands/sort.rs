//! `gloci sort`: external sort, `-m` merge of sorted inputs, `--check`.
use crate::Ctx;
use crate::cli::SortArgs;
use gloci_core::VcfSpan;
use gloci_io::{
    Format, Header, OrderGuard, OutFormat, OutputOptions, ReaderOptions, RecordBatch, RecordReader,
    UsageError, WriterOptions, check_family, chrom_table_for, next_batch_size, open_output,
    open_reader, open_record_writer, parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::{
    Dedup, MIN_MAX_MEM, SortConfig, Sorter, Tie, check_definitions_pair, check_samples_pair,
    merge_sorted, merge_sorted_batches,
};
use std::io::IsTerminal;

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
        threads: ctx.threads,
        format,
        vcf_span: VcfSpan::Ref,
        invalid: if a.invalid.skip_invalid {
            gloci_io::InvalidPolicy::Skip
        } else {
            gloci_io::InvalidPolicy::Fail
        },
    };
    let (format, header, mut chroms, single_reader) = if files.len() == 1 {
        let r0 = open_reader(&files[0], &opts)?;
        let format = r0.format();
        let header = r0.header().clone();
        let chroms = chrom_table_for(
            a.genome.genome.as_deref(),
            &[r0.as_ref()],
            a.chrom_order.chrom_order,
        )?;
        (format, header, chroms, Some(r0))
    } else {
        let r0 = open_reader(&files[0], &opts)?;
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
            chroms = Some(chrom_table_for(
                None,
                &[r0.as_ref()],
                a.chrom_order.chrom_order,
            )?);
        }
        drop(r0);

        for f in &files[1..] {
            let r = open_reader(f, &opts)?;
            if OutFormat::family(r.format()) != OutFormat::family(format) {
                return Err(usage(format!(
                    "{} is {:?} but {} is {:?}; sort needs inputs of one format",
                    files[0],
                    format,
                    r.name(),
                    r.format()
                )));
            }
            if let (Header::Vcf(first), Header::Vcf(other)) = (&header, r.header()) {
                check_samples_pair(&files[0], first, r.name(), other)?;
                if !a.check {
                    check_definitions_pair(&files[0], first, r.name(), other)?;
                }
                if chroms.is_none() && !other.contigs.is_empty() {
                    chroms = Some(chrom_table_for(
                        None,
                        &[r.as_ref()],
                        a.chrom_order.chrom_order,
                    )?);
                }
            }
        }
        let chroms = match chroms {
            Some(c) => c,
            None => chrom_table_for(None, &[], a.chrom_order.chrom_order)?,
        };
        (format, header, chroms, None)
    };

    if a.check {
        let mut batch = RecordBatch::default();
        if let Some(mut r) = single_reader {
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
        } else {
            for f in &files {
                let mut r = open_reader(f, &opts)?;
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
            if let Some(mut r) = single_reader {
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
            } else {
                for f in &files {
                    let mut r = open_reader(f, &opts)?;
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
            if let Some(r) = single_reader {
                let mut readers: [Box<dyn RecordReader>; 1] = [r];
                merge_sorted(&mut readers, &mut chroms, tie, &mut |raw| {
                    w.write_record(raw)
                })
                .map_err(|e| e.with_output_name(&dest))?;
                if let Some(s) = readers[0].diagnostics().summary(readers[0].name()) {
                    log::warn!("{s}");
                }
            } else {
                merge_sorted_batches(
                    files.len(),
                    &mut chroms,
                    tie,
                    a.max_mem,
                    &ctx.tmp_dir,
                    |start, count| {
                        let mut batch = Vec::with_capacity(count);
                        for f in &files[start..start + count] {
                            batch.push(open_reader(f, &opts)?);
                        }
                        Ok(batch)
                    },
                    &mut |raw| w.write_record(raw),
                )
                .map_err(|e| e.with_output_name(&dest))?;
            }
            w.finish().map_err(|e| e.with_output_name(&dest))?;
        }
    }
    Ok(())
}
