//! `gloci view`: stream a file unchanged or re-compressed, convert VCF and
//! BCF to each other, and write only the header or only the records.

use gloci_core::ChromOrder;
use gloci_io::{
    HeaderMode, IoError, OutFormat, OutputOptions, OutputSpec, ReaderOptions, RecordBatch,
    UsageError, WriterOptions, check_family, chrom_table_for, convert_variants, next_batch_size,
    open_headerless_writer, open_output, open_reader, open_record_writer, refuse_binary_to_tty,
};
use std::io::IsTerminal;

#[derive(Debug, Clone, Default)]
pub struct ViewOptions {
    /// Write the header and no records.
    pub header_only: bool,
    /// Write the records and no header.
    pub no_header: bool,
}

fn usage(msg: impl Into<String>) -> IoError {
    UsageError(msg.into()).into()
}

/// View `path` into the output `out_spec`. With no regions the file is
/// streamed in file order; sortedness is not required.
pub fn run_view(
    path: &str,
    opts: &ViewOptions,
    ro: &ReaderOptions,
    out_spec: &OutputSpec,
    oo: &OutputOptions,
    wo: &WriterOptions,
) -> Result<(), IoError> {
    if opts.header_only && opts.no_header {
        return Err(usage("--header-only and --no-header cannot be combined"));
    }
    let mut reader = open_reader(path, ro)?;
    let format = reader.format();
    let family = OutFormat::family(format);
    let out_format = out_spec.format.unwrap_or(family);
    let is_variant = |f: OutFormat| matches!(f, OutFormat::Vcf | OutFormat::Bcf);
    if out_format != family && !(is_variant(family) && is_variant(out_format)) {
        return Err(usage(
            "use gloci convert for format conversion (coming in a later release); \
view only converts between VCF and BCF",
        ));
    }
    if opts.no_header && out_format == OutFormat::Bcf {
        return Err(usage(
            "--no-header cannot be used with BCF output; BCF records need their header",
        ));
    }
    refuse_binary_to_tty(
        out_spec.path.is_none() && std::io::stdout().is_terminal(),
        out_format,
        out_spec.bgzf,
    )?;
    let dest = out_spec
        .path
        .as_ref()
        .map_or_else(|| "<stdout>".to_string(), |p| p.display().to_string());

    let mut wo = wo.clone();
    // View keeps whatever header the file has unless told otherwise.
    wo.write_header = !opts.no_header;
    let out = open_output(out_spec, out_format, oo)?;

    if opts.header_only {
        // No records are read, so nothing is checked either.
        let mode = HeaderMode::Only;
        if out_format != family {
            let mut chroms = chrom_table_for(None, &[], ChromOrder::default())?;
            convert_variants(
                reader.as_mut(),
                &mut chroms,
                out,
                out_format,
                &wo,
                mode,
                None,
            )
            .map_err(|e| e.with_output_name(&dest))?;
        } else {
            let header = reader.header().clone();
            open_record_writer(out, format, out_format, &header, &wo)
                .map_err(|e| e.with_output_name(&dest))?
                .finish()
                .map_err(|e| e.with_output_name(&dest))?;
        }
        return Ok(());
    }

    let refs = [reader.as_ref()];
    let mut chroms = chrom_table_for(None, &refs, ChromOrder::default())?;
    if out_format != family {
        let mode = if opts.no_header {
            HeaderMode::None
        } else {
            HeaderMode::Full
        };
        convert_variants(
            reader.as_mut(),
            &mut chroms,
            out,
            out_format,
            &wo,
            mode,
            None,
        )
        .map_err(|e| e.with_output_name(&dest))?;
    } else {
        check_family(format, out_format)?;
        let header = reader.header().clone();
        let mut w = if opts.no_header {
            open_headerless_writer(out, format, out_format)
        } else {
            open_record_writer(out, format, out_format, &header, &wo)
        }
        .map_err(|e| e.with_output_name(&dest))?;
        let mut batch = RecordBatch::default();
        let mut max = 1;
        loop {
            batch.clear();
            let more = reader.read_batch(&mut chroms, &mut batch, max)?;
            for i in 0..batch.len() {
                w.write_record(batch.raw(i))
                    .map_err(|source| IoError::Write {
                        path: dest.clone(),
                        source,
                    })?;
            }
            max = next_batch_size(&batch, max);
            if !more {
                break;
            }
        }
        w.finish().map_err(|e| e.with_output_name(&dest))?;
    }
    if let Some(s) = reader.diagnostics().summary(reader.name()) {
        log::warn!("{s}");
    }
    Ok(())
}
