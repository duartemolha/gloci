//! `gloci liftover`: move BED and bedGraph records to another assembly through a chain, by
//! CrossMap's rules (D174).
use crate::Ctx;
use crate::cli::LiftoverArgs;
use gloci_core::ChromTable;
use gloci_io::{
    ChainIndex, Format, InvalidPolicy, IoError, OutputOptions, OutputSpec, ReaderOptions,
    UsageError, WriterOptions, input_family_output_spec, open_reader, parse_output_spec,
};
use gloci_ops::liftover::{ChromStyle, LiftoverCounts, LiftoverOptions, liftover};
use gloci_ops::{ReportOutput, open_report};
use std::path::Path;

fn style(s: &str) -> ChromStyle {
    match s {
        "short" => ChromStyle::Short,
        "long" => ChromStyle::Long,
        "chain" => ChromStyle::Chain,
        _ => ChromStyle::AsInput,
    }
}

/// The same file, spelt two ways: by name, through a link (same device and inode), or by where
/// it is (or would be) on disk.
fn same_file(a: &Path, b: &Path) -> bool {
    let real = |p: &Path| {
        let name = p.file_name()?.to_owned();
        let dir = match p.parent() {
            Some(d) if !d.as_os_str().is_empty() => d,
            _ => Path::new("."),
        };
        Some(dir.canonicalize().ok()?.join(name))
    };
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        if let (Ok(x), Ok(y)) = (std::fs::metadata(a), std::fs::metadata(b))
            && x.dev() == y.dev()
            && x.ino() == y.ino()
        {
            return true;
        }
    }
    a == b || (real(a).is_some() && real(a) == real(b))
}

fn check_outputs(mapped: &OutputSpec, unmapped: &OutputSpec) -> Result<(), UsageError> {
    match (&mapped.path, &unmapped.path) {
        (None, None) => Err(UsageError(
            "--unmapped - and the output are both stdout: name a file with -o or --unmapped".into(),
        )),
        (Some(m), Some(u)) if !mapped.direct && !unmapped.direct && same_file(m, u) => {
            Err(UsageError(format!(
                "-o and --unmapped are both {}: the two outputs need different files",
                m.display()
            )))
        }
        _ => Ok(()),
    }
}

fn reasons(c: &LiftoverCounts) -> String {
    c.unmapped
        .iter()
        .map(|(why, n)| format!("{n} {}", why.describe()))
        .collect::<Vec<_>>()
        .join(", ")
}

/// Name the output a write error came from (the core labels them `output` and
/// `unmapped output`).
fn name_write(e: IoError, mapped: &str, unmapped: &str) -> IoError {
    match e {
        IoError::Write { path, source } if path == "unmapped output" => IoError::Write {
            path: unmapped.to_string(),
            source,
        },
        e => e.with_output_name(mapped),
    }
}

pub fn run(ctx: &Ctx, a: &LiftoverArgs) -> anyhow::Result<()> {
    let spec = parse_output_spec(a.out.output.as_deref(), a.out.output_format.as_deref())?;
    let uspec = a
        .unmapped
        .as_ref()
        .map(|p| input_family_output_spec(Some(p)));
    if let Some(u) = &uspec {
        check_outputs(&spec, u)?;
    }
    let chain = ChainIndex::load(&a.chain)?;
    let ro = ReaderOptions {
        accept: gloci_io::Accept::default(),
        chrom_map: a.chrom_map.load()?,
        threads: ctx.threads,
        format: match a.format.format.as_deref() {
            None => None,
            Some(f) => Some(f.parse::<Format>().map_err(UsageError)?),
        },
        vcf_span: gloci_core::VcfSpan::Allele,
        invalid: if a.invalid.skip_invalid {
            InvalidPolicy::Skip
        } else {
            InvalidPolicy::Fail
        },
    };
    let reader = open_reader(&a.file, &ro)?;
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.out.compression_level,
        write_index: a.out.write_index,
    };
    let wo = WriterOptions {
        write_header: a.header.header,
        command_line: Some(ctx.command_line.clone()),
        version: ctx.version,
    };
    let header = a.header.header;
    let mut out = open_report(
        reader.as_ref(),
        &ReportOutput::SameAsA,
        &spec,
        &oo,
        &wo,
        header,
    )?;
    let mut un = match &uspec {
        Some(u) => Some(
            open_report(reader.as_ref(), &ReportOutput::SameAsA, u, &oo, &wo, header)
                .map_err(|e| e.with_output_name(&super::dest_name(u)))?,
        ),
        None => None,
    };
    let (mapped_name, unmapped_name) = (
        super::dest_name(&spec),
        uspec.as_ref().map_or_else(String::new, super::dest_name),
    );
    let o = LiftoverOptions {
        chrom_style: style(&a.chrom_style),
        skip_invalid: a.invalid.skip_invalid,
        ignore_blocks: a.ignore_blocks,
    };
    // Records are written in input order, so the table only gives names ids.
    let mut chroms = ChromTable::first_appearance();
    let counts = liftover(
        reader,
        &mut chroms,
        &chain,
        &o,
        out.as_mut(),
        un.as_mut()
            .map(|w| w.as_mut() as &mut dyn gloci_io::RecordWriter),
    )
    .map_err(|e| name_write(e, &mapped_name, &unmapped_name))?;
    out.finish()
        .map_err(|e| name_write(e, &mapped_name, &unmapped_name))?;
    if let Some(u) = un {
        u.finish().map_err(|e| e.with_output_name(&unmapped_name))?;
    }

    let n = counts.unmapped_total();
    if n > 0 {
        if uspec.is_some() {
            log::info!(
                "{n} of {} records did not lift and are in {unmapped_name}: {}",
                counts.records,
                reasons(&counts)
            );
        } else {
            log::warn!(
                "{n} of {} records did not lift: {}; write them with --unmapped FILE",
                counts.records,
                reasons(&counts)
            );
        }
    }
    if counts.split + counts.shortened > 0 {
        log::warn!(
            "{} records were split into several records at chain gaps and {} lost bases no chain \
block covers (CrossMap does both without a word)",
            counts.split,
            counts.shortened
        );
    }
    if counts.invalid > 0 {
        let n = counts.invalid;
        log::warn!(
            "{n} record{} skipped as invalid",
            if n == 1 { " was" } else { "s were" }
        );
    }
    log::info!(
        "{} chains, {} blocks from {}",
        gloci_core::thousands(chain.chains() as u64),
        gloci_core::thousands(chain.blocks() as u64),
        a.chain
    );
    Ok(())
}
