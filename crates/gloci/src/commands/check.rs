//! `gloci check`: validate files and print a one-screen summary of each.
use crate::Ctx;
use crate::cli::{CheckArgs, VcfSpan};
use gloci_io::{
    Accept, AlignOptions, Format, InvalidPolicy, IoError, ReaderOptions, UsageError, resolve_genome,
};
use gloci_ops::{CheckReport, QuickOptions, check_file, quick_check};

/// Returns true when no file has errors.
pub fn run(ctx: &Ctx, a: &CheckArgs) -> anyhow::Result<bool> {
    let format = match a.format.format.as_deref() {
        None => None,
        Some(s) => Some(s.parse::<Format>().map_err(UsageError)?),
    };
    let ro = ReaderOptions {
        accept: Accept {
            alignments: Some(AlignOptions {
                reference: a.reference.clone(),
                full_records: true,
                ..AlignOptions::default()
            }),
            convert_only: false,
        },
        chrom_map: a.chrom_map.load()?,
        threads: ctx.threads,
        format,
        vcf_span: match a.vcf_span.vcf_span {
            VcfSpan::Allele => gloci_core::VcfSpan::Allele,
            VcfSpan::Ref => gloci_core::VcfSpan::Ref,
        },
        invalid: InvalidPolicy::Skip,
    };
    // A bad -g is one error for the command, not one per file.
    let genome = a.genome.genome.as_deref().map(resolve_genome).transpose()?;
    use std::io::Write;
    let mut stdout = std::io::stdout().lock();
    let mut failed = 0;
    for (i, f) in a.files.iter().enumerate() {
        let report = if a.quick {
            quick_check(
                f,
                &ro,
                &QuickOptions {
                    unaligned: a.unaligned,
                },
            )
        } else {
            match check_file(f, genome.as_ref(), a.chrom_order.chrom_order, &ro) {
                Ok(r) => r,
                // A usage problem applies to the whole command.
                Err(e @ (IoError::Usage(_) | IoError::Gate { .. })) => return Err(e.into()),
                Err(e) => CheckReport::unreadable(f, e.to_string()),
            }
        };
        if a.list_failed {
            // Only the failing paths on stdout; what is wrong with each goes to stderr.
            if report.has_errors() {
                failed += 1;
                writeln!(stdout, "{f}")?;
                stdout.flush()?;
                for line in &report.errors.first {
                    if line.contains(f.as_str()) {
                        eprintln!("{line}");
                    } else {
                        eprintln!("{f}: {line}");
                    }
                }
            }
            continue;
        }
        if i > 0 {
            stdout.write_all(b"\n")?;
        }
        stdout.write_all(report.render().as_bytes())?;
        stdout.flush()?;
        if report.has_errors() {
            failed += 1;
        }
    }
    if a.list_failed {
        // The paths are the whole of stdout.
    } else if failed == 0 {
        writeln!(stdout, "result: OK")?;
    } else {
        writeln!(
            stdout,
            "result: FAILED ({failed} of {} files have errors)",
            a.files.len()
        )?;
    }
    stdout.flush()?;
    Ok(failed == 0)
}
