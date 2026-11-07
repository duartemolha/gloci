//! `gloci check`: validate files and print a one-screen summary of each.
use crate::Ctx;
use crate::cli::{CheckArgs, VcfSpan};
use gloci_io::{Format, InvalidPolicy, ReaderOptions, UsageError};
use gloci_ops::check_file;

/// Returns true when no file has errors.
pub fn run(ctx: &Ctx, a: &CheckArgs) -> anyhow::Result<bool> {
    let format = match a.format.format.as_deref() {
        None => None,
        Some(s) => Some(s.parse::<Format>().map_err(UsageError)?),
    };
    let ro = ReaderOptions {
        threads: ctx.threads,
        format,
        vcf_span: match a.vcf_span.vcf_span {
            VcfSpan::Allele => gloci_core::VcfSpan::Allele,
            VcfSpan::Ref => gloci_core::VcfSpan::Ref,
        },
        invalid: InvalidPolicy::Skip,
    };
    let mut out = String::new();
    let mut failed = 0;
    for (i, f) in a.files.iter().enumerate() {
        let report = check_file(
            f,
            a.genome.genome.as_deref(),
            a.chrom_order.chrom_order,
            &ro,
        )?;
        if i > 0 {
            out.push('\n');
        }
        out.push_str(&report.render());
        if report.has_errors() {
            failed += 1;
        }
    }
    if failed == 0 {
        out.push_str("result: OK\n");
    } else {
        out.push_str(&format!(
            "result: FAILED ({failed} of {} files have errors)\n",
            a.files.len()
        ));
    }
    use std::io::Write;
    std::io::stdout().lock().write_all(out.as_bytes())?;
    Ok(failed == 0)
}
