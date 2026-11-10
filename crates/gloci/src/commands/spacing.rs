//! `gloci spacing`: each record with its distance to the previous one.
use super::cluster::{AppendArgs, run_appending};
use crate::Ctx;
use crate::cli::SpacingArgs;
use gloci_ops::spacing::spacing;

pub fn run(ctx: &Ctx, a: &SpacingArgs) -> anyhow::Result<()> {
    run_appending(
        ctx,
        &AppendArgs {
            file: &a.file,
            genome: &a.genome,
            chrom_order: &a.chrom_order,
            vcf_span: &a.vcf_span,
            chrom_map: &a.chrom_map,
            format: &a.format,
            header: &a.header,
            invalid: &a.invalid,
            sort: &a.sort,
            out: &a.out,
        },
        "spacing",
        |reader, chroms, w| spacing(reader, chroms, w),
    )
}
