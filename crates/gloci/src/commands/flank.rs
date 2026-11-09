//! `gloci flank`: the strips beside each feature.
use super::slop::{Shared, run_pad};
use crate::Ctx;
use crate::cli::FlankArgs;

pub fn run(ctx: &Ctx, a: &FlankArgs) -> anyhow::Result<()> {
    run_pad(
        ctx,
        &Shared {
            name: "flank",
            file: &a.file,
            genome: a.genome.genome.as_deref(),
            format: a.format.format.as_deref(),
            header: a.header.header,
            skip_invalid: a.invalid.skip_invalid,
            chrom_map: &a.chrom_map,
            out: &a.out,
        },
        &a.pad,
        gloci_ops::transform::flank,
    )
}
