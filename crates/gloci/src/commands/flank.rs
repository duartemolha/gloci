//! `gloci flank`: the strips beside each feature.
use super::slop::{Shared, run_shared};
use crate::Ctx;
use crate::cli::FlankArgs;

pub fn run(ctx: &Ctx, a: &FlankArgs) -> anyhow::Result<()> {
    run_shared(
        ctx,
        &Shared {
            name: "flank",
            file: &a.file,
            pad: &a.pad,
            genome: a.genome.genome.as_deref(),
            format: a.format.format.as_deref(),
            header: a.header.header,
            skip_invalid: a.invalid.skip_invalid,
            chrom_map: &a.chrom_map,
            out: &a.out,
        },
        gloci_ops::transform::flank,
    )
}
