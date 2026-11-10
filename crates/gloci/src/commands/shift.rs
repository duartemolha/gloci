//! `gloci shift`: move each feature along its chromosome.
use super::slop::{Shared, run_shared};
use crate::Ctx;
use crate::cli::ShiftArgs;
use gloci_core::{Fraction, FractionError};
use gloci_io::UsageError;
use gloci_ops::shift::{Shift, ShiftBy};

fn parse_by(v: &str, pct: bool) -> Result<ShiftBy, UsageError> {
    if pct {
        let (negative, body) = match v.strip_prefix('-') {
            Some(b) => (true, b),
            None => (false, v.strip_prefix('+').unwrap_or(v)),
        };
        Fraction::parse_ratio(body)
            .map(|fraction| ShiftBy::Fraction { negative, fraction })
            .map_err(|e| match e {
                FractionError::RatioRange(_) => {
                    UsageError(format!("--by: fraction \"{v}\" is too large"))
                }
                e => UsageError(format!("--by: {e}")),
            })
    } else {
        v.strip_prefix('+')
            .unwrap_or(v)
            .parse::<i64>()
            .map(ShiftBy::Bases)
            .map_err(|_| {
                UsageError(format!(
                    "--by: \"{v}\" is not a whole number of bases; use --pct for a fraction"
                ))
            })
    }
}

pub fn run(ctx: &Ctx, a: &ShiftArgs) -> anyhow::Result<()> {
    let shift = Shift {
        by: parse_by(&a.by, a.pct)?,
        strand_aware: a.strand,
    };
    run_shared(
        ctx,
        &Shared {
            name: "shift",
            file: &a.file,
            genome: a.genome.genome.as_deref(),
            format: a.format.format.as_deref(),
            header: a.header.header,
            skip_invalid: a.invalid.skip_invalid,
            chrom_map: &a.chrom_map,
            out: &a.out,
        },
        |reader, chroms, genome, w| gloci_ops::shift::shift(reader, chroms, &shift, genome, w),
    )
}
