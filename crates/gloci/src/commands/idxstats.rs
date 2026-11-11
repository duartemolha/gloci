//! `gloci idxstats`: reads per reference of a BAM, CRAM or SAM (D127).
use crate::Ctx;
use crate::cli::IdxstatsArgs;
use gloci_io::{Accept, AlignOptions, InvalidPolicy, ReaderOptions};
use gloci_ops::idxstats;
use std::io::BufWriter;

pub fn run(ctx: &Ctx, a: &IdxstatsArgs) -> anyhow::Result<()> {
    let ro = ReaderOptions {
        accept: Accept {
            alignments: Some(AlignOptions {
                reference: a.reference.clone(),
                ..AlignOptions::default()
            }),
            convert_only: false,
        },
        chrom_map: None,
        threads: ctx.threads,
        format: None,
        vcf_span: gloci_core::VcfSpan::Allele,
        invalid: InvalidPolicy::Fail,
    };
    // Counted whole before anything is written, so a failing scan prints nothing.
    let mut text = Vec::new();
    idxstats(&a.file, &ro, &mut text)?;
    let mut out = BufWriter::new(std::io::stdout().lock());
    std::io::Write::write_all(&mut out, &text).map_err(|e| gloci_io::IoError::Write {
        path: "stdout".into(),
        source: e,
    })?;
    std::io::Write::flush(&mut out).map_err(|e| gloci_io::IoError::Write {
        path: "stdout".into(),
        source: e,
    })?;
    Ok(())
}
