//! `gloci head`: the header and the first records of any input as text (D134).
use crate::Ctx;
use crate::cli::HeadArgs;
use gloci_io::{Accept, AlignOptions, InvalidPolicy, IoError, ReaderOptions};
use gloci_ops::{HeadOptions, head};
use std::io::{BufWriter, Write};

pub fn run(ctx: &Ctx, a: &HeadArgs) -> anyhow::Result<()> {
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
    let o = HeadOptions {
        headers: a.headers,
        records: a.records,
    };
    let mut out = BufWriter::new(std::io::stdout().lock());
    // What was read before a failure stays written, as `head` does with a cut file.
    let result = head(&a.file, &o, &ro, &mut out).map_err(|e| e.with_output_name("stdout"));
    let flushed = out.flush().map_err(|e| IoError::Write {
        path: "stdout".into(),
        source: e,
    });
    result?;
    flushed?;
    Ok(())
}
