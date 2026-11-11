//! `gloci samples`: the samples named by the read groups of alignment files (D135).
use crate::Ctx;
use crate::cli::SamplesArgs;
use gloci_io::{IoError, OutFormat, OutputOptions, UsageError, open_output, parse_output_spec};
use gloci_ops::{SamplesOptions, samples};

pub fn run(ctx: &Ctx, a: &SamplesArgs) -> anyhow::Result<()> {
    let tag: [u8; 2] = a.tag.as_bytes().try_into().map_err(|_| {
        UsageError(format!(
            "--tag must be two characters, as the tags of a header line are; got {:?}",
            a.tag
        ))
    })?;
    let o = SamplesOptions {
        tag,
        show_index: a.show_index,
        match_reference: a.match_reference.clone(),
        column_names: a.column_names.column_names,
    };
    let spec = parse_output_spec(a.output.as_deref(), Some("tsv"))?;
    let oo = OutputOptions {
        threads: ctx.threads,
        level: 6,
        write_index: false,
    };
    let dest = super::dest_name(&spec);
    let mut out = open_output(&spec, OutFormat::Tsv, &oo)?;
    samples(&a.files, &o, &mut out).map_err(|e: IoError| e.with_output_name(&dest))?;
    out.finish()?;
    Ok(())
}
