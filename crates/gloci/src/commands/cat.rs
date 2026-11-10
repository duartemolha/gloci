//! `gloci cat`: join BAM files by copying their compressed blocks (D123).
use crate::Ctx;
use crate::cli::CatArgs;
use gloci_io::{
    OutFormat, OutputOptions, OutputSpec, UsageError, WriterOptions, open_raw_output,
    parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::cat_bams;
use std::io::IsTerminal;

/// The output spec: BAM to a file named for BAM (or for nothing it knows), or to stdout.
fn bam_spec(a: &CatArgs) -> Result<OutputSpec, UsageError> {
    match parse_output_spec(a.output.as_deref(), None) {
        Ok(s) if s.format.is_none() || s.format == Some(OutFormat::Bam) => Ok(s),
        Ok(s) => Err(UsageError(format!(
            "cat writes BAM, but {} names {}; use an output name ending in .bam",
            a.output
                .as_deref()
                .map_or_else(String::new, |p| p.display().to_string()),
            s.format.map_or("another format", OutFormat::name)
        ))),
        // An output name without a known extension (a device, say) is taken as BAM.
        Err(_) => parse_output_spec(a.output.as_deref(), Some("bam")),
    }
}

pub fn run(ctx: &Ctx, a: &CatArgs) -> anyhow::Result<()> {
    if a.write_index {
        return Err(UsageError(
            "--write-index is not available for cat yet: a BAM index for cat output arrives with \
gloci index for BAM. Until then, index the output with samtools index"
                .into(),
        )
        .into());
    }
    let spec = bam_spec(a)?;
    refuse_binary_to_tty(
        spec.path.is_none() && std::io::stdout().is_terminal(),
        OutFormat::Bam,
        true,
    )?;
    let oo = OutputOptions {
        threads: ctx.threads,
        level: 6,
        write_index: false,
    };
    let wo = WriterOptions {
        write_header: true,
        command_line: Some(ctx.command_line.clone()),
        version: ctx.version,
    };
    let dest = spec
        .path
        .as_ref()
        .map_or_else(|| "<stdout>".to_string(), |p| p.display().to_string());
    let out = open_raw_output(&spec, OutFormat::Bam, &oo)?;
    cat_bams(&a.files, out, &wo).map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}
