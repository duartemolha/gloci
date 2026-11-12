//! `gloci reheader`: give a BAM or CRAM file a new header (D130).
use crate::Ctx;
use crate::cli::ReheaderArgs;
use gloci_io::{
    OutFormat, OutputOptions, OutputSpec, UsageError, WriterOptions, open_output, open_raw_output,
    parse_output_spec, refuse_binary_to_tty,
};
use gloci_ops::reheader::reheader;
use std::io::{IsTerminal, Read};
use std::path::Path;

/// Does `path` start with the CRAM magic? Anything unreadable is left to `reheader` to report.
fn is_cram(path: &Path) -> bool {
    let mut magic = [0u8; 4];
    std::fs::File::open(path)
        .and_then(|mut f| f.read_exact(&mut magic))
        .is_ok()
        && &magic == b"CRAM"
}

/// The output spec for a BAM input: BAM to a file named for BAM (or for nothing it knows), or
/// to stdout.
fn bam_spec(a: &ReheaderArgs) -> Result<OutputSpec, UsageError> {
    match parse_output_spec(a.output.as_deref(), None) {
        Ok(s) if s.format.is_none() || s.format == Some(OutFormat::Bam) => Ok(s),
        Ok(s) => Err(UsageError(format!(
            "{} is BAM, so the output is BAM, but {} names {}; use an output name ending in .bam",
            a.input.display(),
            a.output
                .as_deref()
                .map_or_else(String::new, |p| p.display().to_string()),
            s.format.map_or("another format", OutFormat::name)
        ))),
        Err(_) => parse_output_spec(a.output.as_deref(), Some("bam")),
    }
}

/// The output spec for a CRAM input: a plain file or stdout, in whatever name `-o` has unless
/// it names another format.
fn cram_spec(a: &ReheaderArgs) -> Result<OutputSpec, UsageError> {
    let path = a
        .output
        .as_deref()
        .filter(|p| p.as_os_str() != "-" && *p != Path::new("/dev/stdout"));
    if let Some(p) = path
        && let Some(ext) = p.extension().and_then(|e| e.to_str())
        && matches!(
            ext.to_ascii_lowercase().as_str(),
            "bam" | "sam" | "bcf" | "vcf" | "bed" | "gz"
        )
    {
        return Err(UsageError(format!(
            "{} is CRAM, so the output is CRAM, but {} names another format; use an output name \
ending in .cram",
            a.input.display(),
            p.display()
        )));
    }
    let direct = path.is_some_and(|p| {
        p.starts_with("/dev/fd")
            || p.starts_with("/proc/self/fd")
            || p == Path::new("/dev/stderr")
            || std::fs::metadata(p).is_ok_and(|m| !m.is_file() && !m.is_dir())
    });
    Ok(OutputSpec {
        path: path.map(Path::to_path_buf),
        format: None,
        bgzf: false,
        direct,
    })
}

pub fn run(ctx: &Ctx, a: &ReheaderArgs) -> anyhow::Result<()> {
    let cram = is_cram(&a.input);
    let spec = if cram { cram_spec(a)? } else { bam_spec(a)? };
    if spec.path.is_none() && std::io::stdout().is_terminal() {
        refuse_binary_to_tty(true, OutFormat::Bam, true)?;
    }
    let oo = OutputOptions {
        threads: ctx.threads,
        level: a.compression_level,
        write_index: false,
    };
    let wo = WriterOptions {
        write_header: true,
        command_line: Some(ctx.command_line.clone()),
        version: ctx.version,
    };
    let dest = super::dest_name(&spec);
    let out = if cram {
        open_output(&spec, OutFormat::Tsv, &oo)?
    } else {
        open_raw_output(&spec, OutFormat::Bam, &oo)?
    };
    reheader(&a.header, &a.input, out, &wo, a.compression_level)
        .map_err(|e| e.with_output_name(&dest))?;
    Ok(())
}
