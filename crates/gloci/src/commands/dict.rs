//! `gloci dict`: a SAM sequence dictionary of a FASTA (D132).
use crate::Ctx;
use crate::cli::DictArgs;
use gloci_io::{FastaStream, IoError, OutFormat, OutputOptions, open_output, parse_output_spec};
use gloci_ops::dict::{DictOptions, dict, parse_alt_list};
use std::path::Path;

pub fn run(ctx: &Ctx, a: &DictArgs) -> anyhow::Result<()> {
    // The dictionary is plain text whatever -o is called (`.dict` is not an
    // extension gloci recognises), so the format is fixed.
    let spec = parse_output_spec(a.output.as_deref(), Some("tsv"))?;
    let alt = match &a.alt {
        None => None,
        Some(p) => {
            let text = std::fs::read(p).map_err(|source| IoError::Open {
                path: p.display().to_string(),
                source,
            })?;
            Some(parse_alt_list(&text))
        }
    };
    let o = DictOptions {
        hd_line: !a.no_header,
        assembly: a.assembly.clone(),
        species: a.species.clone(),
        uri: a.uri.clone(),
        alias: a.alias,
        alt,
    };
    let fasta = FastaStream::open(&a.fasta, ctx.threads)?;
    let oo = OutputOptions {
        threads: ctx.threads,
        level: 6,
        write_index: false,
    };
    let dest = super::dest_name(&spec);
    let mut out = open_output(&spec, OutFormat::Tsv, &oo)?;
    let path = (a.fasta != "-").then(|| Path::new(&a.fasta));
    dict(fasta, path, &o, &mut out).map_err(|e| e.with_output_name(&dest))?;
    out.finish()?;
    Ok(())
}
