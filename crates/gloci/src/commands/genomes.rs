use crate::cli::{GenomesAction, GenomesArgs};
use gloci_io::{builtin, builtin_genomes, unknown_genome};
use std::io::Write;

pub fn run(args: &GenomesArgs) -> anyhow::Result<()> {
    let mut out = std::io::BufWriter::new(std::io::stdout().lock());
    match &args.action {
        None | Some(GenomesAction::List) => {
            for b in builtin_genomes() {
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    b.id,
                    b.species,
                    b.assembly,
                    b.naming.as_str(),
                    b.len(),
                    b.total_length()
                )?;
            }
        }
        Some(GenomesAction::Show { name }) => {
            let b = builtin(name).ok_or_else(|| unknown_genome(name))?;
            let g = b.genome();
            for i in 0..g.len() {
                out.write_all(g.name(i))?;
                writeln!(out, "\t{}", g.length(i).unwrap_or(0))?;
            }
        }
    }
    out.flush()?;
    Ok(())
}
