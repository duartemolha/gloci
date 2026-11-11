//! `gloci flags`: SAM FLAGS as numbers and names (D136).
use crate::cli::FlagsArgs;
use gloci_ops::aln_flags::{flags_line, flags_table};
use std::io::Write;

pub fn run(a: &FlagsArgs) -> anyhow::Result<()> {
    // Every argument is checked before anything is written.
    let text = if a.flags.is_empty() {
        flags_table()
    } else {
        a.flags
            .iter()
            .map(|f| flags_line(f))
            .collect::<Result<String, _>>()?
    };
    let mut out = std::io::stdout().lock();
    out.write_all(text.as_bytes())?;
    out.flush()?;
    Ok(())
}
