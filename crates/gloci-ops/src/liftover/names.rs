//! Chromosome names for liftover (D174, R4).

use gloci_io::{ChainIndex, SourceRef};

/// How the target chromosome is named (`--chrom-style`; CrossMap's `--chromid`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ChromStyle {
    /// `chr` prefix as the input has it (`a`).
    #[default]
    AsInput,
    /// Every `chr` removed (`s`).
    Short,
    /// `chr` prepended unless present (`l`).
    Long,
    /// The chain's own name (`n`).
    Chain,
}

/// The chain's source chromosome for an input name: as written, then with every `chr` removed,
/// then with `chr` prepended (R4). `chrM` finds no `M`, and is never sent to `MT`.
pub fn source_name(chain: &ChainIndex, chrom: &str) -> Option<SourceRef> {
    chain
        .source(chrom)
        .or_else(|| chain.source(&chrom.replace("chr", "")))
        .or_else(|| chain.source(&format!("chr{chrom}")))
}

fn long(name: &str) -> String {
    if name.starts_with("chr") {
        name.to_string()
    } else {
        format!("chr{name}")
    }
}

fn short(name: &str) -> String {
    if name.starts_with("chr") {
        name.replace("chr", "")
    } else {
        name.to_string()
    }
}

/// The name written for `target`, a chain's target chromosome, of a record read as `input` (R4).
pub fn target_name(style: ChromStyle, input: &str, target: &str) -> String {
    match style {
        ChromStyle::AsInput if input.starts_with("chr") => long(target),
        ChromStyle::AsInput | ChromStyle::Short => short(target),
        ChromStyle::Long => long(target),
        ChromStyle::Chain => target.to_string(),
    }
}
