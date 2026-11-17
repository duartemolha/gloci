//! Chromosome names for liftover (D174, R4).

use gloci_io::{ChainIndex, MtRule, SourceRef};

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

/// Is `chrom` hg19's mitochondrial chromosome, which the built-in GRCh37-to-GRCh38 chain refuses
/// because it is not GRCh37's `MT` (D174, R13)?
pub fn refused_chrm(chain: &ChainIndex, chrom: &str) -> bool {
    chain.mt_rule() == MtRule::Hg19 && chrom.strip_prefix("chr").unwrap_or(chrom) == "M"
}

/// The chain's source chromosome for an input name: as written, then with every `chr` removed,
/// then with `chr` prepended (R4). `chrM` finds no `M`, and is never sent to `MT`, except by a
/// built-in chain whose source is GRCh38 or CHM13v2, where it finds `MT` (R13).
pub fn source_name(chain: &ChainIndex, chrom: &str) -> Option<SourceRef> {
    if refused_chrm(chain, chrom) {
        return None;
    }
    chain
        .source(chrom)
        .or_else(|| chain.source(&chrom.replace("chr", "")))
        .or_else(|| chain.source(&format!("chr{chrom}")))
        .or_else(|| {
            (chain.mt_rule() == MtRule::ChrM && chrom == "chrM")
                .then(|| chain.source("MT"))
                .flatten()
        })
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

/// A built-in chain's own target name for the mitochondrial chromosome: GRCh38's `MT` is hg38's
/// `chrM`, so `as-input` writes it back as `chrM` for an input with a `chr` prefix (R13).
pub(super) fn mt_target(
    chain: &ChainIndex,
    style: ChromStyle,
    input_chr: bool,
    target: &str,
) -> bool {
    chain.mt_rule() == MtRule::ChrM && style == ChromStyle::AsInput && input_chr && target == "MT"
}
