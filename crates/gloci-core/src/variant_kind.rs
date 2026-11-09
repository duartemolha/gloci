//! Variant kinds per ALT (spec 4.2), sharing the footprint code's ALT splitting.

use crate::footprint::{Alt, Variant, classify, split_alts};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VariantKind {
    Snv,
    Mnv,
    Ins,
    Del,
    /// Symbolic or breakend ALT.
    Sv,
    /// `.`, `<*>`, `<NON_REF>`: no variant allele.
    Ref,
    /// Complex sequence ALTs (and the `*` spanning-deletion allele).
    Other,
}

/// Strip the common prefix, then the common suffix (ASCII case-insensitive).
/// Returns the 0-based start of the trimmed REF and the trimmed REF and ALT.
pub fn trimmed_alleles<'a>(pos0: u64, ref_: &'a [u8], alt: &'a [u8]) -> (u64, &'a [u8], &'a [u8]) {
    let prefix = ref_
        .iter()
        .zip(alt)
        .take_while(|(x, y)| x.eq_ignore_ascii_case(y))
        .count();
    let (r, a) = (&ref_[prefix..], &alt[prefix..]);
    let suffix = r
        .iter()
        .rev()
        .zip(a.iter().rev())
        .take_while(|(x, y)| x.eq_ignore_ascii_case(y))
        .count();
    (
        pos0 + prefix as u64,
        &r[..r.len() - suffix],
        &a[..a.len() - suffix],
    )
}

/// The kind of one ALT (a single comma-split element) of a record whose
/// 0-based start is `pos0`. Does not allocate.
pub fn alt_kind(pos0: u64, reference: &[u8], alt: &[u8]) -> VariantKind {
    if alt == b"*" {
        return VariantKind::Other;
    }
    match classify(alt) {
        Alt::None | Alt::Gvcf => VariantKind::Ref,
        Alt::Ins | Alt::Sv(_) | Alt::Breakend => VariantKind::Sv,
        Alt::Seq(a) => {
            let (_, r, a) = trimmed_alleles(pos0, reference, a);
            match (r.len(), a.len()) {
                (0, 0) => VariantKind::Ref,
                (0, _) => VariantKind::Ins,
                (_, 0) => VariantKind::Del,
                (1, 1) => VariantKind::Snv,
                (x, y) if x == y => VariantKind::Mnv,
                _ => VariantKind::Other,
            }
        }
    }
}

/// The kind of each ALT of `v`, in ALT order.
pub fn alt_kinds(v: &Variant) -> Vec<VariantKind> {
    split_alts(v.alts)
        .into_iter()
        .map(|raw| alt_kind(v.pos.saturating_sub(1), v.reference, raw))
        .collect()
}
