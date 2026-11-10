//! Allele matching (`--match alleles`), shared by `intersect` and
//! `subtract`: which VCF records share a trimmed (position, REF, ALT), the
//! option checks that go with it, and the candidate window that finds the
//! pairs.

use super::RelateOptions;
use gloci_core::{StrandRule, Window, trimmed_alleles};
use gloci_io::{Format, IoError, UsageError, field};

fn usage(msg: impl Into<String>) -> IoError {
    UsageError(msg.into()).into()
}

fn is_variant(f: Format) -> bool {
    matches!(f, Format::Vcf | Format::Bcf)
}

/// Allele matching compares ALT alleles, so A and every B must be VCF/BCF.
pub(crate) fn check_alleles_formats(a: Format, b: &[Format]) -> Result<(), IoError> {
    if !(is_variant(a) && b.iter().all(|f| is_variant(*f))) {
        return Err(usage(
            "--match alleles compares ALT alleles, so -a and every -b must be VCF or BCF files; \
use --match overlap for other formats",
        ));
    }
    Ok(())
}

/// Allele matching decides by the alleles alone, so the overlap criteria
/// are refused.
pub(crate) fn check_alleles_flags(r: &RelateOptions) -> Result<(), IoError> {
    let c = &r.criteria;
    let given = [
        (c.min_frac_a.is_some(), "-f"),
        (c.min_frac_b.is_some(), "-F"),
        (c.reciprocal, "-r"),
        (c.either, "-e"),
        (c.min_bases.is_some(), "--min-bases"),
        (c.strand != StrandRule::Any, "-s/-S"),
        (c.inclusive_points, "--inclusive-points"),
        (r.split, "--split"),
        (r.window.is_some(), "-w/--left/--right"),
    ];
    if let Some((_, flag)) = given.iter().find(|(on, _)| *on) {
        return Err(usage(format!(
            "{flag} cannot be combined with --match alleles: allele matching compares the \
alleles of two records, not how much they overlap; drop {flag} or use --match overlap"
        )));
    }
    Ok(())
}

/// Set up the relation to yield candidates for allele matching: a shared
/// trimmed allele lies at most one base from each record's footprint (an
/// insertion is a point just past an anchor base, and a footprint can end
/// where the other starts), so pad B by one base on each side. The
/// candidates are then accepted by allele equality alone. Points count at
/// the edge of the padded interval (`inclusive_points`): the padding
/// saturates at position 0, so a point at 0 (an insertion at POS 1) would
/// otherwise sit on the edge of its own padded copy and never be a
/// candidate.
pub(crate) fn allele_candidates(relate: &mut RelateOptions) {
    relate.criteria.inclusive_points = true;
    relate.window = Some(Window {
        left: 1,
        right: 1,
        strand_aware: false,
    });
}

/// Calls `f(pos, ref, alt)` with the trimmed (0-based position, REF, ALT)
/// of each usable ALT of a VCF line until it returns true. Symbolic,
/// breakend, `*` and missing ALTs have no sequence to compare and are
/// skipped. Allocates nothing; compare the slices case-insensitively.
fn any_allele(line: &[u8], mut f: impl FnMut(u64, &[u8], &[u8]) -> bool) -> bool {
    let pos = field(line, 1)
        .and_then(|p| std::str::from_utf8(p).ok())
        .and_then(|p| p.parse::<u64>().ok());
    let (Some(pos), Some(reference), Some(alt)) = (pos, field(line, 3), field(line, 4)) else {
        return false;
    };
    alt.split(|&c| c == b',')
        .filter(|a| {
            !a.is_empty()
                && !a
                    .iter()
                    .any(|c| matches!(c, b'<' | b'[' | b']' | b'*' | b'.'))
        })
        .any(|a| {
            let (p, r, a) = trimmed_alleles(pos.saturating_sub(1), reference, a);
            f(p, r, a)
        })
}

/// Do two VCF lines share a trimmed (position, REF, ALT)?
pub(crate) fn share_allele(a: &[u8], b: &[u8]) -> bool {
    any_allele(b, |bp, br, ba| {
        any_allele(a, |ap, ar, aa| {
            ap == bp && ar.eq_ignore_ascii_case(br) && aa.eq_ignore_ascii_case(ba)
        })
    })
}
