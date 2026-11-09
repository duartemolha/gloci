use gloci_core::VariantKind::{self, *};
use gloci_core::{Variant, alt_kinds, trimmed_alleles};

fn kinds(r: &str, a: &str) -> Vec<VariantKind> {
    alt_kinds(&Variant {
        pos: 100,
        reference: r.as_bytes(),
        alts: a.as_bytes(),
        end: None,
        svlen: None,
    })
}

#[test]
fn pos_zero_does_not_underflow() {
    // Task 26 item 8: POS = 0 saturating arithmetic pos.saturating_sub(1)
    let v = Variant {
        pos: 0,
        reference: b"A",
        alts: b"G",
        end: None,
        svlen: None,
    };
    assert_eq!(alt_kinds(&v), [Snv]);
}

#[test]
fn spec_4_2_rows() {
    assert_eq!(kinds("A", "G"), [Snv]);
    assert_eq!(kinds("AC", "GT"), [Mnv]);
    assert_eq!(kinds("ACGT", "A"), [Del]);
    assert_eq!(kinds("A", "ACGT"), [Ins]);
    assert_eq!(kinds("ATG", "AC"), [Other]);
    assert_eq!(kinds("N", "<DEL>"), [Sv]);
    assert_eq!(kinds("G", "G]17:198982]"), [Sv]);
    assert_eq!(kinds("A", "<*>"), [Ref]);
    assert_eq!(kinds("A", "<NON_REF>"), [Ref]);
    assert_eq!(kinds("A", "."), [Ref]);
    assert_eq!(kinds("ACGT", "A,AC"), [Del, Del]);
    assert_eq!(kinds("A", "G,<DEL>,."), [Snv, Sv, Ref]);
}

#[test]
fn star_allele_is_other() {
    assert_eq!(kinds("A", "*"), [Other]);
}

#[test]
fn trimming() {
    assert_eq!(
        trimmed_alleles(99, b"ATG", b"AC"),
        (100, &b"TG"[..], &b"C"[..])
    );
    assert_eq!(
        trimmed_alleles(99, b"A", b"ACGT"),
        (100, &b""[..], &b"CGT"[..])
    );
    assert_eq!(
        trimmed_alleles(99, b"ACGT", b"A"),
        (100, &b"CGT"[..], &b""[..])
    );
    assert_eq!(trimmed_alleles(99, b"A", b"G"), (99, &b"A"[..], &b"G"[..]));
}
