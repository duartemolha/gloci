use gloci_core::{FootprintError, Interval, Variant, VcfSpan, alt_footprints, footprint};

fn iv(s: u64, e: u64) -> Interval {
    Interval::new(s, e).unwrap()
}
fn fp(
    pos: u64,
    r: &str,
    a: &str,
    end: Option<u64>,
    svlen: Option<&str>,
    span: VcfSpan,
) -> Result<Interval, FootprintError> {
    footprint(
        &Variant {
            pos,
            reference: r.as_bytes(),
            alts: a.as_bytes(),
            end,
            svlen: svlen.map(str::as_bytes),
        },
        span,
    )
}

#[test]
fn spec_4_2_table() {
    let p = Interval::point;
    let cases = [
        (100, "A", "G", None, iv(99, 100), iv(99, 100)),
        (100, "AC", "GT", None, iv(99, 101), iv(99, 101)),
        (100, "ACGT", "A", None, iv(99, 103), iv(100, 103)),
        (100, "A", "ACGT", None, iv(99, 100), p(100)),
        (100, "ATG", "AC", None, iv(99, 102), iv(100, 102)),
        (100, "ACGT", "A,AC", None, iv(99, 103), iv(100, 103)),
        (100, "N", "<DEL>", Some(500), iv(99, 500), iv(100, 500)),
        (100, "N", "<INS>", None, iv(99, 100), p(100)),
        (100, "N", "<DUP>", Some(700), iv(99, 700), iv(100, 700)),
        (100, "N", "<INV>", Some(700), iv(99, 700), iv(100, 700)),
        (100, "N", "<CNV>", Some(700), iv(99, 700), iv(100, 700)),
        (100, "G", "G]17:198982]", None, iv(99, 100), iv(99, 100)),
        (100, "A", "<*>", Some(250), iv(99, 250), iv(99, 250)),
        (100, "A", "<NON_REF>", Some(250), iv(99, 250), iv(99, 250)),
        (100, "A", ".", None, iv(99, 100), iv(99, 100)),
    ];
    for (pos, r, a, end, want_ref, want_allele) in cases {
        assert_eq!(
            fp(pos, r, a, end, None, VcfSpan::Ref).unwrap(),
            want_ref,
            "ref {r}>{a}"
        );
        assert_eq!(
            fp(pos, r, a, end, None, VcfSpan::Allele).unwrap(),
            want_allele,
            "allele {r}>{a}"
        );
    }
}

#[test]
fn further_rules() {
    use VcfSpan::*;
    assert_eq!(
        fp(100, "N", "<DEL>", None, Some("-400"), Allele).unwrap(),
        iv(100, 500)
    ); // SVLEN fallback
    assert_eq!(
        fp(100, "N", "<DEL>", Some(300), Some("-400"), Allele).unwrap(),
        iv(100, 300)
    ); // END wins
    assert_eq!(
        fp(100, "N", "<DEL>,<DEL>", None, Some("-100,-300"), Allele).unwrap(),
        iv(100, 400)
    );
    assert!(matches!(
        fp(100, "N", "<DEL>", None, None, Allele),
        Err(FootprintError::MissingSvEnd { .. })
    ));
    assert!(matches!(
        fp(100, "N", "<DEL>", Some(50), None, Allele),
        Err(FootprintError::EndBeforePos { .. })
    ));
    assert_eq!(
        fp(100, "ACGT", "A,*", None, None, Allele).unwrap(),
        iv(100, 103)
    ); // * ignored
    assert_eq!(
        fp(100, "A", "G,<NON_REF>", Some(900), None, Allele).unwrap(),
        iv(99, 100)
    ); // NON_REF ignored
    assert_eq!(fp(100, "A", "*", None, None, Allele).unwrap(), iv(99, 100)); // * alone
    assert_eq!(
        fp(100, "acgt", "A", None, None, Allele).unwrap(),
        iv(100, 103)
    ); // case-insensitive
    assert_eq!(
        fp(100, "AA", "A", None, None, Allele).unwrap(),
        iv(100, 101)
    ); // prefix first
    assert_eq!(
        fp(100, "N", "<INS:ME:ALU>", Some(100), None, Allele).unwrap(),
        Interval::point(100)
    );
    assert_eq!(
        fp(100, "N", "<DEL:ME>", Some(400), None, Allele).unwrap(),
        iv(100, 400)
    );
    assert_eq!(fp(100, "G", "G.", None, None, Allele).unwrap(), iv(99, 100)); // single breakend
    assert_eq!(
        fp(100, "ACGTT", "A", Some(999), None, Ref).unwrap(),
        iv(99, 104)
    ); // END ignored, not symbolic
    assert!(matches!(
        fp(0, "A", "G", None, None, Allele),
        Err(FootprintError::PosZero)
    ));
    let per_alt = alt_footprints(&Variant {
        pos: 100,
        reference: b"ACGT",
        alts: b"A,AC,*",
        end: None,
        svlen: None,
    })
    .unwrap();
    assert_eq!(per_alt, vec![Some(iv(100, 103)), Some(iv(101, 103)), None]);
}

#[test]
fn error_cases() {
    use VcfSpan::*;
    assert!(matches!(
        fp(100, "", "G", None, None, Allele),
        Err(FootprintError::EmptyRef)
    ));
    assert!(matches!(
        fp(100, "N", "<DEL>", None, Some("x"), Allele),
        Err(FootprintError::BadSvlen { .. })
    ));
    assert!(matches!(
        fp(100, "N", "<DEL>,<DEL>,<DUP>", None, Some("-1,-2"), Allele),
        Err(FootprintError::SvlenCount { .. })
    ));
    // Under ref a missing END and SVLEN is not an error: the REF span stands.
    assert_eq!(
        fp(100, "N", "<DEL>", None, Some("."), Ref).unwrap(),
        iv(99, 100)
    );
    assert_eq!(fp(100, "N", "<DEL>", None, None, Ref).unwrap(), iv(99, 100));
    assert_eq!(
        fp(100, "NN", "<DUP>,<CNV>", None, None, Ref).unwrap(),
        iv(99, 101)
    );
    assert!(matches!(
        fp(100, "N", "<DEL>", None, Some("."), Allele),
        Err(FootprintError::MissingSvEnd { .. })
    ));
    // Other symbolic ALTs never need either.
    for alt in ["<INS>", "<*>", "<NON_REF>"] {
        assert!(fp(100, "N", alt, None, None, Allele).is_ok(), "{alt}");
    }
    assert_eq!(
        fp(100, "N", "<DEL>", None, Some("-400"), Ref).unwrap(),
        iv(99, 500)
    );
    assert_eq!(
        fp(100, "N", "<DEL>", None, None, Allele)
            .unwrap_err()
            .to_string(),
        "symbolic ALT <DEL> needs INFO/END or INFO/SVLEN to give it a span; add INFO/END or \
INFO/SVLEN, use --vcf-span ref, or skip invalid records where the command supports it"
    );
}

/// Found by the structured VCF fuzz test: POS near u64::MAX overflowed in debug
/// builds (and wrapped silently in release).
#[test]
fn pos_plus_ref_length_overflow_is_an_error_not_a_panic() {
    for alts in [&b"A"[..], b"<DEL>", b"A]chr1:1]", b"*"] {
        let v = Variant {
            pos: u64::MAX,
            reference: b"ACGT",
            alts,
            end: None,
            svlen: None,
        };
        assert!(matches!(
            footprint(&v, VcfSpan::Allele),
            Err(FootprintError::PosOverflow { .. })
        ));
        assert!(alt_footprints(&v).is_err());
    }
}
