use gloci_core::{Fraction, Interval, OverlapCriteria, Strand, StrandRule, Window, overlaps};

fn iv(s: u64, e: u64) -> Interval {
    Interval::new(s, e).unwrap()
}

#[test]
fn spec_4_3_worked_examples() {
    // (a, b, default, with --inclusive-points)
    let cases = [
        (iv(10, 10), iv(5, 15), true, true), // insertion site inside a feature IS found
        (iv(5, 5), iv(5, 15), false, true),
        (iv(15, 15), iv(5, 15), false, true),
        (iv(10, 10), iv(10, 10), true, true),
        (iv(10, 10), iv(11, 11), false, false),
        (iv(5, 10), iv(10, 15), false, false), // bookended
        (iv(5, 11), iv(10, 15), true, true),
    ];
    for (a, b, dflt, incl) in cases {
        assert_eq!(overlaps(a, b, false), dflt, "{a:?} {b:?}");
        assert_eq!(overlaps(b, a, false), dflt, "symmetric {a:?} {b:?}");
        assert_eq!(overlaps(a, b, true), incl, "inclusive {a:?} {b:?}");
    }
}

#[test]
fn fractions_are_exact() {
    assert!(Fraction::parse("0.3").unwrap().met(3, 10));
    assert!(!Fraction::parse("0.3").unwrap().met(2, 10));
    assert!(Fraction::parse("0.1").unwrap().met(1, 10));
    assert!(Fraction::parse("1").unwrap().met(10, 10) && Fraction::parse("1.0").is_ok());
    assert!(
        Fraction::parse("0.5")
            .unwrap()
            .met(u64::MAX / 2 + 1, u64::MAX)
    );
    for bad in ["0", "1.5", "-0.1", "abc", "5e-1", ""] {
        assert!(Fraction::parse(bad).is_err(), "{bad}");
    }
}

#[test]
fn criteria() {
    let (a, b) = (iv(0, 10), iv(5, 20)); // overlap 5: 0.5 of A, 1/3 of B
    let f = |s| Some(Fraction::parse(s).unwrap());
    let c = |m: OverlapCriteria| m.accepts(a, Strand::Plus, b, Strand::Plus);
    assert!(c(OverlapCriteria {
        min_frac_a: f("0.5"),
        ..Default::default()
    }));
    assert!(!c(OverlapCriteria {
        min_frac_a: f("0.5"),
        reciprocal: true,
        ..Default::default()
    }));
    assert!(c(OverlapCriteria {
        min_frac_a: f("0.9"),
        min_frac_b: f("0.3"),
        either: true,
        ..Default::default()
    }));
    assert!(!c(OverlapCriteria {
        min_bases: Some(6),
        ..Default::default()
    }));
    // Spec 4.3 (amended 2026-11-08, D34, superseding D5): a threshold counts as met
    // when EITHER side is a point, not only the point's own side. See overlap_p1.rs.
    let point_ok = OverlapCriteria {
        min_frac_a: f("0.9"),
        min_bases: Some(6),
        ..Default::default()
    };
    assert!(point_ok.accepts(Interval::point(10), Strand::Plus, b, Strand::Plus));
    let same = OverlapCriteria {
        strand: StrandRule::Same,
        ..Default::default()
    };
    assert!(
        !same.accepts(a, Strand::Plus, b, Strand::Minus)
            && !same.accepts(a, Strand::Unknown, b, Strand::Unknown)
    );
    let opp = OverlapCriteria {
        strand: StrandRule::Opposite,
        ..Default::default()
    };
    assert!(
        opp.accepts(a, Strand::Plus, b, Strand::Minus)
            && !opp.accepts(a, Strand::Unknown, b, Strand::Minus)
    );
}

#[test]
fn windows() {
    assert_eq!(
        Window {
            left: 5,
            right: 5,
            strand_aware: false
        }
        .apply(Interval::point(10), Strand::Plus),
        iv(5, 15)
    );
    assert_eq!(
        Window {
            left: 50,
            right: 0,
            strand_aware: false
        }
        .apply(iv(10, 20), Strand::Plus),
        iv(0, 20)
    );
    assert_eq!(
        Window {
            left: 5,
            right: 0,
            strand_aware: true
        }
        .apply(iv(10, 20), Strand::Minus),
        iv(10, 25)
    );
}
