use gloci_core::{Fraction, Interval, OverlapCriteria, Strand, blocks_overlap, blocks_overlap_len};

fn iv(s: u64, e: u64) -> Interval {
    Interval::new(s, e).unwrap()
}

#[test]
fn d34_thresholds_met_when_either_side_is_point() {
    let u = Strand::Unknown;
    let half = Fraction::parse("0.5").unwrap();
    let c = OverlapCriteria {
        min_frac_b: Some(half),
        ..Default::default()
    };
    assert!(c.accepts(Interval::point(10), u, iv(5, 15), u)); // point A, -F 0.5
    let c = OverlapCriteria {
        min_bases: Some(3),
        ..Default::default()
    };
    assert!(c.accepts(Interval::point(10), u, iv(5, 15), u));
    let c = OverlapCriteria {
        min_frac_a: Some(Fraction::parse("1").unwrap()),
        ..Default::default()
    };
    assert!(c.accepts(iv(5, 15), u, Interval::point(10), u)); // point B, -f 1
    assert!(!c.accepts(Interval::point(5), u, iv(5, 15), u)); // edge: no overlap at all
}

#[test]
fn d34_reciprocal_with_only_f_big() {
    let u = Strand::Unknown;
    let half = Fraction::parse("0.5").unwrap();
    let c = OverlapCriteria {
        min_frac_b: Some(half),
        reciprocal: true,
        ..Default::default()
    };
    assert!(!c.accepts(iv(0, 10), u, iv(5, 100), u)); // 5/10 ok, 5/95 not
    assert!(c.accepts(iv(0, 10), u, iv(5, 10), u));
}

#[test]
fn fraction_equality_by_value() {
    assert_eq!(
        Fraction::parse("1").unwrap(),
        Fraction::parse("1.0").unwrap()
    );
    assert_eq!(
        Fraction::parse("50%").unwrap(),
        Fraction::parse("0.5").unwrap()
    );
    assert_ne!(
        Fraction::parse("0.5").unwrap(),
        Fraction::parse("0.51").unwrap()
    );
    assert!(Fraction::parse("50%").unwrap().met(5, 10));
    assert!(Fraction::parse("0%").is_err() && Fraction::parse("101%").is_err());
}

#[test]
fn blocks() {
    assert_eq!(
        blocks_overlap_len(&[iv(0, 10), iv(20, 30)], &[iv(8, 22)]),
        4
    );
    assert!(!blocks_overlap(&[iv(0, 10)], &[iv(10, 20)], false));
    assert!(blocks_overlap(&[iv(0, 10)], &[iv(9, 20)], false));
    assert!(blocks_overlap(&[iv(0, 10)], &[Interval::point(10)], true));
    assert!(!blocks_overlap(&[iv(0, 10)], &[Interval::point(10)], false));
    assert!(blocks_overlap(
        &[iv(0, 5), iv(20, 30)],
        &[iv(25, 26)],
        false
    ));
    assert!(!blocks_overlap(&[], &[iv(0, 5)], false));
    assert_eq!(blocks_overlap_len(&[iv(0, 5)], &[iv(5, 9)]), 0);
}

#[test]
fn parse_ratio_bounds() {
    use gloci_core::Fraction;
    for ok in ["0", "0.5", "50%", "1", "2", "150%", "18446744073709551615"] {
        assert!(Fraction::parse_ratio(ok).is_ok(), "{ok}");
    }
    for bad in [
        "",
        ".",
        "-1",
        "abc",
        "1e3",
        "18446744073709551616",
        "0.1234567890123456789",
    ] {
        assert!(Fraction::parse_ratio(bad).is_err(), "{bad}");
    }
    // The unit parser keeps its (0, 1] range.
    assert!(Fraction::parse("2").is_err());
    assert!(Fraction::parse("0").is_err());
}

#[test]
fn floor_mul_exact_and_saturating() {
    use gloci_core::Fraction;
    let f = |s| Fraction::parse_ratio(s).unwrap();
    assert_eq!(f("0.5").floor_mul(7), 3);
    assert_eq!(f("0.1").floor_mul(30), 3);
    assert_eq!(f("0").floor_mul(u64::MAX), 0);
    assert_eq!(f("1").floor_mul(u64::MAX), u64::MAX);
    assert_eq!(f("0.5").floor_mul(u64::MAX), u64::MAX / 2);
    assert_eq!(f("2").floor_mul(u64::MAX), u64::MAX, "saturates");
    assert_eq!(f("18446744073709551615").floor_mul(u64::MAX), u64::MAX);
}
