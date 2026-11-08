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
