use gloci_core::{Interval, SignFrame, Strand, distance, signed_distance};

fn iv(s: u64, e: u64) -> Interval {
    Interval::new(s, e).unwrap()
}

#[test]
fn unsigned_distance() {
    assert_eq!(distance(iv(5, 15), iv(10, 20), false), 0);
    assert_eq!(distance(iv(5, 10), iv(10, 15), false), 1); // bookended
    assert_eq!(distance(iv(0, 10), iv(15, 20), false), 6); // gap 5
    assert_eq!(distance(Interval::point(10), iv(5, 15), false), 0);
    assert_eq!(distance(Interval::point(5), iv(5, 15), false), 1);
    assert_eq!(distance(Interval::point(5), iv(5, 15), true), 0);
    assert_eq!(distance(Interval::point(10), Interval::point(12), false), 3);
}

#[test]
fn signed_distances() {
    let (a, up, down) = (iv(100, 200), iv(50, 60), iv(250, 260));
    use SignFrame::*;
    use Strand::*;
    assert_eq!(signed_distance(a, Plus, up, Plus, Ref, false), -41);
    assert_eq!(signed_distance(a, Minus, up, Plus, A, false), 41);
    assert_eq!(signed_distance(a, Plus, up, Plus, B, false), 41);
    assert_eq!(signed_distance(a, Plus, up, Minus, B, false), -41);
    assert_eq!(signed_distance(a, Plus, down, Plus, Ref, false), 51);
    assert_eq!(signed_distance(a, Minus, down, Plus, A, false), -51);
    assert_eq!(signed_distance(a, Plus, down, Plus, B, false), -51);
    assert_eq!(signed_distance(a, Unknown, up, Plus, A, false), -41); // D6
    assert_eq!(signed_distance(a, Plus, iv(150, 160), Minus, B, false), 0);
}
