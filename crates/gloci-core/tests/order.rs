use gloci_core::{ChromTable, Genome, OrderChecker, OrderSource, OrderViolation};

#[test]
fn lexicographic() {
    let mut t = ChromTable::lexicographic();
    let mut c = OrderChecker::new();
    for (n, k) in [("chr1", 5), ("chr1", 5), ("chr10", 0), ("chr2", 3)] {
        let id = t.intern(n.as_bytes()).unwrap();
        c.check(&t, id, k).unwrap();
    }
    let mut c = OrderChecker::new();
    c.check(&t, t.get(b"chr2").unwrap(), 0).unwrap();
    assert!(matches!(
        c.check(&t, t.get(b"chr10").unwrap(), 0),
        Err(OrderViolation::ChromBackwards { .. })
    ));
}
#[test]
fn genome_order_and_reappearance() {
    let g = Genome::parse_sizes(b"chr1\t1000\nchr2\t1000\nchr10\t1000\n").unwrap();
    let t = ChromTable::closed(&g, OrderSource::Lexicographic);
    let id = |n: &[u8]| t.get(n).unwrap();
    let mut c = OrderChecker::new();
    c.check(&t, id(b"chr2"), 0).unwrap();
    c.check(&t, id(b"chr10"), 0).unwrap();
    assert!(c.check(&t, id(b"chr2"), 50).is_err()); // chr2 reappears
    let mut c = OrderChecker::new();
    c.check(&t, id(b"chr1"), 100).unwrap();
    assert!(matches!(
        c.check(&t, id(b"chr1"), 99),
        Err(OrderViolation::StartBackwards {
            prev: 100,
            next: 99,
            ..
        })
    ));
    c.check(&t, id(b"chr1"), 120).unwrap(); // continues after violation
}
#[test]
fn first_appearance() {
    let mut t = ChromTable::first_appearance();
    let mut c = OrderChecker::new();
    for n in ["chr2", "chr1", "chr10"] {
        let id = t.intern(n.as_bytes()).unwrap();
        c.check(&t, id, 0).unwrap();
    }
    assert!(c.check(&t, t.get(b"chr1").unwrap(), 0).is_err());
}

#[test]
fn karyotype_checker() {
    let mut t = ChromTable::karyotype();
    let mut c = OrderChecker::new();
    for (n, k) in [
        ("chr1", 5),
        ("chr2", 0),
        ("chr10", 3),
        ("chrX", 0),
        ("chrM", 0),
    ] {
        let id = t.intern(n.as_bytes()).unwrap();
        c.check(&t, id, k).unwrap();
    }
    let mut t = ChromTable::karyotype();
    let mut c = OrderChecker::new();
    let (c10, c2) = (t.intern(b"chr10").unwrap(), t.intern(b"chr2").unwrap());
    c.check(&t, c10, 0).unwrap();
    assert!(matches!(
        c.check(&t, c2, 0),
        Err(OrderViolation::ChromBackwards { .. })
    ));
}
