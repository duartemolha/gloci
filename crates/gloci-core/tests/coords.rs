use gloci_core::*;
use std::cmp::Ordering;

#[test]
fn coordinate_parsing() {
    assert_eq!(parse_coord(b"007").unwrap(), 7);
    assert_eq!(
        parse_coord(b"12a4").unwrap_err().to_string(),
        r#""12a4" is not an integer"#
    );
    assert_eq!(
        parse_coord(b"-1").unwrap_err().to_string(),
        r#""-1" is negative"#
    );
    assert!(
        parse_coord(b"18446744073709551616").is_err()
            && parse_coord(b"").is_err()
            && parse_coord(b"+5").is_err()
    );
    assert_eq!(
        parse_coord(b"18446744073709551616")
            .unwrap_err()
            .to_string(),
        r#""18446744073709551616" is too large"#
    );
}

#[test]
fn intervals_and_formatting() {
    assert!(Interval::new(10, 5).is_err());
    assert_eq!(
        Interval::new(10, 5).unwrap_err().to_string(),
        "start 10 is greater than end 5"
    );
    assert!(Interval::new(10, 10).unwrap().is_point());
    assert_eq!(Interval::new(3, 10).unwrap().len(), 7);
    assert_eq!(Interval::point(4), Interval::new(4, 4).unwrap());
    assert_eq!(
        one_based_closed(100, 200).unwrap(),
        Interval::new(99, 200).unwrap()
    ); // GFF, region strings
    assert!(one_based_closed(0, 5).is_err());
    assert_eq!(
        one_based_closed(0, 5).unwrap_err().to_string(),
        "start must be >= 1 (1-based coordinates)"
    );
    assert_eq!(thousands(1204331), "1,204,331");
    assert_eq!(thousands(0), "0");
    assert_eq!(thousands(999), "999");
    assert_eq!(Strand::parse(b"?"), Some(Strand::Unknown));
    assert_eq!(Strand::parse(b"+"), Some(Strand::Plus));
    assert_eq!(Strand::parse(b"-"), Some(Strand::Minus));
    assert_eq!(Strand::parse(b"."), Some(Strand::Unknown));
    assert_eq!(Strand::parse(b"x"), None);
}

#[test]
fn region_strings() {
    let none = |_: &[u8]| false;
    assert_eq!(
        Region::parse("chr1:1,000-2,000", none).unwrap(),
        Region {
            chrom: b"chr1".to_vec(),
            start: 999,
            end: Some(2000)
        }
    );
    assert_eq!(
        Region::parse("chr1:100-100", none).unwrap().interval(None),
        Interval::new(99, 100).unwrap()
    );
    assert_eq!(Region::parse("chr1", none).unwrap().end, None);
    assert_eq!(
        Region::parse("chr1:500", none).unwrap(),
        Region {
            chrom: b"chr1".to_vec(),
            start: 499,
            end: None
        }
    );
    assert_eq!(Region::parse("chr1:500-", none).unwrap().start, 499);
    assert!(
        Region::parse("chr1:0-10", none).is_err() && Region::parse("chr1:20-10", none).is_err()
    );
    let hla = |n: &[u8]| n == b"HLA-A*01:01:01:01";
    assert_eq!(
        Region::parse("HLA-A*01:01:01:01", hla).unwrap().chrom,
        b"HLA-A*01:01:01:01"
    );
    assert_eq!(
        Region::parse("{HLA-A*01:01:01:01}:5-10", none)
            .unwrap()
            .start,
        4
    );
    let both = |n: &[u8]| n == b"x:1" || n == b"x";
    assert!(Region::parse("x:1", both).is_err()); // ambiguous
}

#[test]
fn region_interval_open_end() {
    let none = |_: &[u8]| false;
    let r = Region::parse("chr1:5", none).unwrap();
    assert_eq!(r.interval(None), Interval::new(4, u64::MAX).unwrap());
    assert_eq!(r.interval(Some(100)), Interval::new(4, 100).unwrap());
}

#[test]
fn chromosome_tables() {
    let mut t = ChromTable::lexicographic();
    let (c2, c10, c1) = (
        t.intern(b"chr2").unwrap(),
        t.intern(b"chr10").unwrap(),
        t.intern(b"chr1").unwrap(),
    );
    assert_eq!(t.ranks(), vec![2, 1, 0]);
    assert_eq!(t.cmp(c10, c2), Ordering::Less);
    assert_eq!(t.cmp(c1, c10), Ordering::Less); // prefix sorts first
    assert_eq!(t.intern(b"chr2").unwrap(), c2);
    assert_eq!(t.name(c10), b"chr10");
    assert_eq!(t.len(), 3);
    assert_eq!(t.order_source(), &OrderSource::Lexicographic);

    let mut f = ChromTable::first_appearance();
    let (a, b) = (f.intern(b"z").unwrap(), f.intern(b"a").unwrap());
    assert_eq!(f.cmp(a, b), Ordering::Less);
    assert_eq!(f.ranks(), vec![0, 1]);
    assert_eq!(f.get(b"a"), Some(b));
    assert_eq!(f.get(b"q"), None);
}

#[test]
fn genome_and_closed_tables() {
    let g = Genome::parse_sizes(b"chr1\t100\nchr2\t200\nchr10\t50\n").unwrap();
    let mut t = ChromTable::closed(&g, OrderSource::Lexicographic);
    assert!(t.intern(b"chrUn").is_err());
    assert_eq!(
        t.cmp(t.get(b"chr2").unwrap(), t.get(b"chr10").unwrap()),
        Ordering::Less
    );
    assert_eq!(g.len(), 3);
    assert_eq!(g.name(1), b"chr2");
    assert_eq!(g.length(1), Some(200));
    assert!(Genome::parse_sizes(b"chr1\t1\nchr1\t2\n").is_err()); // duplicate
    assert_eq!(
        Genome::parse_dict(b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n")
            .unwrap()
            .length(0),
        Some(248956422)
    );
}

#[test]
fn compare_names_in_table_order() {
    let lex = ChromTable::lexicographic();
    assert_eq!(lex.cmp_names(b"chr10", b"chr2"), Some(Ordering::Less));
    let kar = ChromTable::karyotype();
    assert_eq!(kar.cmp_names(b"chr10", b"chr2"), Some(Ordering::Greater));
    assert_eq!(kar.cmp_names(b"chrX", b"chrX"), Some(Ordering::Equal));
    assert!(kar.is_empty()); // nothing was interned
    let g = Genome::parse_sizes(b"chr2\t100\nchr1\t100\n").unwrap();
    let closed = ChromTable::closed(&g, OrderSource::GenomeFile("g".into()));
    assert_eq!(closed.cmp_names(b"chr2", b"chr1"), Some(Ordering::Less));
    assert_eq!(closed.cmp_names(b"chr2", b"chr3"), None);
    let fa = ChromTable::first_appearance();
    assert_eq!(fa.cmp_names(b"a", b"b"), None);
}

#[test]
fn region_sets_test_each_region() {
    let iv = |s, e| Interval { start: s, end: e };
    let set = RegionSet::new(vec![
        iv(20, 30),
        iv(10, 20),
        iv(10, 20),
        iv(0, 100),
        iv(40, 40),
    ]);
    assert_eq!(
        set.regions(),
        [iv(0, 100), iv(10, 20), iv(20, 30), iv(40, 40)]
    );
    assert_eq!(set.end(), 100);
    assert!(set.matches(iv(99, 120), false));
    assert!(!set.matches(iv(100, 120), false));
    assert!(set.matches(iv(100, 100), true));
    assert!(set.matches(iv(40, 40), false)); // equal points
    // a point between two bookended regions is in neither
    let set = RegionSet::new(vec![iv(10, 20), iv(20, 30)]);
    assert!(!set.matches(iv(20, 20), false));
    assert!(set.matches(iv(20, 20), true));
    assert!(set.matches(iv(19, 21), false));
    // a long region early on is still tested for a late span
    let set = RegionSet::new(vec![iv(0, 1000), iv(5, 6), iv(7, 8)]);
    assert!(set.matches(iv(500, 501), false));
    assert!(!set.matches(iv(1000, 1001), false));
    assert!(RegionSet::default().is_empty());
    assert_eq!(RegionSet::default().end(), 0);
    assert!(!RegionSet::default().matches(iv(0, 5), true));
}

#[test]
fn many_contigs_stay_fast() {
    let mut t = ChromTable::lexicographic();
    for i in 0..100_000 {
        t.intern(format!("scaffold_{i}").as_bytes()).unwrap();
    }
    assert_eq!(t.ranks().len(), 100_000);
}

fn ksort(v: &[&str]) -> Vec<String> {
    let mut v: Vec<_> = v.to_vec();
    v.sort_by(|a, b| karyotype_cmp(a.as_bytes(), b.as_bytes()));
    v.iter().map(|s| s.to_string()).collect()
}

#[test]
fn karyotype_order() {
    assert_eq!(
        ksort(&[
            "chrM",
            "chr10",
            "chr1_KI270706v1_random",
            "chr2",
            "chrY",
            "chr1",
            "chrX",
            "chrUn_GL000220v1",
            "chr22",
            "chrEBV",
            "HLA-A*01:01:01:01"
        ]),
        [
            "chr1",
            "chr2",
            "chr10",
            "chr22",
            "chrX",
            "chrY",
            "chrM",
            "chr1_KI270706v1_random",
            "chrEBV",
            "HLA-A*01:01:01:01",
            "chrUn_GL000220v1"
        ]
    );
    assert_eq!(
        ksort(&["MT", "X", "10", "2", "1", "Y"]),
        ["1", "2", "10", "X", "Y", "MT"]
    ); // Ensembl
    assert_eq!(
        ksort(&[
            "chrX", "chr4", "chr3R", "chr2L", "chrY", "chr3L", "chr2R", "chrM"
        ]),
        [
            "chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY", "chrM"
        ]
    ); // dm6
    assert_eq!(
        ksort(&["chr2B", "chr1", "chr2A", "chr3"]),
        ["chr1", "chr2A", "chr2B", "chr3"]
    ); // panTro
    assert_ne!(karyotype_cmp(b"1", b"chr1"), Ordering::Equal);
    assert_ne!(karyotype_cmp(b"chr01", b"chr1"), Ordering::Equal);
    assert_eq!(
        karyotype_cmp(
            b"chr99999999999999999999999",
            b"chr100000000000000000000000"
        ),
        Ordering::Less
    ); // no overflow
}

proptest::proptest! {
    #[test]
    fn karyotype_cmp_is_a_total_order(a in "[a-zA-Z0-9_]{0,8}", b in "[a-zA-Z0-9_]{0,8}", c in "[a-zA-Z0-9_]{0,8}") {
        let (a, b, c) = (a.as_bytes(), b.as_bytes(), c.as_bytes());
        proptest::prop_assert_eq!(karyotype_cmp(a, b), karyotype_cmp(b, a).reverse());
        proptest::prop_assert_eq!(karyotype_cmp(a, b) == Ordering::Equal, a == b);
        if karyotype_cmp(a, b) != Ordering::Greater && karyotype_cmp(b, c) != Ordering::Greater {
            proptest::prop_assert_ne!(karyotype_cmp(a, c), Ordering::Greater);
        }
    }
}

#[test]
fn karyotype_table_ranks() {
    let mut t = ChromTable::karyotype();
    for n in ["chr2", "chr10", "chr1", "chrX"] {
        t.intern(n.as_bytes()).unwrap();
    }
    assert_eq!(t.ranks(), vec![1, 2, 0, 3]);
    assert_eq!(*t.order_source(), OrderSource::Karyotype);
}

#[test]
fn chrom_order_parses() {
    assert_eq!("karyotype".parse(), Ok(ChromOrder::Karyotype));
    assert_eq!("lex".parse(), Ok(ChromOrder::Lex));
    assert!("natural".parse::<ChromOrder>().is_err());
    assert_eq!(ChromOrder::default(), ChromOrder::Karyotype);
}
