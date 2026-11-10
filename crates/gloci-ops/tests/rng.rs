use gloci_ops::rng::{Stream, seed_or_entropy};

/// The first four words of `Stream::new(42, 0)`. Changing them changes every
/// seeded output of gloci: it is a breaking change for the release notes.
const KNOWN: [u64; 4] = [
    6_424_161_053_832_095_879,
    5_270_208_426_312_333_099,
    9_102_960_255_288_774_902,
    6_786_414_922_848_728_802,
];

#[test]
fn rng_known_answers() {
    let mut s = Stream::new(42, 0);
    let got: Vec<u64> = (0..4).map(|_| s.next_u64()).collect();
    assert_eq!(got, KNOWN);
    assert_ne!(Stream::new(42, 1).next_u64(), KNOWN[0]);
    assert_eq!(Stream::at_word(42, 0, 2).next_u64(), KNOWN[1]);
    assert_eq!(Stream::at_word(42, 0, 4).next_u64(), KNOWN[2]);
}

#[test]
fn rng_below() {
    let mut s = Stream::new(7, 0);
    assert_eq!(s.below(1), 0);
    assert_eq!(s.below(0), 0);
    let _ = s.below(u64::MAX);
    // A bound that forces the rejection path often: 2^63 + 1.
    for _ in 0..1000 {
        assert!(s.below((1 << 63) + 1) < (1 << 63) + 1);
    }
    let mut counts = [0u64; 3];
    let n = 300_000u64;
    for _ in 0..n {
        counts[s.below(3) as usize] += 1;
    }
    let e = n as f64 / 3.0;
    let chi2: f64 = counts.iter().map(|&c| (c as f64 - e).powi(2) / e).sum();
    assert!(chi2 < 13.82, "chi2 = {chi2} for {counts:?}");
}

#[test]
fn rng_unit_and_coin() {
    let mut s = Stream::new(7, 3);
    for _ in 0..10_000 {
        let u = s.unit();
        assert!((0.0..1.0).contains(&u), "{u}");
    }
    let heads = (0..100_000).filter(|_| s.coin()).count() as f64 / 100_000.0;
    assert!((heads - 0.5).abs() < 0.005, "{heads}");
}

#[test]
fn seed_is_kept_when_given() {
    assert_eq!(seed_or_entropy(Some(9)), 9);
}
