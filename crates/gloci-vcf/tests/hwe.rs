use gloci_vcf::hwe_p;

fn lfact(n: u64) -> f64 {
    (2..=n).map(|k| (k as f64).ln()).sum()
}

/// Brute force: exact probability of every het count with the observed
/// allele counts, summing those no larger than the observed one.
fn brute(het: u64, h1: u64, h2: u64) -> f64 {
    let n = het + h1 + h2;
    let homr = h1.min(h2);
    let rare = 2 * homr + het;
    let logp = |h: u64| {
        let r = (rare - h) / 2;
        let c = n - h - r;
        lfact(n) - lfact(r) - lfact(h) - lfact(c)
            + (h as f64) * 2f64.ln()
            + lfact(rare)
            + lfact(2 * n - rare)
            - lfact(2 * n)
    };
    let obs = logp(het);
    let mut total = 0.0;
    let mut p = 0.0;
    for h in (rare % 2..=rare).step_by(2) {
        if h > n || (rare - h) / 2 + h > n {
            continue;
        }
        let v = logp(h).exp();
        total += v;
        if logp(h) <= obs + 1e-12 {
            p += v;
        }
    }
    (p / total).min(1.0)
}

#[test]
fn hwe() {
    assert!((hwe_p(0, 1, 1) - 1.0 / 3.0).abs() < 1e-12);
    assert_eq!(hwe_p(2, 0, 0), 1.0);
    assert_eq!(hwe_p(0, 0, 0), 1.0);
    assert!((hwe_p(57, 14, 50) - brute(57, 14, 50)).abs() < 1e-10);
}

#[test]
fn hwe_matches_brute_force_on_a_grid() {
    for het in 0..12 {
        for h1 in 0..8 {
            for h2 in 0..8 {
                let (a, b) = (hwe_p(het, h1, h2), brute(het, h1, h2));
                assert!((a - b).abs() < 1e-10, "({het},{h1},{h2}): {a} vs {b}");
                assert!((0.0..=1.0).contains(&a));
            }
        }
    }
}
