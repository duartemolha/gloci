//! Hardy-Weinberg exact test (Wigginton, Cutler & Abecasis 2005).

/// Two-sided exact p-value: the total probability of every heterozygote
/// count that is no more likely than the observed one, capped at 1.
/// `hom1` and `hom2` may be given in either order.
pub fn hwe_p(het: u64, hom1: u64, hom2: u64) -> f64 {
    let n = het + hom1 + hom2;
    if n == 0 {
        return 1.0;
    }
    let homr = hom1.min(hom2);
    let rare = 2 * homr + het; // copies of the rarer allele
    let mut probs = vec![0.0f64; rare as usize + 1];

    // Start at the most likely het count and recurse outwards.
    let mut mid = (rare as u128 * (2 * n - rare) as u128 / (2 * n) as u128) as u64;
    if mid % 2 != rare % 2 {
        mid += 1;
    }
    probs[mid as usize] = 1.0;
    let mut sum = 1.0;

    let (mut h, mut r, mut c) = (mid, (rare - mid) / 2, n - mid - (rare - mid) / 2);
    while h >= 2 {
        let p = probs[h as usize] * (h * (h - 1)) as f64 / (4.0 * (r + 1) as f64 * (c + 1) as f64);
        probs[h as usize - 2] = p;
        sum += p;
        h -= 2;
        r += 1;
        c += 1;
    }
    let (mut h, mut r, mut c) = (mid, (rare - mid) / 2, n - mid - (rare - mid) / 2);
    while h + 2 <= rare {
        let p = probs[h as usize] * 4.0 * r as f64 * c as f64 / ((h + 2) * (h + 1)) as f64;
        probs[h as usize + 2] = p;
        sum += p;
        h += 2;
        r -= 1;
        c -= 1;
    }

    let obs = probs[het as usize];
    let tail: f64 = probs.iter().filter(|&&p| p <= obs).sum();
    (tail / sum).min(1.0)
}
