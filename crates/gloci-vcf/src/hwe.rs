//! Hardy-Weinberg exact test (Wigginton, Cutler & Abecasis 2005).

use std::cell::RefCell;

thread_local! {
    static SCRATCH: RefCell<Vec<f64>> = const { RefCell::new(Vec::new()) };
}

/// The exact test's three p-values for one genotype table.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HweTails {
    /// Two-sided: the total probability of every heterozygote count that is
    /// no more likely than the observed one, capped at 1.
    pub two_sided: f64,
    /// One-sided: the probability of at most the observed heterozygotes.
    pub het_deficit: f64,
    /// One-sided: the probability of at least the observed heterozygotes.
    pub het_excess: f64,
}

/// Two-sided exact p-value: the total probability of every heterozygote
/// count that is no more likely than the observed one, capped at 1.
/// `hom1` and `hom2` may be given in either order.
pub fn hwe_p(het: u64, hom1: u64, hom2: u64) -> f64 {
    hwe_tails(het, hom1, hom2).two_sided
}

/// All three exact p-values (two-sided and both one-sided tails).
/// `hom1` and `hom2` may be given in either order.
pub fn hwe_tails(het: u64, hom1: u64, hom2: u64) -> HweTails {
    SCRATCH.with(|s| hwe_with(&mut s.borrow_mut(), het, hom1, hom2))
}

fn hwe_with(probs: &mut Vec<f64>, het: u64, hom1: u64, hom2: u64) -> HweTails {
    let n = het + hom1 + hom2;
    if n == 0 {
        return HweTails {
            two_sided: 1.0,
            het_deficit: 1.0,
            het_excess: 1.0,
        };
    }
    let homr = hom1.min(hom2);
    let rare = 2 * homr + het; // copies of the rarer allele
    probs.clear();
    probs.resize(rare as usize + 1, 0.0);

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
    // A relative tolerance keeps exact ties (equal in theory, unequal in the
    // last bits) on the same side, as plink2 does.
    let limit = obs * (1.0 + 1e-7);
    let tail: f64 = probs.iter().filter(|&&p| p <= limit).sum();
    // Heterozygote counts only take the parity of `rare`, so entries of the
    // other parity are zero and add nothing to either one-sided sum.
    let deficit: f64 = probs[..=het as usize].iter().sum();
    let excess: f64 = probs[het as usize..].iter().sum();
    HweTails {
        two_sided: (tail / sum).min(1.0),
        het_deficit: (deficit / sum).min(1.0),
        het_excess: (excess / sum).min(1.0),
    }
}
