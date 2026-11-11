//! The set of places a shuffled record may go, with exact uniform draws (D152).
//!
//! The allowed set is the `--include` regions (else the whole genome) minus
//! the `--exclude` regions, merged per chromosome into disjoint intervals
//! `[a, b)`. A record of length `L` may start at any `s` with `a <= s <= b - L`
//! inside one interval, so an interval of `n` bases holds `n - L + 1`
//! placements (none when `n < L`). Counting and drawing are exact:
//!
//! - intervals sit in one array sorted by (chromosome, length, start), and a
//!   second view lists their indices by (length, index); each view has prefix
//!   sums of `length + 1` (u128), so the placements of every interval at least
//!   `L` long are `sum(length + 1) - L * count`, found with one binary search;
//! - the *r*-th placement is a second binary search, because the placements an
//!   interval holds grow with its length;
//! - chromosomes are also listed by their longest interval, for `--chrom-first`.
//!
//! Memory is about 52 bytes per allowed interval.

use crate::rng::Stream;
use gloci_core::{ChromId, ChromTable, Genome, Interval};

const NONE: u32 = u32::MAX;

/// Where records may be placed.
pub struct Allowed {
    /// Genome index of each chromosome id (`NONE`: not in the genome).
    gidx: Vec<u32>,
    /// Chromosome id of each genome index.
    ids: Vec<Option<ChromId>>,
    /// `offs[g]..offs[g + 1]` are chromosome `g`'s intervals in `ivs`.
    offs: Vec<u32>,
    /// `(start, length)`, by (chromosome, length, start).
    ivs: Vec<(u64, u64)>,
    /// `pre[k]` is the sum of `length + 1` over `ivs[..k]`.
    pre: Vec<u128>,
    /// Indices into `ivs`, by (length, index).
    order: Vec<u32>,
    /// `gpre[k]` is the sum of `length + 1` over the first `k` of `order`.
    gpre: Vec<u128>,
    /// `(longest interval, genome index)` of the chromosomes that have one,
    /// by (longest, genome index).
    longest: Vec<(u64, u32)>,
}

/// Regions as `(genome index, start, end)`, non-empty, sorted and merged
/// (overlapping and bookended regions join).
fn merged(regions: impl Iterator<Item = (u32, u64, u64)>) -> Vec<(u32, u64, u64)> {
    let mut v: Vec<(u32, u64, u64)> = regions.filter(|&(_, s, e)| s < e).collect();
    v.sort_unstable();
    let mut out: Vec<(u32, u64, u64)> = Vec::with_capacity(v.len());
    for (g, s, e) in v {
        match out.last_mut() {
            Some(last) if last.0 == g && s <= last.2 => last.2 = last.2.max(e),
            _ => out.push((g, s, e)),
        }
    }
    out
}

impl Allowed {
    /// The allowed set of `genome`: `include` (else the whole genome) less
    /// `exclude`. Chromosomes are addressed by their ids in `chroms`; regions on
    /// other chromosomes are ignored, and a region is cut to its chromosome.
    /// A genome chromosome without a length, or missing from `chroms`, has no
    /// allowed intervals. Zero-length regions are ignored.
    pub fn build(
        genome: &Genome,
        chroms: &ChromTable,
        include: Option<&[(ChromId, Interval)]>,
        exclude: &[(ChromId, Interval)],
    ) -> Allowed {
        let n = genome.len();
        let mut gidx = vec![NONE; chroms.len()];
        let mut ids = Vec::with_capacity(n);
        for g in 0..n {
            let id = chroms.get(genome.name(g));
            if let Some(id) = id {
                gidx[id.0 as usize] = g as u32;
            }
            ids.push(id);
        }
        let lens: Vec<u64> = (0..n).map(|g| genome.length(g).unwrap_or(0)).collect();
        let clip = |&(c, iv): &(ChromId, Interval)| -> Option<(u32, u64, u64)> {
            let g = *gidx.get(c.0 as usize)?;
            (g != NONE).then(|| {
                (
                    g,
                    iv.start.min(lens[g as usize]),
                    iv.end.min(lens[g as usize]),
                )
            })
        };
        let base = match include {
            Some(inc) => merged(inc.iter().filter_map(clip)),
            None => merged((0..n).map(|g| (g as u32, 0, lens[g]))),
        };
        let excl = merged(exclude.iter().filter_map(clip));
        // Allowed segments: each base interval less the exclusions on its chromosome.
        let mut segs: Vec<(u32, u64, u64)> = Vec::with_capacity(base.len());
        let mut e = 0usize;
        for &(g, mut s, end) in &base {
            while e < excl.len() && (excl[e].0, excl[e].2) <= (g, s) {
                e += 1;
            }
            let mut k = e;
            while s < end {
                match excl.get(k) {
                    Some(&(eg, es, ee)) if eg == g && es < end => {
                        if es > s {
                            segs.push((g, s, es));
                        }
                        s = s.max(ee);
                        k += 1;
                    }
                    _ => {
                        segs.push((g, s, end));
                        break;
                    }
                }
            }
        }
        let mut keyed: Vec<(u32, u64, u64)> =
            segs.into_iter().map(|(g, s, e)| (g, e - s, s)).collect();
        keyed.sort_unstable();
        let mut offs = vec![0u32; n + 1];
        for &(g, _, _) in &keyed {
            offs[g as usize + 1] += 1;
        }
        for g in 0..n {
            offs[g + 1] += offs[g];
        }
        let ivs: Vec<(u64, u64)> = keyed.iter().map(|&(_, l, s)| (s, l)).collect();
        let mut pre = Vec::with_capacity(ivs.len() + 1);
        let mut sum = 0u128;
        pre.push(0);
        for &(_, l) in &ivs {
            sum += u128::from(l) + 1;
            pre.push(sum);
        }
        let mut order: Vec<u32> = (0..ivs.len() as u32).collect();
        order.sort_unstable_by_key(|&i| (ivs[i as usize].1, i));
        let mut gpre = Vec::with_capacity(ivs.len() + 1);
        let mut sum = 0u128;
        gpre.push(0);
        for &i in &order {
            sum += u128::from(ivs[i as usize].1) + 1;
            gpre.push(sum);
        }
        let mut longest: Vec<(u64, u32)> = (0..n)
            .filter(|&g| offs[g] < offs[g + 1])
            // The last interval of a chromosome is its longest (sorted by length).
            .map(|g| (ivs[offs[g + 1] as usize - 1].1, g as u32))
            .collect();
        longest.sort_unstable();
        Allowed {
            gidx,
            ids,
            offs,
            ivs,
            pre,
            order,
            gpre,
            longest,
        }
    }

    /// Whether `c` is a chromosome of the genome.
    pub fn knows(&self, c: ChromId) -> bool {
        self.gidx.get(c.0 as usize).is_some_and(|&g| g != NONE)
    }

    /// The number of allowed intervals.
    pub fn intervals(&self) -> usize {
        self.ivs.len()
    }

    fn genome_index(&self, c: ChromId) -> Option<usize> {
        self.gidx
            .get(c.0 as usize)
            .filter(|&&g| g != NONE)
            .map(|&g| g as usize)
    }

    /// How many `(chromosome, start)` places a record of `len` bases has: on
    /// chromosome `chrom` only, or anywhere.
    pub fn placements(&self, len: u64, chrom: Option<ChromId>) -> u128 {
        match chrom {
            None => {
                let from = self
                    .order
                    .partition_point(|&i| self.ivs[i as usize].1 < len);
                let count = (self.order.len() - from) as u128;
                self.gpre[self.order.len()] - self.gpre[from] - u128::from(len) * count
            }
            Some(c) => match self.genome_index(c) {
                None => 0,
                Some(g) => {
                    let (lo, hi) = (self.offs[g] as usize, self.offs[g + 1] as usize);
                    let from = lo + self.ivs[lo..hi].partition_point(|iv| iv.1 < len);
                    self.pre[hi] - self.pre[from] - u128::from(len) * (hi - from) as u128
                }
            },
        }
    }

    /// Draw one place for a record of `len` bases uniformly among its
    /// placements (all, or those on `chrom`); `None` when it has none.
    pub fn place(
        &self,
        len: u64,
        chrom: Option<ChromId>,
        rng: &mut Stream,
    ) -> Option<(ChromId, u64)> {
        match chrom {
            None => {
                let total = self.placements(len, None);
                if total == 0 {
                    return None;
                }
                let r = below(rng, total);
                let from = self
                    .order
                    .partition_point(|&i| self.ivs[i as usize].1 < len);
                // Placements of the first k intervals of `order[from..]`.
                let cum = |k: usize| {
                    self.gpre[k] - self.gpre[from] - u128::from(len) * (k - from) as u128
                };
                let k = find(from, self.order.len(), r, cum);
                let i = self.order[k] as usize;
                Some(self.at(i, r - cum(k)))
            }
            Some(c) => {
                let g = self.genome_index(c)?;
                self.place_in(g, len, rng)
            }
        }
    }

    /// Draw a chromosome uniformly among those with a placement for `len`,
    /// then a place on it.
    pub fn place_chrom_first(&self, len: u64, rng: &mut Stream) -> Option<(ChromId, u64)> {
        let from = self.longest.partition_point(|&(l, _)| l < len);
        let n = self.longest.len() - from;
        if n == 0 {
            return None;
        }
        let g = self.longest[from + rng.below(n as u64) as usize].1;
        self.place_in(g as usize, len, rng)
    }

    fn place_in(&self, g: usize, len: u64, rng: &mut Stream) -> Option<(ChromId, u64)> {
        let (lo, hi) = (self.offs[g] as usize, self.offs[g + 1] as usize);
        let from = lo + self.ivs[lo..hi].partition_point(|iv| iv.1 < len);
        let cum = |k: usize| self.pre[k] - self.pre[from] - u128::from(len) * (k - from) as u128;
        let total = cum(hi);
        if total == 0 {
            return None;
        }
        let r = below(rng, total);
        let k = find(from, hi, r, cum);
        Some(self.at(k, r - cum(k)))
    }

    /// The place at `offset` bases into interval `i`.
    fn at(&self, i: usize, offset: u128) -> (ChromId, u64) {
        let g = self.offs.partition_point(|&o| o as usize <= i) - 1;
        let id = self.ids[g].expect("a chromosome with intervals has an id");
        (id, self.ivs[i].0 + offset as u64)
    }
}

/// The smallest `k` in `from..to` whose interval holds placement `r`, that is
/// `cum(k) <= r < cum(k + 1)`.
fn find(from: usize, to: usize, r: u128, cum: impl Fn(usize) -> u128) -> usize {
    let (mut lo, mut hi) = (from, to);
    while hi - lo > 1 {
        let mid = lo + (hi - lo) / 2;
        if cum(mid) <= r {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo
}

/// A uniform integer in `0..n` (`n > 0`), beyond 64 bits too.
fn below(rng: &mut Stream, n: u128) -> u128 {
    if let Ok(n) = u64::try_from(n) {
        return u128::from(rng.below(n));
    }
    let mask = u128::MAX >> n.leading_zeros();
    loop {
        let x = ((u128::from(rng.next_u64()) << 64) | u128::from(rng.next_u64())) & mask;
        if x < n {
            return x;
        }
    }
}

/// Merged regions per chromosome, to test whether an interval touches one.
pub struct Regions {
    by_chrom: Vec<Vec<(u64, u64)>>,
}

impl Regions {
    pub fn new(regions: &[(ChromId, Interval)], chroms: &ChromTable) -> Regions {
        let mut by_chrom = vec![Vec::new(); chroms.len()];
        for (g, s, e) in merged(regions.iter().map(|(c, iv)| (c.0, iv.start, iv.end))) {
            by_chrom[g as usize].push((s, e));
        }
        Regions { by_chrom }
    }

    pub fn is_empty(&self) -> bool {
        self.by_chrom.iter().all(Vec::is_empty)
    }

    /// Whether `[start, end)` overlaps a region. An interval that only touches
    /// a region's edge does not overlap it, nor does a point on an edge.
    pub fn overlaps(&self, c: ChromId, start: u64, end: u64) -> bool {
        let Some(v) = self.by_chrom.get(c.0 as usize) else {
            return false;
        };
        // Regions are disjoint and sorted, so the last one starting before
        // `end` is the one reaching furthest.
        let i = v.partition_point(|&(s, _)| s < end);
        i > 0 && v[i - 1].1 > start
    }
}
