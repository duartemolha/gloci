//! Genotype counts taken from the GT bytes of a BCF record, without decoding
//! one sample at a time.
//!
//! A cohort file has thousands of samples per record and millions of
//! records, so the site statistics are a loop over the GT vector that must
//! cost about what reading the bytes costs. [`count`] gives what
//! [`crate::site_stats`] gives for the genotypes of unmasked samples, or
//! `None` for a record the per-sample decoder reports an error on: the
//! caller then runs that decoder, so the error is its error.

use gloci_io::bcf::typed::Ty;

/// What [`count`] counts besides the calls of each allele.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub(crate) struct GtCounts {
    /// The samples counted.
    pub samples: u64,
    /// Samples with no genotype, or with a missing allele in it.
    pub missing: u64,
    pub called_alleles: u64,
    /// Alleles written, missing ones included; two for a sample with no genotype.
    pub data_alleles: u64,
    pub missing_alleles: u64,
    pub haploid: bool,
    pub polyploid: bool,
    /// A genotype of two or more alleles has a `/`.
    pub unphased: bool,
    /// Fully called two-allele genotypes by the sum of their allele indices:
    /// 0, 1, and 2 or more (0/0, 0/1 and 1/1 where the site has two alleles).
    pub diploid: (u64, u64, u64),
}

/// Element `e` (`W` bytes) of an integer vector, sign-extended.
#[inline(always)]
fn value<const W: usize>(e: &[u8]) -> i32 {
    match W {
        1 => i32::from(e[0] as i8),
        2 => i32::from(i16::from_le_bytes([e[0], e[1]])),
        _ => i32::from_le_bytes([e[0], e[1], e[2], e[3]]),
    }
}

/// Count one sample's genotype: the `W`-byte values of `sample`, up to the
/// end-of-vector value. False for a value the decoder rejects (negative, or
/// an allele `calls` has no entry for).
#[inline(always)]
fn count_sample<const W: usize>(sample: &[u8], calls: &mut [u64], c: &mut GtCounts) -> bool {
    let missing = match W {
        1 => i32::from(i8::MIN),
        2 => i32::from(i16::MIN),
        _ => i32::MIN,
    };
    let (mut ploidy, mut dots, mut bars) = (0u64, 0u64, 0u64);
    let mut index_sum = 0usize;
    for (j, e) in sample.as_chunks::<W>().0.iter().enumerate() {
        let x = value::<W>(e);
        if x == missing + 1 {
            break;
        }
        ploidy += 1;
        if x == missing {
            dots += 1;
            continue;
        }
        if x < 0 {
            return false;
        }
        // The low bit says the allele is phased with the one before it.
        bars += u64::from(j > 0 && x & 1 == 1);
        let Some(allele) = ((x >> 1) as usize).checked_sub(1) else {
            dots += 1;
            continue;
        };
        let Some(n) = calls.get_mut(allele) else {
            return false;
        };
        *n += 1;
        index_sum += allele;
    }
    c.samples += 1;
    if ploidy == 0 {
        // No genotype: a missing sample that stands for two missing alleles.
        c.missing += 1;
        c.data_alleles += 2;
        c.missing_alleles += 2;
        return true;
    }
    c.unphased |= ploidy >= 2 && bars != ploidy - 1;
    c.haploid |= ploidy == 1;
    c.polyploid |= ploidy > 2;
    c.data_alleles += ploidy;
    c.missing_alleles += dots;
    c.called_alleles += ploidy - dots;
    if dots > 0 {
        c.missing += 1;
    } else if ploidy == 2 {
        match index_sum {
            0 => c.diploid.0 += 1,
            1 => c.diploid.1 += 1,
            _ => c.diploid.2 += 1,
        }
    }
    true
}

/// Every sample in turn: `n` values of `W` bytes each.
fn count_each<const W: usize>(
    n: usize,
    data: &[u8],
    keep: Option<&[bool]>,
    calls: &mut [u64],
    c: &mut GtCounts,
) -> bool {
    for (i, sample) in data.chunks_exact(n * W).enumerate() {
        if keep.is_some_and(|k| !k[i]) {
            continue;
        }
        if !count_sample::<W>(sample, calls, c) {
            return false;
        }
    }
    true
}

/// [`count_sample`] for each two-value int8 genotype of `samples`.
#[inline(always)]
fn count_pairs(samples: &[u8], calls: &mut [u64], c: &mut GtCounts) -> bool {
    let (pairs, _) = samples.as_chunks::<2>();
    pairs
        .iter()
        .all(|sample| count_sample::<1>(sample, calls, c))
}

/// The low bit of every byte of a word.
const LOW: u64 = 0x0101_0101_0101_0101;
/// The low bit of the first allele of each of a word's four genotypes, and
/// of the second.
const FIRST: u64 = 0x0001_0001_0001_0001;
const SECOND: u64 = FIRST << 8;
/// Words a block holds: a count per byte of a word must stay below 256.
const BLOCK: usize = 255;

/// The sum of the four 16-bit counts of `lanes`, each below 2^14.
#[inline(always)]
fn sum_lanes(lanes: u64) -> u64 {
    lanes.wrapping_mul(FIRST) >> 48
}

/// Two int8 values per sample, every sample counted, at a site with at least
/// two alleles: the shape of a diploid cohort. Eight bytes (four genotypes)
/// are counted at once when each is 2, 3, 4 or 5, which is allele 0 or 1
/// after a `/` or a `|`; any other word goes through [`count_sample`].
fn count_diploid_int8(data: &[u8], calls: &mut [u64], c: &mut GtCounts) -> bool {
    let (mut words, mut alt, mut het) = (0u64, 0u64, 0u64);
    // The bits every counted word has: its second alleles' low bits say `|`.
    let mut common = SECOND;
    let (chunks, rest) = data.as_chunks::<8>();
    for block in chunks.chunks(BLOCK) {
        // Per byte of a word, the words of this block whose allele there is
        // 1; per genotype (16 bits), those whose two alleles differ.
        let (mut alts, mut hets) = (0u64, 0u64);
        for chunk in block {
            let w = u64::from_le_bytes(*chunk);
            // 2 to 5: nothing above the third bit, and one of the second and third.
            if w & !(7 * LOW) == 0 && ((w >> 1) ^ (w >> 2)) & LOW == LOW {
                // The third bit is set for allele 1.
                let alt1 = (w >> 2) & LOW;
                words += 1;
                alts += alt1;
                hets += (alt1 ^ (alt1 >> 8)) & FIRST;
                common &= w;
            } else if !count_pairs(chunk, calls, c) {
                return false;
            }
        }
        alt += sum_lanes((alts & (0xff * FIRST)) + ((alts >> 8) & (0xff * FIRST)));
        het += sum_lanes(hets);
    }
    if !count_pairs(rest, calls, c) {
        return false;
    }
    calls[0] += 8 * words - alt;
    calls[1] += alt;
    c.samples += 4 * words;
    c.called_alleles += 8 * words;
    c.data_alleles += 8 * words;
    c.unphased |= common != SECOND;
    // Each genotype is 0/0, 0/1 or 1/1, so the alleles 1 not in a
    // heterozygote come in pairs.
    let hom_alt = (alt - het) / 2;
    c.diploid.0 += 4 * words - het - hom_alt;
    c.diploid.1 += het;
    c.diploid.2 += hom_alt;
    true
}

/// Count the genotypes of a record whose GT field holds `n` values of the
/// integer type `ty` per sample in `data`, for the samples `keep` selects
/// (`None` = all; else one flag per sample). `calls[a]` gets the calls of
/// allele `a`. `None` when the field is not one of integers, or holds a value
/// the per-sample decoder rejects: a negative one that is neither missing nor
/// the end of the vector, or an allele `calls` has no entry for. `calls` is
/// then partly counted.
pub(crate) fn count(
    ty: Ty,
    n: usize,
    data: &[u8],
    keep: Option<&[bool]>,
    calls: &mut [u64],
) -> Option<GtCounts> {
    let mut c = GtCounts::default();
    let ok = match ty {
        _ if n == 0 => false,
        Ty::Int8 if n == 2 && keep.is_none() && calls.len() >= 2 => {
            count_diploid_int8(data, calls, &mut c)
        }
        Ty::Int8 => count_each::<1>(n, data, keep, calls, &mut c),
        Ty::Int16 => count_each::<2>(n, data, keep, calls, &mut c),
        Ty::Int32 => count_each::<4>(n, data, keep, calls, &mut c),
        Ty::Missing | Ty::Float | Ty::Char => false,
    };
    ok.then_some(c)
}
