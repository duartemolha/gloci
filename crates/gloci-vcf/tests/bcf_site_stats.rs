//! `BcfGenotypes::site_stats` counts a record's genotypes from the GT bytes.
//! It must give what `site_stats` gives sample by sample, for every way BCF
//! stores a genotype, and what the plain reference in the testkit counts.

use gloci_io::BcfDicts;
use gloci_testkit::naive::{NaiveGt, naive_site_counts};
use gloci_vcf::{GenotypeSource, MaskRules, SiteStats, bcf_genotypes, site_stats};
use proptest::prelude::*;

const GT: u8 = 1;
const DP: u8 = 2;
const GQ: u8 = 3;
const FT: u8 = 4;

const INT8: u8 = 1;
const INT16: u8 = 2;
const INT32: u8 = 3;
const FLOAT: u8 = 5;
const CHAR: u8 = 7;

fn dicts() -> BcfDicts {
    BcfDicts {
        contigs: vec![b"1".to_vec()],
        strings: ["PASS", "GT", "DP", "GQ", "FT"]
            .iter()
            .map(|s| s.as_bytes().to_vec())
            .collect(),
    }
}

/// One FORMAT field: its key, type code, values per sample and the bytes of
/// all samples.
struct Field {
    key: u8,
    ty: u8,
    n: usize,
    data: Vec<u8>,
}

/// A framed BCF record of `n_samples` samples with the FORMAT fields `fields`.
fn record(n_samples: usize, fields: &[Field]) -> Vec<u8> {
    let mut shared = vec![0u8; 24];
    shared[20..24]
        .copy_from_slice(&(((fields.len() as u32) << 24) | n_samples as u32).to_le_bytes());
    shared.extend_from_slice(&[0x07, 0x07]);
    let mut indiv = Vec::new();
    for f in fields {
        assert!(f.n < 15, "the test encoder writes short vectors only");
        indiv.extend_from_slice(&[0x11, f.key, ((f.n as u8) << 4) | f.ty]);
        indiv.extend_from_slice(&f.data);
    }
    let mut raw = Vec::new();
    raw.extend_from_slice(&(shared.len() as u32).to_le_bytes());
    raw.extend_from_slice(&(indiv.len() as u32).to_le_bytes());
    raw.extend_from_slice(&shared);
    raw.extend_from_slice(&indiv);
    raw
}

/// The integer type's missing value; one more is its end-of-vector value.
fn missing(ty: u8) -> i32 {
    match ty {
        INT8 => i32::from(i8::MIN),
        INT16 => i32::from(i16::MIN),
        _ => i32::MIN,
    }
}

fn ints(ty: u8, values: &[i32]) -> Vec<u8> {
    let mut out = Vec::new();
    for &v in values {
        match ty {
            INT8 => out.push(v as i8 as u8),
            INT16 => out.extend_from_slice(&(v as i16).to_le_bytes()),
            _ => out.extend_from_slice(&v.to_le_bytes()),
        }
    }
    out
}

/// A record whose only FORMAT field is GT: `n` values of `ty` per sample.
fn gt_record(ty: u8, n: usize, values: &[i32]) -> Vec<u8> {
    let n_samples = values.len().checked_div(n).unwrap_or(0);
    record(
        n_samples,
        &[Field {
            key: GT,
            ty,
            n,
            data: ints(ty, values),
        }],
    )
}

/// The value BCF stores for allele `a` (`None` is `.`) after a `|` or a `/`.
fn val(a: Option<u32>, bar: bool) -> i32 {
    ((a.map_or(0, |a| a as i32 + 1)) << 1) | i32::from(bar)
}

const NO_MASK: MaskRules = MaskRules {
    min_gq: None,
    min_dp: None,
    max_dp: None,
    mask_filtered: false,
};

/// Both ways of counting `raw`: from the bytes, and sample by sample.
#[allow(clippy::type_complexity)]
fn both(
    raw: &[u8],
    keep: Option<&[bool]>,
    mask: &MaskRules,
    n_alleles: usize,
) -> (
    (Result<SiteStats, String>, Vec<bool>),
    (Result<SiteStats, String>, Vec<bool>),
) {
    let d = dicts();
    let b = bcf_genotypes(raw, &d).unwrap();
    let (mut m1, mut m2) = (vec![true; 3], vec![true; 3]);
    let fast = b.site_stats(keep, mask, n_alleles, &mut m1);
    let slow = site_stats(&b as &dyn GenotypeSource, keep, mask, n_alleles, &mut m2);
    ((fast, m1), (slow, m2))
}

/// The counts of `raw`, which both ways must agree on.
fn counts(raw: &[u8], keep: Option<&[bool]>, n_alleles: usize) -> Result<SiteStats, String> {
    let ((fast, m1), (slow, m2)) = both(raw, keep, &NO_MASK, n_alleles);
    assert_eq!(fast, slow);
    assert_eq!(m1, m2);
    fast
}

#[test]
fn counts_int8_diploid_phased() {
    // 0|0 0|1 1|0 1|1 0|0, the shape of a phased cohort file.
    let raw = gt_record(INT8, 2, &[2, 3, 2, 5, 4, 3, 4, 5, 2, 3]);
    let s = counts(&raw, None, 2).unwrap();
    assert_eq!(s.allele_counts, [6, 4]);
    assert_eq!((s.called_alleles, s.samples, s.missing), (10, 5, 0));
    assert_eq!((s.data_alleles, s.missing_alleles), (10, 0));
    assert_eq!(s.hwe_counts, Some((2, 2, 1)));
    assert!(s.all_phased && !s.haploid && !s.polyploid);
    assert_eq!((s.dp_sum, s.dp_n, s.dp_sumsq, s.filtered), (0, 0, 0, 0));
}

#[test]
fn counts_many_samples_cross_every_word_boundary() {
    // From no sample to several eight-byte words and a remainder, with one
    // odd genotype moved through every position.
    for n_samples in 0..40 {
        for odd in 0..=n_samples {
            let mut v = Vec::new();
            for i in 0..n_samples {
                if i == odd {
                    v.extend_from_slice(&[val(None, false), val(Some(1), false)]); // ./1
                } else {
                    v.extend_from_slice(&[val(Some((i % 2) as u32), false), val(Some(1), true)]);
                }
            }
            let raw = gt_record(INT8, 2, &v);
            let s = counts(&raw, None, 2).unwrap();
            let has_odd = odd < n_samples;
            assert_eq!(s.samples, n_samples as u64);
            assert_eq!(s.missing, u64::from(has_odd));
            assert_eq!(s.called_alleles, 2 * n_samples as u64 - u64::from(has_odd));
            assert_eq!(s.all_phased, !has_odd, "{n_samples} samples, odd at {odd}");
        }
    }
}

#[test]
fn counts_records_of_thousands_of_samples() {
    // Whole words are counted in blocks of 255 (1,020 samples): the counts
    // must hold at a block's last word, past it and over several blocks,
    // when every genotype adds to the same counters.
    // The genotype, its two values, its calls of each allele, and which of
    // 0/0, 0/1 and 1/1 it is.
    let patterns = [
        ("1|1", [4i32, 5], [0u64, 2], (0u64, 0u64, 1u64)),
        ("0|1", [2, 5], [1, 1], (0, 1, 0)),
        ("1|0", [4, 3], [1, 1], (0, 1, 0)),
        ("0|0", [2, 3], [2, 0], (1, 0, 0)),
    ];
    for n_samples in [1016usize, 1019, 1020, 1021, 1024, 1028, 2040, 2044, 5003] {
        for (name, pair, calls, kind) in patterns {
            let n = n_samples as u64;
            let v = pair.repeat(n_samples);
            let s = counts(&gt_record(INT8, 2, &v), None, 2).unwrap();
            assert_eq!(
                s.allele_counts,
                [calls[0] * n, calls[1] * n],
                "{n} x {name}"
            );
            assert_eq!(
                s.hwe_counts,
                Some((kind.0 * n, kind.1 * n, kind.2 * n)),
                "{n} x {name}"
            );
            assert_eq!((s.samples, s.missing, s.called_alleles), (n, 0, 2 * n));
            assert!(s.all_phased);
            // One unphased, half-missing genotype anywhere leaves the rest counted.
            for at in [0, 3, 1015, n_samples - 5, n_samples - 1] {
                let mut v = v.clone();
                v[2 * at] = 0;
                v[2 * at + 1] = 4;
                let s = counts(&gt_record(INT8, 2, &v), None, 2).unwrap();
                assert_eq!(
                    s.allele_counts,
                    [calls[0] * (n - 1), calls[1] * (n - 1) + 1],
                    "{n} x {name}, ./1 at {at}"
                );
                assert_eq!((s.missing, s.all_phased), (1, false));
            }
        }
    }
}

#[test]
fn counts_int16_and_int32_vectors() {
    // 0/70 70|0 1/1: allele 70 needs int16.
    let v = [
        val(Some(0), false),
        val(Some(70), false),
        val(Some(70), false),
        val(Some(0), true),
        val(Some(1), false),
        val(Some(1), false),
    ];
    for ty in [INT16, INT32] {
        let s = counts(&gt_record(ty, 2, &v), None, 71).unwrap();
        assert_eq!((s.allele_counts[0], s.allele_counts[1]), (2, 2));
        assert_eq!(s.allele_counts[70], 2);
        assert_eq!((s.called_alleles, s.missing), (6, 0));
        assert_eq!(s.hwe_counts, None);
        assert!(!s.all_phased);
    }
    // The same genotypes at each width count the same.
    let v = [2, 5, 4, 5, 0, 1, 2, 2];
    let a = counts(&gt_record(INT8, 2, &v), None, 2).unwrap();
    assert_eq!(a, counts(&gt_record(INT16, 2, &v), None, 2).unwrap());
    assert_eq!(a, counts(&gt_record(INT32, 2, &v), None, 2).unwrap());
}

#[test]
fn counts_missing_alleles_and_the_missing_value() {
    for ty in [INT8, INT16, INT32] {
        let m = missing(ty);
        // ./. .|. ./1 0|. and a genotype of the integer missing value twice.
        let v = [0, 0, 0, 1, 0, 4, 2, 1, m, m];
        let s = counts(&gt_record(ty, 2, &v), None, 2).unwrap();
        assert_eq!(s.allele_counts, [1, 1]);
        assert_eq!((s.samples, s.missing), (5, 5));
        assert_eq!((s.data_alleles, s.missing_alleles), (10, 8));
        assert_eq!(s.hwe_counts, Some((0, 0, 0)));
        // `./.`, `./1` and the missing value have no `|`.
        assert!(!s.all_phased);
        // `.|.` alone is phased, called or not.
        assert!(
            counts(&gt_record(ty, 2, &[0, 1, 2, 1]), None, 2)
                .unwrap()
                .all_phased
        );
        // The missing value after an allele is a `.` with no `|` before it.
        let s = counts(&gt_record(ty, 2, &[2, m]), None, 2).unwrap();
        assert_eq!((s.missing, s.missing_alleles, s.all_phased), (1, 1, false));
    }
}

#[test]
fn counts_end_of_vector_mixed_ploidy() {
    for ty in [INT8, INT16, INT32] {
        let e = missing(ty) + 1;
        // 0/1/1, 1 (haploid), no GT at all, 0|1 padded to three.
        let v = [2, 4, 4, 4, e, e, e, e, e, 2, 5, e];
        let s = counts(&gt_record(ty, 3, &v), None, 2).unwrap();
        assert_eq!(s.allele_counts, [2, 4]);
        assert_eq!((s.samples, s.missing), (4, 1));
        // Three, one, two for the sample without GT, two.
        assert_eq!((s.data_alleles, s.missing_alleles), (8, 2));
        assert!(s.haploid && s.polyploid && !s.all_phased);
        assert_eq!(s.hwe_counts, Some((0, 1, 0)));
        // What follows the end of a vector is not read, whatever it is: a
        // value that would be an error, or one that would be an allele.
        let s = counts(&gt_record(ty, 3, &[2, e, -9, 2, 3, e]), None, 2).unwrap();
        assert_eq!(
            (s.allele_counts[0], s.haploid, s.all_phased),
            (3, true, true)
        );
        let s = counts(&gt_record(ty, 3, &[e, 4, 5, 2, e, 4]), None, 2).unwrap();
        assert_eq!(s.allele_counts, [1, 0]);
        assert_eq!((s.missing, s.data_alleles, s.missing_alleles), (1, 3, 2));
    }
}

#[test]
fn counts_unphased_and_phased() {
    let phased = [2, 3, 4, 5, 2, 5, 4, 3, 2, 3];
    assert!(
        counts(&gt_record(INT8, 2, &phased), None, 2)
            .unwrap()
            .all_phased
    );
    // One `/` anywhere, in a word of plain genotypes or in the remainder.
    for at in 0..5 {
        let mut v = phased;
        v[2 * at + 1] &= !1;
        let s = counts(&gt_record(INT8, 2, &v), None, 2).unwrap();
        assert!(!s.all_phased, "sample {at}");
        assert_eq!(s.allele_counts, [6, 4]);
    }
    // The first allele's low bit says nothing.
    let s = counts(&gt_record(INT8, 2, &[3, 3, 5, 5]), None, 2).unwrap();
    assert!(s.all_phased);
    // Three alleles need both separators.
    assert!(
        counts(&gt_record(INT8, 3, &[2, 3, 5]), None, 2)
            .unwrap()
            .all_phased
    );
    assert!(
        !counts(&gt_record(INT8, 3, &[2, 3, 4]), None, 2)
            .unwrap()
            .all_phased
    );
}

#[test]
fn counts_haploid_records() {
    for ty in [INT8, INT16, INT32] {
        // 0 1 . 1
        let s = counts(&gt_record(ty, 1, &[2, 4, 0, 5]), None, 2).unwrap();
        assert_eq!(s.allele_counts, [1, 2]);
        assert_eq!((s.samples, s.missing, s.data_alleles), (4, 1, 4));
        assert!(s.haploid && !s.polyploid && s.all_phased);
        assert_eq!(s.hwe_counts, Some((0, 0, 0)));
    }
}

#[test]
fn counts_more_than_two_alleles() {
    // 0|2 2|1 1|0 3|3 0|0 of a site with four alleles.
    let v = [2, 7, 6, 5, 4, 3, 8, 9, 2, 3];
    let s = counts(&gt_record(INT8, 2, &v), None, 4).unwrap();
    assert_eq!(s.allele_counts, [4, 2, 2, 2]);
    assert_eq!((s.called_alleles, s.missing), (10, 0));
    assert_eq!(s.hwe_counts, None);
    assert!(s.all_phased);
    // A site with no ALT has the reference allele alone.
    let s = counts(&gt_record(INT8, 2, &[2, 3, 2, 2]), None, 1).unwrap();
    assert_eq!(
        (s.allele_counts.as_slice(), s.all_phased),
        (&[4][..], false)
    );
}

#[test]
fn counts_zero_samples() {
    for n_alleles in [1, 2, 3] {
        let s = counts(&gt_record(INT8, 2, &[]), None, n_alleles).unwrap();
        assert_eq!(s.allele_counts, vec![0; n_alleles]);
        assert_eq!((s.samples, s.missing, s.called_alleles), (0, 0, 0));
        assert_eq!(s.hwe_counts, (n_alleles == 2).then_some((0, 0, 0)));
        assert!(s.all_phased);
    }
    // No FORMAT field at all, and a GT of no values.
    assert_eq!(counts(&record(0, &[]), None, 2).unwrap().samples, 0);
    let s = counts(&record(3, &[]), None, 2).unwrap();
    assert_eq!(
        (s.samples, s.missing, s.data_alleles, s.missing_alleles),
        (3, 3, 6, 6)
    );
    let raw = record(
        2,
        &[Field {
            key: GT,
            ty: INT8,
            n: 0,
            data: vec![],
        }],
    );
    assert_eq!(counts(&raw, None, 2).unwrap().missing, 2);
}

#[test]
fn counts_sample_selection() {
    // 0|0 1|1 0|1 ./. 1|1 0/0
    let v = [2, 3, 4, 5, 2, 5, 0, 0, 4, 5, 2, 2];
    let raw = gt_record(INT8, 2, &v);
    let keep = [true, false, true, false, true, false];
    let s = counts(&raw, Some(&keep), 2).unwrap();
    assert_eq!(s.allele_counts, [3, 3]);
    assert_eq!((s.samples, s.missing), (3, 0));
    assert_eq!(s.hwe_counts, Some((1, 1, 1)));
    assert!(s.all_phased);
    // Nobody selected: nothing is read.
    let s = counts(&raw, Some(&[false; 6]), 2).unwrap();
    assert_eq!((s.samples, s.called_alleles), (0, 0));
    // A selection of the wrong length is the error it always was.
    let e = counts(&raw, Some(&[true; 5]), 2).unwrap_err();
    assert_eq!(
        e,
        "the sample selection has 5 entries; the record has 6 samples"
    );
}

#[test]
fn counts_errors_are_those_of_the_first_bad_sample() {
    // An allele the site does not have.
    let e = counts(&gt_record(INT8, 2, &[2, 3, 2, 7]), None, 2).unwrap_err();
    assert_eq!(e, "genotype allele 2 but the site lists 2 alleles");
    // A negative value that is neither missing nor the end of the vector.
    let e = counts(&gt_record(INT8, 2, &[2, 3, -3, 2]), None, 2).unwrap_err();
    assert_eq!(e, "GT value -3 is negative");
    // The earlier sample's error wins, whichever kind it is.
    let e = counts(&gt_record(INT8, 2, &[2, 9, -3, 2]), None, 2).unwrap_err();
    assert_eq!(e, "genotype allele 3 but the site lists 2 alleles");
    let e = counts(&gt_record(INT8, 2, &[-3, 2, 2, 9]), None, 2).unwrap_err();
    assert_eq!(e, "GT value -3 is negative");
    // Within a sample, a value that cannot be read comes before an allele out of range.
    let e = counts(&gt_record(INT8, 2, &[9, -3]), None, 2).unwrap_err();
    assert_eq!(e, "GT value -3 is negative");
    // A sample left out is not read, so its bad values are no error.
    let raw = gt_record(INT8, 2, &[2, 3, -3, 9]);
    assert!(counts(&raw, Some(&[true, false]), 2).is_ok());
    // The same errors at the wider types, and in a record long enough for whole words.
    for ty in [INT8, INT16, INT32] {
        let mut v = [2, 3].repeat(21);
        v[29] = 8;
        let e = counts(&gt_record(ty, 2, &v), None, 3).unwrap_err();
        assert_eq!(e, "genotype allele 3 but the site lists 3 alleles");
        v[11] = -2;
        let e = counts(&gt_record(ty, 2, &v), None, 3).unwrap_err();
        assert_eq!(e, "GT value -2 is negative");
    }
    // A GT that is not stored as integers.
    let raw = record(
        1,
        &[Field {
            key: GT,
            ty: CHAR,
            n: 2,
            data: b"01".to_vec(),
        }],
    );
    assert_eq!(
        counts(&raw, None, 2).unwrap_err(),
        "FORMAT/GT is not an integer field"
    );
}

#[test]
fn counts_with_depths_and_mask_rules() {
    // 0|1:7:30 1|1:.:5 0|0:12:. 1|0:3:40, FT failing for the last.
    let fields = || {
        vec![
            Field {
                key: GT,
                ty: INT8,
                n: 2,
                data: ints(INT8, &[2, 5, 4, 5, 2, 3, 4, 3]),
            },
            Field {
                key: DP,
                ty: INT16,
                n: 1,
                data: ints(INT16, &[7, missing(INT16), 300, 3]),
            },
            Field {
                key: GQ,
                ty: INT8,
                n: 1,
                data: ints(INT8, &[30, 5, missing(INT8), 40]),
            },
            Field {
                key: FT,
                ty: CHAR,
                n: 4,
                data: b"PASSPASS.\0\0\0LowQ".to_vec(),
            },
        ]
    };
    let raw = record(4, &fields());
    let s = counts(&raw, None, 2).unwrap();
    assert_eq!(s.allele_counts, [4, 4]);
    assert_eq!((s.dp_sum, s.dp_n, s.dp_sumsq), (310, 3, 49 + 90_000 + 9));
    let s = counts(&raw, Some(&[true, true, false, true]), 2).unwrap();
    assert_eq!((s.dp_sum, s.dp_n, s.samples), (10, 2, 3));
    // A negative depth and the end-of-vector value are no depth.
    let mut f = fields();
    f[1].data = ints(INT16, &[7, -4, missing(INT16) + 1, 3]);
    let s = counts(&record(4, &f), None, 2).unwrap();
    assert_eq!((s.dp_sum, s.dp_n, s.dp_sumsq), (10, 2, 58));
    // Every mask rule gives what the per-sample count gives, flags included.
    for mask in [
        MaskRules {
            min_dp: Some(5),
            ..NO_MASK
        },
        MaskRules {
            min_gq: Some(20),
            max_dp: Some(100),
            ..NO_MASK
        },
        MaskRules {
            mask_filtered: true,
            ..NO_MASK
        },
    ] {
        let ((fast, m1), (slow, m2)) = both(&raw, None, &mask, 2);
        assert_eq!(fast, slow, "{mask:?}");
        assert_eq!(m1, m2, "{mask:?}");
        assert!(fast.unwrap().filtered > 0, "{mask:?}");
    }
    // A depth that is not an integer is an error once a sample is read.
    let mut f = fields();
    f[1] = Field {
        key: DP,
        ty: FLOAT,
        n: 1,
        data: vec![0; 16],
    };
    let raw = record(4, &f);
    assert_eq!(
        counts(&raw, None, 2).unwrap_err(),
        "FORMAT/DP is not an integer field"
    );
    assert!(counts(&raw, Some(&[false; 4]), 2).is_ok());
}

/// One sample for the reference and its BCF values: the genotype, and how
/// the encoder writes what BCF can write in more than one way.
#[derive(Debug, Clone)]
struct Sample {
    gt: Option<NaiveGt>,
    /// The first allele's low bit, which nothing reads.
    first_bit: bool,
    /// Write an unphased `.` as the integer missing value.
    missing_value: bool,
}

/// A sample of a site with `n_alleles` alleles and room for `max_ploidy`.
/// `odd` is how often, in about 40, a genotype is not a full, called one:
/// 0 for a clean site.
fn sample(max_ploidy: usize, n_alleles: u32, odd: u32) -> BoxedStrategy<Sample> {
    // Mostly the first two alleles, as in real data.
    let called = prop_oneof![
        5 => (0..n_alleles.min(2)).prop_map(Some),
        1 => (0..n_alleles).prop_map(Some),
    ];
    let phase = prop::bool::weighted(0.9);
    let full = prop::collection::vec((called.clone(), phase), max_ploidy);
    let gt = if odd == 0 {
        full.prop_map(|g| Some(NaiveGt(g))).boxed()
    } else {
        // A missing allele, now and then an allele the site does not have,
        // fewer alleles than there is room for, or no genotype at all.
        let allele = prop_oneof![
            36 => called,
            odd => Just(None),
            1 => Just(Some(n_alleles)),
        ];
        let any = prop::collection::vec((allele, any::<bool>()), 1..=max_ploidy);
        prop_oneof![
            40 => full.prop_map(|g| Some(NaiveGt(g))),
            odd => any.prop_map(|g| Some(NaiveGt(g))),
            1 => Just(None),
        ]
        .boxed()
    };
    (gt, any::<bool>(), any::<bool>())
        .prop_map(|(gt, first_bit, missing_value)| Sample {
            gt,
            first_bit,
            missing_value,
        })
        .boxed()
}

#[derive(Debug, Clone)]
struct Site {
    ty: u8,
    n: usize,
    n_alleles: usize,
    samples: Vec<Sample>,
    keep: Option<Vec<bool>>,
}

fn site() -> impl Strategy<Value = Site> {
    let n_alleles = prop::sample::select(vec![1u32, 2, 2, 2, 2, 3, 4, 7, 62, 63, 200]);
    // Clean sites (every genotype plain) and sites with odd genotypes mixed in.
    let odd = prop::sample::select(vec![0u32, 0, 1, 4, 12]);
    let n = prop::sample::select(vec![1usize, 2, 2, 2, 2, 3, 4]);
    (n_alleles, odd, n, 0usize..70)
        .prop_flat_map(|(n_alleles, odd, n, n_samples)| {
            // The narrowest type that holds allele `n_alleles`, or a wider one.
            let fits: Vec<u8> = [(INT8, 62), (INT16, 16_382), (INT32, u32::MAX)]
                .iter()
                .filter(|(_, max)| n_alleles <= *max)
                .map(|(ty, _)| *ty)
                .collect();
            (
                prop::sample::select(fits),
                Just(n),
                Just(n_alleles as usize),
                prop::collection::vec(sample(n, n_alleles, odd), n_samples),
                prop::option::weighted(0.3, prop::collection::vec(any::<bool>(), n_samples)),
            )
        })
        .prop_map(|(ty, n, n_alleles, samples, keep)| Site {
            ty,
            n,
            n_alleles,
            samples,
            keep,
        })
}

/// The BCF values of `site`'s genotypes, each padded to `site.n`.
fn encode(site: &Site) -> Vec<i32> {
    let (m, e) = (missing(site.ty), missing(site.ty) + 1);
    let mut out = Vec::new();
    for s in &site.samples {
        let g = s.gt.as_ref().map_or(&[][..], |g| &g.0[..]);
        for (j, &(a, bar)) in g.iter().enumerate() {
            let bar = if j == 0 { s.first_bit } else { bar };
            out.push(if a.is_none() && !bar && s.missing_value {
                m
            } else {
                val(a, bar)
            });
        }
        out.extend(std::iter::repeat_n(e, site.n - g.len()));
    }
    out
}

/// One value of a GT vector, before its integer type is chosen.
#[derive(Debug, Clone, Copy)]
enum Cell {
    Value(i32),
    Missing,
    End,
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(2000))]

    /// Against the plain reference: every integer type, ploidy up to four
    /// with shorter and absent genotypes, missing alleles written both ways,
    /// phased and unphased, up to 200 alleles, with and without a selection.
    #[test]
    fn counts_match_the_plain_reference(site in site()) {
        let raw = gt_record(site.ty, site.n, &encode(&site));
        let gts: Vec<Option<NaiveGt>> = site.samples.iter().map(|s| s.gt.clone()).collect();
        let keep = site.keep.as_deref();
        let got = counts(&raw, keep, site.n_alleles);
        match naive_site_counts(&gts, keep, site.n_alleles) {
            Err(a) => prop_assert_eq!(
                got.unwrap_err(),
                format!("genotype allele {a} but the site lists {} alleles", site.n_alleles)
            ),
            Ok(want) => {
                let s = got.unwrap();
                prop_assert_eq!(s.called_alleles, want.allele_counts.iter().sum::<u64>());
                prop_assert_eq!(&s.allele_counts, &want.allele_counts);
                prop_assert_eq!(s.samples, want.samples);
                prop_assert_eq!(s.missing, want.missing);
                prop_assert_eq!(s.data_alleles, want.data_alleles);
                prop_assert_eq!(s.missing_alleles, want.missing_alleles);
                prop_assert_eq!(s.haploid, want.haploid);
                prop_assert_eq!(s.polyploid, want.polyploid);
                prop_assert_eq!(s.hwe_counts, want.diploid);
                prop_assert_eq!(s.all_phased, want.all_phased);
                prop_assert_eq!((s.filtered, s.dp_sum, s.dp_n, s.dp_sumsq), (0, 0, 0, 0));
            }
        }
    }

    /// Against the per-sample count on values no encoder writes: stray
    /// negative values, alleles out of range, values after the end of a
    /// vector, a GT of any type, other FORMAT fields of any type, mask rules.
    #[test]
    fn counts_equal_the_per_sample_count_on_any_values(
        ty in prop::sample::select(vec![INT8, INT8, INT8, INT16, INT16, INT32, INT32, FLOAT, CHAR, 0]),
        n in prop::sample::select(vec![0usize, 1, 2, 2, 2, 3]),
        n_alleles in 1usize..5,
        cells in (prop::sample::select(vec![0u32, 0, 1, 6]), 0usize..260).prop_flat_map(|(noise, len)| {
            // Mostly values a genotype has. A stray value ends the count with
            // an error, so `noise` says how many there are: none in a record
            // that is to be counted to its end.
            let odd = prop_oneof![
                6 => (0i32..12).prop_map(Cell::Value),
                3 => Just(Cell::End),
                2 => Just(Cell::Missing),
            ];
            let cell = if noise == 0 {
                prop_oneof![40 => (2i32..6).prop_map(Cell::Value), 11 => odd].boxed()
            } else {
                prop_oneof![
                    40 => (2i32..6).prop_map(Cell::Value),
                    11 => odd,
                    noise => (-200i32..300).prop_map(Cell::Value),
                ]
                .boxed()
            };
            prop::collection::vec(cell, len)
        }),
        keep in prop::option::weighted(0.3, prop::collection::vec(any::<bool>(), 0..66)),
        other in prop::option::weighted(0.3, (
            prop::sample::select(vec![DP, GQ, FT]),
            prop::sample::select(vec![INT8, INT16, INT32, FLOAT, CHAR, 0]),
            0usize..3,
            any::<i8>(),
        )),
        gt_first in any::<bool>(),
        mask in prop::option::weighted(0.2, (
            prop::option::of(0u32..40), prop::option::of(0u32..40), any::<bool>(),
        )),
    ) {
        // Whole samples only; a GT of no values, or of a type with no width, has no bytes.
        let per_sample = if ty == 0 { 0 } else { n };
        let n_samples = cells.len().checked_div(per_sample).unwrap_or(cells.len() % 7);
        let values: Vec<i32> = cells[..n_samples * per_sample]
            .iter()
            .map(|c| match c {
                Cell::Value(v) => *v,
                Cell::Missing => missing(ty),
                Cell::End => missing(ty) + 1,
            })
            .collect();
        let data = match ty {
            FLOAT => ints(INT32, &values),
            CHAR => ints(INT8, &values),
            _ => ints(ty, &values),
        };
        let mut fields = vec![Field { key: GT, ty, n, data }];
        if let Some((key, oty, on, fill)) = other {
            // Small values of either sign, the missing value and the end of
            // the vector (sums of squares of huge depths are not the subject).
            let values: Vec<i32> = (0..n_samples * on)
                .map(|i| match fill.wrapping_mul(i as i8 | 1) {
                    v @ (i8::MIN | -127) => missing(oty) + (i32::from(v) - i32::from(i8::MIN)),
                    v => i32::from(v),
                })
                .collect();
            let data = match oty {
                0 => Vec::new(),
                FLOAT => ints(INT32, &values),
                CHAR => ints(INT8, &values),
                _ => ints(oty, &values),
            };
            let f = Field { key, ty: oty, n: on, data };
            if gt_first { fields.push(f) } else { fields.insert(0, f) }
        }
        let raw = record(n_samples, &fields);
        let keep = keep.map(|mut k| {
            // Mostly the right length; a wrong one is an error both ways.
            if k.len() % 5 != 0 { k.resize(n_samples, true) }
            k
        });
        let mask = mask.map_or(NO_MASK, |(min_dp, min_gq, mask_filtered)| MaskRules {
            min_dp, min_gq, max_dp: None, mask_filtered,
        });
        let ((fast, m1), (slow, m2)) = both(&raw, keep.as_deref(), &mask, n_alleles);
        prop_assert_eq!(&fast, &slow);
        // After an error the flags are whatever the count had reached.
        if fast.is_ok() {
            prop_assert_eq!(m1, m2);
        }
    }
}
