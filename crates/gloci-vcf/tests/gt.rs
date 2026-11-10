mod common;
use common::*;
use gloci_vcf::{
    GenotypeSource, GtFields, SampleCache, SampleFields, bcf_genotypes, text_genotypes,
    text_genotypes_with,
};

type Row = (
    Option<(Vec<Option<u32>>, bool)>,
    Option<u32>,
    Option<u32>,
    bool,
);

fn all(src: &dyn GenotypeSource) -> Vec<Row> {
    let mut sf = SampleFields::default();
    (0..src.n_samples())
        .map(|i| {
            src.sample(i, &mut sf).unwrap();
            (
                sf.gt.as_ref().map(|g| (g.alleles.clone(), g.phased)),
                sf.dp,
                sf.gq,
                sf.ft_pass,
            )
        })
        .collect()
}

fn line(format: &str, samples: &[&str]) -> Vec<u8> {
    let mut l = format!("1\t10\t.\tA\tC,G\t50\tPASS\t.\t{format}");
    for s in samples {
        l.push('\t');
        l.push_str(s);
    }
    l.into_bytes()
}

#[test]
fn gt_text() {
    let l = line("GT", &["0/1", "1|0", "./.", "./1", "2", "0/1/2", "."]);
    let src = text_genotypes(&l).unwrap();
    assert_eq!(src.n_samples(), 7);
    let r = all(&src);
    let gt = |i: usize| r[i].0.clone().unwrap();
    assert_eq!(gt(0), (vec![Some(0), Some(1)], false));
    assert_eq!(gt(1), (vec![Some(1), Some(0)], true));
    assert_eq!(gt(2), (vec![None, None], false));
    assert_eq!(gt(3), (vec![None, Some(1)], false));
    assert_eq!(gt(4), (vec![Some(2)], false));
    assert_eq!(gt(5), (vec![Some(0), Some(1), Some(2)], false));
    assert_eq!(gt(6), (vec![None], false));
    assert!(r.iter().all(|x| x.3 && x.1.is_none() && x.2.is_none()));
}

#[test]
fn gt_text_format_order_and_short_sample() {
    let l = line(
        "DP:GT:GQ:FT",
        &["7:0/1:30:PASS", "5:1|1:.:LowQual", "9:0/0", "3"],
    );
    let r = all(&text_genotypes(&l).unwrap());
    assert_eq!(
        r[0],
        (
            Some((vec![Some(0), Some(1)], false)),
            Some(7),
            Some(30),
            true
        )
    );
    assert_eq!(
        r[1],
        (Some((vec![Some(1), Some(1)], true)), Some(5), None, false)
    );
    // Trailing fields missing: GQ absent, FT absent (passes).
    assert_eq!(
        r[2],
        (Some((vec![Some(0), Some(0)], false)), Some(9), None, true)
    );
    // Only DP present: GT missing entirely.
    assert_eq!(r[3], (None, Some(3), None, true));
}

#[test]
fn gt_text_without_gt_key_and_errors() {
    let l = line("DP", &["4"]);
    let r = all(&text_genotypes(&l).unwrap());
    assert_eq!(r[0], (None, Some(4), None, true));
    let l = line("GT", &["0/x"]);
    let src = text_genotypes(&l).unwrap();
    assert!(src.sample(0, &mut SampleFields::default()).is_err());
    assert!(text_genotypes(b"1\t10\t.\tA\tC").is_err());
}

#[test]
fn gt_text_no_samples() {
    let l = b"1\t10\t.\tA\tC\t50\tPASS\t.".to_vec();
    assert_eq!(text_genotypes(&l).unwrap().n_samples(), 0);
}

#[test]
fn gt_bcf_equals_text() {
    let recs = [
        "1\t10\t.\tA\tC,G\t50\tPASS\t.\tGT:DP:GQ:FT\t0/1:7:30:PASS\t1|0:5:.:LowQual\t./.:.:.:.\t./1:12:99:PASS",
        "1\t20\t.\tA\tC\t50\tPASS\t.\tDP:GT\t9:0/1\t3:2/2\t4:./.\t1:0|1",
        "1\t30\t.\tA\tC\t50\tPASS\t.\tGT\t0/1\t1\t0/0\t1/1",
    ];
    let text = vcf(&["s1", "s2", "s3", "s4"], &recs);
    let bcf = gloci_testkit::fixtures::vcf_to_bcf(text.as_bytes());
    let (dicts, raws) = bcf_records(&bcf);
    assert_eq!(raws.len(), 3);
    for (r, raw) in recs.iter().zip(&raws) {
        let t = text_genotypes(r.as_bytes()).unwrap();
        let b = bcf_genotypes(raw, &dicts).unwrap();
        assert_eq!(b.n_samples(), 4);
        assert_eq!(all(&t), all(&b), "record {r}");
    }
}

/// Mixed ploidy and a missing-GT sample, as hand-encoded BCF (noodles pads
/// mixed-ploidy records inconsistently, so the fixture encoder can't be used).
#[test]
fn gt_bcf_mixed_ploidy_handbuilt() {
    let dicts = gloci_io::BcfDicts {
        contigs: vec![b"1".to_vec()],
        strings: vec![b"PASS".to_vec(), b"GT".to_vec(), b"DP".to_vec()],
    };
    let mut shared = vec![0u8; 24];
    shared[20..24].copy_from_slice(&((2u32 << 24) | 4).to_le_bytes()); // n_fmt=2, n_sample=4
    shared.extend_from_slice(&[0x07, 0x07]); // typed ID, REF placeholders (unused by the decoder)
    let mut indiv = vec![0x11, 1, 0x31]; // key GT, int8 x3 per sample
    //            0/1/2           2 (haploid)         ./.  (as 0,0,End)     1|0|1
    indiv.extend_from_slice(&[2, 4, 6, 6, 0x81, 0x81, 0, 0, 0x81, 4 | 1, 2 | 1, 4 | 1]);
    indiv.extend_from_slice(&[0x11, 2, 0x11, 7, 0x80, 3, 4]); // key DP, int8 x1: 7, missing, 3, 4
    let mut raw = Vec::new();
    raw.extend_from_slice(&(shared.len() as u32).to_le_bytes());
    raw.extend_from_slice(&(indiv.len() as u32).to_le_bytes());
    raw.extend_from_slice(&shared);
    raw.extend_from_slice(&indiv);
    let line = b"1\t10\t.\tA\tC,G\t50\tPASS\t.\tGT:DP\t0/1/2:7\t2:.\t./.:3\t1|0|1:4";
    let b = bcf_genotypes(&raw, &dicts).unwrap();
    assert_eq!(b.n_samples(), 4);
    assert_eq!(all(&text_genotypes(line).unwrap()), all(&b));
}

#[test]
fn gt_text_high_allele_and_mixed_phase() {
    let l = line("GT", &["0/12", "10|11", "0/1|2", "0|1|2"]);
    let r = all(&text_genotypes(&l).unwrap());
    assert_eq!(r[0].0, Some((vec![Some(0), Some(12)], false)));
    assert_eq!(r[1].0, Some((vec![Some(10), Some(11)], true)));
    assert_eq!(r[2].0, Some((vec![Some(0), Some(1), Some(2)], false)));
    assert_eq!(r[3].0, Some((vec![Some(0), Some(1), Some(2)], true)));
}

fn frame(shared: &[u8], indiv: &[u8]) -> Vec<u8> {
    let mut raw = Vec::new();
    raw.extend_from_slice(&(shared.len() as u32).to_le_bytes());
    raw.extend_from_slice(&(indiv.len() as u32).to_le_bytes());
    raw.extend_from_slice(shared);
    raw.extend_from_slice(indiv);
    raw
}

/// int16 GT/DP, the `Int::Missing` GT sentinel, NUL-padded Char FT, and a
/// sample with no value for a FORMAT field.
#[test]
fn gt_bcf_int16_missing_sentinel_ft_and_absent_values() {
    let dicts = gloci_io::BcfDicts {
        contigs: vec![b"1".to_vec()],
        strings: ["PASS", "GT", "DP", "FT", "GQ"]
            .iter()
            .map(|s| s.as_bytes().to_vec())
            .collect(),
    };
    let mut shared = vec![0u8; 24];
    shared[20..24].copy_from_slice(&((4u32 << 24) | 3).to_le_bytes()); // n_fmt=4, n_sample=3
    shared.extend_from_slice(&[0x07, 0x07]);
    let i16s = |v: &[i16]| v.iter().flat_map(|x| x.to_le_bytes()).collect::<Vec<u8>>();
    let mut indiv = vec![0x11, 1, 0x22]; // GT: int16 x2
    // 0/70 (allele 70 -> (71<<1)=142, beyond int8), 70|0 phased, missing sentinel then end
    indiv.extend(i16s(&[2, 142, 142, 3, i16::MIN, -32767]));
    indiv.extend_from_slice(&[0x11, 2, 0x12]); // DP: int16 x1
    indiv.extend(i16s(&[300, i16::MIN, 4000]));
    indiv.extend_from_slice(&[0x11, 3, 0x87]); // FT: char x8
    indiv.extend_from_slice(b"PASS\0\0\0\0LowQual\0.\0\0\0\0\0\0\0");
    indiv.extend_from_slice(&[0x11, 4, 0x11]); // GQ: int8 x1, sample 2 has the missing value
    indiv.extend_from_slice(&[40, 0x80, 41]);
    let raw = frame(&shared, &indiv);
    let b = bcf_genotypes(&raw, &dicts).unwrap();
    let r = all(&b);
    assert_eq!(
        r[0],
        (
            Some((vec![Some(0), Some(70)], false)),
            Some(300),
            Some(40),
            true
        )
    );
    assert_eq!(
        r[1],
        (Some((vec![Some(70), Some(0)], true)), None, None, false)
    );
    assert_eq!(
        r[2],
        (Some((vec![None], false)), Some(4000), Some(41), true)
    );
}

#[test]
fn gt_bcf_field_with_no_value_for_any_sample() {
    // GT present as type Missing (count 0): every sample has no GT.
    let dicts = gloci_io::BcfDicts {
        contigs: vec![b"1".to_vec()],
        strings: ["PASS", "GT"]
            .iter()
            .map(|s| s.as_bytes().to_vec())
            .collect(),
    };
    let mut shared = vec![0u8; 24];
    shared[20..24].copy_from_slice(&((1u32 << 24) | 2).to_le_bytes());
    shared.extend_from_slice(&[0x07, 0x07]);
    let raw = frame(&shared, &[0x11, 1, 0x00]);
    let r = all(&bcf_genotypes(&raw, &dicts).unwrap());
    assert_eq!(r.len(), 2);
    assert!(r.iter().all(|x| x.0.is_none() && x.3));
}

#[test]
fn allele_buffer_survives_across_samples() {
    // Text: sample 1 has no GT, so the buffer is parked and reused.
    let l = line("DP:GT", &["5:0/1", "3", "5:1/1"]);
    let t = text_genotypes(&l).unwrap();
    let mut sf = SampleFields::default();
    t.sample(0, &mut sf).unwrap();
    let (ptr, cap) = {
        let a = &sf.gt.as_ref().unwrap().alleles;
        (a.as_ptr(), a.capacity())
    };
    t.sample(1, &mut sf).unwrap();
    assert!(sf.gt.is_none());
    t.sample(2, &mut sf).unwrap();
    let a = &sf.gt.as_ref().unwrap().alleles;
    assert_eq!((a.as_ptr(), a.capacity()), (ptr, cap));

    // BCF: same, with a GT-less middle sample (typed GT of all-End).
    let dicts = gloci_io::BcfDicts {
        contigs: vec![b"1".to_vec()],
        strings: ["PASS", "GT"]
            .iter()
            .map(|s| s.as_bytes().to_vec())
            .collect(),
    };
    let mut shared = vec![0u8; 24];
    shared[20..24].copy_from_slice(&((1u32 << 24) | 3).to_le_bytes());
    shared.extend_from_slice(&[0x07, 0x07]);
    let raw = frame(&shared, &[0x11, 1, 0x21, 2, 4, 0x81, 0x81, 4, 4]);
    let b = bcf_genotypes(&raw, &dicts).unwrap();
    let mut sf = SampleFields::default();
    b.sample(0, &mut sf).unwrap();
    let (ptr, cap) = {
        let a = &sf.gt.as_ref().unwrap().alleles;
        (a.as_ptr(), a.capacity())
    };
    b.sample(1, &mut sf).unwrap();
    assert!(sf.gt.is_none());
    b.sample(2, &mut sf).unwrap();
    let a = &sf.gt.as_ref().unwrap().alleles;
    assert_eq!((a.as_ptr(), a.capacity()), (ptr, cap));
}

/// Encodes `l` through the native GT path and decodes it again.
fn via_bcf(l: &[u8]) -> Vec<Row> {
    let header =
        gloci_io::VcfHeader::parse(vcf(&["s1", "s2", "s3", "s4", "s5"], &[]).into_bytes()).unwrap();
    let dicts = gloci_io::BcfDicts::from_header(&header.raw).unwrap();
    let mut enc = gloci_io::TextToBcf::new("t.vcf", &header).unwrap();
    let mut rec = Vec::new();
    enc.encode(l, 1, &mut rec).unwrap();
    all(&bcf_genotypes(&rec, &dicts).unwrap())
}

#[test]
fn encode_gt_roundtrip() {
    let alts = (1..=70)
        .map(|i| format!("T{i}"))
        .collect::<Vec<_>>()
        .join(",");
    let cases: Vec<(&str, Vec<&str>)> = vec![
        // A missing haploid next to diploid calls, then a triploid.
        ("GT", vec!["0/1", ".", "1|0", "./.", "0/1/2"]),
        ("GT", vec![".", ".", "1", "0|1", "."]),
        // Alleles past int8 (63) and int16 widths, phased and unphased.
        ("GT", vec!["63/64", "70|1", ".", "0", "65|65"]),
        // GT after another field, and with other fields following.
        (
            "DP:GT:GQ",
            vec!["5:0/1:9", "5:.:9", "5:1|1:9", "5:./1:9", "5:0:9"],
        ),
        // A whole-sample `.` and dropped trailing subfields, GT not first.
        ("DP:GT", vec!["5:0/1", ".", "7:1|0", "8:.", "9:0/1|2"]),
        ("GT:DP:GQ", vec!["0/1:1:2", "1", "0|1:3", ".:4:5", "./."]),
        ("DP:GT:GQ", vec!["5:0/1:9", "6", "7:1", "8:.", "."]),
    ];
    for (fmt, gts) in cases {
        let mut l = format!("1\t10\t.\tA\t{alts}\t50\tPASS\t.\t{fmt}").into_bytes();
        for g in &gts {
            l.push(b'\t');
            l.extend_from_slice(g.as_bytes());
        }
        let mut want = all(&text_genotypes(&l).unwrap());
        // BCF writes a sample with no GT as one missing allele (as htslib
        // does), so it reads back as `.` where the text reader says absent.
        for r in &mut want {
            r.0.get_or_insert((vec![None], false));
        }
        assert_eq!(via_bcf(&l), want, "{fmt} {gts:?}");
    }
}

#[test]
fn encode_gt_roundtrip_without_gt_in_format() {
    let l = b"1\t10\t.\tA\tC\t50\tPASS\t.\tDP\t5\t6\t7\t8\t9".to_vec();
    let want = all(&text_genotypes(&l).unwrap());
    assert!(want.iter().all(|r| r.0.is_none()));
    assert_eq!(via_bcf(&l), want);
}

#[test]
fn bcf_output_refuses_a_leading_phase_separator() {
    let header = gloci_io::VcfHeader::parse(vcf(&["s1"], &[]).into_bytes()).unwrap();
    let mut enc = gloci_io::TextToBcf::new("t.vcf", &header).unwrap();
    for gt in ["|0|1", "/0|1", "|1"] {
        let l = format!("1\t10\t.\tA\tC\t50\tPASS\t.\tGT\t{gt}");
        let e = enc
            .encode(l.as_bytes(), 7, &mut Vec::new())
            .unwrap_err()
            .to_string();
        assert!(e.contains("t.vcf") && e.contains("record 7"), "{e}");
        assert!(
            e.contains("phase separator") && e.contains("remove the leading"),
            "{e}"
        );
    }
}

/// Counts how often each sample is decoded.
struct Counting<'a> {
    inner: &'a dyn GenotypeSource,
    calls: std::cell::RefCell<Vec<u32>>,
}

impl GenotypeSource for Counting<'_> {
    fn n_samples(&self) -> usize {
        self.inner.n_samples()
    }
    fn sample(&self, i: usize, out: &mut SampleFields) -> Result<(), String> {
        self.calls.borrow_mut()[i] += 1;
        self.inner.sample(i, out)
    }
}

#[test]
fn sample_cache_decodes_once_and_reads_alike() {
    let l1 = line(
        "GT:DP:GQ:FT",
        &["0|1:5:20:PASS", "./.:.:.:LowQ", ".", "1/2:9:30:."],
    );
    let l2 = line("GT:DP", &["1/1:7", "0/0:3", "0/1:4", "1|1:8"]);
    let cache = SampleCache::default();
    for l in [&l1, &l2, &l1] {
        let text = text_genotypes(l).unwrap();
        let counting = Counting {
            inner: &text,
            calls: std::cell::RefCell::new(vec![0; text.n_samples()]),
        };
        cache.start_record(text.n_samples());
        let cached = cache.wrap(&counting);
        // Two consumers read every sample, in different orders.
        let first = all(&cached);
        let mut sf = SampleFields::default();
        let mut second = Vec::new();
        for i in (0..cached.n_samples()).rev() {
            cached.sample(i, &mut sf).unwrap();
            second.push(i);
        }
        assert_eq!(first, all(&text));
        assert_eq!(*counting.calls.borrow(), vec![1; text.n_samples()]);
        assert_eq!(second.len(), text.n_samples());
        // A fresh view over the same cache (not started again) reads the same.
        let again = cache.wrap(&text);
        assert_eq!(all(&again), first);
    }
}

#[test]
fn sample_cache_without_start_record_decodes_directly() {
    let l = line("GT:DP", &["0/1:5", "1/1:6"]);
    let text = text_genotypes(&l).unwrap();
    let cache = SampleCache::default();
    assert_eq!(all(&cache.wrap(&text)), all(&text));
}

#[test]
fn gt_only_decoding_skips_optional_fields_and_keeps_gt() {
    let l = line("GT:DP:GQ:FT", &["0|1:5:20:LowQ", "1/1:9:30:PASS"]);
    let full = text_genotypes(&l).unwrap();
    let only = text_genotypes_with(&l, GtFields::GT_ONLY).unwrap();
    let dp_only = text_genotypes_with(
        &l,
        GtFields {
            dp: true,
            ..GtFields::GT_ONLY
        },
    )
    .unwrap();
    for (a, b) in all(&full).into_iter().zip(all(&only)) {
        assert_eq!(a.0, b.0);
        assert_eq!((b.1, b.2, b.3), (None, None, true));
    }
    let rows = all(&dp_only);
    assert_eq!((rows[0].1, rows[0].2, rows[0].3), (Some(5), None, true));
}
