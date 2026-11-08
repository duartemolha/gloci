mod common;
use common::*;
use gloci_vcf::{GenotypeSource, SampleFields, bcf_genotypes, text_genotypes};

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
