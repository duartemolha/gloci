use xtask::synth::{BedOpts, PairOpts, Profile, VcfOpts};

fn bed_bytes(o: &BedOpts) -> Vec<u8> {
    let mut v = Vec::new();
    xtask::synth::bed(o, &mut v).unwrap();
    v
}

fn opts(profile: Profile) -> BedOpts {
    BedOpts {
        records: 1000,
        chroms: 5,
        seed: 7,
        profile,
        ..BedOpts::default()
    }
}

const PROFILES: [Profile; 4] = [
    Profile::Uniform,
    Profile::DeepNesting,
    Profile::WholeChrom,
    Profile::Points,
];

type Rec = (String, u64, u64);

fn parse(bytes: &[u8]) -> Vec<Rec> {
    std::str::from_utf8(bytes)
        .unwrap()
        .lines()
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            assert!(f.len() >= 3, "too few columns: {l}");
            (
                f[0].to_string(),
                f[1].parse().unwrap(),
                f[2].parse().unwrap(),
            )
        })
        .collect()
}

#[test]
fn same_seed_same_bytes() {
    for p in PROFILES {
        assert_eq!(bed_bytes(&opts(p)), bed_bytes(&opts(p)), "{p:?}");
    }
    let other = BedOpts {
        seed: 8,
        ..opts(Profile::Uniform)
    };
    assert_ne!(bed_bytes(&opts(Profile::Uniform)), bed_bytes(&other));
}

#[test]
fn profiles_are_valid_bed() {
    for p in PROFILES {
        let recs = parse(&bed_bytes(&opts(p)));
        assert_eq!(recs.len(), 1000, "{p:?}");
        let mut chroms = std::collections::BTreeSet::new();
        for (c, s, e) in &recs {
            assert!(s <= e, "{p:?}: {c} {s} {e}");
            chroms.insert(c.clone());
        }
        assert!(chroms.len() <= 5 && chroms.len() > 1, "{p:?}");
    }
}

#[test]
fn profile_shapes() {
    let pts = parse(&bed_bytes(&opts(Profile::Points)));
    assert!(pts.iter().all(|(_, s, e)| s == e));
    let whole = parse(&bed_bytes(&opts(Profile::WholeChrom)));
    assert!(whole.iter().all(|(_, s, e)| e - s >= 50_000_000));
    // deep nesting: many records share a centre, so intervals contain one another
    let deep = parse(&bed_bytes(&opts(Profile::DeepNesting)));
    let mids: std::collections::BTreeSet<(String, u64)> = deep
        .iter()
        .map(|(c, s, e)| (c.clone(), (s + e) / 2))
        .collect();
    assert!(mids.len() <= 5 * 100 * 2, "{}", mids.len());
    let sorted = parse(&bed_bytes(&BedOpts {
        sorted: true,
        ..opts(Profile::Uniform)
    }));
    assert!(
        sorted
            .windows(2)
            .all(|w| (&w[0].0, w[0].1) <= (&w[1].0, w[1].1))
    );
}

#[test]
fn length_bounds_respected() {
    let o = BedOpts {
        len_min: 100,
        len_max: 200,
        ..opts(Profile::Uniform)
    };
    for (_, s, e) in parse(&bed_bytes(&o)) {
        assert!((100..=200).contains(&(e - s)));
    }
}

#[test]
fn vcf_shape_and_determinism() {
    let o = VcfOpts {
        records: 500,
        samples: 3,
        seed: 3,
        chroms: 2,
    };
    let (mut a, mut b) = (Vec::new(), Vec::new());
    xtask::synth::vcf(&o, &mut a).unwrap();
    xtask::synth::vcf(&o, &mut b).unwrap();
    assert_eq!(a, b);
    let text = String::from_utf8(a).unwrap();
    let body: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(body.len(), 500);
    for l in &body {
        assert_eq!(l.split('\t').count(), 9 + 3, "{l}");
    }
    assert!(
        body.iter()
            .any(|l| l.contains("<DEL>") && l.contains("END="))
    );
    assert!(
        body.iter()
            .any(|l| l.split('\t').nth(4).unwrap().contains(','))
    );
    assert!(text.lines().any(|l| l.starts_with("#CHROM")));

    // Check VCF header lines for DP, GQ, FT and FILTERs
    assert!(text.contains("##FORMAT=<ID=DP,"));
    assert!(text.contains("##FORMAT=<ID=GQ,"));
    assert!(text.contains("##FORMAT=<ID=FT,"));
    assert!(text.contains("##FILTER=<ID=LowQual,"));

    // Check sorted POS per chromosome
    let mut last_chrom = "";
    let mut last_pos = 0u64;
    for l in &body {
        let cols: Vec<&str> = l.split('\t').collect();
        let chrom = cols[0];
        let pos: u64 = cols[1].parse().unwrap();
        if chrom != last_chrom {
            last_chrom = chrom;
            last_pos = pos;
        } else {
            assert!(pos > last_pos, "POS not sorted: {pos} <= {last_pos}");
            last_pos = pos;
        }
    }

    // Check FORMAT is GT:DP:GQ:FT
    assert!(
        body.iter()
            .all(|l| l.split('\t').nth(8) == Some("GT:DP:GQ:FT"))
    );

    // Check presence of DP/GQ/FT in sample values and mix of QUAL/FILTER
    let has_missing_gt = body
        .iter()
        .any(|l| l.split('\t').skip(9).any(|s| s.starts_with("./.")));
    let has_called_gt = body
        .iter()
        .any(|l| l.split('\t').skip(9).any(|s| !s.starts_with("./.")));
    assert!(has_missing_gt, "should have missing genotypes");
    assert!(has_called_gt, "should have called genotypes");

    let has_pass = body.iter().any(|l| l.split('\t').nth(6) == Some("PASS"));
    let has_non_pass = body.iter().any(|l| {
        let f = l.split('\t').nth(6).unwrap();
        f != "PASS"
    });
    assert!(has_pass, "should have PASS filters");
    assert!(
        has_non_pass,
        "should have non-PASS filters (LowQual, q10, etc.)"
    );

    let has_missing_qual = body.iter().any(|l| l.split('\t').nth(5) == Some("."));
    let has_num_qual = body.iter().any(|l| {
        let q = l.split('\t').nth(5).unwrap();
        q != "." && q.parse::<f64>().is_ok()
    });
    assert!(has_missing_qual, "should have some missing QUAL");
    assert!(has_num_qual, "should have numeric QUAL");
}

#[test]
fn pair_generator_depth_and_determinism() {
    let o = PairOpts {
        records: 500,
        depth: 2.0,
        chroms: 2,
        seed: 42,
        len_min: 50,
        len_max: 500,
        sorted: true,
    };
    let (mut a1, mut b1) = (Vec::new(), Vec::new());
    let (mut a2, mut b2) = (Vec::new(), Vec::new());
    xtask::synth::pair(&o, &mut a1, &mut b1).unwrap();
    xtask::synth::pair(&o, &mut a2, &mut b2).unwrap();
    assert_eq!(a1, a2);
    assert_eq!(b1, b2);

    let recs_a = parse(&a1);
    let recs_b = parse(&b1);
    assert_eq!(recs_a.len(), 500);

    // Verify sortedness
    assert!(
        recs_a
            .windows(2)
            .all(|w| (&w[0].0, w[0].1) <= (&w[1].0, w[1].1))
    );
    assert!(
        recs_b
            .windows(2)
            .all(|w| (&w[0].0, w[0].1) <= (&w[1].0, w[1].1))
    );

    // Measure overlap depth: average number of B intervals overlapping an A interval
    let mut total_overlaps = 0;
    for (ca, sa, ea) in &recs_a {
        for (cb, sb, eb) in &recs_b {
            if ca == cb && *sa < *eb && *sb < *ea {
                total_overlaps += 1;
            }
        }
    }
    let avg_depth = total_overlaps as f64 / recs_a.len() as f64;
    // Expected around 2.0
    assert!(
        (1.2..=2.8).contains(&avg_depth),
        "expected avg depth ~2.0, got {avg_depth}"
    );

    // Test depth = 0
    let o_zero = PairOpts {
        records: 100,
        depth: 0.0,
        chroms: 2,
        seed: 42,
        len_min: 50,
        len_max: 500,
        sorted: true,
    };
    let (mut a_z, mut b_z) = (Vec::new(), Vec::new());
    xtask::synth::pair(&o_zero, &mut a_z, &mut b_z).unwrap();
    let recs_az = parse(&a_z);
    let recs_bz = parse(&b_z);
    let mut overlaps_zero = 0;
    for (ca, sa, ea) in &recs_az {
        for (cb, sb, eb) in &recs_bz {
            if ca == cb && *sa < *eb && *sb < *ea {
                overlaps_zero += 1;
            }
        }
    }
    assert_eq!(overlaps_zero, 0, "depth 0 must yield zero overlaps");
}
