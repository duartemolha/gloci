use xtask::synth::{BedOpts, Profile, VcfOpts};

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
}
