use xtask::synth::{BedOpts, BedpeOpts, PairOpts, Profile, VcfOpts};

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
    assert_eq!(recs_b.len(), 500);

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
    assert_eq!(recs_az.len(), 100);
    assert_eq!(recs_bz.len(), 100);
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

#[test]
fn high_density_pair_generator_distribution_and_depth() {
    let chrom_len = 250_000_000u64;
    // 50,000 records on 1 chromosome with len_max 10,000 => 50,000 * 10,000 = 500,000,000 > CHROM_LEN
    let o = PairOpts {
        records: 50_000,
        depth: 2.0,
        chroms: 1,
        seed: 99,
        len_min: 50,
        len_max: 10_000,
        sorted: true,
    };
    let (mut a, mut b) = (Vec::new(), Vec::new());
    xtask::synth::pair(&o, &mut a, &mut b).unwrap();
    let recs_a = parse(&a);
    let recs_b = parse(&b);

    assert_eq!(recs_a.len(), 50_000);
    assert_eq!(recs_b.len(), 50_000);

    // Verify bounds strictly within 0..CHROM_LEN without piling up at the end
    let min_sa = recs_a.iter().map(|r| r.1).min().unwrap();
    let max_ea = recs_a.iter().map(|r| r.2).max().unwrap();
    let min_sb = recs_b.iter().map(|r| r.1).min().unwrap();
    let max_eb = recs_b.iter().map(|r| r.2).max().unwrap();

    assert!(min_sa > 0, "min_sa must be > 0, got {min_sa}");
    assert!(
        max_ea < chrom_len,
        "max_ea must be < chrom_len, got {max_ea}"
    );
    assert!(min_sb > 0, "min_sb must be > 0, got {min_sb}");
    assert!(
        max_eb < chrom_len,
        "max_eb must be < chrom_len, got {max_eb}"
    );

    // Verify intervals span across the whole chromosome
    assert!(
        min_sa < 100_000,
        "first intervals should start near start of chromosome, got {min_sa}"
    );
    assert!(
        max_ea > 249_000_000,
        "last intervals should reach end of chromosome, got {max_ea}"
    );

    // Verify no clamping pile-up: median distance between consecutive starts should be > 0
    let last_100_starts: Vec<u64> = recs_a.iter().rev().take(100).map(|r| r.1).collect();
    let distinct_last_100: std::collections::HashSet<u64> =
        last_100_starts.iter().copied().collect();
    assert!(
        distinct_last_100.len() >= 90,
        "intervals must not pile up at the end: {} distinct out of 100",
        distinct_last_100.len()
    );

    // Measure overlap depth using two-pointer scan over sorted intervals
    let mut total_overlaps = 0usize;
    let mut b_start_idx = 0;
    for (_ca, sa, ea) in &recs_a {
        while b_start_idx < recs_b.len() && recs_b[b_start_idx].2 <= *sa {
            b_start_idx += 1;
        }
        let mut j = b_start_idx;
        while j < recs_b.len() && recs_b[j].1 < *ea {
            total_overlaps += 1;
            j += 1;
        }
    }
    let avg_depth = total_overlaps as f64 / recs_a.len() as f64;
    assert!(
        (1.5..=2.5).contains(&avg_depth),
        "expected avg depth ~2.0, got {avg_depth}"
    );

    // Test high density with depth = 0
    let o_zero = PairOpts {
        records: 10_000,
        depth: 0.0,
        chroms: 1,
        seed: 99,
        len_min: 50,
        len_max: 10_000,
        sorted: true,
    };
    let (mut a_z, mut b_z) = (Vec::new(), Vec::new());
    xtask::synth::pair(&o_zero, &mut a_z, &mut b_z).unwrap();
    let recs_az = parse(&a_z);
    let recs_bz = parse(&b_z);
    assert_eq!(recs_az.len(), 10_000);
    assert_eq!(recs_bz.len(), 10_000);

    let mut overlaps_zero = 0usize;
    let mut bz_idx = 0;
    for (_ca, sa, ea) in &recs_az {
        while bz_idx < recs_bz.len() && recs_bz[bz_idx].2 <= *sa {
            bz_idx += 1;
        }
        let mut j = bz_idx;
        while j < recs_bz.len() && recs_bz[j].1 < *ea {
            overlaps_zero += 1;
            j += 1;
        }
    }
    assert_eq!(
        overlaps_zero, 0,
        "high-density depth 0 must have zero overlaps"
    );
}

/// The synthetic multi-sample VCF must exercise every `vcf stats` statistic
/// of the benchmark: missing genotypes per site and per sample, hets, both
/// transitions and transversions, and a depth spread around the usual
/// `--min-dp` thresholds.
#[test]
fn vcf_has_the_mix_the_five_statistics_need() {
    let o = VcfOpts {
        records: 400,
        samples: 200,
        seed: 1,
        chroms: 1,
    };
    let mut v = Vec::new();
    xtask::synth::vcf(&o, &mut v).unwrap();
    let text = String::from_utf8(v).unwrap();
    let (mut missing, mut hets, mut homs, mut ts, mut tv) = (0, 0, 0, 0, 0);
    let (mut low_dp, mut high_dp) = (0, 0);
    let mut per_sample_missing = vec![0u32; 200];
    for l in text.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = l.split('\t').collect();
        let (r, a) = (f[3], f[4]);
        if r.len() == 1 && a.len() == 1 {
            match (r, a) {
                ("A", "G") | ("G", "A") | ("C", "T") | ("T", "C") => ts += 1,
                _ => tv += 1,
            }
        }
        for (i, g) in f[9..].iter().enumerate() {
            let p: Vec<&str> = g.split(':').collect();
            match p[0] {
                "./." => {
                    missing += 1;
                    per_sample_missing[i] += 1;
                }
                gt => {
                    let (x, y) = gt.split_once('/').unwrap();
                    if x == y {
                        homs += 1;
                    } else {
                        hets += 1;
                    }
                    let dp: u32 = p[1].parse().unwrap();
                    low_dp += u32::from(dp < 10);
                    high_dp += u32::from(dp >= 50);
                }
            }
        }
    }
    let total = 400 * 200;
    assert!(
        missing > total / 40 && missing < total / 10,
        "missing {missing}"
    );
    assert!(
        hets > total / 10 && homs > total / 10,
        "hets {hets} homs {homs}"
    );
    assert!(ts > 20 && tv > 20, "ts {ts} tv {tv}");
    assert!(low_dp > 0 && high_dp > 0);
    assert!(
        per_sample_missing.iter().all(|&m| m > 0),
        "a sample with no missing genotype"
    );
}

// ---- SAM and FASTA ----

use std::collections::HashMap;
use xtask::synth::{FastaOpts, SamOpts};

fn sam_opts() -> SamOpts {
    SamOpts {
        records: 12000,
        read_len: 100,
        depth: 5.0,
        chroms: 3,
        spliced_frac: 0.2,
        seed: 7,
        ..SamOpts::default()
    }
}

fn sam_bytes(o: &SamOpts) -> (Vec<u8>, Vec<u8>) {
    let (mut s, mut r) = (Vec::new(), Vec::new());
    xtask::synth::sam(o, &mut s).unwrap();
    xtask::synth::sam_reference(o, &mut r).unwrap();
    (s, r)
}

fn read_fasta(bytes: &[u8]) -> HashMap<String, Vec<u8>> {
    let mut m = HashMap::new();
    let mut name = String::new();
    for l in std::str::from_utf8(bytes).unwrap().lines() {
        if let Some(n) = l.strip_prefix('>') {
            name = n.to_string();
            m.insert(name.clone(), Vec::new());
        } else {
            m.get_mut(&name).unwrap().extend_from_slice(l.as_bytes());
        }
    }
    m
}

/// The CIGAR as (length, op) pairs.
fn cigar(c: &str) -> Vec<(usize, char)> {
    let mut out = Vec::new();
    let mut n = 0usize;
    for ch in c.chars() {
        if let Some(d) = ch.to_digit(10) {
            n = n * 10 + d as usize;
        } else {
            out.push((n, ch));
            n = 0;
        }
    }
    out
}

#[test]
fn gen_sam_is_valid() {
    let o = sam_opts();
    let (sam, reference) = sam_bytes(&o);
    let refs = read_fasta(&reference);
    let text = std::str::from_utf8(&sam).unwrap();
    assert!(text.starts_with("@HD\tVN:1.6\tSO:coordinate\n"));
    let len = xtask::synth::sam_chrom_len(&o);
    for c in 1..=3 {
        assert!(text.contains(&format!("@SQ\tSN:chr{c}\tLN:{len}\n")));
        assert_eq!(refs[&format!("chr{c}")].len() as u64, len);
    }
    let (mut n, mut spliced, mut matched) = (0u64, 0u64, 0u64);
    let mut last = (0usize, 0u64, 0u32);
    for l in text.lines().filter(|l| !l.starts_with('@')) {
        let f: Vec<&str> = l.split('\t').collect();
        assert_eq!(f.len(), 11, "{l}");
        let chrom: usize = f[2].strip_prefix("chr").unwrap().parse().unwrap();
        let pos: u64 = f[3].parse().unwrap();
        // Coordinate order: by header chromosome, then position.
        // (samtools order: at one position, a forward read before a reverse one)
        let flag: u32 = f[1].parse().unwrap();
        let key = (chrom, pos, flag & 16);
        assert!(key >= last, "unsorted at {l}");
        last = key;
        // The read is the reference at its blocks (the gaps skipped).
        let rf = &refs[f[2]];
        let mut at = pos as usize - 1;
        let mut want = Vec::new();
        for (k, op) in cigar(f[5]) {
            match op {
                'M' => {
                    want.extend_from_slice(&rf[at..at + k]);
                    matched += k as u64;
                    at += k;
                }
                'N' => {
                    spliced += 1;
                    at += k;
                }
                'D' => at += k,
                other => panic!("unexpected CIGAR op {other}"),
            }
        }
        assert_eq!(
            f[9].as_bytes(),
            &want[..],
            "read {l} differs from the reference"
        );
        assert_eq!(f[9].len(), 100);
        assert_eq!(f[10].len(), 100);
        n += 1;
    }
    assert_eq!(n, 12000);
    let frac = spliced as f64 / n as f64;
    assert!((frac - 0.2).abs() < 0.03, "spliced fraction {frac}");
    let depth = matched as f64 / (3.0 * len as f64);
    assert!((depth - 5.0).abs() < 0.5, "depth {depth}");
}

#[test]
fn gen_sam_shuffled_is_unsorted_and_deterministic() {
    let o = SamOpts {
        shuffled: true,
        ..sam_opts()
    };
    let (a, ra) = sam_bytes(&o);
    let (b, rb) = sam_bytes(&o);
    assert_eq!((&a, &ra), (&b, &rb), "same seed, same bytes");
    let text = std::str::from_utf8(&a).unwrap();
    assert!(text.starts_with("@HD\tVN:1.6\tSO:unsorted\n"));
    let keys: Vec<(usize, u64)> = text
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[2][3..].parse().unwrap(), f[3].parse().unwrap())
        })
        .collect();
    assert_eq!(keys.len(), 12000);
    assert!(keys.windows(2).any(|w| w[0] > w[1]), "output is sorted");
    // The reference does not depend on the order; other seeds differ.
    assert_eq!(ra, sam_bytes(&sam_opts()).1);
    assert_ne!(
        ra,
        sam_bytes(&SamOpts {
            seed: 8,
            ..o.clone()
        })
        .1
    );
    // A part of the same reference holds other reads.
    let p1 = sam_bytes(&SamOpts { part: 1, ..o }).0;
    assert_ne!(a, p1);
}

#[test]
fn gen_fasta_has_masks_and_matching_intervals() {
    let o = FastaOpts {
        records: 500,
        chroms: 3,
        chrom_len: 400_000,
        soft_frac: 0.1,
        n_frac: 0.02,
        seed: 3,
        ..FastaOpts::default()
    };
    let (mut fa, mut bed) = (Vec::new(), Vec::new());
    xtask::synth::fasta(&o, &mut fa).unwrap();
    xtask::synth::fasta_intervals(&o, &mut bed).unwrap();
    let refs = read_fasta(&fa);
    assert_eq!(refs.len(), 3);
    let total: usize = refs.values().map(Vec::len).sum();
    assert_eq!(total, 1_200_000);
    let all: Vec<u8> = refs.values().flatten().copied().collect();
    let lower = all.iter().filter(|b| b.is_ascii_lowercase()).count() as f64 / total as f64;
    let n = all.iter().filter(|b| **b == b'N').count() as f64 / total as f64;
    assert!((lower - 0.1).abs() < 0.05, "soft-masked {lower}");
    assert!(n > 0.002 && n < 0.06, "N {n}");
    assert!(
        all.iter().all(|b| b"ACGTNacgt".contains(b)),
        "unexpected base"
    );
    let recs = parse(&bed);
    assert_eq!(recs.len(), 500);
    for (c, s, e) in &recs {
        assert!(s < e && *e as usize <= refs[c].len());
    }
    // Same seed as a SAM's reference: the A/C/G/T agree where not masked.
    let so = SamOpts {
        records: 20_000,
        read_len: 100,
        depth: 10.0,
        chroms: 3,
        seed: 3,
        ..SamOpts::default()
    };
    let len = xtask::synth::sam_chrom_len(&so) as usize;
    let sref = read_fasta(&sam_bytes(&so).1);
    let agree = refs["chr1"][..len.min(400_000)]
        .iter()
        .zip(&sref["chr1"])
        .filter(|(a, b)| a.to_ascii_uppercase() == **b)
        .count();
    assert!(agree as f64 > 0.8 * len.min(400_000) as f64);
}

fn bedpe_text(o: &BedpeOpts) -> String {
    let mut v = Vec::new();
    xtask::synth::bedpe(o, &mut v).unwrap();
    String::from_utf8(v).unwrap()
}

/// One BEDPE line as `BedpeReader` takes it: ten tab-separated fields, ends
/// with start < end inside the chromosome, a score and two strands.
fn check_bedpe_line(l: &str) -> (String, String) {
    let f: Vec<&str> = l.split('\t').collect();
    assert_eq!(f.len(), 10, "{l}");
    for (c, s, e) in [(0, 1, 2), (3, 4, 5)] {
        let (s, e): (u64, u64) = (f[s].parse().unwrap(), f[e].parse().unwrap());
        assert!(f[c].starts_with("chr"), "{l}");
        assert!(s < e && e <= xtask::synth::CHROM_LEN, "{l}");
    }
    let _: u32 = f[7].parse().unwrap();
    assert!(f[8] == "+" || f[8] == "-", "{l}");
    assert!(f[9] == "+" || f[9] == "-", "{l}");
    (f[0].to_string(), f[3].to_string())
}

#[test]
fn gen_bedpe_is_valid() {
    let o = BedpeOpts {
        pairs: 20_000,
        chroms: 5,
        inter_frac: 0.3,
        seed: 4,
        ..BedpeOpts::default()
    };
    let text = bedpe_text(&o);
    assert_eq!(text, bedpe_text(&o), "same seed, same bytes");
    assert_ne!(
        text,
        bedpe_text(&BedpeOpts {
            seed: 5,
            ..o.clone()
        })
    );
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), 20_000);
    let inter = lines
        .iter()
        .filter(|l| {
            let (a, b) = check_bedpe_line(l);
            a != b
        })
        .count();
    let frac = inter as f64 / 20_000.0;
    assert!(
        (frac - 0.3).abs() < 0.02,
        "inter-chromosomal fraction {frac}"
    );
    // On one chromosome end 2 comes after end 1.
    for l in &lines {
        let f: Vec<&str> = l.split('\t').collect();
        if f[0] == f[3] {
            assert!(
                f[2].parse::<u64>().unwrap() < f[4].parse::<u64>().unwrap(),
                "{l}"
            );
        }
    }
    let one = bedpe_text(&BedpeOpts {
        inter_frac: 0.0,
        ..o.clone()
    });
    assert!(one.lines().all(|l| {
        let (a, b) = check_bedpe_line(l);
        a == b
    }));
    assert!(
        xtask::synth::bedpe(&BedpeOpts { chroms: 1, ..o }, &mut Vec::new()).is_err(),
        "an inter-chromosomal fraction needs two chromosomes"
    );
}

#[test]
fn gen_bedpe_b_files() {
    let o = BedpeOpts {
        pairs: 500,
        chroms: 3,
        seed: 9,
        ..BedpeOpts::default()
    };
    // B as pairs: A's stream, every end moved but still overlapping its original.
    let mut b = Vec::new();
    xtask::synth::bedpe_jittered(&o, xtask::synth::B_JITTER, &mut b).unwrap();
    let (a, b) = (bedpe_text(&o), String::from_utf8(b).unwrap());
    assert_ne!(a, b);
    assert_eq!(a.lines().count(), b.lines().count());
    for (la, lb) in a.lines().zip(b.lines()) {
        let (fa, fb): (Vec<&str>, Vec<&str>) = (la.split('\t').collect(), lb.split('\t').collect());
        assert_eq!((fa[0], fa[3], fa[6]), (fb[0], fb[3], fb[6]));
        for (s, e) in [(1, 2), (4, 5)] {
            let n = |f: &[&str], i: usize| f[i].parse::<u64>().unwrap();
            assert!(
                n(&fa, s) < n(&fb, e) && n(&fb, s) < n(&fa, e),
                "{la} / {lb}"
            );
        }
    }
    assert!(xtask::synth::bedpe_jittered(&o, 50, &mut Vec::new()).is_err());
    // B as BED.
    let mut bed = Vec::new();
    xtask::synth::bedpe_bed(&o, 300, &mut bed).unwrap();
    let recs = parse(&bed);
    assert_eq!(recs.len(), 300);
    assert!(recs.iter().all(|(_, s, e)| s < e));
}

#[test]
fn gen_vcf_pair_differs_only_in_genotypes() {
    let o = VcfOpts {
        records: 2000,
        samples: 50,
        seed: 3,
        chroms: 2,
    };
    let (mut a, mut b) = (Vec::new(), Vec::new());
    xtask::synth::vcf_pair(&o, 0.2, &mut a, &mut b).unwrap();
    let mut plain = Vec::new();
    xtask::synth::vcf(&o, &mut plain).unwrap();
    assert_eq!(a, plain, "A is the single-file output");
    let (a, b) = (String::from_utf8(a).unwrap(), String::from_utf8(b).unwrap());
    assert_eq!(a.lines().count(), b.lines().count());
    let (mut calls, mut differ) = (0u32, 0u32);
    for (la, lb) in a.lines().zip(b.lines()) {
        if la.starts_with('#') {
            assert_eq!(la, lb);
            continue;
        }
        let (fa, fb): (Vec<&str>, Vec<&str>) = (la.split('\t').collect(), lb.split('\t').collect());
        assert_eq!(fa[..9], fb[..9], "sites and INFO are the same");
        for (ca, cb) in fa[9..].iter().zip(&fb[9..]) {
            let (ga, gb) = (
                ca.split(':').collect::<Vec<_>>(),
                cb.split(':').collect::<Vec<_>>(),
            );
            assert_eq!(ga[1..], gb[1..], "only the GT field changes");
            if ga[0] != "./." {
                calls += 1;
                differ += u32::from(ga[0] != gb[0]);
            } else {
                assert_eq!(ga[0], gb[0], "missing stays missing");
            }
        }
    }
    // 0.2 redraws, 3/8 of them equal to the original at biallelic sites.
    let frac = f64::from(differ) / f64::from(calls);
    assert!((0.09..0.14).contains(&frac), "{frac}");
    let (mut a2, mut b2) = (Vec::new(), Vec::new());
    xtask::synth::vcf_pair(&o, 0.0, &mut a2, &mut b2).unwrap();
    assert_eq!(a2, b2, "no discordance, same files");
}

#[test]
fn gen_sam_paired_by_name_is_name_grouped() {
    let o = SamOpts {
        records: 4000,
        depth: 5.0,
        chroms: 3,
        paired: true,
        by_name: true,
        seed: 5,
        ..SamOpts::default()
    };
    let (a, ra) = sam_bytes(&o);
    let (b, rb) = sam_bytes(&o);
    assert_eq!((&a, &ra), (&b, &rb), "same seed, same bytes");
    let text = std::str::from_utf8(&a).unwrap();
    assert!(text.starts_with("@HD\tVN:1.6\tSO:queryname\n"));
    let reads: Vec<Vec<&str>> = text
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(|l| l.split('\t').collect())
        .collect();
    assert_eq!(reads.len(), 4000);
    let refs = read_fasta(&ra);
    let mut last = "";
    for pair in reads.chunks(2) {
        let (r1, r2) = (&pair[0], &pair[1]);
        assert_eq!(r1[0], r2[0], "mates share a name");
        assert!(r1[0] > last, "pairs come in name order (zero-padded)");
        last = r1[0];
        let (f1, f2): (u16, u16) = (r1[1].parse().unwrap(), r2[1].parse().unwrap());
        assert_eq!((f1 & 0x1c0, f2 & 0x1c0), (0x41 & 0x1c0, 0x81 & 0x1c0));
        assert_eq!(f1 & 0x1, 1);
        // Each is on the strand the other's flag says, and they differ.
        assert_eq!((f1 & 0x10) != 0, (f2 & 0x20) != 0);
        assert_eq!((f2 & 0x10) != 0, (f1 & 0x20) != 0);
        assert_ne!(f1 & 0x10, f2 & 0x10, "opposite strands");
        assert_eq!((r1[2], r1[6]), (r2[2], "="));
        assert_eq!(r1[7], r2[3], "PNEXT is the mate's POS");
        assert_eq!(r2[7], r1[3]);
        let (t1, t2): (i64, i64) = (r1[8].parse().unwrap(), r2[8].parse().unwrap());
        assert_eq!(t1, -t2);
        // The reads match the reference they are mapped to.
        for r in [r1, r2] {
            let pos: usize = r[3].parse().unwrap();
            let chrom = &refs[r[2]];
            let cig = cigar(r[5]);
            let (mut rp, mut qp) = (pos - 1, 0);
            for (n, op) in cig {
                match op {
                    'M' => {
                        assert_eq!(&chrom[rp..rp + n], r[9][qp..qp + n].as_bytes());
                        rp += n;
                        qp += n;
                    }
                    'N' | 'D' => rp += n,
                    _ => panic!("{op}"),
                }
            }
        }
    }
    // Only name-grouped pairs, and no odd count.
    for bad in [
        SamOpts {
            by_name: false,
            ..o.clone()
        },
        SamOpts {
            paired: false,
            ..o.clone()
        },
        SamOpts {
            records: 4001,
            ..o.clone()
        },
        SamOpts {
            shuffled: true,
            ..o
        },
    ] {
        assert!(xtask::synth::sam(&bad, &mut Vec::new()).is_err());
    }
}

#[test]
fn run_args_writes_the_two_file_setups() {
    let dir = tempfile::tempdir().unwrap();
    let p = |n: &str| dir.path().join(n).to_string_lossy().into_owned();
    let run = |words: &str| {
        let args: Vec<String> = words.split_whitespace().map(String::from).collect();
        xtask::synth::run_args(&args, &mut Vec::new()).unwrap();
    };
    run(&format!(
        "bedpe --pairs 100 --chroms 3 --b-bed 50 --out-a {} --out-b {}",
        p("a.bedpe"),
        p("b.bed")
    ));
    assert_eq!(
        std::fs::read_to_string(p("a.bedpe"))
            .unwrap()
            .lines()
            .count(),
        100
    );
    assert_eq!(
        std::fs::read_to_string(p("b.bed")).unwrap().lines().count(),
        50
    );
    run(&format!(
        "bedpe --pairs 100 --chroms 3 --b-pairs 100 --out-a {} --out-b {}",
        p("a2.bedpe"),
        p("b2.bedpe")
    ));
    assert_eq!(
        std::fs::read_to_string(p("b2.bedpe"))
            .unwrap()
            .lines()
            .count(),
        100
    );
    run(&format!(
        "vcf-pair --records 50 --samples 4 --discord 0.5 --out-a {} --out-b {}",
        p("a.vcf"),
        p("b.vcf")
    ));
    assert_ne!(
        std::fs::read_to_string(p("a.vcf")).unwrap(),
        std::fs::read_to_string(p("b.vcf")).unwrap()
    );
}
