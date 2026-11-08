mod common;
use common::{bcf_records, vcf};
use gloci_core::{ChromId, Interval, RegionSet, Variant, VariantKind, alt_kinds};
use gloci_vcf::{
    BcfSite, MaskRules, PositionSet, SiteFilters, SiteView, Thinner, VcfType, read_ids,
    read_positions, site_stats, text_genotypes,
};
use std::collections::{HashMap, HashSet};

fn span_of(line: &str) -> Interval {
    let f: Vec<&str> = line.split('\t').collect();
    let pos: u64 = f[1].parse().unwrap();
    Interval {
        start: pos - 1,
        end: pos - 1 + f[3].len() as u64,
    }
}

/// The 1-based POS values of `lines` that pass the static filters.
fn kept(f: &SiteFilters, lines: &[&str]) -> Vec<u64> {
    lines
        .iter()
        .filter(|l| {
            let v = SiteView::from_text(l.as_bytes(), span_of(l)).unwrap();
            f.static_pass(&v)
        })
        .map(|l| l.split('\t').nth(1).unwrap().parse().unwrap())
        .collect()
}

fn set(items: &[&str]) -> HashSet<Vec<u8>> {
    items.iter().map(|s| s.as_bytes().to_vec()).collect()
}
fn list(items: &[&str]) -> Vec<Vec<u8>> {
    items.iter().map(|s| s.as_bytes().to_vec()).collect()
}

const R: &[&str] = &[
    "1\t10\trs1\tA\tC\t50\tPASS\tDB;DP=4",
    "1\t20\t.\tA\tC\t.\t.\tDP=4",
    "1\t30\trs3;rs4\tA\tG,T\t10\tLowQual\tDB",
    "1\t40\t.\tAT\tA\t30\tq10;LowQual\tNS=3",
    "2\t50\t.\tA\tAT\t30\tPASS\t.",
    "2\t60\t.\tAT\tGC\t30\tPASS\t.",
    "2\t70\t.\tA\t<DEL>\t30\tPASS\tEND=90",
    "2\t80\t.\tA\t.\t30\tPASS\t.",
];

#[test]
fn pass_only_keeps_dot() {
    let f = SiteFilters {
        pass_only: true,
        ..Default::default()
    };
    assert_eq!(kept(&f, R), vec![10, 20, 50, 60, 70, 80]);
}

#[test]
fn missing_qual_fails_min_qual() {
    let f = SiteFilters {
        min_qual: Some(30.0),
        ..Default::default()
    };
    // POS 20 has QUAL `.`; POS 30 has 10.
    assert_eq!(kept(&f, R), vec![10, 40, 50, 60, 70, 80]);
    let f = SiteFilters {
        min_qual: Some(0.0),
        ..Default::default()
    };
    assert!(!kept(&f, R).contains(&20));
}

#[test]
fn keep_and_remove_filter() {
    let f = SiteFilters {
        keep_filter: Some(list(&["LowQual"])),
        ..Default::default()
    };
    assert_eq!(kept(&f, R), vec![30, 40]);
    let f = SiteFilters {
        remove_filter: Some(list(&["LowQual", "q10"])),
        ..Default::default()
    };
    assert_eq!(kept(&f, R), vec![10, 20, 50, 60, 70, 80]);
    let f = SiteFilters {
        keep_filter: Some(list(&["PASS"])),
        ..Default::default()
    };
    assert_eq!(kept(&f, R), vec![10, 50, 60, 70, 80]);
}

#[test]
fn types_type_set() {
    let ty = |t: &[VcfType]| SiteFilters {
        types: Some(t.to_vec()),
        ..Default::default()
    };
    // 30 is a multi-allelic SNV, 40 a deletion, 50 an insertion, 60 an MNV.
    assert_eq!(kept(&ty(&[VcfType::Snv]), R), vec![10, 20, 30]);
    assert_eq!(kept(&ty(&[VcfType::Indel]), R), vec![40, 50]);
    assert_eq!(kept(&ty(&[VcfType::Ins]), R), vec![50]);
    assert_eq!(kept(&ty(&[VcfType::Del]), R), vec![40]);
    assert_eq!(kept(&ty(&[VcfType::Mnv]), R), vec![60]);
    assert_eq!(kept(&ty(&[VcfType::Sv]), R), vec![70]);
    assert_eq!(kept(&ty(&[VcfType::Ref]), R), vec![80]);
    assert_eq!(kept(&ty(&[VcfType::Mnv, VcfType::Sv]), R), vec![60, 70]);
    let ex = SiteFilters {
        exclude_types: Some(vec![VcfType::Snv, VcfType::Indel]),
        ..Default::default()
    };
    assert_eq!(kept(&ex, R), vec![60, 70, 80]);
    // Type set of a mixed record: SNV + deletion meets either type.
    let mixed = ["1\t5\t.\tAT\tA,GT\t.\t.\t."];
    assert_eq!(kept(&ty(&[VcfType::Snv]), &mixed), vec![5]);
    assert_eq!(kept(&ty(&[VcfType::Del]), &mixed), vec![5]);
    assert!(kept(&ty(&[VcfType::Mnv]), &mixed).is_empty());
}

#[test]
fn type_classifier_matches_core() {
    let alts = [
        "C",
        "AT",
        "<DEL>",
        "<INS:ME>",
        "<*>",
        "<NON_REF>",
        ".",
        "*",
        "A]1:5]",
        "]1:5]A",
        ".A",
        "A.",
        "ATG",
        "GCTA",
        "a",
        "C,<DEL>",
        "T,*",
        "",
    ];
    for refs in ["A", "AT", "ATCG", "at"] {
        for alt in alts {
            let line = format!("1\t10\t.\t{refs}\t{alt}\t.\t.\t.");
            let v = Variant {
                pos: 10,
                reference: refs.as_bytes(),
                alts: alt.as_bytes(),
                end: None,
                svlen: None,
            };
            for (name, vt, kind) in [
                ("snv", VcfType::Snv, VariantKind::Snv),
                ("mnv", VcfType::Mnv, VariantKind::Mnv),
                ("ins", VcfType::Ins, VariantKind::Ins),
                ("del", VcfType::Del, VariantKind::Del),
                ("sv", VcfType::Sv, VariantKind::Sv),
                ("ref", VcfType::Ref, VariantKind::Ref),
                ("other", VcfType::Other, VariantKind::Other),
            ] {
                let f = SiteFilters {
                    types: Some(vec![vt]),
                    ..Default::default()
                };
                let want = alt_kinds(&v).contains(&kind);
                assert_eq!(
                    !kept(&f, &[&line]).is_empty(),
                    want,
                    "{refs} -> {alt:?} as {name}"
                );
            }
        }
    }
}

#[test]
fn alleles_counts() {
    let f = SiteFilters {
        min_alleles: Some(3),
        ..Default::default()
    };
    assert_eq!(kept(&f, R), vec![30]);
    let f = SiteFilters {
        max_alleles: Some(1),
        ..Default::default()
    };
    assert_eq!(kept(&f, R), vec![80]); // ALT `.` gives 1
    let f = SiteFilters {
        min_alleles: Some(2),
        max_alleles: Some(2),
        ..Default::default()
    };
    assert_eq!(kept(&f, R), vec![10, 20, 40, 50, 60, 70]);
}

#[test]
fn info_flags() {
    let f = SiteFilters {
        require_info: list(&["DB"]),
        ..Default::default()
    };
    assert_eq!(kept(&f, R), vec![10, 30]);
    let f = SiteFilters {
        require_info: list(&["DP"]),
        ..Default::default()
    };
    assert_eq!(kept(&f, R), vec![10, 20]);
    let f = SiteFilters {
        exclude_info: list(&["DB", "NS"]),
        ..Default::default()
    };
    assert_eq!(kept(&f, R), vec![20, 50, 60, 70, 80]);
    // A key is not a prefix match: D is neither DB nor DP.
    let f = SiteFilters {
        require_info: list(&["D"]),
        ..Default::default()
    };
    assert!(kept(&f, R).is_empty());
}

#[test]
fn positions_and_ids() {
    let pos: PositionSet = [
        (b"1".to_vec(), 20),
        (b"2".to_vec(), 60),
        (b"2".to_vec(), 20),
    ]
    .into_iter()
    .collect();
    let f = SiteFilters {
        positions: Some(pos.clone()),
        ..Default::default()
    };
    assert_eq!(kept(&f, R), vec![20, 60]);
    let f = SiteFilters {
        exclude_positions: Some(pos),
        ..Default::default()
    };
    assert_eq!(kept(&f, R), vec![10, 30, 40, 50, 70, 80]);
    // An ID cell is split on `;`; any listed ID matches.
    let f = SiteFilters {
        ids: Some(set(&["rs4", "rsX"])),
        ..Default::default()
    };
    assert_eq!(kept(&f, R), vec![30]);
    let f = SiteFilters {
        exclude_ids: Some(set(&["rs1", "rs3"])),
        ..Default::default()
    };
    assert_eq!(kept(&f, R), vec![20, 40, 50, 60, 70, 80]);
    let f = SiteFilters {
        chroms: Some(set(&["2"])),
        ..Default::default()
    };
    assert_eq!(kept(&f, R), vec![50, 60, 70, 80]);
    let f = SiteFilters {
        exclude_chroms: set(&["2"]),
        ..Default::default()
    };
    assert_eq!(kept(&f, R), vec![10, 20, 30, 40]);
}

fn gt_line(pos: u64, gts: &[&str]) -> String {
    format!(
        "1\t{pos}\t.\tA\tC\t.\t.\t.\tGT:DP:GQ:FT\t{}",
        gts.join("\t")
    )
}

fn stats_kept(f: &SiteFilters, lines: &[String]) -> Vec<u64> {
    lines
        .iter()
        .filter(|l| {
            let src = text_genotypes(l.as_bytes()).unwrap();
            let mut masked = Vec::new();
            let st = site_stats(&src, None, &MaskRules::default(), 2, &mut masked).unwrap();
            f.stats_pass(&st)
        })
        .map(|l| l.split('\t').nth(1).unwrap().parse().unwrap())
        .collect()
}

#[test]
fn maf_mac_missing_meandp_hwe_phased() {
    // 10 samples, with one record per scenario.
    fn rep(s: &'static str, n: usize) -> Vec<&'static str> {
        vec![s; n]
    }
    let mk = |pos, parts: Vec<&str>| gt_line(pos, &parts);
    let lines = vec![
        // 10: 20 alleles, 2 alt: MAF 0.1, MAC 2, DP 10, no missing
        mk(
            10,
            [rep("0/0:10:50:PASS", 8), rep("0/1:10:50:PASS", 2)].concat(),
        ),
        // 20: 4 missing of 10; 5 alt of 12 called alleles, MAC 5
        mk(
            20,
            [
                rep("./.:.:.:.", 4),
                rep("0/1:20:50:PASS", 5),
                rep("0/0:20:50:PASS", 1),
            ]
            .concat(),
        ),
        // 30: all het and phased: MAF 0.5, MAC 10
        mk(30, rep("0|1:30:50:PASS", 10)),
        // 40: all missing
        mk(40, rep("./.:.:.:.", 10)),
    ];
    let one = |f: SiteFilters| stats_kept(&f, &lines);
    assert_eq!(
        one(SiteFilters {
            min_maf: Some(0.2),
            ..Default::default()
        }),
        vec![20, 30]
    );
    assert_eq!(
        one(SiteFilters {
            max_maf: Some(0.2),
            ..Default::default()
        }),
        vec![10]
    );
    assert_eq!(
        one(SiteFilters {
            min_mac: Some(5),
            ..Default::default()
        }),
        vec![20, 30]
    );
    assert_eq!(
        one(SiteFilters {
            max_mac: Some(2),
            ..Default::default()
        }),
        vec![10]
    );
    // missing fraction: 0, 0.4, 0, 1
    assert_eq!(
        one(SiteFilters {
            max_missing: Some(0.4),
            ..Default::default()
        }),
        vec![10, 20, 30]
    );
    assert_eq!(
        one(SiteFilters {
            max_missing: Some(0.0),
            ..Default::default()
        }),
        vec![10, 30]
    );
    assert_eq!(
        one(SiteFilters {
            max_missing_count: Some(4),
            ..Default::default()
        }),
        vec![10, 20, 30]
    );
    assert_eq!(
        one(SiteFilters {
            max_missing_count: Some(3),
            ..Default::default()
        }),
        vec![10, 30]
    );
    // mean DP: 10, 20, 30, none
    assert_eq!(
        one(SiteFilters {
            min_mean_dp: Some(15.0),
            ..Default::default()
        }),
        vec![20, 30]
    );
    assert_eq!(
        one(SiteFilters {
            max_mean_dp: Some(20.0),
            ..Default::default()
        }),
        vec![10, 20]
    );
    // HWE: 10 hets out of 10 is far from equilibrium; the rest pass or have no test.
    let hwe = one(SiteFilters {
        hwe: Some(0.01),
        ..Default::default()
    });
    assert!(!hwe.contains(&30), "{hwe:?}");
    assert!(hwe.contains(&10) && hwe.contains(&40));
    // phased-only: a site with no called genotype has no unphased one
    assert_eq!(
        one(SiteFilters {
            phased_only: true,
            ..Default::default()
        }),
        vec![30, 40]
    );
}

#[test]
fn needs_genotypes_only_for_stats_filters() {
    assert!(!SiteFilters::default().needs_genotypes());
    assert!(
        !SiteFilters {
            pass_only: true,
            min_qual: Some(1.0),
            ..Default::default()
        }
        .needs_genotypes()
    );
    for f in [
        SiteFilters {
            min_maf: Some(0.1),
            ..Default::default()
        },
        SiteFilters {
            max_maf: Some(0.1),
            ..Default::default()
        },
        SiteFilters {
            min_mac: Some(1),
            ..Default::default()
        },
        SiteFilters {
            max_mac: Some(1),
            ..Default::default()
        },
        SiteFilters {
            max_missing: Some(0.1),
            ..Default::default()
        },
        SiteFilters {
            max_missing_count: Some(1),
            ..Default::default()
        },
        SiteFilters {
            min_mean_dp: Some(1.0),
            ..Default::default()
        },
        SiteFilters {
            max_mean_dp: Some(1.0),
            ..Default::default()
        },
        SiteFilters {
            hwe: Some(0.1),
            ..Default::default()
        },
        SiteFilters {
            phased_only: true,
            ..Default::default()
        },
    ] {
        assert!(f.needs_genotypes());
    }
}

#[test]
fn thin_first_wins() {
    let f = SiteFilters {
        thin: Some(10),
        ..Default::default()
    };
    let mut t = f.thinner().unwrap();
    let (c0, c1) = (ChromId(0), ChromId(1));
    let got: Vec<bool> = [
        (c0, 100),
        (c0, 105),
        (c0, 110),
        (c0, 111),
        (c0, 120),
        (c1, 121),
        (c1, 122),
        (c0, 130),
    ]
    .into_iter()
    .map(|(c, p)| t.keep(c, p))
    .collect();
    // A dropped site does not move the reference point; a kept one does.
    assert_eq!(got, [true, false, true, false, true, true, false, true]);
    assert!(SiteFilters::default().thinner().is_none());
    let mut all = Thinner::new(0);
    assert!(all.keep(c0, 5) && all.keep(c0, 5));
}

#[test]
fn exclude_regions_points() {
    let mut m = HashMap::new();
    m.insert(
        b"1".to_vec(),
        RegionSet::new(vec![Interval { start: 19, end: 30 }]),
    );
    let lines = [
        "1\t10\t.\tA\tC\t.\t.\t.",
        "1\t20\t.\tA\tC\t.\t.\t.", // SNV [19,20) shares a base with [19,30)
        "1\t25\t.\tA\tC\t.\t.\t.",
        "1\t31\t.\tA\tC\t.\t.\t.",
        "2\t25\t.\tA\tC\t.\t.\t.",
    ];
    let f = SiteFilters {
        exclude_regions: Some(m.clone()),
        ..Default::default()
    };
    assert_eq!(kept(&f, &lines), vec![10, 31, 25]);
    // A point footprint (an insertion `A -> AT` has the point [p, p)) sits on
    // the region edge: kept by default, excluded with inclusive points.
    let pt = ["1\t30\t.\tA\tAT\t.\t.\t.", "1\t25\t.\tA\tAT\t.\t.\t."];
    let point = |l: &str, inc: bool| {
        let f = SiteFilters {
            exclude_regions: Some(m.clone()),
            inclusive_points: inc,
            ..Default::default()
        };
        let v = SiteView::from_text(l.as_bytes(), Interval::point(30)).unwrap();
        f.static_pass(&v)
    };
    assert!(point(pt[0], false));
    assert!(!point(pt[0], true));
    let inside = SiteView::from_text(pt[1].as_bytes(), Interval::point(25)).unwrap();
    let f = SiteFilters {
        exclude_regions: Some(m),
        ..Default::default()
    };
    assert!(!f.static_pass(&inside));
}

#[test]
fn lists_read_from_files() {
    let dir = std::env::temp_dir().join(format!("gloci-lists-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let p = dir.join("pos.txt");
    std::fs::write(&p, "# header\n1\t100\n\n1\t200\textra\nchrX\t5\n").unwrap();
    let got = read_positions(p.to_str().unwrap()).unwrap();
    let want: PositionSet = [
        (b"1".to_vec(), 100),
        (b"1".to_vec(), 200),
        (b"chrX".to_vec(), 5),
    ]
    .into_iter()
    .collect();
    assert_eq!(got, want);
    std::fs::write(&p, "1\tabc\n").unwrap();
    let e = read_positions(p.to_str().unwrap()).unwrap_err().to_string();
    assert!(e.contains("line 1") && e.contains("abc"), "{e}");
    std::fs::write(&p, "1\t0\n").unwrap();
    assert!(read_positions(p.to_str().unwrap()).is_err());
    std::fs::write(&p, "1\n").unwrap();
    assert!(read_positions(p.to_str().unwrap()).is_err());
    std::fs::write(&p, "# ids\nrs1\n\nrs2\r\nrs1\n").unwrap();
    assert_eq!(read_ids(p.to_str().unwrap()).unwrap(), set(&["rs1", "rs2"]));
    assert!(read_ids(dir.join("nope").to_str().unwrap()).is_err());
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn text_view_errors() {
    assert!(SiteView::from_text(b"1\t10\t.\tA\tC", Interval::point(9)).is_err());
    assert!(SiteView::from_text(b"1\tx\t.\tA\tC\t.\t.\t.", Interval::point(9)).is_err());
    assert!(SiteView::from_text(b"1\t10\t.\tA\tC\tzz\t.\t.", Interval::point(9)).is_err());
    let v = SiteView::from_text(
        b"1\t10\trs\tA\tC\t30.5\tPASS\tDP=1\tGT\t0/1",
        Interval::point(9),
    )
    .unwrap();
    assert_eq!(
        (v.chrom, v.pos, v.id, v.ref_, v.alts),
        (&b"1"[..], 10, &b"rs"[..], &b"A"[..], &b"C"[..])
    );
    assert_eq!(
        (v.qual, v.filter, v.info),
        (Some(30.5), &b"PASS"[..], &b"DP=1"[..])
    );
}

#[test]
fn bcf_view_equals_text_view() {
    let recs = [
        "1\t10\trs1\tA\tC\t50\tPASS\tDP=4",
        "1\t20\t.\tA\tC,G\t.\t.\tDP=4",
        "1\t30\trs3\tAT\tA\t0.7\tLowQual\tDP=2",
        "1\t40\t.\tA\t.\t30\tLowQual;PASS\t.",
    ];
    let text = vcf(
        &["s1"],
        &recs
            .map(|r| format!("{r}\tGT\t0/1"))
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
    );
    let bcf = gloci_testkit::fixtures::vcf_to_bcf(text.as_bytes());
    let (dicts, raws) = bcf_records(&bcf);
    assert_eq!(raws.len(), 4);
    let mut site = BcfSite::default();
    for (raw, line) in raws.iter().zip(recs) {
        let t = SiteView::from_text(line.as_bytes(), span_of(line)).unwrap();
        let b = site.view(raw, &dicts, span_of(line)).unwrap();
        assert_eq!(
            (b.chrom, b.pos, b.id, b.ref_, b.alts),
            (t.chrom, t.pos, t.id, t.ref_, t.alts),
            "{line}"
        );
        assert_eq!(b.qual, t.qual, "{line}");
        let norm = |f: &[u8]| {
            if f.is_empty() {
                b".".to_vec()
            } else {
                f.to_vec()
            }
        };
        assert_eq!(norm(b.filter), norm(t.filter), "{line}");
        // The BCF view's info holds the keys of the record.
        let keys = |i: &[u8]| -> Vec<Vec<u8>> {
            i.split(|&c| c == b';')
                .filter(|k| *k != b"." && !k.is_empty())
                .map(|k| k.split(|&c| c == b'=').next().unwrap().to_vec())
                .collect()
        };
        assert_eq!(keys(b.info), keys(t.info), "{line}");
    }
    // 0.7 is not exactly representable in f32; it must still meet --min-qual 0.7.
    let f = SiteFilters {
        min_qual: Some(0.7),
        ..Default::default()
    };
    let v = site.view(&raws[2], &dicts, Interval::point(29)).unwrap();
    assert!(f.static_pass(&v));
}
