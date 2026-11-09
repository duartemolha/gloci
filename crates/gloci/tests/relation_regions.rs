//! `--region`/`--regions-file` on the relation commands (D75): A is
//! restricted to the records that overlap a region, B is not, and the output
//! is the unrestricted run's output for those A records, at every thread
//! count, for streaming and indexed A alike.
mod common;
use common::*;
use predicates::str::contains;
use std::path::Path;

fn s(p: &Path) -> String {
    p.to_string_lossy().into_owned()
}

fn run(args: &[String]) -> Vec<u8> {
    let o = gloci().args(args).assert().success();
    o.get_output().stdout.clone()
}

fn strs(a: &[&str]) -> Vec<String> {
    a.iter().map(|x| x.to_string()).collect()
}

/// A: records around the edges of `chr2:100-500` (0-based [99, 500)), a
/// point on each edge and one inside, records at a chromosome start, and
/// records on other chromosomes.
fn a_rows() -> Vec<(&'static str, u64, u64, String)> {
    let mut v: Vec<(&str, u64, u64, String)> = Vec::new();
    for (i, (st, en)) in [(0, 40), (10, 60), (70, 90), (500, 520)].iter().enumerate() {
        v.push(("chr1", *st, *en, format!("a1_{i}")));
    }
    for (i, (st, en)) in [
        (10, 99),   // ends where the region starts: no overlap
        (50, 100),  // overlaps by one base
        (99, 99),   // point on the start edge
        (100, 100), // point inside
        (120, 130),
        (400, 900), // reaches far past the region
        (500, 500), // point on the end edge
        (500, 600), // starts where the region ends: no overlap
        (700, 800),
    ]
    .iter()
    .enumerate()
    {
        v.push(("chr2", *st, *en, format!("a2_{i}")));
    }
    for (i, (st, en)) in [(5, 25), (30, 70), (90, 150)].iter().enumerate() {
        v.push(("chr3", *st, *en, format!("a3_{i}")));
    }
    v
}

fn bed(rows: &[(&str, u64, u64, String)]) -> String {
    rows.iter()
        .map(|(c, st, en, n)| format!("{c}\t{st}\t{en}\t{n}\t0\t+\n"))
        .collect()
}

/// B reaches outside every region used here, including a record on chr2
/// far from the region, which `closest` must find.
fn b_text() -> String {
    let mut v: Vec<(&str, u64, u64, String)> = Vec::new();
    for (i, (st, en)) in [(20, 30), (65, 80), (600, 650)].iter().enumerate() {
        v.push(("chr1", *st, *en, format!("b1_{i}")));
    }
    for (i, (st, en)) in [(60, 70), (450, 460), (850, 950), (2000, 2100)]
        .iter()
        .enumerate()
    {
        v.push(("chr2", *st, *en, format!("b2_{i}")));
    }
    for (i, (st, en)) in [(40, 60), (200, 300)].iter().enumerate() {
        v.push(("chr3", *st, *en, format!("b3_{i}")));
    }
    bed(&v)
}

/// The overlap rule of `view --region` (4.3): a point counts only strictly
/// inside a region, or on its edge with `inclusive`.
fn overlaps(st: u64, en: u64, (rs, re): (u64, u64), inclusive: bool) -> bool {
    if st == en {
        if inclusive {
            rs <= st && st <= re
        } else {
            rs < st && st < re
        }
    } else {
        st < re && rs < en
    }
}

struct Fx {
    dir: tempfile::TempDir,
    a_plain: String,
    a_idx: String,
    b_idx: String,
    b_plain: String,
}

impl Fx {
    fn new() -> Fx {
        let dir = tempdir();
        let a = bed(&a_rows());
        let a_plain = s(&plain_file(dir.path(), "a.bed", &a));
        let a_idx = s(&bgzf_file(dir.path(), "a.bed.gz", &a));
        let b_idx = s(&bgzf_file(dir.path(), "b.bed.gz", &b_text()));
        let b_plain = s(&plain_file(dir.path(), "b.bed", &b_text()));
        gloci().args(["index", &a_idx]).assert().success();
        gloci().args(["index", &b_idx]).assert().success();
        Fx {
            dir,
            a_plain,
            a_idx,
            b_idx,
            b_plain,
        }
    }

    /// A file holding only the A records that overlap `regions`.
    fn filtered_a(&self, regions: &[(&str, u64, u64)], inclusive: bool) -> String {
        let rows: Vec<_> = a_rows()
            .into_iter()
            .filter(|(c, st, en, _)| {
                regions
                    .iter()
                    .any(|(rc, rs, re)| rc == c && overlaps(*st, *en, (*rs, *re), inclusive))
            })
            .collect();
        s(&plain_file(self.dir.path(), "a_filtered.bed", &bed(&rows)))
    }

    fn regions_file(&self, regions: &[(&str, u64, u64)]) -> String {
        let t: String = regions
            .iter()
            .map(|(c, a, b)| format!("{c}\t{a}\t{b}\n"))
            .collect();
        s(&plain_file(self.dir.path(), "regions.bed", &t))
    }
}

const COMMANDS: [&str; 5] = ["intersect", "closest", "map", "coverage", "subtract"];

fn cmd(name: &str, a: &str, b: &str, t: &str) -> Vec<String> {
    let mut v = strs(&[name, "-a", a, "-b", b, "-t", t]);
    match name {
        "intersect" => v.extend(strs(&["--mode", "pairs"])),
        "map" => v.extend(strs(&["-c", "4:count"])),
        _ => {}
    }
    v
}

const REGION_SETS: &[&[(&str, u64, u64)]] = &[
    &[("chr2", 99, 500)],
    // at a chromosome start
    &[("chr1", 0, 50)],
    // past the chromosome end (no record reaches it)
    &[("chr3", 100_000, 200_000)],
    // an absent chromosome
    &[("chrZ", 0, 100)],
    // several, on several chromosomes, overlapping each other
    &[
        ("chr1", 0, 50),
        ("chr2", 99, 500),
        ("chr2", 450, 800),
        ("chr3", 0, 60),
    ],
];

#[test]
fn region_restricts_a() {
    let fx = Fx::new();
    for inclusive in [false, true] {
        for regions in REGION_SETS {
            let rf = fx.regions_file(regions);
            let a_filtered = fx.filtered_a(regions, inclusive);
            for name in COMMANDS {
                for (a, b) in [(&fx.a_plain, &fx.b_plain), (&fx.a_idx, &fx.b_idx)] {
                    for t in ["1", "4"] {
                        let mut args = cmd(name, a, b, t);
                        args.extend(strs(&["--regions-file", &rf]));
                        if inclusive {
                            args.push("--inclusive-points".into());
                        }
                        let mut want_args = cmd(name, &a_filtered, b, t);
                        if inclusive {
                            want_args.push("--inclusive-points".into());
                        }
                        assert_eq!(
                            String::from_utf8_lossy(&run(&args)),
                            String::from_utf8_lossy(&run(&want_args)),
                            "{name} {regions:?} inclusive={inclusive} t={t} a={a}"
                        );
                    }
                }
            }
        }
    }
}

/// The same, as the full run's rows filtered by their A columns.
#[test]
fn region_output_is_the_full_runs_rows_for_the_selected_a() {
    let fx = Fx::new();
    for name in ["intersect", "closest", "coverage", "map"] {
        for regions in REGION_SETS {
            let full = String::from_utf8(run(&cmd(name, &fx.a_idx, &fx.b_idx, "1"))).unwrap();
            let want: String = full
                .lines()
                .filter(|l| {
                    let f: Vec<&str> = l.split('\t').collect();
                    let (c, st, en) = (f[0], f[1].parse().unwrap(), f[2].parse().unwrap());
                    regions
                        .iter()
                        .any(|(rc, rs, re)| *rc == c && overlaps(st, en, (*rs, *re), false))
                })
                .map(|l| format!("{l}\n"))
                .collect();
            let mut args = cmd(name, &fx.a_idx, &fx.b_idx, "4");
            args.extend(strs(&["--regions-file", &fx.regions_file(regions)]));
            assert_eq!(
                String::from_utf8(run(&args)).unwrap(),
                want,
                "{name} {regions:?}"
            );
        }
    }
}

#[test]
fn region_flag_parses_one_based_closed() {
    let fx = Fx::new();
    let a_filtered = fx.filtered_a(&[("chr2", 99, 500)], false);
    for a in [&fx.a_plain, &fx.a_idx] {
        let mut args = cmd("coverage", a, &fx.b_idx, "1");
        args.extend(strs(&["--region", "chr2:100-500"]));
        let want = run(&cmd("coverage", &a_filtered, &fx.b_idx, "1"));
        assert_eq!(run(&args), want);
    }
}

#[test]
fn closest_still_finds_b_outside_the_region() {
    let fx = Fx::new();
    // a2_3 [100,100) is a point inside the region [99,101); its nearest B is
    // b2_0 [60,70), which lies outside the region.
    for (a, b) in [(&fx.a_plain, &fx.b_plain), (&fx.a_idx, &fx.b_idx)] {
        for t in ["1", "4"] {
            let mut args = cmd("closest", a, b, t);
            args.extend(strs(&["--region", "chr2:100-101"]));
            let out = String::from_utf8(run(&args)).unwrap();
            assert!(out.contains("a2_3"), "{out}");
            assert!(out.contains("b2_0"), "{out}");
            // a2_5 [400,900) is selected and meets b2_2 [850,950), which
            // starts outside the region [99,101).
            let mut args = cmd("intersect", a, b, t);
            args.extend(strs(&["--region", "chr2:450-460"]));
            let out = String::from_utf8(run(&args)).unwrap();
            assert!(out.contains("a2_5") && out.contains("b2_2"), "{out}");
        }
    }
}

#[test]
fn inclusive_points_at_region_edges() {
    let fx = Fx::new();
    // Region [99, 500): points at 99 (start edge) and 500 (end edge).
    for a in [&fx.a_plain, &fx.a_idx] {
        let names = |inclusive: bool, t: &str| -> Vec<String> {
            let mut args = cmd("coverage", a, &fx.b_idx, t);
            args.extend(strs(&["--region", "chr2:100-500"]));
            if inclusive {
                args.push("--inclusive-points".into());
            }
            String::from_utf8(run(&args))
                .unwrap()
                .lines()
                .map(|l| l.split('\t').nth(3).unwrap().to_string())
                .collect()
        };
        for t in ["1", "4"] {
            let strict = names(false, t);
            assert!(!strict.contains(&"a2_2".to_string()), "{strict:?}");
            assert!(!strict.contains(&"a2_6".to_string()), "{strict:?}");
            assert!(strict.contains(&"a2_3".to_string()), "{strict:?}");
            let incl = names(true, t);
            assert!(incl.contains(&"a2_2".to_string()), "{incl:?}");
            assert!(incl.contains(&"a2_6".to_string()), "{incl:?}");
        }
    }
}

#[test]
fn output_is_identical_at_every_thread_count() {
    let fx = Fx::new();
    for name in COMMANDS {
        for regions in REGION_SETS {
            let rf = fx.regions_file(regions);
            let mut outs = Vec::new();
            for (a, b) in [(&fx.a_plain, &fx.b_plain), (&fx.a_idx, &fx.b_idx)] {
                for t in ["1", "2", "3", "8"] {
                    let mut args = cmd(name, a, b, t);
                    args.extend(strs(&["--regions-file", &rf]));
                    outs.push((format!("{a} t={t}"), run(&args)));
                }
            }
            for (label, o) in &outs[1..] {
                assert_eq!(o, &outs[0].1, "{name} {regions:?} {label} vs {}", outs[0].0);
            }
        }
    }
}

#[test]
fn indexed_a_is_read_through_its_index() {
    let fx = Fx::new();
    for t in ["1", "4"] {
        let mut args = cmd("intersect", &fx.a_idx, &fx.b_idx, t);
        args.extend(strs(&["--region", "chr2:100-500", "-vv"]));
        gloci()
            .args(&args)
            .assert()
            .success()
            .stderr(contains("through its index"));
    }
    // a plain file is scanned, and view's warning says so
    let mut args = cmd("intersect", &fx.a_plain, &fx.b_plain, "1");
    args.extend(strs(&["--region", "chr2:100-500"]));
    gloci()
        .args(&args)
        .assert()
        .success()
        .stderr(contains("scanning the whole file"));
}

#[test]
fn absent_chromosome_warns_as_view_does() {
    let fx = Fx::new();
    for (a, b) in [(&fx.a_plain, &fx.b_plain), (&fx.a_idx, &fx.b_idx)] {
        for t in ["1", "4"] {
            let mut args = cmd("coverage", a, b, t);
            args.extend(strs(&["--region", "chrZ:1-100"]));
            gloci()
                .args(&args)
                .assert()
                .success()
                .stdout("")
                .stderr(contains("chromosome \"chrZ\" not found"));
        }
    }
}

#[test]
fn unsorted_and_sort_keep_working_with_regions() {
    let fx = Fx::new();
    // an unsorted A: chr2 rows before chr1 rows, out of order within chr2
    let mut rows = a_rows();
    rows.reverse();
    let shuffled = s(&plain_file(fx.dir.path(), "shuffled.bed", &bed(&rows)));
    let regions = [("chr2", 99, 500)];
    let rf = fx.regions_file(&regions);
    let filtered_rows: Vec<_> = rows
        .iter()
        .filter(|(c, st, en, _)| *c == "chr2" && overlaps(*st, *en, (99, 500), false))
        .cloned()
        .collect();
    let filtered = s(&plain_file(fx.dir.path(), "f.bed", &bed(&filtered_rows)));
    for flag in ["--unsorted", "--sort"] {
        let mut args = cmd("coverage", &shuffled, &fx.b_plain, "1");
        args.extend(strs(&["--regions-file", &rf, flag]));
        let mut want = cmd("coverage", &filtered, &fx.b_plain, "1");
        want.push(flag.into());
        assert_eq!(run(&args), run(&want), "{flag}");
    }
}

fn map_file(fx: &Fx) -> String {
    s(&plain_file(fx.dir.path(), "map.tsv", "1\tchr1\n2\tchr2\n"))
}

#[test]
fn region_names_are_in_the_renamed_space() {
    let fx = Fx::new();
    let map = map_file(&fx);
    let a = bed(&[
        ("1", 10, 20, "x".into()),
        ("1", 100, 200, "y".into()),
        ("2", 5, 9, "z".into()),
    ]);
    let b = bed(&[("1", 15, 30, "p".into()), ("2", 1, 3, "q".into())]);
    let a_idx = s(&bgzf_file(fx.dir.path(), "m_a.bed.gz", &a));
    let b_idx = s(&bgzf_file(fx.dir.path(), "m_b.bed.gz", &b));
    gloci().args(["index", &a_idx]).assert().success();
    gloci().args(["index", &b_idx]).assert().success();
    for t in ["1", "4"] {
        let o = run(&strs(&[
            "intersect",
            "-a",
            &a_idx,
            "-b",
            &b_idx,
            "-t",
            t,
            "--chrom-map",
            &map,
            "--region",
            "chr1:1-50",
        ]));
        assert_eq!(
            String::from_utf8(o).unwrap(),
            "chr1\t10\t20\tx\t0\t+\n",
            "t={t}"
        );
        // the old name is not a region of the renamed records
        gloci()
            .args([
                "intersect",
                "-a",
                &a_idx,
                "-b",
                &b_idx,
                "-t",
                t,
                "--chrom-map",
                &map,
            ])
            .args(["--region", "1:1-50"])
            .assert()
            .success()
            .stdout("")
            .stderr(contains("not found"));
    }
}

#[test]
fn indexed_collision_under_a_map_is_refused() {
    let fx = Fx::new();
    let map = map_file(&fx);
    let both = bed(&[("1", 10, 20, "x".into()), ("chr1", 30, 40, "y".into())]);
    let a_idx = s(&bgzf_file(fx.dir.path(), "both.bed.gz", &both));
    gloci().args(["index", &a_idx]).assert().success();
    gloci()
        .args([
            "intersect",
            "-a",
            &a_idx,
            "-b",
            &fx.b_idx,
            "-t",
            "2",
            "--chrom-map",
            &map,
        ])
        .assert()
        .code(2)
        .stderr(contains("both.bed.gz has both \"1\" and \"chr1\""));
}

#[test]
fn indexed_parse_error_names_the_renamed_chromosome() {
    let fx = Fx::new();
    let map = map_file(&fx);
    let good = bed(&[
        ("1", 10, 50, "a".into()),
        ("2", 10, 50, "b".into()),
        ("3", 10, 50, "c".into()),
    ]);
    let dir = fx.dir.path().join("bad");
    std::fs::create_dir(&dir).unwrap();
    let a = s(&bgzf_file(&dir, "bad.bed.gz", &good));
    gloci().args(["index", &a]).assert().success();
    // corrupt chr 2 and keep the index
    let bad = good.replace("2\t10\t50", "2\tNOT_A_NUMBER\t50");
    let _ = bgzf_file(&dir, "bad.bed.gz", &bad);
    let o = gloci()
        .args([
            "intersect",
            "-a",
            &a,
            "-b",
            &fx.b_idx,
            "-t",
            "2",
            "--chrom-map",
            &map,
        ])
        .output()
        .unwrap();
    assert!(!o.status.success());
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("record 1 of chr2"), "{err}");
}
