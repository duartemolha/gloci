//! `flagstat` (D128): the loop's counters, samtools' f32 percentages and the three layouts.
use gloci_ops::{FlagStats, FlagstatFormat, percent, render_flagstat};

/// Counter names in the order of `FlagStats`' fields, for the table below.
const NAMES: [&str; 16] = [
    "total",
    "primary",
    "secondary",
    "supplementary",
    "duplicates",
    "primary_duplicates",
    "mapped",
    "primary_mapped",
    "paired",
    "read1",
    "read2",
    "proper_pair",
    "both_mapped",
    "singletons",
    "mate_other_ref",
    "mate_other_ref_q5",
];

fn counters(s: &FlagStats) -> [[u64; 2]; 16] {
    [
        s.total,
        s.primary,
        s.secondary,
        s.supplementary,
        s.duplicates,
        s.primary_duplicates,
        s.mapped,
        s.primary_mapped,
        s.paired,
        s.read1,
        s.read2,
        s.proper_pair,
        s.both_mapped,
        s.singletons,
        s.mate_other_ref,
        s.mate_other_ref_q5,
    ]
}

/// The counters that moved, as (name, column) pairs.
fn moved(flag: u16, r: Option<usize>, m: Option<usize>, q: u8) -> Vec<(&'static str, usize)> {
    let mut s = FlagStats::default();
    s.add(flag, r, m, q);
    let mut v = Vec::new();
    for (i, c) in counters(&s).iter().enumerate() {
        for (w, n) in c.iter().enumerate() {
            if *n > 0 {
                assert_eq!(*n, 1);
                v.push((NAMES[i], w));
            }
        }
    }
    v
}

fn pass(names: &[&'static str]) -> Vec<(&'static str, usize)> {
    names.iter().map(|n| (*n, 0)).collect()
}

#[test]
fn flagstat_loop_table() {
    // Primary, paired, proper, read1, mate reverse, both mapped on one reference.
    assert_eq!(
        moved(0x41 | 0x2 | 0x20, Some(0), Some(0), 60),
        pass(&[
            "total",
            "primary",
            "mapped",
            "primary_mapped",
            "paired",
            "read1",
            "proper_pair",
            "both_mapped"
        ])
    );
    // A duplicate secondary counts as a duplicate, not as a primary duplicate.
    assert_eq!(
        moved(0x100 | 0x400, Some(0), None, 0),
        pass(&["total", "secondary", "duplicates", "mapped"])
    );
    assert_eq!(
        moved(0x800, Some(0), None, 0),
        pass(&["total", "supplementary", "mapped"])
    );
    // A secondary read flagged paired is not counted as paired.
    assert_eq!(
        moved(0x100 | 0x1 | 0x2, Some(0), Some(1), 60),
        pass(&["total", "secondary", "mapped"])
    );
    // A primary duplicate.
    assert_eq!(
        moved(0x400, Some(0), None, 0),
        pass(&["total", "primary", "duplicates", "primary_duplicates", "mapped", "primary_mapped"])
    );
    // Mapped, mate unmapped: a singleton, and not "both mapped".
    assert_eq!(
        moved(0x1 | 0x8, Some(0), Some(0), 60),
        pass(&["total", "primary", "mapped", "primary_mapped", "paired", "singletons"])
    );
    // Unmapped with 0x2 set: paired, not properly paired, not a singleton.
    assert_eq!(
        moved(0x1 | 0x4 | 0x2, None, None, 0),
        pass(&["total", "primary", "paired"])
    );
    // Both unmapped: paired only.
    assert_eq!(
        moved(0x1 | 0x4 | 0x8 | 0x80, None, None, 0),
        pass(&["total", "primary", "paired", "read2"])
    );
    // Mate on another reference: MAPQ 4 is below the threshold, 5 is not.
    let both = ["total", "primary", "mapped", "primary_mapped", "paired", "both_mapped"];
    let mut q4 = both.to_vec();
    q4.push("mate_other_ref");
    assert_eq!(moved(0x1, Some(0), Some(1), 4), pass(&q4));
    q4.push("mate_other_ref_q5");
    assert_eq!(moved(0x1, Some(0), Some(1), 5), pass(&q4));
    // A mate with no reference is another reference.
    assert_eq!(moved(0x1, Some(0), None, 4), pass(&q4[..7]));
    // QC-failed reads move column 1 only.
    let v = moved(0x200 | 0x1 | 0x2 | 0x40, Some(0), Some(0), 60);
    assert!(!v.is_empty() && v.iter().all(|(_, w)| *w == 1), "{v:?}");
    assert_eq!(v.len(), 8);
}

#[test]
fn percent_f32() {
    // f64 gives 95.62 and 4.38 here; samtools' f32 quotient gives 95.63 and 4.37.
    assert_eq!(percent(153, 160, false), "95.63%");
    assert_eq!(percent(7, 160, false), "4.37%");
    assert_eq!(percent(1, 160, false), "0.63%");
    assert_eq!(percent(0, 0, false), "N/A");
    assert_eq!(percent(0, 0, true), "null");
    assert_eq!(percent(153, 160, true), "95.63");
    assert_eq!(percent(5, 5, false), "100.00%");
    assert_eq!(percent(0, 7, true), "0.00");
    // The f64 route the formula avoids.
    assert_eq!(format!("{:.2}", 153f64 / 160f64 * 100.0), "95.62");
}

fn sample() -> FlagStats {
    let mut s = FlagStats::default();
    s.add(0x41 | 0x2, Some(0), Some(0), 60);
    s.add(0x4, None, None, 0);
    s.add(0x200 | 0x10, Some(1), None, 3);
    s
}

#[test]
fn render_layouts() {
    let s = sample();
    let text = render_flagstat(&s, FlagstatFormat::Text);
    assert!(text.starts_with("2 + 1 in total (QC-passed reads + QC-failed reads)\n"));
    assert!(text.contains("\n1 + 1 mapped (50.00% : 100.00%)\n"));
    assert!(text.contains("\n1 + 0 properly paired (100.00% : N/A)\n"));
    assert!(text.ends_with("0 + 0 with mate mapped to a different chr (mapQ>=5)\n"));
    let tsv = render_flagstat(&s, FlagstatFormat::Tsv);
    assert!(tsv.starts_with("2\t1\ttotal (QC-passed reads + QC-failed reads)\n"));
    assert!(tsv.contains("\n50.00%\t100.00%\tmapped %\n"));
    let json = render_flagstat(&s, FlagstatFormat::Json);
    assert!(json.starts_with("{\n \"QC-passed reads\": { \n  \"total\": 2, \n"));
    assert!(json.contains("  \"mapped %\": 50.00, \n"));
    assert!(json.contains("  \"properly paired %\": null, \n"));
    assert!(json.contains(" },\n \"QC-failed reads\": { \n"));
    assert!(json.ends_with(
        "  \"with mate mapped to a different chr (mapQ >= 5)\": 0 \n }\n}\n"
    ));
}
