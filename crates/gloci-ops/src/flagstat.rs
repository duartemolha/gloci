//! `gloci flagstat` (D128): samtools' flag counts, in its text, TSV and JSON layouts.

use gloci_io::{AlnSource, IoError, Loc, aln_ref};

/// The counters of samtools' `flagstat_loop`: index 0 holds the QC-passed reads, index 1 the
/// QC-failed ones (FLAGS bit 0x200).
#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct FlagStats {
    pub total: [u64; 2],
    pub primary: [u64; 2],
    pub secondary: [u64; 2],
    pub supplementary: [u64; 2],
    pub duplicates: [u64; 2],
    pub primary_duplicates: [u64; 2],
    pub mapped: [u64; 2],
    pub primary_mapped: [u64; 2],
    pub paired: [u64; 2],
    pub read1: [u64; 2],
    pub read2: [u64; 2],
    pub proper_pair: [u64; 2],
    pub both_mapped: [u64; 2],
    pub singletons: [u64; 2],
    pub mate_other_ref: [u64; 2],
    pub mate_other_ref_q5: [u64; 2],
}

impl FlagStats {
    /// Count one record.
    pub fn add(&mut self, flag: u16, ref_id: Option<usize>, mate_ref_id: Option<usize>, mapq: u8) {
        let w = (flag & 0x200 != 0) as usize;
        self.total[w] += 1;
        if flag & 0x100 != 0 {
            self.secondary[w] += 1;
        } else if flag & 0x800 != 0 {
            self.supplementary[w] += 1;
        } else {
            self.primary[w] += 1;
            if flag & 0x1 != 0 {
                self.paired[w] += 1;
                if flag & 0x2 != 0 && flag & 0x4 == 0 {
                    self.proper_pair[w] += 1;
                }
                if flag & 0x40 != 0 {
                    self.read1[w] += 1;
                }
                if flag & 0x80 != 0 {
                    self.read2[w] += 1;
                }
                if flag & 0x8 != 0 && flag & 0x4 == 0 {
                    self.singletons[w] += 1;
                }
                if flag & 0xC == 0 {
                    self.both_mapped[w] += 1;
                    if mate_ref_id != ref_id {
                        self.mate_other_ref[w] += 1;
                        if mapq >= 5 {
                            self.mate_other_ref_q5[w] += 1;
                        }
                    }
                }
            }
            if flag & 0x4 == 0 {
                self.primary_mapped[w] += 1;
            }
            if flag & 0x400 != 0 {
                self.primary_duplicates[w] += 1;
            }
        }
        if flag & 0x4 == 0 {
            self.mapped[w] += 1;
        }
        if flag & 0x400 != 0 {
            self.duplicates[w] += 1;
        }
    }
}

/// The three layouts of `flagstat`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlagstatFormat {
    Text,
    Tsv,
    Json,
}

impl FlagstatFormat {
    /// The names `-O` accepts.
    pub const NAMES: &'static str = "text, tsv, json";

    /// The format a `-O` value names (case-insensitive).
    pub fn from_name(s: &str) -> Option<Self> {
        match s.to_ascii_lowercase().as_str() {
            "text" => Some(Self::Text),
            "tsv" => Some(Self::Tsv),
            "json" => Some(Self::Json),
            _ => None,
        }
    }
}

/// `n` as a percentage of `total`, as samtools computes it: a 32-bit float quotient, widened
/// and multiplied by 100 in double precision, printed with two decimals. A total of 0 gives
/// `N/A` (`null` in JSON); text and TSV add a `%` sign. IEEE division is exact in both
/// precisions and `{:.2}` rounds the exact binary value, as C's `%.2f` does, so the result
/// does not depend on the CPU.
pub fn percent(n: u64, total: u64, json: bool) -> String {
    if total == 0 {
        return if json { "null" } else { "N/A" }.to_string();
    }
    let p = (n as f32 / total as f32) as f64 * 100.0;
    if json {
        format!("{p:.2}")
    } else {
        format!("{p:.2}%")
    }
}

/// One row of the report: its label, the counter, and the counter its percentage is of.
type Row = (
    &'static str,
    fn(&FlagStats) -> &[u64; 2],
    Option<fn(&FlagStats) -> &[u64; 2]>,
);

const ROWS: [Row; 16] = [
    ("total", |s| &s.total, None),
    ("primary", |s| &s.primary, None),
    ("secondary", |s| &s.secondary, None),
    ("supplementary", |s| &s.supplementary, None),
    ("duplicates", |s| &s.duplicates, None),
    ("primary duplicates", |s| &s.primary_duplicates, None),
    ("mapped", |s| &s.mapped, Some(|s| &s.total)),
    (
        "primary mapped",
        |s| &s.primary_mapped,
        Some(|s| &s.primary),
    ),
    ("paired in sequencing", |s| &s.paired, None),
    ("read1", |s| &s.read1, None),
    ("read2", |s| &s.read2, None),
    ("properly paired", |s| &s.proper_pair, Some(|s| &s.paired)),
    ("with itself and mate mapped", |s| &s.both_mapped, None),
    ("singletons", |s| &s.singletons, Some(|s| &s.paired)),
    (
        "with mate mapped to a different chr",
        |s| &s.mate_other_ref,
        None,
    ),
    (
        "with mate mapped to a different chr (mapQ>=5)",
        |s| &s.mate_other_ref_q5,
        None,
    ),
];

/// The report, byte for byte as `samtools flagstat -O text|tsv|json` writes it.
pub fn render_flagstat(s: &FlagStats, f: FlagstatFormat) -> String {
    use std::fmt::Write;
    let mut o = String::new();
    match f {
        FlagstatFormat::Text => {
            for (label, get, base) in ROWS {
                let a = get(s);
                let label = if label == "total" {
                    "in total (QC-passed reads + QC-failed reads)"
                } else {
                    label
                };
                let _ = write!(o, "{} + {} {label}", a[0], a[1]);
                if let Some(base) = base {
                    let t = base(s);
                    let _ = write!(
                        o,
                        " ({} : {})",
                        percent(a[0], t[0], false),
                        percent(a[1], t[1], false)
                    );
                }
                o.push('\n');
            }
        }
        FlagstatFormat::Tsv => {
            for (label, get, base) in ROWS {
                let a = get(s);
                let label = if label == "total" {
                    "total (QC-passed reads + QC-failed reads)"
                } else {
                    label
                };
                let _ = writeln!(o, "{}\t{}\t{label}", a[0], a[1]);
                if let Some(base) = base {
                    let t = base(s);
                    let _ = writeln!(
                        o,
                        "{}\t{}\t{label} %",
                        percent(a[0], t[0], false),
                        percent(a[1], t[1], false)
                    );
                }
            }
        }
        FlagstatFormat::Json => {
            o.push_str("{\n");
            for w in 0..2 {
                let head = if w == 0 {
                    "QC-passed reads"
                } else {
                    "QC-failed reads"
                };
                let _ = writeln!(o, " \"{head}\": {{ ");
                for (i, (label, get, base)) in ROWS.iter().enumerate() {
                    let last = i == ROWS.len() - 1;
                    let a = get(s);
                    // The JSON key says `mapQ >= 5`, the text label `mapQ>=5`.
                    let key = label.replace("mapQ>=5", "mapQ >= 5");
                    let sep = if last { "" } else { "," };
                    let _ = writeln!(o, "  \"{key}\": {}{sep} ", a[w]);
                    if let Some(base) = base {
                        let t = base(s);
                        let _ = writeln!(o, "  \"{key} %\": {}, ", percent(a[w], t[w], true));
                    }
                }
                o.push_str(if w == 0 { " },\n" } else { " }\n" });
            }
            o.push_str("}\n");
        }
    }
    o
}

/// Count every record of `src`, unmapped and unplaced ones included.
pub fn flagstat(mut src: Box<dyn AlnSource>) -> Result<FlagStats, IoError> {
    let header = src.header().clone();
    let name = src.name().to_string();
    let layout = src.layout();
    let mut stats = FlagStats::default();
    let mut buf = Vec::new();
    let mut n = 0u64;
    loop {
        buf.clear();
        if src.next_record(&mut buf)?.is_none() {
            break;
        }
        n += 1;
        let r = aln_ref(layout, &buf, &header).map_err(|msg| IoError::Parse {
            path: name.clone(),
            loc: Loc::Record(n),
            msg,
        })?;
        stats.add(r.flag(), r.ref_id(), r.mate_ref_id(), r.mapq());
    }
    Ok(stats)
}
