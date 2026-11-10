//! `--summary` / `--json`: a short overview of the selected records.

use super::engine::{Out, Stat, StatSink};
use super::fmt::{write_g6, write_u64};
use crate::filter::SiteView;
use crate::gt::GenotypeSource;
use crate::stats::SiteStats;
use gloci_core::{VariantKind, alt_kind};
use std::io;

/// The version of the JSON layout; bumped when a key is renamed or removed.
pub const SUMMARY_SCHEMA_VERSION: u64 = 1;

/// Record types in output order.
const TYPES: [&str; 5] = ["snp", "mnp", "indel", "other", "ref"];

fn type_index(kind: VariantKind) -> usize {
    match kind {
        VariantKind::Snv => 0,
        VariantKind::Mnv => 1,
        VariantKind::Ins | VariantKind::Del => 2,
        VariantKind::Sv | VariantKind::Other => 3,
        VariantKind::Ref => 4,
    }
}

/// The type of a record: that of its ALT alleles when they agree (a ref-only
/// ALT such as `<*>` is ignored beside a real one), `other` when they mix,
/// `ref` when there is no ALT.
fn record_type(pos: u64, ref_: &[u8], alts: &[u8]) -> usize {
    let mut seen = [false; 5];
    for alt in alts.split(|&c| c == b',') {
        seen[type_index(alt_kind(pos.saturating_sub(1), ref_, alt))] = true;
    }
    let real: Vec<usize> = (0..4).filter(|&i| seen[i]).collect();
    match real[..] {
        [] => 4,
        [one] => one,
        _ => 3,
    }
}

/// Transition or transversion of a biallelic single-base SNP.
fn ts_or_tv(ref_: &[u8], alts: &[u8]) -> Option<bool> {
    let (&[r], &[a]) = (ref_, alts) else {
        return None;
    };
    let (r, a) = (r.to_ascii_uppercase(), a.to_ascii_uppercase());
    if r == a || !b"ACGT".contains(&r) || !b"ACGT".contains(&a) {
        return None;
    }
    let purine = |b: u8| b == b'A' || b == b'G';
    Some(purine(r) == purine(a))
}

/// The counters behind both summary files.
#[derive(Default)]
struct Totals {
    records: u64,
    multiallelic: u64,
    types: [u64; 5],
    ts: u64,
    tv: u64,
    genotypes: u64,
    missing: u64,
    filtered: u64,
    dp_sum: u64,
    dp_n: u64,
}

fn ratio(num: u64, den: u64) -> Option<f64> {
    (den > 0).then(|| num as f64 / den as f64)
}

fn num_or(b: &mut Vec<u8>, v: Option<f64>, absent: &[u8]) {
    match v {
        Some(v) => write_g6(b, v),
        None => b.extend_from_slice(absent),
    }
}

fn kv_u64(b: &mut Vec<u8>, key: &str, v: u64) {
    b.extend_from_slice(key.as_bytes());
    write_u64(b, v);
}

impl Totals {
    fn text(&self, samples: u64, b: &mut Vec<u8>) {
        kv_u64(b, "Records:        ", self.records);
        b.push(b'\n');
        kv_u64(b, "Samples:        ", samples);
        b.push(b'\n');
        kv_u64(b, "Multiallelic:   ", self.multiallelic);
        b.extend_from_slice(b" records\nCounts by type:\n");
        for (name, n) in TYPES.iter().zip(self.types) {
            b.extend_from_slice(b"  ");
            b.extend_from_slice(name.as_bytes());
            for _ in name.len()..6 {
                b.push(b' ');
            }
            write_u64(b, n);
            b.push(b'\n');
        }
        b.extend_from_slice(b"Ts/Tv:          ");
        num_or(b, ratio(self.ts, self.tv), b"n/a");
        b.extend_from_slice(b" (");
        write_u64(b, self.ts);
        b.extend_from_slice(b" transitions, ");
        write_u64(b, self.tv);
        b.extend_from_slice(b" transversions, biallelic SNPs)\nMean depth:     ");
        num_or(b, ratio(self.dp_sum, self.dp_n), b"n/a");
        b.extend_from_slice(b" (over ");
        write_u64(b, self.dp_n);
        b.extend_from_slice(b" genotypes with a depth)\nMissingness:    ");
        num_or(b, ratio(self.missing, self.genotypes), b"n/a");
        b.extend_from_slice(b" (");
        write_u64(b, self.missing);
        b.extend_from_slice(b" of ");
        write_u64(b, self.genotypes);
        b.extend_from_slice(b" genotypes missing, ");
        write_u64(b, self.filtered);
        b.extend_from_slice(b" of them masked)\n");
    }

    /// Keys are sorted, at every level, so the file is stable.
    fn json(&self, samples: u64, b: &mut Vec<u8>) {
        b.extend_from_slice(b"{\n  \"counts_by_type\": {\n");
        // Alphabetical: indel, mnp, other, ref, snp.
        let order = [2usize, 1, 3, 4, 0];
        for (j, &i) in order.iter().enumerate() {
            b.extend_from_slice(b"    \"");
            b.extend_from_slice(TYPES[i].as_bytes());
            b.extend_from_slice(b"\": ");
            write_u64(b, self.types[i]);
            b.extend_from_slice(if j + 1 < order.len() { b",\n" } else { b"\n" });
        }
        b.extend_from_slice(b"  },\n  \"depth\": {\n    \"mean\": ");
        num_or(b, ratio(self.dp_sum, self.dp_n), b"null");
        b.extend_from_slice(b",\n");
        kv_u64(b, "    \"values\": ", self.dp_n);
        b.extend_from_slice(b"\n  },\n  \"genotypes\": {\n");
        kv_u64(b, "    \"filtered\": ", self.filtered);
        b.extend_from_slice(b",\n    \"missing\": ");
        write_u64(b, self.missing);
        b.extend_from_slice(b",\n    \"missing_fraction\": ");
        num_or(b, ratio(self.missing, self.genotypes), b"null");
        b.extend_from_slice(b",\n");
        kv_u64(b, "    \"total\": ", self.genotypes);
        b.extend_from_slice(b"\n  },\n");
        kv_u64(b, "  \"multiallelic_sites\": ", self.multiallelic);
        b.extend_from_slice(b",\n");
        kv_u64(b, "  \"records\": ", self.records);
        b.extend_from_slice(b",\n");
        kv_u64(b, "  \"samples\": ", samples);
        b.extend_from_slice(b",\n");
        kv_u64(b, "  \"schema_version\": ", SUMMARY_SCHEMA_VERSION);
        b.extend_from_slice(b",\n  \"ts_tv\": {\n    \"ratio\": ");
        num_or(b, ratio(self.ts, self.tv), b"null");
        b.extend_from_slice(b",\n");
        kv_u64(b, "    \"transitions\": ", self.ts);
        b.extend_from_slice(b",\n");
        kv_u64(b, "    \"transversions\": ", self.tv);
        b.extend_from_slice(b"\n  }\n}\n");
    }
}

/// `--summary` (`summary.txt`) and `--json` (`summary.json`): the record
/// count and types, Ts/Tv, mean depth and missingness of the selected
/// records and samples. Masked genotypes count as missing (D50) and are also
/// reported as masked. It never fails on polyploid genotypes.
pub struct Summary {
    samples: u64,
    t: Totals,
    txt: Option<Out>,
    json: Option<Out>,
}

impl Summary {
    /// `samples`: the number of selected samples. `txt` and `json` choose
    /// the files written; at least one must be set.
    pub fn new(samples: u64, txt: bool, json: bool) -> Summary {
        Summary {
            samples,
            t: Totals::default(),
            txt: txt.then(Out::closed),
            json: json.then(Out::closed),
        }
    }
}

impl Stat for Summary {
    fn start(&mut self, sink: &mut dyn StatSink) -> io::Result<()> {
        if let Some(o) = &mut self.txt {
            o.open(sink, "summary.txt")?;
        }
        if let Some(o) = &mut self.json {
            o.open(sink, "summary.json")?;
        }
        Ok(())
    }
    fn site(
        &mut self,
        s: &SiteView,
        _gt: Option<&dyn GenotypeSource>,
        st: Option<&SiteStats>,
    ) -> Result<(), String> {
        let st = st.ok_or("internal error: the summary needs the site statistics")?;
        let t = &mut self.t;
        t.records += 1;
        if s.alts.contains(&b',') {
            t.multiallelic += 1;
        }
        t.types[record_type(s.pos, s.ref_, s.alts)] += 1;
        match ts_or_tv(s.ref_, s.alts) {
            Some(true) => t.ts += 1,
            Some(false) => t.tv += 1,
            None => {}
        }
        t.genotypes += st.samples;
        t.missing += st.missing;
        t.filtered += st.filtered;
        t.dp_sum += st.dp_sum;
        t.dp_n += st.dp_n;
        Ok(())
    }
    fn finish(&mut self) -> io::Result<()> {
        if let Some(o) = &mut self.txt {
            self.t.text(self.samples, &mut o.buf);
            o.row().map_err(io::Error::other)?;
            o.commit()?;
        }
        if let Some(o) = &mut self.json {
            self.t.json(self.samples, &mut o.buf);
            o.row().map_err(io::Error::other)?;
            o.commit()?;
        }
        Ok(())
    }
    fn needs_genotypes(&self) -> bool {
        true
    }
    fn ext(&self) -> &'static str {
        if self.txt.is_some() {
            "summary.txt"
        } else {
            "summary.json"
        }
    }
}
