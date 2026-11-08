//! Resolving specs against a file and accumulating values.

use std::collections::{BTreeSet, HashSet};
use std::io::Write;

use gloci_io::{Format, Header, IoError, Loc, UsageError, field};

use super::num::{Num, fmt_f64, parse_num};
use super::spec::{AggSpec, Col, Op};

/// Where an aggregation reads its value.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Src {
    /// The whole record (count, bases, indicator only).
    Whole,
    /// A 0-based tab column.
    Col(usize),
    /// `INFO/KEY` of a VCF record.
    Info(Vec<u8>),
}

#[derive(Debug, Clone)]
pub struct ResolvedAgg {
    src: Src,
    op: Op,
    label: String,
    /// How the column is named in a data error.
    col_text: String,
}

const VCF_INFO_COL: usize = 7;

/// Bind specs to the columns of a file of `format` with `header`.
pub fn resolve(
    specs: &[AggSpec],
    format: Format,
    header: &Header,
) -> Result<Vec<ResolvedAgg>, UsageError> {
    let is_vcf = matches!(format, Format::Vcf | Format::Bcf) || matches!(header, Header::Vcf(_));
    specs
        .iter()
        .map(|s| {
            let (src, col_text) = match &s.col {
                Col::Whole => (Src::Whole, "0".to_string()),
                Col::Index(n) => (Src::Col(n - 1), n.to_string()),
                Col::Name(name) if is_vcf => (vcf_name(name, &s.label)?, name.clone()),
                Col::Name(name) => (named_column(name, &s.label, header)?, name.clone()),
            };
            Ok(ResolvedAgg {
                src,
                op: s.op,
                label: s.label.clone(),
                col_text,
            })
        })
        .collect()
}

fn vcf_name(name: &str, label: &str) -> Result<Src, UsageError> {
    match name {
        "ID" => Ok(Src::Col(2)),
        "QUAL" => Ok(Src::Col(5)),
        "FILTER" => Ok(Src::Col(6)),
        _ => match name.strip_prefix("INFO/") {
            Some(key) if !key.is_empty() => Ok(Src::Info(key.as_bytes().to_vec())),
            _ => Err(UsageError(format!(
                "-c {label}: unknown VCF column `{name}`; use a column number, \
                 QUAL, ID, FILTER or INFO/KEY"
            ))),
        },
    }
}

fn header_names(header: &Header) -> Vec<String> {
    let Header::Lines(bytes) = header else {
        return Vec::new();
    };
    let last = bytes
        .split(|&b| b == b'\n')
        .map(|l| l.strip_suffix(b"\r").unwrap_or(l))
        .rfind(|l| !l.is_empty())
        .unwrap_or(b"");
    let last = last.strip_prefix(b"#").unwrap_or(last);
    gloci_io::split_fields(last)
        .map(|f| String::from_utf8_lossy(f).into_owned())
        .collect()
}

fn named_column(name: &str, label: &str, header: &Header) -> Result<Src, UsageError> {
    let names = header_names(header);
    match names.iter().position(|n| n == name) {
        Some(i) => Ok(Src::Col(i)),
        None if names.is_empty() => Err(UsageError(format!(
            "-c {label}: column name `{name}` needs a header line; the input has none \
             (use a column number)"
        ))),
        None => Err(UsageError(format!(
            "-c {label}: no column named `{name}`; available names: {}",
            names.join(", ")
        ))),
    }
}

fn is_missing(v: &[u8]) -> bool {
    v.is_empty() || v == b"."
}

/// The value of `INFO` key `key`: a flag gives `1`.
fn info_value<'a>(info: &'a [u8], key: &[u8]) -> Option<&'a [u8]> {
    for kv in info.split(|&b| b == b';') {
        match kv.iter().position(|&b| b == b'=') {
            Some(i) if &kv[..i] == key => return Some(&kv[i + 1..]),
            None if kv == key => return Some(b"1"),
            _ => {}
        }
    }
    None
}

#[derive(Debug, Default)]
struct Numeric {
    n: u64,
    all_int: bool,
    isum: i128,
    fsum: f64,
    imin: i128,
    imax: i128,
    fmin: f64,
    fmax: f64,
    mean: f64,
    m2: f64,
    vals: Vec<f64>,
}

impl Numeric {
    fn reset(&mut self) {
        self.n = 0;
        self.all_int = true;
        self.isum = 0;
        self.fsum = 0.0;
        self.imin = i128::MAX;
        self.imax = i128::MIN;
        self.fmin = f64::INFINITY;
        self.fmax = f64::NEG_INFINITY;
        self.mean = 0.0;
        self.m2 = 0.0;
        self.vals.clear();
    }

    fn push(&mut self, v: Num, keep_vals: bool) {
        self.n += 1;
        let f = v.as_f64();
        match v {
            Num::Int(i) => {
                match self.isum.checked_add(i) {
                    Some(s) => self.isum = s,
                    None => self.all_int = false,
                }
                self.imin = self.imin.min(i);
                self.imax = self.imax.max(i);
            }
            Num::Float(_) => self.all_int = false,
        }
        self.fsum += f;
        self.fmin = self.fmin.min(f);
        self.fmax = self.fmax.max(f);
        let d = f - self.mean;
        self.mean += d / self.n as f64;
        self.m2 += d * (f - self.mean);
        if keep_vals {
            self.vals.push(f);
        }
    }
}

#[derive(Debug)]
enum State {
    Count(u64),
    Bases(u64),
    Indicator(bool),
    Numeric(Numeric),
    CountDistinct(HashSet<Vec<u8>>),
    List(Vec<Vec<u8>>),
    Distinct(HashSet<Vec<u8>>, Vec<Vec<u8>>),
    Sorted(BTreeSet<Vec<u8>>),
    First(Option<Vec<u8>>),
    Last(Option<Vec<u8>>),
}

impl State {
    fn new(op: Op) -> State {
        let mut s = match op {
            Op::Count => State::Count(0),
            Op::Bases => State::Bases(0),
            Op::Indicator => State::Indicator(false),
            Op::Sum | Op::Mean | Op::Median | Op::Min | Op::Max | Op::Stdev | Op::Var => {
                State::Numeric(Numeric::default())
            }
            Op::CountDistinct => State::CountDistinct(HashSet::new()),
            Op::Collapse => State::List(Vec::new()),
            Op::Distinct => State::Distinct(HashSet::new(), Vec::new()),
            Op::DistinctSorted => State::Sorted(BTreeSet::new()),
            Op::First => State::First(None),
            Op::Last => State::Last(None),
        };
        s.reset();
        s
    }

    fn reset(&mut self) {
        match self {
            State::Count(n) | State::Bases(n) => *n = 0,
            State::Indicator(b) => *b = false,
            State::Numeric(n) => n.reset(),
            State::CountDistinct(s) => s.clear(),
            State::List(v) => v.clear(),
            State::Distinct(s, v) => {
                s.clear();
                v.clear();
            }
            State::Sorted(s) => s.clear(),
            State::First(o) | State::Last(o) => *o = None,
        }
    }
}

/// One accumulator per aggregation, reused across A records.
pub struct Accumulators {
    aggs: Vec<ResolvedAgg>,
    states: Vec<State>,
    empty: String,
    delim: String,
}

impl Accumulators {
    pub fn new(aggs: Vec<ResolvedAgg>, empty: String, delim: String) -> Self {
        let states = aggs.iter().map(|a| State::new(a.op)).collect();
        Accumulators {
            aggs,
            states,
            empty,
            delim,
        }
    }

    /// Fold one overlapping record (`raw` is its text line, `overlap_bases`
    /// its overlap with the current A record). Missing values (`.`, empty, an
    /// absent column or INFO key) are skipped by every op except `count`,
    /// which counts records.
    pub fn add(
        &mut self,
        raw: &[u8],
        overlap_bases: u64,
        file: &str,
        line: u64,
    ) -> Result<(), IoError> {
        for (agg, state) in self.aggs.iter().zip(self.states.iter_mut()) {
            match state {
                State::Count(n) => {
                    *n += 1;
                    continue;
                }
                State::Bases(n) => {
                    *n += overlap_bases;
                    continue;
                }
                State::Indicator(b) => {
                    *b = true;
                    continue;
                }
                _ => {}
            }
            let value = match &agg.src {
                Src::Whole => continue,
                Src::Col(i) => field(raw, *i),
                Src::Info(key) => field(raw, VCF_INFO_COL).and_then(|i| info_value(i, key)),
            };
            let Some(v) = value.filter(|v| !is_missing(v)) else {
                continue;
            };
            match state {
                State::Numeric(n) => {
                    let Some(x) = parse_num(v) else {
                        return Err(IoError::Parse {
                            path: file.to_string(),
                            loc: Loc::Line(line),
                            msg: format!(
                                "column {} value \"{}\" is not a number ({})",
                                agg.col_text,
                                String::from_utf8_lossy(v),
                                agg.label
                            ),
                        });
                    };
                    n.push(x, agg.op == Op::Median);
                }
                State::CountDistinct(s) => {
                    if !s.contains(v) {
                        s.insert(v.to_vec());
                    }
                }
                State::List(l) => l.push(v.to_vec()),
                State::Distinct(s, l) => {
                    if !s.contains(v) {
                        s.insert(v.to_vec());
                        l.push(v.to_vec());
                    }
                }
                State::Sorted(s) => {
                    if !s.contains(v) {
                        s.insert(v.to_vec());
                    }
                }
                State::First(o) => {
                    if o.is_none() {
                        *o = Some(v.to_vec());
                    }
                }
                State::Last(o) => *o = Some(v.to_vec()),
                State::Count(_) | State::Bases(_) | State::Indicator(_) => unreachable!(),
            }
        }
        Ok(())
    }

    /// Append `\tv1\tv2…` for the current A record and reset.
    pub fn finish_into(&mut self, out: &mut Vec<u8>) {
        for (agg, state) in self.aggs.iter().zip(self.states.iter_mut()) {
            out.push(b'\t');
            let delim = self.delim.as_bytes();
            let join = |items: &mut dyn Iterator<Item = &[u8]>, out: &mut Vec<u8>, empty: &str| {
                let mut any = false;
                for it in items {
                    if any {
                        out.extend_from_slice(delim);
                    }
                    any = true;
                    out.extend_from_slice(it);
                }
                if !any {
                    out.extend_from_slice(empty.as_bytes());
                }
            };
            match state {
                State::Count(n) | State::Bases(n) => {
                    let _ = write!(out, "{n}");
                }
                State::Indicator(b) => out.push(if *b { b'1' } else { b'0' }),
                State::CountDistinct(s) => {
                    let _ = write!(out, "{}", s.len());
                }
                State::Numeric(n) => finish_numeric(agg.op, n, &self.empty, out),
                State::List(l) => join(&mut l.iter().map(Vec::as_slice), out, &self.empty),
                State::Distinct(_, l) => join(&mut l.iter().map(Vec::as_slice), out, &self.empty),
                State::Sorted(s) => {
                    join(&mut s.iter().map(Vec::as_slice), out, &self.empty);
                }
                State::First(o) | State::Last(o) => match o {
                    Some(v) => out.extend_from_slice(v),
                    None => out.extend_from_slice(self.empty.as_bytes()),
                },
            }
            state.reset();
        }
    }
}

fn finish_numeric(op: Op, n: &mut Numeric, empty: &str, out: &mut Vec<u8>) {
    let none = |out: &mut Vec<u8>| out.extend_from_slice(empty.as_bytes());
    if n.n == 0 {
        return none(out);
    }
    match op {
        Op::Sum if n.all_int => {
            let _ = write!(out, "{}", n.isum);
        }
        Op::Sum => fmt_f64(n.fsum, out),
        Op::Min if n.all_int => {
            let _ = write!(out, "{}", n.imin);
        }
        Op::Min => fmt_f64(n.fmin, out),
        Op::Max if n.all_int => {
            let _ = write!(out, "{}", n.imax);
        }
        Op::Max => fmt_f64(n.fmax, out),
        Op::Mean => {
            let m = if n.all_int {
                n.isum as f64 / n.n as f64
            } else {
                n.fsum / n.n as f64
            };
            fmt_f64(m, out);
        }
        Op::Median => {
            n.vals.sort_unstable_by(|a, b| a.total_cmp(b));
            let k = n.vals.len();
            let m = if k % 2 == 1 {
                n.vals[k / 2]
            } else {
                (n.vals[k / 2 - 1] + n.vals[k / 2]) / 2.0
            };
            fmt_f64(m, out);
        }
        Op::Var | Op::Stdev => {
            if n.n < 2 {
                return none(out);
            }
            let var = n.m2 / (n.n - 1) as f64;
            fmt_f64(if op == Op::Var { var } else { var.sqrt() }, out);
        }
        _ => unreachable!("not a numeric op"),
    }
}
