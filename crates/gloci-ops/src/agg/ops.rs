//! Resolving specs against a file and accumulating values.

use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::io::Write;

use gloci_core::Interval;
use gloci_io::{Format, Header, IoError, Loc, UsageError, field};

use super::num::{Num, fmt_f64, parse_num};
use super::spec::{AggSpec, Col, Op, OpArg};

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
    arg: OpArg,
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
                arg: s.arg,
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

pub(crate) fn header_names(header: &Header) -> Vec<String> {
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
    /// Sum of weights and of weight * value (`wmean`).
    wsum: f64,
    wxsum: f64,
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
        self.wsum = 0.0;
        self.wxsum = 0.0;
    }

    fn push(&mut self, v: Num, keep_vals: bool, weight: u64) {
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
        // strict comparisons keep the first of `-0` and `0`, as bedtools does
        if f < self.fmin {
            self.fmin = f;
        }
        if f > self.fmax {
            self.fmax = f;
        }
        let d = f - self.mean;
        self.mean += d / self.n as f64;
        self.m2 += d * (f - self.mean);
        self.wsum += weight as f64;
        self.wxsum += weight as f64 * f;
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
    Mode(Tally),
    /// `bases_uniq` and `frac_a`: the shared intervals, and A's length.
    Uniq {
        ivs: Vec<Interval>,
        hits: u64,
        a_len: u64,
    },
    FracB {
        hits: u64,
        shared: u64,
        b_len: u64,
    },
    Records(Vec<Vec<u8>>),
    /// `min_record` / `max_record`: the best value so far and its record.
    Best(Option<(f64, Vec<u8>)>),
}

/// Distinct values with counts, in bytewise order (bedtools' tie rule).
#[derive(Debug, Default)]
struct Tally(BTreeMap<Vec<u8>, u64>);

impl Tally {
    fn add(&mut self, v: &[u8]) {
        match self.0.get_mut(v) {
            Some(c) => *c += 1,
            None => {
                self.0.insert(v.to_vec(), 1);
            }
        }
    }

    fn clear(&mut self) {
        self.0.clear();
    }

    /// The most (or least) frequent value; among equals the bytewise smallest.
    fn pick(&self, most: bool) -> Option<&[u8]> {
        let mut best: Option<(&Vec<u8>, u64)> = None;
        for (v, &c) in &self.0 {
            let better = best.is_none_or(|(_, bc)| if most { c > bc } else { c < bc });
            if better {
                best = Some((v, c));
            }
        }
        best.map(|(v, _)| v.as_slice())
    }
}

/// What the overlap of one accepted hit with the current A record looks
/// like, for the ops that need more than the record's text.
#[derive(Debug, Clone, Copy, Default)]
pub struct HitInfo<'a> {
    /// Bases shared by A and the hit (0 for points).
    pub bases: u64,
    /// A's length (0 for a point).
    pub a_len: u64,
    /// The hit's length (0 for a point).
    pub b_len: u64,
    /// The shared intervals, sorted and without empty ones; their union is
    /// what `bases_uniq` counts.
    pub shared: &'a [Interval],
}

impl State {
    fn new(op: Op) -> State {
        let mut s = match op {
            Op::Count => State::Count(0),
            Op::Bases => State::Bases(0),
            Op::Indicator => State::Indicator(false),
            Op::Mode | Op::AntiMode => State::Mode(Tally::default()),
            Op::BasesUniq | Op::FracA => State::Uniq {
                ivs: Vec::new(),
                hits: 0,
                a_len: 0,
            },
            Op::FracB => State::FracB {
                hits: 0,
                shared: 0,
                b_len: 0,
            },
            Op::Records => State::Records(Vec::new()),
            Op::MinRecord | Op::MaxRecord => State::Best(None),
            Op::Sum
            | Op::Mean
            | Op::Cv
            | Op::Mad
            | Op::Quantile
            | Op::Tmean
            | Op::WMean
            | Op::Median
            | Op::Min
            | Op::Max
            | Op::Stdev
            | Op::Var
            | Op::StdevPop
            | Op::VarPop => State::Numeric(Numeric::default()),
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
            State::Mode(t) => t.clear(),
            State::Uniq { ivs, hits, a_len } => {
                ivs.clear();
                *hits = 0;
                *a_len = 0;
            }
            State::FracB {
                hits,
                shared,
                b_len,
            } => {
                *hits = 0;
                *shared = 0;
                *b_len = 0;
            }
            State::Records(v) => v.clear(),
            State::Best(o) => *o = None,
        }
    }
}

fn not_a_number(agg: &ResolvedAgg, v: &[u8], file: &str, line: u64) -> IoError {
    IoError::Parse {
        path: file.to_string(),
        loc: Loc::Line(line),
        msg: format!(
            "column {} value \"{}\" is not a number ({})",
            agg.col_text,
            String::from_utf8_lossy(v),
            agg.label
        ),
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
    /// absent INFO key) are skipped by every op except `count`, which counts
    /// records. A column past the end of the line is a parse error.
    pub fn add(
        &mut self,
        raw: &[u8],
        overlap_bases: u64,
        file: &str,
        line: u64,
    ) -> Result<(), IoError> {
        let info = HitInfo {
            bases: overlap_bases,
            ..HitInfo::default()
        };
        self.add_hit(raw, &info, file, line)
    }

    /// As [`add`](Self::add), with the full overlap facts the two-input ops
    /// (`wmean`, `bases_uniq`, `frac_*`) need.
    pub fn add_hit(
        &mut self,
        raw: &[u8],
        info: &HitInfo<'_>,
        file: &str,
        line: u64,
    ) -> Result<(), IoError> {
        let overlap_bases = info.bases;
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
                State::Uniq { ivs, hits, a_len } => {
                    *hits += 1;
                    *a_len = info.a_len;
                    ivs.extend_from_slice(info.shared);
                    continue;
                }
                State::FracB {
                    hits,
                    shared,
                    b_len,
                } => {
                    *hits += 1;
                    *shared += overlap_bases;
                    *b_len += info.b_len;
                    continue;
                }
                State::Records(v) => {
                    v.push(
                        raw.iter()
                            .map(|&b| if b == b'\t' { b'|' } else { b })
                            .collect(),
                    );
                    continue;
                }
                _ => {}
            }
            let value = match &agg.src {
                Src::Whole => continue,
                Src::Col(i) => match field(raw, *i) {
                    Some(v) => Some(v),
                    None => {
                        return Err(IoError::Parse {
                            path: file.to_string(),
                            loc: Loc::Line(line),
                            msg: format!(
                                "column {} ({}) is missing: the line has {} tab-separated \
                                 column{}; check the column number and that the file is \
                                 tab-separated (an empty value or `.` is skipped, a missing \
                                 column is an error)",
                                agg.col_text,
                                agg.label,
                                gloci_io::column_count(raw),
                                if gloci_io::column_count(raw) == 1 {
                                    ""
                                } else {
                                    "s"
                                }
                            ),
                        });
                    }
                },
                Src::Info(key) => field(raw, VCF_INFO_COL).and_then(|i| info_value(i, key)),
            };
            let Some(v) = value.filter(|v| !is_missing(v)) else {
                continue;
            };
            match state {
                State::Numeric(n) => {
                    let Some(x) = parse_num(v) else {
                        return Err(not_a_number(agg, v, file, line));
                    };
                    let keep = matches!(agg.op, Op::Median | Op::Mad | Op::Quantile | Op::Tmean);
                    n.push(x, keep, overlap_bases);
                }
                State::Mode(t) => t.add(v),
                State::Best(o) => {
                    let Some(x) = parse_num(v) else {
                        return Err(not_a_number(agg, v, file, line));
                    };
                    let x = x.as_f64();
                    let better = o.as_ref().is_none_or(|(b, _)| {
                        if agg.op == Op::MinRecord {
                            x < *b
                        } else {
                            x > *b
                        }
                    });
                    if better {
                        *o = Some((
                            x,
                            raw.iter()
                                .map(|&b| if b == b'\t' { b'|' } else { b })
                                .collect(),
                        ));
                    }
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
                State::Count(_)
                | State::Bases(_)
                | State::Indicator(_)
                | State::Uniq { .. }
                | State::FracB { .. }
                | State::Records(_) => unreachable!(),
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
                State::Numeric(n) => finish_numeric(agg, n, &self.empty, out),
                State::Mode(t) => match t.pick(agg.op == Op::Mode) {
                    Some(v) => out.extend_from_slice(v),
                    None => out.extend_from_slice(self.empty.as_bytes()),
                },
                State::Uniq { ivs, hits, a_len } => {
                    let uniq = union_len(ivs);
                    if agg.op == Op::BasesUniq {
                        let _ = write!(out, "{uniq}");
                    } else if *hits == 0 || *a_len == 0 {
                        out.extend_from_slice(self.empty.as_bytes());
                    } else {
                        fmt_f64(uniq as f64 / *a_len as f64, out);
                    }
                }
                State::FracB {
                    hits,
                    shared,
                    b_len,
                } => {
                    if *hits == 0 || *b_len == 0 {
                        out.extend_from_slice(self.empty.as_bytes());
                    } else {
                        fmt_f64(*shared as f64 / *b_len as f64, out);
                    }
                }
                State::Records(l) => join(&mut l.iter().map(Vec::as_slice), out, &self.empty),
                State::Best(o) => match o {
                    Some((_, rec)) => out.extend_from_slice(rec),
                    None => out.extend_from_slice(self.empty.as_bytes()),
                },
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

/// Total length of the union of `ivs` (sorted in place).
fn union_len(ivs: &mut [Interval]) -> u64 {
    ivs.sort_unstable();
    let (mut total, mut end) = (0u64, 0u64);
    for iv in ivs.iter() {
        let start = iv.start.max(end);
        if iv.end > start {
            total += iv.end - start;
        }
        end = end.max(iv.end);
    }
    total
}

/// Numeric order where `-0` and `0` tie (values are finite), so a stable sort keeps
/// input order between them as bedtools does.
fn cmp_stable(a: &f64, b: &f64) -> std::cmp::Ordering {
    a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
}

/// The median of sorted, non-empty `v`.
fn median_sorted(v: &[f64]) -> f64 {
    let k = v.len();
    if k % 2 == 1 {
        v[k / 2]
    } else {
        (v[k / 2 - 1] + v[k / 2]) / 2.0
    }
}

/// Linear interpolation between order statistics (type 7) of sorted `v`.
fn quantile_sorted(v: &[f64], q: f64) -> f64 {
    let h = (v.len() - 1) as f64 * q;
    let lo = h.floor() as usize;
    let hi = (lo + 1).min(v.len() - 1);
    if h == lo as f64 || v[lo].to_bits() == v[hi].to_bits() {
        // exact rank, or equal neighbours: no arithmetic, so -0 stays -0 as in `median`
        return v[lo];
    }
    v[lo] + (h - lo as f64) * (v[hi] - v[lo])
}

fn finish_numeric(agg: &ResolvedAgg, n: &mut Numeric, empty: &str, out: &mut Vec<u8>) {
    let op = agg.op;
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
            n.vals.sort_by(cmp_stable);
            fmt_f64(median_sorted(&n.vals), out);
        }
        Op::Mad => {
            n.vals.sort_by(cmp_stable);
            let med = median_sorted(&n.vals);
            let mut dev: Vec<f64> = n.vals.iter().map(|x| (x - med).abs()).collect();
            dev.sort_by(cmp_stable);
            fmt_f64(median_sorted(&dev), out);
        }
        Op::Quantile => {
            let OpArg::Q(q) = agg.arg else { unreachable!() };
            n.vals.sort_by(cmp_stable);
            fmt_f64(quantile_sorted(&n.vals, q), out);
        }
        Op::Tmean => {
            let OpArg::Range(lo, hi) = agg.arg else {
                unreachable!()
            };
            n.vals.sort_by(cmp_stable);
            let last = n.vals.len() - 1;
            // rank fraction of the i-th smallest is i / (n - 1); a lone value is kept
            let kept: Vec<f64> = n
                .vals
                .iter()
                .enumerate()
                .filter(|(i, _)| {
                    last == 0 || {
                        let f = *i as f64 / last as f64;
                        lo <= f && f <= hi
                    }
                })
                .map(|(_, x)| *x)
                .collect();
            if kept.is_empty() {
                return none(out);
            }
            fmt_f64(kept.iter().sum::<f64>() / kept.len() as f64, out);
        }
        Op::WMean => {
            if n.wsum == 0.0 {
                return none(out);
            }
            fmt_f64(n.wxsum / n.wsum, out);
        }
        Op::Cv => {
            if n.n < 2 {
                return none(out);
            }
            let mean = if n.all_int {
                n.isum as f64 / n.n as f64
            } else {
                n.fsum / n.n as f64
            };
            if mean == 0.0 {
                return none(out);
            }
            fmt_f64((n.m2 / (n.n - 1) as f64).sqrt() / mean, out);
        }
        Op::VarPop | Op::StdevPop => {
            let var = n.m2 / n.n as f64;
            fmt_f64(if op == Op::VarPop { var } else { var.sqrt() }, out);
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
