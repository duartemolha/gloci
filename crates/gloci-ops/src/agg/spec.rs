//! Parsing `-c COL:OP` arguments.

use gloci_io::UsageError;

/// The P1 aggregation operations (D46).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    Count,
    CountDistinct,
    Sum,
    Mean,
    Median,
    Min,
    Max,
    Stdev,
    Var,
    StdevPop,
    VarPop,
    Collapse,
    Distinct,
    DistinctSorted,
    First,
    Last,
    Bases,
    Indicator,
    Cv,
    Mad,
    /// `quantile:Q` (argument in [`OpArg::Q`]).
    Quantile,
    /// `tmean:LO-HI` (argument in [`OpArg::Range`]).
    Tmean,
    WMean,
    Mode,
    AntiMode,
    BasesUniq,
    FracA,
    FracB,
    MinRecord,
    MaxRecord,
    Records,
}

impl Op {
    /// Ops that read the overlap with a B record, so they have no meaning
    /// where there is no second input (`merge`, `groupby`).
    pub fn needs_second_input(self) -> bool {
        matches!(
            self,
            Op::Bases
                | Op::WMean
                | Op::BasesUniq
                | Op::FracA
                | Op::FracB
                | Op::MinRecord
                | Op::MaxRecord
                | Op::Records
        )
    }

    /// Ops that need the per-hit overlap facts (`HitInfo`).
    pub fn needs_overlap(self) -> bool {
        matches!(
            self,
            Op::Bases | Op::WMean | Op::BasesUniq | Op::FracA | Op::FracB
        )
    }
}

/// The argument of `quantile:Q` and `tmean:LO-HI`.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum OpArg {
    #[default]
    None,
    Q(f64),
    Range(f64, f64),
}

/// A column selector: `0` is the whole record, a number is 1-based, anything
/// else is a header or VCF field name.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Col {
    Whole,
    Index(usize),
    Name(String),
}

#[derive(Debug, Clone, PartialEq)]
pub struct AggSpec {
    pub col: Col,
    pub op: Op,
    pub arg: OpArg,
    /// The text as typed, e.g. `5:mean`.
    pub label: String,
}

fn parse_op(name: &str) -> Option<Op> {
    Some(match name {
        "count" => Op::Count,
        "count_distinct" => Op::CountDistinct,
        "sum" => Op::Sum,
        "mean" => Op::Mean,
        "median" => Op::Median,
        "min" => Op::Min,
        "max" => Op::Max,
        "stdev" => Op::Stdev,
        "var" => Op::Var,
        "stdev_pop" => Op::StdevPop,
        "var_pop" => Op::VarPop,
        "collapse" => Op::Collapse,
        "distinct" => Op::Distinct,
        "distinct_sorted" => Op::DistinctSorted,
        "first" => Op::First,
        "last" => Op::Last,
        "bases" => Op::Bases,
        "indicator" => Op::Indicator,
        "cv" => Op::Cv,
        "mad" => Op::Mad,
        "quantile" => Op::Quantile,
        "tmean" => Op::Tmean,
        "wmean" => Op::WMean,
        "mode" => Op::Mode,
        "antimode" => Op::AntiMode,
        "bases_uniq" => Op::BasesUniq,
        "frac_a" => Op::FracA,
        "frac_b" => Op::FracB,
        "min_record" => Op::MinRecord,
        "max_record" => Op::MaxRecord,
        "records" => Op::Records,
        _ => return None,
    })
}

/// Parse repeatable, comma-separated `COL:OP` arguments.
pub fn parse_aggs(args: &[String]) -> Result<Vec<AggSpec>, UsageError> {
    let mut out = Vec::new();
    for arg in args {
        for item in arg.split(',').map(str::trim).filter(|s| !s.is_empty()) {
            out.push(parse_one(item)?);
        }
    }
    Ok(out)
}

fn parse_one(item: &str) -> Result<AggSpec, UsageError> {
    let Some((col_text, op_text)) = item.split_once(':') else {
        return Err(UsageError(format!(
            "-c {item}: expected COL:OP (for example 5:mean); \
             ops: count count_distinct sum mean median min max stdev var stdev_pop var_pop collapse \
             distinct distinct_sorted first last bases indicator cv mad quantile:Q tmean:LO-HI wmean mode \
             antimode bases_uniq frac_a frac_b min_record max_record records"
        )));
    };
    // `quantile:0.5` style ops carry arguments after their own colon.
    let (op_name, arg_text) = match op_text.split_once(':') {
        Some((n, a)) => (n, Some(a)),
        None => (op_text, None),
    };
    let Some(op) = parse_op(op_name) else {
        return Err(UsageError(format!(
            "-c {item}: unknown op `{op_name}`; ops: count count_distinct sum mean \
             median min max stdev var stdev_pop var_pop collapse distinct distinct_sorted first last \
             bases indicator cv mad quantile:Q tmean:LO-HI wmean mode antimode bases_uniq frac_a \
             frac_b min_record max_record records"
        )));
    };
    let arg = parse_arg(item, op, op_name, arg_text)?;
    let col = match col_text {
        "" => return Err(UsageError(format!("-c {item}: missing column before `:`"))),
        "0" => Col::Whole,
        t if t.bytes().all(|b| b.is_ascii_digit()) => match t.parse::<usize>() {
            Ok(n) => Col::Index(n),
            Err(_) => return Err(UsageError(format!("-c {item}: column `{t}` is too large"))),
        },
        t => Col::Name(t.to_string()),
    };
    if col == Col::Whole
        && !matches!(
            op,
            Op::Count
                | Op::Bases
                | Op::Indicator
                | Op::BasesUniq
                | Op::FracA
                | Op::FracB
                | Op::Records
        )
    {
        return Err(UsageError(format!(
            "-c {item}: column 0 (the whole record) works only with count, bases, indicator, \
             bases_uniq, frac_a, frac_b and records"
        )));
    }
    Ok(AggSpec {
        col,
        op,
        arg,
        label: item.to_string(),
    })
}

fn unit(item: &str, what: &str, text: &str) -> Result<f64, UsageError> {
    match text.parse::<f64>() {
        Ok(v) if (0.0..=1.0).contains(&v) => Ok(v),
        _ => Err(UsageError(format!(
            "-c {item}: {what} `{text}` must be a number from 0 to 1"
        ))),
    }
}

fn parse_arg(item: &str, op: Op, op_name: &str, arg: Option<&str>) -> Result<OpArg, UsageError> {
    match (op, arg) {
        (Op::Quantile, Some(q)) => Ok(OpArg::Q(unit(item, "quantile", q)?)),
        (Op::Quantile, None) => Err(UsageError(format!(
            "-c {item}: quantile needs a quantile, for example 5:quantile:0.9"
        ))),
        (Op::Tmean, Some(r)) => {
            let hint = || {
                UsageError(format!(
                    "-c {item}: tmean needs a range LO-HI, for example 5:tmean:0.1-0.9"
                ))
            };
            let Some((lo, hi)) = r
                .split_once('-')
                .filter(|(l, h)| !l.is_empty() && !h.is_empty())
            else {
                return Err(hint());
            };
            let (lo, hi) = (
                unit(item, "tmean bound", lo)?,
                unit(item, "tmean bound", hi)?,
            );
            if lo >= hi {
                return Err(UsageError(format!(
                    "-c {item}: tmean range must have LO below HI; write for example \
                     5:tmean:0.1-0.9"
                )));
            }
            Ok(OpArg::Range(lo, hi))
        }
        (Op::Tmean, None) => Err(UsageError(format!(
            "-c {item}: tmean needs a range LO-HI, for example 5:tmean:0.1-0.9"
        ))),
        (_, None) => Ok(OpArg::None),
        (_, Some(_)) => Err(UsageError(format!(
            "-c {item}: {op_name} takes no argument"
        ))),
    }
}

/// Refuse the ops that read an overlapping B record in a command that has no
/// second input. `command` is the command's name for the message (`merge`,
/// `groupby`).
pub fn refuse_two_input_ops(aggs: &[AggSpec], command: &str) -> Result<(), UsageError> {
    match aggs.iter().find(|a| a.op.needs_second_input()) {
        Some(a) => Err(UsageError(format!(
            "-c {}: this op measures or copies the overlapping record of another input, and \
             {command} has none; use map for it, or count / collapse / first / last here",
            a.label
        ))),
        None => Ok(()),
    }
}
