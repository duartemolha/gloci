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
    Collapse,
    Distinct,
    DistinctSorted,
    First,
    Last,
    Bases,
    Indicator,
}

/// A column selector: `0` is the whole record, a number is 1-based, anything
/// else is a header or VCF field name.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Col {
    Whole,
    Index(usize),
    Name(String),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AggSpec {
    pub col: Col,
    pub op: Op,
    /// The text as typed, e.g. `5:mean`.
    pub label: String,
}

/// Spec §6.3 ops that parse but are not in P1.
const NOT_YET: &[&str] = &[
    "stdev_pop",
    "var_pop",
    "cv",
    "mad",
    "quantile",
    "tmean",
    "wmean",
    "mode",
    "antimode",
    "bases_uniq",
    "frac_a",
    "frac_b",
    "min_record",
    "max_record",
    "records",
];

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
        "collapse" => Op::Collapse,
        "distinct" => Op::Distinct,
        "distinct_sorted" => Op::DistinctSorted,
        "first" => Op::First,
        "last" => Op::Last,
        "bases" => Op::Bases,
        "indicator" => Op::Indicator,
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
             ops: count count_distinct sum mean median min max stdev var collapse \
             distinct distinct_sorted first last bases indicator"
        )));
    };
    // `quantile:0.5` style ops carry arguments after their own colon.
    let op_name = op_text.split(':').next().unwrap_or(op_text);
    let Some(op) = parse_op(op_name) else {
        if NOT_YET.contains(&op_name) {
            return Err(UsageError(format!(
                "-c {item}: {op_name} is not available yet"
            )));
        }
        return Err(UsageError(format!(
            "-c {item}: unknown op `{op_name}`; ops: count count_distinct sum mean \
             median min max stdev var collapse distinct distinct_sorted first last \
             bases indicator"
        )));
    };
    if op_name != op_text {
        return Err(UsageError(format!(
            "-c {item}: {op_name} takes no argument"
        )));
    }
    let col = match col_text {
        "" => return Err(UsageError(format!("-c {item}: missing column before `:`"))),
        "0" => Col::Whole,
        t if t.bytes().all(|b| b.is_ascii_digit()) => match t.parse::<usize>() {
            Ok(n) => Col::Index(n),
            Err(_) => return Err(UsageError(format!("-c {item}: column `{t}` is too large"))),
        },
        t => Col::Name(t.to_string()),
    };
    if col == Col::Whole && !matches!(op, Op::Count | Op::Bases | Op::Indicator) {
        return Err(UsageError(format!(
            "-c {item}: column 0 (the whole record) works only with count, bases and indicator"
        )));
    }
    Ok(AggSpec {
        col,
        op,
        label: item.to_string(),
    })
}
