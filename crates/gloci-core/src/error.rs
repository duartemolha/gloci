//! Error types for the core crate.

use thiserror::Error;

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum CoordError {
    #[error("\"{0}\" is not an integer")]
    NotInteger(String),
    #[error("\"{0}\" is negative")]
    Negative(String),
    #[error("\"{0}\" is too large")]
    TooLarge(String),
    #[error("start {start} is greater than end {end}")]
    StartAfterEnd { start: u64, end: u64 },
    #[error("start must be >= 1 (1-based coordinates)")]
    StartBelowOne,
}

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum RegionError {
    #[error("empty region")]
    Empty,
    #[error("invalid region \"{region}\": {source}")]
    Coord {
        region: String,
        #[source]
        source: CoordError,
    },
    #[error("invalid region \"{0}\": missing chromosome name")]
    MissingChrom(String),
    #[error("invalid region \"{0}\": unterminated '{{' (expected {{name}}:start-end)")]
    Unterminated(String),
    #[error(
        "region \"{0}\" is ambiguous: both the whole string and its prefix are chromosome names; write it as {{name}}:start-end"
    )]
    Ambiguous(String),
}

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum GenomeError {
    #[error("duplicate chromosome \"{0}\" in genome file")]
    Duplicate(String),
    #[error("genome line {line}: {msg}")]
    Line { line: usize, msg: String },
}

#[derive(Debug, Clone, PartialEq, Eq, Error)]
#[error("chromosome \"{0}\" is not in the genome")]
pub struct UnknownChrom(pub String);

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum FractionError {
    #[error("\"{0}\" is not a plain decimal fraction (use e.g. 0.5 or 1; no signs or exponents)")]
    Syntax(String),
    #[error("fraction \"{0}\" must be greater than 0 and at most 1")]
    OutOfRange(String),
    #[error("fraction \"{0}\" has more than 18 fractional digits")]
    TooPrecise(String),
}
