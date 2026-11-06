use gloci_core::thousands;
use std::fmt;

/// Where in the input a parse problem was found.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Loc {
    Line(u64),
    Record(u64),
}

#[derive(Debug, thiserror::Error)]
pub enum IoError {
    #[error("cannot open {path}: {source}; check that the path exists and is readable")]
    Open {
        path: String,
        #[source]
        source: std::io::Error,
    },
    #[error("error reading {path}: {source}; the file may be corrupt or the stream was cut off")]
    Read {
        path: String,
        #[source]
        source: std::io::Error,
    },
    #[error("{path} {loc}: {msg}")]
    Parse { path: String, loc: Loc, msg: String },
    #[error("{msg}")]
    Detect { path: String, msg: String },
    #[error("error writing {path}: {source}; check free disk space and permissions")]
    Write {
        path: String,
        #[source]
        source: std::io::Error,
    },
    #[error("{0}")]
    Usage(#[from] UsageError),
    #[error("{path}: BGZF end-of-file marker missing; the file is probably truncated")]
    Truncated { path: String },
}

impl fmt::Display for Loc {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Loc::Line(n) => write!(f, "line {}", thousands(*n)),
            Loc::Record(n) => write!(f, "record {}", thousands(*n)),
        }
    }
}

/// A command-line usage problem (exit code 2).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UsageError(pub String);

impl fmt::Display for UsageError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for UsageError {}
