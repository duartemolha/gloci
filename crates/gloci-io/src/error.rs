use gloci_core::thousands;
use std::fmt;

/// Where in the input a parse problem was found.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Loc {
    Line(u64),
    Record(u64),
}

/// Start of the name of the reader over `--sort`ed inputs.
pub const SORTED_INPUTS: &str = "the sorted inputs (";

#[derive(Debug, thiserror::Error)]
pub enum IoError {
    #[error("cannot open {path}: {source}; check that the path exists and is readable")]
    Open {
        path: String,
        #[source]
        source: std::io::Error,
    },
    #[error("cannot create {path}: {source}; check that the directory exists and is writable")]
    Create {
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
    #[error("{0}")]
    Order(Box<OrderError>),
    #[error("internal error: {0}; this is a bug, please report it")]
    Internal(String),
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

/// A record that breaks the input's sort order; both records are quoted.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OrderError {
    pub path: String,
    pub prev_loc: Loc,
    pub prev: String,
    pub loc: Loc,
    pub rec: String,
    pub hint: String,
}

impl fmt::Display for OrderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} is not sorted: {} ({}) comes after {} ({})\nhint: {}",
            self.path, self.loc, self.rec, self.prev_loc, self.prev, self.hint
        )
    }
}

impl std::error::Error for OrderError {}

// Boxed so that `Result<_, IoError>` stays small (clippy::result_large_err).
impl From<OrderError> for IoError {
    fn from(e: OrderError) -> IoError {
        IoError::Order(Box::new(e))
    }
}

impl IoError {
    /// Libraries report output write errors under the placeholder path
    /// `output`; name the real destination, keeping the `io::ErrorKind` so a
    /// broken pipe is still recognised.
    pub fn with_output_name(self, name: &str) -> IoError {
        match self {
            IoError::Write { path, source } if path == "output" => IoError::Write {
                path: name.to_string(),
                source,
            },
            e => e,
        }
    }

    /// A data error found in the stream `--sort` produced names the sorted
    /// inputs, whose line numbers mean nothing to the user: count records.
    pub fn with_sorted_input_loc(self) -> IoError {
        match self {
            IoError::Parse {
                path,
                loc: Loc::Line(n),
                msg,
            } if path.starts_with(SORTED_INPUTS) => IoError::Parse {
                path,
                loc: Loc::Record(n),
                msg,
            },
            e => e,
        }
    }

    /// An order error from a command that can sort its inputs itself gets
    /// `extra` appended to its hint (for example the `--sort` flag).
    pub fn with_order_hint(self, extra: &str) -> IoError {
        match self {
            IoError::Order(mut e) => {
                e.hint.push_str(extra);
                IoError::Order(e)
            }
            e => e,
        }
    }

    /// An order error found while indexing the output of `input` names the
    /// file that is out of order. The output is not published yet, so it
    /// is not the file the user can fix.
    pub fn with_index_input(self, input: &str, output: &str) -> IoError {
        match self {
            IoError::Order(mut e) if e.path == output => {
                e.hint = e.hint.replace(output, input);
                e.path = input.to_string();
                IoError::Order(e)
            }
            e => e,
        }
    }

    /// Like `with_output_name`, for read errors reported under the
    /// placeholder path `input`.
    pub fn with_input_name(self, name: &str) -> IoError {
        match self {
            IoError::Read { path, source } if path == "input" => IoError::Read {
                path: name.to_string(),
                source,
            },
            e => e,
        }
    }
}
