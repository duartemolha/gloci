//! A forward stream over the lines of a FASTA file (plain, gzip or BGZF).

use crate::detect::Format;
use crate::error::IoError;
use crate::input::{InputOptions, open_input};
use memchr::memchr;
use std::io::Read;
use std::num::NonZeroUsize;

/// The most bytes of one line the stream holds; longer sequence lines come in pieces.
const CAP: usize = 1 << 20;

/// One item of a [`FastaStream`]; line breaks (LF or CRLF) are removed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FastaLine<'a> {
    /// A header line without its `>`.
    Header(&'a [u8]),
    /// A sequence line; a blank line is an empty `Seq`. A line longer than
    /// 1 MiB is passed in several consecutive pieces.
    Seq(&'a [u8]),
}

/// Lines of a FASTA, one at a time, in bounded memory (about 1 MiB).
///
/// A line longer than 1 MiB comes in pieces, and the stream may emit an empty
/// `Seq` piece when a line break falls exactly on the 1 MiB boundary.
/// Consumers must tolerate an empty `Seq`.
pub struct FastaStream {
    name: String,
    src: Box<dyn Read + Send>,
    buf: Vec<u8>,
    start: usize,
    end: usize,
    /// Where the item last returned ends in `buf`.
    consumed_to: usize,
    eof: bool,
    /// The last item ended without its line break (a piece of a long line).
    mid_line: bool,
    line: u64,
    crlf: bool,
    terminated: bool,
}

impl FastaStream {
    /// Open `path` (`-` is standard input). The format must be FASTA.
    pub fn open(path: &str, threads: NonZeroUsize) -> Result<Self, IoError> {
        let mut input = open_input(
            path,
            &InputOptions {
                threads,
                format: None,
            },
        )?;
        if input.format != Format::Fasta {
            // An empty file is detected as BED (an empty interval file is
            // valid), but an empty FASTA is a mistake worth naming.
            let mut first = [0u8; 1];
            let empty = loop {
                match input.reader.read(&mut first) {
                    Ok(n) => break n == 0,
                    Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                    Err(_) => break false,
                }
            };
            if empty {
                return Err(IoError::empty_input(&input.name));
            }
            return Err(IoError::Detect {
                path: input.name.clone(),
                msg: format!(
                    "{} is {}, not FASTA; pass the FASTA file the intervals refer to",
                    input.name,
                    input.format.name()
                ),
            });
        }
        Ok(FastaStream {
            name: input.name,
            src: Box::new(input.reader),
            buf: vec![0; CAP],
            start: 0,
            end: 0,
            consumed_to: 0,
            eof: false,
            mid_line: false,
            line: 0,
            crlf: false,
            terminated: false,
        })
    }

    /// The number of the line of the item last returned, from 1.
    pub fn line(&self) -> u64 {
        self.line
    }

    /// Whether the line of the item last returned ended in CRLF.
    pub fn crlf(&self) -> bool {
        self.crlf
    }

    /// Whether the item last returned ended its line with a line feed. False
    /// for a piece of a long line and for a last line that has no line break.
    /// A line break that falls just after a 1 MiB piece comes as an empty
    /// `Seq` for which this is true.
    pub fn line_ended(&self) -> bool {
        self.terminated
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    fn fill(&mut self) -> Result<(), IoError> {
        if self.start > 0 {
            self.buf.copy_within(self.start..self.end, 0);
            self.end -= self.start;
            self.start = 0;
        }
        loop {
            match self.src.read(&mut self.buf[self.end..]) {
                Ok(0) => {
                    self.eof = true;
                    return Ok(());
                }
                Ok(n) => {
                    self.end += n;
                    return Ok(());
                }
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                Err(source) => {
                    return Err(IoError::Read {
                        path: self.name.clone(),
                        source,
                    });
                }
            }
        }
    }

    /// The next line, or `None` at the end.
    // Not `Iterator`: the item borrows the stream's buffer.
    #[allow(clippy::should_implement_trait)]
    pub fn next(&mut self) -> Result<Option<FastaLine<'_>>, IoError> {
        self.start = self.consumed_to;
        let (from, to, ended) = loop {
            let avail = &self.buf[self.start..self.end];
            if let Some(i) = memchr(b'\n', avail) {
                break (self.start, self.start + i, true);
            }
            if self.eof {
                if avail.is_empty() {
                    self.consumed_to = self.end;
                    return Ok(None);
                }
                break (self.start, self.end, false);
            }
            if avail.len() >= CAP {
                let starts_header = !self.mid_line && avail[0] == b'>';
                if starts_header {
                    return Err(IoError::Parse {
                        path: self.name.clone(),
                        loc: crate::error::Loc::Line(self.line + 1),
                        msg: "header line is longer than 1 MiB; this is not a FASTA header".into(),
                    });
                }
                // A piece of a long sequence line; a CR at its end may belong
                // to the line break, so it waits for the next piece.
                let mut to = self.end;
                if self.buf[to - 1] == b'\r' && to - self.start > 1 {
                    to -= 1;
                }
                break (self.start, to, false);
            }
            self.fill()?;
        };
        let content_end = if ended || self.eof {
            // Strip the CR of a CRLF break.
            if to > from && self.buf[to - 1] == b'\r' {
                to - 1
            } else {
                to
            }
        } else {
            to
        };
        let at_line_start = !self.mid_line;
        if at_line_start {
            self.line += 1;
        }
        let whole = ended || self.eof;
        self.crlf = whole && content_end < to;
        self.terminated = ended;
        self.consumed_to = if ended { to + 1 } else { to };
        self.mid_line = !whole;
        let data = &self.buf[from..content_end];
        Ok(Some(if at_line_start && data.first() == Some(&b'>') {
            FastaLine::Header(&data[1..])
        } else {
            FastaLine::Seq(data)
        }))
    }
}
