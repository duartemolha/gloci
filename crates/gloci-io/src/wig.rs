//! WIG reader (D63): `fixedStep` and `variableStep` sections, read as
//! bedGraph records (chrom, start, end, value).
//!
//! WIG positions are 1-based and a data line covers `span` bases (default 1),
//! so a line at position `p` becomes the 0-based interval `[p-1, p-1+span)`.
//! Each data line gives one record; equal neighbours are not merged. The
//! chromosome is named once, on the step line, so `--chrom-map` renames it
//! there and the records carry the new name.

use crate::chrom_map::{ChromMap, Renamer};
use crate::detect::{Compression, Format};
use crate::diag::{Diagnostics, InvalidPolicy};
use crate::error::{IoError, Loc};
use crate::input::Input;
use crate::record::{
    Header, RecordBatch, RecordMeta, RecordReader, coord, intern_chrom, warn_undeclared,
};
use gloci_core::{ChromId, ChromTable, Interval, Strand};
use std::io::{BufRead, Write};
use std::sync::Arc;

/// Is `line` a `fixedStep` or `variableStep` declaration?
pub fn is_step_line(line: &[u8]) -> bool {
    step_kind(line).is_some()
}

fn step_kind(line: &[u8]) -> Option<bool> {
    for (word, fixed) in [(&b"fixedStep"[..], true), (&b"variableStep"[..], false)] {
        if line.starts_with(word) && line.get(word.len()).is_some_and(u8::is_ascii_whitespace) {
            return Some(fixed);
        }
    }
    None
}

fn is_header_line(line: &[u8]) -> bool {
    let word =
        |w: &[u8]| line.starts_with(w) && line.get(w.len()).is_none_or(u8::is_ascii_whitespace);
    line.first() == Some(&b'#') || word(b"track") || word(b"browser")
}

/// The section a data line belongs to.
struct Section {
    chrom: ChromId,
    /// The (renamed) chromosome name, as written in each record.
    name: Vec<u8>,
    span: u64,
    /// `fixedStep`: the 0-based start of the next data line and the step.
    fixed: Option<(u64, u64)>,
}

pub(crate) struct WigReader {
    name: String,
    compression: Compression,
    header: Header,
    diag: Diagnostics,
    policy: InvalidPolicy,
    src: Box<dyn BufRead + Send>,
    line: u64,
    pending: Option<Vec<u8>>,
    eof: bool,
    section: Option<Section>,
    renamer: Option<Renamer>,
}

impl WigReader {
    pub(crate) fn open(
        input: Input,
        policy: InvalidPolicy,
        map: Option<Arc<ChromMap>>,
    ) -> Result<Self, IoError> {
        let mut r = WigReader {
            name: input.name,
            compression: input.compression,
            header: Header::None,
            diag: Diagnostics::default(),
            policy,
            src: input.reader,
            line: 0,
            pending: None,
            eof: false,
            section: None,
            renamer: map.map(Renamer::new),
        };
        let mut head = Vec::new();
        let mut line = Vec::new();
        loop {
            line.clear();
            if !r.read_line(&mut line)? {
                r.eof = true;
                break;
            }
            if line.trim_ascii().is_empty() {
                r.diag.record_blank(r.line);
            } else if is_header_line(&line) {
                head.extend_from_slice(&line);
                head.push(b'\n');
            } else {
                r.line -= 1;
                r.pending = Some(line);
                break;
            }
        }
        if !head.is_empty() {
            r.header = Header::Lines(head);
        }
        Ok(r)
    }

    fn err(&self, msg: String) -> IoError {
        IoError::Parse {
            path: self.name.clone(),
            loc: Loc::Line(self.line),
            msg,
        }
    }

    /// Replace `out` with the next line, terminator stripped.
    fn read_line(&mut self, out: &mut Vec<u8>) -> Result<bool, IoError> {
        if let Some(p) = self.pending.take() {
            *out = p;
            self.line += 1;
            return Ok(true);
        }
        let n = self
            .src
            .read_until(b'\n', out)
            .map_err(|source| IoError::Read {
                path: self.name.clone(),
                source,
            })?;
        if n == 0 {
            return Ok(false);
        }
        self.line += 1;
        if out.last() == Some(&b'\n') {
            out.pop();
        }
        if out.last() == Some(&b'\r') {
            out.pop();
        }
        Ok(true)
    }

    /// Start the section a step line declares.
    fn declare(&mut self, line: &[u8], fixed: bool, chroms: &mut ChromTable) -> Result<(), String> {
        let kind = if fixed { "fixedStep" } else { "variableStep" };
        let (mut chrom, mut start, mut step, mut span) = (None, None, None, None);
        for tok in line.split(|b| b.is_ascii_whitespace()).skip(1) {
            if tok.is_empty() {
                continue;
            }
            let Some(eq) = tok.iter().position(|&b| b == b'=') else {
                return Err(format!(
                    "{kind} line has '{}' where key=value was expected; \
write for example `{kind} chrom=chr1 {}`",
                    String::from_utf8_lossy(tok),
                    if fixed { "start=1 step=1" } else { "span=1" }
                ));
            };
            let val = &tok[eq + 1..];
            match &tok[..eq] {
                b"chrom" => chrom = Some(val),
                b"start" => start = Some(coord("start", val)?),
                b"step" => step = Some(coord("step", val)?),
                b"span" => span = Some(coord("span", val)?),
                // Unknown keys are ignored, as UCSC tools do.
                _ => {}
            }
        }
        let chrom = chrom.ok_or_else(|| {
            format!("{kind} line has no chrom=NAME; add it, for example `{kind} chrom=chr1`")
        })?;
        let span = span.unwrap_or(1);
        if span == 0 {
            return Err("span must be at least 1; use span=1 for single bases".into());
        }
        let fixed = if fixed {
            let start =
                start.ok_or("fixedStep line has no start=POSITION; add start=1 (1-based)")?;
            if start == 0 {
                return Err("start is 0 but WIG positions are 1-based; use start=1".into());
            }
            let step = step.ok_or("fixedStep line has no step=N; add step=1")?;
            if step == 0 {
                return Err("step must be at least 1".into());
            }
            Some((start - 1, step))
        } else {
            None
        };
        // Rename the chromosome once per section.
        let mut name = chrom.to_vec();
        if let Some(r) = &mut self.renamer {
            r.apply(&mut name, 0);
        }
        let id = intern_chrom(chroms, &name)?;
        self.section = Some(Section {
            chrom: id,
            name,
            span,
            fixed,
        });
        Ok(())
    }

    /// Turn one data line into a record's bytes and interval.
    fn data(&mut self, line: &[u8], out: &mut Vec<u8>) -> Result<(ChromId, Interval), String> {
        let sec = self.section.as_mut().expect("section checked by caller");
        let mut toks = line
            .split(|b| b.is_ascii_whitespace())
            .filter(|t| !t.is_empty());
        let (start, value) = if let Some((next, step)) = &mut sec.fixed {
            let value = toks.next().ok_or("empty data line")?;
            if toks.next().is_some() {
                return Err(
                    "a fixedStep data line is one value; a position and value belong to variableStep"
                        .into(),
                );
            }
            let start = *next;
            // Advance even when this line turns out invalid, so later lines
            // keep their positions.
            *next = next.saturating_add(*step);
            (start, value)
        } else {
            let (Some(pos), Some(value)) = (toks.next(), toks.next()) else {
                return Err(
                    "a variableStep data line is `POSITION VALUE`, separated by whitespace".into(),
                );
            };
            if toks.next().is_some() {
                return Err("a variableStep data line has more than POSITION and VALUE".into());
            }
            let pos = coord("position", pos)?;
            if pos == 0 {
                return Err("position is 0 but WIG positions are 1-based".into());
            }
            (pos - 1, value)
        };
        if std::str::from_utf8(value)
            .ok()
            .and_then(|v| v.parse::<f64>().ok())
            .is_none_or(|v| !v.is_finite())
        {
            return Err(format!(
                "value '{}' is not a finite number; WIG values must be finite, \
so remove or replace nan/inf",
                String::from_utf8_lossy(value)
            ));
        }
        let end = start
            .checked_add(sec.span)
            .ok_or("position plus span overflows")?;
        out.extend_from_slice(&sec.name);
        write!(out, "\t{start}\t{end}\t").expect("write to Vec");
        out.extend_from_slice(value);
        Ok((sec.chrom, Interval { start, end }))
    }
}

impl RecordReader for WigReader {
    fn name(&self) -> &str {
        &self.name
    }
    fn format(&self) -> Format {
        Format::Wig
    }
    fn compression(&self) -> Compression {
        self.compression
    }
    fn header(&self) -> &Header {
        &self.header
    }
    fn diagnostics(&self) -> &Diagnostics {
        &self.diag
    }

    fn read_batch(
        &mut self,
        chroms: &mut ChromTable,
        batch: &mut RecordBatch,
        max: usize,
    ) -> Result<bool, IoError> {
        debug_assert!(max > 0, "read_batch max must be > 0");
        let mut n = 0;
        let mut line = Vec::new();
        while n < max {
            if self.eof {
                return Ok(false);
            }
            line.clear();
            if !self.read_line(&mut line)? {
                self.eof = true;
                return Ok(false);
            }
            if line.trim_ascii().is_empty() {
                self.diag.record_blank(self.line);
                continue;
            }
            if is_header_line(&line) {
                continue;
            }
            if let Some(fixed) = step_kind(&line) {
                self.declare(&line, fixed, chroms)
                    .map_err(|m| self.err(m))?;
                continue;
            }
            if self.section.is_none() {
                let shown: String = String::from_utf8_lossy(&line).chars().take(60).collect();
                return Err(self.err(format!(
                    "data line before any fixedStep or variableStep line: {shown}; \
start the data with a line like `fixedStep chrom=chr1 start=1 step=1` or \
`variableStep chrom=chr1`"
                )));
            }
            let start = batch.buf.len();
            match self.data(&line, &mut batch.buf) {
                Ok((chrom, span)) => {
                    let len = batch.buf.len() - start;
                    batch.recs.push(RecordMeta {
                        chrom,
                        key: span.start,
                        span,
                        strand: Strand::Unknown,
                        loc: self.line,
                        off: start,
                        len: len as u32,
                    });
                    n += 1;
                    if chroms.has_undeclared() {
                        warn_undeclared(chroms, &self.name, Some(Loc::Line(self.line)));
                    }
                }
                Err(msg) => {
                    batch.buf.truncate(start);
                    match self.policy {
                        InvalidPolicy::Fail => return Err(self.err(msg)),
                        InvalidPolicy::Skip => self.diag.record(self.line, &msg),
                    }
                }
            }
        }
        Ok(true)
    }
}
