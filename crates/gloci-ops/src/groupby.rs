//! `gloci groupby` (D67): aggregate consecutive lines with equal key columns.
//!
//! Works on any tab-separated text, not only genomic records. Like bedtools,
//! the input needs no sort order: a key that comes back after other keys
//! starts a new group.

use std::io::{BufRead, Write};

use gloci_io::{Format, Header, IoError, Loc, UsageError, split_fields};

use crate::agg::{
    Accumulators, AggSpec, Col, header_names, parse_aggs, refuse_two_input_ops, resolve,
};

/// One key column as typed: a 1-based number or a header name.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KeyCol {
    Index(usize),
    Name(String),
}

/// Parse `-g 1,2` or `-g chrom,start` (names need `--header`).
pub fn parse_keys(text: &str) -> Result<Vec<KeyCol>, UsageError> {
    let mut out = Vec::new();
    for item in text.split(',').map(str::trim).filter(|s| !s.is_empty()) {
        out.push(if item.bytes().all(|b| b.is_ascii_digit()) {
            match item.parse::<usize>() {
                Ok(n) if n >= 1 => KeyCol::Index(n),
                _ => {
                    return Err(UsageError(format!(
                        "-g {item}: key columns are numbered from 1"
                    )));
                }
            }
        } else {
            KeyCol::Name(item.to_string())
        });
    }
    if out.is_empty() {
        return Err(UsageError(
            "-g needs at least one key column, such as -g 1,2".into(),
        ));
    }
    Ok(out)
}

pub struct GroupbyOptions {
    pub keys: Vec<KeyCol>,
    pub aggs: Vec<String>,
    /// The first line names the columns: `-g` and `-c` may use the names, and
    /// the output starts with a line naming its columns.
    pub header: bool,
    pub empty: String,
    pub delim: String,
}

/// Lines `groupby` passes over: blank or whitespace-only lines, `#` comments and UCSC
/// `track`/`browser` lines (all skipped by bedtools too).
fn is_skipped(line: &[u8]) -> bool {
    line.iter().all(u8::is_ascii_whitespace)
        || line[0] == b'#'
        || line.starts_with(b"track")
        || line.starts_with(b"browser")
}

/// Read the next line without its `\n` or `\r\n`; false at end of input.
fn next_line(
    r: &mut dyn BufRead,
    buf: &mut Vec<u8>,
    name: &str,
    n: &mut u64,
) -> Result<bool, IoError> {
    buf.clear();
    let got = r.read_until(b'\n', buf).map_err(|source| IoError::Read {
        path: name.to_string(),
        source,
    })?;
    if got == 0 {
        return Ok(false);
    }
    *n += 1;
    if buf.last() == Some(&b'\n') {
        buf.pop();
    }
    if buf.last() == Some(&b'\r') {
        buf.pop();
    }
    Ok(true)
}

fn write_err(dest: &str) -> impl Fn(std::io::Error) -> IoError + '_ {
    move |source| IoError::Write {
        path: dest.to_string(),
        source,
    }
}

/// Group `input` (named `name` in messages) and write one line per group to `out`.
pub fn groupby(
    input: &mut dyn BufRead,
    name: &str,
    o: &GroupbyOptions,
    out: &mut dyn Write,
    dest: &str,
) -> Result<(), IoError> {
    let specs = parse_aggs(&o.aggs)?;
    if specs.is_empty() {
        return Err(UsageError("-c needs at least one COL:OP, such as -c 5:sum".into()).into());
    }
    refuse_two_input_ops(&specs, "groupby")?;
    if !o.header {
        let named = o.keys.iter().find_map(|k| match k {
            KeyCol::Name(n) => Some(n.clone()),
            KeyCol::Index(_) => None,
        });
        let named = named.or_else(|| {
            specs.iter().find_map(|s| match &s.col {
                Col::Name(n) => Some(n.clone()),
                _ => None,
            })
        });
        if let Some(n) = named {
            return Err(UsageError(format!(
                "column name `{n}` needs --header, which reads the names from the first line; \
                 use a column number, or pass --header"
            ))
            .into());
        }
    }

    let mut buf = Vec::new();
    let mut line_no = 0u64;
    let mut header_line: Option<Vec<u8>> = None;
    if o.header {
        if next_line(input, &mut buf, name, &mut line_no)? {
            header_line = Some(buf.clone());
        } else {
            return Ok(());
        }
    }
    let header = match &header_line {
        Some(h) => {
            let mut v = h.clone();
            v.push(b'\n');
            Header::Lines(v)
        }
        None => Header::None,
    };
    let names = header_names(&header);
    let key_idx: Vec<usize> = o
        .keys
        .iter()
        .map(|k| match k {
            KeyCol::Index(n) => Ok(n - 1),
            KeyCol::Name(n) => names.iter().position(|x| x == n).ok_or_else(|| {
                UsageError(format!(
                    "-g {n}: no column named `{n}` in the header of {name}; available names: {}",
                    names.join(", ")
                ))
            }),
        })
        .collect::<Result<_, _>>()?;
    let resolved = resolve(&specs, Format::Bed, &header)?;
    let mut acc = Accumulators::new(resolved, o.empty.clone(), o.delim.clone());

    if let Some(h) = &header_line {
        out.write_all(&header_row(h, &key_idx, &specs))
            .map_err(write_err(dest))?;
    }

    let need = key_idx.iter().max().map_or(0, |m| m + 1);
    let mut key: Vec<u8> = Vec::new();
    let mut cur: Vec<u8> = Vec::new();
    let mut in_group = false;
    let mut first_data_line = None;
    let mut row: Vec<u8> = Vec::new();
    while next_line(input, &mut buf, name, &mut line_no)? {
        if is_skipped(&buf) {
            continue;
        }
        key.clear();
        let fields: Vec<&[u8]> = split_fields(&buf).collect();
        if fields.len() < need {
            return Err(IoError::Parse {
                path: name.to_string(),
                loc: Loc::Line(line_no),
                msg: format!(
                    "the line has {} tab-separated column{} but -g needs column {need}; \
                     check the -g numbers (columns count from 1) and that the file is \
                     tab-separated",
                    fields.len(),
                    if fields.len() == 1 { "" } else { "s" }
                ),
            });
        }
        for (i, &k) in key_idx.iter().enumerate() {
            if i > 0 {
                key.push(b'\t');
            }
            key.extend_from_slice(fields[k]);
        }
        if in_group && key != cur {
            flush(&mut acc, &cur, &mut row, out, dest)?;
        }
        if !in_group || key != cur {
            std::mem::swap(&mut cur, &mut key);
            in_group = true;
        }
        let first = *first_data_line.get_or_insert(line_no);
        acc.add(&buf, 0, name, line_no).map_err(|e| match e {
            IoError::Parse {
                path,
                loc: Loc::Line(n),
                msg,
            } if n == first && !o.header && msg.contains("is not a number") => IoError::Parse {
                path,
                loc: Loc::Line(n),
                msg: format!("{msg}; if the first line is a header, pass --header"),
            },
            e => e,
        })?;
    }
    if in_group {
        flush(&mut acc, &cur, &mut row, out, dest)?;
    }
    Ok(())
}

fn flush(
    acc: &mut Accumulators,
    key: &[u8],
    row: &mut Vec<u8>,
    out: &mut dyn Write,
    dest: &str,
) -> Result<(), IoError> {
    row.clear();
    row.extend_from_slice(key);
    acc.finish_into(row);
    row.push(b'\n');
    out.write_all(row).map_err(write_err(dest))
}

/// The key names from the header, then the `-c` labels as typed.
fn header_row(header: &[u8], key_idx: &[usize], specs: &[AggSpec]) -> Vec<u8> {
    let hash = header.starts_with(b"#");
    let text = header.strip_prefix(b"#").unwrap_or(header);
    let fields: Vec<&[u8]> = split_fields(text).collect();
    let mut row = Vec::new();
    if hash {
        row.push(b'#');
    }
    for (i, &k) in key_idx.iter().enumerate() {
        if i > 0 {
            row.push(b'\t');
        }
        // an unnamed key column (a header with fewer names) is written as its number
        match fields.get(k) {
            Some(f) => row.extend_from_slice(f),
            None => row.extend_from_slice((k + 1).to_string().as_bytes()),
        }
    }
    for s in specs {
        row.push(b'\t');
        row.extend_from_slice(s.label.as_bytes());
    }
    row.push(b'\n');
    row
}
