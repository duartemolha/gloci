//! The `vcf query` format language (D172): the format string of bcftools 1.21
//! `query -f`, parsed on its own (`parse`) and then bound to a header's fields
//! and declared tags (`resolve`).

mod parse;
mod resolve;

use gloci_io::{UsageError, VcfHeader};

/// A parsed format string: its source text and syntax items, before any header.
#[derive(Clone, Debug)]
pub struct QueryFormat {
    src: String,
    items: Vec<parse::Item>,
}

impl QueryFormat {
    /// Syntax only, no header: errors are usage errors that quote the format
    /// and give a 1-based column.
    pub fn parse(src: &str) -> Result<QueryFormat, UsageError> {
        Ok(QueryFormat {
            src: src.to_string(),
            items: parse::parse(src)?,
        })
    }

    /// Bind names to fields and declared tags; `print_header` is `-H`
    /// (refused with `%LINE`).
    pub fn resolve(&self, header: &VcfHeader, print_header: bool) -> Result<Query, UsageError> {
        resolve::resolve(&self.src, &self.items, header, print_header)
    }
}

/// A resolved format: what to write for each record.
#[derive(Clone, Debug, PartialEq)]
pub struct Query {
    pub pieces: Vec<Piece>,
}

/// One part of a format.
#[derive(Clone, Debug, PartialEq)]
pub enum Piece {
    Text(Vec<u8>),
    Field(Field),
    /// A `[ ]` block, written once per printed sample.
    Block(Vec<Piece>),
}

/// A field a format can name.
#[derive(Clone, Debug, PartialEq)]
pub enum Field {
    Chrom,
    Pos,
    Pos0,
    End,
    End0,
    Id,
    Ref,
    Alt(Option<usize>),
    FirstAlt,
    Qual,
    Filter,
    Type,
    /// The INFO column.
    Info,
    InfoTag(Tag),
    Line,
    Sample,
    Gt,
    Tgt,
    FormatTag(Tag),
}

/// A declared tag with its optional `{N}` subscript.
#[derive(Clone, Debug, PartialEq)]
pub struct Tag {
    pub name: Vec<u8>,
    pub ty: TagType,
    pub sub: Option<usize>,
}

/// A tag's declared type; Character and a missing Type are `String`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TagType {
    Integer,
    Float,
    String,
    Flag,
}
