//! Arrow output (D89, R3): `to_arrow()`, `to_polars()` and `to_pandas()` take the records the
//! output reads back as ([`crate::records`]) as one table, built in memory.
//!
//! The sink ([`ArrowSink`]) buffers each column's text as the records come, on the host
//! stage's thread; once they have all come, it infers each column's type and builds the
//! columns with arrow-array's builders. The table is handed to Python through the Arrow
//! PyCapsule interface ([`ArrowStream`]), so polars needs no pyarrow.
//!
//! **Columns.** `chrom start end` first, for records with coordinates (`start`/`end` int64;
//! null where a record has none), then:
//! - VCF: `pos id ref alt qual filter info format` and the samples (`qual` is float64, `.`
//!   being null);
//! - GFF: `source feature score strand frame attributes`;
//! - SAM: `qname flag pos mapq cigar rnext pnext tlen seq qual`, then the tags;
//! - FASTA: `name sequence` (no coordinates);
//! - BED and others: each column named by a `#` header line with as many fields, else
//!   `col<N>`, N being the column's 1-based number in the record;
//! - a table: named by its header row, else as BED's are (and no coordinates).
//!
//! **Types** are inferred per column: int64 when every value is an integer, else float64 when
//! every value is a number, else string. `.` is a missing value for a number (null), and a
//! value of a string column. A field a row does not have is null.

use crate::records::{Layout, Row, RowSink};
use arrow_array::builder::{Float64Builder, Int64Builder, StringBuilder};
use arrow_array::ffi_stream::FFI_ArrowArrayStream;
use arrow_array::{ArrayRef, RecordBatch, RecordBatchIterator};
use arrow_schema::{DataType, Field, Schema, SchemaRef};
use pyo3::prelude::*;
use pyo3::types::PyCapsule;
use std::io;
use std::sync::{Arc, Mutex};

/// The most text a record batch holds, so that no string column overflows Arrow's 32-bit
/// offsets.
const CHUNK_BYTES: usize = 1 << 30;

/// Column names of VCF's fields 2 to 9 (CHROM is `chrom`).
const VCF: [&str; 8] = [
    "pos", "id", "ref", "alt", "qual", "filter", "info", "format",
];
/// GFF's fields 2, 3 and 6 to 9 (seqid, start and end are `chrom start end`).
const GFF: [(usize, &str); 6] = [
    (1, "source"),
    (2, "feature"),
    (5, "score"),
    (6, "strand"),
    (7, "frame"),
    (8, "attributes"),
];
/// SAM's fields 1, 2 and 4 to 11 (RNAME is `chrom`).
const SAM: [(usize, &str); 10] = [
    (0, "qname"),
    (1, "flag"),
    (3, "pos"),
    (4, "mapq"),
    (5, "cigar"),
    (6, "rnext"),
    (7, "pnext"),
    (8, "tlen"),
    (9, "seq"),
    (10, "qual"),
];

/// The values of one column, as text, row by row.
#[derive(Default)]
struct TextColumn {
    data: Vec<u8>,
    lens: Vec<u32>,
    /// The rows (in order) with no value: the record has no such field.
    absent: Vec<usize>,
    /// Every value seen is an integer, a number; and at least one is (`.` aside).
    all_int: bool,
    all_float: bool,
    any_number: bool,
}

impl TextColumn {
    fn new(rows_before: usize) -> TextColumn {
        TextColumn {
            data: Vec::new(),
            lens: vec![0; rows_before],
            absent: (0..rows_before).collect(),
            all_int: true,
            all_float: true,
            any_number: false,
        }
    }

    fn push(&mut self, v: &[u8]) {
        self.data.extend_from_slice(v);
        self.lens.push(v.len() as u32);
        if v != b"." && (self.all_int || self.all_float) {
            if self.all_int && int_value(v).is_none() {
                self.all_int = false;
            }
            if self.all_float && float_value(v).is_none() {
                self.all_float = false;
            }
            self.any_number = true;
        }
    }

    fn push_absent(&mut self) {
        self.absent.push(self.lens.len());
        self.lens.push(0);
    }

    fn kind(&self) -> Kind {
        if self.any_number && self.all_int {
            Kind::Int
        } else if self.any_number && self.all_float {
            Kind::Float
        } else {
            Kind::Text
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    Int,
    Float,
    Text,
}

/// An integer as text writes it: an optional `-`, then digits with no leading zero.
fn int_value(v: &[u8]) -> Option<i64> {
    let digits = v.strip_prefix(b"-").unwrap_or(v);
    let plain = !digits.is_empty()
        && digits.iter().all(u8::is_ascii_digit)
        && (digits.len() == 1 || digits[0] != b'0');
    if !plain {
        return None;
    }
    std::str::from_utf8(v).ok()?.parse().ok()
}

/// A number: what Rust reads as an f64 (`0.5`, `1e-8`, `inf`, `nan`), but not with a `+` or a
/// zero-padded integer part (`007.5`), which would not be written back the same.
fn float_value(v: &[u8]) -> Option<f64> {
    let body = v.strip_prefix(b"-").unwrap_or(v);
    if body.first() == Some(&b'+') || body.first() == Some(&b'-') {
        return None;
    }
    if body.len() > 1 && body[0] == b'0' && body[1].is_ascii_digit() {
        return None;
    }
    std::str::from_utf8(v).ok()?.parse().ok()
}

/// Reads a column's values in order.
struct Cursor<'c> {
    col: &'c TextColumn,
    row: usize,
    at: usize,
    absent: usize,
}

impl<'c> Cursor<'c> {
    fn new(col: &'c TextColumn) -> Cursor<'c> {
        Cursor {
            col,
            row: 0,
            at: 0,
            absent: 0,
        }
    }

    fn next(&mut self) -> Option<&'c [u8]> {
        let len = self.col.lens[self.row] as usize;
        let is_absent = self.col.absent.get(self.absent) == Some(&self.row);
        self.row += 1;
        if is_absent {
            self.absent += 1;
            return None;
        }
        let v = &self.col.data[self.at..self.at + len];
        self.at += len;
        Some(v)
    }
}

/// Whether field `i` of a record laid out as `layout` is a column of its own (not `chrom`,
/// `start` or `end`).
fn is_column(layout: Layout, i: usize) -> bool {
    !match layout {
        Layout::Bed => i < 3,
        Layout::Vcf => i == 0,
        Layout::Gff => matches!(i, 0 | 3 | 4),
        Layout::Sam => i == 2,
        Layout::Fasta | Layout::Table => false,
    }
}

/// The table being built from the records (R3).
pub(crate) struct ArrowSink {
    layout: Layout,
    /// The last `#` line before the first record, and a table's header row.
    hash_line: Option<Vec<u8>>,
    header_row: Option<Vec<u8>>,
    rows: usize,
    chrom: TextColumn,
    start: Vec<Option<i64>>,
    end: Vec<Option<i64>>,
    /// The record's fields that are columns, by field index (`None`: the field is `chrom`,
    /// `start` or `end`).
    cols: Vec<Option<TextColumn>>,
    /// Where a record batch ends, by row.
    chunks: Vec<usize>,
    chunk_bytes: usize,
    out: Arc<Mutex<Option<Table>>>,
}

impl ArrowSink {
    pub(crate) fn new(out: Arc<Mutex<Option<Table>>>) -> ArrowSink {
        ArrowSink {
            layout: Layout::Table,
            hash_line: None,
            header_row: None,
            rows: 0,
            chrom: TextColumn::new(0),
            start: Vec::new(),
            end: Vec::new(),
            cols: Vec::new(),
            chunks: Vec::new(),
            chunk_bytes: 0,
            out,
        }
    }

    fn coords(&self) -> bool {
        !matches!(self.layout, Layout::Table | Layout::Fasta)
    }

    fn is_column(&self, i: usize) -> bool {
        is_column(self.layout, i)
    }

    /// The column names: the fixed ones of the layout, then those of the header row or `#`
    /// line, else `col<N>`.
    fn column_names(&self) -> Vec<String> {
        let split = |l: &Vec<u8>| -> Vec<String> {
            let l = l.strip_prefix(b"#").unwrap_or(l);
            let l = l.strip_prefix(b"#").unwrap_or(l);
            l.split(|&b| b == b'\t')
                .map(|f| String::from_utf8_lossy(f).into_owned())
                .collect()
        };
        let n = self.cols.len();
        let header: Option<Vec<String>> = match self.layout {
            Layout::Table => self
                .header_row
                .as_ref()
                .map(split)
                .or_else(|| self.hash_line.as_ref().map(split).filter(|h| h.len() == n)),
            Layout::Vcf => self.hash_line.as_ref().map(split),
            Layout::Bed => self.hash_line.as_ref().map(split).filter(|h| h.len() == n),
            _ => None,
        };
        (0..n)
            .filter(|&i| self.is_column(i))
            .map(|i| {
                let fixed = match self.layout {
                    Layout::Vcf if (1..=8).contains(&i) => Some(VCF[i - 1]),
                    Layout::Gff => GFF.iter().find(|(f, _)| *f == i).map(|(_, n)| *n),
                    Layout::Sam => SAM.iter().find(|(f, _)| *f == i).map(|(_, n)| *n),
                    Layout::Fasta => Some(["name", "sequence"][i]),
                    _ => None,
                };
                fixed
                    .map(str::to_string)
                    .or_else(|| {
                        // VCF: the samples' names come from `#CHROM`.
                        let from_header = self.layout != Layout::Vcf || i >= 9;
                        header
                            .as_ref()
                            .filter(|_| from_header)
                            .and_then(|h| h.get(i))
                            .filter(|s| !s.is_empty())
                            .cloned()
                    })
                    .unwrap_or_else(|| format!("col{}", i + 1))
            })
            .collect()
    }

    fn build(&mut self) -> Table {
        // The columns a header names, whether or not any record has them (a VCF with no
        // records still has its columns).
        let named = match self.layout {
            Layout::Vcf => self.hash_line.as_ref(),
            Layout::Table => self.header_row.as_ref(),
            _ => None,
        }
        .map_or(0, |l| l.split(|&b| b == b'\t').count());
        while self.cols.len() < named {
            let i = self.cols.len();
            let rows = self.rows;
            self.cols
                .push(is_column(self.layout, i).then(|| TextColumn::new(rows)));
        }
        let mut fields = Vec::new();
        let names = self.column_names();
        if self.coords() {
            fields.push(Field::new("chrom", DataType::Utf8, true));
            fields.push(Field::new("start", DataType::Int64, true));
            fields.push(Field::new("end", DataType::Int64, true));
        }
        let (pos, qual) = (1, 5); // VCF's POS and QUAL fields
        let kinds: Vec<(usize, Kind)> = (0..self.cols.len())
            .filter(|&i| self.is_column(i))
            .map(|i| {
                let col = self.cols[i].as_ref().expect("a column");
                let kind = match (self.layout, col.kind()) {
                    (Layout::Vcf, _) if i == pos => Kind::Int,
                    // QUAL is a number, `.` being null (D89), whenever every value is one.
                    (Layout::Vcf, k) if i == qual && k != Kind::Text => Kind::Float,
                    (Layout::Vcf, Kind::Text) if i == qual && col.all_float => Kind::Float,
                    (_, k) => k,
                };
                (i, kind)
            })
            .collect();
        for ((_, kind), name) in kinds.iter().zip(&names) {
            let t = match kind {
                Kind::Int => DataType::Int64,
                Kind::Float => DataType::Float64,
                Kind::Text => DataType::Utf8,
            };
            fields.push(Field::new(name, t, true));
        }
        let schema: SchemaRef = Arc::new(Schema::new(fields));

        let mut bounds = std::mem::take(&mut self.chunks);
        if bounds.last() != Some(&self.rows) {
            bounds.push(self.rows);
        }
        let mut chrom = Cursor::new(&self.chrom);
        let mut cursors: Vec<Cursor<'_>> = kinds
            .iter()
            .map(|(i, _)| Cursor::new(self.cols[*i].as_ref().expect("a column")))
            .collect();
        let mut batches = Vec::new();
        let mut from = 0;
        for to in bounds {
            if to == from {
                continue;
            }
            let n = to - from;
            let mut arrays: Vec<ArrayRef> = Vec::with_capacity(schema.fields().len());
            if self.coords() {
                arrays.push(text_array(&mut chrom, n));
                let ints = |v: &[Option<i64>]| -> ArrayRef {
                    let mut b = Int64Builder::with_capacity(n);
                    b.extend(v[from..to].iter().copied());
                    Arc::new(b.finish())
                };
                arrays.push(ints(&self.start));
                arrays.push(ints(&self.end));
            }
            for (c, (_, kind)) in cursors.iter_mut().zip(&kinds) {
                arrays.push(match kind {
                    Kind::Text => text_array(c, n),
                    Kind::Int => {
                        let mut b = Int64Builder::with_capacity(n);
                        for _ in 0..n {
                            b.append_option(c.next().and_then(int_value));
                        }
                        Arc::new(b.finish())
                    }
                    Kind::Float => {
                        let mut b = Float64Builder::with_capacity(n);
                        for _ in 0..n {
                            b.append_option(c.next().and_then(float_value));
                        }
                        Arc::new(b.finish())
                    }
                });
            }
            batches.push(
                RecordBatch::try_new(schema.clone(), arrays).expect("the columns match the schema"),
            );
            from = to;
        }
        Table { schema, batches }
    }
}

/// `n` values of a string column; an absent one is null.
fn text_array(c: &mut Cursor<'_>, n: usize) -> ArrayRef {
    let mut b = StringBuilder::with_capacity(n, 0);
    for _ in 0..n {
        match c.next() {
            None => b.append_null(),
            Some(v) => match std::str::from_utf8(v) {
                Ok(s) => b.append_value(s),
                Err(_) => b.append_value(String::from_utf8_lossy(v)),
            },
        }
    }
    Arc::new(b.finish())
}

impl RowSink for ArrowSink {
    fn layout(&mut self, layout: Layout) {
        self.layout = layout;
    }

    fn names(&mut self, line: &[u8], header_row: bool) {
        if header_row {
            self.header_row = Some(line.to_vec());
        } else {
            self.hash_line = Some(line.to_vec());
        }
    }

    fn row(&mut self, row: &Row<'_>) -> io::Result<()> {
        let r = self.rows;
        let layout = self.layout;
        let cols = &mut self.cols;
        let mut fields = 0;
        let mut chrom = None;
        row.for_each_field(|i, v| {
            fields = i + 1;
            if Some(i) == row.chrom {
                chrom = Some(v.to_vec());
            }
            if i >= cols.len() {
                cols.resize_with(i + 1, || None);
            }
            if is_column(layout, i) {
                cols[i].get_or_insert_with(|| TextColumn::new(r)).push(v);
            }
        });
        // The columns this record has no field for.
        for col in self.cols[fields..].iter_mut().flatten() {
            col.push_absent();
        }
        if self.coords() {
            match chrom {
                Some(c) => self.chrom.push(&c),
                None => self.chrom.push_absent(),
            }
            let span = row.span.map(|(s, e)| (s as i64, e as i64));
            self.start.push(span.map(|s| s.0));
            self.end.push(span.map(|s| s.1));
        }
        self.rows += 1;
        self.chunk_bytes += row.line.len() + 1;
        if self.chunk_bytes >= CHUNK_BYTES {
            self.chunks.push(self.rows);
            self.chunk_bytes = 0;
        }
        Ok(())
    }

    fn finish(&mut self) -> anyhow::Result<()> {
        // A field index no record has (a gap) cannot happen: fields are numbered from 0.
        let table = self.build();
        *self.out.lock().unwrap_or_else(|p| p.into_inner()) = Some(table);
        Ok(())
    }
}

/// A built table.
pub(crate) struct Table {
    schema: SchemaRef,
    batches: Vec<RecordBatch>,
}

impl Table {
    /// The table with its columns named `names`, one per column.
    pub(crate) fn renamed(self, names: &[String]) -> Result<Table, String> {
        let have: Vec<&str> = self
            .schema
            .fields()
            .iter()
            .map(|f| f.name().as_str())
            .collect();
        if names.len() != have.len() {
            return Err(format!(
                "names= has {} names, but the output has {} columns ({}); give one name per \
column",
                names.len(),
                have.len(),
                have.join(", ")
            ));
        }
        let schema: SchemaRef = Arc::new(Schema::new(
            self.schema
                .fields()
                .iter()
                .zip(names)
                .map(|(f, n)| f.as_ref().clone().with_name(n))
                .collect::<Vec<_>>(),
        ));
        let batches = self
            .batches
            .into_iter()
            .map(|b| RecordBatch::try_new(schema.clone(), b.columns().to_vec()))
            .collect::<Result<_, _>>()
            .map_err(|e| e.to_string())?;
        Ok(Table { schema, batches })
    }
}

/// A table of a pipeline's output, built in memory (D89): an Arrow C stream, through the
/// PyCapsule interface (`__arrow_c_stream__`), which pyarrow, polars and pandas take without
/// copying. `Pipeline.to_arrow_stream()` makes it; it can be exported any number of times.
#[pyclass(frozen, module = "gloci", name = "ArrowStream")]
pub(crate) struct ArrowStream {
    table: Table,
}

impl ArrowStream {
    pub(crate) fn new(table: Table) -> ArrowStream {
        ArrowStream { table }
    }
}

#[pymethods]
impl ArrowStream {
    /// A new Arrow C stream of the table, in a PyCapsule named `arrow_array_stream`. A
    /// requested schema is not followed: the table has its own types (D89).
    #[pyo3(signature = (requested_schema = None))]
    fn __arrow_c_stream__<'py>(
        &self,
        py: Python<'py>,
        requested_schema: Option<Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyCapsule>> {
        let _ = requested_schema;
        let reader = RecordBatchIterator::new(
            self.table.batches.clone().into_iter().map(Ok),
            self.table.schema.clone(),
        );
        let stream = FFI_ArrowArrayStream::new(Box::new(reader));
        // The capsule owns the stream: a consumer moves it out (leaving it released), and the
        // capsule's destructor releases one nobody took.
        PyCapsule::new_with_value(py, stream, c"arrow_array_stream")
    }

    /// The number of rows.
    #[getter]
    fn num_rows(&self) -> usize {
        self.table.batches.iter().map(RecordBatch::num_rows).sum()
    }

    /// The column names, in order.
    #[getter]
    fn column_names(&self) -> Vec<String> {
        self.table
            .schema
            .fields()
            .iter()
            .map(|f| f.name().clone())
            .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "<gloci.ArrowStream: {} rows; {}>",
            self.num_rows(),
            self.table
                .schema
                .fields()
                .iter()
                .map(|f| format!("{}: {}", f.name(), f.data_type()))
                .collect::<Vec<_>>()
                .join(", ")
        )
    }
}
