//! Arrow sources (D88): an object with `__arrow_c_stream__` (a pyarrow Table or
//! RecordBatchReader, a polars or pandas DataFrame) is read through the Arrow C stream interface
//! with arrow-rs, so reading it needs no Python package (spec 8.5, R3).
//!
//! Each row is a line: the integer columns `chrom start end` first, then the other columns in
//! their order, tab-separated, with a null written `.`. Values are written as Python's `str()`
//! writes them (`True`, `0.1`, `1e+16`), as a tuple's are. A null `chrom`, `start` or `end`, a
//! negative coordinate, or a tab or line break in a value is a `DataError` naming the column
//! and the 1-based row. pandas' index columns (`__index_level_0__`) are left out.

use crate::errors;
use crate::sources::{Fill, breaks, push_float};
use arrow_array::cast::AsArray;
use arrow_array::ffi_stream::{FFI_ArrowArrayStream, RecordBatchStreamReader};
use arrow_array::types::{
    ArrowDictionaryKeyType, Float16Type, Float32Type, Float64Type, Int8Type, Int16Type, Int32Type,
    Int64Type, UInt8Type, UInt16Type, UInt32Type, UInt64Type,
};
use arrow_array::{Array, ArrowPrimitiveType, RecordBatch, RecordBatchReader, new_empty_array};
use arrow_schema::{DataType, Schema};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::PyCapsule;

/// The columns a row is written from, by index in the schema.
struct Columns {
    chrom: usize,
    start: usize,
    end: usize,
    others: Vec<usize>,
    names: Vec<String>,
}

/// An Arrow stream being read.
pub(crate) struct Frame {
    reader: RecordBatchStreamReader,
    columns: Columns,
    batch: Option<RecordBatch>,
    /// The next row of `batch` to write.
    row: usize,
    /// The rows of the batches before `batch`.
    before: u64,
}

impl Frame {
    /// Take the stream of `obj` (its `__arrow_c_stream__()`), and check its columns: a missing
    /// `chrom`, `start` or `end` is a `UsageError`; a coordinate that is not an integer column,
    /// or a column gloci cannot write as text, a `DataError`.
    pub(crate) fn open(py: Python<'_>, obj: &Bound<'_, PyAny>, display: &str) -> PyResult<Frame> {
        let capsule = obj.call_method0("__arrow_c_stream__")?;
        let capsule = capsule.cast::<PyCapsule>().map_err(|_| {
            PyTypeError::new_err(format!(
                "{display}.__arrow_c_stream__() did not return a PyCapsule"
            ))
        })?;
        let ptr = capsule.pointer_checked(Some(c"arrow_array_stream"))?;
        // SAFETY: the capsule is named "arrow_array_stream", so by the Arrow PyCapsule
        // interface it holds an `ArrowArrayStream`, which lives as long as the capsule (held
        // above). `from_raw` moves the stream out, leaving a released one, so the capsule's
        // destructor does not release it again.
        let reader = unsafe {
            RecordBatchStreamReader::from_raw(ptr.as_ptr().cast::<FFI_ArrowArrayStream>())
        }
        .map_err(|e| {
            errors::data_error(
                py,
                format!("{display}: cannot read its Arrow stream: {e}"),
                display,
            )
        })?;
        let columns = Columns::of(py, &reader.schema(), display)?;
        Ok(Frame {
            reader,
            columns,
            batch: None,
            row: 0,
            before: 0,
        })
    }

    /// Write rows until `buf` holds `limit` bytes or the stream ends; `rows` counts the rows
    /// read.
    pub(crate) fn fill(
        &mut self,
        py: Python<'_>,
        buf: &mut Vec<u8>,
        limit: usize,
        rows: &mut u64,
        display: &str,
    ) -> Fill {
        while buf.len() < limit {
            let batch = match &self.batch {
                Some(b) if self.row < b.num_rows() => b,
                _ => {
                    if let Some(done) = self.batch.take() {
                        self.before += done.num_rows() as u64;
                    }
                    self.row = 0;
                    match self.reader.next() {
                        None => return Fill::End,
                        Some(Ok(b)) => self.batch = Some(b),
                        Some(Err(e)) => {
                            return Fill::Failed(errors::data_error(
                                py,
                                format!("{display}: cannot read its Arrow stream: {e}"),
                                display,
                            ));
                        }
                    }
                    continue;
                }
            };
            match write_rows(&self.columns, batch, self.row, self.before, limit, buf) {
                Ok(next) => {
                    self.row = next;
                    *rows = self.before + next as u64;
                }
                Err((row, msg)) => {
                    *rows = row;
                    return Fill::Failed(errors::data_at(
                        py,
                        format!("{display} row {}: {msg}", gloci_core::thousands(row)),
                        display,
                        row,
                    ));
                }
            }
        }
        Fill::More
    }
}

/// pandas writes a DataFrame's index, unless it is a plain range, as `__index_level_N__`.
fn is_index_column(name: &str) -> bool {
    name.strip_prefix("__index_level_")
        .and_then(|r| r.strip_suffix("__"))
        .is_some_and(|n| !n.is_empty() && n.bytes().all(|c| c.is_ascii_digit()))
}

fn is_coordinate_type(t: &DataType) -> bool {
    t.is_integer()
}

impl Columns {
    fn of(py: Python<'_>, schema: &Schema, display: &str) -> PyResult<Columns> {
        let names: Vec<String> = schema.fields().iter().map(|f| f.name().clone()).collect();
        let find = |n: &str| names.iter().position(|x| x == n);
        let missing: Vec<&str> = ["chrom", "start", "end"]
            .into_iter()
            .filter(|n| find(n).is_none())
            .collect();
        if !missing.is_empty() {
            let has = if names.is_empty() {
                "it has none".to_string()
            } else {
                format!("it has {}", names.join(", "))
            };
            return Err(errors::usage(
                py,
                format!(
                    "{display} has no {} column; an Arrow source needs the columns chrom, start \
and end ({has})",
                    missing.join(" or ")
                ),
            ));
        }
        let (chrom, start, end) = (
            find("chrom").expect("checked"),
            find("start").expect("checked"),
            find("end").expect("checked"),
        );
        for at in [start, end] {
            let t = schema.field(at).data_type();
            if !is_coordinate_type(t) {
                let kind = if t.is_floating() {
                    format!("floating-point numbers ({t})")
                } else {
                    t.to_string()
                };
                return Err(errors::data_error(
                    py,
                    format!(
                        "{display} column {} holds {kind}; start and end must be integer \
columns (int32, int64, uint32 or uint64): cast the column",
                        names[at]
                    ),
                    display,
                ));
            }
        }
        let others: Vec<usize> = (0..names.len())
            .filter(|&i| ![chrom, start, end].contains(&i) && !is_index_column(&names[i]))
            .collect();
        for &at in std::iter::once(&chrom).chain(&others) {
            let t = schema.field(at).data_type();
            if cells(new_empty_array(t).as_ref()).is_none() {
                return Err(errors::data_error(
                    py,
                    format!(
                        "{display} column {} holds {t}, which gloci cannot write as text; cast \
it to strings or numbers",
                        names[at]
                    ),
                    display,
                ));
            }
        }
        Ok(Columns {
            chrom,
            start,
            end,
            others,
            names,
        })
    }
}

/// Write the rows of `batch` from `from` until `buf` holds `limit` bytes; the next row, or the
/// 1-based row (counting `before` rows of earlier batches) and what is wrong with it.
fn write_rows(
    cols: &Columns,
    batch: &RecordBatch,
    from: usize,
    before: u64,
    limit: usize,
    buf: &mut Vec<u8>,
) -> Result<usize, (u64, String)> {
    let chrom = cells(batch.column(cols.chrom).as_ref()).expect("checked when opened");
    let start = coords(batch.column(cols.start).as_ref()).expect("checked when opened");
    let end = coords(batch.column(cols.end).as_ref()).expect("checked when opened");
    let others: Vec<(usize, Box<dyn Cells + '_>)> = cols
        .others
        .iter()
        .map(|&i| {
            (
                i,
                cells(batch.column(i).as_ref()).expect("checked when opened"),
            )
        })
        .collect();
    let n = batch.num_rows();
    let mut r = from;
    while r < n && buf.len() < limit {
        let row = before + r as u64 + 1;
        let mark = buf.len();
        let fail = |buf: &mut Vec<u8>, msg: String| {
            buf.truncate(mark);
            Err((row, msg))
        };
        let at = buf.len();
        if !chrom.write(r, buf) {
            return fail(buf, null_message(&cols.names[cols.chrom]));
        }
        if let Some(what) = breaks(&buf[at..]) {
            let msg = split_message(&cols.names[cols.chrom], &buf[at..], what);
            return fail(buf, msg);
        }
        for (at, c) in [(cols.start, &start), (cols.end, &end)] {
            buf.push(b'\t');
            match c.get(r) {
                None => return fail(buf, null_message(&cols.names[at])),
                Some(Err(v)) => {
                    return fail(
                        buf,
                        format!(
                            "column {} is {v}; coordinates are 0 or more",
                            cols.names[at]
                        ),
                    );
                }
                Some(Ok(v)) => buf.extend_from_slice(itoa::Buffer::new().format(v).as_bytes()),
            }
        }
        for (i, c) in &others {
            buf.push(b'\t');
            let at = buf.len();
            if !c.write(r, buf) {
                buf.push(b'.');
            } else if let Some(what) = breaks(&buf[at..]) {
                let msg = split_message(&cols.names[*i], &buf[at..], what);
                return fail(buf, msg);
            }
        }
        buf.push(b'\n');
        r += 1;
    }
    Ok(r)
}

fn null_message(column: &str) -> String {
    format!(
        "column {column} is null; chrom, start and end must be set in every row: drop or fill it"
    )
}

fn split_message(column: &str, value: &[u8], what: &str) -> String {
    format!(
        "column {column} ({:?}) holds {what}, which would split the line; replace or remove it",
        String::from_utf8_lossy(value)
    )
}

/// A column's values as text.
trait Cells {
    /// Write row `i`'s value; false (writing nothing) when it is null.
    fn write(&self, i: usize, buf: &mut Vec<u8>) -> bool;
}

struct Ints<'a, T: ArrowPrimitiveType>(&'a arrow_array::PrimitiveArray<T>);

impl<T: ArrowPrimitiveType> Cells for Ints<'_, T>
where
    T::Native: itoa::Integer,
{
    fn write(&self, i: usize, buf: &mut Vec<u8>) -> bool {
        if self.0.is_null(i) {
            return false;
        }
        buf.extend_from_slice(itoa::Buffer::new().format(self.0.value(i)).as_bytes());
        true
    }
}

struct Floats<'a, T: ArrowPrimitiveType>(&'a arrow_array::PrimitiveArray<T>, fn(T::Native) -> f64);

impl<T: ArrowPrimitiveType> Cells for Floats<'_, T> {
    fn write(&self, i: usize, buf: &mut Vec<u8>) -> bool {
        if self.0.is_null(i) {
            return false;
        }
        push_float((self.1)(self.0.value(i)), buf);
        true
    }
}

/// Strings and binaries, as their bytes.
struct Bytes<'a, A: ?Sized>(&'a A, fn(&A, usize) -> &[u8]);

impl<A: Array + ?Sized> Cells for Bytes<'_, A> {
    fn write(&self, i: usize, buf: &mut Vec<u8>) -> bool {
        if self.0.is_null(i) {
            return false;
        }
        buf.extend_from_slice((self.1)(self.0, i));
        true
    }
}

struct Bools<'a>(&'a arrow_array::BooleanArray);

impl Cells for Bools<'_> {
    fn write(&self, i: usize, buf: &mut Vec<u8>) -> bool {
        if self.0.is_null(i) {
            return false;
        }
        buf.extend_from_slice(if self.0.value(i) { b"True" } else { b"False" });
        true
    }
}

struct Nulls;

impl Cells for Nulls {
    fn write(&self, _: usize, _: &mut Vec<u8>) -> bool {
        false
    }
}

/// A dictionary-encoded column (pandas' categoricals, polars' Categorical): its values, by key.
struct Dict<'a, K: ArrowDictionaryKeyType> {
    keys: &'a arrow_array::PrimitiveArray<K>,
    len: usize,
    values: Box<dyn Cells + 'a>,
}

impl<K: ArrowDictionaryKeyType> Cells for Dict<'_, K>
where
    K::Native: TryInto<usize>,
{
    fn write(&self, i: usize, buf: &mut Vec<u8>) -> bool {
        if self.keys.is_null(i) {
            return false;
        }
        match self.keys.value(i).try_into() {
            Ok(k) if k < self.len => self.values.write(k, buf),
            _ => false,
        }
    }
}

fn dict<'a, K: ArrowDictionaryKeyType>(a: &'a dyn Array) -> Option<Box<dyn Cells + 'a>>
where
    K::Native: TryInto<usize>,
{
    let d = a.as_dictionary_opt::<K>()?;
    Some(Box::new(Dict {
        keys: d.keys(),
        len: d.values().len(),
        values: cells(d.values().as_ref())?,
    }))
}

/// The writer of a column's values, or None for a type gloci does not write as text.
fn cells(a: &dyn Array) -> Option<Box<dyn Cells + '_>> {
    Some(match a.data_type() {
        DataType::Null => Box::new(Nulls),
        DataType::Boolean => Box::new(Bools(a.as_boolean_opt()?)),
        DataType::Int8 => Box::new(Ints(a.as_primitive_opt::<Int8Type>()?)),
        DataType::Int16 => Box::new(Ints(a.as_primitive_opt::<Int16Type>()?)),
        DataType::Int32 => Box::new(Ints(a.as_primitive_opt::<Int32Type>()?)),
        DataType::Int64 => Box::new(Ints(a.as_primitive_opt::<Int64Type>()?)),
        DataType::UInt8 => Box::new(Ints(a.as_primitive_opt::<UInt8Type>()?)),
        DataType::UInt16 => Box::new(Ints(a.as_primitive_opt::<UInt16Type>()?)),
        DataType::UInt32 => Box::new(Ints(a.as_primitive_opt::<UInt32Type>()?)),
        DataType::UInt64 => Box::new(Ints(a.as_primitive_opt::<UInt64Type>()?)),
        DataType::Float16 => Box::new(Floats(a.as_primitive_opt::<Float16Type>()?, |v| v.to_f64())),
        DataType::Float32 => Box::new(Floats(a.as_primitive_opt::<Float32Type>()?, f64::from)),
        DataType::Float64 => Box::new(Floats(a.as_primitive_opt::<Float64Type>()?, |v| v)),
        DataType::Utf8 => Box::new(Bytes(a.as_string_opt::<i32>()?, |s, i| {
            s.value(i).as_bytes()
        })),
        DataType::LargeUtf8 => Box::new(Bytes(a.as_string_opt::<i64>()?, |s, i| {
            s.value(i).as_bytes()
        })),
        DataType::Utf8View => {
            Box::new(Bytes(a.as_string_view_opt()?, |s, i| s.value(i).as_bytes()))
        }
        DataType::Binary => Box::new(Bytes(a.as_binary_opt::<i32>()?, |s, i| s.value(i))),
        DataType::LargeBinary => Box::new(Bytes(a.as_binary_opt::<i64>()?, |s, i| s.value(i))),
        DataType::BinaryView => Box::new(Bytes(a.as_binary_view_opt()?, |s, i| s.value(i))),
        DataType::Dictionary(k, _) => match k.as_ref() {
            DataType::Int8 => dict::<Int8Type>(a)?,
            DataType::Int16 => dict::<Int16Type>(a)?,
            DataType::Int32 => dict::<Int32Type>(a)?,
            DataType::Int64 => dict::<Int64Type>(a)?,
            DataType::UInt8 => dict::<UInt8Type>(a)?,
            DataType::UInt16 => dict::<UInt16Type>(a)?,
            DataType::UInt32 => dict::<UInt32Type>(a)?,
            DataType::UInt64 => dict::<UInt64Type>(a)?,
            _ => return None,
        },
        _ => return None,
    })
}

/// A coordinate column's values: a value, a negative one (`Err`), or None for a null.
trait Coords {
    fn get(&self, i: usize) -> Option<Result<u64, i64>>;
}

struct Signed<'a, T: ArrowPrimitiveType>(&'a arrow_array::PrimitiveArray<T>);

impl<T: ArrowPrimitiveType> Coords for Signed<'_, T>
where
    T::Native: Into<i64>,
{
    fn get(&self, i: usize) -> Option<Result<u64, i64>> {
        if self.0.is_null(i) {
            return None;
        }
        let v: i64 = self.0.value(i).into();
        Some(u64::try_from(v).map_err(|_| v))
    }
}

struct Unsigned<'a, T: ArrowPrimitiveType>(&'a arrow_array::PrimitiveArray<T>);

impl<T: ArrowPrimitiveType> Coords for Unsigned<'_, T>
where
    T::Native: Into<u64>,
{
    fn get(&self, i: usize) -> Option<Result<u64, i64>> {
        if self.0.is_null(i) {
            return None;
        }
        Some(Ok(self.0.value(i).into()))
    }
}

fn coords(a: &dyn Array) -> Option<Box<dyn Coords + '_>> {
    Some(match a.data_type() {
        DataType::Int8 => Box::new(Signed(a.as_primitive_opt::<Int8Type>()?)),
        DataType::Int16 => Box::new(Signed(a.as_primitive_opt::<Int16Type>()?)),
        DataType::Int32 => Box::new(Signed(a.as_primitive_opt::<Int32Type>()?)),
        DataType::Int64 => Box::new(Signed(a.as_primitive_opt::<Int64Type>()?)),
        DataType::UInt8 => Box::new(Unsigned(a.as_primitive_opt::<UInt8Type>()?)),
        DataType::UInt16 => Box::new(Unsigned(a.as_primitive_opt::<UInt16Type>()?)),
        DataType::UInt32 => Box::new(Unsigned(a.as_primitive_opt::<UInt32Type>()?)),
        DataType::UInt64 => Box::new(Unsigned(a.as_primitive_opt::<UInt64Type>()?)),
        _ => return None,
    })
}
