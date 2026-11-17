//! The bigBed column layout and its autoSql (D171).
//!
//! The first record decides the layout: its standard BED columns are the longest prefix of
//! size 3, 4, 5, 6, 8, 9 or 12 whose fields are valid, and any further columns are extra.
//! Every later record must have the same column count and valid standard columns. An autoSql
//! given with `--autosql`, or a bigBed input's definedFieldCount, says instead how many of the
//! columns are standard (D171).

use crate::bed12_blocks;
use bigtools::bed::autosql::parse::{Field, FieldType, ParseError};

/// The sizes a BED prefix can have.
const STANDARD_SIZES: [usize; 7] = [3, 4, 5, 6, 8, 9, 12];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BedLayout {
    /// The column count every record has.
    pub columns: usize,
    /// How many leading columns are standard BED columns.
    pub standard: usize,
}

fn u32_field(name: &str, v: &[u8]) -> Result<u32, String> {
    std::str::from_utf8(v)
        .ok()
        .filter(|s| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit()))
        .and_then(|s| s.parse::<u32>().ok())
        .ok_or_else(|| {
            format!(
                "{name} '{}' is not an unsigned 32-bit integer",
                String::from_utf8_lossy(v)
            )
        })
}

fn check_rgb(v: &[u8]) -> Result<(), String> {
    let text = String::from_utf8_lossy(v);
    if v.contains(&b',') {
        let parts: Vec<&str> = text.split(',').collect();
        let ok = parts.len() == 3
            && parts.iter().all(|p| {
                !p.is_empty() && p.bytes().all(|b| b.is_ascii_digit()) && p.parse::<u8>().is_ok()
            });
        if ok {
            return Ok(());
        }
        return Err(format!(
            "itemRgb '{text}' is not R,G,B with each value 0 to 255"
        ));
    }
    u32_field("itemRgb", v).map(|_| ()).map_err(|_| {
        format!(
            "itemRgb '{text}' is neither an unsigned integer nor R,G,B with each value 0 to 255"
        )
    })
}

/// Are the standard columns of `f` up to the prefix of size `size` valid?
fn check_prefix(f: &[&[u8]], size: usize) -> Result<(), String> {
    if size >= 4 && f[3].is_empty() {
        return Err("name is empty".into());
    }
    if size >= 5 {
        u32_field("score", f[4])?;
    }
    if size >= 6 && !matches!(f[5], b"+" | b"-" | b".") {
        return Err(format!(
            "strand '{}' is not +, - or .",
            String::from_utf8_lossy(f[5])
        ));
    }
    if size >= 8 {
        u32_field("thickStart", f[6])?;
        u32_field("thickEnd", f[7])?;
    }
    if size >= 9 {
        check_rgb(f[8])?;
    }
    if size >= 12 {
        let head = f[..12].join(&b'\t');
        bed12_blocks(&head)?;
    }
    Ok(())
}

fn fields(line: &[u8]) -> Vec<&[u8]> {
    line.split(|&b| b == b'\t').collect()
}

/// The layout the first record `first` decides.
pub fn bed_layout(first: &[u8]) -> Result<BedLayout, String> {
    let f = fields(first);
    if f.len() < 3 {
        return Err(format!(
            "a BED record needs at least 3 columns, this one has {}",
            f.len()
        ));
    }
    let mut standard = 3;
    for size in STANDARD_SIZES.into_iter().skip(1) {
        if size > f.len() || check_prefix(&f, size).is_err() {
            break;
        }
        standard = size;
    }
    Ok(BedLayout {
        columns: f.len(),
        standard,
    })
}

/// Does `line` fit `layout`: its column count, and valid standard columns?
pub fn check_bed_line(line: &[u8], layout: &BedLayout) -> Result<(), String> {
    let f = fields(line);
    if f.len() != layout.columns {
        return Err(format!(
            "the record has {} columns, but the first record has {} columns; every record of a bigBed \
needs the same columns",
            f.len(),
            layout.columns
        ));
    }
    check_prefix(&f, layout.standard)
}

/// The autoSql of `layout`: the standard columns by their BED names, the others `field1`, ...
pub fn bed_autosql(layout: &BedLayout) -> String {
    const STANDARD: [&str; 12] = [
        "string chrom; \"Chromosome\"",
        "uint chromStart; \"Start, 0-based\"",
        "uint chromEnd; \"End, exclusive\"",
        "string name; \"Name\"",
        "uint score; \"Score\"",
        "char[1] strand; \"Strand: +, - or .\"",
        "uint thickStart; \"Start of the thick part\"",
        "uint thickEnd; \"End of the thick part\"",
        "uint itemRgb; \"Colour, R,G,B\"",
        "int blockCount; \"Number of blocks\"",
        "int[blockCount] blockSizes; \"Block sizes\"",
        "int[blockCount] blockStarts; \"Block starts, relative to chromStart\"",
    ];
    let mut out = String::from("table bed\n\"BED records written by gloci\"\n(\n");
    for row in &STANDARD[..layout.standard] {
        out.push_str("    ");
        out.push_str(row);
        out.push('\n');
    }
    for i in 1..=layout.columns - layout.standard {
        out.push_str(&format!("    string field{i}; \"Extra column {i}\"\n"));
    }
    out.push_str(")\n");
    out
}

/// The fields of the table an autoSql text declares, by bigtools' parser.
fn autosql_fields(text: &str) -> Result<Vec<Field>, String> {
    let mut declarations = bigtools::bed::autosql::parse::parse_autosql(text)
        .map_err(|e| format!("the autoSql does not parse: {}", describe(&e)))?;
    declarations
        .pop()
        .map(|d| d.fields)
        .ok_or_else(|| "the autoSql declares no table".to_string())
}

/// The number of fields an autoSql text declares, by bigtools' parser.
pub fn autosql_field_count(text: &str) -> Result<usize, String> {
    autosql_fields(text).map(|f| f.len())
}

/// The largest size a BED prefix can have that is at most `n` (at least 3).
pub fn standard_size(n: usize) -> usize {
    STANDARD_SIZES
        .into_iter()
        .rev()
        .find(|&s| s <= n)
        .unwrap_or(3)
}

/// The kind of value a standard BED field holds, as an autoSql types it.
#[derive(Clone, Copy)]
enum Kind {
    Text,
    Int,
    Char1,
    IntList,
}

/// The standard BED fields after chrom, chromStart and chromEnd, in order: the names each may
/// have, and its kind.
const BED_FIELDS: [(&[&str], Kind); 9] = [
    (&["name"], Kind::Text),
    (&["score"], Kind::Int),
    (&["strand"], Kind::Char1),
    (&["thickStart"], Kind::Int),
    (&["thickEnd"], Kind::Int),
    (&["itemRgb", "reserved"], Kind::Int),
    (&["blockCount"], Kind::Int),
    (&["blockSizes"], Kind::IntList),
    (&["chromStarts", "blockStarts"], Kind::IntList),
];

fn is_bed_field(f: &Field, names: &[&str], kind: Kind) -> bool {
    let int = matches!(
        f.field_type,
        FieldType::Int
            | FieldType::Uint
            | FieldType::Short
            | FieldType::Ushort
            | FieldType::Byte
            | FieldType::Ubyte
    );
    let typed = match kind {
        Kind::Text => {
            matches!(f.field_type, FieldType::String | FieldType::Lstring) && f.field_size.is_none()
        }
        Kind::Int => int && f.field_size.is_none(),
        Kind::Char1 => {
            matches!(f.field_type, FieldType::Char) && f.field_size.as_deref() == Some("1")
        }
        Kind::IntList => int && f.field_size.is_some(),
    };
    typed && names.contains(&f.name.as_str())
}

/// How many of the leading fields of an autoSql are standard BED columns, by name and type in
/// BED's order: the first three always (every bigBed record starts with chrom, chromStart and
/// chromEnd, whatever the autoSql calls them), then as long as the fields match, rounded down
/// to a size BED has.
pub fn autosql_standard_columns(text: &str) -> Result<usize, String> {
    let fields = autosql_fields(text)?;
    let matched = fields
        .iter()
        .skip(3)
        .zip(BED_FIELDS)
        .take_while(|(f, (names, kind))| is_bed_field(f, names, *kind))
        .count();
    Ok(standard_size(3 + matched))
}

/// bigtools' `ParseError` in words (it has no `Display`).
fn describe(e: &ParseError) -> String {
    match e {
        ParseError::InvalidDeclareType(s) => {
            format!("invalid declaration type '{s}' (expected table, simple or object)")
        }
        ParseError::InvalidDeclareName(s) => format!("invalid declaration name '{s}'"),
        ParseError::InvalidDeclareBrackets(s) => {
            format!("the declaration's brackets are wrong, at '{s}'")
        }
        ParseError::InvalidFieldSizeClose(s) => {
            format!("a field's size bracket is not closed, at '{s}'")
        }
        ParseError::InvalidFieldCommentSeparater(s) => {
            format!("a field's comment is not separated by a semicolon, at '{s}'")
        }
        ParseError::InvalidFieldValuesBrackets(s) => {
            format!("a field's values brackets are wrong, at '{s}'")
        }
        ParseError::InvalidIndexSizeBrackets(s) => {
            format!("an index size bracket is wrong, at '{s}'")
        }
    }
}
