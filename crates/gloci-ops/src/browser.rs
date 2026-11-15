//! `gloci igv` and `gloci links`: scripts and pages that point a genome browser at each feature
//! of a BED, GFF or VCF file (the bedtools commands of the same names).
//!
//! Both read the features in input order and write as they go. They share bedtools' reading of a
//! record: BED and GFF need the same number of columns on every record, a zero-length BED record
//! is shown one base wider on each side, and a VCF record is its REF bases.

use gloci_core::ChromTable;
use gloci_io::{Format, IoError, Loc, RecordBatch, RecordReader, field};
use std::io::Write;

fn write_err(source: std::io::Error) -> IoError {
    IoError::Write {
        path: "output".into(),
        source,
    }
}

/// One feature as the browser commands see it.
struct Feature<'a> {
    chrom: &'a [u8],
    start: i128,
    end: i128,
    /// The name the record gives: BED column 4, the GFF feature type, a VCF record's
    /// `REF/ALT[_ID]`; absent from a BED record of three columns.
    name: Option<Vec<u8>>,
    /// The cells `links` writes after the position.
    cells: Vec<Vec<u8>>,
}

/// The formats both commands read.
fn check_format(command: &str, format: Format) -> Result<(), IoError> {
    match format {
        Format::Bed | Format::BedGraph | Format::Gff3 | Format::Gtf | Format::Vcf => Ok(()),
        other => Err(gloci_io::UsageError(format!(
            "{command} reads BED, GFF or VCF; this input is {}",
            other.name()
        ))
        .into()),
    }
}

/// What `links` writes for a record of `columns` columns: the name, then the score, then the
/// strand. bedtools has these for BED records of 4, 5, 6, 9 and 12 columns and for GFF, and for
/// a VCF record of 9 or 12 columns, whose score is empty and strand `+`; no other record gets
/// a cell.
fn cell_count(columns: usize) -> usize {
    match columns {
        4 => 1,
        5 => 2,
        6 | 9 | 12 => 3,
        _ => 0,
    }
}

fn feature<'a>(
    format: Format,
    raw: &'a [u8],
    span: gloci_core::Interval,
    columns: usize,
) -> Result<Feature<'a>, String> {
    let col = |n: usize| field(raw, n).unwrap_or_default();
    let chrom = col(0);
    match format {
        Format::Vcf => {
            let pos: i128 = std::str::from_utf8(col(1))
                .ok()
                .and_then(|p| p.parse().ok())
                .ok_or("the POS is not a whole number")?;
            let mut name = col(3).to_vec();
            name.push(b'/');
            name.extend_from_slice(col(4));
            if col(2) != b"." {
                name.push(b'_');
                name.extend_from_slice(col(2));
            }
            let cells = match columns {
                9 | 12 => vec![name.clone(), Vec::new(), b"+".to_vec()],
                _ => Vec::new(),
            };
            Ok(Feature {
                chrom,
                start: pos - 1,
                end: pos - 1 + col(3).len() as i128,
                name: Some(name),
                cells,
            })
        }
        Format::Gff3 | Format::Gtf => Ok(Feature {
            chrom,
            start: i128::from(span.start),
            end: i128::from(span.end),
            name: Some(col(2).to_vec()),
            cells: vec![col(2).to_vec(), col(5).to_vec(), col(6).to_vec()],
        }),
        _ => {
            let (mut start, mut end) = (i128::from(span.start), i128::from(span.end));
            if start == end {
                start -= 1;
                end += 1;
            }
            let cells = [3usize, 4, 5]
                .iter()
                .take(cell_count(columns))
                .map(|&n| col(n).to_vec())
                .collect();
            Ok(Feature {
                chrom,
                start,
                end,
                name: (columns >= 4).then(|| col(3).to_vec()),
                cells,
            })
        }
    }
}

/// Call `each` with every feature of `input`, in order, and `line` (the physical line of the
/// record) for an error of its own.
fn for_each_feature(
    command: &str,
    mut input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    mut each: impl FnMut(&Feature, u64) -> Result<(), IoError>,
) -> Result<(), IoError> {
    let format = input.format();
    check_format(command, format)?;
    let path = input.name().to_string();
    let mut first_columns: Option<usize> = None;
    let mut batch = RecordBatch::default();
    loop {
        batch.clear();
        let more = input.read_batch(chroms, &mut batch, 4096)?;
        for (i, m) in batch.recs.iter().enumerate() {
            // bedtools ignores one tab at the end of a line.
            let raw = batch.raw(i);
            let raw = raw.strip_suffix(b"\t").unwrap_or(raw);
            let columns = gloci_io::column_count(raw);
            let parse_error = |msg: String| IoError::Parse {
                path: path.clone(),
                loc: Loc::Line(m.loc),
                msg,
            };
            if format != Format::Vcf {
                let first = *first_columns.get_or_insert(columns);
                if columns != first {
                    return Err(parse_error(format!(
                        "{columns} columns, but the first record has {first}; every record \
must have the same number of columns"
                    )));
                }
            }
            let f = feature(format, raw, m.span, columns).map_err(parse_error)?;
            each(&f, m.loc)?;
        }
        if !more {
            return Ok(());
        }
    }
}

/// The sorting `igv` can ask IGV for.
pub const IGV_SORTS: &[&str] = &["base", "position", "strand", "quality", "sample", "readGroup"];

/// What `igv` writes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IgvOptions {
    /// The directory the snapshots go to, written as given.
    pub path: String,
    /// A session file IGV loads first.
    pub session: Option<String>,
    /// One of [`IGV_SORTS`].
    pub sort: Option<String>,
    pub collapse: bool,
    /// Name each image after the record's name.
    pub name: bool,
    /// Bases of each side added to the region shown.
    pub slop: i64,
    /// The image extension: `png`, `eps` or `svg`.
    pub image: String,
}

/// Write the IGV batch script for the features of `input`.
pub fn igv(
    input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    o: &IgvOptions,
    out: &mut dyn Write,
) -> Result<(), IoError> {
    let path = input.name().to_string();
    out.write_all(format!("snapshotDirectory {}\n", o.path).as_bytes())
        .map_err(write_err)?;
    if let Some(s) = &o.session {
        out.write_all(format!("load {s}\n").as_bytes())
            .map_err(write_err)?;
    }
    let mut line = Vec::new();
    for_each_feature("igv", input, chroms, |f, loc| {
        line.clear();
        line.extend_from_slice(b"goto ");
        line.extend_from_slice(f.chrom);
        line.extend_from_slice(
            format!(":{}-{}\n", f.start - i128::from(o.slop), f.end + i128::from(o.slop)).as_bytes(),
        );
        if let Some(s) = &o.sort {
            line.extend_from_slice(format!("sort {s}\n").as_bytes());
        }
        if o.collapse {
            line.extend_from_slice(b"collapse\n");
        }
        line.extend_from_slice(b"snapshot ");
        line.extend_from_slice(f.chrom);
        line.extend_from_slice(format!("_{}_{}", f.start, f.end).as_bytes());
        if o.name {
            match f.name.as_deref().filter(|n| !n.is_empty()) {
                Some(n) => {
                    line.push(b'_');
                    line.extend_from_slice(n);
                }
                None => {
                    return Err(IoError::Parse {
                        path: path.clone(),
                        loc: Loc::Line(loc),
                        msg: "--name names each image after the name column, which is empty here"
                            .into(),
                    });
                }
            }
        }
        if o.slop > 0 {
            line.extend_from_slice(format!("_slop{}", o.slop).as_bytes());
        }
        line.extend_from_slice(format!(".{}\n", o.image).as_bytes());
        out.write_all(&line).map_err(write_err)
    })
}

/// What `links` writes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LinksOptions {
    pub base: String,
    pub org: String,
    pub db: String,
    /// The page's title: the name of the input.
    pub title: String,
}

/// Write the HTML page of links for the features of `input`.
pub fn links(
    input: Box<dyn RecordReader>,
    chroms: &mut ChromTable,
    o: &LinksOptions,
    out: &mut dyn Write,
) -> Result<(), IoError> {
    let head = format!(
        "<html>\n\t<body>\n<title>{}</title>\n\
<br>Firefox users: Press and hold the \"apple\" or \"alt\" key and click link to open in new tab.\n\
<p style=\"font-family:courier\">\n\
<table border=\"0\" align=\"justify\"\n\
<h3>BED Entries from: stdin </h3>\n",
        o.title
    );
    out.write_all(head.as_bytes()).map_err(write_err)?;
    let mut row = Vec::new();
    for_each_feature("links", input, chroms, |f, _| {
        row.clear();
        row.extend_from_slice(b"<tr>\n\t<td>\n\t\t<a href=");
        row.extend_from_slice(
            format!(
                "{}/cgi-bin/hgTracks?org={}&db={}&position=",
                o.base, o.org, o.db
            )
            .as_bytes(),
        );
        row.extend_from_slice(f.chrom);
        row.extend_from_slice(format!(":{}-{}>", f.start, f.end).as_bytes());
        row.extend_from_slice(f.chrom);
        row.extend_from_slice(format!(":{}-{}</a>\n\t</td>\n", f.start + 1, f.end).as_bytes());
        for cell in &f.cells {
            row.extend_from_slice(b"\t<td>\n");
            row.extend_from_slice(cell);
            row.extend_from_slice(b"\n\t</td>\n");
        }
        row.extend_from_slice(b"</tr>\n");
        out.write_all(&row).map_err(write_err)
    })?;
    out.write_all(b"</table>\n</p>\n\t</body>\n</html>\n")
        .map_err(write_err)
}
