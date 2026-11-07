//! On-the-fly sortedness check shared by every streaming command (spec §4.5).

use crate::detect::Format;
use crate::error::{IoError, Loc, OrderError};
use crate::record::RecordBatch;
use gloci_core::{ChromId, ChromTable, OrderChecker, OrderSource, OrderViolation};

/// Longest record excerpt shown in an error, in characters.
const SHOW_CHARS: usize = 120;
/// Bytes kept per record: enough for `SHOW_CHARS` characters of UTF-8.
const KEEP_BYTES: usize = SHOW_CHARS * 4;

/// The `hint:` line for a sortedness error.
pub fn sort_hint(name: &str, source: &OrderSource) -> String {
    match source {
        OrderSource::GenomeFile(p) => format!("run: gloci sort {name} -g {}", p.display()),
        OrderSource::Lexicographic => format!("run: gloci sort {name} --chrom-order lex"),
        _ => format!("run: gloci sort {name}"),
    }
}

/// Added to the hint when a karyotype-order check fails on a pair of
/// chromosomes that are in lexicographic order: the file was probably
/// sorted by an old tool.
const LEX_SORTED_NOTE: &str = " (this file looks sorted by sort -k1,1, sort-bed or bedtools sort; re-sort it with gloci sort, or pass --chrom-order lex; see docs/migration.md)";

/// Checks that one input's records arrive in `(chrom, key)` order.
pub struct OrderGuard {
    name: String,
    binary: bool,
    checker: OrderChecker,
    /// The last record seen: `(loc, shown bytes)`; the bytes live in `prev`.
    prev_loc: u64,
    prev: Vec<u8>,
    /// Contiguity mode: the chromosomes in the order this input shows them.
    seen: Option<Seen>,
}

/// Chromosomes numbered by first appearance in one input.
struct Seen {
    table: ChromTable,
    /// The last `(id in the caller's table, id here)`.
    last: Option<(ChromId, ChromId)>,
}

impl Seen {
    fn local(&mut self, chroms: &ChromTable, id: ChromId) -> ChromId {
        if let Some((outer, local)) = self.last
            && outer == id
        {
            return local;
        }
        let local = self
            .table
            .intern(chroms.name(id))
            .expect("a first-appearance table is open");
        self.last = Some((id, local));
        local
    }
}

impl OrderGuard {
    /// Records must follow the order of the table passed to `check`.
    pub fn new(name: &str, format: Format) -> OrderGuard {
        OrderGuard {
            name: name.to_string(),
            binary: format == Format::Bcf,
            checker: OrderChecker::new(),
            prev_loc: 0,
            prev: Vec::new(),
            seen: None,
        }
    }

    /// Records of a chromosome must be contiguous and sorted by key, in any
    /// chromosome order (what an index, or a scan for regions, needs).
    pub fn contiguous(name: &str, format: Format) -> OrderGuard {
        OrderGuard {
            seen: Some(Seen {
                table: ChromTable::first_appearance(),
                last: None,
            }),
            ..OrderGuard::new(name, format)
        }
    }

    /// Check every record of `batch`, continuing from the previous batch.
    pub fn check(&mut self, chroms: &ChromTable, batch: &RecordBatch) -> Result<(), IoError> {
        for (i, m) in batch.recs.iter().enumerate() {
            let chrom = match self.seen.as_mut() {
                Some(s) => s.local(chroms, m.chrom),
                None => m.chrom,
            };
            let order = self.seen.as_ref().map_or(chroms, |s| &s.table);
            if let Err(v) = self.checker.check(order, chrom, m.key) {
                let (prev_loc, prev) = if i == 0 {
                    (self.prev_loc, show(&self.prev))
                } else {
                    let p = &batch.recs[i - 1];
                    let mut tmp = Vec::new();
                    self.render(chroms, batch, i - 1, &mut tmp);
                    (p.loc, show(&tmp))
                };
                let mut cur = Vec::new();
                self.render(chroms, batch, i, &mut cur);
                let loc = |n| {
                    if self.binary {
                        Loc::Record(n)
                    } else {
                        Loc::Line(n)
                    }
                };
                return Err(OrderError {
                    path: self.name.clone(),
                    prev_loc: loc(prev_loc),
                    prev,
                    loc: loc(m.loc),
                    rec: show(&cur),
                    hint: {
                        let mut hint = sort_hint(&self.name, chroms.order_source());
                        if let OrderViolation::ChromBackwards { prev, next } = v
                            && self.seen.is_none()
                            && *chroms.order_source() == OrderSource::Karyotype
                            && chroms.name(prev) < chroms.name(next)
                        {
                            hint.push_str(LEX_SORTED_NOTE);
                        }
                        hint
                    },
                }
                .into());
            }
        }
        if let Some(last) = batch.recs.len().checked_sub(1) {
            let mut prev = std::mem::take(&mut self.prev);
            prev.clear();
            self.render(chroms, batch, last, &mut prev);
            self.prev = prev;
            self.prev_loc = batch.recs[last].loc;
        }
        Ok(())
    }

    /// Append record `i` as shown in an error: the (cut) raw line for text
    /// formats, `CHROM POS` for BCF.
    fn render(&self, chroms: &ChromTable, batch: &RecordBatch, i: usize, out: &mut Vec<u8>) {
        let m = &batch.recs[i];
        if self.binary {
            out.extend_from_slice(chroms.name(m.chrom));
            out.push(b' ');
            out.extend_from_slice((m.key + 1).to_string().as_bytes());
        } else {
            let raw = batch.raw(i);
            out.extend_from_slice(&raw[..raw.len().min(KEEP_BYTES)]);
        }
    }
}

/// Lossy UTF-8, cut to `SHOW_CHARS` characters.
fn show(bytes: &[u8]) -> String {
    let s = String::from_utf8_lossy(bytes);
    match s.char_indices().nth(SHOW_CHARS) {
        Some((cut, _)) => format!("{}...", &s[..cut]),
        None => s.into_owned(),
    }
}
