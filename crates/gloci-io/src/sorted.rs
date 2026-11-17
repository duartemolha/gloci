//! On-the-fly sortedness check shared by every streaming command (spec §4.5).

use crate::align::SortOrder;
use crate::detect::Format;
use crate::error::{IoError, Loc, OrderError};
use crate::record::RecordBatch;
use gloci_core::{ChromId, ChromTable, OrderChecker, OrderSource, OrderViolation};

/// Longest record excerpt shown in an error, in characters.
const SHOW_CHARS: usize = 120;
/// Bytes kept per record: enough for `SHOW_CHARS` characters of UTF-8.
const KEEP_BYTES: usize = SHOW_CHARS * 4;

/// `name` as a command line names it: standard input, which a reader calls
/// `stdin`, is `-`.
pub fn as_argument(name: &str) -> &str {
    if name == "stdin" { "-" } else { name }
}

/// The `hint:` line for a sortedness error.
pub fn sort_hint(name: &str, source: &OrderSource) -> String {
    let name = as_argument(name);
    match source {
        OrderSource::GenomeFile(p) => format!("run: gloci sort {name} -g {}", p.display()),
        OrderSource::BuiltinGenome(id) => format!("run: gloci sort {name} -g {id}"),
        OrderSource::AlignmentHeader(p) => format!("run: gloci sort {name} -g {p}"),
        OrderSource::Karyotype => format!("run: gloci sort {name} --chrom-order karyotype"),
        _ => format!("run: gloci sort {name}"),
    }
}

/// The `hint:` line for an alignment input that is out of order: sort it, or let
/// the command sort it (D103). `so` is what its header declares.
pub fn alignment_sort_hint(name: &str, source: &OrderSource, so: SortOrder) -> String {
    let name = as_argument(name);
    let out = match name.strip_suffix(".bam") {
        Some(stem) => format!("{stem}.sorted.bam"),
        None => "sorted.bam".to_string(),
    };
    let genome = match source {
        OrderSource::GenomeFile(p) => format!(" -g {}", p.display()),
        OrderSource::BuiltinGenome(id) => format!(" -g {id}"),
        _ => String::new(),
    };
    let declared = if so == SortOrder::Queryname {
        "; the header says SO:queryname, so it is not coordinate-sorted"
    } else {
        ""
    };
    format!("run: gloci sort {name} -o {out}{genome}, or pass --sort{declared}")
}

/// Added to the hint when a karyotype-order check fails on a pair of
/// chromosomes that are in lexicographic order: the file was probably
/// sorted by sort -k1,1, sort-bed or bedtools sort.
pub const LEX_SORTED_NOTE: &str = " (this file looks sorted by sort -k1,1, sort-bed or bedtools sort; drop --chrom-order karyotype, or re-sort it with gloci sort --chrom-order karyotype; see docs/migration.md)";

/// Added to the hint when a lexicographic-order check fails on a pair of
/// chromosomes that are in karyotype order: the file was probably written by
/// `gloci convert` from a VCF or BAM, or sorted with `--chrom-order karyotype`.
/// This is the note of a command that has `-g` and `--chrom-order`; see
/// [`OrderOptions`] for the others.
pub const KARYOTYPE_SORTED_NOTE: &str = " (this file looks sorted in karyotype order, as a VCF or BAM header lists it; pass -g GENOME or --chrom-order karyotype, or re-sort it with gloci sort)";

/// Which of the options that set the chromosome order the command has, so that
/// the note of an order error suggests only those.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct OrderOptions {
    /// `-g GENOME`.
    pub genome: bool,
    /// `--chrom-order`.
    pub chrom_order: bool,
}

impl OrderOptions {
    /// `-g` and `--chrom-order`, which most commands that check the order have.
    pub const BOTH: OrderOptions = OrderOptions {
        genome: true,
        chrom_order: true,
    };
    /// `--chrom-order` without `-g` (`vcf concat`, `vcf merge`, `vcf annotate`).
    pub const CHROM_ORDER: OrderOptions = OrderOptions {
        genome: false,
        chrom_order: true,
    };

    /// The note for a file in karyotype order read in lexicographic order.
    fn karyotype_note(self) -> &'static str {
        match (self.genome, self.chrom_order) {
            (true, true) => KARYOTYPE_SORTED_NOTE,
            (false, true) => {
                " (this file looks sorted in karyotype order, as a VCF or BAM header lists it; pass --chrom-order karyotype, or re-sort it with gloci sort)"
            }
            (true, false) => {
                " (this file looks sorted in karyotype order, as a VCF or BAM header lists it; pass -g GENOME, or re-sort it with gloci sort)"
            }
            (false, false) => {
                " (this file looks sorted in karyotype order, as a VCF or BAM header lists it; re-sort it with gloci sort)"
            }
        }
    }
}

/// The note for an out-of-order pair of chromosomes (`prev` came first, `next`
/// sorts before it), when the pair is in the other mode's order. `options` are
/// the order options of the command, which the note may suggest.
pub fn order_note(
    source: &OrderSource,
    prev: &[u8],
    next: &[u8],
    options: OrderOptions,
) -> Option<&'static str> {
    match source {
        // Karyotype order comes from --chrom-order karyotype, which the
        // command therefore has.
        OrderSource::Karyotype if prev < next => Some(LEX_SORTED_NOTE),
        OrderSource::Lexicographic
            if gloci_core::karyotype_cmp(prev, next) == std::cmp::Ordering::Less =>
        {
            Some(options.karyotype_note())
        }
        _ => None,
    }
}

/// Checks that one input's records arrive in `(chrom, key)` order.
pub struct OrderGuard {
    name: String,
    binary: bool,
    checker: OrderChecker,
    /// For an alignment input, what its header declares (it picks the hint).
    alignment: Option<SortOrder>,
    /// The last record seen: `(loc, shown bytes)`; the bytes live in `prev`.
    prev_loc: u64,
    prev: Vec<u8>,
    /// Contiguity mode: the chromosomes in the order this input shows them.
    seen: Option<Seen>,
    /// The order options of the command, for the note of an error.
    options: OrderOptions,
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
            binary: format.is_binary(),
            alignment: format.is_alignment().then_some(SortOrder::Unknown),
            checker: OrderChecker::new(),
            prev_loc: 0,
            prev: Vec::new(),
            seen: None,
            options: OrderOptions::BOTH,
        }
    }

    /// Name only `options` in the note of an order error: the order options
    /// the command has (both, unless this is called).
    pub fn with_options(mut self, options: OrderOptions) -> OrderGuard {
        self.options = options;
        self
    }

    /// Note what the alignment header of this input declares, for the hint.
    pub fn with_sort_order(mut self, so: SortOrder) -> OrderGuard {
        if self.alignment.is_some() {
            self.alignment = Some(so);
        }
        self
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
        self.check_range(chroms, batch, 0..batch.recs.len())?;
        self.settle(chroms, batch, batch.recs.len());
        Ok(())
    }

    /// Check `batch.recs[range]` only, in order after whatever was checked
    /// before: for a consumer that validates each record as it takes it, not
    /// a whole read-ahead batch at once. Records of a batch must be checked
    /// in order, and [`OrderGuard::settle`] called before the batch is
    /// replaced.
    pub fn check_range(
        &mut self,
        chroms: &ChromTable,
        batch: &RecordBatch,
        range: std::ops::Range<usize>,
    ) -> Result<(), IoError> {
        for i in range {
            let m = &batch.recs[i];
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
                        let mut hint = match self.alignment {
                            Some(so) => alignment_sort_hint(&self.name, chroms.order_source(), so),
                            None => sort_hint(&self.name, chroms.order_source()),
                        };
                        if let OrderViolation::ChromBackwards { prev, next } = v
                            && self.seen.is_none()
                            && let Some(note) = order_note(
                                chroms.order_source(),
                                chroms.name(prev),
                                chroms.name(next),
                                self.options,
                            )
                        {
                            hint.push_str(note);
                        }
                        hint
                    },
                }
                .into());
            }
        }
        Ok(())
    }

    /// Note the last of the first `checked` records of `batch`, which is
    /// about to be replaced, as the one the next batch's first record
    /// follows in an error.
    pub fn settle(&mut self, chroms: &ChromTable, batch: &RecordBatch, checked: usize) {
        if let Some(last) = checked.checked_sub(1) {
            let mut prev = std::mem::take(&mut self.prev);
            prev.clear();
            self.render(chroms, batch, last, &mut prev);
            self.prev = prev;
            self.prev_loc = batch.recs[last].loc;
        }
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
