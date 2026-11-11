//! `gloci vcf annotate`: the record and header edits of `bcftools annotate`.
//!
//! This module holds everything that does not read files: parsing the
//! options against the headers, the per-record transfers from a VCF record
//! or from tab/BED rows, removals, `--set-id` and `--mark-sites`. The
//! streaming merge of the files is in gloci-ops. Each rule below is checked
//! against bcftools 1.21 (vcfannotate.c, htslib's synced reader).

mod hdr;
mod pair;
mod rec;
mod spec;
mod tab;
mod vcmp;
mod xfer;

pub use hdr::{HeaderEdit, InfoDef, Number, Ty};
pub use pair::{PairRec, pair_rec, pair_rows};
pub use rec::{Rec, format_retain, info_get, info_remove, info_set};
pub use spec::{ColItem, ColKind, ColumnSpec, IdFormat, IdTok, parse_columns, parse_id_format};
pub use tab::{TabLayout, TabRow};

use xfer::Col;

/// Which kind of file `--annotations` names.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SourceKind<'a> {
    /// No `--annotations`.
    None,
    /// A VCF/BCF; its raw header text.
    Vcf(&'a [u8]),
    /// A tab-delimited file; `bed` for a `.bed` name (0-based starts).
    Tab { bed: bool },
}

/// `--mark-sites [+-]TAG`.
#[derive(Clone, Debug)]
pub struct Mark {
    pub tag: Vec<u8>,
    /// `+`/no sign: flag the records listed in the source; `-`: those not.
    pub listed: bool,
}

/// The options of `vcf annotate` that edit records.
#[derive(Clone, Debug, Default)]
pub struct AnnotateOptions {
    pub columns: Option<String>,
    pub remove: Option<String>,
    /// The lines of `--header-lines`.
    pub header_lines: Vec<Vec<u8>>,
    pub set_id: Option<String>,
    pub mark_sites: Option<String>,
    pub min_overlap: Option<String>,
}

/// One removal, run on every record before any transfer.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RmOp {
    Id,
    Qual,
    FilterAll,
    FilterTag(Vec<u8>),
    InfoAll,
    InfoTag(Vec<u8>),
    FormatAll,
    FormatTag(Vec<u8>),
}

/// `--min-overlap`: the fractions of an interval and of the record that must overlap.
#[derive(Clone, Copy, Debug, Default)]
pub struct MinOverlap {
    pub ann: f32,
    pub vcf: f32,
}

/// The compiled options: columns resolved against both headers, the output
/// header built.
pub struct Annotator {
    cols: Vec<Col>,
    rm: Vec<RmOp>,
    keep_pass: bool,
    set_id: Option<IdFormat>,
    mark: Option<Mark>,
    layout: Option<TabLayout>,
    min_overlap: MinOverlap,
    has_source: bool,
    out_header: HeaderEdit,
    vcmp: vcmp::Vcmp,
    scratch: Vec<u8>,
}

fn keyword_match(s: &str, kw: &str) -> bool {
    s.eq_ignore_ascii_case(kw)
}

/// Strip a case-insensitive prefix.
fn strip_ci<'a>(s: &'a str, prefix: &str) -> Option<&'a str> {
    (s.len() >= prefix.len() && s[..prefix.len()].eq_ignore_ascii_case(prefix))
        .then(|| &s[prefix.len()..])
}

/// `-x LIST` compiled: the record operations and the header lines they drop.
fn compile_remove(
    list: &str,
    hdr_in: &HeaderEdit,
    hdr_out: &mut HeaderEdit,
    needs_info: bool,
    set_id_fmt: &str,
) -> Result<(Vec<RmOp>, bool), String> {
    let mut ops: Vec<RmOp> = Vec::new();
    let mut keep_pass = false;
    // keep lists: (kind, id)
    let mut keep: Vec<(&'static str, Vec<u8>)> = Vec::new();
    let (mut keep_info, mut keep_fmt, mut keep_flt) = (false, false, false);
    let clash = |what: &str| {
        format!(
            "--remove {what} runs before --set-id {set_id_fmt}, which reads INFO tags; run two commands (annotate --set-id first, then --remove)"
        )
    };
    for raw in list.split(',') {
        let mut s = raw;
        let mut class: Option<&'static str> = None;
        let mut is_keep = false;
        for (prefix, kind, k) in [
            ("INFO/", "INFO", false),
            ("INF/", "INFO", false),
            ("FORMAT/", "FORMAT", false),
            ("FMT/", "FORMAT", false),
            ("FILTER/", "FILTER", false),
            ("^INFO/", "INFO", true),
            ("^INF/", "INFO", true),
            ("^FORMAT/", "FORMAT", true),
            ("^FMT/", "FORMAT", true),
            ("^FILTER/", "FILTER", true),
        ] {
            if let Some(rest) = strip_ci(s, prefix) {
                s = rest;
                class = Some(kind);
                is_keep = k;
                break;
            }
        }
        if is_keep {
            match class {
                Some("INFO") => keep_info = true,
                Some("FORMAT") => keep_fmt = true,
                _ => keep_flt = true,
            }
        }
        // Once a ^ item of a kind has been seen, the plain items after it
        // of that kind are keep items too (bcftools).
        let is_keep = is_keep
            || match class {
                Some("INFO") => keep_info,
                Some("FORMAT") => keep_fmt,
                Some("FILTER") => keep_flt,
                _ => false,
            };
        let id = s.as_bytes();
        match class {
            Some("FILTER") if !is_keep => {
                keep_pass = true;
                if !hdr_in.has("FILTER", id) {
                    return Err(format!(
                        "--remove FILTER/{s}: the filter `{s}` is not declared in the header (bcftools stops here too); declare it with --header-lines or remove the whole column with -x FILTER"
                    ));
                }
                hdr_out.remove("FILTER", id);
                ops.push(RmOp::FilterTag(id.to_vec()));
            }
            Some("FILTER") => keep.push(("FILTER", id.to_vec())),
            Some(kind) => {
                let declared = hdr_in.has(kind, id);
                if declared && is_keep {
                    keep.push((if kind == "INFO" { "INFO" } else { "FORMAT" }, id.to_vec()));
                } else {
                    // An undeclared tag is still removed from the records
                    // (bcftools warns), but no header line goes.
                    if !declared {
                        log::warn!("vcf annotate: the tag {kind}/{s} is not defined in the header");
                    }
                    if kind == "INFO" {
                        if needs_info {
                            return Err(clash(&format!("INFO/{s}")));
                        }
                        ops.push(RmOp::InfoTag(id.to_vec()));
                    } else {
                        ops.push(RmOp::FormatTag(id.to_vec()));
                    }
                    if declared {
                        hdr_out.remove(kind, id);
                    }
                }
            }
            None => {
                if keyword_match(s, "ID") {
                    ops.push(RmOp::Id);
                } else if keyword_match(s, "FILTER") {
                    ops.push(RmOp::FilterAll);
                    hdr_out.remove_all("FILTER", false);
                } else if keyword_match(s, "QUAL") {
                    ops.push(RmOp::Qual);
                } else if keyword_match(s, "INFO") {
                    if needs_info {
                        return Err(clash("INFO"));
                    }
                    ops.push(RmOp::InfoAll);
                    hdr_out.remove_all("INFO", false);
                } else if keyword_match(s, "FMT") || keyword_match(s, "FORMAT") {
                    ops.push(RmOp::FormatAll);
                    hdr_out.remove_all("FORMAT", true);
                } else if !s.is_empty() {
                    if hdr_in.has("INFO", id) {
                        return Err(format!("--remove {s}: did you mean INFO/{s}?"));
                    }
                    if hdr_in.has("FORMAT", id) {
                        return Err(format!("--remove {s}: did you mean FORMAT/{s}?"));
                    }
                    // `##key` removes the generic lines of that key. A bare key for a
                    // structured line (-x PEDIGREE) removes nothing in bcftools either.
                    if let Some(key) = s.strip_prefix("##") {
                        hdr_out.remove_generic(key.as_bytes());
                    }
                }
            }
        }
    }
    if keep_info || keep_fmt || keep_flt {
        for (kind, on) in [
            ("INFO", keep_info),
            ("FORMAT", keep_fmt),
            ("FILTER", keep_flt),
        ] {
            if !on {
                continue;
            }
            for id in hdr_in.ids(kind) {
                if keep.iter().any(|(k, i)| *k == kind && *i == id) {
                    if kind == "FILTER" && id == b"PASS" {
                        keep_pass = true;
                    }
                    continue;
                }
                ops.push(match kind {
                    "INFO" => RmOp::InfoTag(id.clone()),
                    "FORMAT" => RmOp::FormatTag(id.clone()),
                    _ => RmOp::FilterTag(id.clone()),
                });
                hdr_out.remove(kind, &id);
            }
        }
    }
    if ops.is_empty() {
        return Err(format!(
            "--remove {list}: no matching tag; name ID, QUAL, FILTER, INFO, FORMAT, INFO/TAG, FORMAT/TAG or FILTER/NAME (^ keeps the listed tags)"
        ));
    }
    Ok((ops, keep_pass))
}

impl Annotator {
    /// Compile `opts` against the FILE header `file_header` and the source.
    /// Errors are usage messages.
    pub fn new(
        opts: &AnnotateOptions,
        file_header: &[u8],
        source: SourceKind<'_>,
    ) -> Result<Annotator, String> {
        let hdr_in0 = HeaderEdit::parse(file_header);
        let mut hdr_in = hdr_in0.clone();
        let mut hdr_out = hdr_in0;

        let set_id = match &opts.set_id {
            Some(f) => {
                let fmt = parse_id_format(f)?;
                for tag in fmt.info_tags() {
                    if !hdr_in.has("INFO", tag) {
                        return Err(format!(
                            "--set-id {f}: no such tag defined in the VCF header: INFO/{} (declare it with --header-lines; tags added by --columns are not visible to --set-id)",
                            String::from_utf8_lossy(tag)
                        ));
                    }
                }
                Some(fmt)
            }
            None => None,
        };

        let (rm, keep_pass) = match &opts.remove {
            Some(list) => {
                let needs_info = set_id
                    .as_ref()
                    .is_some_and(|f| f.info_tags().next().is_some());
                compile_remove(
                    list,
                    &hdr_in,
                    &mut hdr_out,
                    needs_info,
                    opts.set_id.as_deref().unwrap_or(""),
                )?
            }
            None => (Vec::new(), false),
        };

        for l in &opts.header_lines {
            hdr_out.append(l)?;
            hdr_in.append(l)?;
        }

        let vcf_src = matches!(source, SourceKind::Vcf(_));
        let src_hdr = match source {
            SourceKind::Vcf(h) => Some(HeaderEdit::parse(h)),
            _ => None,
        };
        let mut cols = Vec::new();
        let mut layout = None;
        if let Some(list) = &opts.columns {
            let spec = parse_columns(list, vcf_src)?;
            match source {
                SourceKind::None => {
                    return Err(
                        "--columns needs --annotations: it names what to transfer from that file"
                            .to_string(),
                    );
                }
                SourceKind::Vcf(_) => {
                    cols = xfer::resolve_vcf(&spec, &mut hdr_out, src_hdr.as_ref().expect("vcf"))?;
                }
                SourceKind::Tab { bed } => {
                    let (c, l) = xfer::resolve_tab(&spec, &mut hdr_out, bed)?;
                    cols = c;
                    layout = Some(l);
                }
            }
        } else if let SourceKind::Tab { .. } = source {
            return Err(
                "--annotations FILE (not a VCF/BCF) needs --columns, e.g. -c CHROM,POS,ID or -c CHROM,FROM,TO,INFO/TAG"
                    .to_string(),
            );
        }

        let min_overlap = match &opts.min_overlap {
            None => MinOverlap::default(),
            Some(s) => {
                if vcf_src || layout.as_ref().is_none_or(|l| l.end.is_none()) {
                    return Err(
                        "--min-overlap applies to intervals: use a BED file or --columns CHROM,FROM,TO,..."
                            .to_string(),
                    );
                }
                parse_min_overlap(s)?
            }
        };

        let mark = match &opts.mark_sites {
            None => None,
            Some(s) => {
                let (listed, tag) = match s.as_bytes().first() {
                    Some(b'+') => (true, &s[1..]),
                    Some(b'-') => (false, &s[1..]),
                    _ => (true, s.as_str()),
                };
                if tag.is_empty() {
                    return Err("--mark-sites needs a tag name, e.g. --mark-sites +KNOWN".into());
                }
                if matches!(source, SourceKind::None) {
                    if !listed {
                        return Err(format!(
                            "--mark-sites -{tag} needs --annotations: without a file every site is listed"
                        ));
                    }
                    log::info!(
                        "vcf annotate: no --annotations given, every site is marked with INFO/{tag}"
                    );
                }
                let desc = if matches!(source, SourceKind::None) {
                    format!("Sites marked with `gloci vcf annotate --mark-sites {s}`")
                } else {
                    format!("Sites {}listed in {tag}", if listed { "" } else { "not " })
                };
                hdr_out.append(
                    format!("##INFO=<ID={tag},Number=0,Type=Flag,Description=\"{desc}\">")
                        .as_bytes(),
                )?;
                Some(Mark {
                    tag: tag.as_bytes().to_vec(),
                    listed,
                })
            }
        };

        Ok(Annotator {
            cols,
            rm,
            keep_pass,
            set_id,
            mark,
            layout,
            min_overlap,
            has_source: !matches!(source, SourceKind::None),
            out_header: hdr_out,
            vcmp: vcmp::Vcmp::default(),
            scratch: Vec::new(),
        })
    }

    /// The output header text.
    pub fn out_header(&self) -> Vec<u8> {
        self.out_header.to_bytes()
    }

    /// How tab/BED lines are read; `None` for a VCF source or no source.
    pub fn layout(&self) -> Option<&TabLayout> {
        self.layout.as_ref()
    }

    pub fn min_overlap(&self) -> MinOverlap {
        self.min_overlap
    }

    /// Whether any edit (removal, transfer, ID, mark) applies to records.
    pub fn edits_records(&self) -> bool {
        !self.rm.is_empty() || !self.cols.is_empty() || self.set_id.is_some() || self.mark.is_some()
    }

    /// Whether a transfer or mark needs a source file read.
    pub fn needs_source(&self) -> bool {
        self.has_source
    }

    /// Run the removals on `rec`.
    pub fn apply_remove(&self, rec: &mut Rec) {
        for op in &self.rm {
            match op {
                RmOp::Id => {
                    if rec.id != b"." {
                        rec.id = b".".to_vec();
                        rec.dirty = true;
                    }
                }
                RmOp::Qual => {
                    if rec.qual != b"." {
                        rec.qual = b".".to_vec();
                        rec.dirty = true;
                    }
                }
                RmOp::FilterAll => {
                    if rec.filter != b"." {
                        rec.filter = b".".to_vec();
                        rec.dirty = true;
                    }
                }
                RmOp::FilterTag(id) => {
                    let parts: Vec<&[u8]> = rec.filter.split(|&b| b == b';').collect();
                    if parts.contains(&id.as_slice()) {
                        let rest: Vec<&[u8]> =
                            parts.into_iter().filter(|p| *p != id.as_slice()).collect();
                        rec.filter = if rest.is_empty() {
                            if self.keep_pass {
                                b"PASS".to_vec()
                            } else {
                                b".".to_vec()
                            }
                        } else {
                            rest.join(&b';')
                        };
                        rec.dirty = true;
                    }
                }
                RmOp::InfoAll => {
                    if rec.info != b"." {
                        rec.info = b".".to_vec();
                        rec.dirty = true;
                    }
                }
                RmOp::InfoTag(id) => {
                    if info_remove(&mut rec.info, id) {
                        rec.dirty = true;
                    }
                }
                RmOp::FormatAll => {
                    if format_retain(&mut rec.tail, |k| k == b"GT") {
                        rec.dirty = true;
                    }
                }
                RmOp::FormatTag(id) => {
                    if format_retain(&mut rec.tail, |k| k != id.as_slice()) {
                        rec.dirty = true;
                    }
                }
            }
        }
    }

    /// Transfer the columns from the paired source record `src`.
    pub fn apply_vcf(&mut self, rec: &mut Rec, src: &Rec) -> Result<(), String> {
        for i in 0..self.cols.len() {
            let col = std::mem::replace(&mut self.cols[i], Col::placeholder());
            let r = xfer::vcf_set(&col, rec, src, &mut self.vcmp);
            self.cols[i] = col;
            r?;
        }
        Ok(())
    }

    /// Transfer the columns from tab/BED `rows` that match the record, in
    /// the order bcftools visits them (a column takes the first row that sets it).
    pub fn apply_rows(&mut self, rec: &mut Rec, rows: &[&TabRow]) -> Result<(), String> {
        let mut done = vec![false; self.cols.len()];
        for row in rows {
            for (i, d) in done.iter_mut().enumerate() {
                if *d {
                    continue;
                }
                let col = std::mem::replace(&mut self.cols[i], Col::placeholder());
                let r = xfer::tab_set(
                    &col,
                    rec,
                    row,
                    self.layout.as_ref().expect("tab layout"),
                    &mut self.vcmp,
                    &hdr_filter_ok(&self.out_header),
                );
                self.cols[i] = col;
                if r? {
                    *d = true;
                }
            }
        }
        Ok(())
    }

    /// Whether the tab line `row` (whose footprint overlaps the record at
    /// 0-based inclusive `[beg, end]`) annotates `rec`: with REF and ALT
    /// columns, REF must agree and one ALT must be among the record's.
    pub fn tab_matches(&mut self, rec: &Rec, beg: i64, row: &TabRow) -> bool {
        let layout = self.layout.as_ref().expect("tab layout");
        let (Some(_), Some(_)) = (layout.ref_, layout.alt) else {
            return true;
        };
        if row.start != beg {
            return false;
        }
        let Some((als, nals)) = row.alleles(layout) else {
            return false;
        };
        if self.vcmp.set_ref(&rec.reference, &als[0]).is_none() {
            return false;
        }
        let alts = rec.alts();
        for al in &als[1..nals] {
            if alts.is_empty() && al.as_slice() == b"." {
                return true;
            }
            if self.vcmp.find_allele(&alts, al).is_some() {
                return true;
            }
        }
        false
    }

    /// The last steps of a record: `--set-id`, then `--mark-sites`.
    pub fn finish(&mut self, rec: &mut Rec, has_overlap: bool) {
        if let Some(fmt) = &self.set_id {
            self.scratch.clear();
            fmt.render(rec, &mut self.scratch);
            if !self.scratch.is_empty() && (!fmt.plus || rec.id == b"." || rec.id.is_empty()) {
                std::mem::swap(&mut rec.id, &mut self.scratch);
                rec.dirty = true;
            }
        }
        if let Some(m) = &self.mark {
            let set = if !self.has_source || has_overlap {
                m.listed
            } else {
                !m.listed
            };
            if set {
                info_set(&mut rec.info, &m.tag, None);
                rec.dirty = true;
            } else if info_remove(&mut rec.info, &m.tag) {
                rec.dirty = true;
            }
        }
    }

    /// Whether the record needs its INFO column or FORMAT parsed (a hint for callers).
    pub fn mark_listed(&self) -> bool {
        self.mark.is_some()
    }
}

/// The FILTER names the output header declares (a tab FILTER column must name one).
fn hdr_filter_ok(h: &HeaderEdit) -> impl Fn(&[u8]) -> bool + '_ {
    move |id| id == b"PASS" || h.has("FILTER", id)
}

fn parse_min_overlap(s: &str) -> Result<MinOverlap, String> {
    let bad = || {
        format!("--min-overlap {s}: expected ANN, ANN:VCF or :VCF with fractions between 0 and 1")
    };
    let mut o = MinOverlap::default();
    let (a, v) = match s.split_once(':') {
        Some((a, v)) => (a, Some(v)),
        None => (s, None),
    };
    if !a.is_empty() {
        o.ann = a.parse().map_err(|_| bad())?;
    }
    if let Some(v) = v {
        o.vcf = v.parse().map_err(|_| bad())?;
    }
    if !(0.0..=1.0).contains(&o.ann) || !(0.0..=1.0).contains(&o.vcf) {
        return Err(bad());
    }
    Ok(o)
}

#[cfg(test)]
mod tests {
    use super::*;

    const FILE_H: &[u8] = b"##fileformat=VCFv4.2\n##INFO=<ID=A1,Number=A,Type=Float,Description=\"a\">\n##INFO=<ID=A2,Number=A,Type=Integer,Description=\"a\">\n##INFO=<ID=S,Number=1,Type=String,Description=\"s\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    const SRC_H: &[u8] = b"##fileformat=VCFv4.2\n##INFO=<ID=A1,Number=A,Type=Float,Description=\"a\">\n##INFO=<ID=A2,Number=A,Type=Integer,Description=\"a\">\n##INFO=<ID=S,Number=1,Type=String,Description=\"s\">\n##INFO=<ID=NEW,Number=1,Type=String,Description=\"n\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    fn rec(line: &str) -> Rec {
        let mut r = Rec::default();
        r.load(line.as_bytes()).unwrap();
        r
    }

    fn run(cols: &str, dst: &str, src: &str) -> String {
        let opts = AnnotateOptions {
            columns: Some(cols.into()),
            ..AnnotateOptions::default()
        };
        let mut a = Annotator::new(&opts, FILE_H, SourceKind::Vcf(SRC_H)).unwrap();
        let mut d = rec(dst);
        a.apply_vcf(&mut d, &rec(src)).unwrap();
        let mut out = Vec::new();
        d.write(&mut out);
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn number_a_follows_the_allele() {
        let got = run(
            "INFO/A1",
            "1\t5\t.\tC\tT,G\t.\t.\t.",
            "1\t5\t.\tC\tG,T\t.\t.\tA1=0.1,0.2",
        );
        assert_eq!(got, "1\t5\t.\tC\tT,G\t.\t.\tA1=0.2,0.1");
        // An allele the source lacks stays missing.
        let got = run(
            "INFO/A1",
            "1\t5\t.\tC\tT,G\t.\t.\t.",
            "1\t5\t.\tC\tG\t.\t.\tA1=0.1",
        );
        assert_eq!(got, "1\t5\t.\tC\tT,G\t.\t.\tA1=.,0.1");
    }

    #[test]
    fn plus_keeps_existing_float_values_but_not_integer_ones() {
        let dst = "1\t5\t.\tC\tT,G\t.\t.\tA1=0.9,.;A2=7,8";
        let src = "1\t5\t.\tC\tT,G\t.\t.\tA1=0.1,0.2;A2=1,2";
        assert_eq!(
            run("+INFO/A1,+INFO/A2", dst, src),
            "1\t5\t.\tC\tT,G\t.\t.\tA1=0.9,0.2;A2=1,2"
        );
        assert_eq!(
            run("INFO/A1,INFO/A2", dst, src),
            "1\t5\t.\tC\tT,G\t.\t.\tA1=0.1,0.2;A2=1,2"
        );
    }

    #[test]
    fn a_wrong_count_and_incompatible_refs_are_errors() {
        let opts = AnnotateOptions {
            columns: Some("INFO/A1".into()),
            ..AnnotateOptions::default()
        };
        let mut a = Annotator::new(&opts, FILE_H, SourceKind::Vcf(SRC_H)).unwrap();
        let mut d = rec("1\t5\t.\tC\tT\t.\t.\t.");
        let e = a
            .apply_vcf(&mut d, &rec("1\t5\t.\tC\tT,G\t.\t.\tA1=0.1"))
            .unwrap_err();
        assert!(
            e.contains("incorrect number of values (1)") && e.contains("1:5"),
            "{e}"
        );
        let e = a
            .apply_vcf(&mut d, &rec("1\t5\t.\tG\tT\t.\t.\tA1=0.1"))
            .unwrap_err();
        assert!(e.contains("REF alleles not compatible"), "{e}");
    }

    #[test]
    fn rename_flag_and_missing_values() {
        assert_eq!(
            run(
                "S:=NEW",
                "1\t5\t.\tC\tT\t.\t.\t.",
                "1\t5\t.\tC\tT\t.\t.\tNEW=x"
            ),
            "1\t5\t.\tC\tT\t.\t.\tS=x"
        );
        // A source value of . is carried over for a VCF source.
        assert_eq!(
            run(
                "INFO/S",
                "1\t5\t.\tC\tT\t.\t.\tS=a",
                "1\t5\t.\tC\tT\t.\t.\tS=."
            ),
            "1\t5\t.\tC\tT\t.\t.\tS=."
        );
        // +TAG replaces a missing value.
        assert_eq!(
            run(
                "+INFO/S",
                "1\t5\t.\tC\tT\t.\t.\tS=.",
                "1\t5\t.\tC\tT\t.\t.\tS=b"
            ),
            "1\t5\t.\tC\tT\t.\t.\tS=b"
        );
    }

    #[test]
    fn tab_rows_match_on_ref_and_an_alt() {
        let opts = AnnotateOptions {
            columns: Some("CHROM,POS,REF,ALT,ID".into()),
            ..AnnotateOptions::default()
        };
        let mut a = Annotator::new(&opts, FILE_H, SourceKind::Tab { bed: false }).unwrap();
        let layout = a.layout().unwrap().clone();
        let row = |l: &str| TabRow::parse(l.as_bytes(), &layout).unwrap().unwrap();
        let r = rec("1\t5\t.\tC\tT,G\t.\t.\t.");
        assert!(a.tab_matches(&r, 4, &row("1\t5\tC\tG\trs1")));
        assert!(!a.tab_matches(&r, 4, &row("1\t5\tC\tA\trs1")));
        assert!(!a.tab_matches(&r, 4, &row("1\t5\tA\tG\trs1")));
        assert!(!a.tab_matches(&r, 4, &row("1\t6\tC\tG\trs1")));
        // A longer REF in the file matches an ALT that carries the extra bases.
        assert!(a.tab_matches(&rec("1\t5\t.\tC\tT\t.\t.\t."), 4, &row("1\t5\tCA\tTA\trs1")));
    }

    #[test]
    fn removals_and_set_id_and_mark() {
        let opts = AnnotateOptions {
            remove: Some("INFO/S,ID".into()),
            set_id: Some("+%CHROM:%POS:%REF:%FIRST_ALT".into()),
            mark_sites: Some("-NOHIT".into()),
            ..AnnotateOptions::default()
        };
        let mut a = Annotator::new(&opts, FILE_H, SourceKind::Vcf(SRC_H)).unwrap();
        let mut d = rec("1\t5\trs9\tC\tT,G\t.\t.\tS=x;A2=1");
        a.apply_remove(&mut d);
        a.finish(&mut d, false);
        let mut out = Vec::new();
        d.write(&mut out);
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "1\t5\t1:5:C:T\tC\tT,G\t.\t.\tA2=1;NOHIT"
        );
        assert!(!String::from_utf8(a.out_header()).unwrap().contains("ID=S,"));
    }
}
