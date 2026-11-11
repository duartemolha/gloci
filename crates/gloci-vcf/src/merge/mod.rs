//! `gloci vcf merge`: the record and header rules of `bcftools merge`
//! (1.21): which records merge is in [`pair`]; this module merges the
//! records of one row into one. The streaming read of the inputs is in
//! gloci-ops.

mod format;
pub mod header;
mod info;
pub mod pair;

use crate::annotate::hdr::{InfoDef, Number, Ty};
use header::{Defs, Merged};
pub use pair::{MergeMode, PairIn, StagedRow, stage_position};

/// Options of one merge.
#[derive(Clone, Debug)]
pub struct MergeOptions {
    pub mode: MergeMode,
    /// A sample absent from a record gets `0/0` instead of `./.`.
    pub missing_to_ref: bool,
    /// `--info-rules`: `None` is bcftools' default, `-` none at all.
    pub info_rules: Option<String>,
    /// `-F x`: a PASS in any input makes the merged FILTER PASS.
    pub filter_x: bool,
}

/// The columns of one input record.
pub(crate) struct P<'a> {
    pub reader: usize,
    pub chrom: &'a [u8],
    pub pos: &'a [u8],
    pub id: &'a [u8],
    pub reference: &'a [u8],
    pub alt: &'a [u8],
    pub qual: &'a [u8],
    pub filter: &'a [u8],
    pub info: &'a [u8],
    pub fmt: Vec<&'a [u8]>,
    pub samples: Vec<&'a [u8]>,
}

impl<'a> P<'a> {
    fn parse(reader: usize, line: &'a [u8]) -> Result<P<'a>, String> {
        let cols: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
        if cols.len() < 8 {
            return Err(format!(
                "the record has {} tab-separated columns; a VCF record needs at least 8",
                cols.len()
            ));
        }
        let fmt = if cols.len() > 8 && cols[8] != b"." {
            cols[8].split(|&b| b == b':').collect()
        } else {
            Vec::new()
        };
        Ok(P {
            reader,
            chrom: cols[0],
            pos: cols[1],
            id: cols[2],
            reference: cols[3],
            alt: cols[4],
            qual: cols[5],
            filter: cols[6],
            info: cols[7],
            fmt,
            samples: cols.get(9..).unwrap_or(&[]).to_vec(),
        })
    }

    pub(crate) fn alts(&self) -> Vec<&'a [u8]> {
        if self.alt == b"." || self.alt.is_empty() {
            Vec::new()
        } else {
            self.alt.split(|&b| b == b',').collect()
        }
    }
}

/// The merged allele list and each input's map into it.
pub struct Alleles {
    pub reference: Vec<u8>,
    pub alts: Vec<Vec<u8>>,
    /// Per row member (same order as the records): input allele index to merged index.
    pub maps: Vec<Vec<usize>>,
    /// Per row member: its alleles are not the merged list as it is.
    pub differs: Vec<bool>,
}

/// What a merge needs of one input.
struct Input {
    name: String,
    samples: usize,
}

/// A configured merge of the inputs whose headers were merged.
pub struct Merger {
    opts: MergeOptions,
    defs: Defs,
    own: Vec<Defs>,
    inputs: Vec<Input>,
    total_samples: usize,
    rules: Vec<info::Rule>,
    /// AC and AN are removed and recomputed from the merged genotypes.
    calc_acan: bool,
    scratch: Vec<u8>,
}

fn pos_label(r: &P<'_>) -> String {
    format!(
        "{}:{}",
        String::from_utf8_lossy(r.chrom),
        String::from_utf8_lossy(r.pos)
    )
}

impl Merger {
    pub fn new(
        opts: MergeOptions,
        merged: &mut Merged,
        names: &[String],
    ) -> Result<Merger, String> {
        let (rules, calc_acan) = info::build_rules(&opts, merged)?;
        // `join` of a tag whose length follows the alleles makes it Number=.,
        // and bcftools moves its header line to the end.
        for r in &rules {
            if r.method == info::Method::Join
                && matches!(r.def.number, Number::A | Number::R | Number::G)
            {
                header::make_variable(merged, &r.tag);
            }
        }
        Ok(Merger {
            opts,
            defs: merged.defs.clone(),
            own: merged.own.clone(),
            inputs: names
                .iter()
                .zip(&merged.counts)
                .map(|(n, &c)| Input {
                    name: n.clone(),
                    samples: c,
                })
                .collect(),
            total_samples: merged.samples.len(),
            rules,
            calc_acan,
            scratch: Vec::new(),
        })
    }

    /// What pairing needs of `line`.
    pub fn pair_in(&self, line: &[u8]) -> Result<PairIn, String> {
        let p = P::parse(0, line)?;
        Ok(PairIn {
            reference: p.reference.to_vec(),
            alts: p.alts().iter().map(|a| a.to_vec()).collect(),
            id: p.id.to_vec(),
            end: info::end_value(p.info),
        })
    }

    pub fn options(&self) -> &MergeOptions {
        &self.opts
    }

    /// Merge the records of one row (`row[i]` is input i's record line, no
    /// terminator) and append the merged line to `out`.
    pub fn merge_row(
        &mut self,
        row: &[Option<&[u8]>],
        al: &Alleles,
        out: &mut Vec<u8>,
    ) -> Result<(), String> {
        let mut recs = Vec::new();
        for (i, l) in row.iter().enumerate() {
            if let Some(l) = l {
                let p = P::parse(i, l).map_err(|e| format!("{}: {e}", self.inputs[i].name))?;
                if p.samples.len() != self.inputs[i].samples {
                    return Err(format!(
                        "{}: the record at {} has {} sample columns but the header has {}",
                        self.inputs[i].name,
                        pos_label(&p),
                        p.samples.len(),
                        self.inputs[i].samples
                    ));
                }
                recs.push(p);
            }
        }
        if recs.is_empty() {
            return Ok(());
        }
        let first = &recs[0];
        out.extend_from_slice(first.chrom);
        out.push(b'\t');
        out.extend_from_slice(first.pos);
        out.push(b'\t');
        merge_ids(&recs, out);
        out.push(b'\t');
        out.extend_from_slice(&al.reference);
        out.push(b'\t');
        if al.alts.is_empty() {
            out.push(b'.');
        } else {
            for (i, a) in al.alts.iter().enumerate() {
                if i > 0 {
                    out.push(b',');
                }
                out.extend_from_slice(a);
            }
        }
        out.push(b'\t');
        merge_qual(&recs, out);
        out.push(b'\t');
        merge_filter(&recs, self.opts.filter_x, out);
        out.push(b'\t');
        let mut fmt = format::Formatted::default();
        if self.total_samples > 0 {
            fmt = format::merge_format(self, &recs, al)?;
        }
        let mut info_out = std::mem::take(&mut self.scratch);
        info_out.clear();
        info::merge_info(self, &recs, al, &fmt.calls, &mut info_out)?;
        if info_out.is_empty() {
            out.push(b'.');
        } else {
            out.extend_from_slice(&info_out);
        }
        self.scratch = info_out;
        if self.total_samples > 0 {
            out.extend_from_slice(&fmt.text);
        }
        Ok(())
    }

    pub(crate) fn out_def(&self, kind_info: bool, key: &[u8]) -> InfoDef {
        let m = if kind_info {
            &self.defs.info
        } else {
            &self.defs.format
        };
        m.get(key).map_or(
            InfoDef {
                ty: Ty::Str,
                number: Number::Var,
            },
            |d| d.1,
        )
    }

    /// An error when input `reader` declares `key` with another type than the merged header.
    pub(crate) fn check_type(
        &self,
        kind_info: bool,
        reader: usize,
        key: &[u8],
        at: &str,
    ) -> Result<(), String> {
        let own = if kind_info {
            self.own[reader].info.get(key)
        } else {
            self.own[reader].format.get(key)
        };
        let Some((_, d)) = own else { return Ok(()) };
        let out = self.out_def(kind_info, key);
        if d.ty != out.ty {
            return Err(format!(
                "{}: {}/{} at {at} is declared with a different type than in the first input that declares it; make the ##{} lines agree (for example with bcftools annotate -h) before merging",
                self.inputs[reader].name,
                if kind_info { "INFO" } else { "FORMAT" },
                String::from_utf8_lossy(key),
                if kind_info { "INFO" } else { "FORMAT" },
            ));
        }
        Ok(())
    }
}

fn merge_ids(recs: &[P<'_>], out: &mut Vec<u8>) {
    let mut seen: Vec<&[u8]> = Vec::new();
    for r in recs {
        if r.id == b"." || r.id.is_empty() || seen.contains(&r.id) {
            continue;
        }
        seen.push(r.id);
    }
    if seen.is_empty() {
        out.push(b'.');
    }
    for (i, s) in seen.iter().enumerate() {
        if i > 0 {
            out.push(b';');
        }
        out.extend_from_slice(s);
    }
}

pub(crate) fn parse_f32(tok: &[u8]) -> Option<f32> {
    std::str::from_utf8(tok)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .map(|v| v as f32)
}

pub(crate) fn put_f32(out: &mut Vec<u8>, v: f32) {
    gloci_core::cfmt::write_g6(out, f64::from(v));
}

fn merge_qual(recs: &[P<'_>], out: &mut Vec<u8>) {
    let mut best: Option<f32> = None;
    for r in recs {
        if let Some(v) = parse_f32(r.qual)
            && best.is_none_or(|b| v > b)
        {
            best = Some(v);
        }
    }
    match best {
        Some(v) => put_f32(out, v),
        None => out.push(b'.'),
    }
}

fn merge_filter(recs: &[P<'_>], filter_x: bool, out: &mut Vec<u8>) {
    // `-F x`: an input without a filter (. or PASS) makes the result PASS.
    if filter_x
        && recs
            .iter()
            .any(|r| r.filter == b"." || r.filter.is_empty() || r.filter == b"PASS")
    {
        out.extend_from_slice(b"PASS");
        return;
    }
    let mut names: Vec<&[u8]> = Vec::new();
    let mut pass = false;
    for r in recs {
        if r.filter == b"." || r.filter.is_empty() {
            continue;
        }
        for f in r.filter.split(|&b| b == b';') {
            if f == b"PASS" {
                pass = true;
            } else if !names.contains(&f) {
                names.push(f);
            }
        }
    }
    if names.is_empty() {
        out.extend_from_slice(if pass { b"PASS" } else { b"." });
        return;
    }
    for (i, f) in names.iter().enumerate() {
        if i > 0 {
            out.push(b';');
        }
        out.extend_from_slice(f);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(items: &[&str]) -> Vec<Vec<u8>> {
        items.iter().map(|s| s.as_bytes().to_vec()).collect()
    }

    #[test]
    fn common_trailing_bases_are_trimmed_but_one_base_stays() {
        let mut a = v(&["AT", "ATT"]);
        pair::trim_common_suffix(&mut a);
        assert_eq!(a, v(&["A", "AT"]));
        let mut a = v(&["AT"]);
        pair::trim_common_suffix(&mut a);
        assert_eq!(a, v(&["A"]));
        // a symbolic ALT or `*` stops the trim
        let mut a = v(&["AT", "<DEL>"]);
        pair::trim_common_suffix(&mut a);
        assert_eq!(a, v(&["AT", "<DEL>"]));
        let mut a = v(&["AT", "*"]);
        pair::trim_common_suffix(&mut a);
        assert_eq!(a, v(&["AT", "*"]));
        // a one-base REF is never trimmed
        let mut a = v(&["A", "AT"]);
        pair::trim_common_suffix(&mut a);
        assert_eq!(a, v(&["A", "AT"]));
    }

    fn merged(lines: &[&str]) -> Result<Vec<u8>, String> {
        let h = |sample: &str| {
            format!(
                "##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
            )
        };
        let hs: Vec<String> = (0..lines.len()).map(|i| h(&format!("S{i}"))).collect();
        let names: Vec<String> = (0..lines.len()).map(|i| format!("f{i}")).collect();
        let ins: Vec<header::HeaderIn<'_>> = hs
            .iter()
            .zip(&names)
            .map(|(h, n)| header::HeaderIn {
                name: n,
                raw: h.as_bytes(),
            })
            .collect();
        let mut m = header::merge_headers(&ins, false)?;
        let mut mg = Merger::new(
            MergeOptions {
                mode: MergeMode::Both,
                missing_to_ref: false,
                info_rules: None,
                filter_x: false,
            },
            &mut m,
            &names,
        )?;
        let ins: Vec<Vec<PairIn>> = lines
            .iter()
            .map(|l| vec![mg.pair_in(l.as_bytes()).unwrap()])
            .collect();
        let staged = stage_position(&ins, MergeMode::All)?;
        let mut out = Vec::new();
        for r in staged {
            let row: Vec<Option<&[u8]>> = r
                .members
                .iter()
                .enumerate()
                .map(|(i, m)| m.map(|_| lines[i].as_bytes()))
                .collect();
            mg.merge_row(&row, &r.alleles, &mut out)?;
        }
        Ok(out)
    }

    #[test]
    fn reference_is_extended_and_alleles_united_in_input_order() {
        let out = merged(&[
            "1\t5\t.\tA\tG\t10\tPASS\tDP=1\tGT:AD\t0/1:3,4",
            "1\t5\t.\tAT\tA\t20\tq1\tDP=2\tGT:AD\t1/1:0,9",
        ])
        .unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "1\t5\t.\tAT\tGT,A\t20\tq1\tDP=3\tGT:AD\t0/1:3,4,.\t2/2:0,.,9"
        );
    }

    #[test]
    fn reference_mismatch_is_an_error() {
        let e = merged(&[
            "1\t5\t.\tA\tG\t.\t.\t.\tGT\t0/1",
            "1\t5\t.\tC\tT\t.\t.\t.\tGT\t0/1",
        ])
        .unwrap_err();
        assert!(e.contains("REF prefixes differ"), "{e}");
    }

    #[test]
    fn genotype_with_an_allele_the_record_lacks_is_an_error() {
        let e = merged(&[
            "1\t5\t.\tA\tG\t.\t.\t.\tGT\t0/2",
            "1\t5\t.\tA\tG\t.\t.\t.\tGT\t0/1",
        ])
        .unwrap_err();
        assert!(e.contains("allele 2"), "{e}");
    }
}
