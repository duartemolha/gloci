//! The merged header of `vcf merge`: bcftools' `bcf_hdr_merge` of the inputs'
//! headers, and the merged sample list.

use crate::annotate::hdr::{InfoDef, Number, Ty, def_of, structured};
use std::collections::HashMap;

/// The standard FILTER line htslib puts after `##fileformat`.
const PASS_LINE: &[u8] = b"##FILTER=<ID=PASS,Description=\"All filters passed\">";

/// One input's header, as text up to and including `#CHROM`.
pub struct HeaderIn<'a> {
    pub name: &'a str,
    pub raw: &'a [u8],
}

/// How the tags of the merged header are declared (the first declaration of
/// each ID wins), by header position.
#[derive(Default, Debug, Clone)]
pub struct Defs {
    pub info: HashMap<Vec<u8>, (usize, InfoDef)>,
    pub format: HashMap<Vec<u8>, (usize, InfoDef)>,
    /// `##INFO` IDs in header order.
    pub info_order: Vec<Vec<u8>>,
    /// `##FORMAT` IDs in header order.
    pub format_order: Vec<Vec<u8>>,
}

/// The result of merging the headers.
pub struct Merged {
    /// The header text, `#CHROM` line included, every line ending in `\n`.
    pub text: Vec<u8>,
    /// The merged sample names, in file order (renamed when forced).
    pub samples: Vec<Vec<u8>>,
    /// Samples per input.
    pub counts: Vec<usize>,
    pub defs: Defs,
    /// Per input, the declared type of each INFO and FORMAT tag in its own header.
    pub own: Vec<Defs>,
}

fn trim(line: &[u8]) -> &[u8] {
    let line = line.strip_suffix(b"\n").unwrap_or(line);
    line.strip_suffix(b"\r").unwrap_or(line)
}

fn lines(raw: &[u8]) -> Vec<&[u8]> {
    raw.split(|&b| b == b'\n')
        .map(trim)
        .filter(|l| !l.is_empty())
        .collect()
}

fn samples_of(line: &[u8]) -> Vec<Vec<u8>> {
    line.split(|&b| b == b'\t')
        .skip(9)
        .map(<[u8]>::to_vec)
        .collect()
}

fn generic_key(line: &[u8]) -> &[u8] {
    let rest = &line[2.min(line.len())..];
    match rest.iter().position(|&b| b == b'=') {
        Some(p) => &rest[..p],
        None => rest,
    }
}

fn record_def(defs: &mut Defs, kind: &[u8], id: &[u8], line: &[u8]) {
    let (map, order) = match kind {
        b"INFO" => (&mut defs.info, &mut defs.info_order),
        b"FORMAT" => (&mut defs.format, &mut defs.format_order),
        _ => return,
    };
    if !map.contains_key(id) {
        map.insert(id.to_vec(), (order.len(), def_of(line)));
        order.push(id.to_vec());
    }
}

fn ty_name(t: Ty) -> &'static str {
    match t {
        Ty::Flag => "Flag",
        Ty::Int => "Integer",
        Ty::Float => "Float",
        Ty::Str => "String",
    }
}

fn number_name(n: Number) -> String {
    match n {
        Number::Fixed(k) => k.to_string(),
        Number::A => "A".into(),
        Number::R => "R".into(),
        Number::G => "G".into(),
        Number::Var => ".".into(),
    }
}

/// Merge the headers. `force_samples` renames a repeated sample name to
/// `N:NAME`, N being the 1-based number of its file; otherwise a repeat is
/// an error naming both files.
pub fn merge_headers(inputs: &[HeaderIn<'_>], force_samples: bool) -> Result<Merged, String> {
    let mut out: Vec<Vec<u8>> = Vec::new();
    let mut seen: HashMap<(Vec<u8>, Vec<u8>), usize> = HashMap::new();
    let mut generic_seen: HashMap<Vec<u8>, ()> = HashMap::new();
    let mut defs = Defs::default();
    let mut own = Vec::new();
    let mut first_decl: HashMap<(Vec<u8>, Vec<u8>), usize> = HashMap::new();
    let mut fileformat: Option<Vec<u8>> = None;
    let mut names: Vec<Vec<u8>> = Vec::new();
    let mut owner: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut counts = Vec::new();

    for (ri, h) in inputs.iter().enumerate() {
        let mut mine = Defs::default();
        let mut chrom_line: &[u8] = b"";
        for line in lines(h.raw) {
            if line.starts_with(b"#CHROM") {
                chrom_line = line;
                continue;
            }
            if !line.starts_with(b"##") {
                continue;
            }
            if line.starts_with(b"##fileformat=") {
                if fileformat.is_none() {
                    fileformat = Some(line.to_vec());
                }
                continue;
            }
            match structured(line) {
                Some((kind, id)) => {
                    record_def(&mut mine, kind, &id, line);
                    if kind == b"FILTER" && id == b"PASS" {
                        continue;
                    }
                    let key = (kind.to_vec(), id.clone());
                    if let Some(&first) = first_decl.get(&key) {
                        // Same ID again: keep the first, warn when its type or length differ.
                        if kind == b"INFO" || kind == b"FORMAT" {
                            let a = &defs_of(&defs, kind).get(&id).map(|x| x.1);
                            let b = def_of(line);
                            if let Some(a) = a {
                                let id_s = String::from_utf8_lossy(&id);
                                let k = String::from_utf8_lossy(kind);
                                if a.ty != b.ty {
                                    log::warn!(
                                        "vcf merge: {k}/{id_s} is declared as {} in {} but as {} in {}; the first declaration is kept, and a record value of the other type is an error",
                                        ty_name(a.ty),
                                        inputs[first].name,
                                        ty_name(b.ty),
                                        h.name
                                    );
                                } else if a.number != b.number {
                                    log::warn!(
                                        "vcf merge: {k}/{id_s} has Number={} in {} but Number={} in {}; the first declaration is kept",
                                        number_name(a.number),
                                        inputs[first].name,
                                        number_name(b.number),
                                        h.name
                                    );
                                }
                            }
                        }
                        continue;
                    }
                    first_decl.insert(key.clone(), ri);
                    seen.insert(key, out.len());
                    record_def(&mut defs, kind, &id, line);
                    out.push(line.to_vec());
                }
                None => {
                    let key = generic_key(line).to_vec();
                    // The first input keeps every line; later ones add a key not yet seen.
                    if ri == 0 || !generic_seen.contains_key(&key) {
                        out.push(line.to_vec());
                    }
                    generic_seen.insert(key, ());
                }
            }
        }
        own.push(mine);
        let ss = samples_of(chrom_line);
        counts.push(ss.len());
        for s in ss {
            if let Some(&prev) = owner.get(&s) {
                if !force_samples {
                    return Err(format!(
                        "sample name {} is in both {} and {}; sample names must be distinct (pass --force-samples to rename the later one to {}:{}, or rename it with bcftools reheader)",
                        String::from_utf8_lossy(&s),
                        inputs[prev].name,
                        h.name,
                        ri + 1,
                        String::from_utf8_lossy(&s)
                    ));
                }
                let mut renamed = format!("{}:", ri + 1).into_bytes();
                renamed.extend_from_slice(&s);
                names.push(renamed);
            } else {
                owner.insert(s.clone(), ri);
                names.push(s);
            }
        }
    }
    let _ = &seen;
    let mut text = Vec::new();
    text.extend_from_slice(fileformat.as_deref().unwrap_or(b"##fileformat=VCFv4.2"));
    text.push(b'\n');
    text.extend_from_slice(PASS_LINE);
    text.push(b'\n');
    for l in &out {
        text.extend_from_slice(l);
        text.push(b'\n');
    }
    text.extend_from_slice(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    if !names.is_empty() {
        text.extend_from_slice(b"\tFORMAT");
        for n in &names {
            text.push(b'\t');
            text.extend_from_slice(n);
        }
    }
    text.push(b'\n');
    Ok(Merged {
        text,
        samples: names,
        counts,
        defs,
        own,
    })
}

fn defs_of<'a>(d: &'a Defs, kind: &[u8]) -> &'a HashMap<Vec<u8>, (usize, InfoDef)> {
    if kind == b"INFO" { &d.info } else { &d.format }
}

/// Declare INFO/`tag` as `Number=.`: its header line moves to the end of the
/// ## lines with Number replaced.
pub(crate) fn make_variable(merged: &mut Merged, tag: &[u8]) {
    let mut lines: Vec<Vec<u8>> = merged
        .text
        .split_inclusive(|&b| b == b'\n')
        .map(<[u8]>::to_vec)
        .collect();
    let prefix = {
        let mut p = b"##INFO=<ID=".to_vec();
        p.extend_from_slice(tag);
        p
    };
    let Some(pos) = lines
        .iter()
        .position(|l| l.starts_with(&prefix) && matches!(l.get(prefix.len()), Some(b',' | b'>')))
    else {
        return;
    };
    let old = lines.remove(pos);
    let text = String::from_utf8_lossy(&old).into_owned();
    let new = match text.find("Number=") {
        Some(i) => {
            let end = text[i..].find(',').map_or(text.len(), |e| i + e);
            format!("{}Number=.{}", &text[..i], &text[end..])
        }
        None => text,
    };
    let chrom = lines
        .iter()
        .position(|l| l.starts_with(b"#CHROM"))
        .unwrap_or(lines.len());
    lines.insert(chrom, new.into_bytes());
    merged.text = lines.concat();
    if let Some((_, d)) = merged.defs.info.get_mut(tag) {
        d.number = Number::Var;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const C: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

    fn h(lines: &[&str], samples: &[&str]) -> Vec<u8> {
        let mut t = String::from("##fileformat=VCFv4.2\n");
        for l in lines {
            t.push_str(l);
            t.push('\n');
        }
        t.push_str(C);
        for s in samples {
            t.push('\t');
            t.push_str(s);
        }
        t.push('\n');
        t.into_bytes()
    }

    #[test]
    fn first_header_wins_and_new_lines_are_added() {
        let a = h(
            &[
                "##source=a",
                "##contig=<ID=1,length=10>",
                "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"a\">",
            ],
            &["A"],
        );
        let b = h(
            &[
                "##source=b",
                "##x=1",
                "##contig=<ID=1,length=20>",
                "##contig=<ID=2>",
                "##INFO=<ID=DP,Number=2,Type=Float,Description=\"b\">",
                "##INFO=<ID=NEW,Number=.,Type=String,Description=\"n\">",
            ],
            &["B"],
        );
        let m = merge_headers(
            &[
                HeaderIn { name: "a", raw: &a },
                HeaderIn { name: "b", raw: &b },
            ],
            false,
        )
        .unwrap();
        let t = String::from_utf8(m.text).unwrap();
        let want = "##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##source=a\n##contig=<ID=1,length=10>\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"a\">\n\
##x=1\n##contig=<ID=2>\n##INFO=<ID=NEW,Number=.,Type=String,Description=\"n\">\n";
        assert!(t.starts_with(want), "{t}");
        assert!(t.ends_with("FORMAT\tA\tB\n"));
        assert_eq!(m.counts, vec![1, 1]);
        assert_eq!(m.defs.info[&b"DP"[..]].1.number, Number::Fixed(1));
    }

    #[test]
    fn repeated_sample_is_an_error_or_renamed_with_its_file_number() {
        let a = h(&[], &["S", "T"]);
        let b = h(&[], &["T", "U"]);
        let ins = [
            HeaderIn {
                name: "a.vcf",
                raw: &a,
            },
            HeaderIn {
                name: "b.vcf",
                raw: &b,
            },
        ];
        let e = merge_headers(&ins, false).err().unwrap();
        assert!(
            e.contains("T") && e.contains("a.vcf") && e.contains("b.vcf"),
            "{e}"
        );
        let m = merge_headers(&ins, true).unwrap();
        let names: Vec<String> = m
            .samples
            .iter()
            .map(|s| String::from_utf8_lossy(s).into_owned())
            .collect();
        assert_eq!(names, ["S", "T", "2:T", "U"]);
    }
}
