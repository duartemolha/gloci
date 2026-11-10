//! Chromosome naming checks, prefix difference hints and genome aliases (D58).

pub use crate::builtin::builtin_alias;
use gloci_core::ChromId;
use std::collections::HashSet;

/// Which chromosome ids a stream has already reported, so a list of the
/// distinct names seen stays linear in the number of contigs. One table gives
/// one id per name.
#[derive(Clone, Debug, Default)]
pub struct SeenChroms(Vec<bool>);

impl SeenChroms {
    /// True the first time `id` is passed.
    pub fn first_time(&mut self, id: ChromId) -> bool {
        let i = id.0 as usize;
        if self.0.len() <= i {
            self.0.resize(i + 1, false);
        }
        !std::mem::replace(&mut self.0[i], true)
    }
}

/// When adding or stripping `"chr"` makes at least one chromosome name match
/// between the two sets, returns the hint clause explaining the mismatch.
pub fn chr_prefix_hint(a_names: &[&[u8]], b_names: &[&[u8]]) -> Option<String> {
    // Hash sets: the sets can hold one name per contig of an assembly.
    let a_set: HashSet<&[u8]> = a_names.iter().copied().collect();
    let b_set: HashSet<&[u8]> = b_names.iter().copied().collect();
    fn stripped(n: &[u8]) -> Option<&[u8]> {
        n.strip_prefix(b"chr").filter(|r| !r.is_empty())
    }
    let explains = a_names
        .iter()
        .any(|&a| stripped(a).is_some_and(|r| b_set.contains(r)))
        || b_names
            .iter()
            .any(|&b| stripped(b).is_some_and(|r| a_set.contains(r)));
    if explains {
        Some(
            "the names differ by a \"chr\" prefix — pass --chrom-map FILE (lines of old name, tab, new name) to rename one file's chromosomes so they match"
                .to_string(),
        )
    } else {
        None
    }
}

/// Format up to 3 distinct chromosome names in seen order, appending `, …`
/// if there are more than 3.
pub fn format_chrom_list(names: &[impl AsRef<[u8]>]) -> String {
    const MAX_SHOWN: usize = 3;
    let shown: Vec<String> = names
        .iter()
        .take(MAX_SHOWN)
        .map(|c| String::from_utf8_lossy(c.as_ref()).into_owned())
        .collect();
    if names.len() > MAX_SHOWN {
        format!("{}, …", shown.join(", "))
    } else {
        shown.join(", ")
    }
}

/// Format the warning message (for `log::warn!`) when two inputs share no chromosome names.
pub fn format_disjoint_chroms_warning(
    name_a: &str,
    names_a: &[impl AsRef<[u8]>],
    name_b: &str,
    names_b: &[impl AsRef<[u8]>],
) -> String {
    let fmt_a = format_chrom_list(names_a);
    let fmt_b = format_chrom_list(names_b);
    let a_bytes: Vec<&[u8]> = names_a.iter().map(|n| n.as_ref()).collect();
    let b_bytes: Vec<&[u8]> = names_b.iter().map(|n| n.as_ref()).collect();
    let suffix = match chr_prefix_hint(&a_bytes, &b_bytes) {
        Some(hint) => format!("; {hint}"),
        None => String::new(),
    };
    format!(
        "{name_a} and {name_b} share no chromosome names ({name_a}: {fmt_a}; {name_b}: {fmt_b}){suffix}"
    )
}

/// For commands over many inputs (`merge`, `complement`): the one warning
/// to print about inputs that share no chromosome names, or `None`. Disjoint
/// chromosome sets are normal there (scatter/gather by chromosome), so only
/// the actionable case warns: two inputs whose names would match with a
/// `chr` prefix added or removed. `names[i]` is input i's name and `seen[i]`
/// the chromosomes its records were on; inputs without records are ignored.
pub fn nary_naming_warning(names: &[String], seen: &[Vec<Vec<u8>>]) -> Option<String> {
    let non_empty: Vec<usize> = (0..names.len().min(seen.len()))
        .filter(|&i| !seen[i].is_empty())
        .collect();
    for (k, &i) in non_empty.iter().enumerate() {
        for &j in &non_empty[k + 1..] {
            let seen_j: HashSet<&Vec<u8>> = seen[j].iter().collect();
            if seen[i].iter().any(|c| seen_j.contains(c)) {
                continue;
            }
            let a: Vec<&[u8]> = seen[i].iter().map(Vec::as_slice).collect();
            let b: Vec<&[u8]> = seen[j].iter().map(Vec::as_slice).collect();
            if chr_prefix_hint(&a, &b).is_some() {
                return Some(format_disjoint_chroms_warning(
                    &names[i], &seen[i], &names[j], &seen[j],
                ));
            }
        }
    }
    None
}

/// Format an unknown-chromosome error message, adding a hint to use the alias
/// genome when `-g` named a built-in genome and its alias contains `chrom`.
pub fn unknown_chrom_err(source: &gloci_core::OrderSource, chrom: &str) -> String {
    let place = match source {
        gloci_core::OrderSource::GenomeFile(p) => format!("the genome file {}", p.display()),
        gloci_core::OrderSource::BuiltinGenome(id) => {
            format!("the built-in genome {id} (gloci genomes show {id} lists its names)")
        }
        gloci_core::OrderSource::VcfHeader(n) => format!("the ##contig lines of {n}"),
        gloci_core::OrderSource::Lexicographic
        | gloci_core::OrderSource::Karyotype
        | gloci_core::OrderSource::FirstAppearance => "the genome".to_string(),
    };
    let hint = match source {
        gloci_core::OrderSource::VcfHeader(_) => format!(
            "; add ##contig=<ID={chrom}> to the header, or pass -g with a genome that lists it"
        ),
        gloci_core::OrderSource::BuiltinGenome(id) => {
            if let Some(alias) = crate::builtin::builtin_alias(id) {
                if let Some(b) = crate::builtin::builtin(alias) {
                    if b.contains_chrom(chrom) {
                        format!("; did you mean -g {alias}?")
                    } else {
                        String::new()
                    }
                } else {
                    String::new()
                }
            } else {
                String::new()
            }
        }
        _ => String::new(),
    };
    format!("chromosome \"{chrom}\" is not in {place}{hint}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_builtin_alias() {
        assert_eq!(builtin_alias("hg38"), Some("GRCh38"));
        assert_eq!(builtin_alias("GRCh38"), Some("hg38"));
        assert_eq!(builtin_alias("HG38"), Some("GRCh38"));
        assert_eq!(builtin_alias("grch38"), Some("hg38"));
        assert_eq!(builtin_alias("mm10"), Some("GRCm38"));
        assert_eq!(builtin_alias("sacCer3"), Some("R64"));
        assert_eq!(builtin_alias("R64"), Some("sacCer3"));
        assert_eq!(builtin_alias("hs1"), None);
        assert_eq!(builtin_alias("foo"), None);
    }

    #[test]
    fn seen_chroms_reports_each_id_once() {
        let mut m = SeenChroms::default();
        assert!(m.first_time(ChromId(3)));
        assert!(!m.first_time(ChromId(3)));
        assert!(m.first_time(ChromId(0)));
        assert!(m.first_time(ChromId(100_000)));
        assert!(!m.first_time(ChromId(0)));
    }

    #[test]
    fn chr_prefix_hint_covers_both_directions_and_bare_chr() {
        assert!(chr_prefix_hint(&[b"chr1"], &[b"chr"]).is_none());
        assert!(chr_prefix_hint(&[b"chr"], &[b"x"]).is_none());
        assert!(chr_prefix_hint(&[b"7", b"8"], &[b"chr9", b"chr8"]).is_some());
    }

    #[test]
    fn test_chr_prefix_hint() {
        assert!(chr_prefix_hint(&[b"chr1"], &[b"1"]).is_some());
        assert!(chr_prefix_hint(&[b"1"], &[b"chr1"]).is_some());
        assert!(chr_prefix_hint(&[b"chr1"], &[b"chr2"]).is_none());
        assert!(chr_prefix_hint(&[b"1"], &[b"2"]).is_none());
        assert!(chr_prefix_hint(&[b"chr1", b"chr2"], &[b"3", b"1"]).is_some());
    }

    #[test]
    fn nary_warns_only_for_the_chr_prefix_case() {
        let names: Vec<String> = (0..24).map(|i| format!("in{i}.bed")).collect();
        let seen: Vec<Vec<Vec<u8>>> = (1..=24)
            .map(|i| vec![format!("chr{i}").into_bytes()])
            .collect();
        assert_eq!(nary_naming_warning(&names, &seen), None);
        let names = vec![
            "a.bed".to_string(),
            "b.bed".to_string(),
            "c.bed".to_string(),
        ];
        let seen = vec![vec![b"chr1".to_vec()], vec![], vec![b"1".to_vec()]];
        let w = nary_naming_warning(&names, &seen).unwrap();
        assert!(
            w.starts_with("a.bed and c.bed share no chromosome names"),
            "{w}"
        );
        assert!(w.contains("\"chr\" prefix"), "{w}");
    }

    #[test]
    fn test_format_chrom_list() {
        assert_eq!(format_chrom_list(&["chr1"]), "chr1");
        assert_eq!(format_chrom_list(&["chr1", "chr2"]), "chr1, chr2");
        assert_eq!(
            format_chrom_list(&["chr1", "chr2", "chr3"]),
            "chr1, chr2, chr3"
        );
        assert_eq!(
            format_chrom_list(&["chr1", "chr2", "chr3", "chr4"]),
            "chr1, chr2, chr3, …"
        );
    }
}
