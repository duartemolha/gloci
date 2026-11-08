//! Genomes built into the binary (`-g hg38`), and how a `-g` value is resolved.

use crate::error::{IoError, UsageError};
use crate::genome::load_genome;
use gloci_core::{Genome, OrderSource};
use std::path::Path;

/// How an assembly names its chromosomes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Naming {
    /// `chr1`, `chrM`.
    Ucsc,
    /// `1`, `MT`.
    Ensembl,
}

impl Naming {
    pub fn as_str(self) -> &'static str {
        match self {
            Naming::Ucsc => "ucsc",
            Naming::Ensembl => "ensembl",
        }
    }
}

/// An assembly whose chromosome names and lengths are compiled in.
#[derive(Debug)]
pub struct BuiltinGenome {
    pub id: &'static str,
    pub species: &'static str,
    pub assembly: &'static str,
    pub naming: Naming,
    data: &'static str,
}

impl BuiltinGenome {
    fn rows(&self) -> impl Iterator<Item = (&'static str, u64)> {
        self.data.lines().filter(|l| !l.is_empty()).map(|l| {
            let (n, len) = l
                .split_once('\t')
                .expect("built-in genome rows have two columns");
            (
                n,
                len.parse().expect("built-in genome lengths are integers"),
            )
        })
    }

    /// The chromosomes in their built-in order.
    pub fn genome(&self) -> Genome {
        Genome::parse_sizes(self.data.as_bytes()).expect("built-in genome data is valid")
    }

    /// The number of sequences.
    pub fn len(&self) -> usize {
        self.rows().count()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// The sum of all sequence lengths.
    pub fn total_length(&self) -> u64 {
        self.rows().map(|r| r.1).sum()
    }
}

macro_rules! genome {
    ($id:literal, $species:literal, $assembly:literal, $naming:ident) => {
        BuiltinGenome {
            id: $id,
            species: $species,
            assembly: $assembly,
            naming: Naming::$naming,
            data: include_str!(concat!("../../../data/genomes/", $id, ".genome")),
        }
    };
}

static GENOMES: [BuiltinGenome; 19] = [
    genome!("hg19", "Homo sapiens", "GRCh37 (UCSC hg19)", Ucsc),
    genome!("hg38", "Homo sapiens", "GRCh38 (UCSC hg38)", Ucsc),
    genome!("hs1", "Homo sapiens", "T2T-CHM13v2.0 (UCSC hs1)", Ucsc),
    genome!("GRCh37", "Homo sapiens", "GRCh37 (Ensembl)", Ensembl),
    genome!("GRCh38", "Homo sapiens", "GRCh38 (Ensembl)", Ensembl),
    genome!("mm10", "Mus musculus", "GRCm38 (UCSC mm10)", Ucsc),
    genome!("mm39", "Mus musculus", "GRCm39 (UCSC mm39)", Ucsc),
    genome!("GRCm38", "Mus musculus", "GRCm38 (Ensembl)", Ensembl),
    genome!("GRCm39", "Mus musculus", "GRCm39 (Ensembl)", Ensembl),
    genome!("rn7", "Rattus norvegicus", "mRatBN7.2 (UCSC rn7)", Ucsc),
    genome!("danRer11", "Danio rerio", "GRCz11 (UCSC danRer11)", Ucsc),
    genome!("GRCz11", "Danio rerio", "GRCz11 (Ensembl)", Ensembl),
    genome!("dm6", "Drosophila melanogaster", "BDGP6 (UCSC dm6)", Ucsc),
    genome!(
        "BDGP6",
        "Drosophila melanogaster",
        "BDGP6 (Ensembl)",
        Ensembl
    ),
    genome!(
        "ce11",
        "Caenorhabditis elegans",
        "WBcel235 (UCSC ce11)",
        Ucsc
    ),
    genome!(
        "WBcel235",
        "Caenorhabditis elegans",
        "WBcel235 (Ensembl)",
        Ensembl
    ),
    genome!(
        "sacCer3",
        "Saccharomyces cerevisiae",
        "R64-1-1 (UCSC sacCer3)",
        Ucsc
    ),
    genome!(
        "R64",
        "Saccharomyces cerevisiae",
        "R64-1-1 (Ensembl)",
        Ensembl
    ),
    genome!(
        "TAIR10",
        "Arabidopsis thaliana",
        "TAIR10 (Ensembl Plants)",
        Ensembl
    ),
];

/// Every built-in genome, in listing order.
pub fn builtin_genomes() -> &'static [BuiltinGenome] {
    &GENOMES
}

/// Look a built-in genome up by id, ignoring case.
pub fn builtin(id: &str) -> Option<&'static BuiltinGenome> {
    GENOMES.iter().find(|g| g.id.eq_ignore_ascii_case(id))
}

fn edit_distance(a: &str, b: &str) -> usize {
    let b: Vec<char> = b.chars().collect();
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    for (i, ca) in a.chars().enumerate() {
        let mut cur = vec![i + 1];
        for (j, cb) in b.iter().enumerate() {
            let sub = prev[j] + usize::from(ca != *cb);
            cur.push(sub.min(prev[j + 1] + 1).min(cur[j] + 1));
        }
        prev = cur;
    }
    prev[b.len()]
}

/// The built-in id closest to `name` (edit distance 2 or less, ignoring case).
fn closest(name: &str) -> Option<&'static str> {
    let name = name.to_lowercase();
    GENOMES
        .iter()
        .map(|g| (edit_distance(&name, &g.id.to_lowercase()), g.id))
        .filter(|(d, _)| *d <= 2)
        .min_by_key(|(d, _)| *d)
        .map(|(_, id)| id)
}

/// The usage error for a `-g` value that is neither a file nor a built-in.
pub fn unknown_genome(arg: &str) -> IoError {
    let ids: Vec<&str> = GENOMES.iter().map(|g| g.id).collect();
    let hint = closest(arg)
        .map(|c| format!(" (did you mean {c}?)"))
        .unwrap_or_default();
    IoError::Usage(UsageError(format!(
        "unknown genome \"{arg}\"; built-in genomes: {}{hint}; or pass a .fai/.genome/.dict file",
        ids.join(", ")
    )))
}

/// Does a `-g` value look like a path rather than a genome id: it has a path
/// separator, or ends in the extension of a genome file?
fn looks_like_path(arg: &str) -> bool {
    const EXTS: [&str; 8] = [
        "fai", "genome", "dict", "sizes", "txt", "bam", "cram", "sam",
    ];
    arg.contains(std::path::MAIN_SEPARATOR)
        || arg.contains('/')
        || Path::new(arg)
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| EXTS.iter().any(|x| e.eq_ignore_ascii_case(x)))
}

/// Resolve a `-g` value. A value with a `/` in it, or one that names an
/// existing file, is a file (a directory does not count); otherwise it is a built-in genome id.
pub fn resolve_genome(arg: &str) -> Result<(Genome, OrderSource), IoError> {
    let path = Path::new(arg);
    let built = builtin(arg);
    if looks_like_path(arg) && !path.exists() {
        return Err(IoError::Usage(UsageError(format!(
            "no such file: {arg}; -g takes a built-in genome (gloci genomes lists them) \
or an existing .fai, .genome, .dict, VCF/BCF or BAM file"
        ))));
    }
    if arg.contains('/') || path.is_file() {
        if built.is_some() {
            log::warn!(
                "-g {arg} names a file in the current directory, so it is used instead of the built-in genome {arg}; to get the built-in genome, rename or move that file; to use the file, pass its path, e.g. -g ./{arg}"
            );
        }
        let g = load_genome(path)?;
        return Ok((g, OrderSource::GenomeFile(path.to_path_buf())));
    }
    match built {
        Some(b) => Ok((b.genome(), OrderSource::BuiltinGenome(b.id.to_string()))),
        None => Err(unknown_genome(arg)),
    }
}
