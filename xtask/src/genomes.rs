//! `cargo xtask genomes refresh`: rebuild `data/genomes/*.genome` from the
//! sources listed in `data/genomes/SOURCES.toml`.

use anyhow::{Context, Result, bail};
use gloci_core::karyotype_cmp;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::{HashMap, HashSet};
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

pub const DEFAULT_DIR: &str = "data/genomes";

/// The Ensembl core-database tables that define the top-level sequences, in
/// the order they are hashed.
const CORE_TABLES: [&str; 4] = [
    "seq_region",
    "seq_region_attrib",
    "attrib_type",
    "coord_system",
];

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum SourceFormat {
    /// A UCSC `*.chrom.sizes` file (or any two-column name/length text).
    ChromSizes,
    /// An Ensembl core database dump directory.
    EnsemblCore,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Entry {
    pub id: String,
    pub format: SourceFormat,
    pub url: String,
    /// `"reference"` keeps the source order; otherwise karyotype order.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub order: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub retrieved: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sha256: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sequences: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub total_length: Option<u64>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Manifest {
    pub genome: Vec<Entry>,
}

const MANIFEST_HEADER: &str = "\
# Provenance of the built-in genomes in this directory.
# `cargo xtask genomes refresh` rewrites the retrieved, sha256, sequences and
# total_length fields. format \"chrom-sizes\": url is the file itself and sha256
# is its hash. format \"ensembl-core\": url is an Ensembl core database dump
# directory; the top-level sequences (the contents of the dna.toplevel FASTA)
# come from seq_region, seq_region_attrib, attrib_type and coord_system (LRG regions
# left out), and sha256 is the hash of those four .txt.gz files concatenated
# in that order.

";

pub struct Args {
    pub dir: PathBuf,
    pub only: Option<String>,
}

/// The reference order of a chromosome named with Roman numerals (yeast):
/// `chrI`..`chrXVI` by number, everything else (mitochondrion, scaffolds)
/// after them in source order. UCSC lists these by size and Ensembl by
/// database id, so the source order is not the reference's.
fn roman_rank(name: &str) -> u32 {
    let n = name.strip_prefix("chr").unwrap_or(name);
    let val = |c| match c {
        'I' => Some(1),
        'V' => Some(5),
        'X' => Some(10),
        'L' => Some(50),
        _ => None,
    };
    let digits: Option<Vec<u32>> = n.chars().map(val).collect();
    let Some(d) = digits.filter(|d| !d.is_empty()) else {
        return u32::MAX;
    };
    let mut total = 0;
    for (i, &v) in d.iter().enumerate() {
        if d.get(i + 1).is_some_and(|&next| next > v) {
            total -= v as i64;
        } else {
            total += v as i64;
        }
    }
    total as u32
}

/// Normalise two-column text: take the first two columns, reject duplicate
/// names, and order with `karyotype_cmp`. With `reference_order` (ids marked
/// `order = "reference"`) the chromosomes keep the reference's own order
/// instead: Roman numerals ascending, then the rest in source order.
pub fn normalize(text: &str, reference_order: bool) -> Result<Vec<(String, u64)>> {
    let mut rows = Vec::new();
    let mut seen = HashSet::new();
    for (i, line) in text.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let mut cols = line.split_whitespace();
        let (Some(name), Some(len)) = (cols.next(), cols.next()) else {
            bail!("line {}: expected a name and a length", i + 1);
        };
        let len: u64 = len
            .parse()
            .with_context(|| format!("line {}: length {len:?} is not an integer", i + 1))?;
        if !seen.insert(name.to_string()) {
            bail!("line {}: duplicate sequence name {name:?}", i + 1);
        }
        rows.push((name.to_string(), len));
    }
    if rows.is_empty() {
        bail!("no sequences");
    }
    if reference_order {
        rows.sort_by_key(|r| roman_rank(&r.0)); // stable
    } else {
        rows.sort_by(|a, b| karyotype_cmp(a.0.as_bytes(), b.0.as_bytes()));
    }
    Ok(rows)
}

/// The top-level sequences of an Ensembl core database as `name<TAB>length`
/// lines, in `seq_region` order. LRG regions are flagged top-level in the
/// database but are not in the `dna.toplevel` FASTA, so they are left out.
pub fn ensembl_toplevel(
    seq_region: &str,
    attrib: &str,
    attrib_type: &str,
    coord_system: &str,
) -> Result<String> {
    let cols = |line: &str| line.split('\t').map(str::to_string).collect::<Vec<_>>();
    let top_id = attrib_type
        .lines()
        .map(cols)
        .find(|c| c.get(1).map(String::as_str) == Some("toplevel"))
        .and_then(|c| c.into_iter().next())
        .context("attrib_type has no \"toplevel\" attribute")?;
    let top: HashSet<String> = attrib
        .lines()
        .map(cols)
        .filter(|c| c.get(1) == Some(&top_id))
        .map(|c| c[0].clone())
        .collect();
    let lrg: HashSet<String> = coord_system
        .lines()
        .map(cols)
        .filter(|c| c.get(2).map(String::as_str) == Some("lrg"))
        .map(|c| c[0].clone())
        .collect();
    let mut out = String::new();
    for line in seq_region.lines() {
        let c = cols(line);
        if c.len() < 4 {
            bail!("malformed seq_region row {line:?}");
        }
        if top.contains(&c[0]) && !lrg.contains(&c[2]) {
            out.push_str(&format!("{}\t{}\n", c[1], c[3]));
        }
    }
    Ok(out)
}

fn fetch(url: &str) -> Result<Vec<u8>> {
    let out = Command::new("curl")
        .args(["-fsSL", "--retry", "3", url])
        .output()
        .context("running curl (is it installed?)")?;
    if !out.status.success() {
        bail!(
            "downloading {url} failed: {}",
            String::from_utf8_lossy(&out.stderr).trim()
        );
    }
    Ok(out.stdout)
}

fn gunzip(bytes: &[u8]) -> Result<String> {
    let mut s = String::new();
    flate2::read::MultiGzDecoder::new(bytes).read_to_string(&mut s)?;
    Ok(s)
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

/// Download one entry: the raw two-column text and the hash of the source.
fn download(e: &Entry) -> Result<(String, String)> {
    match e.format {
        SourceFormat::ChromSizes => {
            let b = fetch(&e.url)?;
            let text = String::from_utf8(b.clone()).context("source is not UTF-8")?;
            Ok((text, hex(&Sha256::digest(&b))))
        }
        SourceFormat::EnsemblCore => {
            let mut h = Sha256::new();
            let mut t = Vec::new();
            for table in CORE_TABLES {
                let url = format!("{}/{table}.txt.gz", e.url.trim_end_matches('/'));
                let b = fetch(&url)?;
                h.update(&b);
                t.push(gunzip(&b).with_context(|| url.clone())?);
            }
            Ok((
                ensembl_toplevel(&t[0], &t[1], &t[2], &t[3])?,
                hex(&h.finalize()),
            ))
        }
    }
}

/// Days since 1970-01-01 to a civil `YYYY-MM-DD` date.
pub fn civil_date(days: i64) -> String {
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    format!("{y:04}-{m:02}-{d:02}")
}

fn today() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_secs());
    civil_date((secs / 86_400) as i64)
}

pub fn render(rows: &[(String, u64)]) -> String {
    let mut s = String::new();
    for (n, l) in rows {
        s.push_str(&format!("{n}\t{l}\n"));
    }
    s
}

pub fn run(args: &Args) -> Result<bool> {
    let path = args.dir.join("SOURCES.toml");
    let text = std::fs::read_to_string(&path).with_context(|| path.display().to_string())?;
    let mut manifest: Manifest =
        toml::from_str(&text).with_context(|| path.display().to_string())?;
    if let Some(id) = &args.only
        && !manifest.genome.iter().any(|e| &e.id == id)
    {
        bail!("no genome {id:?} in {}", path.display());
    }
    let date = today();
    for i in 0..manifest.genome.len() {
        let e = &mut manifest.genome[i];
        if args.only.as_ref().is_some_and(|id| id != &e.id) {
            continue;
        }
        eprintln!("{}: {}", e.id, e.url);
        let (raw, sha) = download(e)?;
        let rows = normalize(&raw, e.order.as_deref() == Some("reference"))
            .with_context(|| e.id.clone())?;
        write_atomic(&args.dir.join(format!("{}.genome", e.id)), &render(&rows))?;
        e.retrieved = Some(date.clone());
        e.sha256 = Some(sha);
        e.sequences = Some(rows.len() as u64);
        e.total_length = Some(rows.iter().map(|r| r.1).sum());
        // Save after every entry so a failure keeps the finished ones.
        write_manifest(&path, &manifest)?;
    }
    Ok(true)
}

fn write_manifest(path: &Path, m: &Manifest) -> Result<()> {
    let body = toml::to_string(m)?;
    write_atomic(path, &format!("{MANIFEST_HEADER}{body}"))
}

fn write_atomic(path: &Path, text: &str) -> Result<()> {
    let tmp = path.with_extension("tmp");
    std::fs::write(&tmp, text)?;
    std::fs::rename(&tmp, path)?;
    Ok(())
}

/// Used by tests: map id to its entry.
pub fn by_id(m: &Manifest) -> HashMap<&str, &Entry> {
    m.genome.iter().map(|e| (e.id.as_str(), e)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(rows: &[(String, u64)]) -> Vec<&str> {
        rows.iter().map(|r| r.0.as_str()).collect()
    }

    #[test]
    fn karyotype_sorting_and_extra_columns() {
        let t = "chrX\t5\tx\nchr10\t4\nchrM\t3\nchr2\t2\nchrUn_1\t1\nchr1\t9\n";
        let rows = normalize(t, false).unwrap();
        assert_eq!(
            names(&rows),
            ["chr1", "chr2", "chr10", "chrX", "chrM", "chrUn_1"]
        );
        assert_eq!(rows[0].1, 9);
    }

    #[test]
    fn reference_order_is_kept() {
        // Roman numerals ascend (not "chrI, chrII, chrIII, chrIV, chrIX, chrV"
        // as karyotype order would give, nor the size order of the source);
        // other names follow in source order.
        let t = "chrM\t9\nchrIX\t5\nchrIV\t4\nchrXVI\t6\nchrI\t1\nplasmid\t2\nchrV\t3\n";
        let rows = normalize(t, true).unwrap();
        assert_eq!(
            names(&rows),
            [
                "chrI", "chrIV", "chrV", "chrIX", "chrXVI", "chrM", "plasmid"
            ]
        );
    }

    #[test]
    fn duplicates_and_bad_lengths_are_rejected() {
        assert!(
            normalize("a\t1\na\t2\n", false)
                .unwrap_err()
                .to_string()
                .contains("duplicate")
        );
        assert!(normalize("a\tx\n", false).is_err());
        assert!(normalize("a\n", false).is_err());
        assert!(normalize("", false).is_err());
    }

    #[test]
    fn ensembl_toplevel_selects_flagged_regions() {
        let sr = "1\tctg1\t1\t10\n2\tI\t2\t100\n";
        let at = "2\t6\tx\n3\t6\t\n3\t9\t1\n";
        let ty = "1\tembl_acc\tx\n6\ttoplevel\tTop Level\n";
        let cs = "1\t1\tcontig\t\\N\t3\n2\t1\tchromosome\tX\t1\n3\t1\tlrg\t\\N\t5\n";
        assert_eq!(ensembl_toplevel(sr, at, ty, cs).unwrap(), "I\t100\n");
        // An LRG region flagged top-level is left out.
        let sr = "2\tI\t2\t100\n3\tLRG_1\t3\t5\n";
        assert_eq!(ensembl_toplevel(sr, at, ty, cs).unwrap(), "I\t100\n");
    }

    #[test]
    fn civil_dates() {
        assert_eq!(civil_date(0), "1970-01-01");
        assert_eq!(civil_date(20_425), "2025-12-03");
        assert_eq!(civil_date(11_016), "2000-02-29");
    }
}
