//! `cargo xtask genomes refresh`: rebuild `crates/gloci-io/data/genomes/*.genome` from the
//! sources listed in `crates/gloci-io/data/genomes/SOURCES.toml`.

use anyhow::{Context, Result, bail};
use gloci_core::karyotype_cmp;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::{HashMap, HashSet};
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

pub const DEFAULT_DIR: &str = "crates/gloci-io/data/genomes";

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
    /// Sequences the source lacks, added on request.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub extra: Vec<Extra>,
}

/// A sequence a genome's source lacks, added on request, with where it comes from.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Extra {
    pub name: String,
    pub length: u64,
    pub source: String,
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
# directory; the primary assembly (chromosomes plus unlocalized and unplaced
# contigs; no haplotypes or patches, which carry the non_ref attribute) comes
# from seq_region, seq_region_attrib, attrib_type and coord_system (LRG regions
# left out), and sha256 is the hash of those four .txt.gz files concatenated
# in that order. extra lists sequences the source lacks, added on request, each
# with its source: a refresh adds them in their karyotype place (sequences and
# total_length count them; sha256 stays the source's), skips one the source has
# with the same length and refuses one it has with another length.

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

/// The primary assembly of an Ensembl core database as `name<TAB>length`
/// lines, in `seq_region` order: the top-level sequences (chromosomes,
/// unlocalized and unplaced contigs) minus the regions that carry the
/// `non_ref` attribute (haplotypes and patches, which Ensembl stores as
/// N-padded copies the length of a chromosome) and minus LRG regions. The
/// dumps' `assembly_exception` table is empty, so `non_ref` is the mechanism.
pub fn ensembl_primary(
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
    let attrib_code = |code: &str| {
        attrib_type
            .lines()
            .map(cols)
            .find(|c| c.get(1).map(String::as_str) == Some(code))
            .and_then(|c| c.into_iter().next())
    };
    let regions_with = |id: &Option<String>| -> HashSet<String> {
        attrib
            .lines()
            .map(cols)
            .filter(|c| c.get(1) == id.as_ref())
            .map(|c| c[0].clone())
            .collect()
    };
    let top = regions_with(&Some(top_id));
    let non_ref = regions_with(&attrib_code("non_ref"));
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
        if top.contains(&c[0]) && !non_ref.contains(&c[0]) && !lrg.contains(&c[2]) {
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
                ensembl_primary(&t[0], &t[1], &t[2], &t[3])?,
                hex(&h.finalize()),
            ))
        }
    }
}

/// The two-column text `raw` with the sequences of `extra` it lacks appended: one it has with
/// the same length is skipped, one it has with another length is an error.
pub fn add_extras(raw: &str, extra: &[Extra]) -> Result<String> {
    let have: HashMap<&str, &str> = raw
        .lines()
        .filter_map(|l| {
            let mut c = l.split_whitespace();
            Some((c.next()?, c.next()?))
        })
        .collect();
    let mut out = raw.to_string();
    if !out.is_empty() && !out.ends_with('\n') {
        out.push('\n');
    }
    for x in extra {
        match have.get(x.name.as_str()) {
            None => out.push_str(&format!("{}\t{}\n", x.name, x.length)),
            Some(len) if *len == x.length.to_string() => {}
            Some(len) => bail!(
                "extra sequence {} is {} long, but the source has it at {len}",
                x.name,
                x.length
            ),
        }
    }
    Ok(out)
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
        let raw = add_extras(&raw, &e.extra).with_context(|| e.id.clone())?;
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
    fn ensembl_primary_selects_flagged_regions() {
        let sr = "1\tctg1\t1\t10\n2\tI\t2\t100\n";
        let at = "2\t6\tx\n3\t6\t\n3\t9\t1\n";
        let ty = "1\tembl_acc\tx\n6\ttoplevel\tTop Level\n";
        let cs = "1\t1\tcontig\t\\N\t3\n2\t1\tchromosome\tX\t1\n3\t1\tlrg\t\\N\t5\n";
        assert_eq!(ensembl_primary(sr, at, ty, cs).unwrap(), "I\t100\n");
        // An LRG region flagged top-level is left out.
        let sr = "2\tI\t2\t100\n3\tLRG_1\t3\t5\n";
        assert_eq!(ensembl_primary(sr, at, ty, cs).unwrap(), "I\t100\n");
        // A non_ref (haplotype or patch) region is left out.
        let sr = "2\tI\t2\t100\n4\tHG1_PATCH\t2\t99\n";
        let at = "2\t6\tx\n4\t6\t\n4\t16\t1\n";
        let ty = "6\ttoplevel\tTop Level\n16\tnon_ref\tNon Reference\n";
        assert_eq!(ensembl_primary(sr, at, ty, cs).unwrap(), "I\t100\n");
    }

    fn extra(name: &str, length: u64) -> Extra {
        Extra {
            name: name.into(),
            length,
            source: "test".into(),
        }
    }

    #[test]
    fn extras_join_in_karyotype_order() {
        let raw = "chrUn_1\t1\nchr2\t2\nchrM\t3\nchr1\t9";
        let rows = normalize(&add_extras(raw, &[extra("chrEBV", 5)]).unwrap(), false).unwrap();
        assert_eq!(names(&rows), ["chr1", "chr2", "chrM", "chrEBV", "chrUn_1"]);
        assert_eq!(rows[3].1, 5);
    }

    #[test]
    fn an_extra_the_source_has_is_kept_once() {
        let raw = "chr1\t9\nchrEBV\t5\n";
        assert_eq!(add_extras(raw, &[extra("chrEBV", 5)]).unwrap(), raw);
    }

    #[test]
    fn an_extra_with_another_length_is_an_error() {
        let e = add_extras("chr1\t9\nchrEBV\t6\n", &[extra("chrEBV", 5)]).unwrap_err();
        assert!(e.to_string().contains("chrEBV is 5 long"), "{e}");
    }

    fn data_dir() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("..")
            .join(DEFAULT_DIR)
    }

    /// The committed manifest is in the form a refresh writes, so its fields (the extras too)
    /// survive a refresh.
    #[test]
    fn the_manifest_round_trips() {
        let path = data_dir().join("SOURCES.toml");
        let text = std::fs::read_to_string(&path).unwrap();
        let m: Manifest = toml::from_str(&text).unwrap();
        assert_eq!(
            format!("{MANIFEST_HEADER}{}", toml::to_string(&m).unwrap()),
            text
        );
    }

    /// hg38.genome is UCSC's table (as committed, without the extras) plus SOURCES.toml's
    /// extras, normalized as a refresh would, and the manifest counts them.
    #[test]
    fn hg38_is_ucscs_table_plus_its_extras() {
        let dir = data_dir();
        let m: Manifest =
            toml::from_str(&std::fs::read_to_string(dir.join("SOURCES.toml")).unwrap()).unwrap();
        let e = by_id(&m)["hg38"];
        assert_eq!(names_of_extras(e), ["chrEBV"]);
        let committed = std::fs::read_to_string(dir.join("hg38.genome")).unwrap();
        let ucsc: String = committed
            .lines()
            .filter(|l| !l.starts_with("chrEBV\t"))
            .map(|l| format!("{l}\n"))
            .collect();
        let rows = normalize(&add_extras(&ucsc, &e.extra).unwrap(), false).unwrap();
        assert_eq!(render(&rows), committed);
        assert_eq!(e.sequences, Some(rows.len() as u64));
        assert_eq!(e.total_length, Some(rows.iter().map(|r| r.1).sum()));
    }

    fn names_of_extras(e: &Entry) -> Vec<&str> {
        e.extra.iter().map(|x| x.name.as_str()).collect()
    }

    #[test]
    fn civil_dates() {
        assert_eq!(civil_date(0), "1970-01-01");
        assert_eq!(civil_date(20_425), "2025-12-03");
        assert_eq!(civil_date(11_016), "2000-02-29");
    }
}
