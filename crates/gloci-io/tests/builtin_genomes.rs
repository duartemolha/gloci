use gloci_core::{OrderSource, karyotype_cmp};
use gloci_io::{builtin, builtin_genomes, resolve_genome};
use gloci_testkit::fixtures;

const IDS: [&str; 19] = [
    "hg19", "hg38", "hs1", "GRCh37", "GRCh38", "mm10", "mm39", "GRCm38", "GRCm39", "rn7",
    "danRer11", "GRCz11", "dm6", "BDGP6", "ce11", "WBcel235", "sacCer3", "R64", "TAIR10",
];

fn len(id: &str, chrom: &str) -> u64 {
    let g = builtin(id).unwrap().genome();
    (0..g.len())
        .find(|&i| g.name(i) == chrom.as_bytes())
        .and_then(|i| g.length(i))
        .unwrap()
}

#[test]
fn known_lengths() {
    for (id, c, l) in [
        ("hg38", "chr1", 248956422),
        ("hg38", "chrX", 156040895),
        ("hg38", "chrM", 16569),
        ("hg38", "chrEBV", 171823),
        ("hg19", "chr1", 249250621),
        ("hg19", "chrM", 16571),
        ("GRCh37", "1", 249250621),
        ("GRCh37", "MT", 16569),
        ("GRCh38", "1", 248956422),
        ("hs1", "chr1", 248387328),
        ("mm10", "chr1", 195471971),
        ("mm39", "chr1", 195154279),
        ("GRCm39", "1", 195154279),
        ("dm6", "chr2L", 23513712),
        ("BDGP6", "2L", 23513712),
        ("ce11", "chrI", 15072434),
        ("WBcel235", "I", 15072434),
        ("sacCer3", "chrI", 230218),
        ("R64", "I", 230218),
        ("danRer11", "chr1", 59578282),
        ("TAIR10", "1", 30427671),
    ] {
        assert_eq!(len(id, c), l, "{id} {c}");
    }
}

#[test]
fn every_id_present_and_case_insensitive() {
    for id in IDS {
        assert!(
            builtin(id).is_some() && builtin(&id.to_uppercase()).is_some(),
            "{id}"
        );
    }
    assert_eq!(builtin_genomes().len(), IDS.len());
}

/// The ids marked `order = "reference"` in SOURCES.toml: the one place the
/// Roman-numeral exception is listed.
fn reference_ordered() -> Vec<String> {
    let manifest = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/data/genomes/SOURCES.toml"
    ))
    .unwrap();
    let t: toml::Table = manifest.parse().unwrap();
    t["genome"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|e| e.get("order").and_then(|o| o.as_str()) == Some("reference"))
        .map(|e| e["id"].as_str().unwrap().to_string())
        .collect()
}

fn names_of(id: &str) -> Vec<String> {
    let g = builtin(id).unwrap().genome();
    (0..g.len())
        .map(|i| String::from_utf8(g.name(i).to_vec()).unwrap())
        .collect()
}

#[test]
fn roman_numeral_genomes_are_in_reference_order() {
    let mut r = reference_ordered();
    r.sort();
    assert_eq!(r, ["R64", "WBcel235", "ce11", "sacCer3"]);
    assert_eq!(
        names_of("ce11"),
        ["chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX", "chrM"]
    );
    assert_eq!(
        names_of("WBcel235"),
        ["I", "II", "III", "IV", "V", "X", "MtDNA"]
    );
    assert_eq!(&names_of("R64")[..5], ["I", "II", "III", "IV", "V"]);
    assert_eq!(names_of("R64")[15..], ["XVI", "Mito"]);
}

#[test]
fn ensembl_tables_are_primary_assemblies() {
    let (human, mouse) = (3_000_000_000..=3_300_000_000, 2_600_000_000..=2_800_000_000);
    for (id, range) in [
        ("GRCh38", &human),
        ("GRCh37", &human),
        ("GRCm38", &mouse),
        ("GRCm39", &mouse),
    ] {
        let t = builtin(id).unwrap().total_length();
        assert!(range.contains(&t), "{id}: {t}");
    }
    for b in builtin_genomes()
        .iter()
        .filter(|b| b.naming == gloci_io::Naming::Ensembl)
    {
        for n in names_of(b.id) {
            assert!(
                !["PATCH", "HSCHR", "CHR_", "LRG"]
                    .iter()
                    .any(|p| n.contains(p)),
                "{}: {n}",
                b.id
            );
        }
        // No N-padded chromosome-length copy: nothing is longer than chr1.
        let g = b.genome();
        let max = (0..g.len()).filter_map(|i| g.length(i)).max().unwrap();
        assert!(max < 400_000_000, "{}: {max}", b.id);
    }
}

#[test]
fn orders() {
    let g = builtin("hg38").unwrap().genome();
    let first: Vec<_> = (0..25)
        .map(|i| String::from_utf8(g.name(i).to_vec()).unwrap())
        .collect();
    assert_eq!(
        first,
        (1..=22)
            .map(|i| format!("chr{i}"))
            .chain(["chrX", "chrY", "chrM"].map(String::from))
            .collect::<Vec<_>>()
    );
    let y = builtin("sacCer3").unwrap().genome();
    assert_eq!(y.name(3), b"chrIV");
    assert_eq!(y.name(8), b"chrIX"); // reference order kept
    let reference = reference_ordered();
    for b in builtin_genomes()
        .iter()
        .filter(|b| !reference.iter().any(|r| r == b.id))
    {
        let g = b.genome();
        let names: Vec<_> = (0..g.len()).map(|i| g.name(i).to_vec()).collect();
        let mut s = names.clone();
        s.sort_by(|a, b| karyotype_cmp(a, b));
        assert_eq!(names, s, "{}", b.id);
    }
}

#[test]
fn resolution() {
    assert!(matches!(
        resolve_genome("HG38").unwrap().1,
        OrderSource::BuiltinGenome(ref id) if id == "hg38"
    ));
    let e = resolve_genome("hg37").unwrap_err().to_string();
    assert!(e.contains("did you mean") && (e.contains("hg38") || e.contains("hg19")));
    assert!(e.contains("built-in genomes: hg19, hg38,"), "{e}");
    assert!(e.contains("gloci genomes list"), "{e}");
    assert!(
        e.contains(
            "-g takes a built-in genome or an existing .fai, .genome, .dict, VCF/BCF or BAM file"
        ),
        "{e}"
    );
    assert!(resolve_genome("./hg38").is_err()); // path form, no such file
    let d = tempfile::tempdir().unwrap();
    let f = fixtures::write(d.path(), "hg38", b"chrA\t10\n");
    assert!(matches!(
        resolve_genome(f.to_str().unwrap()).unwrap().1,
        OrderSource::GenomeFile(_)
    )); // files win
    // A far-off name gets no suggestion.
    let e = resolve_genome("zzzzzz").unwrap_err().to_string();
    assert!(!e.contains("did you mean"), "{e}");
}

#[test]
fn totals_and_counts() {
    let b = builtin("hg38").unwrap();
    assert_eq!(b.len(), 456);
    assert_eq!(b.total_length(), 3209457928);
    assert_eq!(b.genome().len(), b.len());
}

#[test]
fn sources_manifest_matches_data() {
    let manifest = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/data/genomes/SOURCES.toml"
    ))
    .unwrap();
    let t: toml::Table = manifest.parse().unwrap();
    let entries = t["genome"].as_array().unwrap();
    assert_eq!(entries.len(), IDS.len());
    for id in IDS {
        let e = entries
            .iter()
            .find(|e| e["id"].as_str() == Some(id))
            .unwrap_or_else(|| panic!("{id} missing from SOURCES.toml"));
        let url = e["url"].as_str().unwrap();
        assert!(url.starts_with("https://"), "{id}");
        let date = e["retrieved"].as_str().unwrap();
        assert!(
            date.len() == 10 && date.as_bytes()[4] == b'-',
            "{id}: {date}"
        );
        let sha = e["sha256"].as_str().unwrap();
        assert!(
            sha.len() == 64 && sha.bytes().all(|b| b.is_ascii_hexdigit()),
            "{id}"
        );
        let b = builtin(id).unwrap();
        assert_eq!(
            e["sequences"].as_integer().unwrap() as usize,
            b.len(),
            "{id}"
        );
        assert_eq!(
            e["total_length"].as_integer().unwrap() as u64,
            b.total_length(),
            "{id}"
        );
        let file = std::fs::read_to_string(format!(
            "{}/data/genomes/{id}.genome",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        assert_eq!(file.lines().count(), b.len(), "{id}");
    }
}
