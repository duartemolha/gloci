//! The comparison image's Dockerfile: CrossMap comes from PyPI without its dependencies,
//! because every bioconda CrossMap build depends on a UCSC package that the image must not hold.
use std::fs;

#[test]
fn diff_image_pins_crossmap_without_ucsc() {
    let path = xtask::diff::workspace_root().join("xtask/docker/Dockerfile");
    let text = fs::read_to_string(&path).unwrap();
    for needle in [
        "crossmap==0.7.4",
        "--no-deps",
        "--require-hashes",
        "2c976861eb0d64e3afd8cc03413ffe98c496e8fd4e9d60c7913593cff4e62648",
        "pysam=0.23.3",
        "bx-python=0.16.0",
        "pybigwig=",
    ] {
        assert!(text.contains(needle), "Dockerfile lacks `{needle}`");
    }
    // The build-time guard line (`micromamba list ... | grep -qi ucsc`) must be there: it also
    // catches transitive dependencies, which this text check cannot see. Only that exact line may
    // name the word; no other non-comment line may.
    const GUARD: &str = "&& ! micromamba list -n base | grep -qi ucsc";
    assert!(
        text.lines().any(|l| l.trim() == GUARD),
        "Dockerfile lacks its build-time guard line `{GUARD}`"
    );
    let lower = text.to_lowercase();
    let hits: Vec<&str> = lower
        .lines()
        .filter(|l| l.contains("ucsc") && !l.trim_start().starts_with('#') && l.trim() != GUARD)
        .collect();
    assert!(hits.is_empty(), "UCSC named in the Dockerfile: {hits:?}");
}
