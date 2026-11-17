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
    // The guard line in the Dockerfile itself greps for "ucsc", so look only at non-comment words
    // that name a package: no package or tool line may mention it.
    let lower = text.to_lowercase();
    let hits: Vec<&str> = lower
        .lines()
        .filter(|l| {
            l.contains("ucsc") && !l.trim_start().starts_with('#') && !l.contains("grep -qi ucsc")
        })
        .collect();
    assert!(hits.is_empty(), "UCSC named in the Dockerfile: {hits:?}");
}
