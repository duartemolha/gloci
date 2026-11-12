//! Bioconda recipes (D93): they render with the release's version and checksums, the CLI
//! recipe has no run requirements, the library recipe installs no CLI, and every file
//! the recipes test with exists.
use std::fs;
use std::path::{Path, PathBuf};
use xtask::conda::{render, render_all};

const SHA_CLI: &str = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
const SHA_PY: &str = "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb";

fn recipes() -> PathBuf {
    xtask::diff::workspace_root().join("packaging/bioconda")
}

fn cli() -> String {
    render(&recipes().join("gloci/meta.yaml"), "7.8.9", SHA_CLI, None).unwrap()
}

fn py() -> String {
    render(&recipes().join("py-gloci/meta.yaml"), "7.8.9", SHA_PY, None).unwrap()
}

/// The lines of the top-level `key:` block (up to the next top-level key).
fn block<'a>(yml: &'a str, key: &str) -> Vec<&'a str> {
    let mut it = yml.lines().skip_while(|l| *l != format!("{key}:"));
    it.next();
    it.take_while(|l| l.is_empty() || l.starts_with(' ') || l.starts_with('#'))
        .collect()
}

#[test]
fn recipes_render_with_versions() {
    for (name, text, sha) in [("gloci", cli(), SHA_CLI), ("py-gloci", py(), SHA_PY)] {
        assert!(!text.contains('@'), "{name}: placeholder left in\n{text}");
        assert!(
            text.contains("{% set version = \"7.8.9\" %}"),
            "{name}\n{text}"
        );
        assert!(text.contains(&format!("sha256: {sha}")), "{name}\n{text}");
        assert!(
            !text.contains("0.1.0"),
            "{name} keeps the committed version\n{text}"
        );
    }
    assert!(cli().contains("archive/refs/tags/v{{ version }}.tar.gz"));
    // The library recipe builds the sdist attached to the GitHub Release (PyPI is off by default).
    assert!(py().contains("releases/download/v{{ version }}/gloci-{{ version }}.tar.gz"));
    assert!(render(&recipes().join("gloci/meta.yaml"), "7.8.9", "xyz", None).is_err());
    assert!(render(&recipes().join("gloci/meta.yaml"), "latest", SHA_CLI, None).is_err());
}

#[test]
fn local_render_swaps_url_for_path() {
    let src = Path::new("/tmp/some/src");
    let t = render(&recipes().join("gloci/meta.yaml"), "7.8.9", "", Some(src)).unwrap();
    assert!(t.contains("  path: /tmp/some/src\n"), "{t}");
    assert!(!t.contains("\n  url:") && !t.contains("sha256:"), "{t}");
    assert!(!t.contains('@'), "{t}");
}

#[test]
fn render_all_writes_both_recipes_and_their_files() {
    let out = tempfile::tempdir().unwrap();
    render_all(
        &recipes(),
        out.path(),
        "7.8.9",
        Some(SHA_CLI),
        Some(SHA_PY),
        None,
    )
    .unwrap();
    for f in [
        "gloci/meta.yaml",
        "gloci/build.sh",
        "gloci/test/a.bed",
        "gloci/test/b.bed",
        "gloci/test/expected.bed",
        "py-gloci/meta.yaml",
    ] {
        assert!(out.path().join(f).is_file(), "{f} missing");
    }
    let m = fs::read_to_string(out.path().join("gloci/meta.yaml")).unwrap();
    assert!(m.contains("version = \"7.8.9\"") && m.contains(SHA_CLI));
    // Without --local both checksums are required.
    let e = render_all(&recipes(), out.path(), "7.8.9", None, Some(SHA_PY), None).unwrap_err();
    assert!(e.to_string().contains("sha256"), "{e}");
}

#[test]
fn cli_recipe_has_no_run_requirements() {
    let t = cli();
    let req = block(&t, "requirements");
    assert!(req.iter().any(|l| l.trim() == "build:"), "{req:?}");
    for l in &req {
        assert_ne!(
            l.trim(),
            "run:",
            "the static binary needs nothing at run time"
        );
        assert!(!l.contains("host:"), "{l}");
    }
    let build = req.join("\n");
    for need in [
        "{{ compiler('c') }}",
        "{{ compiler('rust') }}",
        "cargo-bundle-licenses",
    ] {
        assert!(build.contains(need), "build requirement {need}\n{build}");
    }
    assert!(t.contains("run_exports:"));
    assert!(t.contains("pin_subpackage('gloci', max_pin=\"x.x\")"));
    assert!(t.contains("license_file:") && t.contains("THIRDPARTY.yml"));
    assert!(t.contains("recipe-maintainers:") && t.contains("- duartemolha"));
    assert!(t.contains("linux-aarch64") && t.contains("osx-arm64"));
    assert!(t.contains("license: MIT"));
    for cmd in [
        "gloci --version",
        "gloci genomes show hg38",
        "gloci intersect -a a.bed -b b.bed | cmp - expected.bed",
        "gloci.1",
        "bash-completion/completions/gloci",
        "site-functions/_gloci",
        "vendor_completions.d/gloci.fish",
    ] {
        assert!(t.contains(cmd), "test command `{cmd}` missing\n{t}");
    }
}

#[test]
fn cli_build_script_follows_the_brief() {
    let b = fs::read_to_string(recipes().join("gloci/build.sh")).unwrap();
    let at = |needle: &str| {
        b.find(needle)
            .unwrap_or_else(|| panic!("build.sh lacks `{needle}`"))
    };
    let steps = [
        at("cargo-bundle-licenses --format yaml --output THIRDPARTY.yml"),
        at("cargo install --no-track --locked --root \"$PREFIX\" --path crates/gloci"),
        at("cargo run --locked -p xtask -- man --out \"$PREFIX/share/man/man1\""),
        at("site-functions"),
        at("vendor_completions.d"),
        at("bash-completion/completions"),
    ];
    assert!(
        steps[..4].windows(2).all(|w| w[0] < w[1]),
        "steps out of order"
    );
}

#[test]
fn py_recipe_installs_no_cli() {
    let t = py();
    assert!(t.contains("test ! -e $PREFIX/bin/gloci"), "{t}");
    assert!(t.contains("import gloci"));
    assert!(t.contains(".count()"));
    let host = block(&t, "requirements").join("\n");
    for need in ["python >=3.10", "pip", "maturin >=1.9.4,<2"] {
        assert!(host.contains(need), "host requirement {need}\n{host}");
    }
    // The wheel bundles every dependency's code, so it ships their licences like the CLI does,
    // and the bundle tool is run (a build requirement nothing runs is dead weight).
    assert!(
        fs::read_to_string(recipes().join("py-gloci/build.sh"))
            .unwrap()
            .contains("cargo-bundle-licenses --format yaml --output THIRDPARTY.yml"),
        "the py-gloci build script must run cargo-bundle-licenses"
    );
    assert!(
        t.contains("- THIRDPARTY.yml") && t.contains("license_file:"),
        "{t}"
    );
    let pyproject =
        fs::read_to_string(xtask::diff::workspace_root().join("crates/gloci-py/pyproject.toml"))
            .unwrap();
    assert!(
        pyproject.contains("maturin>=1.9.4,<2"),
        "the recipe's maturin floor must match pyproject.toml"
    );
    assert!(t.contains("recipe-maintainers:") && t.contains("- duartemolha"));
    assert!(
        !t.contains("cargo install") && !t.contains("entry_points"),
        "{t}"
    );
    assert!(!t.contains("run_exports") || t.contains("python"), "{t}");
}

#[test]
fn recipe_test_files_exist() {
    let d = recipes().join("gloci/test");
    for f in ["a.bed", "b.bed", "expected.bed"] {
        let p = d.join(f);
        assert!(p.is_file(), "{}", p.display());
        assert!(!fs::read_to_string(&p).unwrap().is_empty(), "{f} is empty");
    }
    // The recipe lists them as test files, so conda copies them next to the test commands.
    let t = cli();
    assert!(t.contains("files:"), "{t}");
    for f in ["test/a.bed", "test/b.bed", "test/expected.bed"] {
        assert!(t.contains(f), "{f} not listed\n{t}");
    }
    for r in ["gloci/meta.yaml", "gloci/build.sh", "py-gloci/meta.yaml"] {
        assert!(recipes().join(r).is_file(), "{r}");
    }
    assert!(xtask::diff::workspace_root().join("LICENSE").is_file());
}

#[test]
fn committed_recipes_name_the_workspace_version() {
    let ws = fs::read_to_string(xtask::diff::workspace_root().join("Cargo.toml")).unwrap();
    let v = ws
        .lines()
        .find_map(|l| l.strip_prefix("version = \"")?.strip_suffix('"'))
        .unwrap();
    for r in ["gloci", "py-gloci"] {
        let m = fs::read_to_string(recipes().join(r).join("meta.yaml")).unwrap();
        assert!(m.contains(&format!("{{% set version = \"{v}\" %}}")), "{r}");
    }
}
