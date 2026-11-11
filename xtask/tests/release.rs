use std::fs;
use std::path::Path;
use xtask::release::{changelog_sections, release_check, release_notes};

const CL: &str = "# Changelog\n\n## [Unreleased]\n\n- soon\n\n## [0.2.0] - 2026-12-01\n\n- two\n- more\n\n## [0.1.0] - 2026-11-08\n\n- one\n";

#[test]
fn changelog_parses_keep_a_changelog() {
    let s = changelog_sections(CL).unwrap();
    assert_eq!(s.len(), 3);
    assert_eq!(s[0].version, None);
    assert_eq!(s[0].body, "- soon");
    assert_eq!(s[1].version.as_deref(), Some("0.2.0"));
    assert_eq!(s[1].date.as_deref(), Some("2026-12-01"));
    assert_eq!(s[2].body, "- one");
    // Unreleased must come first, versions must decrease, dates must be ISO.
    assert!(changelog_sections("## [0.1.0] - 2026-11-08\n- a\n## [Unreleased]\n").is_err());
    assert!(changelog_sections("## [0.1.0] - 2026-11-08\n## [0.2.0] - 2026-12-01\n").is_err());
    assert!(changelog_sections("## [0.1.0] - 8 Nov 2026\n").is_err());
}

#[test]
fn footer_links_and_calendar_dates_and_cff_comments() {
    let cl = format!(
        "{CL}\n[Unreleased]: https://example.org/compare/v0.2.0...HEAD\n[0.1.0]: https://example.org/releases/v0.1.0\n"
    );
    let s = changelog_sections(&cl).unwrap();
    assert_eq!(s[2].body, "- one");
    assert_eq!(release_notes(&cl, "0.1.0").unwrap(), "- one");
    assert!(changelog_sections("## [0.1.0] - 2026-02-31\n").is_err());
    assert!(changelog_sections("## [0.1.0] - 2025-02-29\n").is_err());
    assert!(changelog_sections("## [0.1.0] - 2024-02-29\n").is_ok());

    let d = tempfile::tempdir().unwrap();
    fixture(d.path());
    fs::write(
        d.path().join("CITATION.cff"),
        "cff-version: 1.2.0\nversion: 0.2.0 # the release\ndate-released: 2026-12-01 # ISO\n",
    )
    .unwrap();
    assert_eq!(release_check(d.path(), None).unwrap(), Vec::<String>::new());
}

#[test]
fn release_notes_extracts_section() {
    assert_eq!(release_notes(CL, "0.2.0").unwrap(), "- two\n- more");
    let e = release_notes(CL, "0.3.0").unwrap_err().to_string();
    assert!(e.contains("[0.3.0]"), "{e}");
}

fn fixture(root: &Path) {
    let w = |p: &str, c: &str| {
        let p = root.join(p);
        fs::create_dir_all(p.parent().unwrap()).unwrap();
        fs::write(p, c).unwrap();
    };
    w(
        "Cargo.toml",
        "[workspace]\nmembers = [\"crates/*\"]\n[workspace.package]\nversion = \"0.2.0\"\n[workspace.dependencies]\nfoo-core = { path = \"crates/foo-core\", version = \"0.2.0\" }\n",
    );
    w(
        "crates/foo/Cargo.toml",
        "[package]\nname = \"foo\"\nversion.workspace = true\n[dependencies]\nfoo-core = { path = \"../foo-core\", version = \"0.2.0\" }\n",
    );
    w("CHANGELOG.md", CL);
    w(
        "CITATION.cff",
        "cff-version: 1.2.0\nversion: 0.2.0\ndate-released: 2026-12-01\n",
    );
    w("LICENSE", "MIT\n");
}

#[test]
fn release_check_reports_each_mismatch() {
    let d = tempfile::tempdir().unwrap();
    let r = d.path();
    fixture(r);
    assert_eq!(
        release_check(r, Some("v0.2.0")).unwrap(),
        Vec::<String>::new()
    );

    let one = |tag: Option<&str>, needle: &str| {
        let p = release_check(r, tag).unwrap();
        assert_eq!(p.len(), 1, "{p:?}");
        assert!(p[0].contains(needle), "{p:?}");
    };
    one(Some("v0.3.0"), "does not match");

    let cl = fs::read_to_string(r.join("CHANGELOG.md")).unwrap();
    fs::write(
        r.join("CHANGELOG.md"),
        cl.replace("## [0.2.0] - 2026-12-01", "## [0.1.5] - 2026-12-01"),
    )
    .unwrap();
    one(None, "no `## [0.2.0]");
    fs::write(
        r.join("CHANGELOG.md"),
        cl.replace("- two\n- more\n", "prose only\n"),
    )
    .unwrap();
    one(None, "no bullet");
    fs::write(r.join("CHANGELOG.md"), &cl).unwrap();

    let cff = fs::read_to_string(r.join("CITATION.cff")).unwrap();
    fs::write(
        r.join("CITATION.cff"),
        cff.replace("version: 0.2.0", "version: 0.1.0"),
    )
    .unwrap();
    one(None, "CITATION.cff version");
    fs::write(
        r.join("CITATION.cff"),
        cff.replace("2026-12-01", "2026-11-08"),
    )
    .unwrap();
    one(None, "date-released");
    fs::write(r.join("CITATION.cff"), &cff).unwrap();

    fs::remove_file(r.join("LICENSE")).unwrap();
    one(None, "LICENSE");
    fs::write(r.join("LICENSE"), "MIT\n").unwrap();

    let m = fs::read_to_string(r.join("crates/foo/Cargo.toml")).unwrap();
    fs::write(
        r.join("crates/foo/Cargo.toml"),
        m.replace("\"0.2.0\"", "\"0.1.0\""),
    )
    .unwrap();
    one(None, "crates/foo/Cargo.toml");
}

#[test]
fn repo_metadata_is_consistent() {
    let root = xtask::diff::workspace_root();
    let p = release_check(&root, None).unwrap();
    assert!(p.is_empty(), "{p:#?}");
}
