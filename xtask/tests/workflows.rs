use std::fs;

#[test]
fn release_workflow_and_dockerfile_are_valid() {
    let root = xtask::diff::workspace_root();

    // 1. Dockerfile.release checks (D59 / Task 30)
    let dockerfile_path = root.join("xtask/docker/Dockerfile.release");
    assert!(
        dockerfile_path.exists(),
        "xtask/docker/Dockerfile.release must exist"
    );
    let dockerfile = fs::read_to_string(&dockerfile_path).unwrap();
    assert!(
        dockerfile.contains("FROM scratch"),
        "Dockerfile.release must use FROM scratch"
    );
    assert!(
        dockerfile.contains("ARG BIN="),
        "Dockerfile.release must declare ARG BIN"
    );
    assert!(
        dockerfile.contains("COPY ${BIN} /gloci"),
        "Dockerfile.release must COPY ${{BIN}} /gloci"
    );
    assert!(
        dockerfile.contains("ENTRYPOINT [\"/gloci\"]"),
        "Dockerfile.release must have ENTRYPOINT [\"/gloci\"]"
    );

    // 2. release.yml checks (D59 / Task 30)
    let release_path = root.join(".github/workflows/release.yml");
    assert!(
        release_path.exists(),
        ".github/workflows/release.yml must exist"
    );
    let release_yml = fs::read_to_string(&release_path).unwrap();

    // Triggered on v* tag push only
    assert!(
        release_yml.contains("tags:"),
        "release.yml must trigger on tags"
    );
    assert!(
        release_yml.contains("v*"),
        "release.yml must trigger on v* tags"
    );

    // Required targets
    for target in [
        "x86_64-unknown-linux-musl",
        "aarch64-unknown-linux-musl",
        "aarch64-apple-darwin",
    ] {
        assert!(
            release_yml.contains(target),
            "release.yml missing target {target}"
        );
    }

    // Static ELF interpreter check, --version smoke test, SHA256SUMS, GHCR push
    assert!(
        release_yml.contains("readelf -l") && release_yml.contains("INTERP"),
        "release.yml must assert static musl binary via readelf -l / INTERP"
    );
    assert!(
        release_yml.contains("--version"),
        "release.yml must run gloci --version on native builds"
    );
    assert!(
        release_yml.contains("SHA256SUMS"),
        "release.yml must generate SHA256SUMS"
    );
    assert!(
        release_yml.contains("Dockerfile.release"),
        "release.yml must build xtask/docker/Dockerfile.release"
    );
    assert!(
        release_yml.contains("ghcr.io/"),
        "release.yml must push to ghcr.io"
    );
    assert!(
        release_yml.contains("packages: write"),
        "release.yml must declare permissions: packages: write"
    );
    assert!(
        release_yml.contains("contents: write"),
        "release.yml must declare permissions: contents: write"
    );

    // Nothing published to crates.io or bioconda
    assert!(
        !release_yml.contains("cargo publish") && !release_yml.contains("bioconda"),
        "release.yml must not publish to crates.io or bioconda"
    );

    // The tag must equal the workspace version, and the release is created
    // only after the image push (the last step of the publish job).
    assert!(
        release_yml.contains("cargo xtask release-check --tag"),
        "release.yml must check the tag and release metadata with xtask release-check"
    );
    assert!(
        release_yml.contains("cargo xtask release-notes")
            && release_yml.contains("body_path: release-notes.md")
            && !release_yml.contains("generate_release_notes"),
        "release.yml must take the release body from CHANGELOG.md via release-notes"
    );
    let push = release_yml.find("name: Build and push GHCR image").unwrap();
    let release = release_yml.find("name: Create GitHub Release").unwrap();
    assert!(
        push < release,
        "release.yml must create the GitHub release after the image push"
    );

    // 3. ci.yml must include actionlint
    let ci_path = root.join(".github/workflows/ci.yml");
    let ci_yml = fs::read_to_string(&ci_path).unwrap();
    assert!(
        ci_yml.contains("actionlint"),
        "ci.yml must include an actionlint step"
    );

    // 4. Every action in release.yml, ci.yml and nightly.yml must be pinned
    // to a 40-char commit SHA. NOTE: this checks the pins' format only, not
    // that the commits exist or match the `# vX.Y.Z` labels (that needs the
    // network; look a pin up with `gh api repos/OWNER/REPO/commits/TAG`).
    let nightly_yml = fs::read_to_string(root.join(".github/workflows/nightly.yml")).unwrap();
    for (name, content) in [
        ("release.yml", &release_yml),
        ("ci.yml", &ci_yml),
        ("nightly.yml", &nightly_yml),
    ] {
        for (lineno, line) in content.lines().enumerate() {
            let trimmed = line.trim_start();
            let rest = trimmed
                .strip_prefix("- uses:")
                .or_else(|| trimmed.strip_prefix("uses:"));
            if let Some(rest) = rest {
                let value = rest.split('#').next().unwrap().trim();
                let (_, rev) = value
                    .split_once('@')
                    .unwrap_or_else(|| panic!("{name}:{} uses: {value} missing @sha", lineno + 1));
                assert!(
                    rev.len() == 40 && rev.bytes().all(|b| b.is_ascii_hexdigit()),
                    "{name}:{} action {value:?} is not pinned to a 40-char hex commit SHA",
                    lineno + 1
                );
            }
        }
    }
}
