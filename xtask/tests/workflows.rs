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

    // Every job that runs `cargo xtask` installs the pinned toolchain first.
    let mut jobs: Vec<(String, String)> = Vec::new();
    let mut in_jobs = false;
    for line in release_yml.lines() {
        if line == "jobs:" {
            in_jobs = true;
        } else if in_jobs
            && line.starts_with("  ")
            && !line.starts_with("   ")
            && line.trim_end().ends_with(':')
        {
            jobs.push((line.trim().trim_end_matches(':').to_string(), String::new()));
        } else if let Some(j) = jobs.last_mut() {
            j.1.push_str(line);
            j.1.push('\n');
        }
    }
    assert!(!jobs.is_empty(), "release.yml has no jobs");
    for (name, body) in &jobs {
        if body.contains("cargo xtask") {
            let tc = body.find("dtolnay/rust-toolchain@");
            let x = body.find("cargo xtask").unwrap();
            assert!(
                tc.is_some_and(|t| t < x),
                "release.yml job {name} runs cargo xtask without a pinned dtolnay/rust-toolchain step before it"
            );
            assert!(
                body.contains("Swatinem/rust-cache@"),
                "release.yml job {name} runs cargo xtask without the Swatinem/rust-cache step"
            );
        }
    }

    // The release notes are built before the release is created, and the
    // GHCR login and push stay off unless GLOCI_PUBLISH_GHCR is true (D80).
    let notes = release_yml.find("cargo xtask release-notes").unwrap();
    assert!(
        notes < release,
        "release-notes must run before Create GitHub Release"
    );
    for step in [
        "name: Log in to GitHub Container Registry",
        "name: Build and push GHCR image",
    ] {
        let at = release_yml.find(step).unwrap();
        let next = release_yml[at + 1..]
            .find("      - ")
            .map_or(release_yml.len(), |n| at + 1 + n);
        assert!(
            release_yml[at..next].contains("if: vars.GLOCI_PUBLISH_GHCR == 'true'"),
            "{step} must carry `if: vars.GLOCI_PUBLISH_GHCR == 'true'`"
        );
    }

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

/// The `uses:` actions of a workflow, as (action, rev) in order.
fn actions(yml: &str) -> Vec<(String, String)> {
    yml.lines()
        .filter_map(|line| {
            let t = line.trim_start().trim_start_matches("- ").trim_start();
            let spec = t.strip_prefix("uses:")?.split('#').next()?.trim();
            let (action, rev) = spec.split_once('@')?;
            Some((action.to_string(), rev.to_string()))
        })
        .collect()
}

#[test]
fn python_workflow_runs_the_python_tests() {
    let root = xtask::diff::workspace_root();
    let read = |name: &str| fs::read_to_string(root.join(".github/workflows").join(name)).unwrap();
    let py = read("python.yml");

    // On pull requests and on pushes to main.
    assert!(
        py.contains("\n  pull_request:"),
        "python.yml must run on pull requests"
    );
    assert!(
        py.contains("  push:\n    branches: [main]"),
        "python.yml must run on pushes to main"
    );

    // The matrix: (os, python) pairs.
    let mut pairs = Vec::new();
    let mut os: Option<String> = None;
    for line in py.lines() {
        let t = line.trim();
        if let Some(o) = t.strip_prefix("- os: ") {
            os = Some(o.to_string());
        } else if let Some(v) = t.strip_prefix("python: ")
            && let Some(o) = os.take()
        {
            pairs.push((o, v.trim_matches('"').to_string()));
        }
    }
    let want = [
        ("ubuntu-24.04", "3.10"),
        ("ubuntu-24.04", "3.14"),
        ("ubuntu-24.04-arm", "3.12"),
        ("macos-latest", "3.12"),
    ];
    assert_eq!(
        pairs,
        want.map(|(o, v)| (o.to_string(), v.to_string())),
        "python.yml matrix"
    );
    assert!(py.contains("runs-on: ${{ matrix.os }}"));
    assert!(py.contains("python-version: ${{ matrix.python }}"));

    // Every action is pinned to a commit, the same one ci.yml and release.yml pin.
    let pinned: Vec<(String, String)> = [read("ci.yml"), read("release.yml")]
        .iter()
        .flat_map(|y| actions(y))
        .collect();
    let used = actions(&py);
    assert!(used.iter().any(|(a, _)| a == "actions/setup-python"));
    for (action, rev) in &used {
        assert!(
            rev.len() == 40 && rev.bytes().all(|b| b.is_ascii_hexdigit()),
            "python.yml: {action}@{rev} is not pinned to a 40-char hex commit SHA"
        );
        for (a, r) in &pinned {
            assert!(
                a != action || r == rev,
                "python.yml pins {action}@{rev}, but the other workflows pin {a}@{r}"
            );
        }
    }

    // The steps, in order: the CLI (exported as GLOCI_BIN, which the tests require on CI),
    // the test tools, gloci-py with its test hooks, then pytest without the large tests.
    let at = |s: &str| {
        py.find(s)
            .unwrap_or_else(|| panic!("python.yml must contain `{s}`"))
    };
    let steps = [
        "cargo build -p gloci\n",
        "echo \"GLOCI_BIN=$PWD/target/debug/gloci\" >> \"$GITHUB_ENV\"",
        "pip install -r crates/gloci-py/requirements-ci.txt",
        "maturin develop -m crates/gloci-py/Cargo.toml --features test-hooks",
        "pytest crates/gloci-py/tests -m \"not large\"",
    ];
    for pair in steps.windows(2) {
        assert!(
            at(pair[0]) < at(pair[1]),
            "python.yml: `{}` must come before `{}`",
            pair[0],
            pair[1]
        );
    }

    // It never publishes anything.
    for word in [
        "maturin publish",
        "maturin upload",
        "twine",
        "GLOCI_PUBLISH",
        "secrets.",
    ] {
        assert!(!py.contains(word), "python.yml must not publish ({word})");
    }
}
