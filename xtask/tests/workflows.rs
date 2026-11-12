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

    // Nothing goes to bioconda; PyPI and crates.io publishing exist but are gated (below).
    assert!(
        !release_yml.contains("bioconda"),
        "release.yml must not publish to bioconda"
    );
    for needle in [
        "workflow_dispatch",
        "dry_run",
        "man --out",
        "completions --out",
        "release-check",
        "maturin",
        "manylinux: 2014",
        "quay.io/pypa/manylinux2014_x86_64",
        "quay.io/pypa/manylinux2014_aarch64",
        "auditwheel show",
        "manylinux_2_17",
        "gh-action-pypi-publish",
        "skip-existing: true",
        "id-token: write",
        "cargo publish",
        "cargo info",
        "--profile release-py",
        "MACOSX_DEPLOYMENT_TARGET",
        "test_smoke.py",
    ] {
        assert!(
            release_yml.contains(needle),
            "release.yml must contain `{needle}`"
        );
    }

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
    for job in ["publish-pypi:", "publish-crates:"] {
        let at = release_yml
            .find(&format!("\n  {job}"))
            .unwrap_or_else(|| panic!("release.yml must have job {job}"));
        assert!(at < release, "the Release must come after {job}");
    }

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

/// The text of one top-level job of a workflow.
fn job_body<'a>(yml: &'a str, job: &str) -> &'a str {
    let at = yml
        .find(&format!("\n  {job}:\n"))
        .unwrap_or_else(|| panic!("no job {job}"))
        + 1;
    let rest = &yml[at + 1..];
    let end = rest
        .lines()
        .scan(0, |off, l| {
            let o = *off;
            *off += l.len() + 1;
            Some((o, l))
        })
        .skip(1)
        .find(|(_, l)| l.starts_with("  ") && !l.starts_with("   ") && l.trim_end().ends_with(':'))
        .map_or(rest.len(), |(o, _)| o);
    &yml[at..at + 1 + end]
}

/// Every action of `yml` is pinned to a 40-hex commit, the same one ci.yml and release.yml pin.
fn assert_pinned_like_the_others(name: &str, yml: &str) {
    let root = xtask::diff::workspace_root();
    let read = |n: &str| fs::read_to_string(root.join(".github/workflows").join(n)).unwrap();
    let pinned: Vec<(String, String)> = [read("ci.yml"), read("release.yml"), read("python.yml")]
        .iter()
        .flat_map(|y| actions(y))
        .collect();
    for (action, rev) in actions(yml) {
        assert!(
            rev.len() == 40 && rev.bytes().all(|b| b.is_ascii_hexdigit()),
            "{name}: {action}@{rev} is not pinned to a 40-char hex commit SHA"
        );
        for (a, r) in &pinned {
            assert!(
                *a != action || *r == rev,
                "{name} pins {action}@{rev}, but another workflow pins {a}@{r}"
            );
        }
    }
}

/// `steps` appear in `body` in this order.
fn assert_in_order(what: &str, body: &str, steps: &[&str]) {
    let at = |s: &str| {
        body.find(s)
            .unwrap_or_else(|| panic!("{what} must contain `{s}`"))
    };
    for pair in steps.windows(2) {
        assert!(
            at(pair[0]) < at(pair[1]),
            "{what}: `{}` must come before `{}`",
            pair[0],
            pair[1]
        );
    }
}

#[test]
fn python_jobs_export_the_version_the_smoke_test_expects() {
    let root = xtask::diff::workspace_root();
    let py = fs::read_to_string(root.join(".github/workflows/python.yml")).unwrap();
    // test_smoke.py checks gloci.__version__ against GLOCI_EXPECT_VERSION: both jobs that run
    // pytest export the workspace version, as release.yml's version job reads it.
    for job in ["pytest", "polars-only"] {
        assert_in_order(
            &format!("python.yml {job}"),
            job_body(&py, job),
            &[
                "v=\"$(sed -n 's/^version = \"\\(.*\\)\"$/\\1/p' Cargo.toml | head -1)\"",
                "test -n \"$v\"",
                "echo \"GLOCI_EXPECT_VERSION=$v\" >> \"$GITHUB_ENV\"",
                "pytest crates/gloci-py/tests",
            ],
        );
    }
}

#[test]
fn python_workflow_tests_polars_without_pyarrow() {
    let root = xtask::diff::workspace_root();
    let py = fs::read_to_string(root.join(".github/workflows/python.yml")).unwrap();
    // The Arrow output must work with polars alone (D89, P2b Task 8).
    let job = job_body(&py, "polars-only");
    assert_pinned_like_the_others("python.yml", job);
    assert_in_order(
        "python.yml polars-only",
        job,
        &[
            "cargo build -p gloci\n",
            "echo \"GLOCI_BIN=$PWD/target/debug/gloci\" >> \"$GITHUB_ENV\"",
            "pip install -r crates/gloci-py/requirements-polars-only.txt",
            "-c \"import pyarrow\"",
            "exit 1",
            "maturin develop -m crates/gloci-py/Cargo.toml --features test-hooks",
            "pytest crates/gloci-py/tests -m \"not large\"",
        ],
    );
    // Its requirements have polars and no pyarrow (nor pandas, which needs it).
    let reqs =
        fs::read_to_string(root.join("crates/gloci-py/requirements-polars-only.txt")).unwrap();
    let pins: Vec<&str> = reqs
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .collect();
    assert!(pins.iter().any(|l| l.starts_with("polars==")), "{pins:?}");
    assert!(
        !pins
            .iter()
            .any(|l| l.starts_with("pyarrow") || l.starts_with("pandas")),
        "{pins:?}"
    );
}

#[test]
fn nightly_runs_the_large_python_tests() {
    let root = xtask::diff::workspace_root();
    let y = fs::read_to_string(root.join(".github/workflows/nightly.yml")).unwrap();
    assert!(
        y.contains("  schedule:\n"),
        "nightly.yml runs on a schedule"
    );
    let job = job_body(&y, "python-large");
    assert_pinned_like_the_others("nightly.yml", job);
    assert!(
        actions(job)
            .iter()
            .any(|(a, _)| a == "actions/setup-python")
    );
    // Built as the wheels are (release-py), then only the large tests.
    assert_in_order(
        "nightly.yml python-large",
        job,
        &[
            "pip install -r crates/gloci-py/requirements-ci.txt",
            "maturin develop -m crates/gloci-py/Cargo.toml --profile release-py",
            "pytest crates/gloci-py/tests -m large\n",
        ],
    );
    assert!(
        !job.contains("test-hooks"),
        "the large tests run on a wheel-like build"
    );
}

#[test]
fn publish_steps_are_gated() {
    let root = xtask::diff::workspace_root();
    let y = fs::read_to_string(root.join(".github/workflows/release.yml")).unwrap();

    // Each publish job is off unless its repository variable is true (D80), and runs
    // only for a real tag push with a LICENSE in the tree.
    for (job, var) in [
        ("publish-pypi", "GLOCI_PUBLISH_PYPI"),
        ("publish-crates", "GLOCI_PUBLISH_CRATES"),
    ] {
        let body = job_body(&y, job);
        let cond = body
            .lines()
            .find(|l| l.trim_start().starts_with("if:"))
            .unwrap_or_else(|| panic!("{job} has no job-level `if:`"));
        assert!(
            cond.contains(&format!("vars.{var} == 'true'")),
            "{job}: `{cond}` must test vars.{var}"
        );
        assert!(
            cond.contains("github.event_name == 'push'"),
            "{job}: `{cond}` must require a tag push (not workflow_dispatch)"
        );
        assert!(
            cond.contains("needs.version.outputs.license == 'true'"),
            "{job}: `{cond}` must require the LICENSE file (the version job's `license` output)"
        );
    }
    let version = job_body(&y, "version");
    assert!(version.contains("[ -f LICENSE ]") && version.contains("license="));
    assert!(job_body(&y, "publish-pypi").contains("environment: pypi"));
    assert!(job_body(&y, "publish-pypi").contains("gh-action-pypi-publish"));
    assert!(job_body(&y, "publish-crates").contains("CARGO_REGISTRY_TOKEN"));

    // Nothing publishes outside those jobs.
    let within = |job: &str, needle: &str| {
        let body = job_body(&y, job);
        let start = body.as_ptr() as usize - y.as_ptr() as usize;
        for (i, _) in y.match_indices(needle) {
            assert!(
                i >= start && i < start + body.len(),
                "`{needle}` appears outside {job} (at byte {i})"
            );
        }
    };
    within("publish-pypi", "gh-action-pypi-publish");
    within("publish-crates", "cargo publish");
    within("ghcr", "build-push-action");
    let ghcr = job_body(&y, "ghcr");
    assert!(ghcr.contains("if: github.event_name == 'push'"));
    for step in [
        "Log in to GitHub Container Registry",
        "Build and push GHCR image",
    ] {
        let at = ghcr.find(step).unwrap();
        assert!(
            ghcr[at..]
                .lines()
                .nth(1)
                .unwrap()
                .contains("if: vars.GLOCI_PUBLISH_GHCR == 'true'"),
            "{step} must be gated"
        );
    }

    // The Release job runs when a publish job is skipped, but not after a failure,
    // and only on a tag push.
    let rel = job_body(&y, "release");
    assert!(
        rel.contains("if: ${{ !failure() && !cancelled() && github.event_name == 'push' }}"),
        "the release job must run when publish jobs are skipped"
    );
    for need in ["publish-pypi", "publish-crates", "wheels", "sdist", "build"] {
        assert!(rel.contains(need), "release job must need {need}");
    }

    // The smoke test runs against a wheel on the Linux images and on macOS.
    let smoke = job_body(&y, "wheel-smoke");
    assert!(smoke.contains("cp310-cp310") && smoke.contains("--no-index"));
    assert!(job_body(&y, "sdist-smoke").contains("python:3.12"));

    // The new actions are pinned like the others.
    let used = actions(&y);
    assert!(used.iter().any(|(a, _)| a == "PyO3/maturin-action"));
    assert!(used.iter().any(|(a, _)| a == "pypa/gh-action-pypi-publish"));
}

#[test]
fn nightly_builds_and_installs_the_bioconda_recipes() {
    let root = xtask::diff::workspace_root();
    let y = fs::read_to_string(root.join(".github/workflows/nightly.yml")).unwrap();
    let job = job_body(&y, "conda-recipe");
    // The pinned miniforge image: a digest, never a moving tag.
    let image = job
        .split_whitespace()
        .find(|w| w.starts_with("condaforge/miniforge3"))
        .expect("conda-recipe runs in condaforge/miniforge3");
    let (_, digest) = image
        .split_once("@sha256:")
        .expect("image pinned by digest");
    assert!(
        digest.len() == 64 && digest.bytes().all(|b| b.is_ascii_hexdigit()),
        "{image}"
    );
    assert!(
        job.contains("conda render") && job.contains("--local"),
        "render with --local"
    );
    // The in-container steps live in one script that the job and a local run share.
    let script = fs::read_to_string(root.join("packaging/bioconda/ci-build.sh")).unwrap();
    assert!(job.contains("packaging/bioconda/ci-build.sh"));
    for needle in [
        "channel_priority strict",
        "--strict-channel-priority",
        "-c conda-forge -c bioconda",
        "--override-channels",
        "conda build",
        "gloci --version",
        "gloci.__version__",
    ] {
        assert!(script.contains(needle), "ci-build.sh lacks `{needle}`");
    }
    assert!(
        script.contains("gloci py-gloci"),
        "one env with both packages must install"
    );
    assert!(
        !job.contains("secrets."),
        "the recipe check publishes nothing"
    );
}

#[test]
fn release_renders_and_attaches_the_recipes() {
    let root = xtask::diff::workspace_root();
    let y = fs::read_to_string(root.join(".github/workflows/release.yml")).unwrap();
    let body = job_body(&y, "release");
    assert!(
        body.contains("conda render"),
        "release must render the recipes"
    );
    assert!(body.contains("--sha256-cli") && body.contains("--sha256-py"));
    assert!(
        body.contains("conda-recipes"),
        "recipes are attached to the Release"
    );
}

#[test]
fn docs_workflow_is_valid() {
    let root = xtask::diff::workspace_root();
    let read = |name: &str| fs::read_to_string(root.join(".github/workflows").join(name)).unwrap();
    let y = read("docs.yml");

    // Builds on pull requests; deploys only on a v* tag push.
    assert!(
        y.contains("pull_request:"),
        "docs.yml must build on pull requests"
    );
    assert!(
        y.contains("tags:") && y.contains("\"v*\""),
        "docs.yml must run on v* tags"
    );
    assert!(
        !y.contains("branches:"),
        "docs.yml must not run on branch pushes"
    );
    assert!(
        !y.contains("workflow_dispatch"),
        "docs.yml must not be dispatchable"
    );

    let build = job_body(&y, "build");
    assert!(
        build.contains("tool: mdbook@"),
        "mdbook must be pinned with tool: mdbook@<version>"
    );
    assert!(build.contains("cargo xtask book"));
    assert!(build.contains("mdbook build docs/book"));
    assert!(build.contains("upload-pages-artifact"));
    assert!(
        !build.contains("deploy-pages"),
        "the build job must not deploy"
    );

    let deploy = job_body(&y, "deploy");
    assert!(deploy.contains("actions/deploy-pages@"));
    assert!(deploy.contains("pages: write"), "deploy needs pages: write");
    assert!(
        deploy.contains("id-token: write"),
        "deploy needs id-token: write"
    );
    assert!(
        deploy.contains("startsWith(github.ref, 'refs/tags/v')"),
        "deploy must require a v* tag push"
    );
    assert!(
        deploy.contains("github-pages"),
        "deploy uses the github-pages environment"
    );
    assert!(!build.contains("pages: write") && !build.contains("id-token: write"));

    // Every action is pinned to a 40-hex SHA, the same one the other workflows pin.
    let pinned: Vec<(String, String)> = [read("ci.yml"), read("release.yml")]
        .iter()
        .flat_map(|y| actions(y))
        .collect();
    let used = actions(&y);
    for needed in [
        "actions/checkout",
        "taiki-e/install-action",
        "actions/deploy-pages",
    ] {
        assert!(
            used.iter().any(|(a, _)| a == needed),
            "docs.yml must use {needed}"
        );
    }
    for (action, rev) in &used {
        assert!(
            rev.len() == 40 && rev.bytes().all(|b| b.is_ascii_hexdigit()),
            "docs.yml: {action}@{rev} is not pinned to a 40-char hex commit SHA"
        );
        for (a, r) in &pinned {
            assert!(
                a != action || r == rev,
                "docs.yml pins {action}@{rev}, but the other workflows pin {a}@{r}"
            );
        }
    }
}
