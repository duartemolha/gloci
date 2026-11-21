# Releasing

This page is for the maintainers. It records how a release is cut and which first-release steps are
left to the repository owner. The first release is 0.1.0 (D81). It goes to the GitHub Release, PyPI,
crates.io, `ghcr.io` and this site on GitHub Pages; bioconda and Zenodo come later.

## What a release does

A push of a `v*` tag runs `.github/workflows/release.yml`:

1. `cargo xtask release-check --tag <tag>` requires that the tag is `v<workspace version>`, that
   `CHANGELOG.md` has a section for the version, that `CITATION.cff` carries it, and that `LICENSE` exists.
2. Static musl archives for Linux x86_64 and aarch64 and a macOS arm64 archive, each with the man pages,
   shell completions, `SKILL.md`, `docs/migration.md`, `docs/differences.md`, `CHANGELOG.md`, the
   licence and `THIRD-PARTY-LICENSES`.
3. abi3 wheels (manylinux2014 x86_64 and aarch64, macOS arm64) and an sdist, each smoke-tested on its
   platform (the wheels on the oldest and newest supported CPython).
   Each archive, wheel and the sdist gets a signed build-provenance attestation in the job that builds
   it, before anything is published (D176); a manual run attests nothing.
4. The GitHub Release: the archives, wheels, `SHA256SUMS`, the notes from the `CHANGELOG.md` section, and
   `conda-recipes.tar.gz` (the bioconda recipes rendered with the real checksums).
5. Gated publishing, each step off until the owner turns it on (see below): PyPI, crates.io and the
   `ghcr.io` image.

## Repository rules

Two rulesets protect what a release depends on (D176). They apply to everyone, the owner included:

- **`main`:** no force pushes and no deletion. Direct pushes are allowed.
- **Release tags (`v*`):** no deletion and no moving. A published version always points at the commit
  its files were built from; a wrong release is fixed by a new version, never by moving its tag.

They live under the repository's Settings, Rules, Rulesets, where the owner can change them.

`.github/workflows/docs.yml` builds this book on every pull request and deploys it to GitHub Pages only
on a `v*` tag push, and only when the repository variable `GLOCI_PUBLISH_PAGES` is `true`.

On a private repository the release workflow skips the conda-recipe rendering (it downloads the public tag
archive), so the Release carries no `conda-recipes.tar.gz`. The bioconda recipes need a public repository.

A manual run of the release workflow (`workflow_dispatch`) builds and tests everything and publishes
nothing, whatever `dry_run` says.

## Release checklist

1. Bump the version (the first release, 0.1.0, needed no bump) in the workspace
   `Cargo.toml`, the crate dependencies' `version` fields, `CITATION.cff` and the two bioconda recipes'
   `{% set version = ... %}` lines (`packaging/bioconda/gloci/meta.yaml` and `py-gloci/meta.yaml`). The
   Python package takes its version from Cargo, so it needs no edit of its own. Move the
   `[Unreleased]` entries of `CHANGELOG.md` under the new version and date.
2. `cargo xtask release-check` passes.
3. CI is green on `main`, and so is the latest run of the nightly workflow (large files, differential,
   benchmark, fuzz and `conda-recipe`). It runs weekly (D175); dispatch it for a fresh run when the
   week's commits warrant one.
4. Re-measure performance on a quiet machine and compare with the previous release's results (results are
   kept per release).
5. Tag `v<version>` and push the tag. The release workflow and the docs workflow run.
6. Check the Release assets against `SHA256SUMS`, then follow the steps for the channels that are enabled.

## Dependency updates

Dependabot proposes version updates every Monday (`.github/dependabot.yml`, D177): one pull request
per ecosystem (Cargo, GitHub Actions, the Python test tools) for the minor and patch updates, and one
for each major update. A release waits a week after it is published before Dependabot proposes it.

- **Cargo:** the pull request changes crate versions, so CI's licence check fails until
  `THIRD-PARTY-LICENSES` is regenerated on its branch:

  ```sh
  gh pr checkout <number>
  cargo xtask licenses
  git commit -am "deps(cargo): regenerate THIRD-PARTY-LICENSES" && git push
  ```

  CI is not given a token that could do this by itself: the job would compile the updated crates
  with write access to the repository.
- **Major updates:** read the release notes before merging. Some majors change behaviour (for
  example, `actions/download-artifact` v8 stops unzipping files and fails on a checksum mismatch).
- **Left alone on purpose:** `requirements-floors.txt` (it pins the lowest supported versions) and
  `dtolnay/rust-toolchain` (the Rust version follows the procedure below).

## Bumping the Rust toolchain

CI builds with the exact Rust version named in `rust-toolchain.toml`, so a new Rust release (and its new
clippy lints) reaches CI only when a commit asks for it. To bump it:

1. Edit `channel` in `rust-toolchain.toml` and every `toolchain:` and `rust-toolchain:` value in
   `.github/workflows/*.yml` to the new `X.Y.Z` in one commit. The fuzz job stays on `nightly`.
   `cargo test -p xtask` fails with the file and line of anything left out of step.
2. Run `cargo clippy --workspace --all-targets --keep-going -- -D warnings` with the new version and fix
   what it reports in the same commit.

The oldest supported Rust is `rust-version` in the workspace `Cargo.toml`. The `msrv` job in `ci.yml` runs
`cargo check --workspace --locked` with it, the `sdist-smoke` job in `release.yml` builds the sdist with it,
and the README and the install page name it. Raising it means editing all of these together; the same test
checks that they agree.

## THIRD-PARTY-LICENSES

`THIRD-PARTY-LICENSES` at the repository root holds the notices of the crates compiled into the `gloci`
binary and the Python extension. It is generated, not edited: `cargo xtask licenses` rewrites it from
`Cargo.lock` (D80), and CI runs `cargo xtask licenses --check` (after `cargo fetch --locked`), which fails
when the file is out of date. Regenerate it and commit the result whenever a dependency changes. The
release workflow copies it into every archive next to `LICENSE`, and `crates/gloci-py/pyproject.toml` lists
it in `license-files` so the wheels and the sdist carry it. A new dependency whose package has no licence
file needs the standard text of its licence in `xtask/licenses/`, or the command stops and says which.

## First-release steps

These steps are the owner's. The project has deferred them on purpose; each is a checklist item to do once.

- [ ] **PyPI trusted publisher.** On PyPI, add a trusted publisher for the project `gloci` (owner
  `duartemolha`, repository `gloci`, workflow `release.yml`, environment `pypi`). In the repository
  settings, create the `pypi` environment, which the `publish-pypi` job uses.
- [ ] **crates.io token.** Add the repository secret `CARGO_REGISTRY_TOKEN`, a crates.io API token that may
  publish the workspace crates. crates.io rate-limits new crates (HTTP 429): if the publish job stops on
  one, re-run it. Each crate already on crates.io is skipped, so the re-run completes the publish.
- [x] **Repository variables.** Set `GLOCI_PUBLISH_PYPI=true`, `GLOCI_PUBLISH_CRATES=true`,
  `GLOCI_PUBLISH_GHCR=true` and `GLOCI_PUBLISH_PAGES=true` when the matching channel should go live. Each job also needs a real tag push
  and a `LICENSE` file.
- [x] **Dry run.** Dispatch the release workflow once with `dry_run` set (it is the default) and check that
  every build and smoke test passes before the first tag.
- [x] **Bioconda recipe check.** Dispatch the nightly workflow once and check its `conda-recipe` job: it builds
  both recipes (`gloci` and `py-gloci`) and installs them into a fresh environment.
- [ ] **Bioconda submission.** After a tagged release, take `conda-recipes.tar.gz` from the Release and open a
  pull request against `bioconda-recipes` with the `gloci` and `py-gloci` recipes.
- [ ] **Zenodo.** Create the Zenodo record for the release (with the paper's evidence bundle, from
  `cargo xtask paper bundle`) and put its DOI in `CITATION.cff`.
- [x] **Pages.** In the repository settings, enable GitHub Pages with the source "GitHub Actions", add a
  tag rule `v*` to the `github-pages` environment's deployment branches and tags (by default it accepts only
  the default branch, so a tag's deploy would be refused), then set `GLOCI_PUBLISH_PAGES=true`. The docs workflow deploys this site on the next `v*` tag; until the variable is
  set the deploy job is skipped, so the first tag does not fail on a disabled Pages site.
- [ ] **Paper benchmarks.** Run the paper benchmarks (`cargo xtask paper fetch`, then `cargo xtask paper run`)
  on the 8-core x86_64 reference machine. CI never downloads the paper datasets.
- [ ] **Quiet-box performance.** Re-measure the performance numbers (`cargo xtask bench`) on a quiet box
  before the release, not on a shared machine.
- [x] **Version 0.1.0.** The first release keeps the version the workspace has carried since P0 (D81,
  owner 2026-11-20; it was planned as 0.2.0 until then). Its changelog section and the `CITATION.cff`
  date are set on the day of the tag.

The [install page](../install.md) lists crates.io, PyPI and `ghcr.io` as live from 0.1.0, and bioconda
under "Once published" until its recipes are accepted.
