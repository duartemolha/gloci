# Releasing

This page is for the maintainers. It records how a release is cut and which first-release steps are
left to the repository owner. Nothing here has been done yet: no package has been published, no Pages site
deployed, and no release workflow run on GitHub.

## What a release does

A push of a `v*` tag runs `.github/workflows/release.yml`:

1. `cargo xtask release-check --tag <tag>` requires that the tag is `v<workspace version>`, that
   `CHANGELOG.md` has a section for the version, that `CITATION.cff` carries it, and that `LICENSE` exists.
2. Static musl archives for Linux x86_64 and aarch64 and a macOS arm64 archive, each with the man pages,
   shell completions, `SKILL.md`, `docs/migration.md`, `docs/differences.md`, `CHANGELOG.md` and the
   licence.
3. abi3 wheels (manylinux2014 x86_64 and aarch64, macOS arm64) and an sdist, each smoke-tested on its
   platform (the wheels on the oldest and newest supported CPython).
4. The GitHub Release: the archives, wheels, `SHA256SUMS`, the notes from the `CHANGELOG.md` section, and
   `conda-recipes.tar.gz` (the bioconda recipes rendered with the real checksums).
5. Gated publishing, each step off until the owner turns it on (see below): PyPI, crates.io and the
   `ghcr.io` image.

`.github/workflows/docs.yml` builds this book on every pull request and deploys it to GitHub Pages only
on a `v*` tag push, and only when the repository variable `GLOCI_PUBLISH_PAGES` is `true`.

On a private repository the release workflow skips the conda-recipe rendering (it downloads the public tag
archive), so the Release carries no `conda-recipes.tar.gz`. The bioconda recipes need a public repository.

A manual run of the release workflow (`workflow_dispatch`) builds and tests everything and publishes
nothing, whatever `dry_run` says.

## Release checklist

1. Bump the version (see [First-release steps](#first-release-steps) for the first one) in the workspace
   `Cargo.toml`, the crate dependencies' `version` fields, `CITATION.cff` and the two bioconda recipes'
   `{% set version = ... %}` lines (`packaging/bioconda/gloci/meta.yaml` and `py-gloci/meta.yaml`). The
   Python package takes its version from Cargo, so it needs no edit of its own. Move the
   `[Unreleased]` entries of `CHANGELOG.md` under the new version and date.
2. `cargo xtask release-check` passes.
3. CI is green on `main`, including the nightly differential, benchmark and `conda-recipe` jobs.
4. Re-measure performance on a quiet machine and compare with the previous release's results (results are
   kept per release).
5. Tag `v<version>` and push the tag. The release workflow and the docs workflow run.
6. Check the Release assets against `SHA256SUMS`, then follow the steps for the channels that are enabled.

## First-release steps

These steps are the owner's. The project has deferred them on purpose; each is a checklist item to do once.

- [ ] **PyPI trusted publisher.** On PyPI, add a trusted publisher for the project `gloci` (owner
  `duartemolha`, repository `gloci`, workflow `release.yml`, environment `pypi`). In the repository
  settings, create the `pypi` environment, which the `publish-pypi` job uses.
- [ ] **crates.io token.** Add the repository secret `CARGO_REGISTRY_TOKEN`, a crates.io API token that may
  publish the workspace crates. crates.io rate-limits new crates (HTTP 429): if the publish job stops on
  one, re-run it. Each crate already on crates.io is skipped, so the re-run completes the publish.
- [ ] **Repository variables.** Set `GLOCI_PUBLISH_PYPI=true`, `GLOCI_PUBLISH_CRATES=true`,
  `GLOCI_PUBLISH_GHCR=true` and `GLOCI_PUBLISH_PAGES=true` when the matching channel should go live. Each job also needs a real tag push
  and a `LICENSE` file.
- [ ] **Dry run.** Dispatch the release workflow once with `dry_run` set (it is the default) and check that
  every build and smoke test passes before the first tag.
- [ ] **Bioconda recipe check.** Dispatch the nightly workflow once and check its `conda-recipe` job: it builds
  both recipes (`gloci` and `py-gloci`) and installs them into a fresh environment.
- [ ] **Bioconda submission.** After a tagged release, take `conda-recipes.tar.gz` from the Release and open a
  pull request against `bioconda-recipes` with the `gloci` and `py-gloci` recipes.
- [ ] **Zenodo.** Create the Zenodo record for the release (with the paper's evidence bundle, from
  `cargo xtask paper bundle`) and put its DOI in `CITATION.cff`.
- [ ] **Pages.** In the repository settings, enable GitHub Pages with the source "GitHub Actions", then set
  `GLOCI_PUBLISH_PAGES=true`. The docs workflow deploys this site on the next `v*` tag; until the variable is
  set the deploy job is skipped, so the first tag does not fail on a disabled Pages site.
- [ ] **Paper benchmarks.** Run the paper benchmarks (`cargo xtask paper fetch`, then `cargo xtask paper run`)
  on the 8-core x86_64 reference machine. CI never downloads the paper datasets.
- [ ] **Quiet-box performance.** Re-measure the performance numbers (`cargo xtask bench`) on a quiet box
  before the release, not on a shared machine.
- [ ] **Version 0.2.0.** Bump the version to 0.2.0 at the release (D81): workspace `Cargo.toml`, crate
  dependency versions, `CITATION.cff`, the changelog section and the bioconda recipes' version.

Until the owner has done the relevant steps, the [install page](../install.md) lists crates.io, PyPI,
bioconda and `ghcr.io` only under "Once published".
