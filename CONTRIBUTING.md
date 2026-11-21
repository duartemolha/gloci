# Contributing to gloci

Thank you for helping. This page covers reporting problems, building gloci, the checks a change
must pass, and the conventions the project keeps.

## Reporting a problem

Use the [issue forms](https://github.com/duartemolha/gloci/issues/new/choose). If gloci's output
differs from the tool it replaces (bedtools, bedops, bcftools, vcftools, samtools, CrossMap), look
the case up in [docs/differences.md](docs/differences.md) first: the differences gloci makes on
purpose are listed there, with the option that gives the old behaviour where there is one.

Security problems go through a private report instead; see [SECURITY.md](SECURITY.md).

## Setting up

- **Rust:** install [rustup](https://rustup.rs). The toolchain CI uses is pinned in
  `rust-toolchain.toml`, and the first `cargo` command installs it. The oldest Rust gloci builds
  with is `rust-version` in the workspace `Cargo.toml` (1.91).
- **Docker:** only for the differential harness (`cargo xtask diff`) and the benchmarks.
- **Python 3.10 or later** and [maturin](https://www.maturin.rs), for the Python package.

```sh
cargo build -p gloci                       # target/debug/gloci
cargo run -p gloci -- intersect -a a.bed -b b.bed
```

## Where things are

| Path | What it holds |
|---|---|
| `crates/gloci` | The command line: argument parsing, help, the `gloci` binary |
| `crates/gloci-core` | Core types: intervals, genomes, errors, sort orders |
| `crates/gloci-io` | Readers and writers: BED, GFF, VCF, BCF, BAM, CRAM, BGZF, indexes, built-in genomes |
| `crates/gloci-engine` | The interval index, sweep and external-sort engines |
| `crates/gloci-ops` | The interval and VCF operations behind the subcommands |
| `crates/gloci-vcf` | VCF and BCF record operations |
| `crates/gloci-chains`, `gloci-chains-t2t` | The built-in liftover chains |
| `crates/gloci-py` | The Python package (PyO3, built with maturin) |
| `xtask` | Maintenance tasks: `cargo xtask --help` |
| `docs/` | Design, decisions, differences, migration tables and the user guide (`docs/book`) |

The crates other than `gloci` and the Python package are internal: they have no stable API.

## The checks a change must pass

CI runs these on Linux (x86_64 and aarch64) and macOS. Run them before opening a pull request:

```sh
cargo fmt --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
```

While working, run the tests of what you touched (`cargo test -p gloci --test intersect`) and leave
the full run for the end.

Some changes need one more step:

- **A changed dependency:** run `cargo xtask licenses` and commit `THIRD-PARTY-LICENSES` (CI checks it
  is current).
- **A changed command, option or help text:** run `cargo xtask pystub` (the Python type stub) and
  `cargo xtask docs --skill SKILL.md`; tests fail while either is out of date. The book's command
  reference is generated when the book is built and is not committed.
- **A changed workflow** under `.github/workflows`: run actionlint as CI does,
  `docker run --rm -v "$PWD:/repo" -w /repo rhysd/actionlint:1.7.12 -ignore SC2086`.
- **A change to the Python package:** run its tests against a fresh build:

  ```sh
  python -m venv .venv && . .venv/bin/activate
  pip install -r crates/gloci-py/requirements-ci.txt
  maturin develop -m crates/gloci-py/Cargo.toml --features test-hooks
  cargo build -p gloci
  GLOCI_BIN="$PWD/target/debug/gloci" GLOCI_EXPECT_VERSION=0.1.0 \
      pytest crates/gloci-py/tests -m "not large"
  ```

  `GLOCI_BIN` makes the tests compare the library with the command line; `GLOCI_EXPECT_VERSION` is
  the workspace version in `Cargo.toml`. With a conda environment active, `unset CONDA_PREFIX`
  before `maturin develop`.

The large-input tests (multi-GB files, memory ceilings) are slow and run in the weekly nightly
workflow: `cargo test --release -p gloci --test large -- --ignored`, and `pytest -m large` for
Python.

## Matching the tools gloci replaces

gloci aims to give the same output as the tool each command replaces, and checks it:

- Each shared behaviour is a case in `xtask/diff/cases.toml`. Its reference output is committed under
  `crates/gloci/tests/golden/<case>/`, and every `cargo test` compares gloci with it, without Docker.
- To add or update a case, build the image of pinned reference tools and record the output:

  ```sh
  docker build -t gloci-diff:p0 xtask/docker
  cargo xtask diff --only <case> --update-golden
  ```

  Commit the golden files with the change.
- A difference made on purpose gets a section in [docs/differences.md](docs/differences.md) with an
  anchor, and the case names it (`expected_difference = "#anchor"`). gloci's own output for that case
  is committed too, so the difference cannot drift; if it stops differing, the case fails as stale.

## Conventions

- **Decisions:** design decisions are numbered in [docs/decisions.md](docs/decisions.md). A change of
  behaviour that is a choice adds a decision or amends one; [docs/design.md](docs/design.md)
  describes the whole.
- **Changelog:** a change users will notice goes under `[Unreleased]` in
  [CHANGELOG.md](CHANGELOG.md), under Added, Changed, Fixed or Security.
- **Stability:** changing the options or output of a stable-core command is a breaking change; see
  [docs/book/src/stability.md](docs/book/src/stability.md). Experimental commands may change, but the
  changelog says so.
- **Commit messages:** `area: what changed`, for example `fix(py): …`, `ci: …` or `docs: …`, with a
  body that says why. `Fixes #N` closes the issue it fixes.
- **Dependencies:** a new dependency must have a permissive licence (MIT, Apache-2.0, BSD, ISC,
  zlib). Anything else needs the maintainer's agreement first. GPL tools may serve as test-only
  references in the Docker images, never in what ships. No code, binaries or licence from the UCSC
  Genome Browser tools may enter the repository: their licence does not allow free commercial use,
  and gloci must stay free for commercial use.

## Pull requests

Keep a pull request to one change, say how you tested it, and link the issue it addresses. For a new
difference from a reference tool, link its section in `docs/differences.md`. CI must pass before a
pull request is merged.
