# Install

This page lists what works today first, and then, under [Once published](#once-published), the
channels that will exist after the first public release.

## From a GitHub Release archive

Each tagged release attaches one archive per platform to its
[GitHub Release](https://github.com/duartemolha/gloci/releases), named
`gloci-<tag>-<target>.tar.gz`:

| Platform | `<target>` |
|---|---|
| Linux x86_64 (static musl) | `x86_64-unknown-linux-musl` |
| Linux aarch64 (static musl) | `aarch64-unknown-linux-musl` |
| macOS arm64 | `aarch64-apple-darwin` |

```sh
tar -xzf gloci-v0.1.0-x86_64-unknown-linux-musl.tar.gz
install -m 0755 gloci-v0.1.0-x86_64-unknown-linux-musl/gloci ~/.local/bin/gloci
gloci --version
```

The archive also holds `man/` (man pages), `completions/` (bash, zsh and fish), `SKILL.md` (an agent
usage guide), `docs/` (the migration tables and the differences it points to), `CHANGELOG.md` and the
licence. Each release lists a `SHA256SUMS` file; check the archive
with `sha256sum -c --ignore-missing SHA256SUMS`.

## Python wheels from a GitHub Release

The same Release carries the Python wheels (abi3, so one wheel serves every CPython from 3.10: manylinux2014
on x86_64 and aarch64, and macOS arm64). Install one by its URL:

```sh
pip install https://github.com/duartemolha/gloci/releases/download/v<VERSION>/gloci-<VERSION>-cp310-abi3-<platform>.whl
```

Take the exact file name from the Release's asset list. The wheel has no dependencies beyond Python;
`pyarrow` is optional, for the Arrow hand-off. See [gloci for Python](python.md).

## From source

With a Rust toolchain (1.91 or newer):

```sh
cargo install --git https://github.com/duartemolha/gloci gloci
```

This builds the `gloci` binary from the repository's default branch. Add `--tag v<VERSION>` for a release.

## Once published

The channels below are not published yet. The commands are what they will be; until the first public
release they do not work, and nothing here has been uploaded to any of these registries.

- crates.io: `cargo install gloci`
- PyPI: `pip install gloci`
- Bioconda, the command line tool and the Python library:

  ```sh
  conda install -c conda-forge -c bioconda gloci
  conda install -c conda-forge -c bioconda py-gloci
  ```

- A container image on GitHub Container Registry, a `FROM scratch` image holding the static binary:
  `docker run --rm -v "$PWD:/data" -w /data ghcr.io/duartemolha/gloci:<tag> --version`.

The [release checklist](maintainers/releasing.md) lists what has to happen before each of them goes live.
