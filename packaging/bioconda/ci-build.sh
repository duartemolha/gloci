#!/bin/bash
# Build both recipes and install them together, inside a condaforge/miniforge3 container.
#   ci-build.sh RECIPES_DIR OUT_DIR
# RECIPES_DIR holds the rendered recipes (`cargo xtask conda render --local ...`) and the
# source trees they point at. OUT_DIR is writable; the packages land in OUT_DIR/bld.
# Nothing is uploaded anywhere.
set -euxo pipefail
recipes="$1"
out="$2"

# The container may run as an unprivileged user: keep conda's state under OUT_DIR.
export HOME="$out/home" CONDA_PKGS_DIRS="$out/pkgs"
mkdir -p "$HOME" "$out/bld"
# At most 2 jobs: the runner (or a shared dev box) has few cores.
export CARGO_BUILD_JOBS=2 CPU_COUNT=2

conda config --set channel_priority strict
# conda-forge-pinning brings the variant config bioconda builds with: it defines c_stdlib and
# c_stdlib_version, which the recipes' `{{ stdlib('c') }}` needs.
conda create -y -p "$out/build-env" -c conda-forge -c bioconda --override-channels --strict-channel-priority conda-build conda-forge-pinning
# conda-build lives in its own environment; use that environment's conda from here on.
export PATH="$out/build-env/bin:$PATH"

for r in gloci py-gloci; do
  conda build -m "$out/build-env/conda_build_config.yaml" -c conda-forge -c bioconda --override-channels --output-folder "$out/bld" "$recipes/$r"
done

# Both packages in one fresh environment; the local packages win over the channels.
conda create -y -p "$out/test-env" --strict-channel-priority -c "$out/bld" -c conda-forge -c bioconda --override-channels gloci py-gloci

env="$out/test-env"
cli_version="$("$env/bin/gloci" --version | awk '{print $2}')"
py_version="$("$env/bin/python" -c 'import gloci; print(gloci.__version__)')"
echo "gloci --version: $cli_version; gloci.__version__: $py_version"
test "$cli_version" = "$py_version"

# No file is installed by both packages.
"$env/bin/python" - "$env" <<'PY'
import glob, json, sys
env = sys.argv[1]
owners = {}
for meta in glob.glob(f"{env}/conda-meta/*.json"):
    m = json.load(open(meta))
    if m["name"] in ("gloci", "py-gloci"):
        for f in m["files"]:
            owners.setdefault(f, []).append(m["name"])
twice = {f: o for f, o in owners.items() if len(o) > 1}
assert owners and not twice, twice
print(f"{len(owners)} files, none installed twice")
PY
