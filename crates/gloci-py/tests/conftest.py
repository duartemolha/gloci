"""Shared fixtures: the gloci CLI to compare against, and generated inputs."""

from __future__ import annotations

import os
import random
import shutil
import subprocess
from pathlib import Path

import pytest

TESTS = Path(__file__).resolve().parent
REPO = TESTS.parents[2]
# Committed, sorted BED fixtures (1,200 and 450 records on chr1-chr3).
FIXTURES = REPO / "crates" / "gloci" / "tests" / "embed_stages"


@pytest.fixture(scope="session")
def gloci_bin() -> str:
    """The gloci binary named by $GLOCI_BIN; required when CI=true, else its tests skip."""
    path = os.environ.get("GLOCI_BIN")
    if path:
        if not os.access(path, os.X_OK):
            pytest.fail(f"GLOCI_BIN={path} is not an executable file")
        return path
    if os.environ.get("CI") == "true":
        pytest.fail(
            "GLOCI_BIN must name the gloci binary when CI=true: "
            "cargo build -p gloci, then export GLOCI_BIN=$PWD/target/debug/gloci"
        )
    pytest.skip("set GLOCI_BIN to the gloci binary to compare with the CLI")


@pytest.fixture(scope="session")
def gloci_path_dir(gloci_bin: str, tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A directory whose `gloci` is $GLOCI_BIN, to put on PATH for `sh -c`."""
    d = tmp_path_factory.mktemp("bin")
    (d / "gloci").symlink_to(os.path.abspath(gloci_bin))
    return d


@pytest.fixture
def cli(gloci_bin: str):
    """`cli(args, cwd)` runs the CLI with argv0 `gloci`, as a user types it."""

    def run(args: list[str], cwd: Path | str) -> subprocess.CompletedProcess[bytes]:
        return subprocess.run(
            ["gloci", *args], executable=gloci_bin, cwd=cwd, capture_output=True
        )

    return run


@pytest.fixture
def sh(gloci_path_dir: Path):
    """`sh(cmd, cwd)` runs a shell command line with `gloci` on PATH."""

    def run(cmd: str, cwd: Path | str) -> subprocess.CompletedProcess[bytes]:
        env = dict(os.environ, PATH=f"{gloci_path_dir}{os.pathsep}{os.environ['PATH']}")
        return subprocess.run(["sh", "-c", cmd], cwd=cwd, env=env, capture_output=True)

    return run


@pytest.fixture
def ab(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """`a.bed` and `b.bed` (the committed fixtures) in a fresh directory, made the cwd."""
    for name in ("a.bed", "b.bed"):
        shutil.copy(FIXTURES / name, tmp_path / name)
    monkeypatch.chdir(tmp_path)
    return tmp_path


def write_sorted_bed(path: Path, n: int, rng: random.Random) -> None:
    """`n` sorted BED3 records over chr1-chr3, of 1-500 bases each."""
    per = [n // 3 + (1 if i < n % 3 else 0) for i in range(3)]
    with open(path, "w") as f:
        for c, k in enumerate(per, start=1):
            starts = sorted(rng.randrange(0, 50_000_000) for _ in range(k))
            f.writelines(f"chr{c}\t{s}\t{s + rng.randint(1, 500)}\n" for s in starts)


@pytest.fixture
def pair(tmp_path: Path):
    """`pair(n, seed)` writes two sorted BED files of `n` records each and returns their paths."""

    def make(n: int, seed: int = 0) -> tuple[Path, Path]:
        rng = random.Random(seed)
        a, b = tmp_path / f"pair-a-{n}-{seed}.bed", tmp_path / f"pair-b-{n}-{seed}.bed"
        write_sorted_bed(a, n, rng)
        write_sorted_bed(b, n, rng)
        return a, b

    return make
