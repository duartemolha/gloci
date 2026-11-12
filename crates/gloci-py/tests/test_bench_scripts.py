"""The benchmark scripts of xtask/bench/py (D96): the gloci ones write the CLI's bytes."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from conftest import REPO

SCRIPTS = REPO / "xtask" / "bench" / "py"
TOOLS = ["gloci", "pybedtools", "polars_bio", "pyranges1", "bioframe"]
OPS = ["intersect", "merge"]


def run_script(name: str, *args: Path) -> None:
    done = subprocess.run(
        [sys.executable, "-I", str(SCRIPTS / name), *map(str, args)],
        capture_output=True,
    )
    assert done.returncode == 0, done.stderr.decode()


def test_every_tool_has_both_scripts() -> None:
    for tool in TOOLS:
        for op in OPS:
            script = SCRIPTS / f"{tool}_{op}.py"
            assert script.is_file(), f"missing {script}"
            text = script.read_text()
            assert "sys.argv" in text, f"{script.name} must take its paths from argv"


def test_gloci_bench_scripts_match_cli(ab: Path, cli) -> None:
    done = cli(
        ["intersect", "--chrom-order", "lex", "--mode", "any", "-a", "a.bed", "-b", "b.bed", "-o", "i.cli.bed"],
        ab,
    )
    assert done.returncode == 0, done.stderr
    run_script("gloci_intersect.py", ab / "a.bed", ab / "b.bed", ab / "i.py.bed")
    want = (ab / "i.cli.bed").read_bytes()
    assert len(want.splitlines()) > 100
    assert (ab / "i.py.bed").read_bytes() == want

    done = cli(["merge", "--chrom-order", "lex", "-o", "m.cli.bed", "a.bed"], ab)
    assert done.returncode == 0, done.stderr
    run_script("gloci_merge.py", ab / "a.bed", ab / "m.py.bed")
    want = (ab / "m.cli.bed").read_bytes()
    assert len(want.splitlines()) > 100
    assert (ab / "m.py.bed").read_bytes() == want


@pytest.mark.parametrize("tool", ["pybedtools", "polars_bio", "pyranges1", "bioframe"])
@pytest.mark.parametrize("op", OPS)
def test_other_tools_agree_with_gloci_on_records(ab: Path, tool: str, op: str) -> None:
    """Where the tool is installed, it finds the same intervals (order and format may differ)."""
    module = {"pybedtools": "pybedtools", "polars_bio": "polars_bio",
              "pyranges1": "pyranges1", "bioframe": "bioframe"}[tool]
    pytest.importorskip(module)
    if tool == "pybedtools" and shutil.which("bedtools") is None:
        pytest.skip("pybedtools needs the bedtools binary")
    paths = [ab / "a.bed", ab / "b.bed"][: 2 if op == "intersect" else 1]
    run_script(f"gloci_{op}.py", *paths, ab / "want.bed")
    run_script(f"{tool}_{op}.py", *paths, ab / "got.bed")

    def intervals(p: Path) -> list[tuple[str, int, int]]:
        rows = [tuple(line.split("\t")[:3]) for line in p.read_text().splitlines() if line]
        return sorted((c, int(s), int(e)) for c, s, e in rows)

    assert intervals(ab / "got.bed") == intervals(ab / "want.bed")
