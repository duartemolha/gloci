"""The README's examples run: every `python` block, in order, in one session in a scratch
directory (a block tagged `no-run` is skipped), and every gloci idiom of the "From pybedtools"
table."""

from __future__ import annotations

import importlib.util
import os
import re
from pathlib import Path
from typing import Any

import pytest

README = Path(__file__).resolve().parents[1] / "README.md"
TEXT = README.read_text()

FENCE = re.compile(r"^```(\S*)[ \t]*([^\n]*)\n(.*?)^```[ \t]*$", re.M | re.S)


def missing(info: str) -> list[str]:
    """The optional packages a block needs (`requires=polars,pyarrow` in its tag) that are not
    installed; its examples are then skipped, as the tests that need them are."""
    mods = [m for w in info.split() if w.startswith("requires=") for m in w[9:].split(",")]
    return [m for m in mods if importlib.util.find_spec(m) is None]


def blocks() -> list[tuple[int, str, str]]:
    """(line, info, code) of each python block."""
    out = []
    for m in FENCE.finditer(TEXT):
        lang, rest, code = m.groups()
        if lang == "python":
            out.append((TEXT.count("\n", 0, m.start()) + 1, rest.strip(), code))
    return out


def pybedtools_idioms() -> list[tuple[str, str]]:
    section = TEXT.split("## From pybedtools", 1)[1]
    section = section.split("\n## ", 1)[0]
    rows = [
        line for line in section.splitlines() if line.startswith("|") and not set(line) <= set("|- ")
    ][1:]
    out = []
    for row in rows:
        cells = [c.strip() for c in row.strip().strip("|").split("|")]
        spans = re.findall(r"`([^`]+)`", cells[1])
        assert len(cells) == 2 and len(spans) == 1, f"one code span in the gloci cell: {row}"
        out.append((cells[0], spans[0]))
    return out


def test_readme_has_the_guide_sections() -> None:
    for heading in ("Sources", "Steps", "Terminals", "Frames", "Errors", "Threads and the GIL"):
        assert re.search(rf"^## .*{heading}", TEXT, re.M), heading
    assert "## From pybedtools" in TEXT


def test_pybedtools_table_has_ten_rows() -> None:
    assert len(pybedtools_idioms()) == 10


def test_readme_examples(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    ns: dict[str, Any] = {"__name__": "readme"}
    ran = 0
    for line, info, code in blocks():
        if "no-run" in info.split() or missing(info):
            continue
        try:
            exec(compile(code, f"README.md:{line}", "exec"), ns)
        except Exception as e:
            raise AssertionError(f"README.md:{line}: {type(e).__name__}: {e}") from e
        ran += 1
    assert ran >= 10, ran
    needs = {"to_polars": ["polars"], "to_arrow(": ["pyarrow"], "to_pandas": ["pandas", "pyarrow"]}
    ns["a"] = ns["gloci"].open("a.bed")  # the table's `a`
    for what, code in pybedtools_idioms():
        if any(m for k, mods in needs.items() if k in code for m in missing("requires=" + ",".join(mods))):
            continue
        try:
            exec(compile(code, f"README.md ({what})", "exec"), ns)
        except Exception as e:
            raise AssertionError(f"README.md, pybedtools idiom {what!r}: {e}") from e
    assert os.getcwd() == str(tmp_path)


def test_no_run_is_rare() -> None:
    tagged = [1 for _, info, _ in blocks() if "no-run" in info.split()]
    assert len(tagged) <= 3, "keep no-run rare: prefer small inline sources"
