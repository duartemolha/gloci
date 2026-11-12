"""Wheel smoke test: needs only an installed gloci wheel (no pytest, no repository files).

    GLOCI_EXPECT_VERSION=0.1.0 python test_smoke.py

pytest also collects it. The release workflow runs it against each built wheel.
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

import gloci


def test_import_and_version() -> None:
    expect = os.environ.get("GLOCI_EXPECT_VERSION")
    assert expect, "set GLOCI_EXPECT_VERSION to the version the wheel must report"
    assert gloci.__version__ == expect, (gloci.__version__, expect)


def test_genome_data_is_in_the_wheel() -> None:
    # slop with -g hg38 needs the built-in hg38 chromosome sizes.
    out = gloci.open([("chr1", 1, 5)]).slop(b=1, g="hg38").to_bytes()
    assert out == b"chr1\t0\t6\n", out


def test_bad_input_raises_data_error() -> None:
    with tempfile.TemporaryDirectory() as d:
        bad = Path(d) / "bad.bed"
        bad.write_text("chr1\t10\t20\nchr1\tx\t50\n")
        try:
            gloci.open(str(bad)).merge().to_bytes()
        except gloci.DataError as e:
            assert e.line == 2, e.line
        else:
            raise AssertionError("expected gloci.DataError")


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok {name}")
    sys.exit(0)
