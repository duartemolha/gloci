"""Large inputs (the `large` marker, run nightly with a `--profile release-py` build:
`pytest -m large`): iteration keeps memory flat, and stopping early stops the steps."""

from __future__ import annotations

import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import pytest

import gloci

pytestmark = pytest.mark.large

MIB = 1 << 20


def max_rss() -> int:
    """This process's own peak resident memory (RUSAGE_SELF), in bytes. Read before and after
    in the same process, the difference is how much the peak grew in between."""
    unit = 1 if sys.platform == "darwin" else 1024  # bytes on macOS, KiB on Linux
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * unit


def threads_now() -> int | None:
    """The process's OS threads: `Threads:` of /proc/self/status (Linux), `ps -M` (macOS)."""
    try:
        with open("/proc/self/status") as f:
            for line in f:
                if line.startswith("Threads:"):
                    return int(line.split()[1])
    except OSError:
        pass
    if sys.platform == "darwin":
        done = subprocess.run(
            ["ps", "-M", "-p", str(os.getpid())], capture_output=True, text=True, check=True
        )
        return len(done.stdout.splitlines()) - 1
    return None


#: Writes `n` BED6 records over chr1-chr3 (inside hg38's lengths): sorted, or in random order.
WRITE_BED = """
import random, sys
path, n, seed, ordered = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4] == "1"
rng = random.Random(seed)
per = [n // 3 + (1 if i < n % 3 else 0) for i in range(3)]
with open(path, "w") as f:
    for c, k in enumerate(per, start=1):
        starts = [rng.randrange(0, 150_000_000) for _ in range(k)]
        if ordered:
            starts.sort()
        f.writelines(
            f"chr{c}\\t{s}\\t{s + 1 + s % 500}\\tr{i}\\t{s % 1000}\\t{'+-'[s & 1]}\\n"
            for i, s in enumerate(starts)
        )
"""


def write_bed(path: Path, n: int, seed: int, *, ordered: bool) -> None:
    """`n` BED6 records, written by another process, so that this one's peak memory is not
    raised by making them."""
    subprocess.run(
        [sys.executable, "-c", WRITE_BED, str(path), str(n), str(seed), "1" if ordered else "0"],
        check=True,
    )


def test_iteration_memory_flat(tmp_path: Path) -> None:
    n = 5_000_000
    path = tmp_path / "big.bed"
    write_bed(path, n, 1, ordered=True)
    assert path.stat().st_size > 150 * MIB

    before = max_rss()
    began = time.monotonic()
    seen = 0
    last = None
    for r in gloci.open(path).slop(b=10, g="hg38"):
        seen += 1
        last = r  # only the latest record is kept
    took = time.monotonic() - began
    grown = max_rss() - before
    assert seen == n
    assert last is not None and last.chrom == "chr3"
    # 5 million records (about 190 MiB of text) went through Python without the peak growing
    # by more than the ceiling.
    assert grown < 200 * MIB, f"the peak resident memory grew by {grown / MIB:.0f} MiB"
    print(
        f"\n{n:,} records in {took:.1f} s; ru_maxrss {before / MIB:.0f} MiB before, "
        f"grew by {grown / MIB:.1f} MiB"
    )


def test_break_stops_workers(tmp_path: Path) -> None:
    n = 2_000_000
    path = tmp_path / "unsorted.bed"
    write_bed(path, n, 2, ordered=False)
    spill = tmp_path / "spill"
    spill.mkdir()
    before = threads_now()

    p = gloci.open(path).sort(tmp_dir=str(spill), max_mem="16M")
    spilled: list[Path] = []
    for r in p:
        # The first record comes from merging the sort's spilled runs, which are still there.
        assert r.chrom == "chr1"
        spilled = list(spill.rglob("*"))
        break
    stopped = time.monotonic()
    assert spilled, "the sort did not spill, so the test shows nothing"

    # Breaking out dropped the iterator, which cancelled and joined the steps at once.
    if before is not None:
        now = threads_now()
        while now is not None and now > before and time.monotonic() - stopped < 1:
            time.sleep(0.01)
            now = threads_now()
        assert now is not None and now <= before, f"{now} threads after 1 s, {before} before"
    assert list(spill.iterdir()) == []
