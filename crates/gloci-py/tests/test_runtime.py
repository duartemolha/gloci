"""The runtime (D91): the GIL is released, Ctrl-C cancels, panics and logs reach Python."""

from __future__ import annotations

import _thread
import logging
import os
import threading
import time
from pathlib import Path

import pytest

import gloci
from gloci import _native


@pytest.fixture
def here(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    monkeypatch.chdir(tmp_path)
    return tmp_path


def test_gil_released(pair, here: Path) -> None:
    a, b = pair(500_000)
    counter = [0]
    stop = threading.Event()

    def count() -> None:
        # sleep(0) gives the GIL up at every step, so the main thread is never kept from it,
        # and the counter runs whenever the main thread does not hold it.
        while not stop.is_set():
            counter[0] += 1
            time.sleep(0)

    t = threading.Thread(target=count)
    t.start()
    try:
        time.sleep(0.05)
        before = counter[0]
        gloci.open(a).intersect(b=b).save("o.bed")
        advanced = counter[0] - before
    finally:
        stop.set()
        t.join()
    assert (here / "o.bed").stat().st_size > 0
    # A run that held the GIL throughout would let the counter take at most a step or two.
    assert advanced >= 100, advanced


def test_keyboard_interrupt_cancels(here: Path) -> None:
    # A FIFO fed for 10 s: the run cannot end before the interrupt, however fast it is.
    fifo = here / "slow.bed"
    os.mkfifo(fifo)
    stop = threading.Event()

    def feed() -> None:
        try:
            with open(fifo, "w") as f:
                # Enough at once (about 6 MB) to pass the first step's output buffer and the
                # pipe's chunk, so the last step creates its output file ...
                f.writelines(f"chr1\t{i * 10}\t{i * 10 + 5}\n" for i in range(300_000))
                i = 0
                # ... then a trickle, for 10 s at most: a cancel that does not work fails the
                # test rather than hanging it.
                end = time.monotonic() + 10
                while not stop.is_set() and time.monotonic() < end:
                    f.write(f"chr2\t{i * 10}\t{i * 10 + 5}\n")
                    f.flush()
                    i += 1
                    time.sleep(0.002)
        except (BrokenPipeError, OSError):
            pass

    feeder = threading.Thread(target=feed, daemon=True)
    feeder.start()
    partial: list[Path] = []
    interrupted: list[float] = []

    def interrupt() -> None:
        # Interrupt as soon as the output's temporary file exists (within 5 s), however loaded
        # the machine is; interrupt anyway after that, so a missing file fails the test below.
        deadline = time.monotonic() + 5
        while not partial and time.monotonic() < deadline:
            partial.extend(here.glob("*.tmp-*"))
            time.sleep(0.005)
        interrupted.append(time.monotonic())
        _thread.interrupt_main()

    watcher = threading.Thread(target=interrupt, daemon=True)
    watcher.start()
    try:
        with pytest.raises(KeyboardInterrupt):
            gloci.open("slow.bed").merge(max_gap=1).view().save("out.bed")
        raised = time.monotonic()
    finally:
        stop.set()
        feeder.join(5)
        watcher.join(10)
    assert partial, "the output's temporary file did not appear within 5 s"
    # The cancel reaches every step and the run ends well before the feeder's 10 s.
    assert raised - interrupted[0] < 2, raised - interrupted[0]
    assert not list(here.glob("*.tmp-*"))
    assert not (here / "out.bed").exists()
    # The next call works: the cancel was this run's alone.
    (here / "x.bed").write_text("chr1\t1\t2\n")
    assert gloci.open("x.bed").merge().to_bytes() == b"chr1\t1\t2\n"


def test_panic_is_glocierror(here: Path) -> None:
    if not hasattr(_native, "_test_panic"):
        if os.environ.get("CI") == "true":
            pytest.fail("CI builds gloci-py with --features test-hooks")
        pytest.skip("built without the test-hooks feature")
    # In a stage (twice: the first leaves nothing behind), in the driver's own code, and in
    # the driver while a step that needs the GIL to end is running (no deadlock).
    for place in ("stage", "stage", "driver", "unwind"):
        with pytest.raises(gloci.GlociError) as e:
            _native._test_panic(place)
        assert type(e.value) is gloci.GlociError
        assert str(e.value).startswith("internal error: ")
        assert "test panic" in str(e.value) and "this is a bug" in str(e.value)
    # The interpreter survives and the next call works.
    (here / "x.bed").write_text("chr1\t1\t2\n")
    assert gloci.open("x.bed").merge().to_bytes() == b"chr1\t1\t2\n"


def test_warnings_go_to_logging(here: Path, caplog: pytest.LogCaptureFixture) -> None:
    (here / "a.bed").write_text("chr1\t10\t20\nchr1\t50\t60\n")
    (here / "nochr.bed").write_text("1\t10\t20\n1\t50\t60\n")
    p = gloci.open("a.bed").intersect(b="nochr.bed")
    caplog.set_level(logging.WARNING, logger="gloci")
    assert p.to_bytes() == b""
    got = [r for r in caplog.records if r.name == "gloci"]
    assert [r.levelno for r in got] == [logging.WARNING], got
    msg = got[0].getMessage()
    assert msg.startswith("a.bed and nochr.bed share no chromosome names")
    assert 'differ by a "chr" prefix' in msg
    # The level is the logger's, read at each run: at ERROR the warning is not even made.
    caplog.clear()
    caplog.set_level(logging.ERROR, logger="gloci")
    assert p.to_bytes() == b""
    assert not [r for r in caplog.records if r.name == "gloci"]


def test_concurrent_runs(pair, here: Path, caplog: pytest.LogCaptureFixture) -> None:
    a, b = pair(20_000, seed=3)
    (here / "nochr.bed").write_text("1\t10\t20\n")
    p = gloci.open(a).intersect(b=b).merge(max_gap=10)
    want = p.to_bytes()
    warn = gloci.open(a).intersect(b="nochr.bed")
    results: dict[int, bytes] = {}

    def work(i: int) -> None:
        results[i] = (warn if i % 2 else p).to_bytes(threads=2)

    caplog.set_level(logging.WARNING, logger="gloci")
    threads = [threading.Thread(target=work, args=(i,)) for i in range(6)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert [results[i] for i in range(0, 6, 2)] == [want] * 3
    assert [results[i] for i in range(1, 6, 2)] == [b""] * 3
    # Each run's warning is forwarded once, whichever thread's driver forwards it.
    got = [r for r in caplog.records if r.name == "gloci"]
    assert len(got) == 3 and all("share no chromosome names" in r.getMessage() for r in got)
