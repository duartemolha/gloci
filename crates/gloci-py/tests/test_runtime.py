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
        while not stop.is_set():
            counter[0] += 1
            time.sleep(0.0005)

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
    # A FIFO fed forever: the run cannot end before the interrupt, however fast it is.
    fifo = here / "slow.bed"
    os.mkfifo(fifo)
    stop = threading.Event()

    def feed() -> None:
        try:
            with open(fifo, "w") as f:
                # Enough at once for the format to be read and the output to be created ...
                f.writelines(f"chr1\t{i * 10}\t{i * 10 + 5}\n" for i in range(100_000))
                i = 0
                while not stop.is_set():  # ... then a trickle.
                    f.write(f"chr2\t{i * 10}\t{i * 10 + 5}\n")
                    f.flush()
                    i += 1
                    time.sleep(0.002)
        except (BrokenPipeError, OSError):
            pass

    feeder = threading.Thread(target=feed, daemon=True)
    feeder.start()
    partial: list[Path] = []

    def interrupt() -> None:
        partial.extend(here.glob("*.tmp-*"))
        _thread.interrupt_main()

    timer = threading.Timer(0.3, interrupt)
    start = time.monotonic()
    timer.start()
    try:
        with pytest.raises(KeyboardInterrupt):
            gloci.open("slow.bed").merge(max_gap=1).view().save("out.bed")
        elapsed = time.monotonic() - start
    finally:
        timer.cancel()
        stop.set()
        feeder.join(5)
    assert elapsed < 0.3 + 2, elapsed
    assert partial, "the output's temporary file did not exist when the interrupt came"
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
    for _ in range(2):
        with pytest.raises(gloci.GlociError) as e:
            _native._test_panic()
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
