"""Python-object sources (D88): iterables, file objects, Arrow frames and nested pipelines, and
the module functions and sources of D87."""

from __future__ import annotations

import gc
import gzip
import io
import math
import os
import random
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

import pytest

import gloci
from gloci import _pipeline

from conftest import FIXTURES


def no_nul(e: BaseException) -> None:
    for value in (str(e), *(repr(v) for v in vars(e).values())):
        assert "\0" not in value and "gloci:" not in value, value


def os_threads() -> int | None:
    """The process's OS threads (Rust's stage threads included), or None where unknown."""
    if os.path.isdir("/proc/self/task"):
        return len(os.listdir("/proc/self/task"))
    if sys.platform == "darwin":
        done = subprocess.run(
            ["ps", "-M", "-p", str(os.getpid())], capture_output=True, text=True, check=True
        )
        return len(done.stdout.splitlines()) - 1
    return None


def assert_threads_back_to(before: int | None) -> None:
    """Within 1 s, the process has no more OS threads than `before`."""
    if before is None:
        return
    deadline = time.monotonic() + 1
    now = os_threads()
    while now is not None and now > before and time.monotonic() < deadline:
        time.sleep(0.01)
        now = os_threads()
    assert now is not None and now <= before, f"{now} threads after 1 s, {before} before"


def bed_rows(path: Path) -> list[tuple[object, ...]]:
    """A BED file's records as tuples, with integer coordinates."""
    rows: list[tuple[object, ...]] = []
    for line in path.read_text().splitlines():
        c, s, e, *rest = line.split("\t")
        rows.append((c, int(s), int(e), *rest))
    return rows


@pytest.fixture
def here(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    monkeypatch.chdir(tmp_path)
    return tmp_path


@pytest.fixture
def bash(gloci_path_dir: Path):
    """`bash(cmd, cwd)`: a bash command line with `gloci` on PATH (process substitution)."""
    exe = shutil.which("bash")
    if exe is None:
        pytest.skip("bash is not installed")

    def run(cmd: str, cwd: Path | str) -> subprocess.CompletedProcess[bytes]:
        env = dict(os.environ, PATH=f"{gloci_path_dir}{os.pathsep}{os.environ['PATH']}")
        return subprocess.run([exe, "-c", cmd], cwd=cwd, env=env, capture_output=True)

    return run


# ---- classify_source -------------------------------------------------------------------------


def test_classify_source(here: Path) -> None:
    kinds = [
        ("a.bed", _pipeline.Path),
        (here / "a.bed", _pipeline.Path),
        (io.StringIO("x"), _pipeline.FileObject),
        (io.BytesIO(b"x"), _pipeline.FileObject),
        ((r for r in [("chr1", 1, 2)]), _pipeline.Lines),
        ([("chr1", 1, 2)], _pipeline.Lines),
        (["chr1\t1\t2"], _pipeline.Lines),
        (gloci.open("a.bed"), _pipeline.Nested),
    ]
    for obj, kind in kinds:
        assert type(_pipeline.classify_source(obj)) is kind, obj
    for bad in (5, 2.5, None, object()):
        with pytest.raises(TypeError) as e:
            _pipeline.classify_source(bad)
        msg = str(e.value)
        for word in ("path", "file object", "iterable", "__arrow_c_stream__", "Pipeline"):
            assert word in msg, msg
        assert type(bad).__name__ in msg
    with pytest.raises(TypeError, match="BytesIO"):
        gloci.open(b"chr1\t1\t2\n")  # bytes are neither a path nor lines: say how to pass them
    assert _pipeline.classify_source([("chr1", 1, 2)]).display == "<list>"
    assert _pipeline.classify_source(r for r in []).display == "<generator>"
    assert _pipeline.classify_source(io.BytesIO(b"")).display == "<BytesIO>"


# ---- iterables -------------------------------------------------------------------------------


def test_iterable_of_tuples(ab: Path) -> None:
    rows = [("chr1", 15, 25, None), ("chr1", 1000, 1200, "x"), ("chr2", 5, 9, None)]
    text = "chr1\t15\t25\t.\nchr1\t1000\t1200\tx\nchr2\t5\t9\t.\n"
    (ab / "rows.bed").write_text(text)
    # Fed as tab-joined str() values, None as `.`: view writes the input as it is.
    assert gloci.open(rows).to_bytes() == text.encode()
    for src in (rows, iter(rows), (r for r in rows), tuple(rows), [list(r) for r in rows]):
        got = gloci.open(src).intersect(b="b.bed").to_bytes()
        assert got == gloci.open("rows.bed").intersect(b="b.bed").to_bytes()
    # str lines, with or without their line breaks, and bytes lines.
    lines = text.splitlines()
    assert gloci.open(lines).to_bytes() == text.encode()
    assert gloci.open([l + "\n" for l in lines]).to_bytes() == text.encode()
    assert gloci.open([l.encode() for l in lines]).to_bytes() == text.encode()
    # The fixture as tuples equals the fixture as a file, through a chain.
    a = bed_rows(ab / "a.bed")
    p = gloci.open(a).merge(max_gap=50).sort()
    assert p.to_bytes() == gloci.open("a.bed").merge(max_gap=50).sort().to_bytes()
    assert p.cli() == "gloci merge --max-gap 50 <python object 1> | gloci sort -"
    assert gloci.open(a).count() == len(a)


def test_tuple_values(here: Path) -> None:
    floats = [0.1, 1 / 3, 1e-05, 0.0001, 1e16, 1e15, 123456789.0, -0.0, 2.5e-300, 5e-324,
              1.7976931348623157e308, float("inf"), float("-inf"), float("nan"), 100.0, 1e22]
    rng = random.Random(5)
    floats += [rng.uniform(-1e6, 1e6) for _ in range(200)]
    floats += [math.ldexp(rng.random(), rng.randint(-1060, 1020)) for _ in range(200)]

    class Shown:
        def __str__(self) -> str:
            return "shown"

    rows = [("chr1", i, i + 1, f) for i, f in enumerate(floats)]
    rows.append(("chr1", 900, 901, 2**70, True, Shown(), -5))
    got = gloci.open(rows).to_bytes().decode().splitlines()
    want = ["\t".join(str(v) for v in r) for r in rows]
    assert got == want
    # A tab (or a line break) in a value would split it: a DataError naming the item.
    for bad in ("b\tc", "b\nc"):
        rows = [("chr1", 1, 2, "a"), ("chr1", 3, 4, bad), ("chr1", 5, 6, "d")]
        with pytest.raises(gloci.DataError) as e:
            gloci.open(rows).to_bytes()
        assert type(e.value) is gloci.DataError
        assert (e.value.path, e.value.line) == ("<list>", 2)
        assert str(e.value).startswith("<list> item 2: "), str(e.value)
        assert ("tab" if "\t" in bad else "line break") in str(e.value)
        no_nul(e.value)
    # An item that is neither a line nor a tuple.
    with pytest.raises(gloci.DataError) as e:
        gloci.open(["chr1\t1\t2", 7]).to_bytes()
    assert (e.value.path, e.value.line) == ("<list>", 2)
    assert "int" in str(e.value)
    # A str line holding a line break inside it is not one line.
    with pytest.raises(gloci.DataError) as e:
        gloci.open(["chr1\t1\t2\nchr1\t3\t4"]).to_bytes()
    assert (e.value.path, e.value.line) == ("<list>", 1)
    # The parser's own errors name the item as its line.
    with pytest.raises(gloci.DataError) as e:
        gloci.open([("chr1", 1, 2), ("chr1", 9, 3)]).merge().to_bytes()
    assert (e.value.path, e.value.line) == ("<list>", 2)
    no_nul(e.value)


def test_generator_exception_propagates(here: Path) -> None:
    boom = RuntimeError("boom")

    def rows():
        for i in range(999):
            yield ("chr1", i * 10, i * 10 + 5)
        raise boom

    gc.collect()
    before = os_threads()
    with pytest.raises(RuntimeError) as e:
        gloci.open(rows()).merge().view().save("out.bed")
    assert e.value is boom
    assert not (here / "out.bed").exists()
    assert not list(here.glob("*.tmp-*"))
    assert_threads_back_to(before)
    # The same with the error far into a large source, and from an extra source.
    big = RuntimeError("late")

    def many():
        for i in range(300_000):
            yield ("chr1", i * 10, i * 10 + 5)
        raise big

    with pytest.raises(RuntimeError) as e:
        gloci.open([("chr1", 1, 2)]).merge(many()).save("m.bed")
    assert e.value is big
    assert not list(here.glob("*.tmp-*")) and not (here / "m.bed").exists()
    assert_threads_back_to(before)
    # The exception is delivered after the items before it, so a data error well before it
    # (beyond what the step reads ahead) is reported instead, as with a file.
    def bad_then_boom():
        yield ("chr1", 1, 2)
        yield ("chr1", 9, 3)
        for i in range(300_000):
            yield ("chr1", 10 + i * 10, 15 + i * 10)
        raise boom

    with pytest.raises(gloci.DataError) as e:
        gloci.open(bad_then_boom()).merge().to_bytes()
    assert e.value.line == 2
    # The next run works.
    assert gloci.open([("chr1", 1, 2)]).to_bytes() == b"chr1\t1\t2\n"


def test_keyboard_interrupt_in_a_source(here: Path) -> None:
    def rows():
        for i in range(100_000):
            yield ("chr1", i * 10, i * 10 + 5)
        raise KeyboardInterrupt

    before = os_threads()
    with pytest.raises(KeyboardInterrupt):
        gloci.open(rows()).merge().save("out.bed")
    assert not (here / "out.bed").exists() and not list(here.glob("*.tmp-*"))
    assert_threads_back_to(before)


def test_a_one_shot_source_is_read_once(here: Path) -> None:
    gen = (r for r in [("chr1", 1, 2)])
    p = gloci.open(gen)
    assert p.to_bytes() == b"chr1\t1\t2\n"
    with pytest.raises(gloci.UsageError) as e:
        p.to_bytes()
    assert str(e.value).startswith("<generator> was read by an earlier run"), str(e.value)
    # cli() reads nothing.
    gen = (r for r in [("chr1", 1, 2)])
    q = gloci.open(gen).merge()
    assert q.cli() == "gloci merge <python object 1>"
    assert q.to_bytes() == b"chr1\t1\t2\n"
    # A list is read again at each run.
    rows = [("chr1", 1, 2)]
    r = gloci.open(rows).merge()
    assert r.to_bytes() == r.to_bytes() == b"chr1\t1\t2\n"
    # One iterator in two places of one run would interleave.
    it = iter([("chr1", 1, 2)])
    with pytest.raises(gloci.UsageError, match="used twice"):
        gloci.open(it).intersect(b=it).to_bytes()


def test_needs_file_steps_refuse_python_sources(here: Path) -> None:
    with pytest.raises(gloci.UsageError) as e:
        gloci.open(["chr1\t1\t2"]).faidx(region="chr1:1-2")
    assert str(e.value) == (
        "faidx needs a file input, not <list>; write it to a file and open that"
    )


# ---- file objects ----------------------------------------------------------------------------


def test_text_and_binary_file_objects(here: Path) -> None:
    text = (FIXTURES / "a.bed").read_text()
    p = here / "a.bed"
    p.write_text(text)
    gz = here / "a.bed.gz"
    gz.write_bytes(gzip.compress(text.encode()))
    want = gloci.open(p).merge(max_gap=20).to_bytes()
    assert want.count(b"\n") > 100
    makers = {
        "StringIO": lambda: io.StringIO(text),
        "BytesIO": lambda: io.BytesIO(text.encode()),
        "open(p)": lambda: open(p),
        "open(p, 'rb')": lambda: open(p, "rb"),
        "open(gz, 'rb')": lambda: open(gz, "rb"),
        "gzip.open": lambda: gzip.open(gz),
    }
    for name, make in makers.items():
        with make() as f:
            assert gloci.open(f).merge(max_gap=20).to_bytes() == want, name
    # Errors name the file object by f.name, else by its type.
    bad = "chr1\t10\t20\nchr1\tx\t50\n"
    (here / "bad.bed").write_text(bad)
    with pytest.raises(gloci.DataError) as e:
        gloci.open(io.BytesIO(bad.encode())).merge().to_bytes()
    assert (e.value.path, e.value.line) == ("<BytesIO>", 2)
    assert str(e.value).startswith("<BytesIO> line 2: ")
    with open(here / "bad.bed", "rb") as f, pytest.raises(gloci.DataError) as e:
        gloci.open(f).merge().to_bytes()
    assert (e.value.path, e.value.line) == (f.name, 2)
    no_nul(e.value)
    # A file object's exception propagates unchanged.
    oops = OSError("disk on fire")

    class Failing(io.RawIOBase):
        def readable(self) -> bool:
            return True

        def readinto(self, b) -> int:  # type: ignore[no-untyped-def]
            raise oops

    with pytest.raises(OSError) as e:
        gloci.open(Failing()).merge().to_bytes()
    assert e.value is oops


@pytest.mark.xfail(strict=True, reason="Task 8: gloci.Record and iteration")
def test_records_and_lines_roundtrip(ab: Path) -> None:
    records = list(gloci.open("a.bed"))
    assert gloci.open(records).to_bytes() == gloci.open("a.bed").to_bytes()
    assert gloci.open(str(r) for r in records).to_bytes() == gloci.open("a.bed").to_bytes()


# ---- Arrow frames ----------------------------------------------------------------------------

COLS = {
    "name": ["a", "b", None, "d"],
    "chrom": ["chr1", "chr1", "chr2", "chr2"],
    "start": [15, 1000, 5, 40],
    "end": [25, 1200, 9, 50],
}
FRAME_TEXT = "chr1\t15\t25\ta\nchr1\t1000\t1200\tb\nchr2\t5\t9\t.\nchr2\t40\t50\td\n"


def make_frame(lib: str, cols: dict, types: dict[str, str] | None = None):  # type: ignore[no-untyped-def]
    """`cols` as a pyarrow Table, polars DataFrame or pandas DataFrame; `types` maps columns
    to int32/int64/uint32/uint64 (nullable)."""
    types = types or {}
    if lib == "pyarrow":
        pa = pytest.importorskip("pyarrow")
        return pa.table({k: pa.array(v, type=getattr(pa, types[k])()) if k in types
                         else pa.array(v) for k, v in cols.items()})
    if lib == "polars":
        pl = pytest.importorskip("polars")
        names = {"int32": pl.Int32, "int64": pl.Int64, "uint32": pl.UInt32, "uint64": pl.UInt64}
        return pl.DataFrame(cols, schema_overrides={k: names[t] for k, t in types.items()})
    pd = pytest.importorskip("pandas")
    pytest.importorskip("pyarrow")  # pandas exports Arrow through pyarrow
    names = {"int32": "Int32", "int64": "Int64", "uint32": "UInt32", "uint64": "UInt64"}
    df = pd.DataFrame({k: (pd.array(v, dtype=names[types[k]]) if k in types else v)
                       for k, v in cols.items()})
    return df


DISPLAY = {"pyarrow": "<pyarrow.Table>", "polars": "<polars.DataFrame>",
           "pandas": "<pandas.DataFrame>"}


@pytest.mark.parametrize("lib", ["pyarrow", "polars", "pandas"])
def test_frame_sources(lib: str, ab: Path) -> None:
    (ab / "f.bed").write_text(FRAME_TEXT)
    want = gloci.open("f.bed").intersect(b="b.bed").to_bytes()
    for types in ({}, {"start": "int32", "end": "uint64"}, {"start": "uint32", "end": "int64"}):
        frame = make_frame(lib, COLS, types)
        # chrom start end first, then the other columns in order; nulls are `.`.
        assert gloci.open(frame).to_bytes() == FRAME_TEXT.encode(), types
        assert gloci.open(frame).intersect(b="b.bed").to_bytes() == want, types
    assert gloci.open(make_frame(lib, COLS)).cli() == "gloci view <python object 1>"
    shown = DISPLAY[lib]
    # A null coordinate: a DataError naming the column and the 1-based row.
    cols = dict(COLS, start=[15, 1000, 5, None])
    with pytest.raises(gloci.DataError) as e:
        gloci.open(make_frame(lib, cols, {"start": "int64"})).to_bytes()
    assert (e.value.path, e.value.line) == (shown, 4)
    assert str(e.value).startswith(f"{shown} row 4: ") and "start" in str(e.value)
    no_nul(e.value)
    # A tab in a value, likewise.
    cols = dict(COLS, name=["a", "b\tc", None, "d"])
    with pytest.raises(gloci.DataError) as e:
        gloci.open(make_frame(lib, cols)).to_bytes()
    assert (e.value.path, e.value.line) == (shown, 2)
    assert "name" in str(e.value) and "tab" in str(e.value)
    # A negative coordinate.
    cols = dict(COLS, start=[15, -1000, 5, 40])
    with pytest.raises(gloci.DataError) as e:
        gloci.open(make_frame(lib, cols)).to_bytes()
    assert (e.value.path, e.value.line) == (shown, 2)
    assert "start" in str(e.value) and "-1000" in str(e.value)
    # A float coordinate.
    cols = dict(COLS, end=[25.0, 1200.0, 9.0, 50.0])
    with pytest.raises(gloci.DataError) as e:
        gloci.open(make_frame(lib, cols)).to_bytes()
    assert e.value.path == shown and "end" in str(e.value) and "integer" in str(e.value)
    # A missing column is a usage error, before anything runs.
    cols = {k: v for k, v in COLS.items() if k != "end"}
    with pytest.raises(gloci.UsageError) as e:
        gloci.open(make_frame(lib, cols)).to_bytes()
    assert "end" in str(e.value) and shown in str(e.value)


def test_frame_types_and_batches(here: Path) -> None:
    pa = pytest.importorskip("pyarrow")
    batches = [
        pa.record_batch({"chrom": pa.array(["chr1"] * 3).dictionary_encode(),
                         "start": pa.array([1, 2, 3], pa.uint64()),
                         "end": pa.array([5, 6, 7], pa.int32()),
                         "f": pa.array([0.1, None, 1e16]),
                         "b": pa.array([True, False, None]),
                         "s": pa.array(["x", "y", "z"], pa.large_string()),
                         "v": pa.array(["p", None, "r"], pa.string_view()),
                         "i": pa.array([-1, 0, 2**40], pa.int64())}),
    ] * 3
    reader = pa.RecordBatchReader.from_batches(batches[0].schema, batches)
    got = gloci.open(reader).to_bytes().decode().splitlines()
    one = ["chr1\t1\t5\t0.1\tTrue\tx\tp\t-1", "chr1\t2\t6\t.\tFalse\ty\t.\t0",
           f"chr1\t3\t7\t1e+16\t.\tz\tr\t{2**40}"]
    assert got == one * 3
    # Rows count across batches.
    bad = [batches[0], batches[0].set_column(1, "start", pa.array([1, None, 3], pa.uint64()))]
    with pytest.raises(gloci.DataError) as e:
        gloci.open(pa.RecordBatchReader.from_batches(bad[0].schema, bad)).to_bytes()
    assert (e.value.path, e.value.line) == ("<pyarrow.RecordBatchReader>", 5)


# ---- nested pipelines ------------------------------------------------------------------------


def test_pipeline_as_b(ab: Path, sh, bash) -> None:
    inner = gloci.open("b.bed").slop(b=10, g="hg38")
    p = gloci.open("a.bed").intersect(b=inner)
    shell = sh("gloci slop -b 10 -g hg38 b.bed > s.bed && gloci intersect -a a.bed -b s.bed", ab)
    assert shell.returncode == 0, shell.stderr
    got = p.to_bytes()
    assert got == shell.stdout
    assert got != gloci.open("a.bed").intersect(b="b.bed").to_bytes()
    # cli() shows the nested pipeline as bash's process substitution, which runs it.
    assert p.cli() == "gloci intersect -a a.bed -b <(gloci slop -b 10 -g hg38 b.bed)"
    run = bash(p.cli(), ab)
    assert run.returncode == 0, run.stderr
    assert run.stdout == got
    # The nested pipeline is left as it was, and runs on its own.
    assert inner.cli() == "gloci slop -b 10 -g hg38 b.bed"
    assert inner.to_bytes() == (ab / "s.bed").read_bytes()
    # As an extra source of an n-ary step, in a list for b=, and as open()'s source.
    m = gloci.open("a.bed").merge(inner.merge())
    assert m.cli() == "gloci merge a.bed <(gloci slop -b 10 -g hg38 b.bed | gloci merge -)"
    assert m.to_bytes() == bash(m.cli(), ab).stdout
    # (With two B inputs, pairs are labelled with each one's name: names= sets them.)
    both = gloci.open("a.bed").intersect(b=[inner, "b.bed"], mode="pairs", names="s,b")
    assert both.to_bytes() == bash(both.cli(), ab).stdout
    q = gloci.open(inner).merge()
    assert q.cli() == "gloci slop -b 10 -g hg38 b.bed | gloci merge -"
    assert q.to_bytes() == inner.merge().to_bytes()
    # An error in the nested pipeline names its step's output as the stream it wrote.
    (ab / "bad.bed").write_text("chr1\t10\t20\nchr1\tx\t5\n")
    with pytest.raises(gloci.DataError) as e:
        gloci.open("a.bed").intersect(b=gloci.open("bad.bed").view()).to_bytes()
    assert (e.value.path, e.value.line) == ("bad.bed", 2)


def test_two_inline_inputs(ab: Path) -> None:
    a, b = bed_rows(ab / "a.bed"), bed_rows(ab / "b.bed")
    want = gloci.open("a.bed").intersect(b="b.bed").to_bytes()
    assert want.count(b"\n") > 100
    p = gloci.open(a).intersect(b=b)
    assert p.cli() == "gloci intersect -a <python object 1> -b <python object 2>"
    assert p.to_bytes() == want
    assert gloci.open(iter(a)).intersect(b=(r for r in b)).to_bytes() == want
    # A source and a path in one list repeat the flag.
    q = gloci.open(a).intersect(b=[iter(b), "b.bed"], mode="pairs", names="x,b")
    assert q.cli() == (
        "gloci intersect -a <python object 1> -b <python object 2> -b b.bed --mode pairs "
        "--names x,b"
    )
    pairs = gloci.open("a.bed").intersect(b=["b.bed", "b.bed"], mode="pairs", names="x,b")
    pairs = pairs.to_bytes()
    assert q.to_bytes() == pairs
    # Extra positional sources of an n-ary step.
    m = gloci.open(a).merge(b, io.StringIO((ab / "b.bed").read_text()))
    assert m.to_bytes() == gloci.open("a.bed").merge("b.bed", "b.bed").to_bytes()


# ---- memory ----------------------------------------------------------------------------------


def rss() -> int | None:
    """The process's resident memory now, in bytes (Linux), or None."""
    try:
        with open("/proc/self/statm") as f:
            return int(f.read().split()[1]) * os.sysconf("SC_PAGE_SIZE")
    except OSError:
        return None


def test_a_large_generator_is_not_materialised(here: Path) -> None:
    import resource
    import tracemalloc

    n = 10_000_000
    gen = (("chr1", i, i + 1) for i in range(n))
    peak = [0]
    stop = threading.Event()
    base = rss()

    def sample() -> None:
        while not stop.is_set():
            now = rss()
            if now is not None:
                peak[0] = max(peak[0], now)
            time.sleep(0.005)

    maxrss_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    sampler = threading.Thread(target=sample)
    sampler.start()
    tracemalloc.start()
    try:
        assert gloci.open(gen).count() == n
        _, python_peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
        stop.set()
        sampler.join()
    mib = 1 << 20
    # Python-side: the generator's items are converted as they come, never kept.
    assert python_peak < 16 * mib, python_peak
    # Resident memory, sampled during the run (Linux): the 10 million lines (about 190 MiB as
    # text) are never held at once.
    if base is not None:
        assert peak[0] - base < 64 * mib, (peak[0] - base) / mib
    # An absolute ceiling: the process's peak grew by less than that, wherever it was before
    # (ru_maxrss is KiB on Linux, bytes on macOS).
    unit = 1 if sys.platform == "darwin" else 1024
    grown = (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss - maxrss_before) * unit
    assert grown < 64 * mib, grown / mib


# ---- module sources and functions (D87) ------------------------------------------------------

GENOME = "chr1\t1000\nchr2\t500\n"
VCF_X = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\n"
    "chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT:DP\t0/1:10\t0/0:12\t1/1:8\n"
    "chr1\t200\t.\tG\tT\t40\tPASS\t.\tGT:DP\t0/0:9\t0/1:11\t0/1:7\n"
    "chr1\t300\t.\tC\tG\t30\tPASS\t.\tGT:DP\t1/1:15\t./.:.\t0/1:9\n"
    "chr1\t400\t.\tT\tA\t20\tPASS\t.\tGT:DP\t0/1:6\t0/1:14\t0/0:10\n"
)
VCF_Y = VCF_X.replace("0/0:12", "0/1:12").replace("1/1:15", "0/1:15")


def test_module_sources(here: Path, cli, sh) -> None:
    (here / "g.txt").write_text(GENOME)
    (here / "t.bed").write_text("chr1\t50\t250\nchr2\t0\t30\n")
    w = gloci.windows(g="g.txt", w=100)
    assert isinstance(w, gloci.Pipeline)
    assert w.cli() == "gloci windows -g g.txt -w 100"
    done = cli(["windows", "-g", "g.txt", "-w", "100"], here)
    assert done.returncode == 0, done.stderr
    assert w.to_bytes() == done.stdout and done.stdout.count(b"\n") == 15
    r = gloci.random(g="g.txt", l=10, n=20, seed=7)
    assert r.cli() == "gloci random -g g.txt -l 10 -n 20 --seed 7"
    done = cli(["random", "-g", "g.txt", "-l", "10", "-n", "20", "--seed", "7"], here)
    assert done.returncode == 0, done.stderr
    assert r.to_bytes() == done.stdout and done.stdout.count(b"\n") == 20
    # Sources start chains, are saved, counted and nested like any pipeline.
    chain = w.intersect(b="t.bed").merge()
    assert chain.cli() == "gloci windows -g g.txt -w 100 | gloci intersect -a - -b t.bed | gloci merge -"
    shell = sh(chain.cli(), here)
    assert shell.returncode == 0, shell.stderr
    assert chain.to_bytes() == shell.stdout
    w.save("w.bed.gz", write_index=True)
    assert (here / "w.bed.gz.tbi").exists()
    assert gloci.open("w.bed.gz").to_bytes() == w.to_bytes()
    assert r.count() == 20
    counted = gloci.open("t.bed").intersect(b=w, mode="count").to_bytes()
    assert counted == b"chr1\t50\t250\t3\nchr2\t0\t30\t1\n"
    # The step windows() reads its input with -i; the source cuts a genome with -g.
    assert gloci.open("t.bed").windows(w=100).cli() == "gloci windows -i t.bed -w 100"
    with pytest.raises(gloci.UsageError, match="output is set by save"):
        gloci.windows(g="g.txt", w=100, output="x.bed")
    with pytest.raises(gloci.UsageError, match="unexpected keyword argument 'fooo'"):
        gloci.random(g="g.txt", fooo=1)
    with pytest.raises(TypeError):
        gloci.windows("g.txt")  # type: ignore[misc]


def cli_files(cli, here: Path, args: list[str], d: str) -> dict[str, bytes]:
    (here / d).mkdir()
    done = cli(args, here / d)
    assert done.returncode == 0, done.stderr
    return {p.name: p.read_bytes() for p in sorted((here / d).iterdir())}


def written(paths: list[Path], d: Path) -> dict[str, bytes]:
    assert paths == sorted(paths)
    assert all(isinstance(p, Path) for p in paths)
    return {p.name: (d / p.name).read_bytes() for p in paths}


def test_module_functions(here: Path, cli) -> None:
    for d in ("in", "py"):
        (here / d).mkdir()
        (here / d / "x.vcf").write_text(VCF_X)
        (here / d / "y.vcf").write_text(VCF_Y)
    (here / "py" / "old.frq").write_text("left alone\n")
    x, y = str(here / "in" / "x.vcf"), str(here / "in" / "y.vcf")

    # vcf_stats: PREFIX.<ext> files, as `gloci vcf stats --out-prefix PREFIX`.
    want = cli_files(cli, here, ["vcf", "stats", "--out-prefix", "s", "--freq", "--het",
                                 "--summary", x], "c1")
    got = gloci.vcf_stats(x, out=here / "py" / "s", freq=True, het=True, summary=True)
    assert got == sorted(here / "py" / n for n in want)
    assert written(got, here / "py") == want
    # Any source: lines from a generator give the same files.
    got = gloci.vcf_stats((l for l in VCF_X.splitlines()), out=here / "py" / "g",
                          freq=True, het=True, summary=True)
    assert [p.name for p in got] == [n.replace("s.", "g.", 1) for n in want]
    assert [p.read_bytes() for p in got] == list(want.values())
    # A relative prefix gives relative paths.
    os.chdir(here / "py")
    assert gloci.vcf_stats("x.vcf", out="r", freq=True) == [Path("r.frq")]
    os.chdir(here)

    # vcf_extract_matrix: the 012 matrix files.
    want = cli_files(cli, here, ["vcf", "extract", "--genotype-matrix", "012",
                                 "--out-prefix", "m", x], "c2")
    got = gloci.vcf_extract_matrix(x, out=here / "py" / "m")
    assert [p.name for p in got] == ["m.012", "m.012.indv", "m.012.pos"]
    assert written(got, here / "py") == want
    got = gloci.vcf_extract_matrix(io.StringIO(VCF_X), out=here / "py" / "n", biallelic=True)
    assert [p.read_bytes() for p in got] == list(want.values())

    # vcf_compare: two sources.
    want = cli_files(cli, here, ["vcf", "compare", "--out-prefix", "c", x, y], "c3")
    got = gloci.vcf_compare(x, y, out=here / "py" / "c")
    assert written(got, here / "py") == want
    got = gloci.vcf_compare(gloci.open(x).view(), io.BytesIO(VCF_Y.encode()),
                            out=here / "py" / "d")
    assert [p.read_bytes() for p in got] == list(want.values())

    # index: the index next to the file, as `gloci index FILE`.
    gloci.open(FIXTURES / "a.bed").save(here / "py" / "a.bed.gz")
    shutil.copy(here / "py" / "a.bed.gz", here / "in" / "a.bed.gz")
    done = cli(["index", "a.bed.gz"], here / "in")
    assert done.returncode == 0, done.stderr
    got = gloci.index(here / "py" / "a.bed.gz")
    assert got == [here / "py" / "a.bed.gz.tbi"]
    assert got[0].read_bytes() == (here / "in" / "a.bed.gz.tbi").read_bytes()
    # Unrelated files with the prefix are not reported, nor touched.
    assert (here / "py" / "old.frq").read_text() == "left alone\n"

    # index needs a file.
    for src in (io.BytesIO(b""), [("chr1", 1, 2)], gloci.open("a.bed.gz")):
        with pytest.raises(gloci.UsageError) as e:
            gloci.index(src)  # type: ignore[arg-type]
        assert str(e.value).startswith("index needs a file"), str(e.value)
    # Keyword arguments are checked as the command's, before anything runs.
    with pytest.raises(gloci.UsageError, match="unexpected keyword argument 'fooo'"):
        gloci.vcf_stats(x, out="z", fooo=1)
    with pytest.raises(gloci.UsageError, match="file is an argument of the function"):
        gloci.vcf_stats(x, out="z", file="y.vcf")
    with pytest.raises(TypeError):
        gloci.vcf_stats(x, freq=True)  # type: ignore[call-arg]
    # Errors are the CLI's, as exceptions.
    with pytest.raises(gloci.UsageError, match="no statistic requested"):
        gloci.vcf_stats(x, out=here / "py" / "e")
    with pytest.raises(gloci.DataError) as e:
        gloci.vcf_stats(["#CHROM\tPOS", "chr1\tx"], out=here / "py" / "e", freq=True)
    assert e.value.path == "<list>"
    assert not list((here / "py").glob("e.*"))


def test_a_stalled_step_stops_the_pump(here: Path) -> None:
    # intersect reads the generator only as far as its sweep can go, and a FIFO trickles its
    # b= records: the pump must stop when the source's channel is full, not read on.
    fifo = here / "slow.bed"
    os.mkfifo(fifo)
    n = 3_000_000
    read = [0]
    stalled: list[int] = []

    def rows():
        for i in range(n):
            read[0] += 1
            yield ("chr1", i * 10, i * 10 + 5)

    def feed() -> None:
        with open(fifo, "w") as f:
            for i in range(100):
                f.write(f"chr1\t{i * 10}\t{i * 10 + 1}\n")
                f.flush()
                time.sleep(0.01)
            # The sweep still waits for the rest of b= here.
            stalled.append(read[0])

    # A daemon: if the run fails before it opens the FIFO, the feeder must not hang the exit.
    feeder = threading.Thread(target=feed, daemon=True)
    feeder.start()
    try:
        got = gloci.open(rows()).intersect(b="slow.bed").count()
    finally:
        feeder.join(10)
    assert got == 100
    assert read[0] == n
    assert stalled[0] < 1_000_000, f"{stalled[0]} of {n} items read while the step waited"


# ---- fix round 1 -----------------------------------------------------------------------------


def test_module_functions_report_what_the_run_wrote(here: Path) -> None:
    # I1: the files the run itself committed, not what changed in the directory meanwhile.
    (here / "x.vcf").write_text(VCF_X)
    (here / "s.het").write_text("left alone\n")

    def lines():
        for i, line in enumerate(VCF_X.splitlines()):
            if i == 5:
                # Another job writes next to ours while the run reads its source.
                (here / "s.other_job.frq").write_text("another job\n")
                (here / "s.frq.tmp-otherjob").write_text("another writer's temp file\n")
            yield line

    got = gloci.vcf_stats(lines(), out="s", freq=True)
    assert got == [Path("s.frq")]
    assert (here / "s.het").read_text() == "left alone\n"
    # A re-run with the same prefix reports the same files.
    assert gloci.vcf_stats("x.vcf", out="s", freq=True) == got
    assert gloci.vcf_stats("x.vcf", out=here / "s", freq=True, het=True) == [
        here / "s.frq",
        here / "s.het",
    ]
    # index: only the index it wrote.
    gloci.open(FIXTURES / "a.bed").save("a.bed.gz")
    (here / "a.bed.gz.other").write_text("not an index\n")
    assert gloci.index("a.bed.gz") == [Path("a.bed.gz.tbi")]
    assert gloci.index("a.bed.gz") == [Path("a.bed.gz.tbi")]


def test_cli_quotes_what_only_looks_like_a_placeholder(ab: Path, bash) -> None:
    # I2: a placeholder is known by where it came from, never by its text.
    evil = "<(gloci view x; touch pwned)"
    p = gloci.open("a.bed").intersect(b=evil)
    assert p.cli() == "gloci intersect -a a.bed -b '<(gloci view x; touch pwned)'"
    bash(p.cli(), ab)
    assert not (ab / "pwned").exists()
    named = "<python object 1>"
    (ab / named).write_bytes((ab / "a.bed").read_bytes())
    q = gloci.open(named).merge()
    assert q.cli() == "gloci merge '<python object 1>'"
    run = bash(q.cli(), ab)
    assert run.returncode == 0, run.stderr
    assert run.stdout == q.to_bytes()
    # Real placeholders stay as they are, and recorded command lines hold no NUL.
    rows = gloci.open(bed_rows(ab / "a.bed")).intersect(b=gloci.open("b.bed").view())
    assert rows.cli() == "gloci intersect -a <python object 1> -b <(gloci view b.bed)"
    vcf = gloci.open(VCF_X.splitlines()).view().to_bytes().decode()
    want = f"##glociCommand=gloci view <python object 1>; Version={gloci.__version__}"
    assert want in vcf.splitlines()
    assert "\0" not in vcf


def test_a_record_batch_reader_is_read_once(here: Path) -> None:
    # I3: a RecordBatchReader is an iterator: a second run would read nothing.
    pa = pytest.importorskip("pyarrow")
    batch = pa.record_batch({"chrom": ["chr1"], "start": [1], "end": [2]})
    p = gloci.open(pa.RecordBatchReader.from_batches(batch.schema, [batch]))
    assert p.to_bytes() == b"chr1\t1\t2\n"
    for again in (p.to_bytes, p.count):
        with pytest.raises(gloci.UsageError, match="was read by an earlier run"):
            again()
    # A table is read again at each run.
    t = gloci.open(pa.table({"chrom": ["chr1"], "start": [1], "end": [2]}))
    assert t.to_bytes() == t.to_bytes() == b"chr1\t1\t2\n"
