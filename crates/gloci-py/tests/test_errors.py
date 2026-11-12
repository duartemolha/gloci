"""Typed exceptions carrying the CLI's messages (§8.5, D90)."""

from __future__ import annotations

import pickle
from pathlib import Path

import pytest

import gloci

UNSORTED = "chr1\t10\t20\nchr1\t50\t60\nchr1\t30\t40\n"


@pytest.fixture
def here(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    monkeypatch.chdir(tmp_path)
    return tmp_path


def cli_message(done) -> str:
    """The CLI's error message: its stderr without `error: `."""
    assert done.returncode == 1, done.stderr
    text = done.stderr.decode()
    assert text.startswith("error: "), text
    return text[len("error: ") :].rstrip("\n")


def no_nul(e: BaseException) -> None:
    for value in (str(e), *(repr(v) for v in vars(e).values())):
        assert "\0" not in value and "gloci:" not in value, value


def test_usage_error_is_eager(here: Path) -> None:
    p = gloci.open("a.bed")  # never read: the errors come from the step call
    with pytest.raises(gloci.UsageError) as e:
        p.merge(fooo=1)
    assert isinstance(e.value, ValueError) and isinstance(e.value, gloci.GlociError)
    assert "merge() got an unexpected keyword argument 'fooo'" in str(e.value)
    assert "did you mean" in str(e.value)

    with pytest.raises(gloci.UsageError) as e:
        p.merge(max_gap="x")
    assert "max_gap=" in str(e.value)
    assert "--max-gap" not in str(e.value) and "-d" not in str(e.value).split()

    with pytest.raises(gloci.UsageError, match=r"output is set by save\(path\)"):
        p.merge(output="x.bed")
    with pytest.raises(gloci.UsageError, match="threads belongs to the terminal method"):
        p.intersect(b="b.bed", threads=2)
    with pytest.raises(gloci.UsageError, match="b="):
        p.intersect()
    with pytest.raises(TypeError, match="positional"):
        p.intersect("b.bed")


def test_needs_file_step_after_a_step_is_eager(here: Path) -> None:
    gloci.open("ref.fa").faidx(region="chr1:1-10")  # a file input is fine
    with pytest.raises(gloci.UsageError) as e:
        gloci.open("ref.fa").view().faidx(region="chr1:1-10")
    assert str(e.value) == "faidx needs a file input; save() the previous step to a file and open it"


def test_terminal_options_are_checked_before_running(here: Path) -> None:
    (here / "a.bed").write_text("chr1\t10\t20\n")
    with pytest.raises(gloci.UsageError) as e:
        gloci.open("a.bed").merge().save("o.bed.gz", compression_level=99)
    assert "compression_level=" in str(e.value) and "--compression-level" not in str(e.value)
    assert not list(here.glob("o.bed.gz*"))
    with pytest.raises(gloci.UsageError, match="threads="):
        gloci.open("a.bed").merge().to_bytes(threads=0)


def test_data_error_attrs(here: Path, cli) -> None:
    (here / "bad.bed").write_text("chr1\t10\t20\nchr1\t30\t40\nchr1\tx\t50\n")
    with pytest.raises(gloci.DataError) as e:
        gloci.open("bad.bed").merge().to_bytes()
    assert type(e.value) is gloci.DataError
    assert e.value.path == "bad.bed"
    assert e.value.line == 3
    assert str(e.value) == cli_message(cli(["merge", "bad.bed"], here))
    assert str(e.value) == 'bad.bed line 3: start "x" is not an integer'
    no_nul(e.value)


def test_unsorted_error_attrs(here: Path, cli) -> None:
    (here / "unsorted.bed").write_text(UNSORTED)
    with pytest.raises(gloci.UnsortedError) as e:
        gloci.open("unsorted.bed").merge().to_bytes()
    err = e.value
    assert isinstance(err, gloci.DataError)
    assert err.path == "unsorted.bed" and err.line == 3
    assert err.records == ("chr1\t50\t60", "chr1\t30\t40")
    assert err.lines == (2, 3)
    assert "gloci sort" in err.hint
    assert str(err) == cli_message(cli(["merge", "unsorted.bed"], here))
    assert str(err).endswith("\nhint: " + err.hint)
    no_nul(err)
    # The attributes survive pickling (multiprocessing).
    back = pickle.loads(pickle.dumps(err))
    assert type(back) is gloci.UnsortedError and str(back) == str(err)
    assert (back.path, back.line, back.records, back.lines, back.hint) == (
        err.path,
        err.line,
        err.records,
        err.lines,
        err.hint,
    )


def test_error_after_a_pipe_names_the_stream_by_its_display_name(here: Path) -> None:
    (here / "unsorted.bed").write_text(UNSORTED)
    with pytest.raises(gloci.UnsortedError) as e:
        gloci.open("unsorted.bed").view().merge().to_bytes()
    assert e.value.path == "<view output>"
    assert str(e.value).startswith("<view output> is not sorted: line 3 (chr1\t30\t40)")
    no_nul(e.value)
    # The earliest step's own error wins over what it causes downstream.
    (here / "bad.bed").write_text("chr1\t10\t20\nchr1\tx\t50\n")
    with pytest.raises(gloci.DataError) as e:
        gloci.open("bad.bed").view().merge().save("out.bed")
    assert (e.value.path, e.value.line) == ("bad.bed", 2)
    assert not list(here.glob("out.bed*"))
    no_nul(e.value)


def test_open_error(here: Path, cli) -> None:
    with pytest.raises(gloci.DataError) as e:
        gloci.open("missing.bed").merge().to_bytes()
    assert type(e.value) is gloci.DataError
    assert e.value.path == "missing.bed" and e.value.line is None
    assert str(e.value).startswith("cannot open missing.bed: ")
    assert str(e.value) == cli_message(cli(["merge", "missing.bed"], here))


def test_open_takes_paths_only_for_now(here: Path) -> None:
    with pytest.raises(TypeError, match="path"):
        gloci.open(42)  # type: ignore[arg-type]
