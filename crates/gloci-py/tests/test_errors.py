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
    assert str(e.value) == "merge() got an unexpected keyword argument 'fooo'"
    # A near miss is suggested (within two edits); no name is that close to fooo.
    with pytest.raises(gloci.UsageError) as e:
        p.merge(max_gp=5)
    assert str(e.value).endswith("did you mean 'max_gap'?")

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


def test_open_refuses_what_is_not_a_source(here: Path) -> None:
    with pytest.raises(TypeError, match="path"):
        gloci.open(42)  # type: ignore[arg-type]


def test_sort_check_is_reserved(here: Path) -> None:
    # Every terminal writes the step's output, which sort --check does not: it is CLI-only.
    for value in (True, False):
        with pytest.raises(gloci.UsageError) as e:
            gloci.open("a.bed").sort(check=value)
        assert str(e.value) == (
            "sort(check=True) writes no output; check sortedness with the CLI: "
            "gloci sort --check FILE"
        )
    # Never suggested for a near miss either.
    with pytest.raises(gloci.UsageError) as e:
        gloci.open("a.bed").sort(chek=True)
    assert "check" not in str(e.value).split("'chek'")[1]


def test_save_output_format_is_named_as_its_keyword(here: Path) -> None:
    (here / "a.bed").write_text("chr1\t10\t20\n")
    with pytest.raises(gloci.UsageError) as e:
        gloci.open("a.bed").merge().save("o.bed", output_format="xyz")
    assert str(e.value).startswith("unknown output format 'xyz' for output_format=; use one of:")
    assert not list(here.glob("o.bed*"))


def test_hints_in_keyword_terms(here: Path) -> None:
    # A flag's value follows its `=`; a gloci command line in a hint stays a command line.
    (here / "a.bed").write_text("chr1\t1\t5\n")
    with pytest.raises(gloci.UsageError) as e:
        gloci.open("a.bed").sort(max_mem="4M").to_bytes()
    assert str(e.value) == "max_mem=4194304 bytes is too small; use at least 16M"
    bam = (
        Path(__file__).resolve().parents[2]
        / "gloci/tests/golden/view-bam-to-sam-vs-samtools/input.bam"
    )
    with pytest.raises(gloci.UsageError) as e:
        gloci.open(bam).convert(bedpe=True).to_bytes()
    assert str(e.value).endswith(
        f"sort it by name first: gloci sort -n {bam} -o x.byname.bam"
    ), str(e.value)


# The edges of the keyword and terminal arguments.


def test_save_options_hold_no_nul(here: Path) -> None:
    # M1: a NUL in a terminal option is a UsageError (it never reaches the plan, which refuses
    # its own sentinel with a plain ValueError), and no NUL reaches the message.
    (here / "a.bed").write_text("chr1\t1\t5\n")
    p = gloci.open("a.bed").sort()
    for kw in (
        {"output_format": "\0src:0"},
        {"output_format": "bed\0"},
        {"threads": "\0src:0"},
        {"compression_level": "1\0"},
    ):
        with pytest.raises(gloci.UsageError, match="NUL") as e:
            p.save("o.bed", **kw)  # type: ignore[arg-type]
        no_nul(e.value)
    for terminal in (p.to_bytes, p.count, p.to_arrow_stream):
        with pytest.raises(gloci.UsageError, match="NUL") as e:
            terminal(threads="\0src:0")  # type: ignore[operator]
        no_nul(e.value)
    assert not (here / "o.bed").exists()


def test_a_bytes_value_is_refused(here: Path) -> None:
    # M3: bytes would be str()'d into the path "b'x.bed'"; refuse it, as open(b"...") is refused.
    (here / "a.bed").write_text("chr1\t1\t5\n")
    p = gloci.open("a.bed")
    for v in (b"a.bed", bytearray(b"a.bed"), [b"a.bed"]):
        with pytest.raises(TypeError, match="not (bytes|bytearray)"):
            p.intersect(b=v)  # type: ignore[arg-type]
    # A str still works.
    assert p.intersect(b="a.bed").count() == 1


def test_vcf_merge_step_refuses_file_list(here: Path) -> None:
    # M12: the step's inputs are sources.
    with pytest.raises(gloci.UsageError) as e:
        gloci.open("a.vcf").vcf_merge(file_list="list.txt")
    assert str(e.value) == "pass the files as sources: open(a).vcf_merge(b, c)"


def test_usage_messages_name_the_long_keyword(here: Path) -> None:
    # M13: the run's message names `-c`; the keyword the caller can pass is column=.
    (here / "a.bed").write_text("chr1\t1\t5\tx\n")
    for kw in ({"column": "5"}, {"c": "5"}):
        with pytest.raises(gloci.UsageError) as e:
            gloci.open("a.bed").map(b="a.bed", **kw).to_bytes()
        assert "column=5" in str(e.value) and "-c" not in str(e.value), str(e.value)
