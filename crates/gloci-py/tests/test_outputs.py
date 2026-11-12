"""Records, iteration, count() and Arrow output (D89): what a pipeline's output reads back as."""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

import gloci

from conftest import REPO

GOLDEN = REPO / "crates" / "gloci" / "tests" / "golden"
#: A BAM of 14 reads: mapped, secondary, supplementary, and two unmapped (one placed).
BAM = GOLDEN / "view-bam-to-sam-vs-samtools" / "input.bam"

BED6 = "chr1\t10\t20\tn1\t5\t+\nchr1\t30\t40\tn2\t7\t-\n"
INDELS_VCF = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100000>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "chr1\t100\t.\tA\tACGT\t50\tPASS\t.\n"
    "chr1\t200\trs2\tAG\tA\t.\tPASS\t.\n"
)
GFF3 = "##gff-version 3\nchr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=g1\n"
DP_VCF = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100000>\n"
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
    "chr1\t100\trs1\tA\tC\t50\tPASS\tDP=40\tGT\t0/1\t1/1\n"
    "chr1\t200\t.\tA\tACGT\t.\t.\tDP=16\tGT\t0/0\t./.\n"
)


@pytest.fixture
def here(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    monkeypatch.chdir(tmp_path)
    return tmp_path


def records_of(text: bytes, prefixes: tuple[bytes, ...] = (b"#",)) -> list[str]:
    """The lines of `text` that are not blank and do not start with one of `prefixes`."""
    return [
        line.decode()
        for line in text.split(b"\n")
        if line.strip() and not line.startswith(prefixes)
    ]


# ---- Records ---------------------------------------------------------------------------------


def test_iteration_records(here: Path) -> None:
    (here / "x.bed").write_text(BED6)
    recs = list(gloci.open("x.bed"))
    assert [type(r) for r in recs] == [gloci.Record, gloci.Record]
    r = recs[0]
    assert (r.chrom, r.start, r.end) == ("chr1", 10, 20)
    assert r.fields == ("chr1", "10", "20", "n1", "5", "+")
    assert r.line == "chr1\t10\t20\tn1\t5\t+"
    assert str(r) == r.line
    # A sequence over its fields.
    assert len(r) == 6
    assert (r[0], r[3], r[-1]) == ("chr1", "n1", "+")
    assert r[3:5] == ("n1", "5")
    assert list(r) == list(r.fields)
    assert "n1" in r
    with pytest.raises(IndexError):
        r[6]
    # Frozen.
    with pytest.raises(AttributeError):
        r.start = 5  # type: ignore[misc]
    assert recs[1] == list(gloci.open("x.bed"))[1] and recs[0] != recs[1]

    # VCF: the footprint, as --vcf-span says (allele by default): an insertion is a point.
    (here / "indels.vcf").write_text(INDELS_VCF)
    ins, dele = gloci.open("indels.vcf")
    assert (ins.chrom, ins.start, ins.end) == ("chr1", 100, 100)
    assert ins.start == ins.end
    assert ins.fields[:5] == ("chr1", "100", ".", "A", "ACGT")
    assert (dele.start, dele.end) == (200, 201)
    # The step's own --vcf-span: REF's bases.
    ref_ins, ref_del = gloci.open("indels.vcf").view(vcf_span="ref")
    assert (ref_ins.start, ref_ins.end) == (99, 100)
    assert (ref_del.start, ref_del.end) == (199, 201)

    # GFF: column 4 is 1-based.
    (here / "x.gff3").write_text(GFF3)
    (g,) = gloci.open("x.gff3")
    assert (g.chrom, g.start, g.end) == ("chr1", 100, 200)
    assert g.fields[3:5] == ("101", "200")


def test_iteration_alignments(here: Path) -> None:
    # BAM output reads back as SAM text: every record, unmapped ones included.
    recs = list(gloci.open(BAM))
    assert len(recs) == 14
    first = recs[0]
    assert first.line.startswith("r_plain\t0\tchr1\t101\t60\t10M\t")
    assert (first.chrom, first.start, first.end) == ("chr1", 100, 110)
    placed, unplaced = recs[10], recs[13]
    assert placed.fields[0] == "r_unmapped_placed"
    assert (placed.chrom, placed.start, placed.end) == ("chr1", None, None)
    assert unplaced.fields[:3] == ("r_unmapped_unplaced", "4", "*")
    assert (unplaced.chrom, unplaced.start, unplaced.end) == (None, None, None)


def test_iteration_generic_tsv(here: Path) -> None:
    (here / "calls.vcf").write_text(DP_VCF)
    p = gloci.open("calls.vcf").vcf_extract(info="DP")
    text = p.to_bytes().decode()
    assert text.startswith("CHROM\tPOS\tID\tREF\tALT\tDP\n")
    recs = list(p)
    # The header row names the columns; it is not a record.
    assert [r.line for r in recs] == text.splitlines()[1:]
    assert not any(r.line.startswith(("#", "CHROM")) for r in recs)
    assert all((r.chrom, r.start, r.end) == (None, None, None) for r in recs)
    assert recs[0].fields == ("chr1", "100", "rs1", "A", "C", "40")
    assert p.count() == len(recs) == 2

    # A table with no header row: its rows are records with no coordinates either.
    p = gloci.open(BAM).idxstats()
    rows = list(p)
    assert [r.line for r in rows] == ["chr1\t100000\t10\t1", "chr2\t100000\t2\t0", "*\t0\t0\t1"]
    # Read as BED, the first row would be chr1:100000-10, which is no interval.
    assert all((r.chrom, r.start, r.end) == (None, None, None) for r in rows)
    assert p.count() == 3


def test_iteration_errors(here: Path) -> None:
    (here / "bad.bed").write_text("chr1\t1\t5\nchr1\t2\t6\nchr1\tx\t9\n")
    it = iter(gloci.open("bad.bed"))
    with pytest.raises(gloci.DataError) as e:
        list(it)
    assert (e.value.path, e.value.line) == ("bad.bed", 3)
    # The iterator has ended.
    with pytest.raises(StopIteration):
        next(it)

    # An exception a source raises is raised from the iteration, the same object.
    boom = RuntimeError("boom")

    def rows():
        yield ("chr1", 1, 2)
        raise boom

    with pytest.raises(RuntimeError) as raised:
        list(gloci.open(rows()))
    assert raised.value is boom


def test_records_as_a_source(here: Path) -> None:
    (here / "x.bed").write_text(BED6)
    recs = list(gloci.open("x.bed").slop(b=5, g="hg38"))
    assert gloci.open(recs).to_bytes() == gloci.open("x.bed").slop(b=5, g="hg38").to_bytes()
    with pytest.raises(TypeError, match="Record"):
        gloci.open(recs[0])


# ---- count() ---------------------------------------------------------------------------------


def write_fasta(here: Path) -> None:
    (here / "ref.fa").write_text(">chr1\nACGTACGTAC\n>chr2\nGGGGCCCCAA\n")
    (here / "ref.fa.fai").write_text("chr1\t10\t6\t10\t11\nchr2\t10\t23\t10\t11\n")
    (here / "r.bed").write_text("chr1\t0\t5\nchr2\t2\t6\n")


def test_count_matches_iteration(ab: Path) -> None:
    (ab / "calls.vcf").write_text(DP_VCF)
    (ab / "x.gff3").write_text(GFF3)
    (ab / "noisy.bed").write_text(
        "track name=x\n# a comment\nchr1\t1\t5\n\nbrowser position chr1\nchr1\t2\t6\n# late\n"
    )
    write_fasta(ab)
    fisher_genome = ab / "g.txt"
    fisher_genome.write_text("chr1\t100000\nchr2\t100000\nchr3\t100000\n")
    pipelines = {
        "bed": (gloci.open("a.bed").intersect(b="b.bed"), None),
        "noisy bed": (gloci.open("noisy.bed"), 2),
        "vcf": (gloci.open("calls.vcf"), 2),
        "gff": (gloci.open("x.gff3"), 1),
        "fasta": (gloci.open("r.bed").getfasta(fasta="ref.fa"), 2),
        "jaccard": (gloci.open("a.bed").stats_jaccard(b="b.bed"), 1),
        "fisher": (gloci.open("a.bed").stats_fisher(b="b.bed", g=str(fisher_genome)), 1),
        "vcf extract": (gloci.open("calls.vcf").vcf_extract(info="DP"), 2),
        "bam": (gloci.open(BAM), 14),
        "idxstats": (gloci.open(BAM).idxstats(), 3),
        "flagstat": (gloci.open(BAM).flagstat(), None),
    }
    for name, (p, want) in pipelines.items():
        recs = list(p)
        assert p.count() == len(recs), name
        if want is not None:
            assert len(recs) == want, name
    # FASTA: a record is a sequence, its `>` line and its sequence lines.
    seqs = list(pipelines["fasta"][0])
    assert seqs[0].line == ">chr1:0-5\nACGTA"
    assert seqs[0].fields == ("chr1:0-5", "ACGTA")
    assert seqs[0].chrom is None


# ---- Arrow -----------------------------------------------------------------------------------


def schema_of(t) -> list[str]:  # type: ignore[no-untyped-def]
    return [f"{f.name}:{f.type}" for f in t.schema]


def test_to_arrow_schema(ab: Path) -> None:
    pa = pytest.importorskip("pyarrow")
    p = gloci.open("a.bed").intersect(b="b.bed", mode="count")
    t = p.to_arrow()
    assert isinstance(t, pa.Table)
    assert schema_of(t) == [
        "chrom:string",
        "start:int64",
        "end:int64",
        "col4:string",
        "col5:int64",
        "col6:string",
        "col7:int64",
    ]
    lines = p.to_bytes().decode().splitlines()
    assert t.num_rows == len(lines) == 1200
    first = lines[0].split("\t")
    assert t.slice(0, 1).to_pylist()[0] == {
        "chrom": first[0],
        "start": int(first[1]),
        "end": int(first[2]),
        "col4": first[3],
        "col5": int(first[4]),
        "col6": first[5],
        "col7": int(first[6]),
    }

    # VCF columns (D89): qual is float64, with `.` as null.
    (ab / "calls.vcf").write_text(DP_VCF)
    v = gloci.open("calls.vcf").to_arrow()
    assert v.column_names == [
        "chrom", "start", "end", "pos", "id", "ref", "alt", "qual", "filter", "info", "format",
        "s1", "s2",
    ]
    assert v.schema.field("qual").type == pa.float64()
    assert v.schema.field("pos").type == pa.int64()
    assert v.column("qual").to_pylist() == [50.0, None]
    assert v.column("start").to_pylist() == [99, 200]
    assert v.column("end").to_pylist() == [100, 200]
    assert v.column("id").to_pylist() == ["rs1", "."]

    # GFF columns (D89).
    (ab / "x.gff3").write_text(GFF3)
    g = gloci.open("x.gff3").to_arrow()
    assert g.column_names == [
        "chrom", "start", "end", "source", "feature", "score", "strand", "frame", "attributes",
    ]
    assert g.to_pylist()[0]["start"] == 100

    # A header row names a table's columns, and types are inferred: int64, float64, string.
    j = gloci.open("a.bed").stats_jaccard(b="b.bed").to_arrow()
    assert schema_of(j) == [
        "intersection:int64", "union:int64", "jaccard:double", "n_intersections:int64",
    ]
    assert j.num_rows == 1

    # names= overrides them, with exactly one name per column.
    names = ["c", "s", "e", "name", "score", "strand", "n"]
    assert p.to_arrow(names=names).column_names == names
    with pytest.raises(gloci.UsageError, match="names= has 2 names, but the output has 7 columns"):
        p.to_arrow(names=["a", "b"])

    # Threads do not change the table.
    assert p.to_arrow(threads=1).equals(p.to_arrow(threads=4))


def test_to_pandas_and_polars(ab: Path) -> None:
    pytest.importorskip("pyarrow")
    p = gloci.open("a.bed").intersect(b="b.bed", mode="count")
    want = p.to_arrow()
    if (pd := _try_import("pandas")) is not None:
        df = p.to_pandas()
        assert isinstance(df, pd.DataFrame)
        assert list(df.columns) == want.column_names
        assert df["start"].tolist() == want.column("start").to_pylist()
    if (pl := _try_import("polars")) is not None:
        df = p.to_polars()
        assert isinstance(df, pl.DataFrame)
        assert df.columns == want.column_names
        assert df["col7"].to_list() == want.column("col7").to_pylist()


def _try_import(name: str):  # type: ignore[no-untyped-def]
    try:
        return __import__(name)
    except ImportError:
        return None


def test_arrow_capsule_protocol(ab: Path) -> None:
    pa = pytest.importorskip("pyarrow")
    p = gloci.open("a.bed").intersect(b="b.bed", mode="count")
    s = p.to_arrow_stream()
    assert isinstance(s, gloci.ArrowStream)
    capsule = s.__arrow_c_stream__()
    assert type(capsule).__name__ == "PyCapsule"
    del capsule  # an unread stream is released with its capsule
    reader = pa.RecordBatchReader.from_stream(s)
    assert reader.read_all().equals(p.to_arrow())
    # The stream can be exported again (each export is a new stream over the same buffers).
    assert pa.table(s).num_rows == 1200


@pytest.mark.parametrize(
    ("method", "blocked", "extra"),
    [
        ("to_arrow", "pyarrow", "arrow"),
        ("to_polars", "polars", "polars"),
        ("to_pandas", "pandas", "pandas"),
        ("to_pandas", "pyarrow", "pandas"),
    ],
)
def test_missing_extra_message(
    ab: Path, monkeypatch: pytest.MonkeyPatch, method: str, blocked: str, extra: str
) -> None:
    monkeypatch.setitem(sys.modules, blocked, None)  # `import <blocked>` now fails
    with pytest.raises(ImportError, match=re.escape(f"pip install gloci[{extra}]")):
        getattr(gloci.open("a.bed"), method)()


def test_polars_without_pyarrow(ab: Path) -> None:
    pytest.importorskip("polars")
    # A fresh interpreter in which pyarrow cannot be imported, whether or not it is installed.
    code = """
import sys
sys.modules["pyarrow"] = None
import polars as pl
import gloci
p = gloci.open("a.bed").intersect(b="b.bed", mode="count")
df = p.to_polars()
assert df.columns == ["chrom", "start", "end", "col4", "col5", "col6", "col7"], df.columns
assert df.dtypes == [pl.String, pl.Int64, pl.Int64, pl.String, pl.Int64, pl.String, pl.Int64]
assert df.height == 1200, df.height
# polars.from_arrow takes the stream itself (as a struct Series, one field per column).
s = pl.from_arrow(p.to_arrow_stream())
assert s.struct.unnest().equals(df)
assert sys.modules["pyarrow"] is None
print("ok")
"""
    done = subprocess.run(
        [sys.executable, "-c", code], cwd=ab, capture_output=True, text=True, check=False
    )
    assert done.returncode == 0, done.stderr
    assert done.stdout.strip() == "ok"


def test_arrow_of_alignments_and_tables(here: Path) -> None:
    pytest.importorskip("pyarrow")
    t = gloci.open(BAM).to_arrow()
    assert t.column_names[:13] == [
        "chrom", "start", "end", "qname", "flag", "pos", "mapq", "cigar", "rnext", "pnext",
        "tlen", "seq", "qual",
    ]
    assert t.num_rows == 14
    assert t.column("chrom").to_pylist()[-1] is None
    i = gloci.open(BAM).idxstats().to_arrow()
    assert schema_of(i) == ["col1:string", "col2:int64", "col3:int64", "col4:int64"]


# ---- What the readers pass over, empty outputs, and stopping early ---------------------------


def test_lines_a_reader_passes_over_are_records(here: Path) -> None:
    # convert --bedpe writes `.  -1  -1` for an absent end: BED's reader cannot place such a
    # line, but it is still a record, with no coordinates.
    p = gloci.open(GOLDEN / "bamtobed-bedpe-vs-bedtools" / "input.bam").convert(bedpe=True)
    lines = p.to_bytes().decode().splitlines()
    recs = list(p)
    assert [r.line for r in recs] == lines
    assert p.count() == len(lines)
    absent = [r for r in recs if r.line.startswith(".\t-1\t-1")]
    assert absent and all((r.chrom, r.start, r.end) == (None, None, None) for r in absent)
    placed = [r for r in recs if not r.line.startswith(".")]
    assert all(r.start == int(r.fields[1]) for r in placed)
    if _try_import("pyarrow") is not None:
        t = p.to_arrow()
        assert t.num_rows == len(lines)
        assert t.column("chrom").null_count == len(absent)
        assert t.column("col7").to_pylist() == [r.fields[6] for r in recs]


def test_empty_outputs(here: Path) -> None:
    (here / "empty.bed").write_text("")
    (here / "header.vcf").write_text(INDELS_VCF.split("chr1\t100")[0])
    for p in (gloci.open("empty.bed"), gloci.open("header.vcf")):
        assert list(p) == []
        assert p.count() == 0
    if _try_import("pyarrow") is not None:
        assert gloci.open("empty.bed").to_arrow().column_names == ["chrom", "start", "end"]
        v = gloci.open("header.vcf").to_arrow()
        assert v.num_rows == 0
        assert v.column_names[:8] == ["chrom", "start", "end", "pos", "id", "ref", "alt", "qual"]


def threads_now() -> int | None:
    from test_sources import os_threads

    return os_threads()


def test_break_stops_the_steps(here: Path) -> None:
    from test_sources import assert_threads_back_to

    (here / "g.txt").write_text("chr1\t10000000\n")
    before = threads_now()
    # 10 million windows: far more than the pipes and the batches hold.
    p = gloci.windows(g="g.txt", w=1).slop(b=1, g="g.txt")
    for r in p:
        assert (r.chrom, r.start, r.end) == ("chr1", 0, 2)
        break
    assert_threads_back_to(before)

    # A generator source is read no further once the loop is left.
    taken = [0]

    def rows():
        for i in range(10_000_000):
            taken[0] += 1
            yield ("chr1", i, i + 1)

    it = iter(gloci.open(rows()))
    assert next(it).start == 0
    del it
    seen = taken[0]
    assert seen < 10_000_000
    assert_threads_back_to(before)
    assert taken[0] == seen


def test_keyboard_interrupt_in_iteration(here: Path) -> None:
    from test_sources import assert_threads_back_to

    before = threads_now()

    def rows():
        for i in range(50_000):
            yield ("chr1", i, i + 1)
        raise KeyboardInterrupt

    got = 0
    with pytest.raises(KeyboardInterrupt):
        for _ in gloci.open(rows()):
            got += 1
    assert got < 50_000  # the run stopped at once, not after the records before
    assert_threads_back_to(before)


def test_iterators_are_independent_runs(ab: Path) -> None:
    p = gloci.open("a.bed").merge()
    a, b = iter(p), iter(p)
    first_a = [next(a) for _ in range(3)]
    first_b = [next(b) for _ in range(3)]
    assert first_a == first_b
    assert len(list(a)) + 3 == p.count()
    b.close()
    assert list(b) == []


# ---- Shapes that follow the flags given (fix round 1, I1 and I2) ----------------------------


def assert_columns_line_up(t, lines: list[str]) -> None:  # type: ignore[no-untyped-def]
    """Each row of the table `t` (a table's `col<N>` columns) holds its line's fields, in
    order, and nulls past its last field."""
    for row, line in zip(t.to_pylist(), lines, strict=True):
        fields = line.split("\t")
        for k in range(t.num_columns):
            v = row[f"col{k + 1}"]
            if k >= len(fields):
                assert v is None, (line, k)
            elif isinstance(v, float):
                assert v == float(fields[k]), (line, k)
            else:
                assert str(v) == fields[k], (line, k)


def test_tables_do_not_read_back_as_intervals(ab: Path) -> None:
    # coverage --hist: A's lines with histogram columns, then the `all` summary lines.
    p = gloci.open("a.bed").coverage(b="b.bed", hist=True)
    lines = p.to_bytes().decode().splitlines()
    assert any(line.startswith("all\t") for line in lines)
    recs = list(p)
    assert [r.line for r in recs] == lines
    assert p.count() == len(lines)
    assert not any(r.chrom == "all" and r.start is not None for r in recs)
    assert all((r.chrom, r.start, r.end) == (None, None, None) for r in recs)
    if _try_import("pyarrow") is not None:
        t = p.to_arrow()
        assert t.column_names[0] == "col1"
        assert_columns_line_up(t, lines)

    # groupby keyed on other columns than chrom, start and end: `chr2  5310  450` is no
    # interval.
    g = gloci.open("a.bed").groupby(g="1", c=["5:sum", "5:count"])
    lines = g.to_bytes().decode().splitlines()
    rows = list(g)
    assert [r.line for r in rows] == lines and len(rows) == 3
    assert all((r.chrom, r.start, r.end) == (None, None, None) for r in rows)
    # Keyed on chrom, start and end (then more), it writes BED records.
    k = list(gloci.open("a.bed").groupby(g="1,2,3,4", c="5:sum"))
    assert all((r.chrom, r.start, r.end) == (r[0], int(r[1]), int(r[2])) for r in k)


def test_groupby_header_names_the_columns(here: Path) -> None:
    rows = "chr1\t1\t2\ta\t3\nchr1\t5\t6\tb\t4\nchr2\t1\t2\tc\t1\n"
    for header in ("chrom\tstart\tend\tname\tscore", "#chrom\tstart\tend\tname\tscore"):
        (here / "t.tsv").write_text(f"{header}\n{rows}")
        p = gloci.open("t.tsv").groupby(header=True, g="chrom", c="score:sum")
        lines = p.to_bytes().decode().splitlines()
        assert lines[0] == header.split("\t")[0] + "\tscore:sum"
        recs = list(p)
        # The first line is the header, prefixed or not: two groups, two records.
        assert [r.line for r in recs] == lines[1:] == ["chr1\t7", "chr2\t1"]
        assert p.count() == 2
        if _try_import("pyarrow") is not None:
            assert p.to_arrow().column_names == ["chrom", "score:sum"]


def test_ctrl_c_while_iterating(here: Path) -> None:
    # A real Ctrl-C (KeyboardInterrupt in the main thread, as SIGINT gives) while the first
    # record is awaited: sort reads a FIFO, spills runs to tmp_dir, and waits for more input;
    # the interrupt comes once a spill file exists.
    import _thread
    import random
    import threading
    import time

    from test_sources import assert_threads_back_to

    fifo = here / "big.bed"
    os.mkfifo(fifo)
    stop = threading.Event()

    def feed() -> None:
        rng = random.Random(3)
        try:
            with open(fifo, "w") as f:
                # More than sort's 16M budget at once, so it spills ...
                f.writelines(
                    f"chr{rng.randint(1, 3)}\t{s}\t{s + 50}\n"
                    for s in (rng.randrange(0, 100_000_000) for _ in range(1_000_000))
                )
                # ... then a trickle for 10 s at most, so the run cannot end before the
                # interrupt, and a cancel that does not work fails the test rather than hang it.
                end = time.monotonic() + 10
                while not stop.is_set() and time.monotonic() < end:
                    f.write("chr1\t1\t2\n")
                    f.flush()
                    time.sleep(0.002)
        except (BrokenPipeError, OSError):
            pass

    spill = here / "spill"
    spill.mkdir()
    before = threads_now()
    seen: list[Path] = []
    interrupted: list[float] = []

    def interrupt() -> None:
        # As soon as a spill file exists (within 10 s), however loaded the machine is; anyway
        # after that, so that a run that never spills fails the test below.
        deadline = time.monotonic() + 10
        while not seen and time.monotonic() < deadline:
            seen.extend(p for p in spill.rglob("*") if p.is_file())
            time.sleep(0.002)
        interrupted.append(time.monotonic())
        _thread.interrupt_main()

    feeder = threading.Thread(target=feed, daemon=True)
    watcher = threading.Thread(target=interrupt, daemon=True)
    feeder.start()
    watcher.start()
    got = 0
    try:
        with pytest.raises(KeyboardInterrupt):
            for _ in gloci.open("big.bed").sort(tmp_dir=str(spill), max_mem="16M"):
                got += 1
        raised = time.monotonic()
    finally:
        stop.set()
        watcher.join(15)
        feeder.join(15)
    assert seen, "sort did not spill within 10 s"
    assert got == 0
    # The cancel reached the steps at once, well before the feeder's 10 s.
    assert raised - interrupted[0] < 2, raised - interrupted[0]
    assert_threads_back_to(before)
    assert list(spill.iterdir()) == []


def test_iterator_as_a_context_manager(here: Path) -> None:
    from test_sources import assert_threads_back_to

    (here / "g.txt").write_text("chr1\t10000000\n")
    p = gloci.windows(g="g.txt", w=1)
    before = threads_now()
    # `with` closes the iterator when the block ends, however it ends.
    with iter(p) as it:
        assert next(it).start == 0
    assert list(it) == []
    assert_threads_back_to(before)
    with pytest.raises(ZeroDivisionError):
        with iter(p) as it:
            next(it)
            1 / 0  # noqa: B018
    assert list(it) == []
    assert_threads_back_to(before)


def test_track_and_browser_lines_are_not_rows(here: Path) -> None:
    # coverage --hist is a table; --header copies A's track and browser lines into it.
    (here / "a.bed").write_text(
        "track name=x\nbrowser position chr1:1-100\nchr1\t10\t20\nchr1\t30\t40\n"
    )
    (here / "b.bed").write_text("chr1\t15\t35\n")
    p = gloci.open("a.bed").coverage(b="b.bed", hist=True, header=True)
    lines = p.to_bytes().decode().splitlines()
    assert lines[:2] == ["track name=x", "browser position chr1:1-100"]
    assert [r.line for r in p] == lines[2:]
    assert p.count() == len(lines) - 2
