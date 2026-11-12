"""Checked by `mypy --strict` (python.yml), not run by pytest: the stub types real usage."""

from __future__ import annotations

import io
import pathlib
from collections.abc import Iterator

import gloci
from gloci import DataError, GlociError, Pipeline, Record, UnsortedError, UsageError


def sources() -> list[Pipeline]:
    rows: Iterator[tuple[str, int, int]] = iter([("chr1", 1, 5)])
    return [
        gloci.open("a.bed"),
        gloci.open(pathlib.Path("a.bed"), format="bed"),
        gloci.open(io.BytesIO(b"chr1\t1\t5\n")),
        gloci.open(["chr1\t1\t5"]),
        gloci.open([("chr1", 1, 5, None)]),
        gloci.open(rows),
        gloci.open(gloci.open("a.bed").sort()),
        gloci.windows(genome="hg38", size=1_000_000),
        gloci.random(genome="hg38", length=100, count=10, seed=7),
    ]


def pipeline() -> None:
    a = gloci.open("calls.vcf.gz")
    b = gloci.open("targets.bed").slop(genome="hg38", both=10)
    several: list[gloci.Source] = ["targets.bed", b]  # a mixed list needs its type
    p = (
        a.intersect(b=several, min_frac_a=0.5, reciprocal=True, mode="any")
        .merge(max_gap=50)
        .sort()
    )
    p.save("out.bed.gz", write_index=True, threads=2)
    data: bytes = p.to_bytes()
    n: int = p.count(threads=1)
    cli: str = p.cli()
    for r in p:
        c: str | None = r.chrom
        s: int | None = r.start
        f: tuple[str, ...] = r.fields
        line: str = r.line
        _ = (c, s, f, line, r[0], len(r), data, n, cli)
    with iter(p) as records:
        first: Record = next(records)
        records.close()
        _ = first
    p.cat(gloci.open("x.bed"), gloci.open("y.bed"))


def frames() -> None:
    p = gloci.open("a.bed")
    table = p.to_arrow(names=["c", "s", "e"])
    df = p.to_polars()
    pdf = p.to_pandas(threads=2)
    stream: gloci.ArrowStream = p.to_arrow_stream()
    cols: list[str] = stream.column_names
    rows: int = stream.num_rows
    _ = (table.num_rows, df.height, pdf.shape, cols, rows)


def functions() -> None:
    written: list[pathlib.Path] = gloci.vcf_stats("a.vcf", out="stats", freq=True)
    gloci.vcf_extract_matrix("a.vcf", out="m")
    gloci.vcf_compare("a.vcf", gloci.open("b.vcf"), out="cmp", site_discordance=True)
    gloci.index(pathlib.Path("a.bed.gz"))
    _ = written


def errors() -> None:
    try:
        gloci.open("a.bed").merge(max_gap=1).count()
    except UnsortedError as e:
        hint: str | None = e.hint
        _ = hint
    except DataError as e:
        path = e.path
        line = e.line
        _ = (path, line)
    except UsageError as e:
        value: ValueError = e
        _ = value
    except GlociError:
        raise
