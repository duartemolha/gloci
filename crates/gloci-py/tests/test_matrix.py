"""The behaviour matrix (P2b Task 10, ruling R2): every step method, on an input of its kind,
- gives the same bytes at `threads` 1, 2 and 8;
- gives empty or header-only output (no records) for an empty input;
- raises `UnsortedError` for an unsorted input;
- raises `DataError` for a truncated BGZF input;
and with `GLOCI_BIN` set, the output (threads) or the error's `str` is the CLI's for `cli()`.

`MATRIX` names each method's input kind and arguments; a row that cannot apply to a method is
in `EXEMPT`, with the reason. Nothing else skips.
"""

from __future__ import annotations

import os
import random
import struct
import subprocess
import tempfile
import zlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pytest

import gloci
from gloci import _native
from parity.outputs import cli_message, content

ROWS = ("threads", "empty", "unsorted", "truncated")
CHROMS = (("chr1", 100_000), ("chr2", 80_000), ("chr3", 60_000))


# Inputs


def bgzf(data: bytes) -> bytes:
    """`data` as BGZF: gzip members of at most 64 KiB, each with its BC extra field, then the
    empty EOF block."""
    out = bytearray()
    blocks = [data[i : i + 0xFF00] for i in range(0, len(data), 0xFF00)] + [b""]
    for block in blocks:
        c = zlib.compressobj(6, zlib.DEFLATED, -15)
        body = c.compress(block) + c.flush()
        size = 18 + len(body) + 8 - 1
        out += b"\x1f\x8b\x08\x04\0\0\0\0\0\xff\x06\0BC\x02\0" + struct.pack("<H", size)
        out += body + struct.pack("<II", zlib.crc32(block), len(block))
    return bytes(out)


def truncated(data: bytes) -> bytes:
    """A BGZF file cut in the middle of its first block: whatever reads it reads the cut."""
    first = struct.unpack("<H", data[16:18])[0] + 1
    return data[: first // 2]


@dataclass
class Kind:
    """One kind of input, as the files of each variant: `plain` (sorted), `bgzf` (the same,
    BGZF-compressed), `truncated`, `unsorted`, `empty` (no records; the header if it has one)."""

    plain: str
    bgzf: str
    truncated: str
    unsorted: str
    empty: str

    def file(self, variant: str) -> str:
        return getattr(self, variant)


def text(lines: list[str]) -> str:
    return "".join(f"{line}\n" for line in lines)


def _fasta(rng: random.Random) -> dict[str, str]:
    return {c: "".join(rng.choice("ACGT") for _ in range(n)) for c, n in CHROMS}


def _bed(rng: random.Random, n: int) -> list[tuple[str, int, int, str, int, str]]:
    rows = []
    for k in range(n):
        c, size = CHROMS[k % 3]
        s = rng.randrange(0, size - 600)
        rows.append((c, s, s + rng.randint(1, 500), f"r{k}", rng.randint(0, 1000), "+-"[k % 2]))
    order = {c: i for i, (c, _) in enumerate(CHROMS)}
    return sorted(rows, key=lambda r: (order[r[0]], r[1], r[2]))


def _unsorted(lines: list[str]) -> list[str]:
    """`lines` (records, sorted) with two far apart swapped, so the input is unsorted from a
    point near its start."""
    out = list(lines)
    out[2], out[len(out) // 2] = out[len(out) // 2], out[2]
    return out


def _sam_reads(rng: random.Random, n: int) -> list[str]:
    reads = []
    cigars = ["50M", "20M5D30M", "10S40M", "25M100N25M", "30M2I18M"]
    for k in range(n):
        c, size = CHROMS[k % 3]
        pos = rng.randrange(1, size - 500)
        cigar = cigars[k % len(cigars)]
        qlen = sum(int(x) for x, op in _cigar_ops(cigar) if op in "MIS=X")
        seq = "".join(rng.choice("ACGT") for _ in range(qlen))
        flag = 16 if k % 4 == 1 else 0
        reads.append(
            (
                c,
                pos,
                f"q{k}\t{flag}\t{c}\t{pos}\t{rng.choice((0, 20, 60))}\t{cigar}\t*\t0\t0\t"
                f"{seq}\t{'I' * qlen}\tRG:Z:g1",
            )
        )
    order = {c: i for i, (c, _) in enumerate(CHROMS)}
    return [r for _, _, r in sorted(reads, key=lambda r: (order[r[0]], r[1]))]


def _cigar_ops(cigar: str) -> list[tuple[str, str]]:
    out, num = [], ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            out.append((num, ch))
            num = ""
    return out


@dataclass
class Fixtures:
    dir: Path
    kinds: dict[str, Kind] = field(default_factory=dict)

    def text_kind(self, name: str, header: list[str], records: list[str]) -> Kind:
        """A text kind: its records sorted, unsorted, none, and BGZF-compressed."""
        stem, ext = name.split(".", 1)

        def put(file: str, lines: list[str]) -> str:
            (self.dir / file).write_text(text(lines))
            return file

        put(name, header + records)
        put(f"{stem}.unsorted.{ext}", header + _unsorted(records))
        put(f"{stem}.empty.{ext}", header)
        data = bgzf((self.dir / name).read_bytes())
        (self.dir / f"{name}.gz").write_bytes(data)
        (self.dir / f"{stem}.trunc.{ext}.gz").write_bytes(truncated(data))
        return Kind(
            name,
            f"{name}.gz",
            f"{stem}.trunc.{ext}.gz",
            f"{stem}.unsorted.{ext}",
            f"{stem}.empty.{ext}",
        )


def build(d: Path) -> Fixtures:
    rng = random.Random(10)
    fx = Fixtures(d)
    (d / "g.genome").write_text("".join(f"{c}\t{n}\n" for c, n in CHROMS))

    # FASTA, indexed (faidx and getfasta need it), with BGZF and truncated forms
    seqs = _fasta(rng)
    fa = "".join(
        f">{c}\n" + "".join(s[i : i + 60] + "\n" for i in range(0, len(s), 60))
        for c, s in seqs.items()
    )
    (d / "ref.fa").write_text(fa)
    (d / "empty.fa").write_text("")
    data = bgzf(fa.encode())
    (d / "ref.fa.gz").write_bytes(data)
    (d / "ref.trunc.fa.gz").write_bytes(truncated(data))
    for f in ("ref.fa", "ref.fa.gz"):
        gloci.index(d / f)
    (d / "empty.fa.fai").write_text("")  # gloci index refuses a FASTA with no sequence
    # the truncated FASTA keeps the whole file's index, so a reader reads into the cut
    for ext in (".fai", ".gzi"):
        (d / f"ref.trunc.fa.gz{ext}").write_bytes((d / f"ref.fa.gz{ext}").read_bytes())
    fx.kinds["fasta"] = Kind("ref.fa", "ref.fa.gz", "ref.trunc.fa.gz", "ref.fa", "empty.fa")

    bed = ["\t".join(map(str, r)) for r in _bed(rng, 3000)]
    fx.kinds["bed"] = fx.text_kind("a.bed", [], bed)
    (d / "b.bed").write_text("".join("\t".join(map(str, r)) + "\n" for r in _bed(rng, 1500)))

    pairs = []
    for r in _bed(rng, 1000):
        c2, size2 = CHROMS[rng.randrange(3)]
        s2 = rng.randrange(0, size2 - 600)
        pairs.append("\t".join(map(str, (*r[:3], c2, s2, s2 + 300, r[3], r[4], r[5], "+"))))
    fx.kinds["bedpe"] = fx.text_kind("a.bedpe", [], pairs)
    b_pairs = []
    for r in _bed(rng, 500):
        end2 = (r[0], r[2] + 100, r[2] + 400)
        b_pairs.append("\t".join(map(str, (*r[:3], *end2, r[3], r[4], r[5], "-"))))
    (d / "b.bedpe").write_text(text(b_pairs))

    def vcf_header(samples: list[str]) -> list[str]:
        return [
            "##fileformat=VCFv4.2",
            *(f"##contig=<ID={c},length={n}>" for c, n in CHROMS),
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
            "\t".join(
                ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples]
            ),
        ]

    def vcf_records(chroms: list[str], n: int, samples: int) -> list[str]:
        out = []
        for c in chroms:
            seq = seqs[c]
            for pos in sorted(rng.sample(range(10, len(seq) - 10), n)):
                ref = seq[pos - 1]
                alts = [b for b in "ACGT" if b != ref]
                alt = ",".join(alts[:2]) if pos % 7 == 0 else alts[0]
                gts = "\t".join(
                    f"{rng.choice(('0/0', '0/1', '1/1', './.'))}:{rng.randint(1, 60)}"
                    for _ in range(samples)
                )
                out.append(
                    f"{c}\t{pos}\t.\t{ref}\t{alt}\t{rng.randint(1, 99)}\tPASS\t"
                    f"DP={rng.randint(1, 200)}\tGT:DP\t{gts}"
                )
        return out

    fx.kinds["vcf"] = fx.text_kind(
        "calls.vcf", vcf_header(["S1", "S2"]), vcf_records(["chr1", "chr2"], 500, 2)
    )
    (d / "calls3.vcf").write_text(text(vcf_header(["S1", "S2"]) + vcf_records(["chr3"], 200, 2)))
    (d / "other.vcf").write_text(text(vcf_header(["S3"]) + vcf_records(["chr1", "chr2"], 300, 1)))
    # the extra inputs of vcf_concat and vcf_merge, empty, for the empty row
    (d / "calls3.empty.vcf").write_text(text(vcf_header(["S1", "S2"])))
    (d / "other.empty.vcf").write_text(text(vcf_header(["S3"])))

    sam_header = [
        "@HD\tVN:1.6\tSO:coordinate",
        *(f"@SQ\tSN:{c}\tLN:{n}" for c, n in CHROMS),
        "@RG\tID:g1\tSM:s1",
    ]
    reads = _sam_reads(rng, 2000)
    for name, lines in (
        ("reads", sam_header + reads),
        ("reads.unsorted", sam_header + _unsorted(reads)),
        ("reads.empty", sam_header),
    ):
        (d / f"{name}.sam").write_text(text(lines))
        gloci.open(d / f"{name}.sam").save(d / f"{name}.bam")
    (d / "reads.trunc.bam").write_bytes(truncated((d / "reads.bam").read_bytes()))
    gloci.index(d / "reads.bam")  # bedcov reads each region through it
    fx.kinds["bam"] = Kind(
        "reads.bam", "reads.bam", "reads.trunc.bam", "reads.unsorted.bam", "reads.empty.bam"
    )
    (d / "new-header.sam").write_text(text(sam_header + ["@CO\tnew"]))
    return fx


@pytest.fixture(scope="session")
def fixtures(tmp_path_factory: pytest.TempPathFactory) -> Fixtures:
    return build(tmp_path_factory.mktemp("matrix"))


# The table


@dataclass(frozen=True)
class Step:
    """A step method's input kind, its extra positional sources (by kind or file; an extra of
    a kind follows the row's empty variant), and its keyword arguments."""

    kind: str
    extra: tuple[str, ...] = ()
    kwargs: dict[str, Any] = field(default_factory=dict)


def step(kind: str, *extra: str, **kwargs: Any) -> Step:
    return Step(kind, extra, kwargs)


#: Each step method's input and arguments.
MATRIX: dict[str, Step] = {
    "intersect": step("bed", b="b.bed"),
    "subtract": step("bed", b="b.bed"),
    "closest": step("bed", b="b.bed"),
    "map": step("bed", b="b.bed", column="5:sum"),
    "coverage": step("bed", b="b.bed"),
    "pair_intersect": step("bedpe", b="b.bed"),
    "pair_pair": step("bedpe", b="b.bedpe"),
    "stats_jaccard": step("bed", b="b.bed"),
    "stats_fisher": step("bed", b="b.bed", genome="g.genome"),
    "stats_reldist": step("bed", b="b.bed"),
    "merge": step("bed"),
    "partition": step("bed"),
    "complement": step("bed", genome="g.genome"),
    "sort": step("bed"),
    "vcf_concat": step("vcf", "calls3"),
    "vcf_merge": step("vcf", "other"),
    "cat": step("bam", "bam"),
    "depth": step("bam"),
    "samples": step("bam"),
    "view": step("bed"),
    "convert": step("bam"),
    "slop": step("bed", genome="g.genome", both=10),
    "flank": step("bed", genome="g.genome", both=10),
    "shift": step("bed", genome="g.genome", by=5),
    "cluster": step("bed"),
    "spacing": step("bed"),
    "groupby": step("bed", group="1", column="5:sum"),
    "genomecov": step("bed", genome="g.genome"),
    "vcf_filter": step("vcf", include="QUAL>20"),
    "vcf_extract": step("vcf", format="GT"),
    "vcf_split": step("vcf"),
    "vcf_norm": step("vcf", fasta="ref.fa"),
    "vcf_annotate": step("vcf", remove="INFO/DP"),
    "getfasta": step("bed", fasta="ref.fa"),
    "nuc": step("bed", fasta="ref.fa"),
    "sample": step("bed", count=100, seed=1),
    "shuffle": step("bed", genome="g.genome", seed=1),
    "faidx": step("fasta", regions=["chr1:1-100", "chr3:500-700"]),
    "dict": step("fasta"),
    "flagstat": step("bam"),
    "idxstats": step("bam"),
    "head": step("bed", records=100),
    "bedcov": step("bed", files="reads.bam"),
    "maskfasta": step("fasta", b="b.bed"),
    "windows": step("bed", size=100),
    "reheader": step("bam", header="new-header.sam"),
    "igv": step("bed"),
    "links": step("bed"),
    "tag": step("bam", b="b.bed", labels="b"),
}

_EMPTY_FASTA = "an empty FASTA is a DataError (it has no records), on the command line too"
_ANY_ORDER = "reads its input record by record in any order, as the command line does"

#: Rows that cannot apply to a step method, with the reason: (method, row) -> why.
EXEMPT: dict[tuple[str, str], str] = {
    # threads
    ("samples", "threads"): "has no -t: it reads the files' headers only",
    # empty input: output that is not records of the input
    ("complement", "empty"): "the complement of nothing is the whole genome",
    ("genomecov", "empty"): "the genome's coverage: every base at depth 0",
    ("flagstat", "empty"): "a summary: every count is 0",
    ("igv", "empty"): "writes the script's first line whatever the input",
    ("links", "empty"): "writes the head and the end of the page whatever the input",
    ("idxstats", "empty"): "a row per reference of the header, with 0 reads",
    ("samples", "empty"): "lists the header's read groups, which an empty BAM has",
    ("stats_jaccard", "empty"): "a summary row: the statistic of an empty A",
    ("stats_fisher", "empty"): "a summary row: the test of an empty A",
    ("dict", "empty"): _EMPTY_FASTA,
    ("maskfasta", "empty"): _EMPTY_FASTA,
    ("faidx", "empty"): "a region of an empty FASTA names a missing sequence: a DataError, as "
    "on the command line",
    # unsorted input: commands that need no order
    ("sort", "unsorted"): "sorts its input",
    ("view", "unsorted"): _ANY_ORDER,
    ("convert", "unsorted"): _ANY_ORDER,
    ("slop", "unsorted"): _ANY_ORDER,
    ("flank", "unsorted"): _ANY_ORDER,
    ("shift", "unsorted"): _ANY_ORDER,
    ("getfasta", "unsorted"): _ANY_ORDER,
    ("nuc", "unsorted"): _ANY_ORDER,
    ("windows", "unsorted"): _ANY_ORDER,
    ("shuffle", "unsorted"): _ANY_ORDER,
    ("sample", "unsorted"): _ANY_ORDER,
    ("head", "unsorted"): _ANY_ORDER,
    ("igv", "unsorted"): _ANY_ORDER,
    ("links", "unsorted"): _ANY_ORDER,
    ("vcf_filter", "unsorted"): _ANY_ORDER,
    ("vcf_extract", "unsorted"): _ANY_ORDER,
    ("vcf_split", "unsorted"): _ANY_ORDER,
    ("vcf_norm", "unsorted"): _ANY_ORDER,
    ("vcf_annotate", "unsorted"): _ANY_ORDER,
    ("pair_intersect", "unsorted"): "A streams in any order against B held in memory (bedtools "
    "pairtobed)",
    ("pair_pair", "unsorted"): "A streams in any order against B held in memory (bedtools "
    "pairtopair)",
    ("groupby", "unsorted"): "groups runs of adjacent lines with equal keys, in any order "
    "(bedtools groupby)",
    ("bedcov", "unsorted"): "reads each region through the BAM's index, in any order (samtools "
    "bedcov)",
    ("cat", "unsorted"): "copies the BAM files' blocks as they are (samtools cat)",
    ("reheader", "unsorted"): "replaces the header without decoding the records",
    ("tag", "unsorted"): "writes the reads in input order, against features held in memory "
    "(bedtools tag)",
    ("flagstat", "unsorted"): "counts flags, in any order",
    ("idxstats", "unsorted"): "counts reads per reference, in any order",
    ("samples", "unsorted"): "reads the header only",
    ("dict", "unsorted"): "a FASTA has no order to check",
    ("faidx", "unsorted"): "a FASTA has no order to check",
    ("maskfasta", "unsorted"): "a FASTA has no order to check",
}


def pipeline(method: str, fx: Fixtures, variant: str) -> gloci.Pipeline:
    s = MATRIX[method]
    extra = []
    for e in s.extra:
        if e in fx.kinds:
            extra.append(fx.kinds[e].file("empty" if variant == "empty" else "plain"))
        elif "." in e:
            extra.append(e)
        else:  # a file with an empty form: `calls3` is calls3.vcf, or calls3.empty.vcf
            extra.append(f"{e}.empty.vcf" if variant == "empty" else f"{e}.vcf")
    head = gloci.open(fx.kinds[s.kind].file(variant))
    return getattr(head, method)(*extra, **s.kwargs)


def check(method: str, row: str, fx: Fixtures) -> str | None:
    """What is wrong with `row` for `method`, if anything."""
    try:
        return _check(method, row, fx)
    except gloci.GlociError as e:
        return f"{type(e).__name__}: {e}"


def _check(method: str, row: str, fx: Fixtures) -> str | None:
    if row == "threads":
        p = pipeline(method, fx, "bgzf")
        outs = {t: p.to_bytes(threads=t) for t in (1, 2, 8)}
        if not outs[1]:
            return "no output to compare"
        if not outs[1] == outs[2] == outs[8]:
            return "threads 1, 2 and 8 give different bytes"
        return cli_differs(p, out=outs[1])
    if row == "empty":
        p = pipeline(method, fx, "empty")
        n = p.count()
        return None if n == 0 else f"{n} records from an empty input:\n{p.to_bytes()[:300]!r}"
    want = gloci.UnsortedError if row == "unsorted" else gloci.DataError
    variant = "unsorted" if row == "unsorted" else "truncated"
    p = pipeline(method, fx, variant)
    try:
        out = p.to_bytes()
    except want as e:
        return cli_differs(p, message=str(e))
    except gloci.GlociError as e:
        return f"{type(e).__name__} rather than {want.__name__}: {e}"
    return f"no {want.__name__}; output {out[:200]!r}"


def cli_differs(
    p: gloci.Pipeline, out: bytes | None = None, message: str | None = None
) -> str | None:
    """With `GLOCI_BIN` set, how the command line `p.cli()` differs from Python: its output
    from `out`, or its error from `message` (D90: `str()` is the CLI's message without
    `error: `). The lines that record the command line are left out (`parity.outputs`)."""
    binary = os.environ.get("GLOCI_BIN")
    if not binary:
        return None
    with tempfile.TemporaryDirectory() as bin_dir:
        os.symlink(os.path.abspath(binary), os.path.join(bin_dir, "gloci"))
        env = dict(os.environ, PATH=f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
        done = subprocess.run(["sh", "-c", p.cli()], env=env, capture_output=True)
    if message is None and (done.returncode != 0 or content(done.stdout) != content(out or b"")):
        return f"the CLI ({p.cli()}) gives other bytes, or exits {done.returncode}"
    if message is not None and (done.returncode != 1 or cli_message(done.stderr) != message):
        return f"the CLI ({p.cli()}) exits {done.returncode}: {done.stderr.decode()!r}"
    return None


@pytest.mark.parametrize("method", sorted(MATRIX))
def test_matrix(method: str, fixtures: Fixtures, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(fixtures.dir)
    problems = []
    for row in ROWS:
        if (method, row) in EXEMPT:
            continue
        why = check(method, row, fixtures)
        if why is not None:
            problems.append(f"{row}: {why}")
    assert not problems, "\n".join(problems)


def test_matrix_covers_every_binding() -> None:
    methods = {b["method"] for b in _native.bindings()}
    assert set(MATRIX) == methods
    for (method, row), why in EXEMPT.items():
        assert method in methods and row in ROWS and why, (method, row)
    # threads cannot apply exactly where the command has no -t
    no_threads = {b["method"] for b in _native.bindings() if not b["threads"]}
    assert {m for m, r in EXEMPT if r == "threads"} == no_threads
