"""Outputs as the golden harness compares them (`gloci_testkit::differential`), for the
parity tests: the lines that record the command line normalised, BAM decoded, and the CLI's
error message read from its stderr."""

from __future__ import annotations

import gzip
import re
from pathlib import Path

__all__ = ["bam_parts", "cli_message", "content", "normalised", "snapshot"]


def cli_message(stderr: bytes) -> str | None:
    """The CLI's error message: its stderr from the `error: ` line on, without that prefix
    (warnings come before it)."""
    text = stderr.decode("utf-8", "replace")
    at = re.search(r"(?m)^error: ", text)
    return None if at is None else text[at.end() :].rstrip("\n")


def snapshot(d: Path) -> dict[str, bytes]:
    return {str(p.relative_to(d)): p.read_bytes() for p in d.rglob("*") if p.is_file()}


def normalised(text: bytes) -> bytes:
    """Output with the lines that record the command line normalised as the golden harness
    does (`gloci_testkit::differential::observe`): `##glociCommand=` emptied, gloci's `@PG`
    lines dropped. Python names its call in them (`gloci intersect -a a.bed -b b.bed --mode
    pairs`), the CLI the argv as typed."""
    out = []
    for line in text.splitlines(keepends=True):
        if line.startswith(b"@PG\tID:gloci"):
            continue
        if line.startswith(b"##glociCommand="):
            line = b"##glociCommand=<normalised>\n"
        out.append(line)
    return b"".join(out)


def bam_parts(data: bytes) -> tuple[bytes, list[tuple[bytes, int]], bytes]:
    """A BAM file as the golden harness compares it (`bam_parts`): the header text without its
    `@PG` lines, the references, and the record bytes; the compression does not matter."""
    raw = gzip.decompress(data)
    assert raw[:4] == b"BAM\x01", "not BAM"
    at = 4

    def u32() -> int:
        nonlocal at
        v = int.from_bytes(raw[at : at + 4], "little")
        at += 4
        return v

    n = u32()
    text = raw[at : at + n].split(b"\0")[0]
    at += n
    header = b"".join(h for h in text.splitlines(keepends=True) if not h.startswith(b"@PG\t"))
    refs = []
    for _ in range(u32()):
        name = raw[at : at + u32()]
        at += len(name)
        refs.append((name, u32()))
    return header, refs, raw[at:]


def content(data: bytes) -> object:
    """A file's content as compared: BAM as the golden harness reads it, BGZF text
    decompressed, the command-line lines normalised."""
    if data[:2] == b"\x1f\x8b":
        raw = gzip.decompress(data)
        return bam_parts(data) if raw[:4] == b"BAM\x01" else normalised(raw)
    return normalised(data)
