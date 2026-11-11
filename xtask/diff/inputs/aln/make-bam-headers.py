#!/usr/bin/env python3
"""Rebuilds nosq.bam and notext.bam: BAMs whose header text has no @SQ lines (D119).

samtools cannot write them, so they are assembled here. The references are only in the
binary list; htslib (and gloci) write them out as @SQ lines when they read the file.
Run: python3 make-bam-headers.py
"""
import struct
import zlib


def bgzf(data):
    def block(d):
        c = zlib.compressobj(6, zlib.DEFLATED, -15)
        z = c.compress(d) + c.flush()
        return (b"\x1f\x8b\x08\x04\0\0\0\0\0\xff\x06\0BC\x02\0" + struct.pack("<H", len(z) + 25)
                + z + struct.pack("<II", zlib.crc32(d), len(d)))
    return block(data) + block(b"")


def bam(text, refs, records):
    h = b"BAM\x01" + struct.pack("<i", len(text)) + text + struct.pack("<i", len(refs))
    for name, ln in refs:
        h += struct.pack("<i", len(name) + 1) + name + b"\0" + struct.pack("<i", ln)
    return bgzf(h + b"".join(struct.pack("<i", len(r)) + r for r in records))


def record(ref_id, pos, name, cigar_len):
    seq = b"\x12\x48"  # ACGT
    body = struct.pack("<iiBBHHHiiii", ref_id, pos, len(name) + 1, 60, 4681, 1, 0, 4, -1, -1, 0)
    return body + name + b"\0" + struct.pack("<I", cigar_len << 4) + seq + b"\x1e" * 4


refs = [(b"c1", 100), (b"c2", 50)]
records = [record(0, 10, b"r1", 4), record(1, 5, b"r2", 4)]
open("nosq.bam", "wb").write(bam(b"@HD\tVN:1.6\n", refs, records))
open("notext.bam", "wb").write(bam(b"", refs[:1], records[:1]))
