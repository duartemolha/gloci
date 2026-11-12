"""Pipelines over files: the same bytes as the CLI (D83, D85, D87)."""

from __future__ import annotations

from pathlib import Path

import pytest

import gloci
from gloci import _native

CALLS_VCF = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=100000>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "chr1\t100\t.\tA\tC\t.\t.\t.\n"
    "chr1\t200\t.\tA\tG\t.\t.\t.\n"
    "chr1\t500\t.\tA\tG\t.\t.\t.\n"
)


@pytest.fixture
def calls(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """`calls.vcf` and the targets `t.bed`, which hold two of its three records; made the cwd."""
    (tmp_path / "calls.vcf").write_text(CALLS_VCF)
    (tmp_path / "t.bed").write_text("chr1\t90\t210\n")
    monkeypatch.chdir(tmp_path)
    return tmp_path


def test_save_equals_cli(ab: Path, cli) -> None:
    done = cli(["intersect", "-a", "a.bed", "-b", "b.bed", "--mode", "pairs", "-o", "o.tsv"], ab)
    assert done.returncode == 0, done.stderr
    want = (ab / "o.tsv").read_bytes()
    (ab / "o.tsv").unlink()

    p = gloci.open("a.bed").intersect(b="b.bed", mode="pairs")
    assert p.save("o.tsv") is None
    assert (ab / "o.tsv").read_bytes() == want
    assert len(want.splitlines()) > 100
    assert p.cli() == "gloci intersect -a a.bed -b b.bed --mode pairs"


def test_chain_equals_shell_pipeline(ab: Path, sh) -> None:
    p = gloci.open("a.bed").intersect(b="b.bed").merge(max_gap=50)
    assert p.cli() == "gloci intersect -a a.bed -b b.bed | gloci merge --max-gap 50 -"
    shell = sh(p.cli(), ab)
    assert shell.returncode == 0, shell.stderr
    got = p.to_bytes()
    assert got == shell.stdout
    # Both steps did work: the merge joins some of the intersected records.
    assert 100 < got.count(b"\n") < 740


def test_vcf_header_records_command(calls: Path) -> None:
    gloci.open("calls.vcf").intersect(b="t.bed").save("out.vcf")
    lines = (calls / "out.vcf").read_text().splitlines()
    want = f"##glociCommand=gloci intersect -a calls.vcf -b t.bed -o out.vcf; Version={gloci.__version__}"
    assert want in lines
    assert [l for l in lines if not l.startswith("#")] == CALLS_VCF.splitlines()[3:5]
    # Without a file there is no -o, as on the command line; -t never shows.
    text = gloci.open("calls.vcf").intersect(b="t.bed").to_bytes(threads=2).decode()
    want = f"##glociCommand=gloci intersect -a calls.vcf -b t.bed; Version={gloci.__version__}"
    assert want in text.splitlines()


def test_vcf_header_of_a_chained_step_shows_the_stream_input(calls: Path, sh) -> None:
    p = gloci.open("calls.vcf").intersect(b="t.bed").view()
    assert p.cli() == "gloci intersect -a calls.vcf -b t.bed | gloci view -"
    got = p.to_bytes()
    assert f"##glociCommand=gloci view -; Version={gloci.__version__}" in got.decode().splitlines()
    shell = sh(p.cli(), calls)
    assert shell.returncode == 0, shell.stderr
    assert got == shell.stdout


def test_lazy_and_immutable(ab: Path) -> None:
    # Nothing runs until a terminal method: a missing file is not an error yet.
    p = gloci.open("missing.bed")
    q = p.merge(max_gap=5)
    r = p.intersect(b="b.bed")
    assert isinstance(q, gloci.Pipeline) and q is not p
    assert p.cli() == "gloci view missing.bed"
    assert q.cli() == "gloci merge --max-gap 5 missing.bed"
    assert r.cli() == "gloci intersect -a missing.bed -b b.bed"
    # A step returns a new pipeline and leaves its own alone.
    s = q.merge()
    assert q.cli() == "gloci merge --max-gap 5 missing.bed"
    assert s.cli() == "gloci merge --max-gap 5 missing.bed | gloci merge -"
    with pytest.raises(gloci.DataError):
        q.to_bytes()
    for target in (p, q):
        with pytest.raises(AttributeError):
            target.extra = 1  # type: ignore[attr-defined]
        with pytest.raises(AttributeError):
            target._steps = ()  # type: ignore[misc]
    # A pipeline runs again from its sources each time.
    m = gloci.open("a.bed").merge(max_gap=10)
    assert m.to_bytes() == m.to_bytes()
    assert m.to_bytes() == gloci.open("a.bed").merge(max_gap=10).to_bytes()


def test_count(ab: Path, calls: Path) -> None:
    p = gloci.open(ab / "a.bed").intersect(b=ab / "b.bed")
    text = p.to_bytes()
    assert p.count() == text.count(b"\n") > 0
    # Header lines are not records.
    assert gloci.open("calls.vcf").intersect(b="t.bed").count() == 2
    assert gloci.open("calls.vcf").count() == 3
    # FASTA output: a record is a `>` line and its sequence lines.
    (calls / "ref.fa").write_text(">chr1\nACGTACGTAC\n>chr2\nGGGGCCCCAA\n")
    (calls / "ref.fa.fai").write_text("chr1\t10\t6\t10\t11\nchr2\t10\t23\t10\t11\n")
    (calls / "r.bed").write_text("chr1\t0\t5\nchr2\t2\t6\n")
    seqs = gloci.open("r.bed").getfasta(fasta="ref.fa")
    assert seqs.to_bytes() == b">chr1:0-5\nACGTA\n>chr2:2-6\nGGCC\n"
    assert seqs.count() == 2
    region = gloci.open("ref.fa").faidx(region="chr1:1-5")
    assert region.to_bytes() == b">chr1:1-5\nACGTA\n"
    assert region.count() == 1


@pytest.mark.xfail(strict=True, reason="Task 8: count via re-read")
def test_count_summary_tables(ab: Path) -> None:
    # A summary table's header row has no `#`, so count() takes it for a record until Task 8
    # counts what iteration reads.
    p = gloci.open("a.bed").stats_jaccard(b="b.bed")
    assert p.to_bytes().startswith(b"intersection\tunion\tjaccard\tn_intersections\n")
    assert p.count() == 1


def test_threads_deterministic(pair, tmp_path: Path) -> None:
    a, b = pair(20_000, seed=7)
    p = gloci.open(a).intersect(b=b).merge(max_gap=10)
    one = p.to_bytes(threads=1)
    assert one.count(b"\n") > 1000
    for t in (2, 8):
        assert p.to_bytes(threads=t) == one
        out = tmp_path / f"t{t}.bed"
        p.save(out, threads=t)
        assert out.read_bytes() == one
    assert p.count(threads=3) == one.count(b"\n")
    assert "-t" not in p.cli().split()


def test_paths_are_kept_as_given(ab: Path) -> None:
    p = gloci.open(ab / "a.bed").intersect(b=[Path("b.bed"), "b.bed"]).merge("a.bed")
    assert p.cli() == (
        f"gloci intersect -a {ab / 'a.bed'} -b b.bed -b b.bed | gloci merge - a.bed"
    )
    assert p.to_bytes().count(b"\n") > 0
    # A path that needs quoting is quoted as a shell would need it.
    (ab / "my a.bed").write_bytes((ab / "a.bed").read_bytes())
    assert gloci.open("my a.bed").merge().cli() == "gloci merge 'my a.bed'"


def test_save_options(ab: Path) -> None:
    gloci.open("a.bed").merge().save("m.bed.gz", write_index=True, compression_level=1)
    assert (ab / "m.bed.gz").read_bytes()[:2] == b"\x1f\x8b"
    assert (ab / "m.bed.gz.tbi").exists()
    plain = gloci.open("a.bed").merge().to_bytes()
    assert gloci.open("m.bed.gz").to_bytes() == plain
    gloci.open("a.bed").merge().save("m.out", output_format="bed")
    assert (ab / "m.out").read_bytes() == plain
    # The command line recorded in a header spells the option as save()'s keyword does.
    (ab / "c.vcf").write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    gloci.open("c.vcf").view().save("c.out", output_format="vcf")
    want = f"##glociCommand=gloci view -o c.out --output-format vcf c.vcf; Version={gloci.__version__}"
    assert want in (ab / "c.out").read_text().splitlines()


def test_step_methods_cover_bindings() -> None:
    bindings = _native.bindings()
    assert len(bindings) > 40
    for b in bindings:
        method = getattr(gloci.Pipeline, b["method"])
        assert callable(method) and method.__name__ == b["method"]
    assert {"intersect", "merge", "vcf_filter", "pair_intersect", "stats_jaccard"} <= {
        b["method"] for b in bindings
    }
    specs = {tuple(s["path"]) for s in _native.command_specs()}
    assert all(tuple(b["path"]) in specs for b in bindings)


def test_kwargs_to_argv() -> None:
    kw = {"b": Path("x.bed"), "f": 0.5, "s": True, "split": False, "names": None}
    assert _native.kwargs_to_argv(["intersect"], kw) == ["-b", "x.bed", "-f", "0.5", "-s"]
    assert _native.kwargs_to_argv(["merge"], {"max_gap": -2}) == ["--max-gap=-2"]
