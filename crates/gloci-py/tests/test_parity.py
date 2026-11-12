"""Golden parity (P2b Task 10, spec 8.5): every committed golden case of the CLI
(`crates/gloci/tests/golden/<case>/`), replayed through Python, gives the CLI's output.

Each case's argv becomes one Python call (`parity.translate`), run in a scratch copy of the
case's directory:
- with `GLOCI_BIN` set (always on CI), the output equals what that binary writes for the same
  argv, byte for byte, but for the lines that record the command line, which name the call as
  each side spells it (`##glociCommand=`, gloci's `@PG` line): BAM and CRAM output is compared
  as the golden harness compares it, decoded with those lines left out;
- without it, the output equals the committed `expected_gloci` (gloci's own output where it
  differs from the reference tool), or `expected` for an exact case with no difference; other
  cases skip.
A case that exits 1 raises `DataError`, whose `str` is the CLI's message, and one that exits 2
raises `UsageError`.
"""

from __future__ import annotations

import concurrent.futures
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

import pytest

import gloci
from parity.outputs import bam_parts, cli_message, content, normalised, snapshot
from parity.translate import Call, UntranslatableError, split_args, translate

TESTS = Path(__file__).resolve().parent
REPO = TESTS.parents[2]
GOLDEN = REPO / "crates" / "gloci" / "tests" / "golden"
PYTHON_YML = REPO / ".github" / "workflows" / "python.yml"
#: The share of golden cases that must become a Python call (R1).
COVERAGE = 0.95


@dataclass(frozen=True)
class GoldenCase:
    name: str
    dir: Path
    args: str
    compare: str
    expected_difference: str | None
    expected: bytes
    expected_gloci: bytes | None

    @staticmethod
    def load(name: str) -> GoldenCase:
        d = GOLDEN / name

        def line(f: str) -> str | None:
            p = d / f
            return p.read_text().strip() if p.exists() else None

        gloci_out = d / "expected_gloci"
        return GoldenCase(
            name=name,
            dir=d,
            args=line("args") or "",
            compare=line("compare") or "exact",
            expected_difference=line("expected_difference"),
            expected=(d / "expected").read_bytes(),
            expected_gloci=gloci_out.read_bytes() if gloci_out.exists() else None,
        )


CASES = sorted(p.name for p in GOLDEN.iterdir() if p.is_dir())


# What a run gave: the CLI's, or Python's.


@dataclass
class Outcome:
    #: The exit status, or what it would be: 0, 1 for DataError, 2 for UsageError.
    code: int
    #: The output: stdout, or the `-o` file.
    out: bytes = b""
    #: The error message, without `error: ` (the CLI's stderr from its `error:` line).
    message: str | None = None
    #: The files the run wrote (eager module functions), by name.
    files: dict[str, bytes] = field(default_factory=dict)


class CliRuns:
    """`GLOCI_BIN`'s outcome of each case, run once per session in a scratch copy of the case;
    the first `get()` runs them all on a thread pool (R4)."""

    def __init__(self, binary: str, scratch: Path) -> None:
        self.binary = binary
        self.scratch = scratch
        self._runs: dict[str, Outcome] | None = None

    def _one(self, name: str) -> Outcome:
        case = GoldenCase.load(name)
        work = self.scratch / name
        shutil.copytree(case.dir, work)
        before = snapshot(work)
        done = subprocess.run(
            ["gloci", *split_args(case.args)],
            executable=self.binary,
            cwd=work,
            stdin=subprocess.DEVNULL,
            capture_output=True,
        )
        after = snapshot(work)
        files = {n: b for n, b in after.items() if before.get(n) != b}
        try:
            output = translate(case.args).output
        except UntranslatableError:
            output = None
        out = files.get(output, b"") if output is not None else done.stdout
        return Outcome(done.returncode, out, cli_message(done.stderr), files)

    def get(self, name: str) -> Outcome:
        if self._runs is None:
            workers = min(8, os.cpu_count() or 1)
            with concurrent.futures.ThreadPoolExecutor(workers) as pool:
                self._runs = dict(zip(CASES, pool.map(self._one, CASES), strict=True))
        return self._runs[name]

    def sam(self, data: bytes) -> bytes:
        """CRAM (or BAM) as SAM text, by `GLOCI_BIN view -O sam`."""
        done = subprocess.run(
            ["gloci", "view", "-O", "sam", "-"],
            executable=self.binary,
            input=data,
            capture_output=True,
        )
        assert done.returncode == 0, done.stderr.decode()
        return done.stdout


def gloci_bin() -> str | None:
    """$GLOCI_BIN; required when CI=true, where the cases that only the CLI can check must not
    skip."""
    path = os.environ.get("GLOCI_BIN")
    if path:
        if not os.access(path, os.X_OK):
            pytest.fail(f"GLOCI_BIN={path} is not an executable file")
        return path
    if os.environ.get("CI") == "true":
        pytest.fail(
            "GLOCI_BIN must name the gloci binary when CI=true, or the golden cases that only "
            "the CLI can check skip: cargo build -p gloci, then export "
            "GLOCI_BIN=$PWD/target/debug/gloci"
        )
    return None


@pytest.fixture(scope="session")
def cli_runs(tmp_path_factory: pytest.TempPathFactory) -> CliRuns | None:
    binary = gloci_bin()
    return None if binary is None else CliRuns(binary, tmp_path_factory.mktemp("cli-runs"))


def run_python(call: Call) -> Outcome:
    try:
        out = call()
    except gloci.DataError as e:
        return Outcome(1, message=str(e))
    except gloci.UsageError as e:
        return Outcome(2, message=str(e))
    return Outcome(0, out, files={str(p): p.read_bytes() for p in call.written})


def same_output(compare: str, got: bytes, want: bytes, cli: CliRuns | None) -> bool:
    if compare == "bam":
        return bam_parts(got) == bam_parts(want)
    if compare == "cram":
        assert cli is not None
        return normalised(cli.sam(got)) == normalised(cli.sam(want))
    return normalised(got) == normalised(want)


def committed(case: GoldenCase) -> Outcome | None:
    """What the committed golden files say the outcome is, where they can say it: the exit
    status of an exit-status case; gloci's own output (`expected_gloci`); the reference tool's
    output where gloci matches it byte for byte. None otherwise, and for CRAM, which only the
    CLI decodes here."""
    if case.compare == "exit-status":
        return Outcome(int((case.expected_gloci or case.expected).decode()))
    if case.expected_gloci is not None and case.compare != "cram":
        return Outcome(0, case.expected_gloci)
    if case.compare == "exact" and case.expected_difference is None:
        return Outcome(0, case.expected)
    return None


@pytest.mark.parametrize("name", CASES)
def test_golden_parity(
    name: str, cli_runs: CliRuns | None, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    case = GoldenCase.load(name)
    try:
        call = translate(case.args)
    except UntranslatableError as e:
        pytest.skip(f"untranslatable ({e}); test_parity_coverage checks the reason")
    if cli_runs is not None:
        want = cli_runs.get(name)
        # the CLI's output is the golden's, as `cargo test` checks
        assert case.compare != "exit-status" or want.code == committed(case).code
    else:
        want = committed(case)
        if want is None:
            pytest.skip(f"compare mode {case.compare}: set GLOCI_BIN to compare with the CLI")
    work = tmp_path / "case"
    shutil.copytree(case.dir, work)
    monkeypatch.chdir(work)
    got = run_python(call)
    where = f"{call.python()}\n  (gloci {case.args})"
    assert got.code == want.code, (
        f"{where}: exit {want.code} on the command line, {got.code} in Python"
        f" ({got.message or want.message})"
    )
    if got.code == 1 and cli_runs is not None:
        assert got.message == want.message, where
    if got.code != 0 or (case.compare == "exit-status" and cli_runs is None):
        return
    assert same_output(case.compare, got.out, want.out, cli_runs), where
    if call.kind == "function" and cli_runs is not None:
        assert got.files == want.files, where


def test_parity_coverage() -> None:
    """At least 95% of the golden cases are one Python call each, and every case that is not
    is a command that stays on the command line (D87's CLI_ONLY)."""
    untranslated: dict[str, UntranslatableError] = {}
    for name in CASES:
        try:
            translate(GoldenCase.load(name).args)
        except UntranslatableError as e:
            untranslated[name] = e
    share = 1 - len(untranslated) / len(CASES)
    assert share >= COVERAGE, f"{share:.1%} of {len(CASES)} golden cases translate"
    not_cli_only = {n: str(e) for n, e in untranslated.items() if not e.cli_only}
    assert not not_cli_only, not_cli_only


def test_ci_requires_gloci_bin(tmp_path: Path) -> None:
    """CI compares with the CLI: python.yml exports GLOCI_BIN before every job's pytest run
    (the parity run included), and with CI=true and no GLOCI_BIN the parity tests fail rather
    than skip the cases only the CLI can check."""
    yml = PYTHON_YML.read_text()
    jobs = re.split(r"(?m)^  (?=[\w-]+:\n)", yml.split("\njobs:\n", 1)[1])
    pytest_jobs = [j for j in jobs if "pytest crates/gloci-py/tests" in j]
    assert len(pytest_jobs) >= 2, "python.yml: the pytest and polars-only jobs"
    for job in pytest_jobs:
        name = job.split(":", 1)[0]
        export = job.find('echo "GLOCI_BIN=$PWD/target/debug/gloci" >> "$GITHUB_ENV"')
        run = job.find("pytest crates/gloci-py/tests")
        assert 0 <= export < run, f"python.yml: job {name} must export GLOCI_BIN before pytest"
        assert job.find("cargo build -p gloci\n") < export, f"python.yml: job {name}"
        # the run takes the parity tests
        command = job[run:].split("\n", 1)[0]
        assert "parity" not in command and "--ignore" not in command and " -k " not in command

    env = {k: v for k, v in os.environ.items() if k != "GLOCI_BIN"}
    env["CI"] = "true"
    one = f"{TESTS / 'test_parity.py'}::test_golden_parity[{CASES[0]}]"
    done = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", one],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
    )
    assert done.returncode != 0, done.stdout
    assert "GLOCI_BIN must name the gloci binary when CI=true" in done.stdout, done.stdout


# The translation rules (R1), beyond what the golden cases use


@pytest.mark.parametrize(
    ("argv", "python"),
    [
        # canonical keyword names; a flag with several values; a negative value; a cluster
        (
            "intersect -a a.bed -b x.bed y.bed -f 0.5 -s",
            "gloci.open('a.bed').intersect(b=['x.bed', 'y.bed'], min_frac_a='0.5', "
            "same_strand=True).to_bytes()",
        ),
        (
            "cluster -sd -5 a.bed",
            "gloci.open('a.bed').cluster(strand=True, max_gap='-5').to_bytes()",
        ),
        # standard input is the case's file, as a file object
        ("merge -", "gloci.open(open('in.bed', 'rb')).merge().to_bytes()"),
        # -o, and the output options, go to save(); -t to the terminal
        (
            "merge -q a.bed -o out.bed.gz --write-index --compression-level 3 -t 2",
            "gloci.open('a.bed').merge().save('out.bed.gz', write_index=True, "
            "compression_level=3, threads=2)",
        ),
        (
            "view x.sam -O bam",
            "gloci.open('x.sam').view().save('parity-output', output_format='bam')",
        ),
        # n-ary steps, a positional before the input, an attached value
        (
            "sort -m a.bam b.bam",
            "gloci.open('a.bam').sort('b.bam', merge_sorted=True).to_bytes()",
        ),
        ("reheader h.sam in.bam", "gloci.open('in.bam').reheader(header='h.sam').to_bytes()"),
        (
            "vcf norm --split=snps x.vcf --fasta r.fa",
            "gloci.open('x.vcf').vcf_norm(split='snps', fasta='r.fa').to_bytes()",
        ),
        (
            "vcf merge --file-list list.txt",
            "gloci.open(*<the lines of list.txt>).vcf_merge().to_bytes()",
        ),
        # module sources and functions
        ("windows -g g.genome -w 10", "gloci.windows(genome='g.genome', size='10').to_bytes()"),
        ("windows -i a.bed -w 10", "gloci.open('a.bed').windows(size='10').to_bytes()"),
        (
            "random -g g.genome -l 10 -n 5 --seed 1",
            "gloci.random(genome='g.genome', length='10', count='5', seed='1').to_bytes()",
        ),
        (
            "vcf extract x.vcf --genotype-matrix 012 --out-prefix o",
            "gloci.vcf_extract_matrix('x.vcf', genotype_matrix='012', out='o')",
        ),
        (
            "vcf extract x.vcf --format GT",
            "gloci.open('x.vcf').vcf_extract(format='GT').to_bytes()",
        ),
        ("vcf compare a.vcf b.vcf --out-prefix o", "gloci.vcf_compare('a.vcf', 'b.vcf', out='o')"),
        ("index a.bed.gz", "gloci.index('a.bed.gz')"),
    ],
)
def test_translate(argv: str, python: str) -> None:
    assert translate(argv, stdin="in.bed").python() == python


@pytest.mark.parametrize(
    ("argv", "cli_only"),
    [
        ("flags 4", True),
        ("check --quick x.bam", True),
        ("genomes", True),
        # an argument the step never takes, which no module function offers either
        ("sort --check a.bed", True),
        ("merge --no-such-flag a.bed", False),
        ("intersect -b b.bed", False),  # no input for the pipeline
        ("merge -", False),  # standard input, with no file for it
    ],
)
def test_untranslatable(argv: str, cli_only: bool) -> None:
    with pytest.raises(UntranslatableError) as e:
        translate(argv)
    assert e.value.cli_only == cli_only


SAM = (
    "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"
    "r1\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n"
    "r2\t16\tchr1\t20\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n"
)
VCF = (
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    "chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/1\nchr1\t200\t.\tA\tG\t10\tPASS\t.\tGT\t1/1\n"
)


@pytest.mark.parametrize(
    ("argv", "stdin"),
    [
        ("merge -", "a.bed"),
        ("intersect -a - -b b.bed --mode count", "a.bed"),
        ("intersect -a a.bed -b b.bed --mode pairs -o out.tsv", None),
        ("sort a.bed -o out.bed.gz --write-index -t 2", None),
        ("view calls.vcf -o out.vcf.gz --compression-level 1", None),
        ("view x.sam -O bam", None),
        ("vcf stats calls.vcf --out-prefix st --freq", None),
        ("index sorted.bed.gz", None),
    ],
)
def test_translated_calls_equal_the_cli(
    argv: str, stdin: str | None, ab: Path, gloci_bin: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Standard input, -o, output options, threads and module functions: the same output
    files and bytes as the CLI."""
    (ab / "x.sam").write_text(SAM)
    (ab / "calls.vcf").write_text(VCF)
    gloci.open("a.bed").save(ab / "sorted.bed.gz")
    inputs = ("a.bed", "b.bed", "x.sam", "calls.vcf", "sorted.bed.gz")
    there, here = ab / "cli", ab / "py"
    for d in (there, here):
        d.mkdir()
        for f in inputs:
            shutil.copy(ab / f, d / f)

    before = snapshot(there)
    with open(there / stdin if stdin else os.devnull, "rb") as f:
        done = subprocess.run(
            ["gloci", *split_args(argv)],
            executable=gloci_bin,
            cwd=there,
            stdin=f,
            capture_output=True,
        )
    assert done.returncode == 0, done.stderr.decode()
    wrote = {n: b for n, b in snapshot(there).items() if before.get(n) != b}

    monkeypatch.chdir(here)
    before = snapshot(here)
    call = translate(argv, stdin=here / stdin if stdin else None)
    got = call()
    written = {n: b for n, b in snapshot(here).items() if before.get(n) != b}
    written.pop("parity-output", None)

    want = wrote.get(call.output, b"") if call.output else done.stdout
    assert content(got) == content(want), call.python()
    assert sorted(written) == sorted(wrote), call.python()
    for name, data in wrote.items():
        if name.endswith((".tbi", ".csi")) and "vcf" in name:
            continue  # its offsets follow the header, which records the command line
        assert content(written[name]) == content(data), name
    if call.kind == "function":
        assert sorted(str(p) for p in call.written) == sorted(wrote)
