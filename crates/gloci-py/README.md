# gloci for Python

The [gloci](https://github.com/duartemolha/gloci) commands as lazy, streaming Python pipelines.
Steps run in-process, in the same code as the `gloci` command line, so a pipeline gives the same
bytes as the shell pipeline it stands for. Nothing is run as a subprocess, and no temporary files
are passed between steps.

```
pip install gloci                  # no dependencies
pip install gloci[arrow]           # + pyarrow, for to_arrow()
pip install gloci[polars]          # + polars, for to_polars() (no pyarrow needed)
pip install gloci[pandas]          # + pandas and pyarrow, for to_pandas()
```

Every example below runs: the README is a test (`tests/test_readme_examples.py` runs the code
blocks in order, in one scratch directory). The first block writes the small files the others
read.

```python
import pathlib

pathlib.Path("a.bed").write_text(
    "chr1\t100\t200\tgene1\t0\t+\nchr1\t300\t400\tgene2\t0\t-\nchr2\t50\t150\tgene3\t0\t+\n"
)
pathlib.Path("b.bed").write_text("chr1\t150\t350\tpeakA\nchr3\t10\t20\tpeakB\n")
pathlib.Path("calls.vcf").write_text(
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
    "chr1\t100\trs1\tA\tC\t50\tPASS\t.\tGT\t0/1\t1/1\n"
    "chr1\t300\trs2\tG\tT\t10\tq10\t.\tGT\t0/0\t0/1\n"
)
```

```python
import gloci

(gloci.open("a.bed")                  # nothing is read yet
      .intersect(b="b.bed")           # gloci intersect -a a.bed -b b.bed
      .merge(max_gap=50)              # | gloci merge --max-gap 50 -
      .save("hits.bed"))
```

A pipeline is a source and a list of steps, checked as you write them and run only by a terminal
method. `pipeline.cli()` prints the shell pipeline it stands for.

```python
p = gloci.open("a.bed").intersect(b="b.bed").merge(max_gap=50)
print(p.cli())   # gloci intersect -a a.bed -b b.bed | gloci merge --max-gap 50 -
```

The package is typed (`py.typed`, with a stub generated from the command definitions):
`mypy --strict` checks every keyword argument, and `help(gloci.Pipeline.intersect)` lists them
with the command line's help.

## Sources

`gloci.open(src)` starts a pipeline. `format=` overrides format detection. `src` may be:

- a path (`str` or `os.PathLike`; `-` is standard input);
- a file object: binary, or text (encoded as UTF-8), such as `io.BytesIO`, `open(p, "rb")` or
  `gzip.open(p)`;
- an iterable of lines (`str`), or of tuples of values: `("chr1", 15, 25, None)` is the line
  `chr1\t15\t25\t.` (values are written as `str()`, and `None` as `.`);
- an Arrow table or stream: a pyarrow `Table` or `RecordBatchReader`, a polars or pandas
  `DataFrame`. The integer columns `chrom`, `start` and `end` come first, then the others; nulls
  are `.`. It is read through the Arrow C stream interface, so no extra package is needed;
- another pipeline.

```python
import io

rows = [("chr1", 100, 200, "x"), ("chr1", 500, 600, "y")]
print(gloci.open(rows).intersect(b="b.bed").count())                   # a list of tuples
print(gloci.open(io.BytesIO(b"chr1\t1\t5\n")).count())                  # a file object
print(gloci.open("a.bed").intersect(b=gloci.open("b.bed")).count())     # a pipeline
```

```python requires=polars
import polars as pl

df = pl.DataFrame({"chrom": ["chr1"], "start": [100], "end": [200], "name": ["x"]})
print(gloci.open(df).intersect(b="b.bed").count())                      # a data frame
```

`gloci.windows(genome=..., size=...)` and `gloci.random(genome=..., length=..., count=...)`
start a pipeline with windows or random intervals over a genome, with no input:

```python
print(gloci.windows(genome="hg38", size=100_000_000).count())
print(gloci.random(genome="hg38", length=100, count=3, seed=7).count())
```

An input argument such as `b=` takes any source above, or a list of them (`b=["x.bed", rows]`):
`intersect(b=gloci.open("b.bed").slop(genome="hg38", both=10))` runs both pipelines at once.

Python objects are read as the run goes, never whole. A tab in a value, or a null `chrom`,
`start` or `end`, raises `DataError` with `.line` set to the item or row.

An iterator, a file object or a pyarrow `RecordBatchReader` can be read by one run: a second run
raises `UsageError` rather than read nothing. gloci remembers a generator, a file object or a
`RecordBatchReader` that a run has read, whatever pipeline it is given to next. An iterator that
Python cannot hold a weak reference to (`iter()` of a list or tuple, `map()`, `zip()`) is known
only to the pipelines made from the `gloci.open()` (or step) it was given to: pass a list, or a
generator, to be safe. Lists, tuples, tables and data frames are read again at each run.

## Steps

Each step is a method named after its command: `intersect`, `merge`, `vcf_filter` for
`gloci vcf filter`, and so on. Keyword arguments are the command's flags, with `-` written as `_`
(`max_gap=50` is `--max-gap 50`); a flag with no long name keeps its letter (`s=True` on `slop`, `b=` on `intersect`). A one-letter
short flag (`f=0.5`, `g="hg38"`) and a visible alias (`diff_indv_map=`) are keyword arguments
too, typed like the long form and documented as "alias of ...". A name that is not a Python
identifier (`genomecov`'s `--5p` and `--3p`) is not a parameter: pass it as `**{"5p": True}`
(mypy wants a `# type: ignore[arg-type]` there).

- `True` gives a switch, and `False` or `None` leaves it out.
- A list repeats the flag (`b=["x.bed", "y.bed"]`, `region=["chr1:1-100", "chr2"]`). Its items may be
  any values (`column=[4, 5]`); a list that mixes sources and paths needs a type for mypy:
  `several: list[gloci.Source] = ["x.bed", pipeline]`.
- Any other value is passed as `str(value)`.
- A step whose command takes several inputs takes the others as positional sources:
  `a.cat("b.bed", "c.bed")`.
- The output's options belong to the terminal (`save(path, output_format=...)`), the threads
  too (`threads=`), and the log level to the logger (below); a step refuses them with a message
  that says where they go.
- A pipeline is immutable: each step returns a new one, so a pipeline can be kept, extended in
  two ways, and run again.

```python
print(gloci.open("a.bed").intersect(b="b.bed", f=0.5, r=True).count())      # f= is min_frac_a=
print(gloci.windows(g="hg38", w=100_000_000).count())                      # g= is genome=
```

```python
base = gloci.open("a.bed").sort()
print(base.merge().count(), base.slop(genome="hg38", both=100).merge().count())
```

Mistakes are caught when the step is called, as `UsageError` (a `ValueError`), with the keyword
argument in the message, not the flag:

```python
try:
    gloci.open("a.bed").merge(max_gapp=50)
except gloci.UsageError as e:
    print(e)      # ... unexpected keyword argument 'max_gapp'; did you mean 'max_gap'?
```

## Terminals

A terminal method runs the pipeline: every step at once, each on its own thread, chained by
bounded in-memory pipes.

- `save(path, *, output_format=None, write_index=False, compression_level=None, threads=None)`
  writes the file `path`, whose extension gives the format (`.bed.gz` is BGZF-compressed BED).
  The file appears only once the run has succeeded.
- `to_bytes(*, threads=None)` returns the bytes the command line writes to stdout.
- `count(*, threads=None)` returns the number of records: as many as iteration gives. Header
  lines are not records.
- Iteration (`for r in pipeline`) yields `gloci.Record`s, as they come.
  Iteration takes no `threads=`: each step uses the command line's default, every core or
  `$GLOCI_THREADS`; set that variable to limit a loop run from a thread pool.
- `to_arrow()`, `to_polars()`, `to_pandas()` and `to_arrow_stream()` return the whole output as
  a table; see Frames.
- `cli()` is the equivalent shell pipeline.

```python
p = gloci.open("a.bed").intersect(b="b.bed")
p.save("out.bed.gz", write_index=True)
print(p.to_bytes().decode().splitlines()[0])
print(p.count())
for r in p:
    print(r.chrom, r.start, r.end, r.fields[3])
```

Breaking out of the loop (or dropping the iterator) stops the steps, and their temporary files
are removed. Where the iterator may outlive the loop (a loop whose body may raise, the
exception's traceback keeping it), close it: `with iter(pipeline) as records:`, or
`records.close()`.

The module functions run at once, as the command line does, and return the paths written,
sorted:

```python
print(gloci.vcf_stats("calls.vcf", out="stats", freq=True))        # stats.frq
print(gloci.vcf_extract_matrix("calls.vcf", out="matrix"))         # matrix.012, .012.indv, .012.pos
print(gloci.vcf_compare("calls.vcf", "calls.vcf", out="cmp", site_discordance=True))
print(gloci.vcf_popgen("calls.vcf", out="pg", site_pi=True, relatedness2=True))
gloci.open("a.bed").sort().save("a.sorted.bed.gz")
print(gloci.index("a.sorted.bed.gz"))                              # the index next to the file
```

## Frames

Iterating over a pipeline runs it and yields its records as `gloci.Record`s, read back from its
output with gloci's own readers. `count()` and the tables read the output the same way, so they
agree on what a record is.

```python
for r in gloci.open("calls.vcf").vcf_filter(pass_only=True):
    print(r.chrom, r.start, r.end)    # the footprint, 0-based and half-open: chr1 99 100
    print(r.fields[:5])               # the tab-separated fields
    print(r[2], len(r))               # a sequence over its fields; str(r) is the line
```

- `chrom`, `start` and `end` are the record's footprint as gloci reads it: BED's columns, GFF's
  1-based columns 4 and 5 as `start - 1` and `end`, a VCF record's span by the step's
  `vcf_span=` (an insertion is a point: `start == end`), a read's aligned span. BAM output comes
  back as SAM text, BCF output as VCF text, FASTA output as one record per sequence (its fields
  are the name and the sequence).
- They are `None` where the output has no coordinates: a table such as `vcf_extract()`,
  `stats_jaccard()`, `idxstats()`, `depth()`, `coverage(hist=True)` or `view(count=True)`
  writes, or `groupby()` keyed on other columns than 1, 2 and 3; an unmapped read; a line that
  cannot be read as a record.
- Header lines are not records: `#` lines, SAM's `@` lines, `track` and `browser` lines, and a
  table's header row (`vcf_extract()`'s `CHROM  POS  ...`, `groupby(header=True)`'s first
  line), which names its columns.
- Records are handed over in batches of about 4 MiB, so memory stays flat however large the
  output. An error ends the iteration (see Errors); the records just before it may not have been
  given, as a failing step stops its output where it fails.
- A list of records is a source again: `gloci.open(records)` reads their lines.

`to_arrow()` returns a `pyarrow.Table`, `to_polars()` a `polars.DataFrame` and `to_pandas()` a
`pandas.DataFrame`, of the whole output, built in memory:

```python requires=polars
df = gloci.open("a.bed").intersect(b="b.bed", mode="pairs").to_polars()
print(df.columns)    # chrom, start, end, col4, ..., then B's columns
```

```python requires=pyarrow
table = gloci.open("a.bed").to_arrow(names=["chrom", "start", "end", "name", "score", "strand"])
print(table.schema.names)
```

- The columns are `chrom`, `start` and `end` (int64), then VCF's `pos id ref alt qual filter
  info format` and the samples, GFF's `source feature score strand frame attributes`, SAM's
  `qname flag pos mapq cigar rnext pnext tlen seq qual`, or else names from a `#` header line
  (or a table's header row), else `col4`, `col5`, ... by the column's number. A table has no
  `chrom start end` columns of its own.
- Each column's type is inferred: int64 if every value is an integer, else float64 if every
  value is a number, else string. In a number column `.` is null (VCF's `qual` is always
  float64).
- `names=` renames the columns: one name per column, or `UsageError`.
- The table goes to pyarrow and polars through the Arrow PyCapsule interface, without copying:
  `to_arrow_stream()` returns that `gloci.ArrowStream`, which any Arrow consumer takes
  (`pyarrow.table(s)`, `polars.DataFrame(s)`).
- The return types name `pyarrow.Table`, `polars.DataFrame` and `pandas.DataFrame`. mypy knows
  them only where their types are installed (polars ships them; `pandas-stubs` and
  `pyarrow-stubs` do for the others): without those, `to_pandas()` and `to_arrow()` are `Any`
  to mypy.
- pyarrow, polars and pandas are optional: `gloci[arrow]` (pyarrow 23.0.1 or later), `gloci[polars]`
  (polars 1.5 or later alone; no pyarrow needed) or `gloci[pandas]` (pandas 2.2 or later, whose
  DataFrames speak the Arrow stream interface, and pyarrow 23.0.1 or later). Without the one a method needs, it raises `ImportError` naming the extra.

## Errors

| Exception | When |
|---|---|
| `gloci.UsageError` (also a `ValueError`) | a step's keyword arguments are wrong (raised when the step is called) |
| `gloci.DataError` | an input or output cannot be used; `.path` and `.line` say where, when known |
| `gloci.UnsortedError` (a `DataError`) | an input is not sorted; `.records`, `.lines` and `.hint` |
| `gloci.GlociError` | the base class; raised on its own for an internal error (a bug) |

`str(error)` is the command line's message without `error: `.

```python
pathlib.Path("unsorted.bed").write_text("chr1\t5\t6\nchr1\t1\t2\n")
try:
    gloci.open("unsorted.bed").merge().count()
except gloci.UnsortedError as e:
    print(e.path, e.lines)       # unsorted.bed (1, 2)
    print(gloci.open("unsorted.bed").merge(sort=True).count())    # the fix: sort first
try:
    gloci.open("missing.bed").count()
except gloci.DataError as e:
    print(e.path)                # missing.bed
```

An exception a Python source raises (its iterator, or its file object's `read()`) is raised as
it is, the same object, when the step that reads the source reads that far: the run fails there,
as at a bad line of a file. A step that stops early (`head`) may never read that far, and the run
then succeeds. An Arrow stream is different: pyarrow, polars and pandas hand gloci their data
through the Arrow C stream interface, which carries an error as text only. So when a stream that
pyarrow builds from Python code (`pa.RecordBatchReader.from_batches(schema, generator)`) fails,
gloci raises `DataError` with the producer's message (which names the original exception), not
the exception itself.

## Threads and the GIL

- The GIL is released while the steps run, so other Python threads keep going, and several
  pipelines can run in parallel from a thread pool. Python sources are the exception (D91): an
  iterator, a file object or an Arrow stream is pumped with the GIL held, a chunk of at most
  1 MiB at a time, so a slow Python source slows the other Python threads (a path does not).
- `threads=` (on a terminal) is the number of threads per step; it defaults to the command
  line's default (every core, or `$GLOCI_THREADS`). The output is the same at every thread
  count.
- Ctrl-C (`KeyboardInterrupt`) stops the run, removes partial output files and is raised again.
- Warnings and progress go to the `logging` logger `gloci`; its level decides what is made
  (`logging.getLogger("gloci").setLevel(logging.ERROR)` is `-q`).

```python
from concurrent.futures import ThreadPoolExecutor

with ThreadPoolExecutor(2) as pool:
    counts = list(pool.map(lambda f: gloci.open(f).count(threads=1), ["a.bed", "b.bed"]))
print(counts)    # [3, 2]
```

## From pybedtools

The same questions, asked of gloci. pybedtools runs `bedtools` as a subprocess and passes
temporary files between steps; gloci does not. `bedtools intersect` with no mode flag writes the
overlapping part, which gloci's `mode="clip"` does; its default writes each overlapping record of
A once (`-u`). The book's migration page lists every flag.

| pybedtools | gloci |
|---|---|
| `a.intersect(b)` | `a.intersect(b="b.bed", mode="clip").count()` |
| `a.intersect(b, u=True)` | `a.intersect(b="b.bed").count()` |
| `a.intersect(b, v=True)` | `a.intersect(b="b.bed", mode="none").count()` |
| `a.intersect(b, wa=True, wb=True)` | `a.intersect(b="b.bed", mode="pairs").count()` |
| `a.merge(d=100)` | `a.merge(max_gap=100).count()` |
| `a.slop(genome="hg38", b=100)` | `a.slop(genome="hg38", both=100).count()` |
| `a.closest(b, d=True)` | `a.closest(b="b.bed", distance=True).count()` |
| `len(a)` | `a.count()` |
| `a.saveas("out.bed")` | `a.save("out.bed")` |
| `a.to_dataframe()` | `a.to_pandas()` |

Here `a` is `gloci.open("a.bed")`. For a polars data frame, use `a.to_polars()`.

Linux (x86_64, aarch64) and macOS (arm64) are supported. On Windows, use WSL.
