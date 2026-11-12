# gloci for Python

The [gloci](https://github.com/duartemolha/gloci) commands as lazy, streaming Python pipelines.
Steps run in-process, in the same code as the `gloci` command line, so a pipeline gives the same
bytes as the shell pipeline it stands for. Nothing is run as a subprocess, and no temporary files
are passed between steps.

```python
import gloci

(gloci.open("calls.vcf.gz")          # nothing is read yet
      .intersect(b="targets.bed")     # gloci intersect -a calls.vcf.gz -b targets.bed
      .merge(max_gap=50)              # | gloci merge --max-gap 50 -
      .save("on_target.bed.gz", write_index=True))
```

## Pipelines

- `gloci.open(src)` starts a pipeline. `format=` overrides format detection. `src` may be:
  - a path (`-` is standard input);
  - a file object: binary, or text (encoded as UTF-8), such as `io.BytesIO`, `open(p, "rb")`
    or `gzip.open(p)`;
  - an iterable of lines (`str`), or of tuples of values: `("chr1", 15, 25, None)` is the line
    `chr1\t15\t25\t.` (values are written as `str()`, and `None` as `.`);
  - an Arrow table or stream: a pyarrow `Table` or `RecordBatchReader`, a polars or pandas
    `DataFrame`. The integer columns `chrom`, `start` and `end` come first, then the others;
    nulls are `.`. It is read through the Arrow C stream interface, so no extra package is
    needed;
  - another pipeline.

  Python objects are read as the run goes, never whole. A tab in a value, or a null `chrom`,
  `start` or `end`, raises `DataError` with `.line` set to the item or row.
- An iterator, a file object or a pyarrow `RecordBatchReader` can be read by one run: a second
  run raises `UsageError` rather than read nothing. gloci remembers a generator, a file object
  or a `RecordBatchReader` that a run has read, whatever pipeline it is given to next. An
  iterator that Python cannot hold a weak reference to (`iter()` of a list or tuple, `map()`,
  `zip()`) is known only to the pipelines made from the `gloci.open()` (or step) it was given
  to: pass a list, or a generator, to be safe. Lists, tuples, tables and data frames are read
  again at each run.
- `gloci.windows(g=..., w=...)` and `gloci.random(g=..., l=..., n=...)` start a pipeline with
  windows or random intervals over a genome.
- Each step is a method named after its command: `intersect`, `merge`, `vcf_filter` for
  `gloci vcf filter`, and so on. Keyword arguments are the command's flags, with `-` written as
  `_` (`max_gap=50` is `--max-gap 50`); one-letter flags work too (`f=0.5`, `s=True`).
  - `True` gives a switch, and `False` or `None` leaves it out.
  - A list repeats the flag (`b=["x.bed", "y.bed"]`).
  - An input such as `b=` also takes any source above, or a list of them:
    `intersect(b=gloci.open("b.bed").slop(b=10, g="hg38"))` runs both pipelines at once.
  - Any other value is passed as `str(value)`.
- Steps are checked when they are called, and run only when a terminal method is called:
  - `save(path, *, output_format=None, write_index=False, compression_level=None, threads=None)`;
  - `to_bytes(*, threads=None)`, the bytes the command line writes to stdout;
  - `count(*, threads=None)`, the number of records: as many as iteration gives;
  - iteration (`for r in pipeline`), its records, below;
  - `to_arrow()`, `to_polars()`, `to_pandas()` and `to_arrow_stream()`, its records as a table,
    below;
  - `cli()`, the equivalent shell pipeline.
- A pipeline is immutable: each step returns a new one, and a pipeline can be run again.
- `threads=` defaults to the command line's default (every core, or `$GLOCI_THREADS`).
  The output is the same at every thread count.

## Records and tables

Iterating over a pipeline runs it and yields its records as `gloci.Record`s, read back from its
output with gloci's own readers. `count()` and the tables read the output the same way, so they
agree on what a record is.

```python
for r in gloci.open("calls.vcf.gz").intersect(b="targets.bed"):
    r.chrom, r.start, r.end   # "chr1", 99, 100: the footprint, 0-based and half-open
    r.fields                  # ("chr1", "100", "rs1", "A", "C", ...): the tab-separated fields
    r[3], len(r), str(r)      # a sequence over its fields; str() is the line
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
  table's header row (`vcf_extract()`'s `CHROM  POS  ...`), which names its columns.
- Records are handed over in batches of about 4 MiB, so memory stays flat however large the
  output. Breaking out of the loop (or dropping the iterator) stops the steps, and their
  temporary files are removed. An error is raised once the records before it have been given.
- A list of records is a source again: `gloci.open(records)` reads their lines.

`to_arrow()` returns a `pyarrow.Table`, `to_polars()` a `polars.DataFrame` and `to_pandas()` a
`pandas.DataFrame`, of the whole output, built in memory:

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
- pyarrow, polars and pandas are optional: `pip install gloci[arrow]` (pyarrow 15 or later),
  `gloci[polars]` (polars 1.3 or later alone; no pyarrow needed) or `gloci[pandas]` (pandas and
  pyarrow 15 or later).
  Without the one a method needs, it raises `ImportError` naming the extra.

## Functions

These write files at once, as the command line does, and return the paths written, sorted:

- `gloci.vcf_stats(src, *, out, **kw)`: `gloci vcf stats --out-prefix out` (`freq=True`, ...);
- `gloci.vcf_extract_matrix(src, *, out, **kw)`: the 012 genotype matrix;
- `gloci.vcf_compare(a, b, *, out, **kw)`: genotype concordance;
- `gloci.index(path, **kw)`: the index next to the file (a file only).

## Errors

| Exception | When |
|---|---|
| `gloci.UsageError` (also a `ValueError`) | a step's keyword arguments are wrong (raised when the step is called) |
| `gloci.DataError` | an input or output cannot be used; `.path` and `.line` say where, when known |
| `gloci.UnsortedError` (a `DataError`) | an input is not sorted; `.records`, `.lines` and `.hint` |
| `gloci.GlociError` | the base class; raised on its own for an internal error (a bug) |

`str(error)` is the command line's message without `error: `.

An exception a Python source raises (its iterator, or its file object's `read()`) is raised as
it is, the same object, when the step that reads the source reads that far: the run fails there,
as at a bad line of a file. A step that stops early (`head`) may never read that far, and the run
then succeeds. An Arrow stream is different: pyarrow, polars and pandas hand gloci their
data through the Arrow C stream interface, which carries an error as text only. So when a stream
that pyarrow builds from Python code (`pa.RecordBatchReader.from_batches(schema, generator)`)
fails, gloci raises `DataError` with the producer's message (which names the original exception),
not the exception itself.

## Runtime

- The GIL is released while the steps run, so other Python threads keep going.
- Ctrl-C (`KeyboardInterrupt`) stops the run, removes partial output files and is raised again.
- Warnings and progress go to the `logging` logger `gloci`; its level decides what is made.

Linux (x86_64, aarch64) and macOS (arm64) are supported. On Windows, use WSL.
