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

- `gloci.open(path)` starts a pipeline from a file. `format=` overrides format detection.
- Each step is a method named after its command: `intersect`, `merge`, `vcf_filter` for
  `gloci vcf filter`, and so on. Keyword arguments are the command's flags, with `-` written as
  `_` (`max_gap=50` is `--max-gap 50`); one-letter flags work too (`f=0.5`, `s=True`).
  - `True` gives a switch, and `False` or `None` leaves it out.
  - A list repeats the flag (`b=["x.bed", "y.bed"]`).
  - Any other value is passed as `str(value)`.
- Steps are checked when they are called, and run only when a terminal method is called:
  - `save(path, *, output_format=None, write_index=False, compression_level=None, threads=None)`;
  - `to_bytes(*, threads=None)`, the bytes the command line writes to stdout;
  - `count(*, threads=None)`, the number of records (header lines are not counted);
  - `cli()`, the equivalent shell pipeline.
- A pipeline is immutable: each step returns a new one, and a pipeline can be run again.
- `threads=` defaults to the command line's default (every core, or `$GLOCI_THREADS`).
  The output is the same at every thread count.

## Errors

| Exception | When |
|---|---|
| `gloci.UsageError` (also a `ValueError`) | a step's keyword arguments are wrong (raised when the step is called) |
| `gloci.DataError` | an input or output cannot be used; `.path` and `.line` say where, when known |
| `gloci.UnsortedError` (a `DataError`) | an input is not sorted; `.records`, `.lines` and `.hint` |
| `gloci.GlociError` | the base class; raised on its own for an internal error (a bug) |

`str(error)` is the command line's message without `error: `.

## Runtime

- The GIL is released while the steps run, so other Python threads keep going.
- Ctrl-C (`KeyboardInterrupt`) stops the run, removes partial output files and is raised again.
- Warnings and progress go to the `logging` logger `gloci`; its level decides what is made.

Linux (x86_64, aarch64) and macOS (arm64) are supported. On Windows, use WSL.
