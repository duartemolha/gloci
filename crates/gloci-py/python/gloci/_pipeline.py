"""Lazy, immutable pipelines over the gloci commands (D85-D88).

A pipeline is a source and a list of steps. Each step records a command and its keyword
arguments, checked when the step is called; nothing runs until a terminal method. A terminal
runs every step at once, each on its own thread, chained by bounded in-memory pipes, in the
same code as the command line: the output is the bytes the shell pipeline `cli()` prints would
write.

A source (D88) is a path, a file object, an iterable of lines or tuples, an Arrow table or
stream, or another pipeline; `classify_source()` tells which. Python objects are read as the
run goes, never whole.
"""

from __future__ import annotations

import os
from collections.abc import Iterable, Mapping
from typing import Any

from . import _native
from .errors import UsageError

__all__ = ["Pipeline", "classify_source", "open"]

#: Each step method's binding, by method name (`_native.bindings()`).
_BINDINGS: dict[str, dict[str, Any]] = {b["method"]: b for b in _native.bindings()}

# Plan tokens (`_native.run_plan`): the pipe stage `i` writes, and the terminal's output.
_OUT = "\0out"

#: What a source may be, for messages.
_ACCEPTED = (
    "a path (str or os.PathLike), a file object, an iterable of lines (str) or tuples of "
    "values, an Arrow table or stream (an object with __arrow_c_stream__, such as a pyarrow "
    "Table or a polars or pandas DataFrame) or a gloci Pipeline"
)


def _pipe(i: int) -> str:
    return f"\0pipe:{i}"


def _placeholder(text: str) -> str:
    """`text` marked as a placeholder for `_native.command_line`, which shows it as it is. The
    marks hold a NUL, which no value a user gives can, so a path that only looks like a
    placeholder is quoted like any other word."""
    return f"\0{{{text}\0}}"


# Sources (D88)


class Source:
    """What a pipeline reads: what `classify_source()` makes of an object.

    Attributes:
        obj: The object (for a path, the path as a `str`).
        display: Its name in messages: the path, `f.name`, `<list>`, `<pyarrow.Table>`.
        once: It can be read only once (an iterator, a file object).
    """

    __slots__ = ("display", "obj", "once", "used")
    #: The kind the plan gives `_native.run_plan` for a Python object.
    kind = ""

    def __init__(self, obj: Any, display: str, once: bool = False) -> None:
        self.obj = obj
        self.display = display
        self.once = once
        self.used = False

    def __repr__(self) -> str:
        return f"<{type(self).__name__} source {self.display}>"


class Path(Source):
    """A file, by its path (`-` is standard input)."""

    __slots__ = ()
    kind = "path"


class FileObject(Source):
    """An object with `read()`: it gives bytes, or text, which is encoded as UTF-8."""

    __slots__ = ()
    kind = "file"


class Lines(Source):
    """An iterable of lines (`str` or `bytes`, one per item, with or without its line break),
    or of tuples of values, written tab-joined as their `str()`, with None as `.`."""

    __slots__ = ()
    kind = "lines"


class Frame(Source):
    """An object with `__arrow_c_stream__` (a pyarrow Table or RecordBatchReader, a polars or
    pandas DataFrame): each row a line, the integer columns `chrom start end` first."""

    __slots__ = ()
    kind = "frame"


class Nested(Source):
    """A `Pipeline`, run as its own chain of steps, whose output the reading step takes."""

    __slots__ = ()
    kind = "nested"


class _Command(Source):
    """A module source (`gloci.windows(g=...)`, `gloci.random(...)`): its command is the
    pipeline's first stage. `obj` is its `_native.module_functions()` entry and its keyword
    arguments."""

    __slots__ = ()
    kind = "command"


def _file_name(f: object) -> str:
    """A file object's display name: `f.name` when it is a path, else `<BytesIO>`."""
    name = getattr(f, "name", None)
    if isinstance(name, os.PathLike):
        name = os.fspath(name)
    if isinstance(name, bytes):
        name = os.fsdecode(name)
    if isinstance(name, str) and name and "\0" not in name:
        return name
    return f"<{type(f).__name__}>"


def _frame_name(obj: object) -> str:
    """An Arrow object's display name: `<pyarrow.Table>`, `<polars.DataFrame>`."""
    t = type(obj)
    return f"<{t.__module__.split('.')[0]}.{t.__qualname__}>"


def classify_source(obj: object, what: str = "gloci.open()") -> Source:
    """What `obj` is as a source (D88): a `Path`, `FileObject`, `Lines`, `Frame` or `Nested`.

    Raises:
        TypeError: `obj` is none of them; the message lists what a source may be.
        UsageError: A path holds a NUL character.
    """
    if isinstance(obj, Pipeline):
        return Nested(obj, "<pipeline>")
    if isinstance(obj, (str, os.PathLike)):
        path = obj if isinstance(obj, str) else os.fsdecode(os.fspath(obj))
        if "\0" in path:
            raise UsageError(f"{what}: the path {path!r} holds a NUL character")
        return Path(path, path)
    if isinstance(obj, (bytes, bytearray, memoryview)):
        raise TypeError(
            f"{what} takes {_ACCEPTED}, not {type(obj).__name__}; for data in memory, pass "
            "io.BytesIO(data)"
        )
    if hasattr(type(obj), "__arrow_c_stream__"):
        return Frame(obj, _frame_name(obj))
    if callable(getattr(obj, "read", None)):
        return FileObject(obj, _file_name(obj), once=True)
    if isinstance(obj, Iterable) and not isinstance(obj, Mapping):
        return Lines(obj, f"<{type(obj).__name__}>", once=iter(obj) is obj)
    raise TypeError(f"{what} takes {_ACCEPTED}, not {type(obj).__name__}")


def _is_source(value: object) -> bool:
    """Whether a keyword argument's value (or an item of a list of them) is a source rather
    than a value (D86: anything else is its `str()`)."""
    if value is None or isinstance(value, (str, bytes, bytearray, os.PathLike, int, float)):
        return False
    if isinstance(value, Pipeline) or hasattr(type(value), "__arrow_c_stream__"):
        return True
    if callable(getattr(value, "read", None)):
        return True
    return isinstance(value, Iterable) and not isinstance(value, Mapping)


def _kwarg(where: str, name: str, value: Any) -> Any:
    """A keyword argument's value, with the sources in it made `Source`s (D87: `b=` takes any
    source, or a list of them). A list or tuple whose items are all tuples (or lists) is one
    source, of rows; any other list repeats the flag, each item a value or a source."""
    what = f"{where}: {name}="

    def one(v: Any) -> Any:
        if _is_source(v):
            return classify_source(v, what)
        if isinstance(v, (str, os.PathLike)) and "\0" in os.fsdecode(os.fspath(v)):
            raise UsageError(f"{name}= holds a NUL character, which no argument can hold")
        return v

    if isinstance(value, (list, tuple)):
        if value and all(isinstance(v, (tuple, list)) for v in value):
            return classify_source(value, what)
        items = [one(v) for v in value]
        return items if any(isinstance(i, Source) for i in items) else value
    return one(value)


def _placeholders(kwargs: dict[str, Any]) -> dict[str, Any]:
    """`kwargs` with each source as a placeholder value, to check them before anything runs."""

    def one(v: Any) -> Any:
        return "<python object>" if isinstance(v, Source) else v

    return {k: [one(i) for i in v] if isinstance(v, list) else one(v) for k, v in kwargs.items()}


class _Step:
    """One step: the binding of its method, its keyword arguments and its extra sources."""

    __slots__ = ("binding", "extra", "kwargs")

    def __init__(
        self, binding: dict[str, Any], kwargs: dict[str, Any], extra: tuple[Source, ...]
    ) -> None:
        self.binding = binding
        self.kwargs = kwargs
        self.extra = extra


class _Plan:
    """A run's stages and Python-object sources (`_native.run_plan`), built from a pipeline and
    the pipelines nested in it, each nested chain's stages before the stage that reads it.

    Each stage has its run argv (`argv`, when `run`) and its command line as the shell shows it
    (`command_line`): a Python object is `<python object N>`, numbered in the order the run
    reads them, and a nested pipeline is bash's process substitution, `<(gloci ...)`.
    """

    def __init__(self, threads: object = None, run: bool = True) -> None:
        self.threads = threads
        self.run = run
        self.stages: list[dict[str, Any]] = []
        self.sources: dict[str, dict[str, Any]] = {}
        self.lines: list[str] = []
        self._objects = 0
        self._once: list[Source] = []

    def input(self, src: Source) -> tuple[str, str]:
        """A source read by a stage: the word its argv takes, and the word shown."""
        if isinstance(src, Path):
            return src.obj, src.obj
        if isinstance(src, Nested):
            j, lines = self.chain(src.obj, None, [])
            return _pipe(j), _placeholder(f"<({' | '.join(lines)})")
        self._objects += 1
        shown = _placeholder(f"<python object {self._objects}>")
        if not self.run:
            return shown, shown
        if src.once:
            if any(s.obj is src.obj for s in self._once):
                raise UsageError(
                    f"{src.display} is used twice in one run; an iterator or file object can be "
                    "read only once, so pass a list, or a path, for each use"
                )
            if src.used:
                raise UsageError(
                    f"{src.display} was read by an earlier run, and an iterator or file object "
                    "can be read only once; pass a list (or a path) to run a pipeline again"
                )
            self._once.append(src)
        token = f"\0src:{len(self.sources)}"
        self.sources[token] = {"kind": src.kind, "obj": src.obj, "display": src.display}
        return token, shown

    def kwargs(self, kwargs: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
        """Keyword arguments with their sources as inputs: for the argv, and as shown."""
        real: dict[str, Any] = {}
        shown: dict[str, Any] = {}
        for name, v in kwargs.items():
            if isinstance(v, Source):
                real[name], shown[name] = self.input(v)
            elif isinstance(v, list) and any(isinstance(i, Source) for i in v):
                pairs = [self.input(i) if isinstance(i, Source) else (i, i) for i in v]
                real[name] = [r for r, _ in pairs]
                shown[name] = [s for _, s in pairs]
            else:
                real[name] = shown[name] = v
        return real, shown

    def add(self, argv: list[str] | None, shown: list[str]) -> int:
        """Add a stage; its index."""
        stage: dict[str, Any] = {"command_line": _native.command_line(shown)}
        if argv is not None:
            stage["argv"] = argv
        self.stages.append(stage)
        return len(self.stages) - 1

    def _to(
        self, b: dict[str, Any], idx: int, last: bool, out: list[str] | None, shown: list[str]
    ) -> tuple[list[str], list[str]]:
        """The terminal words of stage `idx`, for its argv and as shown: `out` and `shown` for
        the last stage of the chain, else its pipe; `-t` for every stage that has it."""
        if last and out is not None:
            real, seen = list(out), list(shown)
        else:
            real, seen = ["-o", _pipe(idx)], []
        if self.threads is not None and b["threads"]:
            real += ["-t", str(self.threads)]
        return real, seen

    def chain(
        self, p: Pipeline, out: list[str] | None, shown_out: list[str]
    ) -> tuple[int, list[str]]:
        """Add `p`'s stages. The last writes `out` (its pipe when None) and shows `shown_out`.
        Returns its index and the command lines of `p`'s own chain, as the shell shows them."""
        lines: list[str] = []
        steps = p._steps_to_run()
        root = p._source
        prev: tuple[str, str] | None = None
        idx = -1
        if isinstance(root, _Command):
            info, kw = root.obj
            real_kw, shown_kw = self.kwargs(kw)
            idx = len(self.stages)
            to, to_shown = self._to(info, idx, not steps, out, shown_out)
            shown = _native.source_argv(info["path"], shown_kw, to_shown, check=False)
            argv = _native.source_argv(info["path"], real_kw, to) if self.run else None
            self.add(argv, shown)
            lines.append(self.stages[idx]["command_line"])
            prev = (_pipe(idx), "-")
        elif isinstance(root, Nested):
            j, nested = self.chain(root.obj, None, [])
            lines += nested
            prev = (_pipe(j), "-")
        for i, step in enumerate(steps):
            b = step.binding
            first = prev if prev is not None else self.input(root)
            extras = [self.input(s) for s in step.extra]
            real_kw, shown_kw = self.kwargs(step.kwargs)
            idx = len(self.stages)
            to, to_shown = self._to(b, idx, i == len(steps) - 1, out, shown_out)
            shown = _native.step_argv(
                b["path"], [first[1], *(s for _, s in extras)], shown_kw, to_shown, check=False
            )
            argv = None
            if self.run:
                ins = [first[0], *(r for r, _ in extras)]
                argv = _native.step_argv(b["path"], ins, real_kw, to)
            self.add(argv, shown)
            lines.append(self.stages[idx]["command_line"])
            prev = (_pipe(idx), "-")
        return idx, lines

    def run_plan(self, terminal: str) -> Any:
        """Run the plan; a one-shot source is spent from now on, whatever the result."""
        for s in self._once:
            s.used = True
        plan = {"stages": self.stages, "sources": self.sources, "terminal": terminal, "names": None}
        return _native.run_plan(plan)


class Pipeline:
    """A lazy pipeline of gloci steps. Start one with `gloci.open()` (or `gloci.windows()`,
    `gloci.random()`).

    Each step method (`intersect`, `merge`, `vcf_filter`, ...) returns a new pipeline; the
    keyword arguments are the command's flags with `-` written as `_`. A keyword argument that
    takes a file (`b=`) also takes any source: a file object, an iterable, an Arrow table, or
    another pipeline. Terminal methods run it: `save()`, `to_bytes()` and `count()`. `cli()` is
    the equivalent shell pipeline.
    """

    __slots__ = ("_format", "_source", "_steps")

    _source: Source
    _format: str | None
    _steps: tuple[_Step, ...]

    def __init__(self) -> None:
        raise TypeError("start a pipeline with gloci.open()")

    @classmethod
    def _make(cls, source: Source, fmt: str | None, steps: tuple[_Step, ...]) -> Pipeline:
        p = object.__new__(cls)
        object.__setattr__(p, "_source", source)
        object.__setattr__(p, "_format", fmt)
        object.__setattr__(p, "_steps", steps)
        return p

    def __setattr__(self, name: str, value: object) -> None:
        raise AttributeError(
            f"a Pipeline is immutable, so {name!r} cannot be set; each step returns a new one"
        )

    def __delattr__(self, name: str) -> None:
        raise AttributeError(f"a Pipeline is immutable, so {name!r} cannot be deleted")

    def __repr__(self) -> str:
        return f"<gloci.Pipeline: {self.cli()}>"

    # Steps

    def _then(self, b: dict[str, Any], sources: tuple[object, ...], kwargs: dict[str, Any]):
        method = b["method"]
        if sources and not b["nary"]:
            raise TypeError(
                f"{method}() takes no positional arguments: pass its options as keyword arguments"
            )
        if b["needs_file"]:
            if self._steps or isinstance(self._source, (Nested, _Command)):
                raise UsageError(
                    f"{method} needs a file input; save() the previous step to a file and open it"
                )
            if not isinstance(self._source, Path):
                raise UsageError(
                    f"{method} needs a file input, not {self._source.display}; write it to a "
                    "file and open that"
                )
        extra = tuple(classify_source(s, f"{method}()") for s in sources)
        kw = {name: _kwarg(f"{method}()", name, value) for name, value in kwargs.items()}
        if not self._steps and self._format is not None:
            kw = self._with_format(b, kw)
        _native.kwargs_to_argv(b["path"], _placeholders(kw))  # raises UsageError at once (D85)
        return Pipeline._make(self._source, self._format, (*self._steps, _Step(b, kw, extra)))

    def _with_format(self, b: dict[str, Any], kwargs: dict[str, Any]) -> dict[str, Any]:
        """`kwargs` with the format given to `open()`, for the step that reads the source."""
        name = b["input_format"]
        if name is None:
            raise UsageError(
                f"{b['method']}() cannot take open(format=...): gloci {' '.join(b['path'])} has "
                "no option for its input's format"
            )
        if name in kwargs:
            raise UsageError(
                f"open(format=...) and {name}= both set the input's format; pass one of them"
            )
        return {name: self._format, **kwargs}

    # Terminals

    def _steps_to_run(self) -> tuple[_Step, ...]:
        """The steps; a bare source is read with `view`, which writes it as it is, and a module
        source (`gloci.windows()`) is its own command."""
        if self._steps:
            return self._steps
        if isinstance(self._source, _Command):
            return ()
        view = _BINDINGS["view"]
        kw = self._with_format(view, {}) if self._format is not None else {}
        return (_Step(view, kw, ()),)

    def _last(self) -> dict[str, Any]:
        """The binding (or module source) of the command whose output the terminal takes."""
        steps = self._steps_to_run()
        if steps:
            return steps[-1].binding
        info: dict[str, Any] = self._source.obj[0]
        return info

    def _plan(
        self, out: str | None, save_opts: list[str], threads: object, run: bool = True
    ) -> _Plan:
        """The run's plan: the last stage writes the file `out` with `save_opts` (`save()`), or
        the terminal's output (`out` None)."""
        plan = _Plan(threads, run)
        if out is None:
            _, plan.lines = plan.chain(self, ["-o", _OUT], [])
        else:
            words = ["-o", out, *save_opts]
            _, plan.lines = plan.chain(self, words, words)
        return plan

    def _stages(
        self, out: str | None, save_opts: list[str], threads: object, run: bool = True
    ) -> list[dict[str, Any]]:
        """Each stage's argv (when `run`) and command line, nested pipelines' first."""
        return self._plan(out, save_opts, threads, run).stages

    def _run(self, terminal: str, out: str | None, save_opts: list[str], threads: object) -> Any:
        return self._plan(out, save_opts, threads).run_plan(terminal)

    def save(
        self,
        path: str | os.PathLike[str],
        *,
        output_format: str | None = None,
        write_index: bool = False,
        compression_level: int | None = None,
        threads: int | None = None,
    ) -> None:
        """Run the pipeline into the file `path`, whose extension gives the format (`.bed.gz`
        is BGZF-compressed BED), as `-o path` does. The file appears only once the run has
        succeeded.

        Args:
            output_format: The format, overriding the extension (`--output-format`).
            write_index: Also write an index next to the file (`--write-index`).
            compression_level: The BGZF compression level, 0-12 (`--compression-level`).
            threads: Threads per step (`-t`); the command line's default when None.
        """
        if not isinstance(path, (str, os.PathLike)):
            raise TypeError(f"save() takes a path (str or os.PathLike), not {type(path).__name__}")
        out = classify_source(path, "save()")
        last = self._last()
        given: list[tuple[str, list[str]]] = []
        if output_format is not None:
            given.append(("output_format", ["--output-format", str(output_format)]))
        if write_index:
            given.append(("write_index", ["--write-index"]))
        if compression_level is not None:
            given.append(("compression_level", ["--compression-level", str(compression_level)]))
        opts: list[str] = []
        for name, words in given:
            if name not in last["outputs"]:
                raise UsageError(
                    f"save(): {last['method']}() has no {name}= (gloci "
                    f"{' '.join(last['path'])} has no {words[0]} option)"
                )
            opts += words
        self._run("save", out.obj, opts, threads)

    def to_bytes(self, *, threads: int | None = None) -> bytes:
        """Run the pipeline and return what the command line writes to stdout."""
        result: bytes = self._run("bytes", None, [], threads)
        return result

    def count(self, *, threads: int | None = None) -> int:
        """Run the pipeline and return the number of records it writes (header lines, such as
        `#` lines and SAM's `@` lines, are not records)."""
        result: int = self._run("count", None, [], threads)
        return result

    def cli(self) -> str:
        """The equivalent shell pipeline, for stdout: `gloci intersect -a a.bed -b b.bed |
        gloci merge -`. A step reading the previous one's output reads `-`, a Python object is
        `<python object N>`, and a pipeline given as a source is bash's `<(gloci ...)`."""
        return " | ".join(self._plan(None, [], None, run=False).lines)


def _step_method(b: dict[str, Any]):  # type: ignore[no-untyped-def]
    method = b["method"]

    def step(self: Pipeline, *sources: object, **kwargs: Any) -> Pipeline:
        return self._then(b, sources, kwargs)

    step.__name__ = method
    step.__qualname__ = f"Pipeline.{method}"
    step.__doc__ = f"{b['about']}\n\nCLI: gloci {' '.join(b['path'])}"
    return step


for _b in _BINDINGS.values():
    if hasattr(Pipeline, _b["method"]):
        raise RuntimeError(f"the step method {_b['method']}() clashes with a Pipeline attribute")
    setattr(Pipeline, _b["method"], _step_method(_b))
del _b


def open(src: object, *, format: str | None = None) -> Pipeline:  # noqa: A001
    """Start a pipeline that reads `src`. Nothing is read until a terminal method runs.

    Args:
        src: A path (`-` is standard input); a file object (bytes, or text read as UTF-8); an
            iterable of lines (str) or of tuples of values (tab-joined `str()` values, None as
            `.`); an Arrow table or stream (a pyarrow Table, a polars or pandas DataFrame:
            the integer columns `chrom start end` first, then the others); or a Pipeline.
            Python objects are read as the run goes; an iterator or a file object can be read
            by one run only.
        format: The input's format (bed, vcf, ...), overriding detection; the first step
            takes it as its input-format option (`--format`, `--a-format`).
    """
    source = classify_source(src, "gloci.open()")
    if isinstance(source, Nested) and format is None:
        return source.obj  # type: ignore[no-any-return]
    return Pipeline._make(source, None if format is None else str(format), ())
