"""Lazy, immutable pipelines over the gloci commands (D85-D87).

A pipeline is a source and a list of steps. Each step records a command and its keyword
arguments, checked when the step is called; nothing runs until a terminal method. A terminal
runs every step at once, each on its own thread, chained by bounded in-memory pipes, in the
same code as the command line: the output is the bytes the shell pipeline `cli()` prints would
write.
"""

from __future__ import annotations

import os
from typing import Any

from . import _native
from .errors import UsageError

__all__ = ["Pipeline", "open"]

#: Each step method's binding, by method name (`_native.bindings()`).
_BINDINGS: dict[str, dict[str, Any]] = {b["method"]: b for b in _native.bindings()}

# Plan tokens (`_native.run_plan`): the pipe stage `i` writes, and the terminal's output.
_OUT = "\0out"


def _pipe(i: int) -> str:
    return f"\0pipe:{i}"


def _path(value: object, what: str) -> str:
    """`value` as a path: a `str`, or an `os.PathLike`."""
    if isinstance(value, str):
        path = value
    elif isinstance(value, os.PathLike):
        path = os.fsdecode(os.fspath(value))
    else:
        raise TypeError(f"{what} takes a path (str or os.PathLike), not {type(value).__name__}")
    if "\0" in path:
        raise UsageError(f"{what}: the path {path!r} holds a NUL character")
    return path


def _check_value(method: str, name: str, value: object) -> None:
    items = value if isinstance(value, (list, tuple)) else (value,)
    for item in items:
        if isinstance(item, Pipeline):
            raise TypeError(
                f"{method}(): {name}= takes paths and values; a Pipeline as an input is not "
                "supported yet, so save() it to a file and pass the file"
            )


class _Step:
    """One step: the binding of its method, its keyword arguments and its extra sources."""

    __slots__ = ("binding", "extra", "kwargs")

    def __init__(self, binding: dict[str, Any], kwargs: dict[str, Any], extra: tuple[str, ...]):
        self.binding = binding
        self.kwargs = kwargs
        self.extra = extra


class Pipeline:
    """A lazy pipeline of gloci steps. Start one with `gloci.open()`.

    Each step method (`intersect`, `merge`, `vcf_filter`, ...) returns a new pipeline; the
    keyword arguments are the command's flags with `-` written as `_`. Terminal methods run
    it: `save()`, `to_bytes()` and `count()`. `cli()` is the equivalent shell pipeline.
    """

    __slots__ = ("_format", "_source", "_steps")

    _source: str
    _format: str | None
    _steps: tuple[_Step, ...]

    def __init__(self) -> None:
        raise TypeError("start a pipeline with gloci.open()")

    @classmethod
    def _make(cls, source: str, fmt: str | None, steps: tuple[_Step, ...]) -> Pipeline:
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
        if b["needs_file"] and self._steps:
            raise UsageError(
                f"{method} needs a file input; save() the previous step to a file and open it"
            )
        extra = tuple(_path(s, f"{method}()") for s in sources)
        for name, value in kwargs.items():
            _check_value(method, name, value)
        kw = dict(kwargs)
        if not self._steps and self._format is not None:
            kw = self._with_format(b, kw)
        _native.kwargs_to_argv(b["path"], kw)  # raises UsageError at once (D85)
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
        """The steps; a bare source is read with `view`, which writes it as it is."""
        if self._steps:
            return self._steps
        view = _BINDINGS["view"]
        kw = self._with_format(view, {}) if self._format is not None else {}
        return (_Step(view, kw, ()),)

    def _stages(
        self,
        out: str | None,
        save_opts: list[str],
        threads: object,
        run: bool = True,
    ) -> list[dict[str, Any]]:
        """Each step's argv (when `run`) and command line. `out` is the file `save()` writes
        (`None` for the other terminals), and `save_opts` its other options."""
        steps = self._steps_to_run()
        stages = []
        for i, step in enumerate(steps):
            b = step.binding
            last = i == len(steps) - 1
            shown_in = self._source if i == 0 else "-"
            shown_out = ["-o", out, *save_opts] if last and out is not None else []
            shown = _native.step_argv(
                b["path"], [shown_in, *step.extra], step.kwargs, shown_out, check=False
            )
            stage: dict[str, Any] = {"command_line": _native.command_line(b["path"], shown)}
            if run:
                real_in = self._source if i == 0 else _pipe(i - 1)
                if not last:
                    to = ["-o", _pipe(i)]
                elif out is None:
                    to = ["-o", _OUT]
                else:
                    to = ["-o", out, *save_opts]
                if threads is not None and b["threads"]:
                    to += ["-t", str(threads)]
                stage["argv"] = _native.step_argv(
                    b["path"], [real_in, *step.extra], step.kwargs, to
                )
            stages.append(stage)
        return stages

    def _run(self, terminal: str, out: str | None, save_opts: list[str], threads: object) -> Any:
        plan = {
            "stages": self._stages(out, save_opts, threads),
            "sources": {},
            "terminal": terminal,
            "names": None,
        }
        return _native.run_plan(plan)

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
            output_format: The format, overriding the extension (`-O`).
            write_index: Also write an index next to the file (`--write-index`).
            compression_level: The BGZF compression level, 0-12 (`--compression-level`).
            threads: Threads per step (`-t`); the command line's default when None.
        """
        out = _path(path, "save()")
        last = self._steps_to_run()[-1].binding
        opts: list[str] = []
        for name, flag, value in (
            ("output_format", "-O", output_format),
            ("write_index", "--write-index", write_index or None),
            ("compression_level", "--compression-level", compression_level),
        ):
            if value is None:
                continue
            if name not in last["outputs"]:
                raise UsageError(
                    f"save(): {last['method']}() has no {name}= (gloci "
                    f"{' '.join(last['path'])} has no {flag} option)"
                )
            opts += [flag] if value is True else [flag, str(value)]
        self._run("save", out, opts, threads)

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
        gloci merge -`. A step reading the previous one's output reads `-`."""
        return " | ".join(s["command_line"] for s in self._stages(None, [], None, run=False))


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


def open(src: str | os.PathLike[str], *, format: str | None = None) -> Pipeline:  # noqa: A001
    """Start a pipeline that reads the file `src`. Nothing is read until a terminal method runs.

    Args:
        src: The path; `-` is standard input.
        format: The input's format (bed, vcf, ...), overriding detection; the first step
            takes it as its input-format option (`--format`, `--a-format`).
    """
    return Pipeline._make(_path(src, "gloci.open()"), None if format is None else str(format), ())
