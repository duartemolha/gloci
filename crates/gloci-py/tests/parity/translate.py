"""A CLI argv as the Python call that asks for the same output (P2b Task 10, ruling R1).

`translate(args)` reads a gloci argv (without `gloci`) with the command's own spec
(`_native.command_specs()`): which words are options, which take values and how many, which
are positionals. It splits them into the pipeline's input (the argument the step binds,
`_native.bindings()`), the other inputs (an n-ary step's extra sources, a module function's
positionals) and the keyword arguments, named canonically (`min_frac_a`, never `f`). The
result is a `Call`, which runs the call when called and returns the output's bytes:

- a step: `gloci.open(input).<method>(*extra, **kwargs)`, then `to_bytes(threads=)`, or
  `save(path, output_format=..., ...)` when the argv has `-o path` or an output option
  (`-O sam` to stdout: saved to a scratch file, whose bytes are the output);
- a module source (`windows -g`, `random`): `gloci.windows(**kwargs)` and the same terminals;
- an eager module function (`vcf stats`, `vcf extract --genotype-matrix`, `vcf compare`,
  `vcf popgen`, `index`): the function, with `--out-prefix` as `out=`; it writes files, whose
  paths are kept in `Call.written`, and the output (its stdout) is empty.

Standard input (`-`) is the file given as `stdin`, as a file object. `vcf merge --file-list
F` passes F's files as the step's sources, which is how a step takes several inputs.

What stays on the command line raises `UntranslatableError` with `cli_only` set: a command
D87 lists as `CLI_ONLY` (neither a step nor a module function), and an argument a step never
takes that no module function offers instead (`sort --check`: it writes no output for a
terminal to take, P2b Task 7 R2). Anything else untranslatable (a word the spec does not
know, a step with no input, standard input with no file for it) has `cli_only` unset.
"""

from __future__ import annotations

import os
import pathlib
from collections.abc import Callable
from typing import Any

import gloci
from gloci import _native

__all__ = ["Call", "UntranslatableError", "split_args", "translate"]

#: Each command's spec, by path.
SPECS: dict[tuple[str, ...], dict[str, Any]] = {
    tuple(s["path"]): s for s in _native.command_specs()
}
#: The step methods, by command path.
BINDINGS: dict[tuple[str, ...], dict[str, Any]] = {tuple(b["path"]): b for b in _native.bindings()}
#: The module functions and sources, by command path.
FUNCTIONS: dict[tuple[str, ...], dict[str, Any]] = {
    tuple(f["path"]): f for f in _native.module_functions()
}
_SIGNATURES: dict[str, Any] = _native.signatures()

#: The options a step leaves to its terminal method (`save()`, `to_bytes()`), by clap id.
_TERMINAL = ("output", "output_format", "write_index", "compression_level", "threads")
#: The global log-level switches: they change what the CLI writes to stderr only, and no
#: keyword sets them (the level is the `gloci` logger's).
_LOG_LEVEL = ("quiet", "verbose")
#: The scratch file a step's output goes to when the argv has an output option but no `-o`.
SCRATCH_OUTPUT = "parity-output"


class UntranslatableError(Exception):
    """An argv with no Python call that asks for the same thing.

    Attributes:
        cli_only: The command stays on the command line (D87's `CLI_ONLY`).
    """

    def __init__(self, message: str, *, cli_only: bool = False) -> None:
        super().__init__(message)
        self.cli_only = cli_only


def split_args(line: str) -> list[str]:
    """A golden `args` line as words, as `gloci_testkit::differential::split_args` splits it:
    whitespace separates, and a single-quoted run is text, quotes dropped."""
    out: list[str] = []
    cur: list[str] = []
    in_word = quoted = False
    for c in line:
        if c == "'":
            quoted = not quoted
            in_word = True
        elif c.isspace() and not quoted:
            if in_word:
                out.append("".join(cur))
                cur = []
                in_word = False
        else:
            cur.append(c)
            in_word = True
    if quoted:
        raise ValueError(f"unbalanced quote in {line!r}")
    if in_word:
        out.append("".join(cur))
    return out


def is_cli_only(path: tuple[str, ...]) -> bool:
    """Whether the command at `path` is neither a step nor a module function: CLI_ONLY (D87)
    and the commands not placed yet (`embed::UNBOUND`, which
    `every_cli_command_is_bound_or_excluded` checks)."""
    return path in SPECS and path not in BINDINGS and path not in FUNCTIONS


class _Parsed:
    """An argv read with its command's spec: each option's values by clap id, in the order
    first given (a switch's value is True), and the positionals."""

    def __init__(self, path: tuple[str, ...], words: list[str]) -> None:
        self.path = path
        self.spec = SPECS[path]
        self.options: dict[str, list[str] | bool] = {}
        self.positionals: list[str] = []
        self._read(words)

    def arg(self, arg_id: str) -> dict[str, Any]:
        return next(a for a in self.spec["args"] if a["id"] == arg_id)

    def _put(self, a: dict[str, Any], values: list[str] | None) -> None:
        if values is None:
            self.options[a["id"]] = True
            return
        old = self.options.get(a["id"])
        if a["multiple"] and isinstance(old, list):
            old.extend(values)
        else:
            # clap keeps the last of an option given twice that is not repeatable
            self.options[a["id"]] = list(values)

    def _read(self, words: list[str]) -> None:
        cmd = " ".join(self.path)
        longs: dict[str, dict[str, Any]] = {}
        shorts: dict[str, dict[str, Any]] = {}
        for a in self.spec["args"]:
            if a["positional"]:
                continue
            for name in [a["long"], *a["aliases"]]:
                if name:
                    longs[name] = a
            if a["short"]:
                shorts[a["short"]] = a
        i = 0
        only_positionals = False
        while i < len(words):
            w = words[i]
            i += 1
            if only_positionals or w == "-" or not w.startswith("-"):
                self.positionals.append(w)
                continue
            if w == "--":
                only_positionals = True
                continue
            attached: str | None = None
            if w.startswith("--"):
                name, eq, value = w[2:].partition("=")
                a = longs.get(name)
                if a is None:
                    raise UntranslatableError(
                        f"--{name} is not an argument gloci {cmd} offers to a host (hidden?)"
                    )
                if eq:
                    attached = value
            else:
                # -x, -xVALUE, -x=VALUE, or a cluster of switches -xyz
                j = 1
                while True:
                    a = shorts.get(w[j])
                    if a is None:
                        raise UntranslatableError(
                            f"-{w[j]} (in {w}) is not an argument gloci {cmd} offers to a host"
                        )
                    rest = w[j + 1 :]
                    if a["takes_value"] or not rest:
                        if rest:
                            attached = rest[1:] if rest.startswith("=") else rest
                        break
                    self._put(a, None)
                    j += 1
            if not a["takes_value"]:
                if attached is not None:
                    raise UntranslatableError(f"{w}: a switch with a value")
                self._put(a, None)
                continue
            if attached is not None:
                self._put(a, [attached])
                continue
            if a["require_equals"]:
                # `--split` alone: its default-missing value, which `split=True` gives too
                self._put(a, None)
                continue
            if i >= len(words):
                raise UntranslatableError(f"{w} has no value")
            # The CLI took this argv, so the next word is the value, whatever it starts with
            # (`--by -100`, `--strand -`).
            values = [words[i]]
            i += 1
            if a["max_values"] is None:
                # `-b x y`: clap takes values until the next word that looks like a flag
                while i < len(words) and (words[i] == "-" or not words[i].startswith("-")):
                    values.append(words[i])
                    i += 1
            self._put(a, values)

    def positional_values(self) -> dict[str, list[str]]:
        """The positionals by clap id, as clap hands them out: in order, one each to a single
        positional, all that are left to a multiple one (less one for each single positional
        after it)."""
        args = [a for a in self.spec["args"] if a["positional"]]
        words = list(self.positionals)
        out: dict[str, list[str]] = {}
        for k, a in enumerate(args):
            if not words:
                break
            if a["multiple"]:
                after = sum(1 for b in args[k + 1 :] if not b["multiple"])
                n = max(len(words) - after, 1)
            else:
                n = 1
            out[a["id"]], words = words[:n], words[n:]
        if words:
            raise UntranslatableError(
                f"gloci {' '.join(self.path)} takes no positional {words[0]!r}"
            )
        return out


def _value(a: dict[str, Any], values: list[str] | bool) -> Any:
    """A keyword argument's value: True for a switch, a str, or a list of str when the
    argument was given several values."""
    if values is True:
        return True
    assert isinstance(values, list)
    return values[0] if len(values) == 1 else list(values)


def _command(words: list[str]) -> tuple[tuple[str, ...], list[str]]:
    for n in (2, 1):
        path = tuple(words[:n])
        if len(path) == n and path in SPECS:
            return path, words[n:]
    raise UntranslatableError(f"no gloci command in {' '.join(words[:2])!r}")


class Stdin:
    """Standard input (`-`): the file standing for it, opened as a file object when the call
    runs (R1)."""

    def __init__(self, path: str | os.PathLike[str]) -> None:
        self.path = path

    def __repr__(self) -> str:
        return f"open({os.fspath(self.path)!r}, 'rb')"

    def open(self) -> Any:
        return open(self.path, "rb")  # noqa: SIM115 - the run reads it to its end


def _input(i: Any) -> Any:
    return i.open() if isinstance(i, Stdin) else i


class Call:
    """One Python call (`translate()`), run by calling it: the output's bytes."""

    def __init__(
        self,
        argv: list[str],
        kind: str,
        name: str,
        inputs: list[Any],
        kwargs: dict[str, Any],
        terminal: dict[str, Any],
        output: str | None,
    ) -> None:
        #: The argv translated (without `gloci`).
        self.argv = argv
        #: "step", "source" or "function".
        self.kind = kind
        #: The step method or module function.
        self.name = name
        #: The inputs: the pipeline's input and a step's extra sources, or a function's
        #: positionals. A path, or `Stdin`.
        self.inputs = inputs
        self.kwargs = kwargs
        #: `save()`'s options (`output_format`, `write_index`, `compression_level`) and
        #: `threads`.
        self.terminal = terminal
        #: The file the output goes to (`-o`), else None (stdout).
        self.output = output
        #: The files an eager function reported writing, once it has run.
        self.written: list[pathlib.Path] = []
        #: A file listing the inputs (`--file-list`), read when the call runs, else None.
        self.file_list: str | None = None

    def __repr__(self) -> str:
        return self.python()

    def python(self) -> str:
        """The call as Python source, for messages."""
        kw = ", ".join(f"{k}={v!r}" for k, v in self.kwargs.items())
        ins = [repr(i) for i in self.inputs]
        if self.file_list is not None:
            ins = [f"*<the lines of {self.file_list}>"]
        if self.kind == "function":
            return f"gloci.{self.name}({', '.join([*ins, *([kw] if kw else [])])})"
        if self.kind == "source":
            head = f"gloci.{self.name}({kw})"
        else:
            args = ", ".join([*ins[1:], *([kw] if kw else [])])
            head = f"gloci.open({ins[0]}).{self.name}({args})"
        return f"{head}.{self._terminal_text()}"

    def _terminal_text(self) -> str:
        threads = self.terminal.get("threads")
        save = {k: v for k, v in self.terminal.items() if k != "threads"}
        if self.output is None and not save:
            return f"to_bytes({'' if threads is None else f'threads={threads}'})"
        out = self.output or SCRATCH_OUTPUT
        opts = "".join(f", {k}={v!r}" for k, v in self.terminal.items())
        return f"save({out!r}{opts})"

    def pipeline(self) -> gloci.Pipeline:
        """The pipeline of a step or module source."""
        if self.kind == "source":
            source: Callable[..., gloci.Pipeline] = getattr(gloci, self.name)
            return source(**self.kwargs)
        inputs = [_input(i) for i in self.inputs]
        if self.file_list is not None:
            text = pathlib.Path(self.file_list).read_text()
            inputs = [line.removesuffix("\r") for line in text.split("\n") if line.strip("\r")]
        head = gloci.open(inputs[0])
        return getattr(head, self.name)(*inputs[1:], **self.kwargs)

    def __call__(self) -> bytes:
        if self.kind == "function":
            f: Callable[..., list[pathlib.Path]] = getattr(gloci, self.name)
            self.written = f(*map(_input, self.inputs), **self.kwargs)
            return b""
        p = self.pipeline()
        threads = self.terminal.get("threads")
        save = {k: v for k, v in self.terminal.items() if k != "threads"}
        if self.output is None and not save:
            return p.to_bytes(threads=threads)
        out = pathlib.Path(self.output or SCRATCH_OUTPUT)
        p.save(out, threads=threads, **save)
        return out.read_bytes()


def _stdin(word: str, stdin: str | os.PathLike[str] | None) -> Any:
    if word != "-":
        return word
    if stdin is None:
        raise UntranslatableError("reads standard input, and the case gives none")
    return Stdin(stdin)


def _bound_input(b: dict[str, Any], parsed: _Parsed) -> str | None:
    """The clap id of the argument a step's input fills (`embed::input_arg`)."""
    args = parsed.spec["args"]
    if b["input"] == "flag":
        f = b["flag"]
        return next(
            a["id"]
            for a in args
            if not a["positional"] and f in (a["short"], a["long"], a["kwarg"])
        )
    positionals = [a["id"] for a in args if a["positional"]]
    if b["input"] == "positional":
        return positionals[0]
    before = positionals.index(b["flag"])
    return positionals[before + 1]


def _uses_function(path: tuple[str, ...], parsed: _Parsed, given: set[str]) -> bool:
    """Whether the argv asks for the module function of a command that is also a step: it
    gives an argument the function takes and the step does not (`vcf extract
    --genotype-matrix`), or it lacks the step's input (`windows -g`)."""
    f = FUNCTIONS.get(path)
    b = BINDINGS.get(path)
    if f is None:
        return False
    if b is None:
        return True
    step = {a["id"] for a in _SIGNATURES[b["method"]]["args"]}
    function = {a["id"] for a in _SIGNATURES[f["function"]]["args"]}
    if any(i in function and i not in step and i not in _TERMINAL for i in given):
        return True
    return _bound_input(b, parsed) not in given


def _check_not_reserved(b: dict[str, Any], parsed: _Parsed, bound: str | None) -> None:
    """Raise `UntranslatableError` (`cli_only`) for an argument the step never takes and no
    module function offers instead (`sort --check`, which writes no output for a terminal to
    take; P2b Task 7 R2): it stays on the command line."""
    sig = _SIGNATURES[b["method"]]
    params = {a["id"] for a in sig["args"]}
    unnamed = {s["name"] for s in sig["skipped"]}  # `5p`: not identifiers, still kwargs
    for arg_id in parsed.options:
        a = parsed.arg(arg_id)
        if arg_id in params or a["kwarg"] in unnamed or arg_id == bound:
            continue
        if arg_id in _TERMINAL or arg_id in _LOG_LEVEL:
            continue
        if arg_id == "file_list" and b["nary"]:
            continue  # the inputs: the translation reads the list, and passes them as sources
        try:
            _native.kwargs_to_argv(list(parsed.path), {a["kwarg"]: True})
            why = "not a parameter of the step"
        except gloci.UsageError as e:
            why = str(e)
        flag = f"--{a['long']}" if a["long"] else f"-{a['short']}"
        raise UntranslatableError(
            f"gloci {' '.join(parsed.path)} {flag} stays on the command line: {why}",
            cli_only=True,
        )


def _terminal(parsed: _Parsed, options: dict[str, Any]) -> dict[str, Any]:
    """The terminal method's options: `threads`, and `save()`'s output options."""
    out = {}
    for t in _TERMINAL:
        if t != "output" and t in options:
            v = _value(parsed.arg(t), options[t])
            # `-t 4`, `--compression-level 1`: the terminals take numbers
            out[t] = int(v) if isinstance(v, str) and v.isdigit() else v
    return out


def _output(parsed: _Parsed, options: dict[str, Any]) -> str | None:
    """The file `-o` names, if any."""
    return _value(parsed.arg("output"), options["output"]) if "output" in options else None


def translate(args: list[str] | str, *, stdin: str | os.PathLike[str] | None = None) -> Call:
    """The Python call for the gloci argv `args` (words, or a golden `args` line), which reads
    its inputs by the same relative paths. `stdin` is the file standing for standard input
    (`-`). Raises `UntranslatableError`."""
    words = split_args(args) if isinstance(args, str) else list(args)
    path, rest = _command(words)
    if is_cli_only(path):
        raise UntranslatableError(
            f"gloci {' '.join(path)} stays on the command line (CLI_ONLY)", cli_only=True
        )
    parsed = _Parsed(path, rest)
    positional = parsed.positional_values()
    given = set(parsed.options) | set(positional)
    options = dict(parsed.options)
    for level in _LOG_LEVEL:
        options.pop(level, None)

    def kwargs(skip: set[str]) -> dict[str, Any]:
        kw: dict[str, Any] = {}
        for arg_id, values in [*options.items(), *positional.items()]:
            if arg_id in skip:
                continue
            a = parsed.arg(arg_id)
            kw[a["kwarg"]] = _value(a, values)
        return kw

    if _uses_function(path, parsed, given):
        f = FUNCTIONS[path]
        if f["source"]:
            return Call(
                words,
                "source",
                f["function"],
                [],
                kwargs(set(_TERMINAL)),
                _terminal(parsed, options),
                _output(parsed, options),
            )
        inputs = [_stdin(w, stdin) for i in f["inputs"] for w in positional.get(i, [])]
        kw = kwargs(set(f["inputs"]))
        if "out_prefix" in kw:
            kw["out"] = kw.pop("out_prefix")
        return Call(words, "function", f["function"], inputs, kw, {}, None)

    b = BINDINGS[path]
    bound = _bound_input(b, parsed)
    _check_not_reserved(b, parsed, bound)
    values = options.get(bound) if b["input"] == "flag" else positional.get(bound)
    if values is None and b["nary"] and "file_list" in options:
        # `vcf merge --file-list list.txt`: a step's inputs are the pipeline and the sources
        # passed to it, so the list's files are those, read when the call runs.
        listed = _value(parsed.arg("file_list"), options["file_list"])
        terminal = _terminal(parsed, options)
        skip = {bound, "file_list", *_TERMINAL}
        output = _output(parsed, options)
        call = Call(words, "step", b["method"], [], kwargs(skip), terminal, output)
        call.file_list = listed
        return call
    if not isinstance(values, list) or not values:
        raise UntranslatableError(
            f"gloci {' '.join(path)} without its input ({bound}): a step reads the pipeline it "
            "is called on"
        )
    if b["input"] == "flag" and len(values) != 1:
        raise UntranslatableError(f"{bound} is given {len(values)} times")
    inputs = [_stdin(w, stdin) for w in values]
    terminal = _terminal(parsed, options)
    output = _output(parsed, options)
    if output is None and set(terminal) - {"threads"} and "output_format" not in terminal:
        # Saved to the scratch file, which has no extension: name the format the command
        # writes to standard output when a flag rather than `-O` chooses it.
        if path == ("convert",) and options.get("bam") is True:
            terminal["output_format"] = "bam"
    return Call(
        words,
        "step",
        b["method"],
        inputs,
        kwargs({bound, *_TERMINAL}),
        terminal,
        output,
    )
