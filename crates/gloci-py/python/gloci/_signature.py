"""Runtime `__signature__` and `__doc__` of the generated step methods and module functions
(P2b Task 9, R3).

They are built from the same specs as the type stub (`cargo xtask pystub` writes it from
`gloci::embed`), which `_native.signatures()` and `_native.command_specs()` report: the parameter
names and kinds, and the docstrings, are the stub's. `test_api_surface.py` compares them.
"""

from __future__ import annotations

import inspect
from typing import Any

from . import _native

__all__ = ["make_doc", "make_signature"]

_P = inspect.Parameter

_SPECS = _native.command_specs()

#: Each command's arguments by clap id: `{path: {id: spec}}`.
_ARGS: dict[tuple[str, ...], dict[str, dict[str, Any]]] = {
    tuple(spec["path"]): {a["id"]: a for a in spec["args"]} for spec in _SPECS
}
_ABOUT: dict[tuple[str, ...], str] = {tuple(spec["path"]): spec["about"] for spec in _SPECS}

#: Each step method's and function's parameters (`_native.signatures()`).
SIGNATURES: dict[str, dict[str, Any]] = _native.signatures()


def _squash(text: str) -> str:
    return " ".join(text.split())


def make_doc(path: tuple[str, ...], sig: dict[str, Any]) -> str:
    """The docstring: the command's about text, `Args:` with each argument's help, and the
    `CLI: gloci <path>` line."""
    docs: list[tuple[str, str]] = [(name, help_) for name, help_ in sig["lead"]]
    if sig["out"] is not None:
        docs.append(("out", sig["out"]))
    for k in sig["args"]:
        help_ = f"alias of {k['alias_of']}" if k["alias_of"] else _ARGS[path][k["id"]]["help"]
        docs.append((k["name"], help_))
    text = _squash(_ABOUT[path])
    if docs:
        text += "\n\nArgs:"
        for name, help_ in docs:
            help_ = _squash(help_)
            text += f"\n    {name}: {help_}" if help_ else f"\n    {name}:"
    if sig["skipped"]:
        first = sig["skipped"][0]
        names = ", ".join(k["name"] for k in sig["skipped"])
        value = "True" if first["switch"] else '"..."'
        text += (
            f"\n\nNot Python identifiers, so not parameters: {names}. "
            f'Pass one as **{{"{first["name"]}": {value}}}.'
        )
    return f"{text}\n\nCLI: gloci {' '.join(path)}"


def make_signature(path: tuple[str, ...], sig: dict[str, Any], method: bool) -> inspect.Signature:
    """The signature: `self` for a method; a function's leading positional parameters; the extra
    sources of an n-ary step as `*sources`; then the keyword-only parameters (`out`, required, for
    a function that writes files; each argument of the command that a keyword may set, then its
    short flag and aliases). A switch defaults to False, an argument that is not required to
    None."""
    params: list[inspect.Parameter] = []
    if method:
        params.append(_P("self", _P.POSITIONAL_OR_KEYWORD))
    params += [_P(name, _P.POSITIONAL_OR_KEYWORD) for name, _ in sig["lead"]]
    if sig["nary"]:
        params.append(_P("sources", _P.VAR_POSITIONAL))
    if sig["out"] is not None:
        params.append(_P("out", _P.KEYWORD_ONLY))
    for k in sig["args"]:
        if not _ARGS[path][k["id"]]["takes_value"]:
            default: Any = False
        elif k["required"]:
            default = _P.empty
        else:
            default = None
        params.append(_P(k["name"], _P.KEYWORD_ONLY, default=default))
    return inspect.Signature(params)
