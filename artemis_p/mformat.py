"""MATLAB-compatible sprintf.

MATLAB's sprintf/fprintf differ from C's printf in ways that show up in the
PHREEQC input files that getphreeqcrun.m writes:

* arguments are flattened (column-major) and the format is recycled until
  they are used up; output stops at the first conversion with no data left;
* ``%d`` with a non-integer value prints in ``%e`` form (``5.500000e+00``);
* NaN and Inf print as ``NaN`` and ``Inf``.
"""

import math
import re

import numpy as np

_SPEC = re.compile(r"%(?P<flags>[-+ 0#]*)(?P<width>\d+)?(?:\.(?P<prec>\d+))?(?P<conv>[diouxXfeEgGcs%])")


def _flatten(args):
    flat = []
    for a in args:
        if isinstance(a, (str, bytes)):
            flat.append(a)
        elif isinstance(a, np.ndarray):
            flat.extend(a.flatten(order="F").tolist())
        elif isinstance(a, (list, tuple)):
            flat.extend(_flatten(a))
        else:
            flat.append(a)
    return flat


def _num(v):
    if isinstance(v, (bool, np.bool_)):
        return int(v)
    if isinstance(v, np.generic):
        return v.item()
    return v


def _convert(m, v):
    flags, width, prec, conv = m.group("flags"), m.group("width") or "", m.group("prec"), m.group("conv")
    p = f".{prec}" if prec is not None else ""
    if conv == "s":
        if not isinstance(v, str):
            v = _num(v)
            v = ("%d" % v) if float(v).is_integer() else ("%e" % v)  # MATLAB prints numbers given to %s compactly
        return f"%{flags}{width}{p}s" % v
    v = _num(v)
    if conv == "c":
        return f"%{flags}{width}c" % (v if isinstance(v, str) else chr(int(v)))
    x = float(v)
    if math.isnan(x) or math.isinf(x):
        txt = "NaN" if math.isnan(x) else ("Inf" if x > 0 else "-Inf")
        return f"%{flags.replace('0', '')}{width}s" % txt
    if conv in "diouxX":
        if x.is_integer():
            return f"%{flags}{width}{p}{conv}" % int(x)
        return f"%{flags}{width}e" % x  # MATLAB: %d of a non-integer uses %e
    return f"%{flags}{width}{p}{conv}" % x


def sprintf(fmt, *args):
    """Format like MATLAB's sprintf (escape sequences are taken from the Python string)."""
    specs = [m for m in _SPEC.finditer(fmt) if m.group("conv") != "%"]
    flat = _flatten(args)
    out, idx = [], 0
    while True:
        pos = 0
        for m in _SPEC.finditer(fmt):
            out.append(fmt[pos:m.start()])
            pos = m.end()
            if m.group("conv") == "%":
                out.append("%")
                continue
            if idx >= len(flat):
                if flat or idx > 0:
                    return "".join(out)
                pos = len(fmt)  # no arguments at all: MATLAB prints the literal text only
                break
            out.append(_convert(m, flat[idx]))
            idx += 1
        out.append(fmt[pos:])
        if not specs or idx >= len(flat):
            return "".join(out)
