"""Load the authors' saved MATLAB run structures as reference data for the port.

getphreeqcrun.m saves the run struct (named after the run, e.g. ``rBcomp``)
to ``<fname>.mat``. Its ``ds`` field is the full output of setupdatastruc.m,
so it gives intermediate values to test the Python preprocessor against.
"""

import numpy as np
import pandas as pd

try:  # mat-io decodes MATLAB tables, strings and datetimes; scipy cannot
    import matio
except ImportError:  # pragma: no cover
    matio = None
import scipy.io as sio
from scipy.io.matlab import mat_struct


def _convert_scipy(obj):
    if isinstance(obj, mat_struct):
        return {f: _convert_scipy(getattr(obj, f)) for f in obj._fieldnames}
    if isinstance(obj, np.ndarray) and obj.dtype == object:
        return [_convert_scipy(o) for o in obj.ravel()]
    return _leaf(obj)


def _leaf(obj):
    if isinstance(obj, np.ndarray):
        if obj.dtype.kind == "U":
            if obj.size == 0:
                return ""  # MATLAB empty char
            if obj.size == 1:
                return str(obj.ravel()[0])
            return [str(x) for x in obj.ravel()]
        if obj.dtype.kind in "iu":
            # MAT-files store integer-valued doubles in the smallest integer type
            obj = obj.astype(float)
        obj = np.squeeze(obj)
        if obj.ndim == 0 and obj.dtype.kind in "fb":
            return obj.item()
        return obj
    return obj


def _convert(obj):
    """Turn mat-io output into dicts, lists of strings and squeezed arrays."""
    if isinstance(obj, pd.DataFrame):
        out = {c: _convert_column(obj[c]) for c in obj.columns}
        out["_rownames"] = [str(i) for i in obj.index]
        return out
    if isinstance(obj, np.ndarray) and obj.dtype.names:
        if obj.size == 1:
            rec = obj.ravel()[0]
            return {f: _convert(rec[f]) for f in obj.dtype.names}
        return [{f: _convert(r[f]) for f in obj.dtype.names} for r in obj.ravel()]
    if isinstance(obj, np.ndarray) and obj.dtype == object:
        items = [_convert(o) for o in obj.ravel()]
        return items[0] if obj.size == 1 else items  # 1x1 cell unwrapped, as scipy's squeeze_me does
    return _leaf(obj)


def _convert_column(col):
    if col.dtype == object:
        return [_convert(v) for v in col]
    return col.to_numpy()


def load_run(path, name=None):
    """Return the run struct in ``path`` as nested dicts.

    ``name`` is the MATLAB variable name; by default the first non-internal
    one. MATLAB tables become dicts of columns (plus ``_rownames``).
    """
    if matio is not None:
        m = matio.load_from_mat(str(path))
        if name is None:
            name = next(k for k in m if not k.startswith("__"))
        return _convert(m[name])
    m = sio.loadmat(path, squeeze_me=True, struct_as_record=False)
    if name is None:
        name = next(k for k in m if not k.startswith("__"))
    return _convert_scipy(m[name])


def differences(py, ref, path="", rtol=1e-12, atol=0.0):
    """List the places where a Python structure differs from a MATLAB reference.

    Dicts are compared key by key (keys missing on either side are reported),
    numeric arrays with ``rtol``/``atol`` (NaNs equal), strings exactly.
    A MATLAB scalar may match a length-1 sequence.
    """
    out = []
    if isinstance(ref, dict):
        if not isinstance(py, dict):
            return [f"{path}: expected dict, got {type(py).__name__}"]
        for k in sorted(set(ref) | set(py)):
            if k not in py:
                out.append(f"{path}.{k}: missing in Python")
            elif k not in ref:
                out.append(f"{path}.{k}: not in MATLAB reference")
            else:
                out += differences(py[k], ref[k], f"{path}.{k}", rtol, atol)
        return out
    if isinstance(ref, str) or isinstance(py, str):
        return [] if py == ref else [f"{path}: {py!r} != {ref!r}"]
    if isinstance(ref, list) or isinstance(py, list):
        r = ref if isinstance(ref, list) else [ref]
        p = list(py) if isinstance(py, (list, np.ndarray)) else [py]
        if len(r) != len(p):
            return [f"{path}: length {len(p)} != {len(r)}"]
        for i, (a, b) in enumerate(zip(p, r)):
            out += differences(a, b, f"{path}[{i}]", rtol, atol)
        return out
    try:
        a = np.atleast_1d(np.asarray(py, dtype=float))
        b = np.atleast_1d(np.asarray(ref, dtype=float))
    except (TypeError, ValueError):
        return [f"{path}: cannot compare {type(py).__name__} with {type(ref).__name__}"]
    if a.shape != b.shape and a.size == b.size:
        a = a.reshape(b.shape)
    if a.shape != b.shape:
        return [f"{path}: shape {a.shape} != {b.shape}"]
    if not np.allclose(a, b, rtol=rtol, atol=atol, equal_nan=True):
        return [f"{path}: values differ (max abs {np.nanmax(np.abs(a - b)):.3g})"]
    return []
