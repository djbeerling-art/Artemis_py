"""Read PHREEQC SELECTED_OUTPUT (.sel) files and compare them.

Python counterpart of the reading part of getselectedoutput.m.
"""

from pathlib import Path

import numpy as np
import pandas as pd


def read_sel(path):
    """Return a .sel file as a DataFrame.

    PHREEQC writes a tab-separated header line followed by rows of
    whitespace-padded values. Text columns (``state``) are kept as strings,
    everything else becomes float.
    """
    df = pd.read_csv(Path(path), sep="\t", skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    df = df.loc[:, [c for c in df.columns if c]]  # trailing tab gives an empty column
    for c in df.columns:
        if not pd.api.types.is_numeric_dtype(df[c]):
            df[c] = df[c].str.strip()
            converted = pd.to_numeric(df[c], errors="coerce")
            if converted.notna().all():
                df[c] = converted
    return df


def compare_sel(test, ref, rtol=1e-6, atol=1e-12):
    """Compare two .sel files column by column.

    Returns a DataFrame with one row per numeric column: the maximum absolute
    difference, the maximum relative difference, and the number of values
    outside ``rtol``/``atol``. Raises ValueError if the shapes or headers differ.
    """
    a, b = read_sel(test), read_sel(ref)
    if list(a.columns) != list(b.columns):
        raise ValueError("selected-output headers differ")
    if len(a) != len(b):
        raise ValueError(f"row counts differ: {len(a)} vs {len(b)}")
    rows = []
    for c in a.columns:
        if not (pd.api.types.is_numeric_dtype(a[c]) and pd.api.types.is_numeric_dtype(b[c])):
            if not (a[c] == b[c]).all():
                rows.append((c, np.nan, np.nan, int((a[c] != b[c]).sum())))
            continue
        x, y = a[c].to_numpy(float), b[c].to_numpy(float)
        d = np.abs(x - y)
        rel = d / np.maximum(np.abs(y), atol)
        nbad = int((d > atol + rtol * np.abs(y)).sum())
        rows.append((c, d.max(), rel.max(), nbad))
    return pd.DataFrame(rows, columns=["column", "max_abs", "max_rel", "n_outside"]).set_index("column")


def matlab_field(header):
    """getselectedoutput.m's struct field name for a .sel column heading:
    brackets removed, and a trailing charge cut off (``m_HCO3-`` -> ``m_HCO3``)
    or an inner +/- replaced by an underscore."""
    fld = header.strip().replace("(", "").replace(")", "")
    letters = [j for j, c in enumerate(fld) if c.isalpha()]
    for ch in "+-":
        if ch in fld:
            w = fld.index(ch)
            fld = fld.replace(ch, "_") if letters and max(letters) > w else fld[:w]
    return fld


def get_selected_output(path):
    """getselectedoutput.m: a .sel file as a dict of columns keyed by MATLAB field names.

    The leading text columns (``sim`` and ``state``, as MATLAB's importdata reads
    them) are lists of strings, the rest are float arrays. Repeated headings keep
    the last column, as MATLAB struct assignment does. ``statesim`` labels
    repeated state/sim pairs with a, b, c, ...
    """
    lines = Path(path).read_text().splitlines()
    heads = [h.strip() for h in lines[0].split("\t")]
    rows = [[v.strip() for v in l.split("\t")] for l in lines[1:] if l.strip()]
    ntext = 0
    for j in range(len(heads)):
        if heads[j] and any(_is_text(r[j]) for r in rows if j < len(r)):
            ntext = j + 1
    s = {"filename": str(path), "label": {}}
    for j, h in enumerate(heads):
        if not h:
            continue
        fld = matlab_field(h)
        s["label"][fld] = h.replace("_", " ")
        if j < ntext:
            s[fld] = [r[j] for r in rows]
        else:
            s[fld] = np.array([float(r[j]) for r in rows])
    if "state" in s and "sim" in s:
        statesim = [f"{a}_{b}" for a, b in zip(s["state"], s["sim"])]
        out = list(statesim)
        for k, ss in enumerate(statesim):
            w = [j for j, t in enumerate(statesim) if t == ss]
            if len(w) > 1:
                for n, j in enumerate(w):
                    out[j] = statesim[j] + chr(ord("a") + n)
        s["statesim"] = out
    return s


def _is_text(v):
    try:
        float(v)
        return False
    except ValueError:
        return True


def _checkfld(m, f, nsim):
    v = m[f]
    return isinstance(v, (dict, str, list)) or np.size(v) != nsim


def reshape_selected_output(m, out, iter, ileach=None):
    """reshapeselectedoutput.m: store one day's selected output in ``out``.

    Soil cells go to ``out['soil'][fld][cell-1, iter-1]``, solution 0 to
    ``out['influent']`` and the leachate to ``out['leachate']``. Daily changes
    (``dk_*``, ``d_*``) are summed over the day's steps; everything else takes
    the last step. On the first iteration all steps are kept in ``out['init']``.
    ``iter`` is the day index (getphreeqcrun passes ``'iter',i``).
    """
    soln = np.asarray(m["soln"])
    nsim = len(soln)
    usoln = np.unique(soln)
    flds = list(m)
    out.setdefault("label", m["label"])
    nrun = out["_nrun"]
    if "soil" in out and "pH" in out["soil"]:
        ileach = out["soil"]["pH"].shape[0] + 1
    if ileach is None:
        ileach = int(usoln.max())
    if not (soln == ileach).any():
        for f in flds:
            _put(out, "leachate", f, 0, iter, np.nan, 1, nrun)
    for sval in usoln:
        w = np.flatnonzero(soln == sval)
        if sval == 0:
            j, ostruc, nrow = 0, "influent", 1
        elif sval == ileach:
            j, ostruc, nrow = 0, "leachate", 1
        else:
            j, ostruc, nrow = int(sval) - 1, "soil", out["_ncells"]
        for f in flds:
            if iter == 1:
                init = out.setdefault("init", {}).setdefault(ostruc, {})
                if f == "sim":
                    init.setdefault("sim", {})[j] = np.array([float(m["sim"][x]) for x in w])
                if f == "state":
                    init.setdefault("state", {})[j] = [m["state"][x] for x in w]
                if "statesim" in m:
                    init.setdefault("statesim", {})[j] = [m["statesim"][x] for x in w]
                if not _checkfld(m, f, nsim):
                    init.setdefault(f, {})[j] = np.asarray(m[f])[w]
            if _checkfld(m, f, nsim):
                continue
            v = np.asarray(m[f])
            val = v[w].sum() if (f.startswith("dk_") or f.startswith("d_")) else v[w[-1]]
            _put(out, ostruc, f, j, iter, val, nrow, nrun)
    return out


def _put(out, ostruc, f, j, iter, val, nrow, nrun):
    d = out.setdefault(ostruc, {})
    if f not in d:
        d[f] = np.zeros((nrow, nrun))  # MATLAB grows arrays with zeros
    d[f][j, iter - 1] = val
