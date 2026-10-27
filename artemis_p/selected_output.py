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
