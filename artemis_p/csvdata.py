"""Read ARTEMIS "field"-format CSV files.

Python counterpart of readcsvdata.m. These files have a first column whose
header is ``field``. Rows whose first cell is ``data`` hold values; every other
row is metadata (``label``, ``units``, ``format``, ...) and is stored as a
sub-dict keyed by column name, e.g. ``s["units"]["molmass"]``.

MATLAB's readtable decides each column's type by content. Here a column is
numeric when every non-empty data value parses as a number, with empty cells
read as NaN; otherwise it is a list of strings with empty cells as ''.
"""

import csv
import re
from pathlib import Path

import numpy as np

# readcsvdata.m only lower-cases these metadata row names.
_LOWERED = {"field", "label", "units"}


def valid_name(name):
    """Approximate matlab.lang.makeValidName for column headers."""
    name = name.strip()
    name = re.sub(r"\s+(\w)", lambda m: m.group(1).upper(), name)
    name = re.sub(r"\W", "_", name)
    if not name or not name[0].isalpha():
        name = "x" + name
    return name


def _to_float(v):
    try:
        return float(v)
    except ValueError:
        return None


def _column(values):
    if all(v == "" or _to_float(v) is not None for v in values) and any(v != "" for v in values):
        return np.array([np.nan if v == "" else float(v) for v in values])
    return list(values)


def read_csv_data(path):
    """Return the columns and metadata of a field-format CSV file as a dict."""
    with open(Path(path), newline="", encoding="utf-8-sig") as f:
        rows = [r for r in csv.reader(f)]
    rows = [r for r in rows if any(c.strip() for c in r)]
    header = [valid_name(h) for h in rows[0]]
    if header[0].lower() != "field":
        raise ValueError(f"{path}: first column header is not 'field'")
    ncol = len(header)
    body = [(r + [""] * ncol)[:ncol] for r in rows[1:]]
    data = [r for r in body if r[0].strip().lower() == "data"]
    meta = [r for r in body if r[0].strip().lower() != "data"]

    s = {}
    for j, name in enumerate(header[1:], start=1):
        s[name] = _column([r[j].strip() for r in data])
    for r in meta:
        key = r[0].strip()
        if key.lower() in _LOWERED:
            key = key.lower()
        s[key] = {name: r[j].strip() for j, name in enumerate(header[1:], start=1)}
    return s
