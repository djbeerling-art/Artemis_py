"""Periodic table data. Python counterpart of getperiodictable.m."""

import os
from pathlib import Path

from .csvdata import read_csv_data


def default_file():
    return Path(os.environ["ARTEMIS"]) / "phreeqcdbincludefiles" / "periodictable.csv"


def get_periodic_table(path=None):
    """Return ``(pt, pt2)``.

    ``pt`` holds one array per column, indexed like ``pt["element"]``.
    ``pt2`` maps each element symbol to a dict of its values.
    """
    pt = read_csv_data(path or default_file())
    pt2 = {}
    for i, elt in enumerate(pt["element"]):
        pt2[elt] = {k: v[i] for k, v in pt.items() if k != "element" and not isinstance(v, dict)}
    return pt, pt2
