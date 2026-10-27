"""Mineral phase tables. Python counterparts of getphasetable.m and checkphasetype.m.

A phase table is a pandas DataFrame with one row per phase, indexed by
``MATLABphase`` (the PHREEQC phase name with non-word characters removed).
Column names are normalised as getphasetable.m does: ``Phase`` becomes
``PHREEQCphase``, ``kEq`` becomes ``log10k_Eq``, ``k_Neutral`` becomes
``log10k_Neutral``, ``E_Neutral`` becomes ``Eapp_Neutral``, and so on. The
units row (row 2 of each CSV) is kept in ``df.attrs["units"]``.
"""

import csv
import re
from pathlib import Path

import numpy as np
import pandas as pd

# Canonical column name -> alternative header patterns (anchored at the start,
# case-insensitive), in the order getphasetable.m tries them.
VPHASE = {
    "PHREEQCphase": ["PHREEQCphase", "phase", "mineral"],
    "PHREEQCphasetype": ["PHREEQCphasetype", "phaseflag"],
    "formula": ["formula"],
    "reaction": ["reaction"],
    "InitialMoles": [r"Initial[_\s]*Moles"],
    "frac": ["frac", "mfrac", r"^Mass[_\s]*Fraction", "wgtpct", "weight"],
    "G": ["G", "Density"],
    "reversible": ["reversible", "precipallowed"],
    "log10k_Eq": [r"log10k[_\s]*Eq", r"k[_\s]*Eq", r"log[_\s]*k_eq", r"log[_\s]*k"],
    "deltaH": ["deltaH", "delta_H", "delta H"],
}
VSCALE = {"frac": [1, 1, 1, 0.01, 0.01], "G": [1, 1000]}
VTESTUNITS = {"G": [r"g[/\s]cm[^{-]*3}", r"kg[/\s]m.*3"]}
VGOODUNITS = {"G": "g cm^{-3}"}
VMECH = {
    "log10k": [r"log10[_\s]*k[_\s]", r"log[_\s]*k[_\s]*", r"k[_\s]*"],
    "Eapp": [r"Eapp[_\s]*", r"Ea[_\s]*", r"E[_\s]*"],
    "n": [r"n[_\s]*", r"order[_\s]*", r"rorder[_\s]*"],
    "q": [r"q[_\s]*"],
    "p": [r"p[_\s]*"],
}


def _sub(pattern, repl, names):
    rx = re.compile(pattern, re.IGNORECASE)
    return [rx.sub(repl, n) for n in names]


def _read_csv(path):
    """readtable with VariableUnitsLine=2 and preserved names: header, units, typed columns."""
    with open(Path(path), newline="", encoding="utf-8-sig") as f:
        rows = [r for r in csv.reader(f)]
    names = [h.strip() for h in rows[0]]
    n = len(names)
    units = [u.strip() for u in (rows[1] + [""] * n)[:n]]
    body = [(r + [""] * n)[:n] for r in rows[2:] if any(c.strip() for c in r)]
    cols = {}
    for j, name in enumerate(names):
        vals = [r[j].strip() for r in body]
        try:
            cols[j] = np.array([np.nan if v == "" else float(v) for v in vals])
        except ValueError:
            cols[j] = vals
    return names, units, cols


def _check_varnames(names, units, cols):
    """checkvarnames in getphasetable.m: rename columns to canonical names, rescale."""
    tvar = list(names)
    for canon, alts in VPHASE.items():
        if sum(v == canon for v in tvar) == 1:
            continue
        found = False
        for ialt, alt in enumerate(alts):
            if found:
                continue
            tvar = _sub("^" + alt, canon, tvar)
            wu = [j for j, v in enumerate(tvar) if v.lower() == canon.lower()]
            if not wu:
                continue
            wu = wu[0]
            if canon in VTESTUNITS:
                if re.match("^" + VTESTUNITS[canon][ialt], units[wu], re.IGNORECASE):
                    cols[wu] = cols[wu] * VSCALE[canon][ialt]
                    units[wu] = VGOODUNITS[canon]
                    found = True
            elif canon in VSCALE:
                # getphasetable.m compares the column name with the fields of
                # vscale, so this is true whenever the column is named frac or G.
                if tvar[wu].lower() in (k.lower() for k in VSCALE):
                    cols[wu] = cols[wu] * VSCALE[canon][ialt]
                    found = True
    good = "log10k_"
    for pat in VMECH["log10k"]:
        tvar = _sub("^" + pat + r"[_\s]*", good, tvar)
    mechs = [v[len(good):] for v in tvar if good in v and "log10k_Eq" not in v]
    for mv, pats in VMECH.items():
        if mv == "log10k":
            continue
        for mech in mechs:
            for pat in pats:
                tvar = _sub(pat + r"[_\s]*" + re.escape(mech), f"{mv}_{mech}", tvar)
    return tvar


def check_mech(t):
    """checkmech in getphasetable.m: count rate mechanisms per phase."""
    t["nkinmechs"] = 0.0
    t["nprecipmechs"] = 0.0
    mechs = [c[len("log10k_"):] for c in t.columns if "log10k_" in c and "_Eq" not in c]
    for mech in mechs:
        k = t["log10k_" + mech].astype(float)
        active = k.notna() & (k != 0)
        t.loc[active, "nkinmechs"] += 1
        if "precip" in mech:
            t.loc[active, "nprecipmechs"] += 1
    t.attrs["RATESmechanisms"] = mechs
    return t


def check_phase_type(t):
    """checkphasetype.m called with a table and no options (as getphasetable.m does).

    Phases with no PHREEQCphasetype become ``ignore`` if log10k_Eq is not
    finite, ``equilibrium`` if they have no kinetic mechanisms, otherwise
    ``kinetic``.
    """
    if "PHREEQCphasetype" not in t:
        t["PHREEQCphasetype"] = ""
    pt = t["PHREEQCphasetype"].astype(object).copy()
    empty = pt == ""
    logk = pd.to_numeric(t["log10k_Eq"], errors="coerce")
    pt[empty & ~np.isfinite(logk)] = "ignore"
    empty = pt == ""
    pt[empty & (t["nkinmechs"] <= 0)] = "equilibrium"
    empty = pt == ""
    pt[empty & (t["nkinmechs"] > 0)] = "kinetic"
    t["PHREEQCphasetype"] = pt
    return t


def _pad_value(col):
    return np.nan if pd.api.types.is_numeric_dtype(col) else ""


def get_phase_table(csvfile, table=None, feedstock=None):
    """Read a phases CSV and append it to ``table`` (if given), as getphasetable.m does."""
    names, units, cols = _read_csv(csvfile)
    tvar = _check_varnames(names, units, cols)
    tj = pd.DataFrame({tvar[j]: cols[j] for j in range(len(tvar))})
    unit_map = dict(zip(tvar, units))
    if "PHREEQCphase" not in tj:
        if "Phase" in tj:
            tj["PHREEQCphase"] = tj["Phase"]
        else:
            print(f"(getphasetable) No PHREEQCphase or Phase variable found in {csvfile} (skipping)")
            return table
    tj = tj[tj["PHREEQCphase"] != ""].reset_index(drop=True)
    if feedstock is not None:
        tj["feedstock"] = feedstock
    elif "feedstock" not in tj:
        tj["feedstock"] = ""
    elif pd.api.types.is_numeric_dtype(tj["feedstock"]):
        # an all-empty feedstock column reads as NaN; it is text in MATLAB
        tj["feedstock"] = ""
    tj["csvfile"] = str(csvfile)
    tj["MATLABphase"] = [re.sub(r"\W", "", p) for p in tj["PHREEQCphase"]]
    if table is not None:
        for c in table.columns:
            if c not in tj:
                tj[c] = _pad_value(table[c])
        for c in tj.columns:
            if c not in table:
                table[c] = _pad_value(tj[c])
        units_all = {**table.attrs.get("units", {}), **unit_map}
        table = pd.concat([table.reset_index(drop=True), tj[table.columns]], ignore_index=True)
        unit_map = units_all
    else:
        table = tj
    table = check_mech(table)
    table = check_phase_type(table)
    table.index = table["MATLABphase"].tolist()
    table.attrs["units"] = unit_map
    return table
