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


PARMFLAGS = ["allkin", "alleq", "allreveq", "noncarbkin", "nonCaCO3kin", "CaCO3eq", "carbeq",
             "ignorefeedstock", "ignorecategory", "ignorenative", "ignoreeq"]
PARMLISTS = ["seteq", "setkin", "ignore", "implicit", "ignorefeedstock", "ignorecategory"]


def _as_list(v):
    if v is None:
        return []
    return list(v) if isinstance(v, (list, tuple)) else [v]


def _find_phase(t, name):
    """First row whose MATLABphase (else PHREEQCphase) matches ``name``, case-insensitively."""
    for col in ("MATLABphase", "PHREEQCphase"):
        hits = [i for i, v in enumerate(t[col]) if str(v).lower() == name.lower()]
        if hits:
            return hits[0]
    return None


def check_phase_type(t, ds=None, quiet=True, **opts):
    """checkphasetype.m.

    Fills empty PHREEQCphasetype values (``ignore`` if log10k_Eq is not
    finite, ``equilibrium`` with no kinetic mechanisms, else ``kinetic``),
    then applies the global flags (``allkin``, ``carbeq``, ...) and the
    per-phase lists (``seteq``, ``setkin``, ``ignore``, ``implicit``, ...).
    Options come from ``ds["parm"]`` when ``ds`` is given, overridden by
    keyword arguments.

    Returns the table, or ``ds`` (with ``ds["phasetable"]`` updated) when ``ds``
    is given. implicitcalc is called when the set of implicit phases changes.
    """
    parm = {k: [] for k in PARMLISTS}
    flags = {k: 0 for k in PARMFLAGS}
    if ds is not None:
        t = ds["phasetable"]
        for k in PARMLISTS:
            if k in ds.get("parm", {}):
                parm[k] = _as_list(ds["parm"][k])
        for k in PARMFLAGS:
            if k in ds.get("parm", {}):
                flags[k] = ds["parm"][k]
    for k, v in opts.items():
        if k in PARMLISTS:
            parm[k] = _as_list(v)
        elif k in PARMFLAGS:
            flags[k] = v
    t = t.copy()
    if "PHREEQCphasetype" not in t:
        t["PHREEQCphasetype"] = ""
    pt = t["PHREEQCphasetype"].astype(object).tolist()
    wimp0 = [p.lower() == "implicit" for p in pt]
    nwimp0 = sum(wimp0) if "impmolrat" in t else 0

    logk = pd.to_numeric(t["log10k_Eq"], errors="coerce").to_numpy()
    nkin = t["nkinmechs"].to_numpy(dtype=float)
    for i in range(len(pt)):
        if pt[i] == "" and not np.isfinite(logk[i]):
            pt[i] = "ignore"
    for i in range(len(pt)):
        if pt[i] == "" and nkin[i] <= 0:
            pt[i] = "equilibrium"
    for i in range(len(pt)):
        if pt[i] == "" and nkin[i] > 0:
            pt[i] = "kinetic"

    low = lambda v: str(v).lower()
    formula = t["formula"].tolist() if "formula" in t else [""] * len(pt)
    if flags["allkin"]:
        for i in range(len(pt)):
            if nkin[i] > 0 and low(pt[i]) == "equilibrium":
                pt[i] = "kinetic"
    elif flags["noncarbkin"] or flags["nonCaCO3kin"]:
        key = "CO3" if flags["noncarbkin"] else "CaCO3"
        bad = [key in str(f) for f in formula]
        for i in range(len(pt)):
            if not bad[i] and nkin[i] > 0 and low(pt[i]) == "equilibrium":
                pt[i] = "kinetic"
        for i in range(len(pt)):
            if bad[i]:
                pt[i] = "equilibrium"
    elif flags["carbeq"] or flags["CaCO3eq"]:
        key = "CO3" if flags["carbeq"] else "CaCO3"
        for i in range(len(pt)):
            if key in str(formula[i]) and low(pt[i]) == "kinetic":
                pt[i] = "equilibrium"
    elif flags["alleq"]:
        pt = ["equilibrium" if low(p) == "kinetic" else p for p in pt]
    elif flags["allreveq"]:
        rev = t["reversible"].to_numpy(dtype=float)
        for i in range(len(pt)):
            if rev[i] > 0 and low(pt[i]) == "kinetic":
                pt[i] = "equilibrium"

    def set_named(names, value):
        for name in names:
            i = _find_phase(t, name)
            if i is None:
                print(f"{name} not found in phase table")
            else:
                pt[i] = value

    set_named(parm["seteq"], "equilibrium")
    set_named(parm["setkin"], "kinetic")
    cats = t["category"].tolist() if "category" in t else [""] * len(pt)
    feeds = t["feedstock"].tolist() if "feedstock" in t else [""] * len(pt)
    for cat in parm["ignorecategory"]:
        if cat.lower() == "all":
            continue
        for i in range(len(pt)):
            if low(cats[i]) == cat.lower():
                pt[i] = "ignore"
    igfs = parm["ignorefeedstock"]
    if igfs and igfs[0].lower() == "all":
        igfs = sorted({feeds[i] for i in range(len(pt)) if low(cats[i]) == "feedstock"})
    for fs in igfs:
        for i in range(len(pt)):
            if low(feeds[i]) == fs.lower():
                pt[i] = "ignore"
    if flags["ignorenative"]:
        for i in range(len(pt)):
            if low(cats[i]) in ("soil", "native"):
                pt[i] = "ignore"
    if flags["ignoreeq"]:
        pt = ["ignore" if low(p) == "equilibrium" else p for p in pt]
    set_named(parm["ignore"], "ignore")
    set_named(parm["implicit"], "implicit")

    t["PHREEQCphasetype"] = pt
    wimp = [p.lower() == "implicit" for p in pt]
    needimpcalc = (nwimp0 != sum(wimp)) or any(a != b for a, b in zip(wimp0, wimp))
    if ds is None:
        if needimpcalc:
            print("WARNING: implicit phases have changed; need to run implicitcalc")
        return t
    ds["phasetable"] = t
    if needimpcalc:
        ds["parm"]["implicit"] = parm["implicit"]
        from .implicit import implicit_calc
        ds = implicit_calc(ds, quiet=quiet)
    return ds


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
