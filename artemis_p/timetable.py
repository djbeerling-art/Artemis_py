"""Management timetable (crops, tillage, fertiliser and feedstock treatments).

Python counterparts of readtimetable.m, getscreenedtable.m and getdecyear.m.
"""

import csv
from pathlib import Path

import numpy as np

from .csvdata import valid_name
from .fixes import enabled as fixed

# Cumulative day of year at the start of each month (non-leap year), indexed by month 1..12
_CMON = np.cumsum([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])


def is_leap(year):
    """Gregorian leap years."""
    year = np.asarray(year, dtype=float)
    return (year % 4 == 0) & ((year % 100 != 0) | (year % 400 == 0))


def is_leap_matlab(year):
    """Leap-year rule exactly as written in getdecyear.m.

    getdecyear.m uses ``year/4 == floor(year/4) & year/400 < ceil(year/400)``:
    divisible by 4 and *not* divisible by 400. That is not the Gregorian rule
    (2000 is a leap year, 1900 and 2100 are not), but it gives the same answer
    for every year from 1901 to 2099 except 2000.
    """
    year = np.asarray(year, dtype=float)
    return (year / 4 == np.floor(year / 4)) & (year / 400 < np.ceil(year / 400))


def get_decyear(t, nomday=15, gregorian=False):
    """Add ``leapy``, ``leapd``, ``doy``, ``decyear`` and labels to a dict of columns."""
    year = np.asarray(t["year"], dtype=float)
    month = np.asarray(t["month"], dtype=float)
    leapy = is_leap(year) if gregorian else is_leap_matlab(year)
    t["leapy"] = leapy
    t["leapd"] = leapy & (month > 2)
    if "day" not in t:
        t["day"] = month * 0 + nomday
    day = np.asarray(t["day"], dtype=float)
    day = np.where(np.isnan(day), nomday, day)
    t["day"] = day
    t["doy"] = _CMON[month.astype(int) - 1] + day + t["leapd"]
    t["decyear"] = year + (t["doy"] - 0.5) / (365 + leapy)
    t.setdefault("label", {})
    t["label"]["decyear"] = "Decimal year"
    t["label"]["doy"] = "Day of year"
    return t


def read_table(path):
    """Read a plain CSV the way MATLAB's readtable types columns.

    Columns where every non-empty value is numeric become float arrays (empty
    cells NaN); other columns become lists of strings (empty cells '').
    """
    with open(Path(path), newline="", encoding="utf-8-sig") as f:
        rows = [r for r in csv.reader(f)]
    header = [valid_name(h) for h in rows[0]]
    body = [(r + [""] * len(header))[:len(header)] for r in rows[1:] if any(c.strip() for c in r)]
    t = {}
    for j, name in enumerate(header):
        vals = [r[j].strip() for r in body]
        try:
            t[name] = np.array([np.nan if v == "" else float(v) for v in vals])
        except ValueError:
            t[name] = vals
    return t


def screen_table(t, ds, idfld="expID", keepnull=False):
    """getscreenedtable.m: keep rows whose ``idfld`` matches ``ds[idfld]``.

    Returns ``t`` unchanged when the table or ``ds`` has no such field.
    """
    if idfld not in t or idfld not in ds:
        return t
    target = str(ds[idfld]).lower()
    ids = [str(v).lower() for v in t[idfld]]
    keep = [i for i, v in enumerate(ids) if v == target or (keepnull and v == "")]
    out = {}
    for k, v in t.items():
        out[k] = v[keep] if isinstance(v, np.ndarray) else [v[i] for i in keep]
    return out


def read_timetable(ds, fname, baskg=None):
    """readtimetable.m, including the long-run options ``repeattimetable``,
    ``ceasefeedstockyear`` and ``nofeedstockyears`` from ``ds['parm']``."""
    parm = ds.get("parm", {})
    t = screen_table(read_table(fname), ds)
    if "inputfiles" in ds:
        ds["inputfiles"].append(str(fname))
    gregorian = fixed(ds, "leap_year")
    tt = ds["timetable"] = get_decyear(dict(t), gregorian=gregorian)
    if baskg is not None and "baskg" in tt:
        b = tt["baskg"]
        b[b > 0] = baskg
    if "repeattimetable" in parm:
        tt = ds["timetable"] = _repeat_timetable(tt, parm, gregorian)
    if "ceasefeedstockyear" in parm:
        fs = (tt["decyear"] > parm["ceasefeedstockyear"]) & _is_feedstock(tt)
        if fs.any():
            print(f"Reducing feedstock dosages to zero for {fs.sum()} treatments starting from "
                  f"{tt['decyear'][fs][0]:g}")
            _no_feedstock(tt, fs)
    elif "nofeedstockyears" in parm:
        bad = np.isin(tt["year"], np.atleast_1d(parm["nofeedstockyears"])) & _is_feedstock(tt)
        print(f"(readtimetable) Reduced feedstock dosages to zero for {bad.sum()} entries")
        _no_feedstock(tt, bad)
    if "nrun" not in ds:
        raise NotImplementedError("readtimetable without ds['nrun'] (julian dates) is not ported yet")
    return ds


def _is_feedstock(tt):
    return np.array([str(f).lower() == "feedstock" for f in tt["doseform"]])


def _no_feedstock(tt, w):
    for i in np.flatnonzero(w):
        tt["doseform"][i] = ""
    tt["dose"][w] = 0


def _repeat_timetable(tt, parm, gregorian):
    """Extend the timetable by repeating the entries from ``repeattimetable[0]`` (decimal
    year; -1 for the last calendar year) until ``repeattimetable[1]``."""
    rep = [float(v) for v in parm["repeattimetable"]]
    last = tt["decyear"][-1]
    if rep[0] < 0:
        rep[0] = parm["repeattimetable"][0] = tt["year"][-1]
    wrep = np.flatnonzero(tt["decyear"] >= rep[0])
    if not len(wrep):
        print(f"Unable to extend timetable to {rep[1]:g} because no entries found after {rep[0]:g}")
        return tt
    nrepyears = len(np.unique(tt["year"][wrep]))
    nyr = int(np.ceil(rep[1] - last))
    # MATLAB round: halves away from zero
    ncycles = int(np.floor(nyr / nrepyears + 0.5))
    print(f"Last timetable entry: {last:g} ({len(tt['decyear'])}) Adding {nyr} years ({ncycles:g} rep cycles)")
    out = {}
    for k, v in tt.items():
        if isinstance(v, dict):
            out[k] = v
        elif isinstance(v, np.ndarray):
            out[k] = np.concatenate([v] + [v[wrep]] * ncycles)
        else:
            out[k] = list(v) + [v[i] for i in wrep] * ncycles
    n = len(tt["decyear"])
    year = np.asarray(out["year"], float)
    for j in range(1, ncycles + 1):
        year[n + (j - 1) * len(wrep):n + j * len(wrep)] += j * nrepyears
    out["year"] = year.astype(np.asarray(tt["year"]).dtype)
    # the decimal years are recomputed, not copied, so that leap years count
    out = get_decyear(out, gregorian=gregorian)
    print(f"Extended timetable to {out['decyear'][-1]:g} (now {len(out['decyear'])} entries) using {len(wrep)} "
          f"entries from {tt['decyear'][wrep[0]]:g} to {tt['decyear'][wrep[-1]]:g}")
    return out
