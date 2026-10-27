"""Daily run calendar and climate forcing time series.

Python counterparts of setruntimearrays (setupdatastruc.m) and of
readtimeseries, checktsdepths, settimeseries and testgetrundays
(readinputfiles.m).

Timeseries CSVs have a units line, a header line (``year,month,...``) and
monthly or daily rows. Depth-resolved files have ``nlay=N`` and one depth per
layer column in the units line; they are interpolated to the soil layers'
middle depths.
"""

import csv
import re
from datetime import date, datetime, timedelta
from pathlib import Path

import numpy as np


def _days_in_year(y):
    return 366 if (y % 4 == 0 and (y % 100 != 0 or y % 400 == 0)) else 365


def decyear(dt):
    """MATLAB decyear for a datetime: year + (day of year - 1 + fraction of day) / days in year."""
    start = datetime(dt.year, 1, 1)
    return dt.year + (dt - start).total_seconds() / 86400 / _days_in_year(dt.year)


def set_run_time_arrays(ds, t1, t2):
    """setruntimearrays: one entry per day from ``t1`` to ``t2`` inclusive (ISO date strings)."""
    datastep = ds.setdefault("datastep", 1)
    ds.setdefault("rollstep", datastep)
    d1, d2 = date.fromisoformat(t1), date.fromisoformat(t2)
    tt = [datetime(d1.year, d1.month, d1.day) + timedelta(days=i)
          for i in range(0, (d2 - d1).days + 1, datastep)]
    run = ds.setdefault("run", {})
    run["tt"] = np.array(tt, dtype="datetime64[ns]")
    ds["nrun"] = len(tt)
    run["year"] = np.array([t.year for t in tt], float)
    run["month"] = np.array([t.month for t in tt], float)
    run["day"] = np.array([t.day for t in tt], float)
    run["decyear"] = np.array([decyear(t) for t in tt])
    run["doy"] = np.array([t.timetuple().tm_yday for t in tt], float)
    run["exptday"] = np.arange(1, len(tt) + 1) / datastep
    return ds


def _interp1(x, y, xq):
    """MATLAB interp1 (linear): NaN outside the data range."""
    out = np.interp(xq, x, y)
    out[(xq < x[0]) | (xq > x[-1])] = np.nan
    return out


def _fill_ends(a, nrun):
    """Fill leading/trailing NaNs with the first/last good value, as settimeseries does."""
    good = np.flatnonzero(~np.isnan(a))
    if good.size:
        a[: good[0] + 1] = a[good[0]]
        a[good[-1]: nrun] = a[good[-1]]
    return a


def _interp2_rows(x, depth, array, xq, yq):
    """MATLAB interp2(x, depth, array, xq, yq) for row vector xq and column vector yq.

    Linear in time, then linear in depth; NaN outside either range.
    """
    in_time = np.vstack([_interp1(x, row, xq) for row in array])  # ndepth x nq
    out = np.full((len(yq), len(xq)), np.nan)
    for i, y in enumerate(yq):
        if y < depth[0] or y > depth[-1]:
            continue
        k = min(np.searchsorted(depth, y, side="right") - 1, len(depth) - 2)
        w = (y - depth[k]) / (depth[k + 1] - depth[k])
        out[i] = (1 - w) * in_time[k] + w * in_time[k + 1]
    return out


def _read(path):
    with open(Path(path), newline="", encoding="utf-8-sig") as f:
        rows = [r for r in csv.reader(f)]
    ucell = [c.strip() for c in rows[0]]
    names = [c.strip() for c in rows[1]]
    body = [r for r in rows[2:] if any(c.strip() for c in r)]
    cols = {}
    for j, n in enumerate(names):
        vals = [(r + [""] * len(names))[j].strip() for r in body]
        cols[n] = np.array([np.nan if v == "" else float(v) for v in vals])
    return rows[0], ucell, names, cols


def read_timeseries(ds, fname, keepmonthly=False):
    """readtimeseries for files with ``year,month[,day]`` columns."""
    if ds.get("parm", {}).get("usemeanforcings"):
        raise NotImplementedError("usemeanforcings is not ported yet")
    raw_units, ucell, names, cols = _read(fname)
    run = ds["run"]
    run_units = run.setdefault("units", {})
    if names[0].lower() != "year" or names[1].lower() != "month":
        raise NotImplementedError(f"{fname}: only year/month timeseries files are ported")
    for n in names[2:]:
        if not n.startswith("lay"):
            run_units[n] = raw_units[0]  # MATLAB splits without trimming

    # testgetrundays
    monthly = "day" not in cols
    day = np.full(len(cols["year"]), 15.0) if monthly else cols["day"]
    dates = [datetime(int(y), int(m), int(d), 12) for y, m, d in zip(cols["year"], cols["month"], day)]
    tdec = np.array([decyear(d) for d in dates])
    dexptday = 14 if monthly else 0
    exptday = np.array([(d - dates[0]).days + 1 + dexptday for d in dates], float)
    dtday = exptday[1] - exptday[0]
    is_monthly = 28 <= dtday <= 31

    # checktsdepths
    wnlay = [i for i, c in enumerate(ucell) if "nlay" in c]
    provided = ds.setdefault("provided", {})
    if wnlay:
        i0 = wnlay[0]
        nlay = int(re.match(r"nlay=(\d+)", ucell[i0]).group(1))
        # MATLAB tests numel(contains(ucell,'cm')), which is always true, so depths are always taken as cm
        depth = np.array([float(re.match(r"\s*([-\d.eE+]+)", ucell[i0 + 1 + i]).group(1)) for i in range(nlay)])
        array = np.vstack([cols[names[i0 + 1 + i]] for i in range(nlay)])
        stem = Path(fname).stem
        if "_timeseries_" in stem:
            stem = stem[stem.index("timeseries_"):]
        fld = stem.replace("timeseries", "").replace("_", "")
        run_units[fld] = ucell[0]
        fscl = 1
        if "month" in run_units[fld]:
            raise NotImplementedError("per-month depth-resolved units are not ported (MATLAB uses a stale index here)")
        val = _interp2_rows(tdec, depth, array, run["decyear"], ds["layer"]["middepth"]) / fscl
        for j in range(ds["nlayers"]):
            _fill_ends(val[j], ds["nrun"])
        run[fld] = val
        provided[fld] = str(fname)
    else:
        extra = {"juldate": None}  # testgetrundays adds juldate, which settimeseries then resamples
        for n in names + list(extra):
            if n in run or n in ("year", "month", "day"):
                continue
            if n == "juldate":
                y = np.array([_julian_date(d) for d in dates])
            else:
                y = cols[n]
            fscl = 1
            if is_monthly and not keepmonthly:
                run_units.setdefault(n, n)
                if "month" in run_units[n]:
                    fscl = 30
                    run_units[n] = run_units[n] + "/30"
            if dtday == 1:
                val = y.astype(float).copy()
            else:
                val = _interp1(tdec, y, run["decyear"]) / fscl
            run[n] = _fill_ends(val, ds["nrun"])
            provided[n] = str(fname)
    ds.setdefault("inputfiles", []).append(str(fname))
    return ds


def _julian_date(dt):
    """MATLAB juliandate for a datetime."""
    return (dt - datetime(2000, 1, 1, 12)).total_seconds() / 86400 + 2451545.0
