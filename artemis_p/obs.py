"""Observations: readobsdata.m, setruncellindices.m, assessrun2.m and checkuselog.m.

Observation files are CSVs with variable names on the first line and units on
the second. Each becomes a DataFrame in ``ds['obs'][<structure>][<name>]``,
with the units in ``df.attrs['units']`` (column -> unit) and the file name in
``df.attrs['description']``, standing in for MATLAB's table properties.
``assessrun2`` compares them with a finished run and stores the matched model
values and fit statistics in ``out['obs']``.
"""

import csv
import datetime as _dt
import re
from pathlib import Path

import numpy as np
import pandas as pd

__all__ = ["read_obs_data", "set_run_cell_indices", "assessrun2", "checkuselog"]


def _valid_name(h):
    """matlab.lang.makeValidName as readtable applies it to CSV headers."""
    h = h.strip()
    h = re.sub(r"\s+(\w)", lambda m: m.group(1).upper(), h)
    h = re.sub(r"[^A-Za-z0-9_]", "_", h)
    if not h or not h[0].isalpha():
        h = "x" + h
    return h


def _isnum(v):
    try:
        float(v)
        return True
    except ValueError:
        return False


def _parse_dates(vals, unit):
    dayfirst = unit.lower().startswith("dd/")
    out = []
    for v in vals:
        v = v.strip()
        if not v:
            out.append(pd.NaT)
        elif "/" in v:
            a, b, y = (int(x) for x in v.split("/"))
            d, m = (a, b) if dayfirst else (b, a)
            out.append(pd.Timestamp(y, m, d))
        else:
            out.append(pd.Timestamp(v))
    return pd.to_datetime(pd.Series(out))


def _decyear(d):
    ndays = np.where(d.dt.is_leap_year, 366, 365)
    return d.dt.year + (d.dt.day_of_year - 1) / ndays


def _read_table(fname):
    """readtable with detectImportOptions: names on line 1, units on line 2.

    A column is numeric when at least as many of its non-empty entries are
    numbers as are not (other entries become NaN); otherwise it is text.
    """
    with open(fname, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.reader(f))
    names = [_valid_name(h) for h in rows[0]]
    units = rows[1] + [""] * (len(names) - len(rows[1]))
    data = [r + [""] * (len(names) - len(r)) for r in rows[2:] if any(c.strip() for c in r)]
    cols, umap = {}, {}
    for j, (n, u) in enumerate(zip(names, units)):
        vals = [r[j] for r in data]
        u = u.strip()
        umap[n] = u
        if n == "datecollected" or "yyyy" in u.lower():
            cols[n] = _parse_dates(vals, u)
            continue
        filled = [v.strip() for v in vals if v.strip()]
        nnum = sum(_isnum(v) for v in filled)
        if nnum >= len(filled) - nnum:
            cols[n] = np.array([float(v) if _isnum(v) else np.nan for v in vals])
        else:
            cols[n] = np.array([v.strip() for v in vals], dtype=object)
    t = pd.DataFrame(cols)
    t.attrs["units"] = umap
    t.attrs["description"] = str(fname)
    return t


def _table_name(fname):
    f = Path(fname).stem
    tname = f.replace("-", "")  # MATLAB discards its strrep(f,'obs_','') here
    if "obs_rain_" in f or "obs_influent" in f or "obs_precip" in f:
        ostruc = "influent"
    elif any(s in f for s in ("obs_stream", "obs_leach", "obs_drain", "obs_effluent")):
        ostruc = "leachate"
    elif any(s in f for s in ("obs_pore", "obs_lys", "obs_rhizon", "obs_ex", "obs_avail")):
        ostruc = "soil"
    elif any(s in f for s in ("obs_bio", "obs_veg", "obs_crop")):
        ostruc = "plant"
    elif "total" in f:
        ostruc = "total"
    elif "soil" in f:
        ostruc = "soil"
    else:
        raise ValueError(f"(readobsdata) cannot tell which structure {fname} belongs to")
    return ostruc, tname


def _run_index(t, run):
    year, month, day = (np.ravel(run[k]) for k in ("year", "month", "day"))
    idx = np.zeros(len(t))
    for i in range(len(t)):
        w = np.flatnonzero((year == t["year"].iat[i]) & (month == t["month"].iat[i]) & (day == t["day"].iat[i]))
        if w.size:
            idx[i] = w[0] + 1
    return idx


def _middepth(t, ds, nan_bottom):
    names = list(t.columns)
    if "topdepth" in names and "bottomdepth" in names:
        top, bot = t["topdepth"].to_numpy(float), t["bottomdepth"].to_numpy(float)
        mid = (bot - top) / 2 + top
        if nan_bottom == "mid":  # setruncellindices/assessrun2: missing bottom -> top + 1
            w = np.isnan(bot)
            mid[w] = top[w] + 1
        else:  # readobsdata sets bottom = top + 1 first, so the mid-point is top + 0.5
            w = np.isnan(bot) & ~np.isnan(top)
            mid[w] = top[w] + 0.5
        return mid
    if "depth" in names:
        return t["depth"].to_numpy(float)
    if "middepth" in names:
        return t["middepth"].to_numpy(float)
    if any("layer" in n for n in names):
        lay = t[[n for n in names if "layer" in n][0]].to_numpy(int) - 1
        return (np.asarray(ds["layer"]["topdepth"])[lay] + np.asarray(ds["layer"]["bottomdepth"])[lay]) / 2
    return None


def _cell_index(t, ds, none_value=0, nan_bottom="mid"):
    """Cell for each row from its mid-point depth (depth units other than cm are not
    rescaled: MATLAB computes the scale but never applies it)."""
    if not any("depth" in n for n in t.columns):
        return np.full(len(t), float(none_value))
    idx = np.zeros(len(t))
    mid = _middepth(t, ds, nan_bottom)
    if mid is None:
        return idx
    top = np.asarray(ds["cell"]["topdepth"], float).reshape(len(ds["cell"]["topdepth"]), -1)[:, 0]
    bot = np.asarray(ds["cell"]["bottomdepth"], float).reshape(len(ds["cell"]["bottomdepth"]), -1)[:, 0]
    for j in range(int(ds["ncells"])):
        idx[(mid >= top[j]) & (mid < bot[j])] = j + 1
    return idx


def read_obs_data(ds, fname, name=None):
    """readobsdata.m: add one observation file to ``ds['obs']``.

    ``name`` is the file name recorded in ``ds['inputfiles']`` (default ``fname``).
    """
    if "obsID" in ds.get("parm", {}):
        raise NotImplementedError("getscreenedtable (parm.obsID) is not ported yet")
    ostruc, tname = _table_name(fname)
    t = _read_table(fname)
    if "datecollected" in t:
        d = t["datecollected"]
        t["year"], t["month"], t["day"] = d.dt.year.astype(float), d.dt.month.astype(float), d.dt.day.astype(float)
        t["doy"] = d.dt.day_of_year.astype(float)
        t["decyear"] = _decyear(d)
    if "run" in ds:
        t["runindex"] = _run_index(t, ds["run"])
        if "exptday" not in t:
            t["exptday"] = t["runindex"]
        if "ncells" in ds:
            t["cellindex"] = _cell_index(t, ds, nan_bottom="readobsdata")
    if "inputfiles" in ds:
        ds["inputfiles"].append(name if name is not None else str(fname))
    ds.setdefault("obs", {}).setdefault(ostruc, {})[tname] = t
    return ds


def set_run_cell_indices(ds):
    """setruncellindices(ds,'obs'): run day and cell for every observation row."""
    if "obs" not in ds:
        print("(setruncellindices) Data structure has no observations")
        return ds
    for tabs in ds["obs"].values():
        for t in tabs.values():
            if "year" in ds.get("run", {}):
                t["runindex"] = _run_index(t, ds["run"])
                if "exptday" not in t:
                    t["exptday"] = t["runindex"]
            t["cellindex"] = _cell_index(t, ds)
    return ds


def checkuselog(fld, *opts):
    """checkuselog.m: 1 if the field should be compared (or plotted) as log10."""
    nologkinetics = any(o.lower() == "nologkinetics" for o in opts)
    low = fld.lower()
    if fld == "pH" or low in ("alk", "mass_h2o", "mmh2o", "charge", "pe") or fld == "pct_err":
        return 0
    if any(s in fld for s in ("release", "BSpct", "eff", "exch")):
        return 0
    uselog = 1
    if len(fld) > 2:
        if low[:3] == "si_" or low[:2] == "d_":
            uselog = 0
        if low[:2] == "k_" or low[:3] == "dk_":
            if nologkinetics:
                uselog = 0
            elif not any(s in fld for s in ("uptake", "itrification", "Ureolysis", "Nfix")):
                uselog = 0
    return uselog


def _convert_units(t, ds):
    """assessrun2's convertunits: element columns to mol/L (or mol/cell for mg/kg)."""
    pt2, solute = ds["pt2"], ds["solute"]
    units = dict(t.attrs.get("units", {}))
    for f in list(t.columns):
        u = units.get(f, "")
        if f in pt2:
            xsp = psp = f
            molmass = pt2[xsp]["molmass"]
        elif f == "SO4":
            xsp = psp = "S"
            molmass = pt2["S"]["molmass"] + 4 * pt2["O"]["molmass"]
        elif f == "NO3":
            xsp = psp = "N"
            molmass = pt2["N"]["molmass"] + 3 * pt2["O"]["molmass"]
        elif f == "NH4" or "Amm" in f:
            xsp, psp = "N", "Amm"
            molmass = pt2["N"]["molmass"] + 4 * pt2["H"]["molmass"]
        elif "DIC" in f or "TIC" in f:
            xsp, psp = "C", "DIC"
            molmass = pt2["C"]["molmass"]
        elif "TOC" in f or "DOC" in f:
            xsp, psp = "C", "DOC"
            molmass = pt2["C"]["molmass"]
        elif f in solute:
            xsp = psp = f
            molmass = pt2[xsp]["molmass"]
        elif "sys" in f:
            xsp, psp = f.replace("sys", ""), f
            molmass = pt2[xsp]["molmass"]
            if xsp not in solute:
                continue
        else:
            continue
        col = t[f]
        if not pd.api.types.is_numeric_dtype(col):  # e.g. 'bdl' entries: str2double gives NaN
            col = pd.to_numeric(col, errors="coerce").astype(float)
            t[f] = col
        t[psp] = col
        if u.lower() in ("mg/l", "ppm"):
            t[psp] = col / molmass / 1000
            units[f] = "mol/L"
        elif "mg/kg" in u:
            cell = t["cellindex"].to_numpy(int)
            if (cell < 1).any():
                raise IndexError(f"(assessrun2/convertunits) {f}: rows without a soil cell")
            gsoil = np.asarray(ds["cell"]["gsoil"], float).reshape(int(ds["ncells"]), -1)[:, 0]
            t[psp] = col / molmass * gsoil[cell - 1] / 1e6
            units[f] = "mol/cell"
        elif "mmol" in u or "mM" in u:
            t[psp] = col / 1000
            units[f] = "mol/L"
        elif "umol" in u or "uM" in u:
            t[psp] = col / 1e6
            units[f] = "mol/L"
        elif "mol" in u or "{\\mu}" in u:
            t[psp] = col / 1e6
            units[f] = "mol/L"
        elif not u:
            print(f"(assessrun2/convertunits) Assuming units of mol/L for {f} in {t.attrs.get('description', '')}")
    t.attrs["units"] = units
    return t


def assessrun2(out, ds):
    """assessrun2.m: model values on observation days and cells, and fit statistics."""
    if "obs" not in ds:
        print("(assessrun2) No observed data")
        return out
    out.pop("obs", None)
    md = out["metadata"]
    runlim1 = int(md.get("istart0", md["istart"]))
    runlim2 = int(md["ilast"])
    ncells = int(ds["ncells"])
    for ostruc, tabs in ds["obs"].items():
        if ostruc == "units" or not isinstance(tabs, dict):
            continue
        modstruct = {"porewater": "soil", "exchanger": "soil"}.get(ostruc.lower(), ostruc)
        mod = out[modstruct]
        if "runindex" not in mod and "runday" in mod:
            mod["runindex"] = np.array(mod["runday"], copy=True)
        for tname, tobs in tabs.items():
            if not isinstance(tobs, pd.DataFrame):
                continue
            t = tobs.copy()
            t.attrs = {k: (dict(v) if isinstance(v, dict) else v) for k, v in tobs.attrs.items()}
            t["cellindex"] = _cell_index(t, ds, none_value=1)
            t = _convert_units(t, ds)
            ri = t["runindex"].to_numpy(float)
            good = (ri >= runlim1) & (ri <= runlim2)
            if modstruct == "soil":
                ci = t["cellindex"].to_numpy(float)
                good &= (ci > 0) & (ci <= ncells)
            if not good.any():
                continue
            attrs = t.attrs
            t = t[good].reset_index(drop=True)
            t.attrs = attrs
            ob = out.setdefault("obs", {}).setdefault(modstruct, {})
            ob.setdefault("yobs", {})[tname] = t
            for f in list(t.columns):
                if _checkrunfld(out, ds, modstruct, f, t.attrs["units"].get(f, "")):
                    _rmseobs(out, t, f, modstruct, tname)
    return out


def _checkrunfld(out, ds, modstruct, fld, fldunits):
    mod = out[modstruct]
    if fld in mod:
        return True
    if fld.lower() == "camgrat":
        mod["CaMgrat"] = mod["Ca"] / mod["Mg"]
        out["label"]["CaMgrat"], out["units"]["CaMgrat"] = "Ca/Mg", "molar ratio"
        return True
    if fld.lower() == "cax2mgx2rat" and "m_CaX2" in mod:
        mod["CaX2MgX2rat"] = mod["m_CaX2"] / mod["m_MgX2"]
        out["label"]["CaX2MgX2rat"], out["units"]["CaX2MgX2rat"] = "CaX2/MgX2", "molar ratio"
        return True
    if "_exch" in fld:
        elt = fld.replace("_exch", "")
        if elt not in ds["pt2"]:
            if elt.lower() not in ("nh4", "ammh"):
                return False
            mfld, molmass = "AmmHX", ds["pt2"]["N"]["molmass"]
        else:
            sp = ds["pt2"][elt]["exchangespecies"]
            mfld, molmass = "m_" + (sp[0] if isinstance(sp, (list, tuple)) else sp), ds["pt2"][elt]["molmass"]
        if mfld not in mod:
            return False
        if fldunits.lower() == "ppm":
            gsoil = np.asarray(ds["cell"]["gsoil"], float).reshape(int(ds["ncells"]), -1)[:, 0]
            if modstruct == "total":
                gsoil = np.sum(gsoil)
            else:
                gsoil = gsoil[:, None]
            mod[fld] = mod[mfld] * mod["mass_H2O"] * molmass * 1e6 / gsoil
            out["units"][fld], out["label"][fld] = fldunits, "Exchangeable " + elt
            return True
    return False


def _rmseobs(out, t, fld, modstruct, tname):
    if fld.lower() in ("year", "month", "day", "doy", "decyear", "exptday", "runindex") or "ind" in fld:
        return
    mod = out[modstruct]
    if fld not in mod:
        return
    uselog = checkuselog(fld, "nologkinetics", "nologuptake")
    wcell = t["cellindex"].to_numpy(int).copy()
    wrun = t["runindex"].to_numpy(int)
    wcell[wcell == 0] = 1
    ob = out["obs"][modstruct]
    ym = ob.setdefault("ymod", {})
    if tname not in ym:
        ym[tname] = pd.DataFrame({"runindex": wrun.astype(float), "cellindex": wcell.astype(float)})
    arr = np.asarray(mod[fld], float)
    if arr.ndim == 1:  # total fields are rows in MATLAB
        arr = arr[None, :]
    ymod = arr[wcell - 1, wrun - 1]
    ym[tname][fld] = ymod
    yobs = t[fld].to_numpy(float).copy()
    yobs[yobs < 0] = np.nan
    with np.errstate(divide="ignore", invalid="ignore"):
        if uselog:
            yobs, ymod = np.log10(yobs), np.log10(ymod)
        diff = ymod - yobs
        w = np.isfinite(diff)
        m = {"RMSE": np.sqrt(np.mean(diff[w] ** 2)) if w.any() else np.nan, "n": float(w.sum())}
        mo, mm = (np.mean(yobs[w]), np.mean(ymod[w])) if w.any() else (np.nan, np.nan)
        m.update(meanobs=mo, meanmod=mm, uselog=float(uselog), bias=mm - mo)
        n = w.sum()
        m["crmsd"] = np.sqrt(np.sum(((ymod[w] - mm) - (yobs[w] - mo)) ** 2) / n) if n else np.nan
        # MATLAB: sqrt of the sum, then divided by n (not the usual RMS difference)
        m["rmsd"] = np.sqrt(np.sum((ymod[w] - yobs[w]) ** 2)) / n if n else np.nan
        m["ccoef"] = np.corrcoef(ymod[w], yobs[w])[0, 1] if n > 1 else np.nan
        m["stdevmod"] = np.std(ymod[w]) if n else np.nan
        m["stdevobs"] = np.std(yobs[w]) if n else np.nan
    tass = pd.DataFrame({fld: list(m.values())}, index=list(m))
    asg = ob.setdefault("assess", {})
    asg[tname] = tass if tname not in asg else pd.concat([asg[tname], tass], axis=1)
