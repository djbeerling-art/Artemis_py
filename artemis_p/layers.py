"""Initial soil layer data. Python counterpart of readinitlayerdata.m.

An initlayers CSV has one row per parameter (``parameter``, ``label``,
``units`` columns) and one column per soil layer (any header containing
"lay"). Rows ending in ``_soln``/``_hotwater`` (and ``pH``) go to
``layer["soilsoln"]``; rows ending in ``_exch``/``_basesat``/``_avail`` go to
``layer["exchangeable"]``/``["basesat"]``/``["available"]``; everything else is
stored directly as ``layer[parameter]``.
"""

import csv
from pathlib import Path

import numpy as np

from .csvdata import valid_name
from .periodictable import get_periodic_table
from .species import set_model_species


def fill_nan(a):
    """fillNaN in readinitlayerdata.m: replace NaNs with the nearest non-NaN layer value.

    With one good value, every NaN takes it. Otherwise each NaN takes the
    closest good layer (the shallower one on a tie).

    MATLAB note: the original loops ``for inan=wnan`` over a column vector,
    which runs once with the whole vector. It only works as intended when
    there is exactly one NaN; with more it errors or fills wrongly. This
    version fills each NaN separately.
    """
    a = np.array(a, dtype=float)
    bad = np.flatnonzero(np.isnan(a))
    good = np.flatnonzero(~np.isnan(a))
    if bad.size == 0 or good.size == 0:
        return a
    out = a.copy()
    for i in bad:
        d = np.abs(good - i)
        out[i] = a[good[np.argmin(d)]]
    return out


def _numeric(vals):
    try:
        return np.array([np.nan if v == "" else float(v) for v in vals])
    except ValueError:
        return None


def read_init_layer_data(ds, fname):
    """Read an initlayers CSV into ``ds["layer"]`` (and set ``ds["nlayers"]``)."""
    parm = ds.setdefault("parm", {})
    if "soil_depth" in parm:
        raise NotImplementedError("readinitlayerdata soil_depth truncation is not ported yet")
    with open(Path(fname), newline="", encoding="utf-8-sig") as f:
        rows = [r for r in csv.reader(f)]
    header = [valid_name(h) for h in rows[0]]
    body = [(r + [""] * len(header))[:len(header)] for r in rows[1:]]
    body = [r for r in body if r[0].strip()]
    lower = [h.lower() for h in header]
    wlay = [j for j, h in enumerate(lower) if "lay" in h] or [j for j, h in enumerate(lower) if "cell" in h]
    if not wlay:
        return ds
    ds["nlayers"] = len(wlay)
    wunits = [j for j, h in enumerate(lower) if "unit" in h]
    wlabel = [j for j, h in enumerate(lower) if "label" in h]
    if "pt" not in ds:
        ds["pt"], ds["pt2"] = get_periodic_table()
    if "solute" not in ds:
        ds["solute"] = set_model_species(ds["pt2"])
    if "expID" in lower:
        raise NotImplementedError("readinitlayerdata expID screening is not ported yet")

    layer = ds.setdefault("layer", {})

    # MATLAB's readtable takes leading rows with text in the layer columns
    # (e.g. "horizon") as extra header lines, so they never reach the main
    # loop. readinitlayerdata.m picks them up again from readcell, storing
    # the raw values without units, label or source.
    nskip = 0
    while nskip < len(body) and _numeric([body[nskip][j].strip() for j in wlay]) is None:
        nskip += 1
    for r in body[:nskip]:
        name = r[0].strip()
        if name.lower() != "parameter":
            layer[name] = [r[j].strip() for j in wlay]
    body = body[nskip:]

    params = [r[0].strip() for r in body]
    units = [r[wunits[0]].strip() if wunits else "" for r in body]

    def values(i):
        return [body[i][j].strip() for j in wlay]

    def soilsoln(i, name):
        sp = name.split("_")[0]
        if sp not in ("H", "pH") and sp not in ds["solute"]:
            print(f"(readinitlayerdata/getsoilsolndata) {sp} is not in the model")
            return
        tval = fill_nan(_numeric(values(i)))
        ss = layer.get("soilsoln")
        if sp != "pH":
            ss = layer.setdefault("soilsoln", {})
            if "_hotwater" in name:
                ss.setdefault("label", {})["msp"] = "porewater (hot water extraction)"
            elif "_soln" in name:
                ss.setdefault("label", {})["msp"] = "porewater"
            ss.setdefault("units", {})["msp"] = units[i]
        if ss is None:
            ss = layer.setdefault("soilsoln", {})
            ss.setdefault("label", {})["msp"] = "porewater"
            ss.setdefault("units", {})["msp"] = units[i]
        ss.setdefault("msp", {})[sp] = tval
        spas = ""
        key = sp
        if sp.lower() == "nh4":
            key, spas = "Amm", " as AmmH+"
        elif sp.lower() == "no3":
            key, spas = "N", " as NO3-"
        elif sp.lower() == "so4":
            key, spas = "S", " as SO4-2"
        ss.setdefault("spas", {})[key] = spas

    def exch(i, name):
        sp = name.split("_")[0]
        sol = ds["solute"].get(sp)
        if sp != "H" and sol is None:
            print(f"(readinitlayerdata/getexchdata) {sp} is not in the model")
            return
        if sp != "H" and "exchangespecies" not in sol:
            print(f"(readinitlayerdata/getexchdata) {sp} is not an exchangeable species in the model")
            return
        exchsp = sol["exchangespecies"]
        tval = fill_nan(_numeric(values(i)))
        u = units[i]
        lunits = layer.setdefault("units", {})
        if "avail" in name:
            layer.setdefault("available", {})[sp] = tval
            lunits["available"] = u
        if u.lower() == "%cec":
            layer.setdefault("basesat", {})[sp] = tval
            if "X" in exchsp:
                layer.setdefault("exchangeable", {})[sp] = tval
                lunits["exchangeable"] = "%CEC"
        elif ("meq" in u or "cmol" in u) and "X" in exchsp:
            layer.setdefault("exchangeable", {})[sp] = tval
            lunits["exchangeable"] = u
        elif "mg" in u and "kg" in u and "X" in exchsp:
            layer.setdefault("exchangeable", {})[sp] = tval / sol["molmass"]
            layer.setdefault("charge", {})[exchsp] = sol["charge"]
            lunits["exchangeable"] = u
        elif "X" in exchsp:
            layer.setdefault("charge", {})[exchsp] = sol["charge"]
            layer.setdefault("exchangeable", {})[sp] = tval
            lunits["exchangeable"] = u

    for i, name in enumerate(params):
        if name in layer:
            continue  # already set by an earlier file (no expID match to override it)
        if "_soln" in name or "_hotwater" in name or name.lower() == "ph":
            soilsoln(i, name)
        elif "_basesat" in name or "_exch" in name or "_avail" in name:
            exch(i, name)
        else:
            num = _numeric(values(i))
            layer[name] = fill_nan(num) if num is not None else values(i)
        layer.setdefault("source", {})[name] = str(fname)

    if wunits:
        layer.setdefault("units", {}).update({p: r[wunits[0]].strip() for p, r in zip(params, body)})
    if wlabel:
        layer.setdefault("label", {}).update({p: r[wlabel[0]].strip() for p, r in zip(params, body)})
    layer.setdefault("layername", [header[j] for j in wlay])
    layer["middepth"] = (layer["topdepth"] + layer["bottomdepth"]) / 2
    ds.setdefault("inputfiles", []).append(str(fname))
    return ds
