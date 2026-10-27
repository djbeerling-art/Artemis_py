"""Feedstock and secondary mineral phases. Python counterpart of setupphases.m
and setnucSA.m.

``setup_phases`` builds ``ds["feedstocktable"]`` (one row per feedstock that
the timetable applies), the per-phase run arrays ``ds["run"]["moles"]`` and
``ds["run"]["RSA"]`` (cells x days), the parsed stoichiometries
``ds["stoich"]``, and finally re-checks the phase types, which builds the
implicit pseudo-phase (implicitcalc).

Ported for the paths the EF runs use: a "vintage" particle size file and the
BET surface area. Rosin-Rammler and Rock-MIP size distributions, timetabled
"phase" treatments and molar masses computed from formulas raise
NotImplementedError.
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd

from .csvdata import read_csv_data
from .formula import hydrate_to_formula, parse_formula
from .phasetable import check_phase_type

PRECIPDISSRAT = {"calcite": 1e-7, "sio2a": 1e-6, "kaolinite": 0.0498, "goethite": 1e-6}


def _units(t):
    u = t.attrs.setdefault("units", {})
    return u


def setup_phases(ds, basedir=None):
    """setupphases.m. ``basedir`` resolves a relative ``ds["feedstockPSDfile"]``."""
    parm = ds["parm"]
    parm.setdefault("RSAscale", 1)
    parm.setdefault("useBET", 1)
    afracfld = "afrac"
    if "feedstock" in parm:
        raise NotImplementedError("parm.feedstock (one feedstock for the whole run) is not ported yet")
    ds.pop("feedstock", None)
    psdfile = ds.get("feedstockPSDfile", "")
    if not psdfile:
        print("(setupphases) No PSD filename given")
    if "feedstock" not in ds["phasetable"].columns:
        print("(setupphases) No feedstock names given, returning")
        return ds
    _init_phase_mol_arrays(ds)
    _set_stoichiometry(ds)
    _prepare_tables(ds, afracfld)
    _check_timetable(ds)
    _get_psd_data(ds, psdfile, basedir)
    _set_rsa(ds, afracfld)
    if any(str(d).lower() == "phase" for d in ds["timetable"].get("doseform", [])):
        raise NotImplementedError("timetabled phase treatments (checktimetablephase) are not ported yet")
    ds["cell"]["nucSA"], ds["cell"]["nucN"] = set_nuc_sa(ds)
    t = ds["phasetable"]
    t["precipdissrat"] = 1.0
    for i, name in enumerate(t["MATLABphase"]):
        if str(name).lower() in PRECIPDISSRAT:
            t.iloc[i, t.columns.get_loc("precipdissrat")] = PRECIPDISSRAT[str(name).lower()]
    ds = check_phase_type(None, ds, quiet=True)
    if "feedstocktable" not in ds:
        print("(setupphases) No feedstock table")
    return ds


def _init_phase_mol_arrays(ds):
    """initphasemolarrays: zeroed cells x days arrays, with day 1 from the soil cells."""
    t, cell, run = ds["phasetable"], ds["cell"], ds.setdefault("run", {})
    shape = (ds["ncells"], ds["nrun"])
    moles, rsa = run.setdefault("moles", {}), run.setdefault("RSA", {})
    initmol = pd.DataFrame(0.0, index=range(ds["ncells"]), columns=list(t["MATLABphase"]))
    initrsa = initmol.copy()
    cmoles = cell.get("moles")
    for i, pname in enumerate(t["MATLABphase"]):
        z = np.zeros(shape)
        rsa[pname] = np.zeros(shape)
        if cmoles is not None:
            if pname in cmoles:
                initmol[pname] = np.asarray(cmoles[pname], float)
                initrsa[pname] = np.asarray(cell["RSA"][pname], float)
                z[:, 0] = initmol[pname]
        elif "InitialMoles" in t and np.isfinite(t["InitialMoles"].iloc[i]):
            z[:, 0] = t["InitialMoles"].iloc[i] * np.asarray(cell["frac"], float)
            initmol[pname] = z[:, 0]
        moles[pname] = z
    initmol.attrs["units"] = "mol/cell"
    initrsa.attrs["units"] = "m2/cell"
    ds["initphasemoltable"], ds["initphaseRSAtable"] = initmol, initrsa


def _set_stoichiometry(ds):
    """setstoichiometry: element counts per phase (StoichTools parse_formula)."""
    t = ds["phasetable"]
    if "formula" not in t:
        print("(setupphases/setstoichiometry) Cannot determine phase stoichiometries without chemical formulas, returning")
    formulas = [hydrate_to_formula(str(f)) for f in t["formula"]]
    if "molmass" not in t or np.any(~np.isfinite(t["molmass"]) | (t["molmass"] <= 0)):
        raise NotImplementedError("molar masses from formulas (StoichTools molweight) are not ported yet")
    stoich = {}
    for name, f in zip(t["MATLABphase"], formulas):
        try:
            stoich[name] = parse_formula(f)
        except ValueError as err:
            print(f"(setupphases/setstoichiometry) Error determining stoichiometry for {name}")
            print(err)
    ds["stoich"] = stoich


def _phase_fractions(t, afracfld):
    """getphasefractions: grain volume and sphere-equivalent area per gram of feedstock."""
    if afracfld in t:
        return t
    if "frac" not in t:
        if "wgtpct" in t:
            t["frac"] = t["wgtpct"] / 100
        elif "weight" in t:
            t["frac"] = t["weight"] / 100
    if afracfld in t:
        return t
    if "G" not in t and "density" in t:
        t["G"] = t["density"]
    if "G" in t and "frac" in t:
        t["vol_cm3pg"] = t["frac"] / t["G"]
        t["area_cm2pg"] = 4 * np.pi * (t["vol_cm3pg"] * 3 / 4 / np.pi) ** (2 / 3)
        _units(t).update(vol_cm3pg="cm3phase/gfeedstock", area_cm2pg="cm2phase/gfeedstock")
    return t


def _prepare_tables(ds, afracfld):
    """preparetables: feedstock table of summed fractions, densities and molar masses."""
    t = ds["phasetable"]
    if "affinity" not in t:
        t["affinity"] = 1.0
    if "category" not in t:
        t["category"] = ""
    if "PHREEQCphasetype" not in t:
        raise NotImplementedError("phase tables without PHREEQCphasetype are not ported yet")
    fs = t["feedstock"].fillna("").astype(str)
    wfs = fs != ""
    if not wfs.any():
        return
    t.loc[wfs, "category"] = "feedstock"
    ufeedstock = sorted(set(fs[wfs]))
    fstable = pd.DataFrame({"feedstock": ufeedstock}, index=ufeedstock)
    fstable.attrs["units"] = {}
    t = _phase_fractions(t, afracfld)
    tunits = _units(t)
    sumfld = ["frac", "vol_cm3pg", "area_cm2pg", "G", "molmass"]
    fracfld = ["frac", "vfrac", "afrac", "G", "molmass"]
    for ifs, name in enumerate(ufeedstock):
        rows = (fs == name).to_numpy()
        for sf, ff in zip(sumfld, fracfld):
            if sf not in t:
                continue
            tot = "tot" + sf
            if tot not in fstable:
                fstable[tot] = 0.0
            if sf == "frac":
                fstable.iloc[ifs, fstable.columns.get_loc(tot)] = t.loc[rows, sf].sum(skipna=False)
            elif "pg" in sf:
                total = t.loc[rows, sf].sum(skipna=False)
                fstable.iloc[ifs, fstable.columns.get_loc(tot)] = total
                if ff not in t:
                    t[ff] = 0.0
                t.loc[rows, ff] = t.loc[rows, sf] / total
                totfr = "tot" + ff
                if totfr not in fstable:
                    fstable[totfr] = 0.0
                fstable.iloc[ifs, fstable.columns.get_loc(totfr)] = t.loc[rows, ff].sum(skipna=False)
                fstable.attrs["units"][totfr] = tunits.get(sf, "")
            else:
                fstable.iloc[ifs, fstable.columns.get_loc(tot)] = (t.loc[rows, sf] * t.loc[rows, "frac"]).sum(skipna=False)
            fstable.attrs["units"][tot] = tunits.get(sf, "")
    gunits = fstable.attrs["units"].get("totG")
    if gunits is not None:
        if re.search(r"kg[/\s]*m", gunits):
            t["rho_kgm3"] = t["G"]
            fstable["rho_kgm3"] = fstable["totG"]
        elif re.search(r"g[/\s]*cm", gunits):
            # MATLAB divides by 1000 here; g/cm3 to kg/m3 is x1000 (see PORTING_PLAN)
            t["rho_kgm3"] = t["G"] / 1000
            fstable["rho_kgm3"] = fstable["totG"] / 1000
    ds["phasetable"], ds["feedstocktable"] = t, fstable


def _check_timetable(ds):
    """checktimetable: feedstock doses per cell and day, and moles added on application days."""
    if "feedstocktable" not in ds:
        return
    fstable, t, tt = ds["feedstocktable"], ds["phasetable"], ds["timetable"]
    run, cell = ds["run"], ds["cell"]
    if "treatment" not in tt and "basalt" in tt:
        tt["treatment"] = tt["basalt"]
    if "dosekg" not in tt:
        tt["dosekg"] = np.array(tt["baskg"], float) if "baskg" in tt else np.zeros(len(tt["runindex"]))
    run["gdust"] = np.zeros((ds["ncells"], ds["nrun"]))
    run.setdefault("units", {})["treatment"] = "g treatment cell^{-1}"
    run.setdefault("label", {})["treatment"] = "Treatment"
    treat = [str(x).lower() for x in tt["treatment"]]
    runindex = np.asarray(tt["runindex"])
    tfs = t["feedstock"].fillna("").astype(str).str.lower().to_numpy()
    tillage = np.asarray(cell["tillagefrac"], float)
    bad = []
    for fsname in list(fstable["feedstock"]):
        wtt = [j for j, x in enumerate(treat) if x == fsname.lower() and runindex[j] > 0]
        if not wtt or not np.any(tfs == fsname.lower()):
            bad.append(fsname)
            continue
        ri = runindex[wtt].astype(int) - 1
        rt = run.setdefault("treatment", {})
        rt[fsname] = np.zeros((ds["ncells"], ds["nrun"]))
        if "doseunits" in tt:
            for j in wtt:
                du = tt["doseunits"][j]
                if du.lower() == "t/ha":
                    tt["dosekg"][j] = tt["dosage"][j] / 10
                elif "kg/ha" in du:
                    tt["dosekg"][j] = tt["dosage"][j] / 10000
                elif "kg/m2" in du:
                    tt["dosekg"][j] = tt["dosage"][j]
        dosekg = np.asarray(tt["dosekg"], float)[wtt]
        rt[fsname][:, ri] = tillage[:, None] * dosekg[None, :] * 1000
        cell.setdefault("treatment", {})[fsname] = rt[fsname][:, 0].copy()
        for i in np.flatnonzero(tfs == fsname.lower()):
            pname = t["MATLABphase"].iloc[i]
            m = run["moles"].setdefault(pname, np.zeros((ds["ncells"], ds["nrun"])))
            if "InitialMoles" in t and np.isfinite(t["InitialMoles"].iloc[i]):
                m[:, 0] = t["InitialMoles"].iloc[i] * np.asarray(cell["frac"], float)
            if np.isfinite(t["frac"].iloc[i]):
                m[:, ri] = m[:, [0]] + rt[fsname][:, ri] * t["frac"].iloc[i] / t["molmass"].iloc[i]
            cell.setdefault("moles", {})[pname] = m[:, 0].copy()
    if bad:
        print(f"(setupphases/checktimetable) Removing superfluous feedstocks: {' '.join(bad)}")
        fstable = fstable.drop(index=bad)
        badl = {b.lower() for b in bad}
        drop = [i for i, f in enumerate(tfs) if f in badl]
        for pname in t["MATLABphase"].iloc[drop]:
            for d in (run["RSA"], run["moles"], cell.get("RSA", {}), cell.get("moles", {})):
                d.pop(pname, None)
        t = t.drop(index=t.index[drop])
    ds["feedstocktable"], ds["phasetable"] = fstable, t


def _get_psd_data(ds, psdfile, basedir):
    """getPSDdata: particle size distributions and mean radii per feedstock."""
    if "feedstocktable" not in ds:
        return
    fstable = ds["feedstocktable"]
    if not len(fstable):
        print("(setupphases/getPSDdata) No feedstocks kept!!")
        return
    fstable["totPSD"] = 0.0
    fstable["meanrad"] = 0.0
    if "PSDtype" not in fstable:
        fstable["PSDtype"] = ""
    if "usePSD" in ds:
        fstable["PSDname"] = ds["usePSD"]
    if "PSDname" not in fstable:
        fstable["PSDname"] = fstable["feedstock"]
    if "RRspread" in ds["parm"] and "P80" in ds["parm"]:
        raise NotImplementedError("Rosin-Rammler particle size distributions are not ported yet")
    for f in [psdfile] if isinstance(psdfile, str) else psdfile:
        if not f:
            continue
        path = Path(f)
        if basedir is not None and not path.is_absolute():
            path = Path(basedir) / path
        with open(path, encoding="utf-8-sig") as fh:
            first = fh.readline().split(",")[0].strip()
        if first.lower() != "field":
            raise NotImplementedError("Rock-MIP particle size files are not ported yet")
        _get_vintage_psd(ds, path, f)


def _get_vintage_psd(ds, path, fname):
    """getvintagePSD: a readcsvdata-format file with one PSD column per feedstock."""
    fstable = ds["feedstocktable"]
    psd = read_csv_data(path)
    units = psd["units"]
    if "radius" not in psd and "diameter" in psd:
        psd["radius"] = psd["diameter"] / 2
        units["radius"] = units["diameter"]
    radunits = units["radius"]
    rscl = 1e6 if ("micro" in radunits.lower() or "mum" in radunits.lower()) else 1
    names = {str(n).lower() for n in fstable["PSDname"]}
    feedstock = ds.setdefault("feedstock", {})
    meta = {"units", "label", "format", "field"}
    for fld, col in psd.items():
        if fld in meta or isinstance(col, dict) or fld.lower() in ("radius", "diameter"):
            continue
        if fld.lower() not in names and not (ds.get("defaultPSD", "").lower() == fld.lower()):
            continue
        col = np.asarray(col, float)
        cdf = np.cumsum(col)
        table = {}
        if "diameter" in psd:
            table["diameter"] = np.asarray(psd["diameter"], float)
        table["radius"] = np.asarray(psd["radius"], float)
        table["r_m"] = table["radius"] / rscl
        table["PSD"], table["CDF"], table["cdf_frac"] = col, cdf, cdf / cdf.max()
        pct = "%" if col.sum() > 10 else "fraction"
        tunits = {"radius": radunits, "r_m": "m", "PSD": pct, "CDF": pct, "cdf_frac": "fraction"}
        if "diameter" in psd:
            tunits["diameter"] = units["diameter"]
        df = pd.DataFrame(table)
        df.attrs["units"] = tunits
        feedstock[fld] = {"PSDtable": df}
    for name in fstable.index:
        pname = fstable.at[name, "PSDname"]
        if pname not in feedstock:
            if ds.get("defaultPSD") in feedstock:
                pname = fstable.at[name, "PSDname"] = ds["defaultPSD"]
            else:
                continue
        tab = feedstock[pname]["PSDtable"]
        fstable.at[name, "totPSD"] = np.nansum(tab["PSD"])
        meanrad = np.nansum(tab["PSD"] * psd["radius"]) / np.sum(tab["PSD"])
        ru = tab.attrs["units"]["radius"]
        scale = {"micro": 1e-4, "milli": 0.1, "centi": 1}
        conv = next((v for k, v in scale.items() if k in ru), None)
        if conv is None:
            conv = {"mm": 0.1, "cm": 1, "m": 100}.get(ru.lower())
        if conv is None:
            print(f"(setupphases/getvintagePSD) Cannot convert {ru} to cm for feedstock {name}")
            conv = 1
        else:
            fstable.attrs["units"]["meanrad"] = "cm"
        fstable.at[name, "meanrad"] = meanrad * conv
        fstable.at[name, "PSDtype"] = "vintage"
    # Metadata rows (BET area, lambda, ...) become feedstock table columns.
    # MATLAB readcsvdata lets detectImportOptions read the all-numeric rows as
    # data, which lower-cases their names; the header pass then adds them again
    # with their original case. Both copies end up as columns (betn2_m2g and
    # BETN2_m2g), lower-case ones first.
    rows = [(f, c) for f, c in psd.items() if f not in meta and isinstance(c, dict)]

    def numeric(c):
        try:
            [float(v) for v in c.values()]
            return True
        except ValueError:
            return False

    fields = [(f.lower(), c) for f, c in rows if numeric(c)]
    fields += [(f, c) for f, c in rows if f not in {g for g, _ in fields}]
    for fld, col in fields:
        for name in fstable.index:
            pname = fstable.at[name, "PSDname"]
            if pname not in col:
                continue
            v = col[pname]
            try:
                v = float(v)
                if fld not in fstable:
                    fstable[fld] = 0.0
            except ValueError:
                if fld not in fstable:
                    fstable[fld] = ""
            fstable.at[name, fld] = v
    ds.setdefault("inputfiles", []).append(str(fname))


def _set_rsa(ds, afracfld):
    """setRSA: feedstock BET or geometric surface areas and per-phase RSA arrays."""
    if "feedstocktable" not in ds:
        return
    fstable, t, run, parm = ds["feedstocktable"], ds["phasetable"], ds["run"], ds["parm"]
    ptype = t["PHREEQCphasetype"].astype(str)
    wk = ptype.str.contains("eq") | ptype.str.contains("kin")
    if not wk.any():
        return
    for k in t.loc[wk, "MATLABphase"]:
        z = np.zeros((ds["ncells"], ds["nrun"]))
        z[:, 0] = ds["initphaseRSAtable"][k]
        run["RSA"][k] = z
    nfs = len(fstable)
    # parm.RSAscale always exists here (set at the top of setupphases)
    rsascale = np.full(nfs, float(parm["RSAscale"]))
    fstable["RSAscale"] = rsascale
    geo = bet = None
    if "particle_size" in parm:
        fstable["meanrad"] = parm["particle_size"] / 2 / 1e4
        fstable.attrs["units"]["meanrad"] = "cm"
        geo = 3 / fstable["meanrad"].to_numpy() / fstable["totG"].to_numpy() / 1e4
    if "meanrad" in fstable:
        fstable["meanrad_m"] = fstable["meanrad"] / 100
        fstable.attrs["units"]["meanrad_m"] = "m"
    cols = list(fstable.columns)
    # lambda always ends up defined (zero if not given), as in MATLAB
    lam = np.full(nfs, float(parm["lambda"])) if "lambda" in parm else (
        fstable["lambda"].to_numpy(float) if "lambda" in cols else np.zeros(nfs))
    if "geo_m2g" in parm:
        geo = np.full(nfs, float(parm["geo_m2g"]))
    elif "geo_m2g" in cols:
        geo = fstable["geo_m2g"].to_numpy(float)
    wbet = [c for c in cols if "BET" in c and "_m2g" in c]
    if "BET_m2g" in parm:
        bet = np.full(nfs, float(parm["BET_m2g"]))
    elif "BETN2_m2g" in parm:
        bet = np.full(nfs, float(parm["BETN2_m2g"]))
    elif wbet:
        bet = fstable[wbet[0]].to_numpy(float)
    if bet is None and geo is not None:
        bet = lam * geo
    elif bet is not None and geo is None:
        geo = lam * bet  # MATLAB: lambda = BET/geo, so this should be BET/lambda
    if bet is None:
        bet = np.ones(nfs)
    if geo is None:
        geo = np.full(nfs, 0.01)
    rsa_m2g = (bet if parm["useBET"] else geo) * rsascale
    fstable["lambda"], fstable["geo_m2g"], fstable["BET_m2g"] = lam, geo, bet
    fstable["RSA_m2g"] = rsa_m2g
    for c in ("RSA_m2g", "geo_m2g", "BET_m2g", "lambda", "RSAscale"):
        t[c] = 0.0
    tfs = t["feedstock"].fillna("").astype(str).str.lower().to_numpy()
    for ifs, fsname in enumerate(fstable["feedstock"]):
        rows = tfs == fsname.lower()
        if not rows.any():
            continue
        af = t.loc[rows, afracfld].to_numpy(float)
        t.loc[rows, "RSA_m2g"] = rsa_m2g[ifs] * af
        t.loc[rows, "geo_m2g"] = geo[ifs] * af
        t.loc[rows, "BET_m2g"] = bet[ifs] * af
        t.loc[rows, "lambda"] = lam[ifs]
        t.loc[rows, "RSAscale"] = rsascale[ifs]
        for kf in t.loc[rows, "MATLABphase"]:
            run["RSA"][kf] = run["treatment"][t.at[kf, "feedstock"]] * t.at[kf, "RSA_m2g"]
    parm["RSAweightingfield"] = afracfld
    run.setdefault("units", {})["RSA"] = "m^2 cell^{-1}"


def set_nuc_sa(ds, Ri=(1.025 / 20, 0.026 / 20, 0.001 / 20)):
    """setnucSA.m: mineral surface area (m2/cell) and particle count available for nucleation."""
    cell = ds["cell"]
    phi = np.asarray(cell["por"], float)[:, None]
    bd = np.asarray(cell["bulkden"] if "bulkden" in cell else cell["BD"], float)[:, None]
    gsoil = np.asarray(cell["gsoil"], float)
    if "sandpct" in cell:
        wi = np.column_stack([cell["sandpct"], cell["siltpct"], cell["claypct"]]) / 100
    else:
        wi = np.column_stack([cell["sandfrac"], cell["siltfrac"], cell["clayfrac"]])
    Ri = np.asarray(Ri)
    ni = 3 * wi / (4 * np.pi * 2.65 * Ri ** 3)
    rin = np.sqrt(0.0717 * phi * (wi / bd) / ni ** (4 / 3) / Ri)
    sai = ni * 4 * np.pi * Ri ** 2
    small = rin < 1e-5
    sai[small] = 0
    ni[small] = 0
    return gsoil * sai.sum(axis=1) / 1e4, gsoil * ni.sum(axis=1)
