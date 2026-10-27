"""Run set-up after the soil cells and phases are defined. Python counterpart of
setuprun.m (with setmixfracs2.m, getsoilCO2g.m) and setupnutuptake.m.

Only the paths used by the published runs are ported: Richards hydrology with
top-down mixing, and soil CO2 from a soil-respiration time series (Cerling
1991) or the atmosphere.
"""

import math

import numpy as np

from .csvdata import read_csv_data
from .fixes import enabled as fixed
from .richards import get_rich_quick
from .swat import day_swat

ET_METHOD_LABELS = {
    "MIXflow": "Removing ET using MIX without uptake from transpiration stream",
    "MIXflownut": "Removing ET using MIX allowing uptake from transpiration stream",
    "MIXpure": "Removing ET using MIX with negative pure water",
    "REACTION": "Removing ET from cell water in REACTION datablock as in Ex. 4 PHREEQC documentation",
    "KINETICS": "Removing ET from cell water in KINETICS datablock",
}


def setup_run(ds, et_method="MIXflow", co2scl=1.0, legacy_airpor=False, **swat_kw):
    """setuprun.m as called by setupdatastruc (``'cellwatermethod','fracET'``).

    ``legacy_airpor`` reproduces the code that made the published runs (before
    LLT's 21 June 2026 change): getrichquick did not set ``run.airpor``, so
    getsoilCO2g used max(por - fc, por - Sws) for the CO2 diffusion.
    """
    et_method = ds.get("ETmethod", et_method)
    if et_method not in ET_METHOD_LABELS:
        raise ValueError(f"(setuprun/setETmethod) Unknown ETmethod: {et_method}")
    ds["ETmethod"], ds["ETmethodlabel"] = et_method, ET_METHOD_LABELS[et_method]
    run, parm = ds["run"], ds["parm"]
    scale = ds.get("scale", {})
    for k, v in list(run.items()):
        if not isinstance(v, np.ndarray) or v.ndim != 1 or not np.issubdtype(v.dtype, np.number):
            continue
        if k in parm:
            run[k] = v * 0 + parm[k]
        elif k in scale:
            run[k] = v * scale[k]
    ds = day_swat(ds, **swat_kw)
    ds["cellwatermethod"] = "fracET"
    ds = set_mix_fracs(ds)
    if legacy_airpor:
        del ds["run"]["airpor"]
    run = ds["run"]
    run["sat"] = run["Sws"] / np.asarray(ds["cell"]["por"], float)[:, None]
    run["units"]["sat"], run["label"]["sat"] = "fraction", "Saturation"
    soilco2opt = parm.get("soilco2opt")
    if soilco2opt is None:
        soilco2opt = "atmospheric"
        crop = parm.get("crop", ds.get("crop", ""))
        if "soilresp" in run:
            if np.nanmax(run["soilresp"]) > 0:
                soilco2opt = "timeseries"
            elif str(crop).lower() in ("sorghum", "oil palm"):
                raise NotImplementedError(f"soil respiration for {crop} is not ported")
        elif "totmassinc" in run and np.nanmax(run["totmassinc"]) > 0:
            raise NotImplementedError("Raich soil respiration is not ported")
        parm["soilco2opt"] = soilco2opt
    if "pCO2g" in scale:
        co2scl = scale["pCO2g"]
    set_soil_co2(ds, soilco2opt, co2scl)
    return ds


# ---------------------------------------------------------------- mixing
def set_mix_fracs(ds):
    """setmixfracs2.m for Richards hydrology and top-down mixing."""
    if ds.get("HYDROmethod", "bucket").lower() != "richards":
        raise NotImplementedError("bucket hydrology is not ported")
    if ds.get("MIXdir", "topdown") != "topdown":
        raise NotImplementedError("bottom-up mixing is not ported")
    run = ds["run"]
    if "infil" not in run:
        raise NotImplementedError("infiltration must come from the forcing data")
    ds = get_rich_quick(ds)
    run = ds["run"]
    nc, n = ds["ncells"], ds["nrun"]
    mmscl = ds["units"]["mmscl"]
    ds.setdefault("nmix", nc + 1)
    ds.setdefault("nmixsoln", nc + 2)
    ds["MIXture"] = [f"MIX {i}" for i in range(1, ds["nmix"] + 1)]
    ds["MIXsoln"] = np.arange(ds["nmixsoln"], dtype=float)
    ds["MIXsolnstr"] = [f"soln {i}" for i in range(ds["nmixsoln"])]
    run.setdefault("rhizon", np.zeros((nc, n)))
    run["bypassin"] = np.zeros((nc, n))
    run.setdefault("qdrain", np.zeros((nc, n)))
    inperc, outperc, targ, infil = run["inperc"], run["outperc"], run["targH2O"], run["infil"]
    inflow = np.fmax(0, inperc) + run["bypassin"]
    outflow = np.fmax(0, outperc) + run["qdrain"]
    upflow = -np.fmin(0, outperc)
    run["cellFC"] = cellFC = targ + inflow
    wleach = np.flatnonzero(outperc[-1] > 0)
    wdrain = np.flatnonzero(run["qdrain"].sum(axis=0) > 0)
    nleach = nc + 1
    if len(wdrain) > 0:
        ds["nmix"] = ndrain = nleach + 1
        ds["MIXture"].append(f"MIX {ndrain}")
    else:
        ds["nmix"], ndrain = nleach, 0
        ds["MIXture"] = ds["MIXture"][:nleach]
    ds["nmixleach"], ds["nmixdrain"] = nleach, ndrain
    mf = np.zeros((ds["nmix"], ds["nmixsoln"], n))
    run.setdefault("exptday", np.arange(1, n + 1, dtype=float))

    def m(i, j):  # MATLAB mixfrac(i,j,:)
        return (i - 1, j - 1)

    with np.errstate(divide="ignore", invalid="ignore"):
        mf[m(1, 1)] = inperc[0] / infil
        postmix = np.zeros((nc, n))
        mf[m(1, 2)] = 1
        mf[m(1, 3)] = upflow[1] / targ[1]
        for i in range(2, nc):
            mf[m(i, i)] = np.fmax(0, inperc[i - 1]) / cellFC[i - 2]
            mf[m(i, i + 1)] = 1
            mf[m(i, i + 2)] = inflow[i] / targ[i]
        i = nc
        mf[m(i, i)] = np.fmax(0, inperc[i - 1]) / cellFC[i - 2]
        mf[m(i, i + 1)] = (targ[i - 1] + upflow[i - 1]) / targ[i - 1]
        # overwrites the values above with upflow(i)/targH2O(i) (see PORTING_PLAN.md)
        if fixed(ds, "upflow_fractions"):  # upflow(i) comes from cell i+1's water
            for i in range(1, nc):
                mf[m(i, i + 2)] = upflow[i - 1] / targ[i]
        else:
            for i in range(1, nc - 1):
                mf[m(i, i + 2)] = upflow[i - 1] / targ[i - 1]
        for i in range(1, nc + 1):
            k = i - 1
            mf[k] = mf[k] * (cellFC[k] + upflow[k] - outflow[k]) / cellFC[k]
            postmix[k] = infil * mf[k, 0] + targ[k] * mf[k, i]
            if i > 1:
                postmix[k] = postmix[k] + cellFC[k - 1] * mf[k, i - 1]
        mf[nleach - 1, nc] = np.fmax(0, outperc[nc - 1]) / (cellFC[nc - 1] - outflow[nc - 1] + upflow[nc - 1])
        if ndrain > 0:
            for j in range(1, nc + 1):
                mf[ndrain - 1, j] = run["qdrain"][j - 1] / cellFC[j - 1]
    mf[mf < 0] = 0
    run["mixfrac"], run["postmix"] = mf, postmix
    leachmix = np.zeros(n)
    for k in wleach:
        leachmix[k] = infil[k] * mf[ds["nmix"] - 1, 0, k] + np.sum(cellFC[:, k] * mf[ds["nmix"] - 1, 1:nc + 1, k])
    run["leachmix"] = leachmix
    # Days with cellFC != targH2O outside the irrigation days (all days here) need no fix
    w = cellFC == 0
    cellFC[w] = targ[w]
    dET = np.fmax(0.01, cellFC - run["cellET"] - run["rhizon"])
    if "MIXflow" in ds["ETmethod"]:
        run["ETmix"] = dET / cellFC
        run["postmixET"] = cellFC * run["ETmix"]
    run["postET"] = dET.sum(axis=0)
    label, units = run["label"], run["units"]
    label["postirr"] = "Post-irrigation moisture minus drainage"
    label["postET"] = "Post-ET moisture"
    units["postirr"] = units["postET"] = mmscl
    label["Sws"], units["Sws"] = "Soil moisture", "V_{H2O} V_{soil}^{-1}"
    label["cellET"], units["inperc"] = "Evapotranspiration", f"{mmscl} cell^{{-1}}"
    return ds


# ---------------------------------------------------------------- soil gases
def _set_label(run, cell, k, label, units):
    run["label"][k], run["units"][k] = label, units
    if cell is not None:
        cell.setdefault("label", {})[k] = label
        cell.setdefault("units", {})[k] = units


def set_soil_co2(ds, soilco2opt, co2scl=1.0):
    """setuprun/setsoilco2 (with prescribesoilpO2g and prescribesoilpCO2g)."""
    run, cell, parm = ds["run"], ds["cell"], ds["parm"]
    nc, n = ds["ncells"], ds["nrun"]
    if "CO2atm" in parm:
        run["CO2atm"] = np.zeros(n) + parm["CO2atm"]
    elif "CO2atm" not in run:
        parm["CO2atm"] = parm.get("pCO2g", 420)
        run["CO2atm"] = np.zeros(n) + parm["CO2atm"]
    run["pCO2g"] = np.zeros((nc, n))
    run["units"].setdefault("pCO2g", "ppmV")
    if "pO2g" in parm:
        cell["pO2g"] = np.zeros(nc) + parm["pO2g"]
        run["pO2g"] = np.zeros((nc, n)) + parm["pO2g"]
        parm["soilo2opt"] = "prescribed"
        _set_label(run, cell, "pO2g", "Soil pO2g", "atm")
    else:
        if "pO2g" not in cell:
            parm["soilo2opt"] = "default"
        # an existing cell.pO2g is not prescribed: setsoilco2 resets it to 0.21 anyway
        cell["pO2g"] = np.zeros(nc) + 0.21
        run["pO2g"] = np.zeros((nc, n)) + 0.21
        _set_label(run, cell, "pO2g", "Atmospheric O_2(g) throughout soil", "atm")
    cell["si_O2g"] = np.log10(cell["pO2g"])
    if "pCO2g" in parm or "pCO2g" in cell:
        raise NotImplementedError("prescribed soil pCO2 is not ported")
    parm["soilco2opt"] = soilco2opt
    if soilco2opt == "atmospheric":
        run["pCO2g"] = np.tile(run["CO2atm"], (nc, 1))
        if "ppm" in run["units"].get("CO2atm", ""):
            run["pCO2g"] = run["pCO2g"] / 1e6
            run["units"]["pCO2g"] = "atm"
        run["label"]["pCO2g"] = "Atmospheric pCO_2 throughout soil"
    elif soilco2opt == "timeseries":
        get_soil_co2g(ds)
        run["label"]["pCO2g"] = "Soil CO_2 gas"
    else:
        raise NotImplementedError(f"soilco2opt {soilco2opt} is not ported")
    run["pCO2g"] = run["pCO2g"] * co2scl
    ds.setdefault("scale", {})["pCO2g"] = co2scl
    if co2scl != 1:
        run["label"]["pCO2g"] = f"{run['label']['pCO2g']}, X {co2scl:g}"
    run["si_CO2g"] = np.log10(run["pCO2g"])
    run["units"]["si_CO2g"] = f"log_{{10}}({run['units']['pCO2g']})"
    run["label"]["si_CO2g"] = run["label"]["pCO2g"]
    run["si_O2g"] = np.log10(run["pO2g"])
    run["units"]["si_O2g"] = f"log_{{10}}({run['units']['pO2g']})"
    run["label"]["si_O2g"] = run["label"]["pO2g"]
    return ds


def get_soil_co2g(ds):
    """getsoilCO2g.m: soil pCO2 from soil respiration after Cerling (1991) eqn 11."""
    run, cell, parm = ds["run"], ds["cell"], ds["parm"]
    nc, n = ds["ncells"], ds["nrun"]
    por = np.asarray(cell["por"], float)[:, None]
    if "airpor" not in run and "Sws" in run:
        if "fc" in cell:
            run["airpor"] = np.fmax(por - np.asarray(cell["fc"], float)[:, None], por - run["Sws"])
        else:
            run["airpor"] = np.fmax(0.09, por - run["Sws"])
        run["trueairpor"] = np.fmax(0, por - run["Sws"])
    if "airpor" in run and "watertabledepth" not in run:
        raise NotImplementedError("water table from air porosity is not ported")
    if "impermeabledepth" not in run and "watertabledepth" in run:
        run["impermeabledepth"] = run["watertabledepth"]
    co2atm_units = run["units"].get("CO2atm", "ppmV")
    airpressure, tortuosity, Dair0, Dexp, R = 1.0, 0.6, 0.144, 1.823, 82.057338
    P0, TK0 = 1.0, 273.15
    soilresp = np.asarray(run["soilresp"], float)
    TKsoil = TK0 + run["soiltemp"]
    airpor = run["airpor"]
    tprat = P0 * TKsoil / airpressure / TK0
    Dair = Dair0 * tprat ** Dexp
    Dsoil = Dair * tortuosity * airpor
    aircm3mol = R * TKsoil / airpressure
    molcmdepth = airpor / aircm3mol
    # getcdepth: heterotrophic respiration depth from the SOC profile
    tcdepth = np.asarray(run["tcdepth"], float)
    if "hetdepth" in parm:
        hetdepth = parm["hetdepth"]
    elif "gSOC" in cell:
        gSOCcm = np.asarray(cell["gSOC"], float) / np.asarray(cell["thickness"], float)
        hetdepth = np.sum(gSOCcm * np.asarray(cell["middepth"], float)) / np.sum(gSOCcm)
    else:
        wgs = tcdepth > 0
        hetdepth = np.nanmean(tcdepth[wgs])
    with np.errstate(divide="ignore", invalid="ignore"):
        run["cdepth"] = (run["rootresp"] * tcdepth + run["hetresp"] * hetdepth) / soilresp
    cdepth = np.zeros(n) + parm["cdepth"] if "cdepth" in parm else run["cdepth"]
    impermeabledepth = np.asarray(run["impermeabledepth"], float)
    soildepth = np.asarray(cell["middepth"], float)[:, None]
    SR = soilresp / 10000 / 86400 * 1e6
    phi0 = SR / cdepth / (1 - np.exp(-impermeabledepth / cdepth))
    Sz11 = phi0 / Dsoil * cdepth ** 2 * (1 - np.exp(-soildepth / cdepth))
    # p.SL is overwritten on every day of the loop, so only the last day's terms survive
    SL = phi0 * soildepth * cdepth[-1] / Dsoil[:, -1:] * np.exp(-impermeabledepth / cdepth[-1])
    C1 = Sz11 / molcmdepth
    run["Sz11"], run["SL"], run["C1"] = Sz11, SL, C1
    if co2atm_units.lower() == "atm":
        run["pCO2g"] = C1 / 1e6 + run["CO2atm"]
    else:
        run["pCO2g"] = (C1 + run["CO2atm"]) / 1e6
    run["units"]["pCO2g"] = "atm"
    run["soilCO2gppmV"] = run["pCO2g"] * 1e6
    run["phi0"] = phi0
    return ds


# ---------------------------------------------------------------- nutrients
def setup_nut_uptake(ds, usecrop=None):
    """setupnutuptake.m: Michaelis-Menten parameter row for each nutrient and day."""
    from pathlib import Path

    run = ds["run"]
    path = Path(ds.get("csvdir", "")) / "MMuptakeparams.csv"
    if "csvdir" not in ds or not path.exists():
        path = Path(ds["dbdir"]) / "MMuptakeparams.csv"
    mm = ds["MMuptake"] = read_csv_data(path)
    ds["nut"] = [k for k in run["totuptake"] if k not in ("label", "units")]
    ds["uptakeformula"] = {}
    crops = list(mm["crop"])
    nutrients = list(mm["nutrient"])
    for nut in ds["nut"]:
        formula = f"{nut} 1"
        if nut in ds["pt2"]:
            charge = ds["pt2"][nut]["charge"]
            if charge < 0:
                formula += " OH- -1"
            elif charge > 0:
                formula += " H+ -1"
        ds["uptakeformula"][nut] = formula
        idx = (np.asarray(run["plant"]) > 0).astype(float)
        for i, crop in enumerate(run["crop"]):
            if crop == "fallow":
                continue
            test = usecrop if usecrop is not None else crop
            if fixed(ds, "crop_case"):
                w = [r for r, (c, nu) in enumerate(zip(crops, nutrients)) if c.lower() == test.lower() and nu == nut]
            else:
                w = [r for r, (c, nu) in enumerate(zip(crops, nutrients)) if c == test and nu == nut]
            if w:
                idx[i] = w[0] + 1
        run[f"{nut}index"] = idx
    return ds


def set_residues(ds):
    """setupdatastruc.m after setupnutuptake: crop residues and root residues left at harvest."""
    run, cell, parm = ds["run"], ds["cell"], ds["parm"]
    till = np.asarray(cell["tillagefrac"])
    if parm.get("usemeanforcings") == 1 and "peaktot" in run:
        zero = np.zeros((ds["ncells"], ds["nrun"]))
        run["rootresidue"] = zero + np.mean(run["peakroot"]) * run["cellup"]
        run["residue"] = zero + max(0, np.mean(run["peaktot"] - run["harvest"] - run["peakroot"])) * till[:, None]
    else:
        w = np.flatnonzero(np.asarray(run["harvest"]) > 0)
        if len(w):
            run["rootresidue"][:, w] = run["rootmass"][w] * run["cellup"][:, w]
            run["residue"][:, w] = np.fmax(0, run["totmass"][w] - run["harvest"][w] - run["rootmass"][w]) * till[:, None]
    return ds


def complete_setup(ds, legacy_airpor=False):
    """The end of setupdatastruc.m: setuprun, setupnutuptake and the harvest residues."""
    ds.setdefault("units", {"mmscl": "L", "scl": "m^{-2}land"})  # defaultsettings: ds.mmscl, ds.scl
    ds = setup_run(ds, legacy_airpor=legacy_airpor)
    return set_residues(setup_nut_uptake(ds))
