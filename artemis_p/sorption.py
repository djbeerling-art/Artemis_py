"""Cation exchange and sorbing surfaces. Python counterpart of CECsorb.m and of
the setsoiltexture nested function in setupdatastruc.m.

Only the branches reachable from a ``ds.cell.CEC`` profile (with or without
``ds.cell.exchangeable`` data) and from a single ``parm.CEC`` value are
ported. NH4-acetate extraction data and base saturation tables raise
NotImplementedError.
"""

import numpy as np

from .fixes import enabled as fixed

CEC_OPTS = ("CECSURF", "CECWHAM", "CECcarboxylic")
BASE_CATIONS = ("ca", "mg", "na", "k")

WHAM = {
    "carbsite": ["H_a", "H_b", "H_c", "H_d"],
    "phensite": ["H_e", "H_f", "H_g", "H_h"],
    "dipsite": ["H_ab", "H_ad", "H_af", "H_ah", "H_bc", "H_be", "H_bg", "H_cd", "H_cf", "H_ch", "H_de", "H_dg"],
}


def set_soil_texture(ds):
    """setsoiltexture: fraction of the tillage depth in each cell."""
    parm, cell = ds["parm"], ds["cell"]
    if "tillagedepth" not in parm:
        parm["tillagedepth"] = parm.get("basalt_depth", np.asarray(cell["bottomdepth"]).ravel()[0])
    parm.setdefault("units", {})["tillagedepth"] = "cm"
    top = np.asarray(cell["topdepth"], float)
    cell["tillagefrac"] = np.fmax(np.fmin(parm["tillagedepth"] - top, cell["thickness"]) / parm["tillagedepth"], 0)
    cell.setdefault("label", {})["tillagefrac"] = "Fraction of tillage depth in each cell"
    return ds


def _per_cec_units(u):
    u = u.lower()
    return ("meq" in u and "100g" in u) or ("cmol" in u and "kg" in u)


def _exacidity_species(ds):
    """Exchangeable acidity species and its charge (default HX)."""
    parm, has_data = ds["parm"], False
    sp, chg = "HX", 1
    if "setexaciditysp" in parm:
        name = parm["setexaciditysp"]
        if name in ds["pt2"]:
            sp, chg = ds["pt2"][name]["exchangespecies"], ds["pt2"][name]["charge"]
        else:
            ex = [str(x).lower() for x in ds["pt"]["exchangespecies"]]
            if name.lower() in ex:
                k = ex.index(name.lower())
                sp, chg = ds["pt"]["exchangespecies"][k], ds["pt"]["charge"][k]
        has_data = name in ds.get("layer", {}).get("exchangeable", {})
    return sp, chg, has_data


def cec_sorb(ds, extractable_fe=None):
    """CECsorb: PHREEQC EXCHANGE species and SURFACE site amounts per cell."""
    parm, cell = ds["parm"], ds["cell"]
    units, label = cell.setdefault("units", {}), cell.setdefault("label", {})
    cecopt = parm.get("CECopt", "")
    gsoil = np.asarray(cell["gsoil"], float)
    n = ds["ncells"]
    exsp, exchg, exdata = _exacidity_species(ds)

    scl = 1
    # MATLAB: ds.scale takes precedence over parm.CECscl even when it has no CEC field
    if "scale" in ds:
        scl = ds["scale"].get("CEC", 1)
    elif "CECscl" in parm:
        scl = ds.setdefault("scale", {})["CEC"] = parm["CECscl"]

    cecsum = 0
    if "CEC" in cell:
        cec = np.asarray(cell["CEC"], float)
        if _per_cec_units(units["CEC"]):
            cell["X"] = cec * gsoil / 1000 / 100 * scl  # eq/cell
        if "exchangeable" in units:
            exu = units["exchangeable"]
            exch_in = cell["exchangeable"]
            if "mg" in exu and "kg" in exu:
                for sp in exch_in:
                    exch_in[sp] = exch_in[sp] / ds["pt2"][sp]["molmass"] / 10 * ds["pt2"][sp]["charge"]
                exu = units["exchangeable"] = exu.replace("mg", "cmolc")
            elif exu.lower() == "%cec":
                for sp in exch_in:
                    exch_in[sp] = exch_in[sp] / 100 * cec
                exu = units["CEC"]
            if _per_cec_units(exu):
                cell["basesatsum"] = np.zeros(n)
                cecsum = np.zeros(n)
                exchange = cell.setdefault("exchange", {})
                for sp, v in exch_in.items():
                    v = np.asarray(v, float)
                    if sp.lower() in BASE_CATIONS:
                        cell["basesatsum"] = cell["basesatsum"] + v
                    cecsum = cecsum + v
                    p = ds["pt2"][sp]
                    exchange[p["exchangespecies"]] = v / p["charge"] * gsoil / 1000 / 100 * scl  # mol/cell
                units["exchange"] = "mol cell^{-1}"
            cell["basesatsum"] = cell["basesatsum"] / cec * 100
            units["basesatsum"], label["basesatsum"] = "percent", "Base saturation"
            cell["cecsum"] = cecsum / cec * 100
            units["cecsum"], label["cecsum"] = "percent", "Sum contributions to CEC"
            if np.any(cecsum < cec) and not exdata:
                cell["exchange"][exsp] = np.fmax(1e-13, cec - cecsum) / cec * cell["X"] / exchg
    elif "exchangeable" in cell or "extractiondata" in parm or "basesat" in cell:
        raise NotImplementedError("CECsorb without a CEC profile is not ported yet")
    elif "CEC" in parm:
        cell["X"] = parm["CEC"] * gsoil / 1000 / 100
    else:
        print("(CECsorb) No CEC info given!!")
    _set_base_sat(ds, exsp)

    units["X"], label["X"] = "eq cell^{-1}", "Cation exchange capacity"
    units["exchange"] = "mol cell^{-1}"
    cell["CECcmolckgsoil"] = cell["X"] / gsoil * 1000 * 100
    units["CEC"], label["CEC"] = "cmolc kgsoil^{-1}", "Cation exchange capacity"
    # MATLAB: cecsum is still in meq/100g here, so eqCEC is 100x the eq/cell it is labelled as
    cell["eqCEC"] = cecsum * gsoil / 1000 / (100 if fixed(ds, "eqcec_units") else 1)
    units["eqCEC"], label["eqCEC"] = "eq cell^{-1}", "Cation exchange capacity"

    if "gSOC" not in cell and "SOC" in cell:
        cell["gSOC"] = np.asarray(cell["SOC"])
    surf = cell.setdefault("surface", {})
    su = units.setdefault("surface", {})
    surf["CARBOXYLICSITES"] = 2.84 / 4 / 1000 * np.asarray(cell["gSOC"], float) * parm["SOCactfrac"]
    su["CARBOXYLICSITES"] = "mol monoprotic carboxylic organic sites cell^{-1}"
    surf["PHENOLICSITES"] = surf["CARBOXYLICSITES"] / 2
    su["PHENOLICSITES"] = "mol monoprotic phenolic organic sites cell^{-1}"
    surf["DIPROTICSITES"] = surf["CARBOXYLICSITES"] / 3
    su["DIPROTICSITES"] = "mol diprotic organic sites cell^{-1}"
    wham = parm.setdefault("WHAM", {})
    for k in ("carb", "phen", "dip"):
        wham[f"{k}site"] = list(WHAM[f"{k}site"])
        wham[f"n{k}sites"] = len(WHAM[f"{k}site"])

    surf["Fe"] = _surface_fe(ds, extractable_fe)
    su["Fe"] = "g Fe g^{-1} soil"
    fe = gsoil * surf["Fe"] * parm["Feactfrac"]
    surf["Hfo_w"] = fe / 55.85 * 0.2
    surf["Hfo_s"] = fe / 55.85 * 0.005
    su["Hfo_s"], su["Hfo_w"] = "mol HFO strong sites cell^{-1}", "mol HFO weak sites cell^{-1}"
    surf["gHFO"] = fe * 88.858 / 55.85
    su["gHFO"] = "g HFO cell^{-1}"
    surf["HfoSA"] = 5.33e4 / 88.858
    cell["SOCeq"] = 7.1 * np.asarray(cell["gSOC"], float) * parm["SOCactfrac"] / 1000
    units["SOCeq"], label["SOCeq"] = "eq cell^{-1}", "SOC CEC"

    X = cell["X"]
    if cecopt.lower() == "cecsurf":
        frac = (X - cell["SOCeq"] - surf["Hfo_w"] - surf["Hfo_s"]) / X
    elif cecopt.lower() == "cecwham":
        frac = (X - cell["SOCeq"]) / X
    elif cecopt.lower() == "ceccarboxylic":
        frac = (X - surf["CARBOXYLICSITES"]) / X
    else:
        frac = np.ones(n)
    cell["XCECfrac"] = frac
    label["XCECfrac"] = "Fraction of CEC on PHREEQC exchangers"
    cell["X"] = X * frac
    if "exchange" in cell:
        for k in cell["exchange"]:
            cell["exchange"][k] = cell["exchange"][k] * frac
        label["exchange"] = "PHREEQC exchange species"
        units["exchange"] = "mol cell^{-1}"
    elif "comp" in ds.get("EXCHmethod", ""):
        print(f"(CECsorb) *****WARNING: No exchange data; PHREEQC may fail with ds.EXCHmethod = {ds['EXCHmethod']}")
    return ds


def _surface_fe(ds, extractable_fe):
    """Fe (g Fe/g soil) available for HFO surface complexation."""
    n, cell = ds["ncells"], ds["cell"]
    fe = np.zeros(n)
    if extractable_fe is not None:
        fe = extractable_fe / 1e3
    else:
        t = ds["phasetable"]
        formula = t["formula"].astype(str)
        moles = cell.get("moles", {})
        for idx in t.index[formula.str.contains("FeOOH", regex=False) | formula.str.contains("Fe(OH)3", regex=False)]:
            ph = t.loc[idx, "PHREEQCphase"]
            if ph in moles:
                # MATLAB divides moles by 100 as if they were a mass percent
                fe = fe + np.asarray(moles[ph]) / 100 * 55.85 / t.loc[idx, "molmass"]
                print(f"(CECsorb) Fe for surface complexation from {ph}")
    if np.sum(fe) == 0:
        fe = fe + 2.79 / 1000
        print(f"(CECsorb) Assuming Fe available for surface complexation: {2.79 / 1000:g} gFe/gsoil")
    return fe


def _set_base_sat(ds, exsp):
    """setbasesat: optionally overwrite exchange with a fixed base saturation."""
    parm, cell = ds["parm"], ds["cell"]
    if "setbasesatpct" not in parm:
        return
    bsfrac = parm["setbasesatpct"] / 100
    bsexchsp = ds["pt2"][parm["setbasesatsp"]]["exchangespecies"] if "setbasesatsp" in parm else "NaX"
    cell["exchange"] = {bsexchsp: cell["X"] * bsfrac, exsp: cell["X"] * (1 - bsfrac)}
