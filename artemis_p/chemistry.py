"""Initial solution chemistry. Python counterpart of readinitialchemistry and
setinitialchemistry in readinputfiles.m.

An initialChemistry CSV has a header row (``species, rainwater, ...``), a
units row, then one row per species. The influent (default column
``rainwater``) becomes ``ds["influent"]``. The soil solution (default column
``porewater``) becomes ``ds["layer"]["soilsoln"]``; if there is no such
column the influent is used instead, without its pH.
"""

import csv
from pathlib import Path

import numpy as np

ALLOWED_UNITS = ("mol/L", "mg/L")


def set_initial_chemistry(ds, species, conc, units, chargefld, nrows):
    """setinitialchemistry: concentrations per species, plus charge equivalents.

    ``nrows`` is the number of rows each species array gets (1 for the
    influent, the number of layers for the soil solution).
    """
    if units not in ALLOWED_UNITS:
        msg = f"Unacceptable units ({units}) for PHREEQC solutions, returning from setinitialchemistry"
        print(msg)
        return {"errormsg": msg, "units": units}
    ispscl = 1 / 1000 if "mg" in units else 1
    soln = {"msp": {}, "spas": {}, "eq": {}}
    for sp, c in zip(species, conc):
        spas, spch = "", sp
        if c < -7000:
            continue
        if sp.lower() == "nh4":
            sp, spas, spch = "Amm", " as AmmH+", "N"
        elif sp.lower() == "no3":
            sp, spas, spch = "N", " as NO3-", "N"
        elif sp.lower() == "so4":
            sp, spas, spch = "S", " as SO4-2", "S"
        if sp.lower() == "ph" and c > 0:
            soln["msp"][sp] = np.full(nrows, c) if nrows > 1 else c
            soln["spas"][sp] = ""
            continue
        if sp.lower() == "pe":
            soln["msp"]["pe"] = np.full(nrows, c) if nrows > 1 else c
            soln["spas"][sp] = ""
            continue
        if sp not in ds["solute"]:
            continue
        soln["msp"][sp] = np.full(nrows, c) if nrows > 1 else c
        soln["spas"][sp] = spas
        if spch in ds["pt2"]:
            scl = ispscl * ds["pt2"][spch]["charge"]
            if ispscl < 1:
                scl = scl / ds["pt2"][spch]["molmass"]
        else:
            scl = 0
        soln["eq"][sp] = c * scl
    soln["irrchemscale"] = 1.0
    soln["units"] = {"msp": units}
    soln["chargefld"] = chargefld or "Cl"
    return soln


def read_initial_chemistry(ds, fname, influentfld="rainwater", soilsolnfld="porewater",
                           influentchargefld="", soilsolnchargefld=""):
    """readinitialchemistry: fill ds["init"], ds["influent"] and ds["layer"]["soilsoln"].

    MATLAB note: this replaces the whole ``ds.layer.soilsoln`` struct, so a pH
    profile read earlier from the initlayers file (readinitlayerdata) is
    discarded. When there is no porewater column the influent pH is then
    removed too, so soil cells start without a measured pH.
    """
    with open(Path(fname), newline="", encoding="utf-8-sig") as f:
        rows = [r for r in csv.reader(f)]
    h = [c.replace(" ", "").lower() for c in rows[0]]
    units = rows[1]
    body = [r for r in rows[2:] if any(c.strip() for c in r)]
    ds.setdefault("inputfiles", []).append(str(fname))
    init = ds.setdefault("init", {})
    init.setdefault("units", {})
    for j, name in enumerate(h):
        vals = [(r + [""] * len(h))[j].strip() for r in body]
        try:
            init[name] = np.array([np.nan if v == "" else float(v) for v in vals])
        except ValueError:
            init[name] = vals
        init["units"][name] = units[j].strip() if j < len(units) else ""
    if influentfld in init:
        inf = set_initial_chemistry(ds, init["species"], init[influentfld],
                                    init["units"][influentfld], influentchargefld, 1)
        inf.setdefault("label", {})
        inf["label"]["msp"] = influentfld
        inf["label"]["influent"] = influentfld
        ds["influent"] = inf
    rwflag = False
    if soilsolnfld not in init:
        soilsolnfld, rwflag = influentfld, True
    if soilsolnfld in init:
        ss = set_initial_chemistry(ds, init["species"], init[soilsolnfld],
                                   init["units"][soilsolnfld], soilsolnchargefld, ds["nlayers"])
        ss.setdefault("label", {})["msp"] = soilsolnfld
        layer = ds.setdefault("layer", {})
        layer["soilsoln"] = ss
        layer.setdefault("source", {})["soilsoln"] = str(fname)
        if rwflag and "pH" in ss["msp"]:
            print("(readinputfiles/readinitialchemistry) Warning: Removing influent pH from soil solution")
            del ss["msp"]["pH"]
    return ds
