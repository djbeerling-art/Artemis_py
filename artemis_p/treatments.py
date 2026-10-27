"""Fertiliser and other timetable treatments. Python counterpart of getfert.m
and of the processtimetable nested function in setupdatastruc.m.

``process_timetable`` matches timetable rows to run days, converts doses to
mol/m2 land and builds ``ds["ferttable"]``. The PHREEQC RATES include files
that MATLAB writes with ``sed`` as a side effect are written separately by
``write_treatment_rates``, so the preprocessor itself touches no files.
"""

import re
from pathlib import Path

import numpy as np

DOSETYPES = ("Pellets", "liquid", "phase")


def get_fert(treatment, units="g/m2land", convmassarea=False, verbose=False):
    """getfert: name, molar mass, PHREEQC formula and dose scale for a treatment.

    The scale converts a dose in ``units`` (per element, e.g. kgN/ha) into
    mass of the whole treatment.
    """
    units2gpm2, convunits = 1.0, units
    if convmassarea:
        if "kg" in units:
            units2gpm2 *= 1000
        if "ha" in units:
            units2gpm2 /= 1e4
        convunits = convunits.replace("kg", "g").replace("ha", "m2")
    if "uan" in treatment.lower():
        longname = "Urea Ammonium Nitrate"
        pct_n = float(re.match(r"UAN([-+\d.eE]+)", treatment).group(1))
        pct_urea = pct_n / 2 * 60.0556 / 2 / 14.0067
        mol_urea = pct_urea / 100 / 60.0556
        pct_nh4 = pct_n / 4 * 18.0385 / 14.0067
        mol_nh4 = pct_nh4 / 100 / 18.0385
        pct_no3 = pct_n / 4 * 62.0049 / 14.0067
        mol_no3 = pct_no3 / 100 / 62.0049
        pct_h2o = max(0, 100 - pct_urea - pct_nh4 - pct_no3)
        mol_h2o = pct_h2o / 100 / 18.0153
        molmass = (pct_nh4 * 18.0385 + pct_no3 * 62.0049 + pct_urea * 60.0556 + pct_h2o * 18.0153) / 100
        if "N" in units and "UAN" not in units:
            units2gpm2 = units2gpm2 / pct_n * 100
            convunits = convunits.replace("N", treatment)
        formula = "Urea%.2f (AmmH)%.2f (NO3)%.2f (H2O)%.2f" % tuple(
            m * molmass for m in (mol_urea, mol_nh4, mol_no3, mol_h2o))
    elif treatment.lower() == "mop":
        longname, formula, molmass = "Muriate of potash", "K Cl", 74.5510
        if "K" in units:
            units2gpm2 = units2gpm2 / 39.0983 * molmass
        convunits = convunits.replace("K", treatment)
    elif treatment.lower() == "soa":
        longname, formula, molmass = "Sulphate of ammonia", "(AmmH)2 SO4", 132.1406
        if treatment not in units:
            for tag, mm in (("N", 14.0067), ("SO4", 96.0636), ("S", 32.0660)):
                if tag in units:
                    units2gpm2 = units2gpm2 / mm * molmass
                    convunits = convunits.replace(tag, treatment)
                    break
    elif treatment.lower() in ("dap", "map"):
        if treatment.lower() == "dap":
            longname, formula, molmass, nn = "Diammonium phosphate", "(AmmH)2 (HPO4)", 132.0562, 2
        else:
            longname, formula, molmass, nn = "Monoammonium phosphate", "(AmmH) (H2PO4)", 115.0256, 1
        if "N" in units:
            units2gpm2 = units2gpm2 / nn / 14.0067 * molmass
            convunits = convunits.replace("N", treatment)
        elif "P" in units:
            units2gpm2 = units2gpm2 / 30.9736 * molmass
            convunits = convunits.replace("P", treatment)
    else:
        # anything else is taken as a chemical formula. MATLAB uses StoichTools molweight; the
        # molar masses here come from ARTEMIS's periodictable.csv, which may differ from
        # StoichTools' in the last digits (not checked: no published run doses one in its first days)
        from .formula import parse_formula
        from .periodictable import get_periodic_table

        mm = get_periodic_table()[1]
        spp = parse_formula(treatment)
        longname = formula = treatment
        molmass = sum(float(mm[el]["molmass"]) * n for el, n in spp.items())
        for el, n in spp.items():
            if el in units:
                units2gpm2 = units2gpm2 / float(mm[el]["molmass"]) / n * molmass
                convunits = convunits.replace(el, treatment)
                break
        treatment = treatment + "_trmt"
    if verbose:
        print(f"(getfert) Treatment {treatment} ({longname}) molar mass {molmass:g} g/mol, formula {formula}")
        print(f"(getfert) {treatment} units {units} (= {units2gpm2:g} {convunits})")
    return treatment, molmass, formula, units2gpm2


def process_timetable(ds):
    """processtimetable: run indices (1-based), doses in mol/m2 and the fertiliser table."""
    tt, run, parm = ds["timetable"], ds["run"], ds["parm"]
    ntt = len(tt["decyear"])
    tt["runindex"] = np.zeros(ntt)
    tt["molpm2"] = np.zeros(ntt)
    tt.setdefault("molmass", np.zeros(ntt))
    tt.setdefault("formula", [""] * ntt)
    if "dose" in tt and "dosage" not in tt:
        tt["dosage"] = np.array(tt["dose"], float)
    nfsset = 0
    for j in range(ntt):
        w = np.flatnonzero((run["year"] == tt["year"][j]) & (run["month"] == tt["month"][j])
                           & (run["day"] == tt["day"][j]))
        if not w.size:
            continue
        tt["runindex"][j] = w[0] + 1
        if tt["doseform"][j].lower() == "feedstock":
            if "fsname" in parm:
                tt["treatment"][j] = parm["fsname"]
                nfsset += 1
            continue
        dose = tt["dosage"][j]
        if np.isnan(dose) or dose <= 0:
            continue
        du = tt["doseunits"][j]
        if "ha" in du:
            dose = dose / 1e4
        if "t" in du:
            dose = dose * 1e6
        elif "kg" in du:
            dose = dose * 1e3
        if "mol/m" in du:
            tt["molpm2"][j] = dose
            continue
        if np.isnan(tt["molmass"][j]):
            tt["molmass"][j] = 0
        if tt["molmass"][j] == 0:
            first = next(i for i in range(ntt) if tt["treatment"][i].lower() == tt["treatment"][j].lower()
                         and tt["doseunits"][i].lower() == du.lower() and (tt["runindex"][i] > 0 or i == j))
            name, mm, formula, scl = get_fert(tt["treatment"][j], du, verbose=(first == j))
            if mm > 0:
                tt["treatment"][j], tt["molmass"][j], tt["formula"][j] = name, mm, formula
                dose = dose * scl
        if tt["molmass"][j] == 0:
            tt["molmass"][j] = 1
        tt["molpm2"][j] = dose / tt["molmass"][j]
    if nfsset:
        print(f"(setupdatastruc/processtimetable) {nfsset} feedstock treatments set to {parm['fsname']}")
    nout = int(np.sum(tt["runindex"] == 0))
    if nout:
        print(f"(setupdatastruc/processtimetable) {nout} of {ntt} timetable entries are outside the run timeline")

    ds["treatments"] = []
    if "doseform" not in tt:
        return ds
    iall, formulae, molmass = [], [], []
    rates = ds.setdefault("treatmentrates", {})
    for dt in DOSETYPES:
        wdose = [j for j in range(ntt) if tt["doseform"][j].lower() == dt.lower() and tt["runindex"][j] > 0]
        if not wdose:
            continue
        if dt != "phase":
            ratesfile = f"{dt}.RATES" if "expID" not in ds else f"{ds['expID']}_{dt}.RATES"
            ds.setdefault("INCLUDE", []).append(ratesfile)
        # MATLAB unique: sorted, index of first occurrence
        utr = sorted({tt["treatment"][j] for j in wdose})
        first = [next(j for j in wdose if tt["treatment"][j] == u) for u in utr]
        iall += first
        formulae += [tt["formula"][j] for j in first]
        molmass += [tt["molmass"][j] for j in first]
        ds["treatments"] += utr
        if dt != "phase":
            rates[ratesfile] = (dt, utr)
    if not ds["treatments"]:
        return ds
    ds["ferttable"] = {"treatment": list(ds["treatments"]), "PHREEQCformula": formulae,
                       "molmass": np.array(molmass), "ttindex": np.array(iall, float) + 1}
    return ds


def treatment_rates_text(dbdir, dosetype, treatments):
    """Text of a treatment RATES include file, as processtimetable builds it with sed.

    The first treatment gets a copy of ``<dbdir>/<dosetype>.dat`` with every
    occurrence of the dose type replaced; later ones replace only the first
    occurrence on each line and drop the RATES keyword and comment lines.
    """
    src = (Path(dbdir) / f"{dosetype}.dat").read_text().splitlines(keepends=True)
    out = [line.replace(dosetype, treatments[0]) for line in src]
    for tr in treatments[1:]:
        for line in src:
            line = line.replace(dosetype, tr, 1)
            if line.startswith("RATES") or line.startswith("#"):
                continue
            out.append(line)
    return "".join(out)


def write_treatment_rates(ds, workdir, overwrite=None):
    """Write the RATES include files listed by process_timetable into ``workdir``.

    As in MATLAB, an existing file is kept unless ``parm.writetrRATESfiles`` is
    set (or ``overwrite`` is True), even if it was written for other treatments.
    """
    if overwrite is None:
        overwrite = bool(ds["parm"].get("writetrRATESfiles", 0))
    written = []
    for fname, (dosetype, treatments) in ds.get("treatmentrates", {}).items():
        path = Path(workdir) / fname
        if path.is_file() and not overwrite:
            continue
        path.write_text(treatment_rates_text(ds["dbdir"], dosetype, treatments))
        written.append(path)
    return written
