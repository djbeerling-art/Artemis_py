"""Implicitly dissolving phases. Python counterpart of implicitcalc.m.

Phases flagged ``implicit`` (e.g. apatite and titanite in the Blue Ridge
feedstock) are lumped into one pseudo-phase per feedstock,
``Implicit_<feedstock>``, that dissolves in step with the feedstock's
explicitly dissolving phases.
"""

import numpy as np
import pandas as pd


def _g(x):
    """MATLAB %g."""
    return "%g" % x


def implicit_calc(ds, quiet=False):
    """implicitcalc.m.

    Existing ``Implicit_*`` rows (and the ``impmolrat`` column) are removed
    first. The pseudo-phase is only built once ``ds["stoich"]`` exists, i.e.
    after setupphases; before that implicitcalc returns. The optional RATES
    file writer (``parm.writeimplicitcalcfile``) is not ported.
    """
    parm = ds["parm"]
    parm.setdefault("writeimplicitcalcfile", 0)
    parm.setdefault("writeimplicitcalcfile_eq", 0)
    t = ds["phasetable"]
    imp_rows = ["Implicit_" in str(p) for p in t["PHREEQCphase"]]
    if any(imp_rows):
        t = t[[not r for r in imp_rows]]
        t = t.drop(columns=["impmolrat"], errors="ignore")
        ds["phasetable"] = t
    wimpall = [str(p).lower() == "implicit" for p in t["PHREEQCphasetype"]]
    if not any(wimpall):
        return ds
    if "stoich" not in ds:
        if not quiet:
            print("(implicitcalc) Unable to proceed without ds.stoich and ds.run.treatment")
        return ds
    for fs in sorted(set(t.loc[wimpall, "feedstock"])):
        _setup_implicit_phase(ds, fs, quiet)
        _get_explicit_grams(ds, fs)
        if parm["writeimplicitcalcfile"]:
            raise NotImplementedError("writing implicit-phase RATES files is not ported yet")
    return ds


def _setup_implicit_phase(ds, feedstock, quiet):
    t = ds["phasetable"]
    impname = f"Implicit_{feedstock}"
    ptype = t["PHREEQCphasetype"].astype(str).str.lower()
    wimp = np.flatnonzero((ptype == "implicit").to_numpy()
                          & (t["feedstock"].astype(str).str.lower() == feedstock.lower()).to_numpy())
    frac, molmass = t["frac"].to_numpy(float), t["molmass"].to_numpy(float)
    impstoich = np.zeros(len(t))
    impstoich[wimp] = frac[wimp] / molmass[wimp]
    impstoich = impstoich / impstoich[wimp[0]]
    impformula = ""
    for i in wimp:
        impf = f"({t['formula'].iloc[i]}){_g(impstoich[i])} "
        impformula += impf
        if not quiet:
            print(f"{t['MATLABphase'].iloc[i]} ({i + 1}) impf={impf}")
    # New row: copy of the first implicit phase with text emptied and numbers zeroed
    row = {}
    for c in t.columns:
        v = t[c].iloc[wimp[0]]
        if pd.api.types.is_numeric_dtype(t[c]) and not pd.api.types.is_bool_dtype(t[c]):
            row[c] = 0.0
        elif isinstance(v, str) or v is None or pd.api.types.is_string_dtype(t[c]):
            row[c] = ""
        else:
            row[c] = v
    row.update(PHREEQCphase=impname, MATLABphase=impname, feedstock=feedstock, PHREEQCphasetype="kinetic",
               frac=frac[wimp].sum(), molmass=np.sum(molmass[wimp] * impstoich[wimp]),
               G=np.sum(t["G"].to_numpy(float)[wimp] * impstoich[wimp]), formula=impformula)
    attrs = dict(t.attrs)
    t = pd.concat([t, pd.DataFrame([row], index=[impname])])
    t.attrs = attrs
    ds["phasetable"] = t
    if not quiet:
        print(f"(implicitcalc/setupimplicitphase) {impname} combined formula: {impformula}")
    _implicit_stoich(ds, impname, wimp, impstoich)
    _implicit_molrat(ds, feedstock, impname)
    run = ds["run"]
    if "treatment" not in run:
        return
    i0 = wimp[0]
    m0 = run["treatment"][feedstock] * frac[i0] / molmass[i0]
    run["moles"][impname] = m0
    run["RSA"][impname] = m0 * 0
    ds["cell"]["moles"][impname] = m0[:, 0].copy()
    ds["cell"]["RSA"][impname] = m0[:, 0] * 0


def _implicit_stoich(ds, impname, wimp, impstoich):
    """getimplicitstoich: element counts of the pseudo-phase (per mole of the first implicit phase)."""
    names = ds["phasetable"]["MATLABphase"]
    st = ds["stoich"]
    st[impname] = dict(st[names.iloc[wimp[0]]])
    for i in wimp[1:]:
        for el, n in st[names.iloc[i]].items():
            st[impname][el] = st[impname].get(el, 0) + n * impstoich[i]


def _implicit_molrat(ds, feedstock, impname):
    """getimplicitmolrat: implicit moles per mole of each explicit phase of the feedstock."""
    t = ds["phasetable"]
    if "impmolrat" not in t:
        t["molfrac"] = t["frac"] / t["molmass"]
        t["impmolrat"] = t["frac"] * 0
    ptype = t["PHREEQCphasetype"].astype(str)
    wexp = (t["feedstock"].astype(str) == feedstock) & (ptype.str.contains("kin") | ptype.str.contains("eq"))
    t.loc[wexp, "impmolrat"] = t.at[impname, "molfrac"] / t.loc[wexp, "molfrac"]


def _get_explicit_grams(ds, feedstock):
    """getexplicitgrams: grams of explicitly dissolving feedstock phases per cell and day."""
    t, run = ds["phasetable"], ds["run"]
    impname = f"Implicit_{feedstock}"
    treat = run["treatment"][feedstock]
    grams = treat * 0
    ptype = t["PHREEQCphasetype"].astype(str).str.lower()
    fs = t["feedstock"].astype(str).str.lower() == feedstock.lower()
    notimp = t["PHREEQCphase"].astype(str).str.lower() != impname.lower()
    for kind in ("kinetic", "equilibrium"):
        for f in t.loc[(ptype == kind) & fs & notimp, "frac"]:
            grams = grams + treat * f
    run.setdefault("grams", {})[f"Explicit_{feedstock}"] = grams
