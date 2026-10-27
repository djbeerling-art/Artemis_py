"""End-of-run totals: gettotalfields.m, getbasesatpct.m, assessuptake.m and etaRH.m.

``gettotalfields(out, ds)`` adds ``out['total']`` (soil-column totals, leachate
fluxes, CDR fields and elemental release) and the extra per-cell fields that
MATLAB adds to ``out['soil']``. It runs after the daily loop, from
getphreeqcrun's ``cleanup``.

Soil arrays are (ncells, nrun); total and leachate time series are 1-D here
(MATLAB's (1, nrun) rows).

The MATLAB bugs that change results are reproduced unless ``ds['fixes']`` turns
their fixes on (see artemis_p.fixes and docs/PORTING_PLAN.md):
native-mineral release goes to the ``c<elt>release`` fields, where the feedstock
step overwrites it, so ``c<elt>natrelease`` stays zero; ``cCDR`` adds the soil
carbonate stock to the daily CDR; ``dt_*uptake`` adds the cumulative passive
uptake instead of the daily active uptake; ``ccateqCECrat`` for the soil cells
uses the first cell's release only.
"""

import numpy as np

from .fixes import enabled as fixed
from .mformat import sprintf

__all__ = ["gettotalfields", "getbasesatpct", "assessuptake", "etaRH"]


def etaRH(salinity=35.0, temperature=17.0, pCO2=400.0):
    """Marine CO2 sequestration efficiency (Renforth and Henderson 2017, Eq. 11)."""
    S, T, p = salinity, temperature, np.asarray(pCO2, float)
    return ((S * 10 ** -3.009 + 10 ** -1.519) * np.log(p) - S * 10 ** -2.100
            - T * p * (S * 10 ** -7.501 - 10 ** -5.598) - T * 10 ** -2.337 + 10 ** -0.102)


def _row(a):
    return np.ravel(np.asarray(a, float))


def getbasesatpct(out, ds):
    soil, total = out["soil"], out["total"]
    if "m_CaX2" not in soil:
        return out
    with np.errstate(divide="ignore", invalid="ignore"):
        soil["BSpct"] = (2 * (soil["m_CaX2"] + soil["m_MgX2"]) + soil["m_NaX"] + soil["m_KX"]) / soil["X"] * 100
    X = np.asarray(ds["cell"]["X"], float).reshape(-1)
    total["BSpct"] = np.sum(soil["BSpct"] * X[:, None], axis=0) / np.sum(X)
    out["units"]["BSpct"] = "percent"
    out["label"]["BSpct"] = "Base saturation"
    total["units"]["BSpct"] = "percent"
    total["label"]["BSpct"] = "Base saturation"
    return out


def gettotalfields(out, ds):
    md, soil, leach = out["metadata"], out["soil"], out["leachate"]
    total = out.setdefault("total", {})
    tunits, tlabel = total.setdefault("units", {}), total.setdefault("label", {})
    units, label = out.setdefault("units", {}), out.setdefault("label", {})
    run = ds["run"]
    ncells = int(ds["ncells"])

    istart, ilast = int(md["istart"]), int(md["ilast"])
    nstop = np.size(out["run"]["stopwatch"])
    if "istart0" in md:
        istart = int(md["istart0"])
    elif ilast - istart + 1 < nstop:
        istart = 1
        md["istart0"] = 1
    tzero = np.zeros(np.shape(leach["pH"])[-1])
    szero = np.zeros(np.shape(soil["pH"]))
    scl = ds.get("units", {}).get("scl", "")
    timelab = " day^{-1}" if ds["datastep"] == 1 else " month^{-1}"

    for f in ("decyear", "year", "doy", "month", "day", "exptday"):
        if f in run:
            total[f] = _row(leach[f])
    total["elapsedyear"] = total["decyear"] - total["decyear"][0]
    tlabel["elapsedyear"] = "Decimal year"
    soil["elapsedyear"] = total["elapsedyear"].copy()
    label["elapsedyear"] = "Decimal year"
    label["decyear"] = "Decimal year"

    H = np.asarray(md["H2Oscale"], float).reshape(-1)[:, None]
    mmscl = ds.get("mmscl", 1)
    tCO2m2 = (12.01 + 16 + 16) / 1e6 / mmscl
    tCm2 = 12.01 / 1e6 / mmscl
    CO2lab, Clab = "tCO_2 ha^{-1}", "tC ha^{-1}"
    tCO2ha, tCha = tCO2m2 * 1e4, tCm2 * 1e4
    md.update(tCO2ha=tCO2ha, tCO2m2=tCO2m2, tCha=tCha, tCm2=tCm2)
    md.setdefault("units", {}).update(tCO2ha=CO2lab, tCO2m2="tCO_2 m^{-2}land", tCm2="tC m^{-2}land", tCha=Clab)
    cfluxlab = "mol " + scl
    soilfld = list(soil)

    # ---- getfluxes
    total["mass_H2O"] = np.sum(soil["mass_H2O"], axis=0)
    tunits["mass_H2O"] = "kg water " + scl
    tlabel["mass_H2O"] = "Total soil water"
    outperc = np.asarray(run["outperc"], float)[ncells - 1]
    total["q"] = outperc.copy()
    tunits["q"] = "kg water " + scl + timelab
    total["qm2"] = outperc / mmscl
    tunits["qm2"] = "kg water m^{-2}land" + timelab
    total["qha"] = total["qm2"] * 1e4
    tunits["qha"] = "kg water ha^{-1}" + timelab
    lm = leach["mass_H2O"]
    nh2o = lm.size
    lm[..., outperc[:nh2o] < 1e-6] = 0
    solute = list(ds["solute"]) + [s for s in ("m_DIN", "m_DIC", "m_HCO3", "m_CO3", "cateq", "Urea") if s in leach]
    noflux = np.isnan(_row(leach["pH"]))
    for s in solute:
        if s not in leach:
            continue
        nom = s.replace("m_", "")
        ff, cff = nom + "flux", "c" + nom + "flux"
        flux = _row(leach[s]) * _row(lm) * mmscl
        total[ff] = flux.copy()
        tunits[ff] = "mol m^2_{land} " + timelab
        tlabel[ff] = nom + " export"
        flux[noflux] = 0
        total[cff] = np.cumsum(flux)
        tunits[cff] = "mol m^2_{land} "
        tlabel[cff] = "Cumulative " + nom + " export"

    # ---- gettoteqphases
    for d in soilfld:
        if "d_" not in d or "soilCarbonate" in d:
            continue
        e = d.replace("d_", "")
        total[e] = np.sum(soil[e] * H, axis=0)
        total[d] = np.sum(soil[d] * H, axis=0)
        tunits[e] = "mol " + scl
        tlabel[e] = "Total " + e
        tunits[d] = tunits[e] + " " + timelab
        tlabel[d] = d.replace("d_", "\\Delta")

    # ---- gettotkinphases
    for d in soilfld:
        if "dk_" not in d:
            continue
        k = d.replace("dk_", "k_")
        total[k] = np.sum(soil[k] * H, axis=0)
        total[d] = np.sum(soil[d] * H, axis=0)
        tunits[k] = "mol " + scl
        tlabel[k] = "Total " + k
        tunits[d] = tunits[k] + " " + timelab
        tlabel[d] = d.replace("d_", "\\Delta")
        if "Nuptake" in k:
            tlabel[k] = tlabel[k].replace("Nuptake", "NO_3^- uptake")
            tlabel[d] = tlabel[d].replace("Nuptake", "NO_3^- uptake")

    # ---- gettotupphases
    if ds["ETmethod"].lower() == "mixflownut":
        try:
            for dp in soilfld:
                if "dp_" not in dp:
                    continue
                p = dp.replace("dp_", "p_")
                plab = p.replace("p_", "").replace("uptake", " uptake")
                dt = dp.replace("dp_", "dt_")
                t = p.replace("p_", "t_")
                k = p.replace("p_", "k_")
                # MATLAB bug: p has no 'dp_', so this is the passive p_ field
                dk = dp.replace("dp_", "dk_") if fixed(ds, "uptake_totals") else p.replace("dp_", "dk_")
                total[p] = np.sum(soil[p] * H, axis=0)
                total[dp] = np.sum(soil[dp] * H, axis=0)
                tunits[p] = "mol " + scl
                tlabel[p] = "Passive " + plab
                tunits[dp] = tunits[p] + " " + timelab
                tlabel[dp] = "Daily passive " + plab
                if k not in soil:
                    continue
                total[t] = total[p] + total[k]
                total[dt] = total[dp] + total[dk]
                tunits[t] = "mol " + scl
                tlabel[t] = "Passive+active " + plab
                tunits[dt] = tunits[t] + " " + timelab
                tlabel[dt] = "Daily passive+active " + plab
        except KeyError as e:  # MATLAB catches and prints the error
            print(f"(gettotalfields/gettotupphases) {e}")

    getbasesatpct(out, ds)

    # ---- gettotsorbed
    for x in soilfld:
        if "X" in x:
            total[x] = np.sum(soil[x] * soil["mass_H2O"], axis=0) * mmscl
            tunits[x] = "mol " + scl
            tlabel[x] = "Total " + x.replace("m_", "")
    for s in soilfld:
        if "surf" in s:
            total[s] = np.sum(soil[s] * H, axis=0)
            tunits[s] = "mol " + scl
            tlabel[s] = "Total " + s
    for s in ("EDL", "WHAM", "HFO", "CEC"):
        if s in soil:
            total[s] = np.sum(soil[s] * H, axis=0)
            tunits[s] = "mol " + scl
            tlabel[s] = s

    # ---- gettotsys
    for s in soilfld:
        if "sys" in s:
            total[s] = np.sum(soil[s] * H, axis=0)
            tunits[s] = "mol " + scl
            tlabel[s] = "Total " + s

    _gettotalCDR(out, ds, istart, ilast, tzero, szero, H, tCO2ha, CO2lab, timelab)
    _getnatsoilrelease(out, ds, ilast, tzero, szero, tCO2ha, CO2lab, cfluxlab)
    _getfsrelease(out, ds, ilast, tzero, szero, tCO2ha, CO2lab, cfluxlab)
    assessuptake(out, ds)
    return out


def _gettotalCDR(out, ds, istart, ilast, tzero, szero, H, tCO2ha, CO2lab, timelab):
    md, soil, total = out["metadata"], out["soil"], out["total"]
    tunits, tlabel = total["units"], total["label"]
    units, label = out["units"], out["label"]
    run = ds["run"]

    total["CO2cons"] = (total["HCO3flux"] + total["CO3flux"]) * tCO2ha
    tunits["CO2cons"] = CO2lab + timelab
    tlabel["CO2cons"] = "CO_2 consumption"
    cC = np.nan_to_num(total["CO2cons"], nan=0.0, posinf=np.inf, neginf=-np.inf)
    total["cCO2cons"] = np.cumsum(cC)
    tunits["cCO2cons"] = CO2lab
    tlabel["cCO2cons"] = "Cumulative CO_2 consumption"

    parm = ds.setdefault("parm", {})
    parm.setdefault("ocean_salinity", 35)
    parm.setdefault("ocean_tmp", 17)
    nday = np.shape(soil["pH"])[1]
    if ilast - istart + 1 > nday:
        ilast = nday
        md["ilast"] = nday
    eta = etaRH(parm["ocean_salinity"], parm["ocean_tmp"], _row(run["CO2atm"])[:ilast])
    total["CDRrh"] = total["cateqflux"] * eta * tCO2ha
    tunits["CDRrh"] = CO2lab + timelab
    tlabel["CDRrh"] = sprintf("Marine CDR {\\eta}=%5.3f", np.mean(eta))
    total["cCDRrh"] = total["ccateqflux"] * eta * tCO2ha
    tunits["cCDRrh"] = CO2lab
    tlabel["cCDRrh"] = sprintf("Cumulative marine CDR {\\eta}=%5.3f", np.mean(eta))
    md["etaRH"] = eta
    md.setdefault("label", {})["etaRH"] = "{\\\\eta} (Renforth and Henderson 2017)"
    total["CDR"] = total["CDRrh"].copy()
    tunits["CDR"] = CO2lab + timelab
    tlabel["CDR"] = "Total soil+marine CDR"
    total["cCDR"] = total["cCDRrh"].copy()
    tunits["cCDR"] = CO2lab
    tlabel["cCDR"] = "Total cumulative soil+marine CDR"

    sinks = []
    pt = ds.get("phasetable")
    if pt is not None:
        w = pt["formula"].astype(str).str.contains("CO3") & ~pt["PHREEQCphasetype"].astype(str).str.contains("ignore")
        for name, ptype in zip(pt.loc[w, "MATLABphase"], pt.loc[w, "PHREEQCphasetype"]):
            sinks.append("k_" + name if "kin" in ptype else name)

    total["soilCarbonate"] = tzero.copy()
    total["d_soilCarbonate"] = tzero.copy()
    total["cd_soilCarbonate"] = tzero.copy()
    tlabel["soilCarbonate"] = "Soil carbonates"
    tlabel["d_soilCarbonate"] = "Change in soil carbonates"
    tlabel["cd_soilCarbonate"] = "Cumulative change in soil carbonates"
    tunits["soilCarbonate"] = CO2lab
    tunits["d_soilCarbonate"] = CO2lab + " " + timelab
    tunits["cd_soilCarbonate"] = CO2lab
    soil["soilCarbonate"] = szero.copy()
    soil["d_soilCarbonate"] = szero.copy()
    soil["cd_soilCarbonate"] = szero.copy()
    label["soilCarbonate"] = "Soil carbonates"
    label["d_soilCarbonate"] = "Change in soil carbonates"
    label["cd_soilCarbonate"] = "Cumulative change in soil carbonates"
    units["soilCarbonate"] = "tCO_2 cell^{-1}"
    units["d_soilCarbonate"] = "tCO_2 cell^{-1} day^{-1}"
    units["cd_soilCarbonate"] = "tCO_2 cell^{-1}"
    for co3 in sinks:
        cdr = co3 + "CDR"
        if "k_" in co3:
            dfld, dcdr = "d" + co3, "d" + cdr
        else:
            dfld, dcdr = "d_" + co3, "d_" + co3 + "CDR"
        soil[cdr] = soil[co3] * H * tCO2ha
        soil[dcdr] = soil[dfld] * H * tCO2ha
        soil["soilCarbonate"] = soil["soilCarbonate"] + soil[cdr]
        soil["d_soilCarbonate"] = soil["d_soilCarbonate"] + soil[dcdr]
        total[cdr] = np.sum(soil[co3][:, istart - 1:ilast], axis=0) * tCO2ha
        tunits[cdr] = CO2lab
        tlabel[cdr] = "Soil " + co3
        total[dcdr] = total[dfld] * tCO2ha
        tunits[dcdr] = CO2lab + " " + timelab
        tlabel[dcdr] = tlabel[dfld]
    total["soilCarbonate"] = np.sum(soil["soilCarbonate"], axis=0)
    total["d_soilCarbonate"] = np.sum(soil["d_soilCarbonate"], axis=0)
    soil["cd_soilCarbonate"] = np.cumsum(soil["d_soilCarbonate"], axis=1)
    total["cd_soilCarbonate"] = np.sum(soil["cd_soilCarbonate"], axis=0)
    total["CDR"] = total["CDR"] + total["d_soilCarbonate"]
    if fixed(ds, "cdr_total"):
        total["cCDR"] = total["cCDRrh"] + total["cd_soilCarbonate"]
    else:  # MATLAB bug: adds the daily CDR to the soil carbonate stock, not the cumulative change
        total["cCDR"] = total["CDR"] + total["soilCarbonate"]

    total["CDRDIC"] = total["DICflux"] * tCO2ha
    tunits["CDRDIC"] = CO2lab + timelab
    tlabel["CDRDIC"] = "Exported DIC"
    total["cCDRDIC"] = total["cDICflux"] * tCO2ha
    tunits["cCDRDIC"] = CO2lab
    tlabel["cCDRDIC"] = "Exported DIC"
    total["CO2cateq"] = total["cateqflux"] * tCO2ha
    tunits["CO2cateq"] = CO2lab + timelab
    tlabel["cCO2cateq"] = "Cationic CO_2 consumption"
    total["cCO2cateq"] = total["ccateqflux"] * tCO2ha
    tunits["cCO2cateq"] = CO2lab
    tlabel["cCO2cateq"] = "Cumulative cationic CO_2 consumption"
    ys = total.setdefault("yalescale", {})
    ys.update(cCDRrh=eta * 2, cCDRcateq=eta * 2, cCDRDIC=eta * 2, cCO2cons=2, CalciteCDR=1, cCDR=1)


def _phases(pt, mask):
    return list(pt.loc[mask, "MATLABphase"])


def _getnatsoilrelease(out, ds, ilast, tzero, szero, tCO2ha, CO2lab, cfluxlab):
    pt = ds.get("phasetable")
    if pt is None or "stoich" not in ds:
        return
    soil, total = out["soil"], out["total"]
    tunits, tlabel = total["units"], total["label"]
    units, label = out["units"], out["label"]
    cat = pt["category"].astype(str).str.lower()
    flds = _phases(pt, cat.isin(["soil", "native", "parent"]))
    elt0 = ["Ca", "Mg", "Na", "K"]
    cfld = [f"c{e}natrelease" for e in elt0]
    for e, c in zip(elt0, cfld):
        soil[c] = szero.copy()
        total[c] = tzero.copy()
        units[c] = cfluxlab
        label[c] = e + " from native soil phases"
        tunits[c] = cfluxlab
        tlabel[c] = e + " from native soil phases"
    moles = ds["run"]["moles"]
    for p in flds:
        if str(pt.loc[p, "PHREEQCphasetype"]).lower() != "kinetic":
            print(f"(gettotalfields/getfsrelease) {p} is not a kinetic phase so will be excluded from native release arrays")
            continue
        kf = "k_" + p
        for e, st in ds["stoich"][p].items():
            if e.lower() in ("h", "o"):
                continue
            # MATLAB bug: looks for c<elt>release in a list of c<elt>natrelease names, so
            # native release goes to the feedstock release fields (and c<elt>natrelease stays 0)
            cf = f"c{e}natrelease" if fixed(ds, "native_release") else f"c{e}release"
            if cf.lower() not in (c.lower() for c in cfld):
                cfld.append(cf)
                fs = e + "native"
                soil[cf] = szero.copy()
                total[cf] = tzero.copy()
                soil[fs] = szero.copy()
                total[fs] = tzero.copy()
                units[cf] = cfluxlab
                label[cf] = e + " from native phases"
                tunits[cf] = cfluxlab
                tlabel[cf] = e + " from native phases"
            diff = np.cumsum(np.asarray(moles[p], float)[:, :ilast], axis=1) - soil[kf]
            soil[cf] = soil[cf] + diff * st
            total[cf] = total[cf] + np.sum(diff, axis=0) * st

    soil["cCaMgNanatrelease"] = (soil["cCanatrelease"] * 2 + soil["cMgnatrelease"] * 2 + soil["cNanatrelease"]) * tCO2ha
    total["cCaMgNanatrelease"] = (total["cCanatrelease"] * 2 + total["cMgnatrelease"] * 2
                                  + total["cNanatrelease"]) * tCO2ha
    ys = total.setdefault("yalescale", {})
    ys["cCaMgNanatrelease"] = 2
    tunits["cCaMgNanatrelease"] = CO2lab
    tlabel["cCaMgNanatrelease"] = "Ca+Mg+Na from native phases"
    label["cCaMgNanatrelease"] = "Ca+Mg+Na from native phases"
    total["ccateqnatrelease"] = total["cCaMgNanatrelease"] + total["cKnatrelease"] * tCO2ha
    soil["ccateqnatrelease"] = soil["cCaMgNanatrelease"] + soil["cKnatrelease"] * tCO2ha
    tunits["ccateqnatrelease"] = CO2lab
    tlabel["ccateqnatrelease"] = "Ca+Mg+Na+K released from native phases"
    ys["ccateqnatrelease"] = 2
    eqCEC = np.asarray(ds["cell"]["eqCEC"], float).reshape(-1)
    total["ccateqnatCECrat"] = total["ccateqnatrelease"] / tCO2ha / np.sum(eqCEC)
    # MATLAB uses the first cell's release for every cell
    rel = soil["ccateqnatrelease"] if fixed(ds, "cec_ratio_cells") else soil["ccateqnatrelease"][0][None, :]
    soil["ccateqnatCECrat"] = rel / tCO2ha / eqCEC[:, None]
    units["ccateqnatCECrat"] = "eq eq^{-1}CEC"
    label["ccateqnatCECrat"] = "Ca+Mg+Na+K released from native phases"
    tunits["ccateqnatCECrat"] = "eq eq^{-1}CEC"
    tlabel["ccateqnatCECrat"] = "Ca+Mg+Na+K released from native phases"


def _getfsrelease(out, ds, ilast, tzero, szero, tCO2ha, CO2lab, cfluxlab):
    soil, total = out["soil"], out["total"]
    if "cCarelease" not in total and "dust" not in out.get("runinfo", {}).get("including", ""):
        return
    if "feedstocktable" not in ds or "stoich" not in ds:
        return
    tunits, tlabel = total["units"], total["label"]
    units, label = out["units"], out["label"]
    pt = ds["phasetable"]
    fs = pt["feedstock"].astype(str).str.lower()
    flds = _phases(pt, ~fs.isin(["", "soil", "native", "parent"]))
    elt0 = ["Ca", "Mg", "Na", "K"]
    cfld = [f"c{e}release" for e in elt0]
    fsfld = [f"{e}feedstock" for e in elt0]
    for e, c, f in zip(elt0, cfld, fsfld):
        soil[c] = szero.copy()
        total[c] = tzero.copy()
        soil[f] = szero.copy()
        total[f] = tzero.copy()
        units[c] = cfluxlab
        label[c] = e + " from feedstock"
        tunits[c] = cfluxlab
        tlabel[c] = e + " from feedstock"
        tunits[f] = cfluxlab
        tlabel[f] = e + " in remaining feedstock"
    M, M0 = out["run"]["M"], out["run"]["M0"]
    for p in flds:
        ptype = str(pt.loc[p, "PHREEQCphasetype"]).lower()
        if ptype in ("ignore", "implicit"):
            continue
        if ptype != "kinetic":
            print(f"(gettotalfields/getfsrelease) {p} is not a kinetic phase so will be excluded from feedstock release arrays")
            continue
        kf = "k_" + p
        klast = soil[kf].shape[1]
        wz = np.flatnonzero(M[p][0, :klast] != soil[kf][0, :klast])
        if wz.size:
            M[p][:, wz] = soil[kf][:, wz]
        if pt.loc[p, "reversible"]:
            diff = np.cumsum(M0[p][:, :ilast] - M[p][:, :ilast], axis=1)
        else:
            diff = np.cumsum(np.maximum(0, M0[p][:, :ilast] - M[p][:, :ilast]), axis=1)
        for e, st in ds["stoich"][p].items():
            if e.lower() in ("h", "o"):
                continue
            cf = f"c{e}release"
            low = [c.lower() for c in cfld]
            if cf.lower() in low:
                k = low.index(cf.lower())
            else:
                k = len(cfld)
                cfld.append(cf)
                fsfld.append(e + "feedstock")
                soil[cf] = szero.copy()
                total[cf] = tzero.copy()
                soil[fsfld[k]] = szero.copy()
                total[fsfld[k]] = tzero.copy()
                units[cf] = cfluxlab
                label[cf] = e + " from feedstock"
                tunits[cf] = cfluxlab
                tlabel[cf] = e + " from feedstock"
                tunits[fsfld[k]] = cfluxlab
                tlabel[fsfld[k]] = e + " in remaining feedstock"
            c, f = cfld[k], fsfld[k]
            soil[c] = soil[c] + diff * st
            total[c] = total[c] + np.sum(diff, axis=0) * st
            soil[f] = soil[f] + soil[kf] * st
            if fixed(ds, "feedstock_remaining"):
                total[f] = total[f] + np.sum(soil[kf] * st, axis=0)
            else:  # MATLAB adds the running cell sum to the total each time
                total[f] = total[f] + np.sum(soil[f], axis=0)

    # MATLAB converts these totals to tCO2/ha here and again in cCaMgNarelease below
    for e, fac in (("Ca", 2), ("Mg", 2), ("Na", 1), ("K", 1)):
        c = f"c{e}release"
        total[c] = total[c] * fac * tCO2ha
        tlabel[c] = e + " from feedstock"
        tunits[c] = CO2lab
    soil["cCaMgNarelease"] = (soil["cCarelease"] * 2 + soil["cMgrelease"] * 2 + soil["cNarelease"]) * tCO2ha
    if fixed(ds, "release_units"):  # the totals are already in tCO2/ha
        total["cCaMgNarelease"] = total["cCarelease"] + total["cMgrelease"] + total["cNarelease"]
    else:
        total["cCaMgNarelease"] = (total["cCarelease"] * 2 + total["cMgrelease"] * 2 + total["cNarelease"]) * tCO2ha
    ys = total.setdefault("yalescale", {})
    ys["cCaMgNarelease"] = 2
    tunits["cCaMgNarelease"] = CO2lab
    tlabel["cCaMgNarelease"] = "Ca+Mg+Na from feedstock"
    label["cCaMgNarelease"] = "Ca+Mg+Na from feedstock"
    total["ccateqrelease"] = total["cCaMgNarelease"] + total["cKrelease"] * (1 if fixed(ds, "release_units") else tCO2ha)
    soil["ccateqrelease"] = soil["cCaMgNarelease"] + soil["cKrelease"] * tCO2ha
    tunits["ccateqrelease"] = CO2lab
    tlabel["ccateqrelease"] = "Ca+Mg+Na+K released from feedstock"
    ys["ccateqrelease"] = 2
    eqCEC = np.asarray(ds["cell"]["eqCEC"], float).reshape(-1)
    total["ccateqCECrat"] = total["ccateqrelease"] / tCO2ha / np.sum(eqCEC)
    # MATLAB uses the first cell's release for every cell
    rel = soil["ccateqrelease"] if fixed(ds, "cec_ratio_cells") else soil["ccateqrelease"][0][None, :]
    soil["ccateqCECrat"] = rel / tCO2ha / eqCEC[:, None]
    units["ccateqCECrat"] = "eq eq^{-1}CEC"
    label["ccateqCECrat"] = "Ca+Mg+Na+K released from feedstock"
    tunits["ccateqCECrat"] = "eq eq^{-1}CEC"
    tlabel["ccateqCECrat"] = "Ca+Mg+Na+K released from feedstock"


def assessuptake(out, ds):
    """assessuptake.m: daily uptake from the cumulative k_ fields, uptake efficiency and annual values."""
    md, soil, total = out["metadata"], out["soil"], out["total"]
    istart, ilast = int(md["istart"]), int(md["ilast"])
    if "istart0" in md:
        istart = int(md["istart0"])
    nrun = ilast - istart + 1
    run = ds["run"]
    plant = _row(run["plant"])[istart - 1:ilast]
    ncells = int(ds["ncells"])
    if "k_Ammuptake" not in soil:
        return out

    def daily(kup):
        d = kup - np.roll(kup, 1, axis=1)
        d[:, 0] = kup[:, 0]
        return d * plant

    w = wds = None

    def efficiency(upnut, dsnut):
        kup, eff = f"k_{upnut}uptake", f"eff_{upnut}uptake"
        soil[eff] = np.zeros((ncells, nrun))
        with np.errstate(divide="ignore", invalid="ignore"):  # MATLAB gives NaN/Inf where demand is 0
            soil[eff][:, w] = soil[kup][:, w] / np.asarray(run["cellcuptake"][dsnut], float)[:, wds]
        out["units"][eff] = "fraction of demand"
        out["label"][eff] = f"Efficiency of {upnut} uptake"
        total[eff] = np.zeros(nrun)
        total[eff][w] = np.sum(soil[kup][:, w], axis=0) / _row(run["totcuptake"][dsnut])[wds]
        total["units"][eff] = "fraction of demand"
        total["label"][eff] = f"Efficiency of {upnut} uptake"

    def annual(upnut):
        ann = out.setdefault("annual", {})
        if "year" not in ann:
            ann["year"] = np.unique(total["year"])
        eff, up = f"eff_{upnut}uptake", f"k_{upnut}uptake"
        ann[up] = np.zeros(len(ann["year"]))
        ann[eff] = np.zeros(len(ann["year"]))
        for yi, y in enumerate(ann["year"]):
            wy = np.flatnonzero((total["year"] == y) & (plant > 0))
            if not wy.size:
                continue
            ann[up][yi] = total[up][wy[-1]]
            ann[eff][yi] = total[eff][wy[-1]]
            ann.setdefault("units", {})[up] = total["units"][up] + " yr^{-1}"
            ann["units"][eff] = total["units"][eff] + " yr^{-1}"
            ann.setdefault("label", {})[eff] = f"Efficiency of {upnut} uptake"
            ann["label"][up] = total["label"][up]

    upeq = soil["k_Ammuptake"].copy()
    dupeq = daily(soil["k_Ammuptake"])
    soil["dk_Ammuptake"] = daily(soil["k_Ammuptake"])
    soil["k_DINuptake"] = soil["k_Ammuptake"] + soil["k_Nuptake"]
    total["k_DINuptake"] = total["k_Ammuptake"] + total["k_Nuptake"]
    total["label"]["k_DINuptake"] = "Total DIN uptake"
    soil.setdefault("label", {})["k_DINuptake"] = "DIN uptake"
    soil["dk_DINuptake"] = daily(soil["k_DINuptake"])
    soil["DINsys"] = soil["Ammsys"] + soil["Nsys"]
    total["DINsys"] = total["Ammsys"] + total["Nsys"]
    for n in ds["nut"]:
        up = f"k_{n}uptake"
        if up not in soil:
            continue
        upi = soil[up]
        dupi = daily(upi)
        soil[f"dk_{n}uptake"] = dupi
        w = np.flatnonzero(plant > 0)
        wds = w + istart - 1
        efficiency(n, n)
        annual(n)
        chg = ds["pt2"][n]["charge"]
        dupeq = dupeq + dupi * chg
        upeq = upeq + upi * chg
    efficiency("Amm", "N")
    annual("Amm")
    efficiency("DIN", "N")
    annual("DIN")
    soil["upeq"] = upeq
    out["units"]["upeq"] = "eq cell^{-1}"
    out["label"]["upeq"] = "Cumulative uptake"
    soil["dupeq"] = dupeq
    out["units"]["dupeq"] = "eq cell^{-1} day^{-1}"
    out["label"]["dupeq"] = "Uptake"
    total["upeq"] = np.sum(upeq, axis=0)
    total["dupeq"] = np.sum(dupeq, axis=0)
    total["units"]["upeq"] = "eq"
    total["units"]["dupeq"] = "eq day^{-1}"
    return out
