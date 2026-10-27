"""Layer-wise CO2 consumption and net CDR. Python counterparts of getlayerfluxes.m
and getnetCDRfluxes.m.

The manuscript's CDR metrics are the bicarbonate export (``cCO2cons``) and
cation export (``cCO2cateq``) at the bottom of layer 5 (76 cm), net of the
control run: ``soil["netcCO2cons"][4]`` and ``soil["netcCO2cateq"][4]`` after
``get_net_cdr_fluxes(treated, control, ds_treated, ds_control)``.
"""

import numpy as np

from .fixes import enabled as fixed

_ADDSOLUTES = ("m_DIN", "m_DIC", "m_HCO3", "m_CO3", "cateq", "Urea")


def _cumsum_nonan(x):
    x = np.array(x, float)
    x[np.isnan(x)] = 0
    return np.cumsum(x, axis=-1)


def get_layer_fluxes(out, ds):
    """getlayerfluxes: per-layer daily and cumulative CO2 consumption (HCO3+CO3) and
    cationic CO2 consumption (2Ca+2Mg+Na+K) in t CO2/ha, and solute exports in
    mol/m2 land, all at the bottom of each cell. Updates and returns ``out``."""
    soil, units, label = out["soil"], out.setdefault("units", {}), out.setdefault("label", {})
    md = out["metadata"]
    tCO2ha, ilast = md["tCO2ha"], int(md["ilast"])
    outperc = np.asarray(ds["run"]["outperc"], float)[:, :ilast]

    def a(f):
        return np.asarray(soil[f], float)

    soil["CO2cons"] = (a("m_HCO3") + a("m_CO3")) * outperc * tCO2ha
    units["CO2cons"], label["CO2cons"] = "t CO_2 ha^{-1} day^{-1}", "CO_2 consumption"
    soil["cCO2cons"] = _cumsum_nonan(soil["CO2cons"])
    units["cCO2cons"], label["cCO2cons"] = "t CO_2 ha^{-1}", "Cumulative CO_2 consumption"
    soil["CO2cateq"] = (a("Ca") * 2 + a("Mg") * 2 + a("Na") + a("K")) * outperc * tCO2ha
    units["CO2cateq"], label["CO2cateq"] = "t CO_2 ha^{-1} day^{-1}", "Cationic CO_2 consumption"
    soil["cCO2cateq"] = _cumsum_nonan(soil["CO2cateq"])
    units["cCO2cateq"], label["cCO2cateq"] = "t CO_2 ha^{-1}", "Cumulative cationic CO_2 consumption"

    solutes = list(ds["solute"]) + [s for s in _ADDSOLUTES if s in soil]
    nanpH = np.isnan(a("pH"))
    mmscl = ds.get("mmscl", 1)
    for s in solutes:
        if s not in soil:
            continue
        nom = s.replace("m_", "")
        fluxfld = nom + "flux"
        conc = a(s)
        # MATLAB uses the top cell's concentration for every layer: run.soil.(solute)(1,:)
        if not fixed(ds, "layer_solute_fluxes"):
            conc = conc[0]
        flux = conc * outperc * mmscl
        soil[fluxfld] = flux
        units[fluxfld], label[fluxfld] = "mol m^2_{land} day^{-1}", f"{nom} export"
        flux = flux.copy()
        flux[nanpH] = 0
        soil["c" + fluxfld] = np.cumsum(flux, axis=-1)
        units["c" + fluxfld], label["c" + fluxfld] = "mol m^2_{land} ", f"Cumulative {nom} export"
    return out


def get_net_cdr_fluxes(tr, ctl, ds_tr, ds_ctl):
    """getnetCDRfluxes: the treated run's CO2 consumption and cation export net of the
    control run, per layer and in total, plus CDR efficiencies (Kanzaki et al. 2025).
    Updates and returns ``tr``."""
    mt, mc = tr["metadata"], ctl["metadata"]
    istart = max(int(mt.get("istart0", mt["istart"])), int(mc.get("istart0", mc["istart"])))
    ilast = min(int(mt["ilast"]), int(mc["ilast"]))
    w = slice(istart - 1, ilast)
    if "cCO2cons" not in tr["soil"]:
        get_layer_fluxes(tr, ds_tr)
    if "cCO2cons" not in ctl["soil"]:
        get_layer_fluxes(ctl, ds_ctl)
    soil, total, units, label = tr["soil"], tr["total"], tr["units"], tr["label"]
    tunits, tlabel = total.setdefault("units", {}), total.setdefault("label", {})
    cs, ct = ctl["soil"], ctl["total"]

    def t(d, f):
        return np.asarray(d[f], float)

    label["cCO2cons"], label["cCO2cateq"] = "HCO_3^- export", "Cation export"
    total["netcCO2cons"] = t(total, "cCO2cons")[..., w] - t(ct, "cCO2cons")[..., w]
    total["netcCO2cateq"] = t(total, "cCO2cateq")[..., w] - t(ct, "cCO2cateq")[..., w]
    tr["leachate"]["netcCO2cons"] = total["netcCO2cons"]
    tr["leachate"]["netcCO2cateq"] = total["netcCO2cateq"]
    soil["netcCO2cons"] = t(soil, "cCO2cons")[:, w] - t(cs, "cCO2cons")[:, w]
    soil["netcCO2cateq"] = t(soil, "cCO2cateq")[:, w] - t(cs, "cCO2cateq")[:, w]
    soil["deltapH"] = t(soil, "pH")[:, w] - t(cs, "pH")[:, w]
    label.update(netcCO2cons="Net HCO_3^- export", netcCO2cateq="Net cation export", deltapH="{\\Delta}pH")
    units.update(netcCO2cons="t CO_2 ha^{-1}", netcCO2cateq="t CO_2 ha^{-1}")

    total["CDRundiss"] = ((t(total, "Cafeedstock") + t(total, "Mgfeedstock")) * 2 + t(total, "Nafeedstock")
                          + t(total, "Kfeedstock")) * mt["tCO2ha"]
    tlabel["CDRundiss"], tunits["CDRundiss"] = "Undissolved feedstock", "t CO_2 ha^{-1}"
    added = t(total, "ccateqrelease") + total["CDRundiss"]
    with np.errstate(divide="ignore", invalid="ignore"):
        total["CDReffdiss"] = np.ravel(total["ccateqrelease"]) / np.ravel(added) * 100
        soil["CDReffcat"] = soil["netcCO2cateq"] / np.ravel(added) * 100
        soil["CDReffcons"] = soil["netcCO2cons"] / np.ravel(added) * 100
    tlabel["CDReffdiss"], tunits["CDReffdiss"] = "Dissolution efficiency", "%"
    label.update(CDReffcat="Cationic export efficiency", CDReffcons="HCO_3^- export efficiency")
    units.update(CDReffcat="%", CDReffcons="%")
    total["cCaMgrelease"] = t(total, "cCarelease") + t(total, "cMgrelease")
    tlabel["cCaMgrelease"] = "Ca+Mg release from feedstock"
    tunits["cCaMgrelease"] = tunits.get("cCarelease", "")
    return tr
