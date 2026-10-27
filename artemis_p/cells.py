"""PHREEQC soil cells. Python counterparts of getcellsfromlayers.m,
setupsoilphases.m and pedotransfer.m.

Index-like fields copied from MATLAB (``cell["ilayer"]``, ``layer["wcell"]``)
keep MATLAB's 1-based values.
"""

import copy

import numpy as np

from .phasetable import check_phase_type

MODNAME = {
    "fc": "fc", "field_capacity": "fc", "wp": "wp", "wilting_point": "wp",
    "por": "por", "porosity": "por", "bulkden": "bulkden", "BD": "bulkden",
    "bulk_density": "bulkden", "SOC": "gSOC", "gSOC": "gSOC", "SON": "gSON",
    "gSON": "gSON", "SOCpct": "SOCpct", "SONpct": "SONpct", "TOC": "SOCpct",
    "TOCpct": "SOCpct", "fSOM": "fSOM", "SOMfrac": "fSOM",
}
MODUNITS = {
    "fc": "volumetric fraction", "wp": "volumetric fraction", "por": "volumetric fraction",
    "bulkden": "g cm^{-3}", "gSOC": "gOC cell^{-1}", "gSON": "gON cell^{-1}",
    "fSOM": "mass fraction", "SOCpct": "mass percent", "SONpct": "mass percent",
}
MODLABEL = {
    "fc": "field capacity", "wp": "wilting point", "por": "porosity",
    "bulkden": "bulk density", "gSOC": "soil organic carbon", "gSON": "soil organic nitrogen",
    "fSOM": "soil organic matter", "SOCpct": "soil organic carbon", "SONpct": "soil organic nitrogen",
}


def _modname(name):
    return MODNAME.get(name, name)


def _apply_layer_parm_scaling(ds):
    """applylayerparmscaling: ds.parm / ds.scale values override or rescale layer fields."""
    parm, layer, nl = ds["parm"], ds["layer"], ds["nlayers"]
    scale = ds.setdefault("scale", {"description": "Scale factors for cell or run arrays"})
    for p in list(parm):
        g = _modname(p)
        if g not in parm:
            parm[g] = parm[p]
        if g in MODNAME:
            layer[g] = np.zeros(nl) + np.ravel(parm[g])
            layer.setdefault("units", {})[g] = MODUNITS[g]
            layer.setdefault("label", {})[g] = MODLABEL[g]
    for p in list(scale):
        g = _modname(p)
        if g not in scale:
            scale[g] = scale[p]
    for f in list(layer):
        if isinstance(layer[f], dict):
            continue
        g = _modname(f)
        if g not in layer:
            layer[g] = layer[f]
            if g in layer.get("units", {}):
                layer["units"][g] = layer["units"][f]
            if g in layer.get("label", {}):
                layer["label"][g] = layer["label"][f]
        if f in parm:
            layer[g] = np.zeros(nl) + parm[f]
            layer["source"][g] = f"Reset as per ds.parm.{f}"
        elif g in parm:
            layer[g] = np.zeros(nl) + parm[g]
            layer["source"][g] = f"Reset as per ds.parm.{g}"
        elif f in scale:
            layer[g] = layer[g] * scale[f]
            layer["source"][g] += f" RESCALED as per ds.scale.{f}"
        elif g in scale:
            layer[g] = layer[g] * scale[g]
            layer["source"][g] += f" RESCALED as per ds.scale.{f}"


def get_cells_from_layers(ds, verbose=False):
    """getcellsfromlayers.m: build ``ds["cell"]`` from ``ds["layer"]``."""
    _apply_layer_parm_scaling(ds)
    parm = ds["parm"]
    if any(k in parm for k in ("cell_thickness", "num_cells", "cell_bottoms")):
        raise NotImplementedError("cells finer than layers are not ported yet")
    ds["cell"] = copy.deepcopy(ds["layer"])
    ds["ncells"] = ds["nlayers"]
    cell = ds["cell"]
    cell["ilayer"] = np.arange(1, ds["nlayers"] + 1, dtype=float)
    ds["layer"]["wcell"] = [float(i) for i in range(1, ds["nlayers"] + 1)]
    _get_cell_phases(ds)
    _get_cell_texture(ds)
    return ds


def _get_cell_phases(ds):
    cell, parm, n = ds["cell"], ds["parm"], ds["ncells"]
    units, label = cell.setdefault("units", {}), cell.setdefault("label", {})
    if "topdepth" in cell and "bottomdepth" in cell:
        cell["middepth"] = (cell["topdepth"] + cell["bottomdepth"]) / 2
        cell["thickness"] = cell["bottomdepth"] - cell["topdepth"]
        cell["frac"] = cell["thickness"] / cell["bottomdepth"][n - 1]
        label["middepth"] = "Cell midpoint"
        units["middepth"] = units["topdepth"]
        label["thickness"] = "Cell thickness"
        units["thickness"] = units["topdepth"]
        label["frac"] = "Amount of model soil in cell"
        # MATLAB sets units.thickness (not units.frac) to 'fraction' here
        units["thickness"] = "fraction"
        cell["vsoil"] = cell["thickness"] * 10  # L/cell (1 L = 1 mm over 1 m2)
        label["vsoil"] = "Soil cell volume"
        units["vsoil"] = "L cell^{-1}"
        cell["depthlab"] = [f"{_g(a)}-{_g(b)}cm" for a, b in zip(cell["topdepth"], cell["bottomdepth"])]
        if "bulkden" not in cell and "bulkden" in parm:
            cell["bulkden"] = np.zeros(n) + parm["bulkden"]
        if "bulkden" in cell:
            cell["gsoil"] = cell["bulkden"] * cell["vsoil"] * 1000  # g soil per cell
    if "fcgrav" in cell and "fc" not in cell and "bulkden" in cell:
        raise NotImplementedError("gravimetric field capacity is not ported yet")
    if "H2Ograv" in cell and "targH2O" not in cell and "bulkden" in cell:
        raise NotImplementedError("gravimetric water content is not ported yet")
    for f in ("por", "fc", "wp"):
        if f in cell and units.get(f, "").lower() in ("percent", "pct", "%"):
            cell[f] = cell[f] / 100
            units[f] = "fraction"
    setup_soil_phases(ds)


def _g(x):
    """MATLAB %g."""
    return f"{x:g}"


def _check_texture_parms(ds):
    cell, parm, n = ds["cell"], ds["parm"], ds["ncells"]
    for k in ("bulkden", "TOC", "Sand", "Silt", "Clay", "field_capacity"):
        if k in parm:
            raise NotImplementedError(f"ds.parm.{k} texture override is not ported yet")
    if "claypct" in cell and "sandpct" in cell:
        cell["siltpct"] = 100 - cell["claypct"] - cell["sandpct"]
    elif "siltpct" in cell and "sandpct" in cell:
        cell["claypct"] = 100 - cell["siltpct"] - cell["sandpct"]
    elif "siltpct" in cell and "claypct" in cell:
        cell["sandpct"] = 100 - cell["claypct"] - cell["siltpct"]
    if "fc" not in cell and "claypct" in cell and "sandpct" in cell:
        pedotransfer(ds, sandpct=cell["sandpct"], claypct=cell["claypct"])


def _get_cell_texture(ds):
    _check_texture_parms(ds)
    cell, n = ds["cell"], ds["ncells"]
    units, label = cell["units"], cell["label"]
    if "gSOC" not in cell:
        if "SOCpct" in cell:
            cell["gSOC"] = cell["SOCpct"] / 100 * cell["gsoil"]
        elif "SOCfrac" in cell:
            cell["gSOC"] = cell["SOCfrac"] * cell["gsoil"]
        elif "gSON" in cell:
            cell["gSOC"] = cell["gSON"] / 0.1
        elif "TC" in cell:
            cell["gSOC"] = cell["TC"]
        else:
            cell["gSOC"] = np.zeros(n)
        label["gSOC"], units["gSOC"] = MODLABEL["gSOC"], MODUNITS["gSOC"]
    if "SOCpct" not in cell:
        if "SOCfrac" in cell:
            cell["SOCpct"] = 100 * cell["SOCfrac"]
        elif "gSOC" in cell:
            cell["SOCpct"] = 100 * cell["gSOC"] / cell["gsoil"]
        else:
            cell["SOCpct"] = np.zeros(n)
        label["SOCpct"], units["SOCpct"] = MODLABEL["SOCpct"], MODUNITS["SOCpct"]
    if "gSON" not in cell:
        if "SONpct" in cell:
            cell["gSON"] = cell["SONpct"] / 100 * cell["gsoil"]
        elif "SONfrac" in cell:
            cell["gSON"] = cell["SONfrac"] * cell["gsoil"]
        elif "gSOC" in cell:
            cell["gSON"] = cell["gSOC"] * 0.1
            print("(getcellsfromlayers) No soil organic nitrogen given; assuming C/N=10")
        else:
            cell["gSON"] = np.zeros(n)
        label["gSON"], units["gSON"] = MODLABEL["gSON"], MODUNITS["gSON"]
    if "fSOM" not in cell:
        cell["fSOM"] = np.zeros(n) + cell["gSOC"] / 0.58 / cell["gsoil"]
        units["fSOM"], label["fSOM"] = MODUNITS["fSOM"], MODLABEL["fSOM"]
    if not all(k in cell for k in ("fc", "wp", "por")):
        old = {k: cell[k] for k in ("fc", "wp", "por") if k in cell}
        pedotransfer(ds)
        if "fc" in old:
            cell["fc"] = old["fc"]
            cell["FC"] = cell["fc"] * cell["vsoil"]
            units["FC"], label["FC"] = "L cell^{-1}", label["fc"]
        if "wp" in old:
            cell["wp"] = old["wp"]
            cell["WP"] = cell["wp"] * cell["vsoil"]
            units["WP"], label["WP"] = "L cell^{-1}", label["wp"]
        if "por" in old:
            cell["por"] = old["por"]
            cell["SAT"] = cell["por"] * cell["vsoil"]
            units["SAT"], label["SAT"] = "L cell^{-1}", label["por"]
    if "FC" not in cell:
        cell["FC"] = cell["fc"] * cell["vsoil"]
        units["FC"], label["FC"] = "L cell^{-1}", label["fc"]
    if "WP" not in cell:
        cell["WP"] = cell["wp"] * cell["vsoil"]
        units["WP"], label["WP"] = "L cell^{-1}", label["wp"]
    if "SAT" not in cell:
        cell["SAT"] = cell["por"] * cell["vsoil"]
        units["SAT"], label["SAT"] = "L cell^{-1}", label["por"]


def pedotransfer(ds, sandpct=None, claypct=None, ptran="Cosby5", include_om=True):
    """pedotransfer.m: soil hydraulic properties from texture (CLM5 / Cosby et al. 1984).

    MATLAB notes: ``if (om_frac>pcalpha)`` on a vector is true only when every
    cell exceeds the percolation threshold, and the organic-matter blend of
    ``bsw`` always uses the Cosby Table 5 formula. Both are kept.
    """
    cell, n = ds["cell"], ds["ncells"]
    sand = np.asarray(cell["sandpct"] if sandpct is None else sandpct, float)
    clay = np.asarray(cell["claypct"] if claypct is None else claypct, float)
    label, units = cell.setdefault("label", {}), cell.setdefault("units", {})
    label["por"] = "Volumetric soil water at saturation watsat"
    units["por"] = "fraction"
    label["bsw"] = "Clapp-Hornberger shape parameter b"
    label["sucsat"], units["sucsat"] = "Minimum soil suction", "mm"
    label["xksat"], units["xksat"] = "Maximum saturated hydraulic conductivity", "mm s^{-1}"
    label["ksatom"], units["ksatom"] = "Saturated hydraulic conductivity organic matter", "mm s^{-1}"
    label["Ksatmmday"], units["Ksatmmday"] = "Saturated hydraulic conductivity", "mm day^{-1}"
    if ptran == "Cosby5":
        label["ptran"] = "Table 5, Cosby et al 1984 Water ResRes 20:682-690"
        watsat = 0.489 - 0.00126 * sand
        sucsat = 10 * (10 ** (1.88 - 0.0131 * sand))
        bsw = 2.91 + 0.159 * clay
        xksat = 0.0070556 * (10 ** (-0.884 + 0.0153 * sand))
    elif ptran == "Cosby4":
        label["ptran"] = "Table 4, Cosby et al 1984 Water ResRes 20:682-690"
        watsat = 0.505 - 0.00142 * sand - 0.00037 * clay
        sucsat = 10 * (10 ** (1.54 - 0.0095 * sand + 0.0063 * (100 - sand - clay)))
        bsw = 3.10 + 0.157 * clay - 0.003 * sand
        xksat = 0.0070556 * (10 ** (-0.60 + 0.0126 * sand - 0.0064 * clay))
    elif ptran == "NoilhanLacarrere1995":
        label["ptran"] = "Noilhan and Lacarrere 1995"
        watsat = -0.00108 * sand + 0.494305
        sucsat = 10 ** (-0.0088 * sand + 2.85)
        bsw = 0.137 * clay + 3.501
        xksat = 10 ** (-0.0582 * clay - 0.00091 * sand + 0.000529 * clay ** 2 - 0.0001203 * sand ** 2 - 1.38)
    else:
        raise ValueError(ptran)
    if "fSOM" in cell and include_om:
        pcbeta, pcalpha, zsapric = 0.139, 0.5, 50
        zsaprat = cell["middepth"] / zsapric
        om_watsat = np.maximum(0.93 - 0.1 * zsaprat, 0.83)
        om_b = np.minimum(2.7 + 9.3 * zsaprat, 12.0)
        om_sucsat = np.minimum(10.3 - 0.2 * zsaprat, 10.1)
        om_hksat = np.maximum(0.28 - 0.2799 * zsaprat, xksat)
        perc_frac = 0.0
        om_frac = cell["fSOM"] * cell["gsoil"] / cell["vsoil"] / 130  # volumetric OM fraction
        if np.all(om_frac > pcalpha):
            perc_norm = (1 - pcalpha) ** (-pcbeta)
            perc_frac = perc_norm * (om_frac - pcalpha) ** pcbeta
        uncon_frac = (1 - om_frac) + (1 - perc_frac) * om_frac
        if np.all(om_frac < 1):
            uncon_hksat = uncon_frac / ((1 - om_frac) / xksat + ((1 - perc_frac) * om_frac) / om_hksat)
        else:
            uncon_hksat = 0
        hksat = uncon_frac * uncon_hksat + (perc_frac * om_frac) * om_hksat
        watsat = (1 - om_frac) * watsat + om_watsat * om_frac
        bsw = (1 - om_frac) * (2.91 + 0.159 * clay) + om_frac * om_b
        sucsat = (1 - om_frac) * sucsat + om_sucsat * om_frac
        zero = np.zeros(n)
        cell["om_hksat"] = zero + om_hksat
        cell["om_watsat"] = zero + om_watsat
        cell["om_bsw"] = zero + om_b
        cell["om_sucsat"] = zero + om_sucsat
        cell["om_frac"] = om_frac
    else:
        hksat = xksat
    zero = np.zeros(n)
    watfc = watsat * (0.1 / (hksat * 86400)) ** (1 / (2 * bsw + 3))
    smpmax = 2.57e5
    wilt = ((smpmax / sucsat) ** (-1 / bsw)) * watsat
    cell["watfc"] = zero + watfc
    cell["wilt"] = zero + wilt
    cell["fc"] = cell["watfc"]
    cell["wp"] = cell["wilt"]
    cell["sucsat"] = zero + sucsat
    cell["bsw"] = zero + bsw
    cell["watsat"] = zero + watsat
    cell["por"] = zero + watsat
    if "vsoil" in cell:
        cell["FC"] = cell["fc"] * cell["vsoil"]
        cell["WP"] = cell["wp"] * cell["vsoil"]
        cell["SAT"] = cell["por"] * cell["vsoil"]
    cell["xksat"] = zero + xksat
    cell["hksat"] = zero + hksat
    cell["pedotransferfunction"] = ptran
    cell["Ksatmmday"] = hksat * 86400
    return ds


def setup_soil_phases(ds):
    """setupsoilphases.m: native mineral moles and reactive surface areas per cell.

    MATLAB note (bug kept for now): the RSA of each mineral uses the mass of
    all soil minerals, not the mineral's own share. See docs/PORTING_PLAN.md.
    """
    t = ds["phasetable"]
    if "category" not in t:
        t = t.copy()
        t["category"] = ""
        ds["phasetable"] = t
    cell = ds["cell"]
    rows = {str(r).lower(): i for i, r in enumerate(t.index)}
    wmin, wcellmin = [], []
    for f in list(cell):
        cfld = f.replace("_native", "").replace("_soil", "").replace("_parent", "")
        eqflag = "_eq" in f.lower()
        if eqflag:
            cfld = cfld[: cfld.index("_eq")]
        kinflag = "_kin" in f.lower()
        if kinflag:
            cfld = cfld[: cfld.index("_kin")]
        if cfld.lower() not in rows:
            continue
        w = rows[cfld.lower()]
        t = ds["phasetable"]
        if t["category"].iloc[w] == "":
            t = t.copy()
            t.iloc[w, t.columns.get_loc("category")] = "native"
            ds["phasetable"] = t
        wmin.append(w)
        wcellmin.append(f)
        if eqflag:
            ds = check_phase_type(None, ds, seteq=[cfld])
        if kinflag:
            ds = check_phase_type(None, ds, setkin=[cfld])
        rows = {str(r).lower(): i for i, r in enumerate(ds["phasetable"].index)}
        cell[cfld] = cell[f]
    if not wmin:
        raise NotImplementedError("native phases from the phase table's feedstock column are not ported yet")
    t = ds["phasetable"]
    phase = [t["MATLABphase"].iloc[w] for w in wmin]
    G = [float(t["G"].iloc[w]) for w in wmin]
    molmass = [float(t["molmass"].iloc[w]) for w in wmin]
    prop = {}
    for ph, f in zip(phase, wcellmin):
        v = np.asarray(cell[f], float)
        u = cell.get("units", {}).get(f, "")
        if "%" in u or "pct" in u or "percent" in u:
            v = v / 100
        prop[ph] = v
    n = ds["ncells"]
    parm = ds["parm"]
    parm.setdefault("nativeRSAscale", 1)
    native_rsa_scale = parm["nativeRSAscale"]
    parm.setdefault("label", {})["nativeRSAscale"] = "Scale factor for native mineral RSA"
    c2dwt = parm.get("C2DWT", 0.58)
    clay, silt, sand = cell["claypct"], cell["siltpct"], cell["sandpct"]
    units, label = cell["units"], cell["label"]
    if "SOCpct" in cell and "orgpct" in cell:
        pass
    elif "SOCpct" in cell or "TOC" in cell:
        raise NotImplementedError("orgpct from SOCpct/TOC is not ported yet")
    elif "orgpct" in cell:
        cell["SOCpct"] = cell["orgpct"] * c2dwt
        label["SOCpct"] = "Organic carbon content"
        units["SOCpct"] = "%"
    pebbles = cell.get("pebbles", cell.get("gravel", np.zeros(n)))
    cobbles = cell.get("cobbles", cell.get("boulders", np.zeros(n)))
    psd = np.column_stack([cobbles, pebbles, sand, silt, clay])
    psdthresh = np.array([1000, 64, 2, 0.05, 0.002])  # size-class thresholds, mm
    totorg = cell["orgpct"] / 100
    minfrac = 1 - totorg
    pmd = cell["gsoil"][:, None] * psd / 100 * minfrac[:, None]
    moles, rsa = cell.setdefault("moles", {}), cell.setdefault("RSA", {})
    for m, ph in enumerate(phase):
        rsad = pmd / psdthresh / 1000 / G[m] * 6
        mass = pmd.sum(axis=1) * prop[ph]
        moles[ph] = mass / molmass[m]
        rsa[ph] = rsad.sum(axis=1) * native_rsa_scale
    return ds
