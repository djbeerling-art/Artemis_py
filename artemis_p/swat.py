"""Plant growth, evapotranspiration partitioning and nutrient uptake. Python
counterpart of daySWAT7.m, which follows the SWAT 2009 crop equations.

Arrays are 1-D (days) for whole-column quantities and 2-D (cells x days) for
per-cell ones; MATLAB's ``circshift`` differences (which wrap the last day
round to the first) are reproduced with ``np.roll``.
"""

import numpy as np

from .fixes import enabled as fixed

from .csvdata import read_csv_data


def _diff_wrap(a):
    """MATLAB a - circshift(a,1) along days."""
    return a - np.roll(a, 1, axis=-1)


def day_swat(ds, zrmxfrPHU=0.4, bwater_default=10, ClKmassrat=0.3):
    """daySWAT7.m for the default options (no potbound or root-tip kludges)."""
    parm, run, cell = ds["parm"], ds["run"], ds["cell"]
    parm.setdefault("nutaqfrac", 1)
    parm.setdefault("nutlist", ["P", "K", "Ca", "Mg", "S", "Na", "Si"])
    parm["zrmxfrPHU"] = zrmxfrPHU
    sclunits = ds["units"]["scl"] if "sclunits" in ds else "m^{-2}land"
    incunits = " day^{-1}" if ds["datastep"] == 1 else (
        "month^{-1}" if ds["datastep"] == 30.5 else " timestep^{-1}")
    massunits = f"g {sclunits}"
    n, nc = ds["nrun"], ds["ncells"]
    zero = np.zeros(n)
    uzero = np.zeros((nc, n))
    ds.setdefault("dbdir", "")
    units, label = run.setdefault("units", {}), run.setdefault("label", {})

    # initarrays
    if "ET" not in run:
        run["ET"] = zero + parm.get("ETday", parm.get("ETyear", 0) / 365 if "ETyear" in parm else parm.get("ET", 0))
        label["ET"], units["ET"] = "Total evapotranspiration", "mm day^{-1}"
    if "zroot" not in run:
        run["zroot"], run["dzroot"] = zero.copy(), zero.copy()
        label["zroot"], units["zroot"] = "Max rooting depth", "cm"
        label["dzroot"], units["dzroot"] = "Max rooting depth increment", "cm"
    elif "zrootmx" in parm:
        run["zroot"] = np.fmin(parm["zrootmx"], run["zroot"])
    if "rootmass" not in run:
        run["rootmass"] = zero.copy()
        label["rootmass"], units["rootmass"] = "Root mass", massunits
    if "totmass" not in run:
        run["totmass"] = zero.copy()
        label["totmass"], units["totmass"] = "Total live biomass", massunits
    if "LAI" not in run:
        run["LAI"] = zero.copy()
        label["LAI"], units["LAI"] = "Leaf area index", f"m^{{2}} {sclunits}"
    if "plant" not in run:
        run["plant"] = zero.copy()
        label["plant"] = "Presence-absence flag for plant"
    label["totmassinc"], units["totmassinc"] = "Total live mass increment", massunits + incunits
    label["rootmassinc"], units["rootmassinc"] = "Root mass increment", massunits + incunits
    label["LAIinc"], units["LAIinc"] = "LAI increment", f"m^{{2}} {sclunits}{incunits}"
    up = run["uptake"] = {
        "fracET": uzero.copy(), "label": {"fracET": "Evapotranspiration", "transpiration": "Transpiration",
                                          "fractransp": "Fraction of transpiration in ET from cell",
                                          "ctranspiration": "Cumulative transpiration",
                                          "fracroot": "Fraction of roots in soil layers"},
        "units": {"fracET": "Fraction from cell", "transpiration": "mm cell^{-1}",
                  "ctranspiration": "mm cell^{-1}"},
        "fractransp": uzero.copy(), "transpiration": uzero.copy(), "ctranspiration": uzero.copy(),
        "fracroot": uzero.copy()}
    run["netET"] = uzero.copy()
    run["rootmassfrac"] = uzero.copy()
    label["rootmassfrac"], units["rootmassfrac"] = "Fraction of final total root biomass", "g g^{-1} layer^{-1}"
    run["rootmassdayfrac"] = uzero.copy()
    label["rootmassdayfrac"], units["rootmassdayfrac"] = "Fraction of current total root biomass", "g g^{-1} layer^{-1}"

    if "snowmelt" not in run:
        raise NotImplementedError("SWATsnow (snow from air temperature) is not ported yet")
    _canopy_interception(ds)
    run["totET"] = np.asarray(run["ET"], float).copy()
    label["totET"] = "Total evapotranspiration"
    if "ET" in units:
        units["totET"] = units["ET"]
    run["ET"] = run["ET"] - run["canevap"]
    top, bot = (np.asarray(cell[k], float) * 10 for k in ("topdepth", "bottomdepth"))
    cell["evapfrac"] = bot / (bot + np.exp(2.374 - 0.00713 * bot)) - top / (top + np.exp(2.374 - 0.00713 * top))
    cell["label"]["evapfrac"] = "Fraction of evaporation from cells on bare ground"
    up["fracET"] = uzero + cell["evapfrac"][:, None]
    run["cellET"] = cell["evapfrac"][:, None] * run["ET"][None, :]

    wr, wrun, wswat = _find_swat_crop(ds)
    _growth_curves(ds, wr, wrun, wswat, zero, incunits, sclunits)
    parm.setdefault("bwater", bwater_default)
    bwater = parm["bwater"]
    cover = np.fmin(run["LAI"], 1)
    bare = 1 - cover
    run["cover"] = cover
    _char_depth(ds, cover, bwater)
    _layer_et(ds, wr, cover, bare, bwater)
    _leg_stage(ds, wr, wswat, zero)
    if "cropstoichfile" not in ds:
        ds["cropstoichfile"] = f"{ds['dbdir']}/cropstoichiometry.csv"
    cs = ds["cropstoich"] = read_csv_data(ds["cropstoichfile"])
    nut = cs["nutrient"]
    wN, wP = nut.index("N"), nut.index("P")
    stoichN = np.zeros((len(nut), n))
    _np_uptake(ds, wr, wswat, zero, stoichN, wN, wP, sclunits, incunits)
    _cell_up(ds, uzero, zero, sclunits, incunits)
    if len(wr):
        _nut_uptake(ds, stoichN, wN, zero, uzero)
    if parm["usemeanforcings"] == 1:
        _use_means(ds["run"], zero)
    return ds


def _mean_row(x):
    """MATLAB x.*0+mean(x,2): every day set to the run mean (NaN propagates)."""
    x = np.asarray(x, float)
    return x * 0 + x.mean(axis=-1, keepdims=True)


def _use_means(run, zero):
    """Spin-ups (usemeanforcings): mean LAI, masses, root depth and uptake."""
    print("daySWAT7 Using mean values for LAI, masses, root depth")
    for f in ("LAI", "totmass", "rootmass", "zroot"):
        run[f] = zero + np.asarray(run[f], float).mean(axis=-1, keepdims=True)
    run["plant"] = (run["totmass"] > 0).astype(float)
    run["cellup"] = _mean_row(run["cellup"])
    for f, v in run["uptake"].items():
        if isinstance(v, np.ndarray):
            run["uptake"][f] = _mean_row(v)
    for el in list(run["totuptake"]):
        for grp in ("totuptake", "totcuptake", "cellcuptake", "celluptake"):
            run[grp][el] = _mean_row(run[grp][el])


def _canopy_interception(ds):
    """getcanopyinterception. daySWAT7 calls it before the crop is looked up,
    so ds.swat does not exist yet and only zero arrays are set."""
    run = ds["run"]
    if "canint" in run:
        return
    run["canint"] = np.zeros(ds["nrun"])
    run["label"]["canint"] = "Canopy interception"
    if "canevap" in run:
        return
    run["canevap"] = np.zeros(ds["nrun"])
    run["label"]["canevap"] = "Canopy evaporation"
    if "precip" not in run:
        return
    run["units"]["canint"] = run["units"]["precip"]
    if "ET" in run["units"]:
        run["units"]["canevap"] = run["units"]["ET"]
    if "swat" in ds:
        raise NotImplementedError("canopy storage with a preloaded SWAT table is not ported yet")


def _find_swat_crop(ds):
    """findSWATcrop: index into the SWAT crop table for each run day."""
    run = ds["run"]
    n = ds["nrun"]
    wr = np.flatnonzero(run["plant"])
    if "swat" not in ds:
        if "SWATdbfile" not in ds:
            ds["SWATdbfile"] = f"{ds['dbdir']}/SWAT_AppendixA_crop_params4.csv"
            print(f"(daySWAT7/findSWATcrop) ds.SWATdb not defined, using {ds['SWATdbfile']}")
        ds["swat"] = read_csv_data(ds["SWATdbfile"])
    sw = ds["swat"]
    keys = [[str(x).lower() for x in sw[c]] for c in ("crop", "altcrop", "code")]
    wrun = np.zeros(n, dtype=int)
    cache = {}
    for j, c in enumerate(run["crop"]):
        c = c.lower()
        if c not in cache:
            cache[c] = next((k.index(c) + 1 for k in keys if c in k), 0)
        wrun[j] = cache[c]
    if wrun.sum() == 0:
        print("(daySWAT7/findSWATcrop) No matching plant params in SWAT crop database")
        return wr, wrun, np.zeros(0, dtype=int)
    wr = np.flatnonzero(wrun > 0)
    wswat = wrun[wr] - 1  # 0-based rows of the SWAT table
    sw["wswat"] = wswat[0] + 1
    return wr, wrun, wswat


def _growth_curves(ds, wr, wrun, wswat, zero, incunits, sclunits):
    """growthcurves: heat units, root depth, biomass and LAI through each season."""
    run, parm, sw = ds["run"], ds["parm"], ds["swat"]
    units, label = run["units"], run["label"]
    if not len(wr):
        return
    if np.sum(run["LAI"] > 0) and np.sum(run["totmass"]) > 0:
        return
    run["PHU"], run["frPHU"] = zero.copy(), zero.copy()
    label["PHU"], units["PHU"] = "Phenological heat units attained", "degree-days"
    label["frPHU"], units["frPHU"] = "Fraction of total phenological heat units attained", ""
    run["vegmass"] = zero.copy()
    label["vegmass"], units["vegmass"] = "Vegetative mass", f"g {sclunits}"
    run["grainmass"] = zero.copy()
    label["grainmass"], units["grainmass"] = "Grain mass", units["vegmass"]
    run["frLAI"] = zero.copy()
    label["frLAI"], units["frLAI"] = "Fraction of LAI attained", ""
    label["legstage"], units["legstage"] = "Growth stage factor for legumes", ""
    label["vegmassinc"], units["vegmassinc"] = "Vegetative mass increment", units["vegmass"] + incunits
    label["grainmassinc"], units["grainmassinc"] = "Grain mass increment", units["vegmassinc"]
    run["canhgt"] = zero.copy()
    label["canhgt"] = "Canopy height"
    if "peaktot" not in run:
        run["peaktot"], run["peakroot"] = zero.copy(), zero.copy()
    airtemp = np.asarray(run["airtemp"], float)
    grow = zero.copy()
    grow[wr] = airtemp[wr] > sw["Tbase"][wswat]
    run["grow"] = grow
    label["grow"], units["grow"] = "Air temperature high enough for crop growth", "0=no, 1=yes"
    if "zrootmx" not in parm:
        parm["zrootmx"] = sw["zrootmx"][sw["wswat"] - 1] * 100
        print(f"(daySWAT7/growthcurves) Setting ds.parm.zrootmx to {parm['zrootmx']:f} cm "
              f"using value for {sw['crop'][sw['wswat'] - 1]}")
    zrootmx = zero.copy()
    zrootmx[wr] = parm["zrootmx"]
    PHUsen = zero.copy()
    PHUsen[wr] = sw["frPHUsen"][wswat]
    plant, peaktot, harvest = run["plant"], run["peaktot"], run["harvest"]
    isow = ihar = ipeak = -1  # 0-based; -1 means unset
    for i in range(ds["nrun"]):
        if isow < 0 and plant[i] > 0:
            isow = i
        if peaktot[i] > 0:
            ipeak = i
        if harvest[i] > 0:
            ihar = i
        elif isow >= 0 and plant[i] == 0:
            ihar = i - 1
        if ihar >= 0:
            if isow < 0:
                raise ValueError("harvest before any sowing (MATLAB would fail indexing day 0)")
            s = slice(isow, ihar + 1)
            rows = wrun[s] - 1
            if np.any(rows < 0):
                raise ValueError("fallow day inside a growing season (MATLAB would fail indexing row 0)")
            run["PHU"][s] = np.cumsum(airtemp[s] * grow[s] * plant[s])
            run["frPHU"][s] = run["PHU"][s] / run["PHU"][ihar] * sw["frPHUhar"][rows]
            run["zroot"][s] = np.fmax(run["sowdepth"][isow],
                                         np.fmin(1, run["frPHU"][s] / parm["zrmxfrPHU"]) * zrootmx[s])
            run["totmass"][s] = peaktot[ipeak] if ipeak >= 0 else _bad_index("peaktot")
            run["rootmass"][s] = run["peakroot"][ipeak]
            gf = (run["frPHU"][s] >= PHUsen[s]) & (run["frPHU"][s] <= 1.0)
            ngf = gf.sum()
            if ngf > 0:
                run["grainmass"][s] = harvest[ihar] * np.cumsum(gf) / ngf
            else:
                run["grainmass"][ipeak + 1: ihar + 1] = harvest[ihar]
            run["grainmass"][s] = run["grainmass"][s] * sw["HIopt"][rows]
            isow = ihar = -1
    if "LAIshp1" in parm:
        sw["LAIshp1"][wswat] = parm["LAIshp1"]
    if "LAIshp2" in parm:
        sw["LAIshp2"][wswat] = parm["LAIshp2"]
    run["dzroot"] = _diff_wrap(run["zroot"])
    frPHU = run["frPHU"]
    run["frLAI"], run["canhgt"] = zero.copy(), zero.copy()
    f = frPHU[wr]
    run["frLAI"][wr] = f / (f + np.exp(sw["LAIshp1"][wswat] - sw["LAIshp2"][wswat] * f))
    run["totmass"] = run["totmass"] * run["frLAI"]
    run["rootmass"] = run["rootmass"] * run["frLAI"]
    run["vegmass"] = run["totmass"] - run["grainmass"]
    run["vegmassinc"] = _diff_wrap(run["vegmass"])
    run["rootmassinc"] = _diff_wrap(run["rootmass"])
    run["totmassinc"] = _diff_wrap(run["totmass"])
    run["canhgt"][wr] = sw["hcmx"][wswat] * np.sqrt(run["frLAI"][wr])
    units["canhgt"] = sw["units"]["hcmx"]
    run["LAI"], run["LAIinc"] = zero.copy(), zero.copy()
    w = (f > sw["frPHUsen"][wswat]) & (f <= 1)
    run["frLAI"][wr[w]] = run["frLAI"][wr[w]] * (1 - f[w]) / (1 - sw["frPHUsen"][wswat[w]])
    run["frLAI"][frPHU > 1] = 0
    if "column" in ds["units"]["scl"] and "hydrocol" in ds:
        raise NotImplementedError("column-scale LAI is not ported yet")
    run["LAI"][wr] = sw["LAImx"][wswat] * run["frLAI"][wr]
    run["LAIinc"] = _diff_wrap(run["LAI"])


def _bad_index(name):
    raise ValueError(f"no {name} before harvest (MATLAB would fail indexing day 0)")


def _char_depth(ds, cover, bwater):
    """getchardepth: depth (cm) whose share of root water uptake is closest to the mean."""
    run, cell = ds["run"], ds["cell"]
    if "liveroot" not in run:
        run["liveroot"] = run["zroot"].copy()
    dmax = int(np.ceil(np.asarray(cell["bottomdepth"])[-1]))
    run["tcdepth"] = np.zeros(ds["nrun"])
    wlr = np.flatnonzero(run["liveroot"] > 0)
    lr = run["liveroot"][wlr]
    idep = np.arange(1, dmax + 1)[:, None]
    dfr = np.zeros((dmax, ds["nrun"]))
    dfr[:, wlr] = cover[wlr] / (1 - np.exp(-bwater)) * (
        (1 - np.exp(-bwater * idep / lr)) - (1 - np.exp(-bwater * (idep - 1) / lr)))
    meandfr = dfr.mean(axis=0)
    idx = np.argsort(np.abs(dfr - meandfr), axis=0, kind="stable")
    run["tcdepth"][wlr] = idx[0, wlr] + 1


def _layer_et(ds, wr, cover, bare, bwater):
    """getlayerET: share of ET and transpiration taken from each cell."""
    run, cell, up = ds["run"], ds["cell"], ds["run"]["uptake"]
    top_c, bot_c = np.asarray(cell["topdepth"], float), np.asarray(cell["bottomdepth"], float)
    lr = run["liveroot"][wr]
    with np.errstate(divide="ignore", invalid="ignore"):
        for i in range(ds["ncells"]):
            up["fracET"][i] = cell["evapfrac"][i] * bare
            bottom = np.fmin(run["zroot"][wr], bot_c[i])
            topd = np.fmin(top_c[i], bottom)
            ft = cover[wr] / (1 - np.exp(-bwater)) * (
                (1 - np.exp(-bwater * bottom / lr)) - (1 - np.exp(-bwater * topd / lr)))
            up["fractransp"][i, wr] = ft
            up["transpiration"][i, wr] = run["ET"][wr] * ft
            up["ctranspiration"][i, wr] = np.cumsum(up["transpiration"][i, wr])
            up["fracroot"][i, wr] = ft
    bad = np.isnan(up["fracroot"])
    up["fracroot"][bad] = 0
    up["fractransp"][bad] = 0
    s = up["fracroot"][:, wr].sum(axis=0)
    wrr = wr[s > 0]
    up["fracroot"][:, wrr] = up["fracroot"][:, wrr] / up["fracroot"][:, wrr].sum(axis=0)
    up["fracET"][:, wr] = up["fracET"][:, wr] + up["fractransp"][:, wr]
    run["cellET"] = np.fmax(0, run["ET"])[None, :] * up["fracET"]
    run["label"]["cellET"] = "Evapotranspiration"
    run["units"]["cellET"] = f"{ds['units']['mmscl']} cell^{{-1}}"


def _leg_stage(ds, wr, wswat, zero):
    """getlegstage: SWAT legume growth-stage factor (eqns 5:2.3.10-14)."""
    run = ds["run"]
    run["legstage"] = zero.copy()
    if ds["parm"]["usemeanforcings"] == 1 or "frPHU" not in run:
        return
    f = run["frPHU"]
    ls = run["legstage"]
    w = (f > 0.15) & (f <= 0.3)
    ls[w] = 6.67 * f[w] - 1
    ls[(f > 0.3) & (f <= 0.55)] = 1
    w = (f > 0.55) & (f <= 0.75)
    ls[w] = 3.75 - 5 * f[w]
    ls[wr] = ls[wr] * ds["swat"]["legume"][wswat]


def _np_uptake(ds, wr, wswat, zero, stoichN, wN, wP, sclunits, incunits):
    """getNPuptake: N and P mass fractions and masses in the crop (SWAT eqn 5:2.3.1)."""
    run, sw, cs = ds["run"], ds["swat"], ds["cropstoich"]
    units, label = run["units"], run["label"]
    run["frN"], run["frP"] = zero.copy(), zero.copy()
    label["frN"], label["frP"] = "Fraction of N in biomass each day", "Fraction of P in biomass each day"
    if "frPHU" in run:
        f = run["frPHU"][wr]
        for el, w, mm in (("N", wN, 14.01), ("P", wP, 30.97)):
            a1, a3 = sw[f"fr{el}1"][wswat], sw[f"fr{el}3"][wswat]
            s1, s2 = sw[f"{el}shp1"][wswat], sw[f"{el}shp2"][wswat]
            run[f"fr{el}"][wr] = (a1 - a3) * (1 - f / (f + np.exp(s1 - s2 * f))) + a3
            cs["conc"][w] = sw[f"fr{el}3"][wswat[0]] / mm * 1e6
            cs["massfrac"][w] = sw[f"fr{el}3"][wswat[0]]
            cs["source"][w] = f"SWAT2009 crop database fr{el}3"
            if el == "N":
                stoichN[:, wr] = 14.01 / sw["frN3"][wswat]
    else:
        run["frN"][wr] += cs["massfrac"][wN]
        run["frP"][wr] += cs["massfrac"][wP]
        stoichN[:, wr] = 14.01 / np.nanmax(run["frN"][wr])
    for el in ("N", "P"):
        run[f"{el}mass"] = run[f"fr{el}"] * run["totmass"]
        run[f"{el}massinc"] = np.fmax(0, _diff_wrap(run[f"{el}mass"]))
        units[f"{el}mass"] = f"g{el} {sclunits}"
        units[f"{el}massinc"] = units[f"{el}mass"] + incunits
    label["Nmass"], label["Nmassinc"] = "Plant N mass", "Plant N mass increment"
    label["Pmass"], label["Pmassinc"] = "Plant P mass", "Plant P mass increment"


def _cell_up(ds, uzero, zero, sclunits, incunits):
    """getcellup: share of nutrient demand met from each cell, and N and P uptake
    rescaled so each season's uptake adds up to its peak plant N and P."""
    run, cell, parm = ds["run"], ds["cell"], ds["parm"]
    parm.setdefault("nupdis", 20)
    nupdis = parm["nupdis"]
    top_c, bot_c = np.asarray(cell["topdepth"], float), np.asarray(cell["bottomdepth"], float)
    run["cellup"] = uzero.copy()
    for i in range(ds["ncells"]):
        w = np.flatnonzero(run["zroot"] > top_c[i])
        if not len(w):
            continue
        lr = run["liveroot"][w]
        bottom = lr if i == ds["ncells"] - 1 else np.fmin(lr, bot_c[i])
        topup = 1 - np.exp(-nupdis * top_c[i] / lr)
        bottomup = 1 - np.exp(-nupdis * bottom / lr)
        run["cellup"][i, w] = (bottomup - topup) / (1 - np.exp(-nupdis))
    units, label = run["units"], run["label"]
    label["cellup"] = "Layer fraction of nutrient demand"
    tot = run["totuptake"] = {"N": run["Nmassinc"] / 14.01, "P": run["Pmassinc"] / 30.97}
    cu = run["celluptake"] = {"N": tot["N"][None, :] * run["cellup"], "P": tot["P"][None, :] * run["cellup"]}
    ccu = run["cellcuptake"] = {"N": uzero.copy(), "P": uzero.copy()}
    label["celluptake"], units["celluptake"] = "Uptake from each cell", f"mol cell^{{-1}} {incunits}"
    label["cellcuptake"], units["cellcuptake"] = "Cumulative uptake from each cell", "mol cell^{-1}"
    label["totuptake"], units["totuptake"] = "Uptake from soil", f"mol {sclunits}{incunits}"
    label["totcuptake"], units["totcuptake"] = "Cumulative uptake from soil", f"mol {sclunits}"
    plant = run["plant"]
    run["GS"] = plant * 0
    if not np.any(plant > 0):
        print("(daySWAT7/getcellup) No plants detected")
        return
    dplant = _diff_wrap(plant)
    dplant[0] = plant[0]
    wstart = np.flatnonzero(dplant > 0)
    tcu = run["totcuptake"] = {"N": zero.copy(), "P": zero.copy()}
    if not len(wstart):
        raise NotImplementedError("non-seasonal plant growth is not ported yet")
    for s in wstart:
        run["GS"][s:] = run["GS"][s:] + plant[s:]
    print(f"(daySWAT7/getcellup) {len(wstart)} growing seasons detected")
    for igs in range(1, len(wstart) + 1):
        w = np.flatnonzero(run["GS"] == igs)
        tm = run["totmass"][w]
        # MATLAB max skips NaN (a one-day season after a harvest gives 0/0 in frPHU)
        wmax = w[np.flatnonzero(tm == np.nanmax(tm))[0]]
        for el in ("N", "P"):
            scl = run[f"{el}mass"][wmax] / np.sum(run[f"{el}massinc"][w])
            run[f"{el}massinc"][w] *= scl
            tot[el][w] *= scl
            if el == "N":
                cu[el][:, w] *= scl
            else:
                # MATLAB rescales only the top cell's P uptake here (celluptake.P(1,w))
                if fixed(ds, "p_uptake_all_cells"):
                    cu[el][:, w] *= scl
                else:
                    cu[el][0, w] *= scl
            run[f"{el}mass"][w] = np.cumsum(run[f"{el}massinc"][w])
            tcu[el][w] = np.cumsum(tot[el][w])
            ccu[el][:, w] = np.cumsum(cu[el][:, w], axis=1)


def _nut_uptake(ds, stoichN, wN, zero, uzero):
    """getnutuptake: other nutrients taken up in fixed ratios to N."""
    run, cs = ds["run"], ds["cropstoich"]
    nut = cs["nutrient"]
    nutflag = np.ones(len(nut), dtype=bool)
    tot = run["totuptake"]
    wNgood = np.flatnonzero(tot["N"] > 0)
    nutrat = run.setdefault("nutrat", {})
    if len(wNgood):
        nutrat["P"] = np.zeros(ds["nrun"]) if "P" not in nutrat else nutrat["P"]
        nutrat["P"][wNgood] = tot["P"][wNgood] / tot["N"][wNgood]
    cs.setdefault("units", {})["stoichN"] = "mol mol^{-1} N"
    cs.setdefault("label", {})["stoichN"] = "Nutrient content at maturity"
    cs["stoichN"] = np.zeros(len(nut))
    for i, el in enumerate(nut):
        cs["stoichN"][i] = cs["conc"][i] / cs["conc"][wN]
        if el.lower() in ("c", "h", "o", "n", "mo") or el not in ds["solute"]:
            continue
        if cs["crop"][i].lower() == "default":
            wcrop = np.flatnonzero(run["totmass"] > 0)
        else:
            # MATLAB contains() is case-sensitive: 'Maize' rows never match crop 'maize'
            if fixed(ds, "crop_case"):
                wcrop = np.flatnonzero([cs["crop"][i].lower() in c.lower() for c in run["crop"]])
            else:
                wcrop = np.flatnonzero([cs["crop"][i] in c for c in run["crop"]])
            if not len(wcrop):
                continue
        stoichN[i, wcrop] = cs["conc"][i] / 1e6 * stoichN[i, wcrop]
        if nutflag[i]:
            tot[el] = zero.copy()
            run["totcuptake"][el] = zero.copy()
            run["celluptake"][el] = uzero.copy()
            run["cellcuptake"][el] = uzero.copy()
            nutrat[el] = uzero.copy()
            nutflag[[j for j, x in enumerate(nut) if x.lower() == el.lower()]] = False
        tot[el][wcrop] = tot["N"][wcrop] * stoichN[i, wcrop]
        run["totcuptake"][el][wcrop] = run["totcuptake"]["N"][wcrop] * stoichN[i, wcrop]
        run["celluptake"][el][:, wcrop] = run["celluptake"]["N"][:, wcrop] * stoichN[i, wcrop]
        run["cellcuptake"][el][:, wcrop] = run["cellcuptake"]["N"][:, wcrop] * stoichN[i, wcrop]
        nutrat[el] = zero.copy()
        if len(wNgood):
            nutrat[el][wNgood] = tot[el][wNgood] / tot["N"][wNgood]
    run["units"]["nutrat"] = "mol mol^{-1}N"
