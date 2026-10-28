"""Rice paddy mode: the SCEPTER paddy chemistry (paddy/paddy.dat) in each ARTEMIS cell.

Turned on by ``ds['paddy']``, a dict with

``flooded``  1 on the days whose cells start the day saturated (one value per run day, or a
             function of ds that returns them)
``q``        percolation under flood, mm/day (default 0.5 m/yr)
``pools``    initial solid pools, mol per m3 of bulk soil (``Ferrihydrite``, ``Goethite_k``, ...),
             a number or one value per cell
``omin``     organic C input, mol C m-3 yr-1 per cell or a function of depth in m (default: 450 g C m-2 yr-1 spread as
             exp(-z/0.3), SCEPTER's zsupp_plant)
``uin``      urea input, mol m-3 yr-1 per cell, or a function of depth in m (default 0)
``n2c``      mol NH3 released per mol CH2O decomposed (default 0.1)
``keb``      CH4 ebullition rate above 0.4 atm, yr-1 (default 1000, spec 6.1)
``kfw``      floodwater gas transfer velocity, m/day (default 0.1, spec 6.2)

The process spec is paddy/docs/riceland_process_spec.md on SCEPTER branch paddy/fe3-reduction.
In paddy mode the SWAT organic matter and N kinetics are off: the spec's acceptor ladder
replaces them. Crop uptake stays.
"""

from pathlib import Path

import numpy as np

DAT = Path(__file__).resolve().parent.parent / "paddy" / "paddy.dat"

N2C = 0.1  # mol NH3 per mol CH2O (SCEPTER default; the spec's production value is 0.03)

FORMULA = {  # net change to the solution per mol of extent
    "g2o2": "C 1 O 3 H 2 Oxg -1",
    "g2fe3": "C 1 Fe 4 O 13 H 14",
    "g2gr": "C 1 Fe 6 O 17 H 18",
    "g2gt": "C 1 Fe 4 O 9 H 6",
    "g2so4": "C 1 O 1 H 1.5 S -0.5 Hsg 0.5",
    "g2ch4": "C 0.5 O 1 Mtg 0.5",
    "ch4o2": "Mtg -1 Oxg -2 C 1 O 4 H 4",
    "fe2o2": "Fe -1 Oxg -0.25 H -3 O -2.5",
    "fe2gr": "Fe -1 H -2 O -2",
    "groxy": "Oxg -0.25 H -1 O -0.5",
    "fe2fs": "Fe -1 Hsg -1 H 1",
    "fesox": "Oxg -2 Fe 1 S 1 O 4",
    "amo2o": "Amm -1 Oxg -2 N 1 O 4 H 3",
    "g2n21": "C 1 None 1 O 1 H 2 N -1",
    "g2n0": "C 1 Nzero 0.8 O 1 H 2 N -0.8",
    "g2n2p": "C 1 O 1 H 2 None -4 Nzero 4",
    "ureah": "C 1 O 2 Amm 2",
    "ch4eb": "Mtg -1",  # CH4 ebullition
    "fwo2": "Oxg -1",  # floodwater exchange (negative: uptake from the atmosphere)
    "fwco2": "C -1 O -2",
    "fwch4": "Mtg -1",
}
OM_REACTIONS = {"g2o2", "g2fe3", "g2gr", "g2gt", "g2so4", "g2ch4", "g2n21", "g2n0", "g2n2p"}
POOLS = ("Om", "Ferrihydrite", "Green_rust", "Goethite_k", "FeS_k", "Urea_k")
GASES = ("Oxg(g)", "Mtg(g)")
OUTPUT_REACTIONS = ("g2ch4", "ch4o2", "ch4eb", "fwch4", "fwco2", "fwo2", "g2n21", "g2n2p")  # dk_ is the day's extent
TOTALS = ("Oxg", "Mtg", "Hsg")

PCH4_ATM = 1.9e-6
SI_NONE = 3  # a gas phase that never forms
PE = 4  # pe of new solutions: O2 is in Oxg, so the coupled O(0) stays negligible
OFF_SWAT = ("includeNkinetics", "includehumus", "includeresidues", "includerootresidues", "includemanure")


def config(ds):
    """``ds['paddy']`` with defaults filled in, or None outside paddy mode."""
    p = ds.get("paddy")
    if p is None:
        return None
    nc, n = ds["ncells"], ds["nrun"]
    z = np.asarray(ds["cell"]["middepth"], float) / 100
    p.setdefault("q", 500 / 365)
    fl = p.get("flooded", np.ones(n))
    if callable(fl):  # e.g. a function of ds['run']['doy']
        fl = fl(ds)
    p["flooded"] = np.asarray(fl, float)[:n] > 0
    p.setdefault("omin", 450 / 12 * np.exp(-z / 0.3) / 0.3)
    p.setdefault("uin", 0.0)
    p.setdefault("n2c", N2C)
    p.setdefault("keb", 1000.0)
    p.setdefault("kfw", 0.1)
    p.setdefault("pools", {})
    for k in ("omin", "uin"):
        if callable(p[k]):  # a function of depth (m)
            p[k] = p[k](z)
        p[k] = np.zeros(nc) + p[k]
    p["pools0"] = {q: np.zeros(nc) + p["pools"].get(q, 0.0) for q in POOLS}
    return p


def apply_hydrology(ds):
    """Saturate the cells on flooded days and percolate ``q`` mm/day out of the column.

    Called after getrichquick: the Richards solution stands on drained days, and the
    water that flooding adds or the drain releases is carried by the cell fluxes, so
    that each cell still ends the day with the next day's water.
    """
    p = config(ds)
    run, cell = ds["run"], ds["cell"]
    nc, n = ds["ncells"], ds["nrun"]
    vsoil = np.asarray(cell["vsoil"], float)[:, None]
    por = np.asarray(cell["por"], float)[:, None]
    fl = p["flooded"]
    old = run["targH2O"].copy()
    new = old.copy()
    new[:, fl] = (por * vsoil + 0 * old)[:, fl]
    d = new - old
    dnext = np.hstack([d[:, 1:], d[:, -1:]])  # change in the water each day ends with
    corr = dnext - d
    inperc, outperc = run["inperc"], run["outperc"]
    infil = np.asarray(run["infil"], float).copy()
    out_new = outperc - np.cumsum(corr, axis=0)
    flnext = np.append(fl[1:], fl[-1])
    extra = np.where(flnext, p["q"] - out_new[-1], 0.0)
    extra = np.fmax(extra, -infil)  # standing water above what the soil takes runs off
    outperc[:] = out_new + extra
    inperc[0] = inperc[0] + extra
    inperc[1:] = outperc[:-1]
    run["infil"] = infil + extra
    run["paddyirr"] = np.fmax(0, extra)
    run["targH2O"] = new
    run["Sws"] = new / vsoil
    if "Swstot" in run:
        run["Swstot"][:, fl] = run["Sws"][:, fl]
    run["airpor"] = np.fmax(0, por - run["Sws"])
    run["flooded"] = np.tile(fl, (nc, 1))
    return ds


def apply_gases(ds):
    """Gas phases for paddy cells: O2 decoupled (Oxg), CH4 vented when drained.

    Flooded cells have no gas phase (si 3, no moles): their gases move by the kinetic
    ebullition and floodwater exchange in Paddy_rates. A gas phase that switched on and
    off inside a day's integration stalled the solver. Drained cells keep ARTEMIS's soil
    CO2 and take O2 as Oxg(g), so Fe(II) and CH4 oxidation stay kinetic; CH4 above the
    atmosphere's leaves.
    """
    config(ds)
    run = ds["run"]
    fl = run["flooded"]
    run["si_CH4g"] = np.where(fl, SI_NONE, np.log10(PCH4_ATM))
    por = np.asarray(ds["cell"]["por"], float)[:, None]
    run["airpor"] = np.fmax(0, por - run["Sws"])  # gasmoles: no gas in flooded cells
    return ds


def co2_phase(ds, k, i):
    """CO2(g) (si, moles) for cell k on day i: ARTEMIS's soil gas when drained, none when flooded."""
    run = ds["run"]
    K, I = k - 1, i - 1
    if run["flooded"][K, I]:
        return SI_NONE, 0.0
    return run["si_CO2g"][K, I], run["CO2gmol"][K, I]


def eqphase_moles(ds, targ, k, i):
    """Oxg(g) moles per kg of the cell's solution on day i (1-based): the air-filled
    pores' O2 in drained cells, none in flooded cells."""
    run = ds["run"]
    K, I = k - 1, i - 1
    return 0.0 if run["flooded"][K, I] else run["O2gmol"][K, I] / targ


def kinetics_parms(ds, k, i, targ):
    """Paddy_rates -parms: f (kg water per m3 bulk), saturation, OM and urea inputs,
    temperature, floodwater weight (top cell, spec 6.2), k_fw (m/yr), thickness (m), ebullition."""
    run, cell, p = ds["run"], ds["cell"], ds["paddy"]
    K, I = k - 1, i - 1
    f = 1000 * run["Sws"][K, I]
    sw = run["Sws"][K, I] / cell["por"][K]
    wfw = min(1, max(0, (sw - 0.9) / 0.1)) if k == 1 else 0
    return (f, sw, p["omin"][K], p["uin"][K], run["soiltemp"][K, I], wfw, p["kfw"] * 365,
            cell["thickness"][K] / 100, p["keb"])


def init_kinetics(w, ds, k, i, targ, tol):
    """The paddy entries of KINETICS k (written inside the cell's block)."""
    p = ds["paddy"]
    w("Paddy_rates\n    -formula H 0\n    -m 1\n    -m0 1\n")
    w("    -parms" + " %.10g" * 9 + "\n", *kinetics_parms(ds, k, i, targ))
    w("    %10s %g \n", "-tol", tol)
    for r, form in FORMULA.items():
        if r in OM_REACTIONS and p["n2c"]:
            form = f"{form} Amm {p['n2c']}"
        w("%s\n    -formula %s\n    -m 1e10\n    -m0 1e10\n", r, form)
    vsoil = ds["cell"]["vsoil"][k - 1]
    for q in POOLS:
        mol = p["pools0"][q][k - 1] * vsoil / 1000 / targ
        w("%s\n    -formula H 0\n    -m %.12g\n    -m0 %.12g\n", q, mol, mol)


def modify_kinetics(w, ds, out, k, i, targ):
    """The paddy entries of KINETICS_MODIFY k: today's parms, and the pools carried over
    in moles per cell and rescaled to today's water."""
    w("    -component Paddy_rates\n")
    w("    %14s" + " %.10g" * 9 + "\n", "-d_params", *kinetics_parms(ds, k, i, targ))
    M = out["run"]["paddy"]
    for q in POOLS:
        mol = M[q][k - 1, i - 2] / targ
        w("    -component %s\n    %14s %.12g\n    %14s %.12g\n", q, "-m", mol, "-m0", mol)


def after_day(run_obj, i):
    """Pools at the end of day i in moles per cell (per m2)."""
    ds, out = run_obj.ds, run_obj.out
    soil = out["soil"]
    M = out["run"].setdefault("paddy", {})
    for q in POOLS:
        M.setdefault(q, np.zeros((ds["ncells"], ds["nrun"])))[:, i - 1] = (
            soil["k_" + q][:, i - 1] * run_obj.targH2O[:, i - 1])
