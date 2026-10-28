"""Soil gas diffusion: a daily 1-D diffusion step for CO2 (and optionally O2) between cells.

Turned on by ``ds['gasdiff']``, a dict with

``gases``   the gases to diffuse, from ``("CO2", "O2")`` (default ``("CO2",)``)
``source``  respiration that enters the gas: ``"soilresp"`` (default, the prescribed total,
            as getsoilCO2g) or ``"rootresp"`` (roots only, for runs where decomposition
            releases its own CO2 in PHREEQC, ``releaseCO2`` = 1)
``top``     ``"open"`` (default, the atmosphere above the top cell) or ``"closed"`` (no flux)
``model``   gas-phase tortuosity: ``"MQ"`` (default, Millington and Quirk 1961,
            ta^(10/3)/phi^2), ``"Moldrup"`` (Moldrup et al. 2000, ta^2.5/phi, for undisturbed
            soils; less steep in wet soil) or ``"SCEPTER"`` (SCEPTER's default, Archie-type
            ta^3.4/phi, for the benchmark)
``icegap``  air-filled porosity left open through ice (default 0.02): frozen soil is not
            sealed, but the ice-filled pores otherwise count as closed to gas
``feedback`` 1 (default) to add the gas each cell's chemistry takes up or releases
            (PHREEQC's ``d_CO2g``, ``d_O2g``) to the sources, as a running mean over ``tau``
            days (default 30) that leaves out the first day, when the starting solutions
            equilibrate. Fed back day by day, the uptake made the profile oscillate: a
            calcite or exchange buffer answers a pressure change at once, so yesterday's
            uptake overshoots today's. The running mean keeps the sustained sink (weathering,
            decomposition) and averages out the buffer's swings.

This replaces Cerling's (1991) analytic steady state (getsoilCO2g) with a backward-Euler
finite-volume step of one day on the cell grid. The effective diffusivity is
Millington and Quirk (1961) in both phases,

    D = Da(T) ta^(10/3)/phi^2 + beta(T) Dw(T) tw^(10/3)/phi^2

and the storage is (ta + beta tw) per unit soil volume, where ta and tw are the air- and
water-filled porosities (ice fills pores, so ta is from the total water), phi the porosity and beta the dimensionless Henry solubility, so
concentrations are gas-equivalent (mol per cm3 of air). Respiration is spread over the
depth with Cerling's exponential profile (cdepth) above the impermeable depth. Each day's
partial pressure is passed to PHREEQC as the CO2(g) (O2(g)) saturation index; the gas
moles stay as set by gasmoles (``fixCO2gmol``), so within a day the gas phase holds the
diffused pressure and the chemistry's draw-down enters the next day's step.
"""

import numpy as np

R = 82.057338  # cm3 atm mol-1 K-1
DAIR0 = {"CO2": 0.144, "O2": 0.178}  # cm2 s-1 at 273.15 K and 1 atm (getsoilCO2g's CO2 value)
DEXP = 1.823
DW25 = {"CO2": 1.91e-5, "O2": 2.10e-5}  # cm2 s-1 in water at 25 C
KH25 = {"CO2": 3.4e-2, "O2": 1.3e-3}  # mol L-1 atm-1 at 25 C
KHT = {"CO2": 2400.0, "O2": 1700.0}  # van 't Hoff d(ln KH)/d(1/T), K (Sander 2015)
PMIN = 1e-8  # atm: a cell that the chemistry would drive below zero holds this
SEC = 86400.0


def config(ds):
    """``ds['gasdiff']`` with defaults filled in, or None when gas diffusion is off."""
    g = ds.get("gasdiff")
    if not g:
        return None
    if g is True or g == 1:
        g = ds["gasdiff"] = {}
    g.setdefault("gases", ("CO2",))
    g.setdefault("source", "soilresp")
    g.setdefault("top", "open")
    g.setdefault("feedback", 1)
    g.setdefault("tau", 30.0)
    g.setdefault("model", "MQ")
    g.setdefault("icegap", 0.02)
    return g


def _water_viscosity(TK):
    """Pa s (Vogel equation)."""
    return 2.414e-5 * 10 ** (247.8 / (TK - 140))


def properties(gas, TK, ta, tw, por, model="MQ"):
    """Effective diffusivity (cm2 s-1) and gas-equivalent capacity (-) per cell."""
    beta = KH25[gas] * np.exp(KHT[gas] * (1 / TK - 1 / 298.15)) * R / 1000 * TK
    da = DAIR0[gas] * (TK / 273.15) ** DEXP
    dw = DW25[gas] * TK / 298.15 * _water_viscosity(298.15) / _water_viscosity(TK)
    ta, tw = np.fmax(ta, 0), np.fmax(tw, 0)
    fa = {"moldrup": ta ** 2.5 / por, "scepter": ta ** 3.4 / por}.get(model.lower(), ta ** (10 / 3) / por ** 2)
    D = da * fa + beta * dw * tw ** (10 / 3) / por ** 2
    return D, ta + beta * tw


def respiration_profile(ds, I, source="soilresp"):
    """Respiration entering each cell on day I (0-based), mol cm-2 s-1 of land.

    Cerling's source phi0 exp(-z/cdepth) between the surface and the impermeable depth,
    integrated over each cell, so that the column total is the day's respiration."""
    run, cell, parm = ds["run"], ds["cell"], ds["parm"]
    total = float(np.asarray(run[source], float)[I]) / 1e4 / SEC
    if total <= 0:
        return np.zeros(ds["ncells"])
    if source == "rootresp":
        c = float(np.asarray(run["tcdepth"], float)[I])
        c = c if c > 0 else float(np.asarray(run["cdepth"], float)[I])
    else:
        c = parm["cdepth"] if "cdepth" in parm else float(np.asarray(run["cdepth"], float)[I])
    L = float(np.asarray(run["impermeabledepth"], float)[I])
    top = np.fmin(np.asarray(cell["topdepth"], float), L)
    bot = np.fmin(np.asarray(cell["bottomdepth"], float), L)
    w = np.exp(-top / c) - np.exp(-bot / c)
    return total * w / w.sum() if w.sum() > 0 else np.zeros(ds["ncells"])


def step(C0, D, cap, dz, S, Ctop, dt=SEC, closed=False):
    """One backward-Euler step: concentrations (mol cm-3 of air-equivalent) at the end.

    dz in cm, S in mol cm-2 s-1 per cell (positive into the gas), Ctop the atmosphere."""
    n = len(C0)
    with np.errstate(divide="ignore"):
        half = np.where(D > 0, dz / (2 * D), np.inf)  # s cm-1, resistance of each half cell
    g = 1 / (half[:-1] + half[1:])  # conductance between neighbours
    g0 = 0.0 if closed else (1 / half[0] if np.isfinite(half[0]) else 0.0)
    store = cap * dz / dt
    A = np.diag(store.copy())
    A[0, 0] += g0
    for j in range(n - 1):
        A[j, j] += g[j]
        A[j + 1, j + 1] += g[j]
        A[j, j + 1] -= g[j]
        A[j + 1, j] -= g[j]
    b = store * C0 + S
    b[0] += g0 * Ctop
    # a cell with no air and no water (frozen solid) has no storage or conductance
    dead = np.diag(A) <= 0
    A[dead, dead] = 1
    b[dead] = C0[dead]
    return np.linalg.solve(A, b)


def _atmosphere(ds, gas, I):
    run = ds["run"]
    if gas == "O2":
        return 0.20946
    units = run["units"].get("CO2atm", "ppmV")
    a = float(np.asarray(run["CO2atm"], float).reshape(-1)[I])
    return a if units.lower() == "atm" else a / 1e6


def update_day(ds, out, i, targ, fixmol=None):
    """Set day i's (1-based) soil pCO2 (pO2), si and, for closed runs, gas moles.

    ``targ`` is PhreeqcRun.targH2O (PHREEQC moles times targ = moles per cell);
    ``fixmol`` maps a gas to its fixed gas moles (fixCO2gmol, fixO2gmol; 0 = from PV = nRT)."""
    g = config(ds)
    run, cell = ds["run"], ds["cell"]
    I = i - 1
    TK = 273.15 + np.asarray(run["soiltemp"], float)[:, I]
    por = np.asarray(cell["por"], float)
    tw = np.asarray(run["Sws"], float)[:, I]
    tot = np.asarray(run.get("Swstot", run["Sws"]), float)[:, I]
    ta = np.fmax(0, por - tot)
    ta = np.fmax(ta, np.fmin(g["icegap"], tot - tw))  # a little air through the ice
    dz = np.asarray(cell["thickness"], float)
    resp = respiration_profile(ds, I, g["source"])
    state = run.setdefault("gasdiff", {})
    for gas in g["gases"]:
        key = "p" + gas + "g"
        D, cap = properties(gas, TK, ta, tw, por, g["model"])
        S = resp if gas == "CO2" else -resp  # respiration quotient 1
        mean = state.setdefault("S" + gas, np.zeros(len(dz)))
        if g["feedback"] and i > 2 and ("d_" + gas + "g") in out["soil"]:
            dchem = out["soil"]["d_" + gas + "g"][:, I - 1] * np.asarray(targ)[:, I - 1]  # mol per cell (m2)
            mean += (np.nan_to_num(dchem) / 1e4 / SEC - mean) / g["tau"]
        S = S + mean
        patm = _atmosphere(ds, gas, I)
        if gas not in state:  # start from the atmosphere and let the first day settle
            C0 = np.full(len(dz), patm) / (R * TK)
            for _ in range(30):
                C0 = step(C0, D, cap, dz, S, patm / (R * TK[0]), closed=g["top"] == "closed")
        else:
            C0 = state[gas] / (R * TK)  # yesterday's pressure at today's temperature
        C = step(C0, D, cap, dz, S, patm / (R * TK[0]), closed=g["top"] == "closed")
        p = np.fmax(C * R * TK, PMIN)
        state[gas] = p
        run[key][:, I] = p
        with np.errstate(divide="ignore"):
            run["si_" + gas + "g"][:, I] = np.log10(p)
        if gas == "CO2":
            run["soilCO2gppmV"][:, I] = p * 1e6
        if (fixmol or {}).get(gas, 1) == 0 and (gas + "gmol") in run:
            airvol = run["airpor"][:, I] * np.asarray(cell["vsoil"], float)
            run[gas + "gmol"][:, I] = p * airvol / 0.08207 / TK
    return ds
