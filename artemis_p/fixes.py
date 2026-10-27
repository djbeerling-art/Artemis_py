"""Switches for fixes to bugs found in the MATLAB code.

The port reproduces the MATLAB code as written by default, so that it can be
checked against the published runs. ``ds['fixes']`` turns fixes on:

* ``"all"`` applies every fix in ``CLEAR`` (coding errors whose intent is plain);
* a collection of names applies just those (names from ``CLEAR`` or ``REVIEW``);
* ``("all", "native_rsa_weighting")`` adds review fixes to the clear ones.

``REVIEW`` fixes change model behaviour in ways the authors should confirm
(they may have tuned parameters around the current behaviour), so ``"all"``
leaves them out. docs/PORTING_PLAN.md describes each bug.
"""

CLEAR = {
    # setup and forcing
    "leap_year": "getdecyear: 2000 is a leap year",
    "depth_units": "checktsdepths and observation depths: convert m and mm to cm",
    "eqcec_units": "CECsorb: eqCEC from CEC in eq/g, not meq/100g",
    "density_units": "setupphases: rho_kgm3 = G*1000 for g/cm3",
    "geo_from_lambda": "setRSA: geometric area = BET/lambda",
    "hornblende_formula": "HornblendeRLBT formula Ca2Mg5Si8O22(OH)2 (input file has Mg52)",
    # uptake
    "crop_case": "setupnutuptake and daySWAT7: match crop names case-insensitively",
    "p_uptake_all_cells": "daySWAT7 getcellup: rescale P uptake in every cell",
    "amm_uptake_mm_row": "Ammonium uptake uses its own Michaelis-Menten row (N), not the last nutrient's",
    "passive_uptake_flag": "Passive uptake tests each nutrient's trnutflag, not the cell index",
    # PHREEQC input and daily loop
    "pctorgc_update": "Update percent organic C for denitrification each day (out.soil.k_Nhuma)",
    "influent_redox": "Influent pe and O2(g) from ds.influent are used",
    "fractal_per_phase": "updateRSA: use per-phase fractal exponents",
    "leachate_ph": "No-infiltration days set the influent pH to NaN (not a stray out.pH)",
    "din_uptake": "k_DINuptake is nitrate plus ammonium uptake during the run",
    # hydrology
    "frost_table": "WaterTable: set k_frz (not kfrz)",
    "clm5_layer_overlap": "getds: overlap fix and initial water for the first added CLM5 layer",
    # end-of-run totals
    "cdr_total": "cCDR = cCDRrh + cumulative change in soil carbonate",
    "native_release": "Record release from native minerals in c<elt>natrelease",
    "release_units": "Convert feedstock Ca/Mg/Na/K release to tCO2/ha once",
    "feedstock_remaining": "Remaining feedstock totals counted once",
    "uptake_totals": "dt_<nut>uptake = daily passive + daily active uptake",
    "cec_ratio_cells": "Per-cell release/CEC ratios use each cell's release",
}

REVIEW = {
    "native_rsa_weighting": "setupsoilphases: weight native RSA by mineral proportion (nativeRSAscale was tuned without it)",
    "upflow_fractions": "setmixfracs2 top-down: upflow fraction relative to the source cell's water",
    "day1_soil_moisture": "getrichquick: fill Sws on day 1",
    "bottom_layer_inflow": "soilwater_zengdecker2009: bottom layer keeps its inflow. NOT VALID ALONE: in rBcomp the water table rises to the surface and drainage stops from 2019",
    "nfix_flag": "includeNfix=0 (the default) removes Nfix from the rates and KINETICS, so legumes stop fixing N",
}


def active(ds):
    """Names of the fixes on for this run."""
    return sorted(n for n in list(CLEAR) + list(REVIEW) if enabled(ds, n))


def enabled(ds, name):
    """True if fix ``name`` is on for this run."""
    if name not in CLEAR and name not in REVIEW:
        raise KeyError(f"unknown fix {name!r}")
    f = ds.get("fixes", ()) if isinstance(ds, dict) else ds
    if isinstance(f, str):
        f = (f,)
    return name in f or ("all" in f and name in CLEAR)
