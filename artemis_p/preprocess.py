"""The preprocessor. Python counterpart of setupdatastruc.m and readinputfiles.m.

``setup_data_struc`` runs the MATLAB steps in the same order, up to the end of
``setrunarrays``. The hydrology and nutrient-uptake setup that follows
(setuprun, setupnutuptake) is stage 4 of the port.

Input files are given as in MATLAB's ``ds.usefiles``. ``resolve`` maps each
name to a local path (for the published runs: relative to examples/EF2, with
the authors' include directory replaced).
"""

import copy
from pathlib import Path

import numpy as np

from .cells import get_cells_from_layers
from .chemistry import read_initial_chemistry
from .forcing import read_timeseries, set_run_time_arrays
from .layers import read_init_layer_data
from .periodictable import get_periodic_table
from .phases import setup_phases
from .obs import read_obs_data, set_run_cell_indices
from .phasetable import get_phase_table
from .sorption import cec_sorb, set_soil_texture
from .species import set_model_species
from .timetable import read_timetable
from .treatments import process_timetable


def read_input_files(ds, files, resolve=Path):
    """readinputfiles.m for initlayers, phases, timetable, timeseries, PSD and
    initialChemistry and observation files.

    File names are recorded as given (in ``inputfiles``, ``provided``, ...),
    not as resolved.
    """
    parm = ds["parm"]
    ds["inputfiles"] = []

    def read(fn, f, *args, **kw):
        n = len(ds["inputfiles"])
        out = fn(ds, resolve(f), *args, **kw)
        out["inputfiles"][n:] = [f] * (len(out["inputfiles"]) - n)
        src = out.get("layer", {}).get("source", {})
        for k, v in src.items():
            if v == str(resolve(f)):
                src[k] = f
        return out

    phases = [f for f in files if "phases" in f]
    for f in files:
        if "parameters" in f:
            raise NotImplementedError("parameters files are not ported yet")
        if "obs" in f:
            continue
        if "initlayers" in f:
            ds = read(read_init_layer_data, f)
        elif "phases" in f:
            ds["phasetable"] = get_phase_table(resolve(f), None if f == phases[0] else ds["phasetable"])
            ds["inputfiles"].append(f)
    for f in files:
        if "timetable" in f:
            ds = read(read_timetable, f)
    for f in files:
        if "timeseries" in f:
            ds = read(read_timeseries, f)
            for k, v in ds["provided"].items():
                if v == str(resolve(f)):
                    ds["provided"][k] = f
        elif "PSD" in f:
            if "feedstockPSDfile" not in ds:
                ds["feedstockPSDfile"] = [f]
            else:
                raise NotImplementedError("more than one PSD file is not ported yet")
        elif "initialChemistry" in f:
            ds = read(read_initial_chemistry, f, influentfld=parm.get("influentfld", "rainwater"),
                      soilsolnfld=parm.get("soilsolnfld", "porewater"))
            ds["layer"]["source"]["soilsoln"] = f
    for f in files:  # observations last; they need run indices and/or cell indices assigned later
        if "obs_" in f:
            ds = read_obs_data(ds, resolve(f), name=f)
    return ds


def setup_data_struc(ds, files, resolve=Path, basedir=None):
    """setupdatastruc.m up to the end of setrunarrays.

    ``ds`` holds the settings (defaultsettings.m and the run's arguments):
    ``parm`` (including ``t1``/``t2``), ``dbdir``, ``EXCHmethod`` and so on.
    ``basedir`` resolves the relative PSD file name.
    """
    ds = copy.deepcopy(ds)
    ds.setdefault("datastep", 1)
    ds["pt"], ds["pt2"] = get_periodic_table()
    ds["solute"] = set_model_species(ds["pt2"])
    ds.setdefault("scale", {"description": "Scale factors for cell or run arrays"})
    ds.setdefault("provided", {"description": "Fields for which input data were provided"})
    parm = ds["parm"]
    parm.setdefault("usemeanforcings", 0)
    if "t1" in parm and "t2" in parm:
        ds = set_run_time_arrays(ds, parm["t1"], parm["t2"])
    else:
        raise NotImplementedError("runs without parm.t1 and parm.t2 are not ported yet")

    ds = read_input_files(ds, files, resolve)
    nrun = ds["nrun"]

    if "layer" in ds:
        ds = get_cells_from_layers(ds)
    else:
        print("(setupdatastruc) No layer data provided")
    if "cell" not in ds:
        print("(setupdatastruc) No PHREEQC cells; quitting")
        return ds
    ds = set_soil_texture(ds)
    ds = set_run_cell_indices(ds)
    if "exchange" in ds["cell"] or "available" in ds["cell"]:
        raise NotImplementedError("CECsorb with existing exchange/available data is not ported yet")
    ds = cec_sorb(ds)

    if "timetable" not in ds:
        raise NotImplementedError("maketimetable is not ported yet")
    ds = process_timetable(ds)
    if parm.get("writephasesratesfile") or ds.get("writephasesratesfile"):
        raise NotImplementedError("writePHASESRATES is not ported yet")
    ds = setup_phases(ds, basedir=basedir)
    _set_feedstock(ds)
    _set_run_arrays(ds)
    soilsoln = ds["cell"].get("soilsoln")
    if soilsoln is not None and "gsoil" in soilsoln["units"]["msp"]:
        raise NotImplementedError("soil solutions per gram of soil are not ported yet")
    ds["nmix"] = ds["ncells"] + 1
    ds["nmixsoln"] = ds["ncells"] + 2
    ds["run"]["mixfrac"] = np.zeros((ds["nmix"], ds["nmixsoln"], nrun))
    ds["MIXture"] = [f"MIX {i}" for i in range(1, ds["nmix"] + 1)]
    ds["MIXsoln"] = np.arange(ds["nmixsoln"], dtype=float)
    ds["MIXsolnstr"] = [f"soln {i}" for i in range(ds["nmixsoln"])]
    return ds


def _set_feedstock(ds):
    """setfeedstock: dose in kg/m2, feedstock name per timetable entry, treated days."""
    tt, run = ds["timetable"], ds["run"]
    if "baskg" in tt and "dosekg" not in tt:
        tt["dosekg"] = np.array(tt["baskg"], float)
    if "dosekg" not in tt and "dosage" in tt:
        # every entry with recognised units, not just feedstocks (as in MATLAB)
        print("(setupdatastruc/setfeedstock) creating dosekg from dosage")
        dosage = np.asarray(tt["dosage"], float)
        units = [str(u).lower() for u in tt["doseunits"]]
        scale = {"kg/ha": 1e-4, "kg/m2": 1.0, "t/ha": 0.1}
        tt["dosekg"] = np.array([dosage[j] * scale[u] if u in scale else dosage[j] * 0
                                 for j, u in enumerate(units)])
    tt["dosekg"] = np.where(np.isnan(tt["dosekg"]), 0, tt["dosekg"])
    run["treated"] = np.zeros(ds["nrun"])
    if "treatment" in tt:
        tt["feedstock"] = [""] * len(tt["decyear"])
        if "doseform" in tt:
            for j, form in enumerate(tt["doseform"]):
                if form.lower() == "feedstock" and tt["runindex"][j] > 0:
                    tt["feedstock"][j] = tt["treatment"][j]
                    run["treated"][int(tt["runindex"][j]) - 1] = 1
    run["tillagefrac"] = np.zeros((ds["ncells"], ds["nrun"])) + np.asarray(ds["cell"]["tillagefrac"])[:, None]
    # markPSDinit (particle size tracking for SAmethod 'discrete'). It matches
    # timetable treatments to feedstocks case-sensitively, so for EF ('Blueridge'
    # vs 'blueridge') it builds nothing and ds.mark holds only its defaults.
    ptype = ds["phasetable"]["PHREEQCphasetype"].astype(str).str.lower()
    if ptype.str.contains("kin").any():
        names = set(ds["feedstocktable"]["feedstock"]) if "feedstocktable" in ds else set()
        if any(t in names for t in tt.get("treatment", [])):
            # the PSD structures need Mark Lomas's routines (PSD_range_2_P, ...), which are not
            # in the archive; setupdatastruc catches the error and leaves ds.mark unset. Only
            # SAmethod 'discrete' (updateRSA) uses them.
            if str(ds.get("SAmethod", "")).lower() == "discrete":
                raise NotImplementedError("markPSDinit PSD tracking (SAmethod 'discrete') is not ported yet")
            print("****Unable to set up Lomas PSD-tracking entities****")
            print("****Avoid SAmethod=\"discrete\"****")
        else:
            ds["mark"] = {"binmethod": "SA-mean", "binspec": 10, "wholebasaltflag": 0}


def _set_crop_from_timetable(ds):
    """setcropfromtimetable: crop, planting and plant-growth arrays per run day.

    Timetable decimal years use the day's midpoint and run decimal years its
    start, so a crop entry takes effect from the following day.
    """
    tt, run, parm = ds["timetable"], ds["run"], ds["parm"]
    if "crop" not in tt:
        return
    n = ds["nrun"]
    run["crop"] = ["fallow"] * n
    for k in ("wtt", "plant", "sowdepth", "harvest", "peaktot", "peakroot"):
        run[k] = np.zeros(n)
    run["units"].update(harvest="g m^{-2}", sowdepth="cm", peaktot="g m^{-2}", peakroot="g m^{-2}")
    if "crop" in parm:
        raise NotImplementedError("parm.crop with a timetable crop is not ported (MATLAB assigns the whole crop array)")
    ri = np.asarray(tt["runindex"])
    for f in ("peaktot", "peakroot", "sowdepth", "harvest"):
        if f not in tt:
            continue
        v = np.asarray(tt[f], float)
        w = np.flatnonzero((v > 0) & (ri > 0))
        run[f] = np.zeros(n)
        run[f][ri[w].astype(int) - 1] = v[w]
    crop = "fallow"
    for itt, c in enumerate(tt["crop"]):
        w = np.flatnonzero(run["decyear"] >= tt["decyear"][itt])
        if c.lower() != "fallow" and len(c) > 0:
            crop = c
            for i in w:
                run["crop"][i] = crop
            run["plant"][w] = 1
            if "tillagedepth" in tt and tt["tillagedepth"][itt] > 0 and len(w) > 1:
                w = w[1:]
                run["plant"][w] = 0
                for i in w:
                    run["crop"][i] = "fallow"
        elif c.lower() == "fallow":
            for i in w:
                run["crop"][i] = "fallow"
            run["plant"][w] = 0
        if crop.lower() == "fallow":
            run["plant"][w] = 0
    if "plantdensity" not in ds and np.sum(run["peaktot"]) > 0:
        ds["plantdensity"] = 1


def _set_run_arrays(ds):
    """setrunarrays: crop, residue, weather defaults and other per-day arrays."""
    run, parm, cell = ds["run"], ds["parm"], ds["cell"]
    n, nc = ds["nrun"], ds["ncells"]
    runzero = np.zeros(n)
    units, label = run.setdefault("units", {}), run.setdefault("label", {})
    ds["nut"] = []
    run["crop"] = ["fallow"] * n
    run["plant"] = runzero.copy()
    if "crop" in parm:
        run["crop"] = [parm["crop"]] * n
        run["plant"][:] = 1
    else:
        _set_crop_from_timetable(ds)
    run["residue"] = np.zeros((nc, n))
    label["residue"], units["residue"] = "Residue tilled into soil", "g cell^{-1}"
    run["rootresidue"] = np.zeros((nc, n))
    label["rootresidue"], units["rootresidue"] = "Root residue", "g cell^{-1}"
    if "ET" not in run:
        if "soilevap" in run and "transpiration" in run:
            run["ET"] = run["soilevap"] + run["transpiration"]
        elif "potET" in run:
            print("(setupdatastruc/setrunarrays) No ET given; using given potential ET")
            run["ET"] = run["potET"]
    if "precip" not in run:
        run["precip"] = run["exptday"] * 0 + parm.get("precip", 5)  # 5 mm/day default
        units["precip"] = "mm"
    if "airtemp" not in run:
        if "airtemp" in parm:
            run["airtemp"] = run["exptday"] * 0 + parm["airtemp"]
        elif "soiltemp" in run:
            run["airtemp"] = np.asarray(run["soiltemp"])[0].copy()
        else:
            run["airtemp"] = run["exptday"] * 0 + 25
    if "targH2O" not in cell:
        cell["targH2O"] = np.asarray(cell["FC"], float).copy()
        cell["units"]["targH2O"] = "kgw cell^{-1}"
        cell["label"]["targH2O"] = "Mean Day 0 cell moisture"
    run["targH2O"] = np.zeros((nc, n)) + np.asarray(cell["targH2O"])[:, None]
    if "CO2atm" in parm:
        run["CO2atm"] = runzero + parm["CO2atm"]
    if "CO2atm" in run:
        run["si_CO2atm"] = np.log10(np.asarray(run["CO2atm"]) / 1e6)
        units["CO2atm"], label["CO2atm"] = "ppmV", "Atmospheric CO_2"
        units["si_CO2atm"] = "log_10(atm)"
    if "soiltemp" not in run:
        run["soiltemp"] = np.tile(run["airtemp"], (nc, 1))
    units["airtemp"], label["airtemp"] = "^{\\circ}C", "Air temperature at soil surface"
    run.setdefault("harvest", runzero.copy())
    if "soilresp" not in run:
        run["soilresp"] = runzero.copy()
        if "soilrespmolCm2day" in parm:
            run["soilresp"] += parm["soilrespmolCm2day"]
        elif "soilrespmumolCm2s" in parm:
            run["soilresp"] += parm["soilrespmumolCm2s"] * 86400 / 1e6
        elif "soilrespkgCO2hahr" in parm:
            run["soilresp"] += parm["soilrespkgCO2hahr"] * 24 / 44.01 * 1000 / 10000
    if "peaktot" in parm:
        run["peaktot"] = runzero + parm["peaktot"]
    if "totmassgm2" in parm:
        run["totmass"] = runzero + parm["totmassgm2"]
    elif "totmasstha" in parm:
        run["totmass"] = runzero + parm["totmasstha"] / 1e4 * 1e6
    if "zrootcm" in parm:
        run["zroot"] = runzero + parm["zrootcm"]
    elif "zrootm" in parm:
        run["zroot"] = runzero + parm["zrootm"] * 100
    if "LAI" in parm:
        run["LAI"] = runzero + parm["LAI"]
