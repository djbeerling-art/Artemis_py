"""Input-reading functions against the authors' saved rBcomp ``ds`` struct."""

import os
from functools import lru_cache
from pathlib import Path

import pytest

from artemis_p import matref
from artemis_p.periodictable import get_periodic_table

RBCOMP = "rBcomp_Blueridge_Jan2015_Dec2020_03-Jun-2026"

pytestmark = pytest.mark.skipif(
    not (os.environ.get("ARTEMIS") and os.environ.get("ARTEMIS_PHREEQC")),
    reason="reference data not configured",
)


@lru_cache(maxsize=1)
def ref_ds():
    return matref.load_run(Path(os.environ["ARTEMIS_PHREEQC"]) / f"{RBCOMP}.mat")["ds"]


def test_periodic_table():
    pt, pt2 = get_periodic_table()
    assert matref.differences(pt, ref_ds()["pt"], "pt") == []
    assert pt2["Ca"]["molmass"] == 40.08  # value in periodictable.csv


def test_timetable():
    from artemis_p.timetable import read_timetable

    fname = Path(os.environ["ARTEMIS"]) / "examples" / "EF_data" / "EF_maize_B0_timetable.csv"
    ds = read_timetable({"parm": {}, "nrun": 2192}, fname)
    ref = ref_ds()["timetable"]
    # Later preprocessing steps add more columns to ds.timetable; check those readtimetable sets.
    # setupdatastruc/processtimetable later renames feedstock treatments to ds.parm.fsname.
    ours = {k: v for k, v in ds["timetable"].items() if k != "treatment"}
    assert matref.differences(ours, {k: ref[k] for k in ours}, "timetable") == []


def test_model_species():
    from artemis_p.species import set_model_species

    assert matref.differences(set_model_species(), ref_ds()["solute"], "solute") == []


def test_init_layer_data():
    from artemis_p.layers import read_init_layer_data

    fname = Path(os.environ["ARTEMIS"]) / "examples" / "EF_data" / "initlayers_Flanagan_EF_0-183cm.csv"
    ds = read_init_layer_data({"parm": {}}, fname)
    ref = ref_ds()["layer"]
    assert ds["nlayers"] == ref_ds()["nlayers"]
    # soilsoln is completed later from the initialChemistry file and wcell by
    # getcellsfromlayers; source holds the authors' relative paths.
    skip = {"soilsoln", "wcell", "source"}
    ours = {k: v for k, v in ds["layer"].items() if k not in skip}
    assert matref.differences(ours, {k: v for k, v in ref.items() if k not in skip}, "layer") == []
    assert set(ds["layer"]["source"]) == set(ref["source"]) - {"soilsoln"}


def test_phase_table():
    from artemis_p.phasetable import get_phase_table

    root = Path(os.environ["ARTEMIS"])
    files = [root / "examples" / "EF_data" / "Richland_loess_Batestown_till_phases_combined.csv"] + [
        root / "phreeqcdbincludefiles" / f
        for f in ("secondary_phases.csv", "feedstock_phases_pure.csv", "feedstock_phases_Lewis_THERMODDEM.csv")
    ]
    t = None
    for f in files:
        t = get_phase_table(f, t)
    ref = ref_ds()["phasetable"]
    rows = ref["_rownames"]
    # setupphases later keeps only the phases in use, adds columns and the
    # implicit phase, and marks Titanite/Hydroxyapatite implicit (parm=implicit).
    later = {"Implicit_blueridge"}
    implicit = {"Titanite_blueridge", "Hydroxyapatite_blueridge"}
    assert set(rows) - set(t.index) == later
    for c in t.columns:
        if c == "csvfile":  # authors' paths
            continue
        keep = [i for i, r in enumerate(rows) if r not in later and not (c == "PHREEQCphasetype" and r in implicit)]
        ours = [t.loc[rows[i], c] for i in keep]
        assert matref.differences(ours, [ref[c][i] for i in keep], c) == []
    assert (t.loc[list(implicit), "PHREEQCphasetype"] == "ignore").all()


def test_initial_chemistry():
    from artemis_p.chemistry import read_initial_chemistry
    from artemis_p.layers import read_init_layer_data

    data = Path(os.environ["ARTEMIS"]) / "examples" / "EF_data"
    ds = read_init_layer_data({"parm": {}}, data / "initlayers_Flanagan_EF_0-183cm.csv")
    ds = read_initial_chemistry(ds, data / "rainwater_Champaign_NTN_2021_initialChemistry.csv")
    ref = ref_ds()
    assert matref.differences(ds["init"], ref["init"], "init") == []
    assert matref.differences(ds["influent"], ref["influent"], "influent") == []
    assert matref.differences(ds["layer"]["soilsoln"], ref["layer"]["soilsoln"], "soilsoln") == []


def test_climate_forcing():
    from artemis_p.forcing import read_timeseries, set_run_time_arrays
    from artemis_p.layers import read_init_layer_data

    data = Path(os.environ["ARTEMIS"]) / "examples" / "EF_data"
    ds = read_init_layer_data({"parm": {}}, data / "initlayers_Flanagan_EF_0-183cm.csv")
    ds = set_run_time_arrays(ds, "2015-01-01", "2020-12-31")
    ref = ref_ds()
    run, rrun = ds["run"], ref["run"]
    assert ds["nrun"] == ref["nrun"]
    for k in ("year", "month", "day", "decyear", "doy", "exptday"):
        assert matref.differences(run[k], rrun[k], k) == []
    assert (run["tt"] == rrun["tt"]).all()
    for f in ("soiltemp", "airtemp", "soilresp", "CO2atm", "snow", "infil_ET"):
        ds = read_timeseries(ds, data / f"EF_CLM5_2015-2070_timeseries_{f}.csv")
    forced = ["soiltemp", "airtemp", "soilresp", "hetresp", "rootresp", "CO2atm",
              "snowfall", "snowmelt", "snowpack", "sublimation", "ET", "infil"]
    for k in forced:
        assert matref.differences(run[k], rrun[k], k, rtol=1e-12, atol=1e-12) == [], k
        assert run["units"][k] == rrun["units"][k]
        assert Path(ds["provided"][k]).name == Path(ref["provided"][k]).name


def _ef_phase_table():
    from artemis_p.phasetable import get_phase_table

    root = Path(os.environ["ARTEMIS"])
    files = [root / "examples" / "EF_data" / "Richland_loess_Batestown_till_phases_combined.csv"] + [
        root / "phreeqcdbincludefiles" / f
        for f in ("secondary_phases.csv", "feedstock_phases_pure.csv", "feedstock_phases_Lewis_THERMODDEM.csv")
    ]
    t = None
    for f in files:
        t = get_phase_table(f, t)
    return t


def test_cells_from_layers():
    import numpy as np

    from artemis_p.cells import get_cells_from_layers
    from artemis_p.chemistry import read_initial_chemistry
    from artemis_p.layers import read_init_layer_data

    data = Path(os.environ["ARTEMIS"]) / "examples" / "EF_data"
    ds = {"parm": {"nativeRSAscale": 1e-6, "implicit": ["Titanite_blueridge", "Hydroxyapatite_blueridge"]},
          "phasetable": _ef_phase_table()}
    ds = read_init_layer_data(ds, data / "initlayers_Flanagan_EF_0-183cm.csv")
    ds = read_initial_chemistry(ds, data / "rainwater_Champaign_NTN_2021_initialChemistry.csv")
    ds = get_cells_from_layers(ds)
    ref = ref_ds()
    cell, rcell = ds["cell"], ref["cell"]
    # Changed later: soilsoln (setupphases adds species), CEC (CECsorb converts
    # units), moles/RSA (moved into ds.run and zeroed), labels/units/source.
    later = {"soilsoln", "moles", "RSA", "CEC", "label", "units", "source"}
    for k in cell:
        if k not in later:
            assert matref.differences(cell[k], rcell[k], k, rtol=1e-10) == []
    for ph, v in cell["moles"].items():
        assert np.allclose(v, np.asarray(ref["run"]["moles"][ph])[:, 0], rtol=1e-10), ph
        assert np.allclose(cell["RSA"][ph], np.asarray(ref["run"]["RSA"][ph])[:, 0], rtol=1e-10), ph
    pt, rpt = ds["phasetable"], ref["phasetable"]
    for i, r in enumerate(rpt["_rownames"]):
        if r in pt.index:
            assert pt.loc[r, "PHREEQCphasetype"] == rpt["PHREEQCphasetype"][i], r
