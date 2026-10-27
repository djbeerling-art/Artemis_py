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
