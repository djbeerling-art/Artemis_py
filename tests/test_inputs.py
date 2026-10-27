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
