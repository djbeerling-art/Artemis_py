"""Regression tests against the published ARTEMIS v1.1 PHREEQC files.

Needs:
  ARTEMIS          v1.1 soilmodel directory
  ARTEMIS_PHREEQC  directory with the published PHREEQCfiles_GMD_R1 files
  ARTEMIS_FLANAGAN path to Flanagan_BlueRidge_PHASESRATES_withimplicitphases
                   (missing from v1.1; take it from v1.0 examples/EF_CLM5)
  PHREEQC          PHREEQC 3.7.3 executable (the version the authors used)
"""

import os
from pathlib import Path

import pytest

from artemis_p import matref, reference

RBCOMP = "rBcomp_Blueridge_Jan2015_Dec2020_03-Jun-2026"

needs_data = pytest.mark.skipif(
    not all(os.environ.get(v) for v in ("ARTEMIS", "ARTEMIS_PHREEQC", "PHREEQC")),
    reason="reference data or PHREEQC not configured",
)


@needs_data
def test_rbcomp_first_iteration_matches_published(tmp_path):
    res, cmp = reference.rerun_init(
        RBCOMP, os.environ["ARTEMIS_PHREEQC"], tmp_path,
        flanagan=os.environ.get("ARTEMIS_FLANAGAN"),
    )
    assert res.ok, res.errors
    # Anatase sits at saturation and flips between present and absent with
    # rounding, so its amounts are excluded; its SI still has to match.
    flips = {"Anatase", "d_Anatase", "AnataseRLBT", "d_AnataseRLBT"}
    close = cmp.drop(index=[c for c in flips if c in cmp.index])
    assert close.loc["pH", "max_abs"] < 1e-5
    assert close.loc["pe", "max_abs"] < 1e-4
    assert close["max_abs"].max() < 1e-3


@needs_data
def test_rbcomp_mat_has_preprocessor_struct():
    path = Path(os.environ["ARTEMIS_PHREEQC"]) / f"{RBCOMP}.mat"
    run = matref.load_run(path)
    assert {"ds", "metadata", "soil"} <= set(run)
    assert run["ds"]["EXCHmethod"] == "comp"
