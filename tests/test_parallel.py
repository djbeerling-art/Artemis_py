"""Cells' kinetics run in parallel (phreeqc.run_cells_parallel) agree with a serial run."""

import copy
import os
from pathlib import Path

import numpy as np
import pytest

from artemis_p import phreeqc

needs_data = pytest.mark.skipif(
    not all(os.environ.get(v) for v in ("ARTEMIS", "ARTEMIS_PHREEQC", "PHREEQC")),
    reason="reference data or PHREEQC not configured",
)


def _run(tmp, ndays, **opts):
    from test_inputs import _rbcomp_setup

    from artemis_p.phreeqc_input import PhreeqcRun
    from artemis_p.reference import rbcomp_run_settings, stage_ef2_dir

    artemis = Path(os.environ["ARTEMIS"])
    ef2 = stage_ef2_dir(tmp, artemis, os.environ.get("ARTEMIS_FLANAGAN"))
    ds = copy.deepcopy(_rbcomp_setup())
    ds.update(rbcomp_run_settings(dbdir=f"{artemis / 'phreeqcdbincludefiles'}/", outdir=f"{ef2}/"))
    r = PhreeqcRun(ds, **opts)
    for i in range(1, ndays + 1):
        assert r.run_day(i).ok
    return r, ef2


@needs_data
def test_parallel_cells_match_serial(tmp_path):
    ndays = 3
    a, _ = _run(tmp_path / "serial", ndays)
    b, _ = _run(tmp_path / "parallel", ndays, cellworkers=4)
    parts = phreeqc._split_day(Path(b.fname + ".in").read_text())
    assert parts is not None and [k for k, _ in parts[3]] == list(range(8, 0, -1))
    # cell 2 (frozen, 0.01 mm of water) can need smaller pe steps when started from a dump
    assert all(b.out["run"]["nwarnings"][i] >= n for i, n in a.out["run"]["nwarnings"].items())
    # Anatase flips between present and absent with rounding; charge and pct_err are ~1e-8 residuals.
    # Cell 2 holds 0.01 mm of (frozen) water, so its concentrations are ill-conditioned; the
    # daily gas exchanges d_CO2g and d_O2g are small differences of large amounts.
    skip = {"Anatase", "d_Anatase", "AnataseRLBT", "d_AnataseRLBT", "si_Anatase", "si_AnataseRLBT", "si_Goethite",
            "charge", "pct_err"}
    for f, x in a.out["soil"].items():
        if f in skip or not isinstance(x, np.ndarray) or x.dtype.kind != "f":
            continue
        x, y = x[:, :ndays], b.out["soil"][f][:, :ndays]
        scale = np.nanmax(np.abs(x), initial=0)
        rtol = 1e-2 if f in ("d_CO2g", "d_O2g") else 1e-4
        for cells, tol in ((np.r_[0, 2:len(x)], rtol), (1, max(10 * rtol, 1e-3))):
            assert np.allclose(x[cells], y[cells], rtol=tol, atol=1e-6 * scale, equal_nan=True), f
