"""The in-memory IPhreeqc engine (artemis_p.iphreeqc) agrees with the PHREEQC executable."""

import ctypes
import os

import numpy as np
import pytest

from test_parallel import _run, needs_data


def _have_lib():
    try:
        ctypes.CDLL(os.environ.get("IPHREEQC_LIB") or "libIPhreeqc.so")
    except OSError:
        return False
    return True


needs_lib = pytest.mark.skipif(not _have_lib(), reason="libIPhreeqc not found (set IPHREEQC_LIB)")

# Anatase flips between present and absent with rounding; charge and pct_err are ~1e-8 residuals
SKIP = {"Anatase", "d_Anatase", "AnataseRLBT", "d_AnataseRLBT", "si_Anatase", "si_AnataseRLBT", "si_Goethite",
        "charge", "pct_err"}


def _fields(r, ndays):
    return {f: x[:, :ndays] for f, x in r.out["soil"].items()
            if f not in SKIP and isinstance(x, np.ndarray) and x.dtype.kind == "f"}


@needs_data
@needs_lib
def test_iphreeqc_matches_executable(tmp_path):
    ndays = 4
    exe, _ = _run(tmp_path / "exe", ndays)
    one, _ = _run(tmp_path / "one", ndays, engine="iphreeqc", printoutput=0)
    par, _ = _run(tmp_path / "par", ndays, engine="iphreeqc", printoutput=0, cellworkers=3)
    a, b, c = _fields(exe, ndays), _fields(one, ndays), _fields(par, ndays)
    for f, x in a.items():
        # the executable passes the state between days through a text dump; small differences
        # grow most in cells with little water
        scale = np.nanmax(np.abs(x), initial=0)
        assert np.allclose(x, b[f], rtol=0, atol=1e-3 * scale, equal_nan=True), f
        assert np.allclose(x, c[f], rtol=0, atol=1e-3 * scale, equal_nan=True), f
    # each cell has its own instance, so results don't depend on the number of workers
    four, _ = _run(tmp_path / "four", ndays, engine="iphreeqc", printoutput=0, cellworkers=4)
    for f, x in _fields(four, ndays).items():
        assert np.array_equal(x, c[f], equal_nan=True), f
    assert exe.out["run"]["nwarnings"] == one.out["run"]["nwarnings"]


def test_split_dump():
    from artemis_p.iphreeqc import _split_dump, _water

    d = ("SOLUTION_RAW 2 x\n  -total_h 1.5\nEQUILIBRIUM_PHASES_RAW 2\n  -x 1\n"
         "SOLUTION_RAW 10\n  -total_h 900\nKINETICS_RAW 10\n")
    s = _split_dump(d)
    assert sorted(s) == [2, 10]
    assert s[2].startswith("SOLUTION_RAW 2") and "EQUILIBRIUM_PHASES_RAW 2" in s[2]
    assert _water(s[2]) == 1.5 and _water(s[10]) == 900 and _water("") == float("inf")
