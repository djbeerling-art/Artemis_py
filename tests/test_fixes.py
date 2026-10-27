"""Fixes to MATLAB bugs (artemis_p.fixes): each switch changes what it should.

The legacy behaviour (no ``ds['fixes']``) is checked against the published run
by the other tests.
"""

import os

import numpy as np
import pytest

from artemis_p import fixes

needs_data = pytest.mark.skipif(not os.environ.get("ARTEMIS_PHREEQC"), reason="needs ARTEMIS_PHREEQC")


def test_registry():
    assert not fixes.enabled({}, "cdr_total")
    assert fixes.enabled({"fixes": "all"}, "cdr_total")
    assert not fixes.enabled({"fixes": "all"}, "native_rsa_weighting")
    assert fixes.enabled({"fixes": ("all", "native_rsa_weighting")}, "native_rsa_weighting")
    assert fixes.enabled({"fixes": ["crop_case"]}, "crop_case")
    assert not set(fixes.CLEAR) & set(fixes.REVIEW)
    with pytest.raises(KeyError):
        fixes.enabled({}, "no_such_fix")
    assert fixes.active({"fixes": "all"}) == sorted(fixes.CLEAR)


@needs_data
def test_totals_fixes():
    """gettotalfields with every clear fix on the saved run."""
    from test_inputs import saved_totals_inputs

    from artemis_p.totals import gettotalfields

    run, out, ds = saved_totals_inputs("all")
    gettotalfields(out, ds)
    t = out["total"]
    # cCDR is the daily CDR plus the change in soil carbonate, not the whole carbonate stock
    assert np.allclose(t["cCDR"], t["cCDRrh"] + t["cd_soilCarbonate"])
    assert t["cCDR"][-1] < 20 < run["total"]["cCDR"][-1]
    # feedstock Ca+Mg+Na release converted to tCO2/ha once
    assert t["cCaMgNarelease"][-1] == pytest.approx(16.313, rel=1e-4)
    assert run["total"]["cCaMgNarelease"][-1] == pytest.approx(14.357, rel=1e-4)
    # release from native minerals is recorded
    assert t["cCanatrelease"][-1] > 1 and t["cMgnatrelease"][-1] > 1
    # daily passive + active uptake
    assert np.allclose(t["dt_Nuptake"], t["dp_Nuptake"] + t["dk_Nuptake"])


@needs_data
def test_setup_fixes():
    """Setup-time fixes on the rBcomp inputs change only what they should."""
    from test_inputs import _rbcomp_preprocessed

    from artemis_p.preprocess import setup_data_struc
    from artemis_p.reference import EF2_FILES, ef2_resolver, rbcomp_settings

    legacy = _rbcomp_preprocessed()
    settings = rbcomp_settings()
    settings["fixes"] = ("hornblende_formula", "density_units", "eqcec_units", "native_rsa_weighting")
    resolve, ef2 = ef2_resolver()
    ds = setup_data_struc(settings, EF2_FILES, resolve, basedir=ef2)

    pt, lpt = ds["phasetable"], legacy["phasetable"]
    assert lpt.loc["HornblendeRLBT", "formula"] == "Ca2Mg52Si8O22(OH)2"
    assert pt.loc["HornblendeRLBT", "formula"] == "Ca2Mg5Si8O22(OH)2"
    assert ds["stoich"]["HornblendeRLBT"]["Mg"] == 5
    assert np.allclose(pt["rho_kgm3"], lpt["rho_kgm3"] * 1e6, equal_nan=True)
    assert np.allclose(ds["cell"]["eqCEC"], legacy["cell"]["eqCEC"] / 100)
    # native RSA scales with each mineral's proportion, so quartz keeps far more area than magnetite
    for ph in ("QuartzRLBT", "MagnetiteRLBT"):
        assert np.all(ds["cell"]["RSA"][ph] <= legacy["cell"]["RSA"][ph] + 1e-30)
    q = ds["cell"]["RSA"]["QuartzRLBT"][0] / ds["cell"]["RSA"]["MagnetiteRLBT"][0]
    assert q > 100 * legacy["cell"]["RSA"]["QuartzRLBT"][0] / legacy["cell"]["RSA"]["MagnetiteRLBT"][0]


@needs_data
def test_crop_case_fix():
    """Crop-specific Michaelis-Menten rows are found when names differ in case."""
    from test_inputs import _rbcomp_setup

    from artemis_p.preprocess import setup_data_struc
    from artemis_p.reference import EF2_FILES, ef2_resolver, rbcomp_settings
    from artemis_p.setuprun import complete_setup

    legacy = _rbcomp_setup()
    settings = rbcomp_settings()
    settings["fixes"] = ("crop_case",)
    resolve, ef2 = ef2_resolver()
    ds = complete_setup(setup_data_struc(settings, EF2_FILES, resolve, basedir=ef2), legacy_airpor=True)
    planted = np.asarray(legacy["run"]["plant"]) > 0
    assert np.all(legacy["run"]["Nindex"][planted] == 1)  # the default row only
    assert np.any(ds["run"]["Nindex"][planted] > 1)
