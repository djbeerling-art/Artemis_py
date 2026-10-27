"""MATLAB sprintf behaviour that shows up in the PHREEQC input files."""

import numpy as np

from artemis_p.mformat import sprintf


def test_integer_d_with_fraction_prints_e():
    assert sprintf("%10s %d # kg", "water", 31.6) == "     water 3.160000e+01 # kg"
    assert sprintf("%d", 12.0) == "12"


def test_nan_and_recycling():
    assert sprintf("%g;", np.array([1.5, np.nan])) == "1.5;NaN;"
    assert sprintf("%s Nhuma %g Phum %g", "x", [0.1, 0.2]) == "x Nhuma 0.1 Phum 0.2"
    assert sprintf("%13.4g in %d steps", 86400, 3) == "     8.64e+04 in 3 steps"
