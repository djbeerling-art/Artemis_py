"""The paddy reactions (paddy/paddy.dat, paddy/kinetics_block.py) are balanced: each formula,
added to a porewater, gives its products without moving electrons between redox states."""

import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "paddy"))
from kinetics_block import FORMULA  # noqa: E402

DB = Path(os.environ.get("ARTEMIS", "")) / "phreeqcdbincludefiles" / "THAmmOrgAlP.dat"
EXE = os.environ.get("PHREEQC")
pytestmark = pytest.mark.skipif(not (EXE and DB.exists()), reason="needs $PHREEQC and $ARTEMIS")

# what each reaction should change, per 1e-5 mol of extent (mol/kgw)
EXPECT = {
    "g2o2": {"C(4)": 1, "Oxg": -1}, "g2fe3": {"C(4)": 1, "Fe(2)": 4}, "g2gr": {"C(4)": 1, "Fe(2)": 6},
    "g2gt": {"C(4)": 1, "Fe(2)": 4}, "g2so4": {"C(4)": 1, "S(6)": -0.5, "Hsg": 0.5},
    "g2ch4": {"C(4)": 0.5, "Mtg": 0.5}, "ch4o2": {"Mtg": -1, "Oxg": -2, "C(4)": 1},
    "fe2o2": {"Fe(2)": -1, "Oxg": -0.25}, "fe2gr": {"Fe(2)": -1}, "groxy": {"Oxg": -0.25},
    "fe2fs": {"Fe(2)": -1, "Hsg": -1}, "fesox": {"Oxg": -2, "Fe(2)": 1, "S(6)": 1},
    "amo2o": {"Amm": -1, "Oxg": -2, "N": 1}, "g2n21": {"C(4)": 1, "N": -1, "None": 1},
    "g2n0": {"C(4)": 1, "N": -0.8, "Nzero": 0.8}, "g2n2p": {"C(4)": 1, "None": -4, "Nzero": 4},
    "ureah": {"C(4)": 1, "Amm": 2}, "ch4eb": {"Mtg": -1}, "fwo2": {"Oxg": -1}, "fwco2": {"C(4)": -1},
    "fwch4": {"Mtg": -1},
}
TOTALS = ["C(4)", "C(-4)", "S(6)", "S(-2)", "Fe(2)", "Fe(3)", "O(0)", "Mtg", "Oxg", "Hsg", "N", "Amm", "None", "Nzero"]


def test_formulas_balanced(tmp_path):
    txt = f"DATABASE {DB}\nINCLUDE$ {ROOT / 'paddy' / 'paddy.dat'}\n"
    txt += """SELECTED_OUTPUT
    -file out.sel
    -reset false
    -high_precision true
    -charge_balance true
    -totals """ + " ".join(TOTALS) + """
SOLUTION 1
    pH 7 charge
    Ca 1; Na 1; Cl 1; S(6) 0.5; N 0.2; Fe 0.1; C 2; Amm 0.5; Oxg 0.1; Mtg 0.1; Hsg 0.05; None 0.1; Nzero 0.01
SAVE solution 1
END
"""
    for r, f in FORMULA.items():
        txt += f"USE solution 1\nREACTION 1\n    {f}; 1e-5\nEND\n"
    (tmp_path / "t.in").write_text(txt)
    subprocess.run([EXE, "t.in", "t.out", str(DB)], cwd=tmp_path, check=True, capture_output=True)
    rows = [l.split() for l in (tmp_path / "out.sel").read_text().splitlines()]
    head, rows = rows[0], [list(map(float, r)) for r in rows[1:]]
    first = dict(zip(head, rows[0]))  # the initial solution
    assert len(rows) == 1 + len(FORMULA)
    for r, row in zip(FORMULA, rows[1:]):
        got = dict(zip(head, row))
        assert abs(got["charge"]) < 1e-12, r
        for t in ("C(-4)", "S(-2)", "O(0)"):
            assert got[t] < 1e-15, (r, t)
        for t in TOTALS:
            want = first[t] + EXPECT[r].get(t, 0) * 1e-5
            assert got[t] == pytest.approx(want, rel=1e-6, abs=1e-12), (r, t)
