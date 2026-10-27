"""The day-1 PHREEQC input of each published EF2 scenario (scripts/check_day1.py)."""

import os
import sys
from pathlib import Path

import pytest

from artemis_p import ef2
from test_parallel import needs_data

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

NOT_PORTED = {"rBbucket": "bucket hydrology (not needed: David, 27 October 2026)"}


@needs_data
@pytest.mark.parametrize("name", list(ef2.SCRIPTS))
def test_day1_input(name):
    import check_day1

    if name in NOT_PORTED:
        pytest.xfail(NOT_PORTED[name] + " is not ported")
    res, _ = check_day1.check(name, os.environ["ARTEMIS_PHREEQC"])
    assert res.startswith("OK"), res
