"""The KINETICS block of a paddy cell: formulas of paddy.dat's reactions and pools."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from artemis_p.paddy import FORMULA, N2C, OM_REACTIONS, POOLS  # noqa: E402


def kinetics(n, pools, f=500.0, sw=1.0, omin=0.0, uin=0.0, tc=25.0, n2c=N2C, steps="86400 in 1", extra=""):
    """KINETICS n; pools in mol m-3 bulk, converted to mol per kgw with f."""
    L = [f"KINETICS {n}", "Paddy_rates", "    -formula H 0", "    -m 1", f"    -parms {f} {sw} {omin} {uin} {tc} 0 36.5 0.01 0"]  # no gas exchange
    for r, form in FORMULA.items():
        if r in OM_REACTIONS and n2c:
            form += f" Amm {n2c}"
        L += [r, f"    -formula {form}", "    -m 1e10"]
    for p in POOLS:
        L += [p, "    -formula H 0", f"    -m {pools.get(p, 0) / f:.12g}"]
    L += [f"    -steps {steps}", "    -cvode true" if "cvode" in extra else "", "    -bad_step_max 1000"]
    return "\n".join(l for l in L if l) + "\n"
