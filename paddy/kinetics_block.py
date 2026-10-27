"""The KINETICS block of a paddy cell: formulas of paddy.dat's reactions and pools."""

N2C = 0.1  # mol NH3 per mol CH2O (SCEPTER default; the spec's production value is 0.03)

FORMULA = {  # net change to the solution per mol of extent
    "g2o2": "C 1 O 3 H 2 Oxg -1",
    "g2fe3": "C 1 Fe 4 O 13 H 14",
    "g2gr": "C 1 Fe 6 O 17 H 18",
    "g2gt": "C 1 Fe 4 O 9 H 6",
    "g2so4": "C 1 O 1 H 1.5 S -0.5 Hsg 0.5",
    "g2ch4": "C 0.5 O 1 Mtg 0.5",
    "ch4o2": "Mtg -1 Oxg -2 C 1 O 4 H 4",
    "fe2o2": "Fe -1 Oxg -0.25 H -3 O -2.5",
    "fe2gr": "Fe -1 H -2 O -2",
    "groxy": "Oxg -0.25 H -1 O -0.5",
    "fe2fs": "Fe -1 Hsg -1 H 1",
    "fesox": "Oxg -2 Fe 1 S 1 O 4",
    "amo2o": "Amm -1 Oxg -2 N 1 O 4 H 3",
    "g2n21": "C 1 None 1 O 1 H 2 N -1",
    "g2n0": "C 1 Nzero 0.8 O 1 H 2 N -0.8",
    "g2n2p": "C 1 O 1 H 2 None -4 Nzero 4",
    "ureah": "C 1 O 2 Amm 2",
}
OM_REACTIONS = {"g2o2", "g2fe3", "g2gr", "g2gt", "g2so4", "g2ch4", "g2n21", "g2n0", "g2n2p"}
POOLS = ("Om", "Ferrihydrite", "Green_rust", "Goethite_k", "FeS_k", "Urea_k")


def kinetics(n, pools, f=500.0, sw=1.0, omin=0.0, uin=0.0, tc=25.0, n2c=N2C, steps="86400 in 1", extra=""):
    """KINETICS n; pools in mol m-3 bulk, converted to mol per kgw with f."""
    L = [f"KINETICS {n}", "Paddy_rates", "    -formula H 0", "    -m 1", f"    -parms {f} {sw} {omin} {uin} {tc}"]
    for r, form in FORMULA.items():
        if r in OM_REACTIONS and n2c:
            form += f" Amm {n2c}"
        L += [r, f"    -formula {form}", "    -m 1e10"]
    for p in POOLS:
        L += [p, "    -formula H 0", f"    -m {pools.get(p, 0) / f:.12g}"]
    L += [f"    -steps {steps}", "    -cvode true" if "cvode" in extra else "", "    -bad_step_max 1000"]
    return "\n".join(l for l in L if l) + "\n"
