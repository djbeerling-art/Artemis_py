"""PHREEQC batch cell for comparison with a SCEPTER paddy run (cell 1 of a 15-cell column, z = 0.01 m).

Usage: python paddy/batch.py OUT [noammx] [scx]

Writes OUT.in, runs PHREEQC ($PHREEQC) and leaves OUT.sel. ``scx`` uses SCEPTER's
exchange selectivities (no NH4+ exchange), ``noammx`` the database's exchanger
without NH4+ exchange; the default is the database's exchanger as it is.
The SCEPTER reference was run with paddy/run_flooded_column.py (branch
paddy/fe3-reduction):
  batch3 fe2o2,g2fe3,g2gr,g2so4,g2ch4,ch4o2,fe2gr,fe2fs,fesox,amo2o,g2n21,g2n2p,g2gt,g2o2,groxy,sidox 0.31
  --ch4 --n --kaol=0.05 --refsoil --gt=0.02 --nz=15 --q=1e-6 --kfw=0 --nfert=120 --om=450 --basalt=0
  --o2inh=2.5e-7 --n2oy=6.2:1:0.1 --k=paddy/krxn_calibrated_v15.in
"""
import subprocess, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from kinetics_block import kinetics

DB = os.path.join(os.environ.get("ARTEMIS", ""), "phreeqcdbincludefiles", "THAmmOrgAlP.dat")
HERE = os.path.dirname(os.path.abspath(__file__))
F = 500.0               # L pore water per m3 bulk (porosity 0.5, flooded)
Z = 0.01                # cell depth, m
OMIN = 450 / 12 * np.exp(-Z / 0.3) / 0.3          # mol C m-3 yr-1 (SCEPTER zsupp_plant 0.3 m)
UIN = 450 / 60.06 * 0.18441163442856148 * np.exp(-Z / 0.3) / 0.3  # mol urea m-3 yr-1
CEC = 188.0             # eq m-3 bulk, SCEPTER's initial exchanger total in this cell
POOLS = {"Ferrihydrite": 446.17382521816779, "Goethite_k": 297.449}
DAYS = 113


# SCEPTER's exchange selectivities (cec.in: log K vs H for Na 5.9, K 4.8, Ca 10.47, Mg 10.786, Al 16.47;
# Fe2+ as Mg + 0.158), rewritten relative to NaX for PHREEQC; no NH4+ exchange, as in SCEPTER
SCEPTER_X = """EXCHANGE_SPECIES
    H+ + X- = HX; log_k 5.9
    K+ + X- = KX; log_k 1.1
    Ca+2 + 2X- = CaX2; log_k 1.33
    Mg+2 + 2X- = MgX2; log_k 1.014
    Fe+2 + 2X- = FeX2; log_k 0.856
    Al+3 + 3X- = AlX3; log_k 1.23
    AmmH+ + X- = AmmHX; log_k -20
"""


def write(path, days=DAYS, ndt=None, ammx=True, scx=False):
    ndt = ndt or days
    noammx = SCEPTER_X if scx else "" if ammx else "EXCHANGE_SPECIES\n    AmmH+ + X- = AmmHX\n    log_k -20\n"
    txt = f"""DATABASE {DB}
INCLUDE$ {HERE}/paddy.dat
{noammx}SOLUTION 1 irrigation water
    temp 25
    pH 7.6 charge
    units mol/kgw
    Ca 5e-4; Mg 2e-4; Na 3e-4; K 5e-5; Si 2e-4; S(6) 2e-4; N 1e-6; Cl 1.15e-3
    Fe 1e-10; Al 1e-15; Amm 1e-12; Mtg 2.7e-9; Hsg 1e-15; None 1.58e-8; Nzero 1e-12
EQUILIBRIUM_PHASES 1
    CO2(g) {np.log10(3.16e-4):.4f} 10
    Oxg(g) {np.log10(0.21):.4f} 10
SAVE solution 1
END
USE solution 1
EXCHANGE 1
    X {CEC / F:.6g}
    -equilibrate with solution 1
SAVE exchange 1
END
INCREMENTAL_REACTIONS true
USE solution 1
USE exchange 1
{kinetics(1, POOLS, f=F, omin=OMIN, uin=UIN, steps=f"{days * 86400} in {ndt}")}
SELECTED_OUTPUT 1
    -file {path}.sel
    -reset false
    -time true
    -pH true
    -pe true
    -totals Oxg N S(6) Hsg Mtg C(4) Fe(2) Fe(3) Amm None Ca Mg
    -molalities Fe+2 FeX2 CaX2 MgX2 HX AmmH+ Hsg-
    -kinetic_reactants Om Ferrihydrite Green_rust Goethite_k FeS_k Urea_k g2o2 g2fe3 g2gr g2gt g2so4 g2ch4 ch4o2 fe2o2 fe2gr groxy fe2fs fesox amo2o g2n21 g2n0 g2n2p ureah
    -saturation_indices Siderite Calcite Mackinawite
END
"""
    open(path + ".in", "w").write(txt)


def run(path, **kw):
    write(path, **kw)
    exe = os.environ.get("PHREEQC", "phreeqc")
    r = subprocess.run([exe, path + ".in", path + ".out", DB], capture_output=True, text=True)
    err = [l for l in open(path + ".out", errors="replace") if "ERROR" in l]
    if err:
        print("".join(err[:10]))
    return read_sel(path + ".sel")


def read_sel(f):
    L = [l.split() for l in open(f)]
    h = L[0]
    a = np.array([[float(x) for x in r] for r in L[1:]])
    return {k: a[:, j] for j, k in enumerate(h)}


if __name__ == "__main__":
    import time
    t = time.time()
    d = run(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "batch"), ammx="noammx" not in sys.argv[2:], scx="scx" in sys.argv[2:])
    print(f"{time.time() - t:.1f} s, {len(d['time'])} rows")
