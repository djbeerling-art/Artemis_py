"""Batch-cell comparison figure: SCEPTER paddy module vs its PHREEQC port.

Usage: python paddy/figure.py SCEPTER_RUNDIR SELFILE OUT.png
"""
import sys
sys.path.insert(0, __import__("os").path.dirname(__file__))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from read_scepter import read
from batch import read_sel, F

S1, S2 = "#2a78d6", "#eb6834"   # categorical slots 1 and 2
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
s = read(sys.argv[1])
p = read_sel(sys.argv[2])
ts, tp = s["time"] * 365, p["time"] / 86400

panels = [
    ("Dissolved O$_2$ (µM)", s["po2"][:, 0] * 1.26e-3 * 1e6, p["Oxg"] * 1e6, False),
    ("Nitrate (µM)", s["no3"][:, 0] * 1e6, p["N"] * 1e6, False),
    ("Sulfate (µM)", s["so4"][:, 0] * 1e6, p["S(6)"] * 1e6, False),
    ("CH$_4$ (atm equivalent)", s["pch4"][:, 0], p["Mtg"] / 1.4e-3, False),
    ("Dissolved Fe(II) (µM)", s["fe2(tot)"][:, 0] * 1e6, p["Fe(2)"] * 1e6, False),
    ("Green rust (mol m$^{-3}$)", s["gr"][:, 0], p["k_Green_rust"] * F, False),
    ("FeS (mol m$^{-3}$)", s["fes"][:, 0], p["k_FeS_k"] * F, False),
    ("pH", s["ph"][:, 0], p["pH"], False),
]
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "font.family": "DejaVu Sans"})
fig, axes = plt.subplots(2, 4, figsize=(12, 5.6), sharex=True)
for ax, (title, a, b, logy) in zip(axes.flat, panels):
    ax.plot(ts, a, color=S1, lw=2, label="SCEPTER (Fortran)")
    ax.plot(tp, b, color=S2, lw=2, ls=(0, (4, 2)), label="ARTEMIS rates in PHREEQC")
    ax.set_title(title, fontsize=9.5, color=INK, loc="left")
    ax.grid(color=GRID, lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(0, 113)
    if logy:
        ax.set_yscale("log")
for ax in axes[1]:
    ax.set_xlabel("Days since flooding")
h, l = axes[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.995))
fig.text(0.5, 0.005, "Flooded cell at 1 cm, reference paddy soil, 450 g C m$^{-2}$ yr$^{-1}$ OM and 120 kg urea-N ha$^{-1}$, "
         "no percolation or gas exchange. SCEPTER: 15-cell column (cell 1). PHREEQC: one cell, SCEPTER's exchange selectivities.",
         ha="center", fontsize=7.5, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 0.95))
out = sys.argv[3]
fig.savefig(out, dpi=160, facecolor="#fcfcfb")
print(out)
