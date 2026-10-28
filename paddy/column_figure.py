"""Paddy step 2 figure: an ARTEMIS column through a flood season, three depths.

Usage: python paddy/column_figure.py WORKDIR OUT.png [SCEPTER_RUNDIR]

WORKDIR holds out.pkl and ds.pkl from scripts/run_paddy.py. With SCEPTER_RUNDIR, the
SCEPTER reference run's top cell (a closed 113-day flood) is drawn in grey from the
day the flood starts, as a guide only: that run has no percolation or floodwater exchange.
"""
import os
import pickle
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

DEPTH = ("#9ec5f0", "#4a90dc", "#174f8f")  # one hue, light to dark with depth
INK, INK2, GRID, REF = "#0b0b0b", "#52514e", "#e4e3df", "#a3a29d"


def main(work, png, scepter=None):
    out = pickle.load(open(os.path.join(work, "out.pkl"), "rb"))
    ds = pickle.load(open(os.path.join(work, "ds.pkl"), "rb"))
    soil, cell = out["soil"], ds["cell"]
    n = int(np.max(np.flatnonzero(soil["pH"][0] != 0))) + 1
    t = np.arange(1, n + 1)
    fl = ds["paddy"]["flooded"][:n]
    vsoil = np.asarray(cell["vsoil"], float)[:, None]
    M = {q: v[:, :n] / vsoil * 1000 for q, v in out["run"]["paddy"].items()}  # mol m-3
    cells = [0, 1, 2]
    panels = [
        ("Dissolved O$_2$ (µM)", soil["Oxg"] * 1e6),
        ("Nitrate (µM)", soil["N"] * 1e6),
        ("Sulfate (µM)", soil["S"] * 1e6),
        ("CH$_4$ (atm equivalent)", soil["Mtg"] / 10 ** -2.86),
        ("Dissolved Fe(II) (mM)", soil["Fe2"] * 1e3),
        ("Green rust (mol m$^{-3}$)", M["Green_rust"]),
        ("FeS (mol m$^{-3}$)", M["FeS_k"]),
        ("pH", soil["pH"]),
    ]
    ref = None
    if scepter:
        from read_scepter import read
        s = read(scepter)
        t0 = np.flatnonzero(fl)[0] + 1
        ts = t0 + s["time"] * 365
        ref = [s["po2"][:, 0] * 1.26e-3 * 1e6, s["no3"][:, 0] * 1e6, s["so4"][:, 0] * 1e6, s["pch4"][:, 0],
               s["fe2(tot)"][:, 0] * 1e3, s["gr"][:, 0], s["fes"][:, 0], s["ph"][:, 0]]
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2, "font.family": "DejaVu Sans"})
    fig, axes = plt.subplots(2, 4, figsize=(12, 5.8), sharex=True)
    on = np.flatnonzero(fl)
    for j, (ax, (title, v)) in enumerate(zip(axes.flat, panels)):
        if len(on):
            ax.axvspan(on[0] + 1, on[-1] + 1, color="#f1f0ec", lw=0)
        if ref is not None:
            ax.plot(ts, ref[j], color=REF, lw=1.5, label="SCEPTER top cell (closed flood)")
        for c, col in zip(cells, DEPTH):
            ax.plot(t, v[c, :n], color=col, lw=2, label=f"{cell['middepth'][c]:.0f} cm")
        ax.set_title(title, fontsize=9.5, color=INK, loc="left")
        ax.grid(color=GRID, lw=0.6)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlim(max(1, on[0] - 9) if len(on) else 1, n)  # from just before the flood
    for ax in axes[1]:
        ax.set_xlabel("Run day (shaded: flooded)")
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=len(lab), frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(png, dpi=150)


if __name__ == "__main__":
    main(*sys.argv[1:4])
