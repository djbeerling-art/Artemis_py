"""Run the rBcomp reference case and compare it with the authors' saved MATLAB run.

Usage: python scripts/run_rbcomp.py NDAYS WORKDIR

Needs ARTEMIS, ARTEMIS_PHREEQC (and ARTEMIS_FLANAGAN if used) as for the tests.
The day timestamps come from the saved run so that the last day's ``.in`` can be
compared byte for byte with the published one. A full run (2,192 days) takes
about 45 minutes.
"""

import copy
import datetime
import difflib
import os
import pickle
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))

from test_inputs import RBCOMP, _rbcomp_setup  # noqa: E402

from artemis_p import matref, phreeqc  # noqa: E402
from artemis_p.phreeqc_input import PhreeqcRun  # noqa: E402
from artemis_p.reference import AUTHOR_DBDIR, AUTHOR_OUTDIR, rbcomp_run_settings, stage_ef2_dir  # noqa: E402


def compare(grp, ours, ref, ndays, tol=1e-3):
    for f in sorted(ours):
        if f not in ref:
            continue
        b = np.asarray(ref[f], float)
        a = np.asarray(ours[f], float).reshape(b.shape)[..., :ndays]
        b = b[..., :ndays]
        scale = np.nanmax(np.abs(b)) if b.size else 0.0
        rel = np.abs(a - b) / (np.abs(b) + 1e-6 * scale + 1e-300)
        rel[np.isnan(a) & np.isnan(b)] = 0
        if np.nanmax(rel, initial=0) > tol:
            print(f"{grp:9s} {f:24s} max rel diff {np.nanmax(rel):.3g}, {(rel > tol).sum()}/{rel.size} values")


def main(ndays, work):
    saved = matref.load_run(Path(os.environ["ARTEMIS_PHREEQC"]) / f"{RBCOMP}.mat")
    artemis = Path(os.environ["ARTEMIS"])
    dbdir = f"{artemis / 'phreeqcdbincludefiles'}/"
    ef2 = stage_ef2_dir(work, artemis, os.environ.get("ARTEMIS_FLANAGAN"))
    ds = copy.deepcopy(_rbcomp_setup())
    ds.update(rbcomp_run_settings(dbdir=dbdir, outdir=f"{ef2}/"))

    def stamp(i):
        s = saved["run"]["current"][i - 1].split("starting: ")[1]
        return datetime.datetime.strptime(s, "%d-%b-%Y %H:%M:%S")

    r = PhreeqcRun(ds, now=stamp(1))
    t0 = time.time()
    for i in range(1, ndays + 1):
        res = r.run_day(i, now=stamp(i))
        if not res.ok:
            sys.exit(f"PHREEQC failed on day {i}: {res.errors[:3]}")
        if i % 50 == 0:
            print(f"day {i}: {time.time() - t0:.0f} s", flush=True)
    print(f"{ndays} days in {time.time() - t0:.0f} s")
    pickle.dump(r.out, open(Path(work) / "out.pkl", "wb"))

    for grp in ("soil", "influent", "leachate"):
        compare(grp, r.out[grp], saved[grp], ndays)
    for k in ("RSA", "M0", "grams"):
        for p, v in saved["run"][k].items():
            if p in r.out["run"][k]:
                compare(k, {p: r.out["run"][k][p]}, {p: v}, ndays)
    if ndays == len(saved["run"]["current"]):
        name = f"{RBCOMP}.in"
        pub = (Path(os.environ["ARTEMIS_PHREEQC"]) / name).read_text()
        pub = phreeqc.relocate_input(pub, {AUTHOR_OUTDIR[:-1]: f"{ef2}/", AUTHOR_DBDIR: dbdir})
        diff = list(difflib.unified_diff(pub.splitlines(), (ef2 / name).read_text().splitlines(), lineterm="", n=0))
        print("last-day .in identical to the published one" if not diff else "\n".join(diff[:40]))


if __name__ == "__main__":
    main(int(sys.argv[1]), Path(sys.argv[2]))
