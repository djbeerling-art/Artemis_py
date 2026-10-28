"""Run a published EF2 scenario with the port.

Usage: python scripts/run_scenario.py NAME WORKDIR [NDAYS [OPTIONS]]

NAME is a key of ``artemis_p.ef2.SCRIPTS`` (e.g. rCcomp). NDAYS defaults to
the whole run. OPTIONS are PhreeqcRun options as for run_rbcomp.py, e.g.
``engine=iphreeqc,cellworkers=4,printoutput=0``; ``fixes=all+native_rsa_weighting``
sets ``ds['fixes']`` (names joined by +); ``noappend=PATH`` restarts from the port's own
spin-up dump (PATH.dmp, or PATH.mat beside it) instead of the authors'. ``gasdiff=CO2+O2``
turns on soil gas diffusion for those gases (artemis_p.gasdiff); ``gastop=``, ``gassource=``,
``gasmodel=``, ``gasfeedback=`` and ``gastau=`` set its other settings. The output structure is
pickled to WORKDIR/out.pkl and the data structure to WORKDIR/ds.pkl.
"""

import os
import pickle
import sys
import time
from pathlib import Path

from check_day1 import build, run_options

from artemis_p.phreeqc_input import PhreeqcRun
from artemis_p.reference import stage_ef2_dir


def main(name, work, ndays=None, options=None):
    work = Path(work)
    artemis = Path(os.environ["ARTEMIS"])
    ef2 = stage_ef2_dir(work, artemis, os.environ.get("ARTEMIS_FLANAGAN"))
    opts = {}
    for kv in (options or "").split(","):
        if kv:
            k, v = kv.split("=")
            opts[k] = int(v) if v.isdigit() else v
    fixes = opts.pop("fixes", "")
    ds, runarg, _ = build(name, fixes=[f for f in fixes.split("+") if f])
    if "gasdiff" in opts:
        ds["gasdiff"] = {"gases": tuple(str(opts.pop("gasdiff")).split("+"))}
        for k in ("top", "source", "model", "feedback", "tau"):
            if "gas" + k in opts:
                ds["gasdiff"][k] = opts.pop("gas" + k)
    dbdir = f"{artemis / 'phreeqcdbincludefiles'}/"
    ds.update(dbdir=dbdir, plantratesfile=dbdir + "plantRATES5", schemeratesfile=dbdir + "SWATRATESPHASEStunable",
              outdir=f"{ef2}/")
    if "INCLUDE" in ds:
        ds["INCLUDE"] = [f.replace(f.split("phreeqcdbincludefiles")[0], f"{artemis}/") if "phreeqcdbincludefiles" in f
                         else f for f in ds["INCLUDE"]]
    opts = {**run_options(ds, runarg), **opts}
    r = PhreeqcRun(ds, **opts)
    ndays = int(ndays) if ndays else r.ilast - r.istart + 1
    t0 = time.time()
    for i in range(1, ndays + 1):
        res = r.run_day(i)
        if not res.ok:
            sys.exit(f"PHREEQC failed on day {i}: {res.errors[:3]}")
        if i % 100 == 0:
            print(f"day {i}: {time.time() - t0:.0f} s", flush=True)
    print(f"{ndays} days in {time.time() - t0:.0f} s")
    if ndays == r.ilast - r.istart + 1:
        r.finish()
    elif r.engine is not None:  # a partial run still leaves the dump a restart reads
        r.engine.write_dump(r.fname + ".dmp")
    print(f"dump: {r.fname}.dmp")
    pickle.dump(r.out, open(work / "out.pkl", "wb"))
    pickle.dump(r.ds, open(work / "ds.pkl", "wb"))


if __name__ == "__main__":
    main(*sys.argv[1:5])
