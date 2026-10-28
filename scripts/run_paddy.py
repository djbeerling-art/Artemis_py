"""Run an EF2 scenario as a rice paddy (paddy step 2: the paddy chemistry in each cell,
with a flood and drain schedule).

Usage: python scripts/run_paddy.py NAME WORKDIR NDAYS [OPTIONS]

NAME is a key of ``artemis_p.ef2.SCRIPTS`` (rCcomp is the control). The cells are
flooded for ``flooddays`` days from day of year ``floodday`` every year (defaults 152
and 113, the SCEPTER reference season) and drained otherwise. The soil starts with the
SCEPTER reference site's ferrihydrite and goethite; OM and urea enter as in the batch
cell (paddy/batch.py). OPTIONS are as for run_scenario.py, plus floodday=, flooddays=.
"""

import os
import pickle
import sys
import time
from pathlib import Path

import numpy as np
from check_day1 import run_options

from artemis_p import ef2
from artemis_p.phreeqc_input import PhreeqcRun
from artemis_p.preprocess import setup_data_struc
from artemis_p.reference import EF2_FILES, ef2_resolver, rbcomp_run_settings, stage_ef2_dir
from artemis_p.setuprun import complete_setup

REFERENCE_POOLS = {"Ferrihydrite": 446.17382521816779, "Goethite_k": 297.449}  # mol m-3, SCEPTER --refsoil


def paddy_settings(floodday=152, flooddays=113):
    def flooded(ds):
        doy = np.asarray(ds["run"]["doy"], float)
        return (doy >= floodday) & (doy < floodday + flooddays)

    def urea(z):  # mol m-3 yr-1, as paddy/batch.py
        return 450 / 60.06 * 0.18441163442856148 * np.exp(-z / 0.3) / 0.3

    return {"flooded": flooded, "pools": dict(REFERENCE_POOLS), "uin": urea}


def build(name, floodday, flooddays):
    settings, runarg = ef2.scenario(name)
    settings["parm"].update(ef2.PUBLISHED_PARMS.get(name, {}))
    settings["paddy"] = paddy_settings(floodday, flooddays)
    resolve, ef2dir = ef2_resolver()
    ds = setup_data_struc(settings, ef2.USEFILES.get(name, EF2_FILES), resolve, basedir=ef2dir)
    ds = complete_setup(ds, legacy_airpor=True)
    rs = rbcomp_run_settings()
    rs["runname"] = name + "paddy"
    rs["runlab"] = settings.get("runlab", rs["runlab"]) + " paddy"
    rs.update(ef2.POST_DS.get(name, {}))
    ds.update(rs)
    return ds, runarg


def main(name, work, ndays, options=None):
    work = Path(work)
    artemis = Path(os.environ["ARTEMIS"])
    ef2dir = stage_ef2_dir(work, artemis, os.environ.get("ARTEMIS_FLANAGAN"))
    opts = {}
    for kv in (options or "").split(","):
        if kv:
            k, v = kv.split("=")
            opts[k] = int(v) if v.isdigit() else v
    ds, runarg = build(name, opts.pop("floodday", 152), opts.pop("flooddays", 113))
    dbdir = f"{artemis / 'phreeqcdbincludefiles'}/"
    ds.update(dbdir=dbdir, plantratesfile=dbdir + "plantRATES5", schemeratesfile=dbdir + "SWATRATESPHASEStunable",
              outdir=f"{ef2dir}/")
    if "INCLUDE" in ds:
        ds["INCLUDE"] = [f.replace(f.split("phreeqcdbincludefiles")[0], f"{artemis}/") if "phreeqcdbincludefiles" in f
                         else f for f in ds["INCLUDE"]]
    opts = {**run_options(ds, runarg), **opts}
    r = PhreeqcRun(ds, **opts)
    t0 = time.time()
    for i in range(1, int(ndays) + 1):
        res = r.run_day(i)
        if not res.ok:
            if r.engine is not None:  # the state when the day failed, for a look with the executable
                r.engine.write_dump(r.fname + ".dmp")
            sys.exit(f"PHREEQC failed on day {i}: {res.errors[:3]}")
        if i % 20 == 0:
            print(f"day {i}: {time.time() - t0:.0f} s", flush=True)
            pickle.dump(r.out, open(work / "out.pkl", "wb"))  # so a long run can be looked at as it goes
            pickle.dump(r.ds, open(work / "ds.pkl", "wb"))
    print(f"{ndays} days in {time.time() - t0:.0f} s")
    pickle.dump(r.out, open(work / "out.pkl", "wb"))
    pickle.dump(r.ds, open(work / "ds.pkl", "wb"))


if __name__ == "__main__":
    main(*sys.argv[1:5])
