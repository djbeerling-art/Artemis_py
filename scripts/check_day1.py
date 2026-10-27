"""Check the port's day-1 PHREEQC input for each published EF2 scenario.

Usage: python scripts/check_day1.py [NAME ...]

For each run script in ``artemis_p.ef2.SCRIPTS`` this builds the settings as the
script does (EFrunsettings2, setupdatastruc, getphreeqcrun arguments), writes
the day-1 input and compares it with the published ``.in_init``. It prints one
line per scenario: OK, the first differing lines, or the step that is not
ported. Needs ARTEMIS, ARTEMIS_PHREEQC and ARTEMIS_FLANAGAN as for the tests.
"""

import copy
import datetime
import difflib
import os
import re
import sys
import traceback
from pathlib import Path

from artemis_p import ef2
from artemis_p.phreeqc_input import PhreeqcRun
from artemis_p.preprocess import setup_data_struc
from artemis_p.reference import EF2_FILES, ef2_resolver, rbcomp_run_settings
from artemis_p.setuprun import complete_setup

# getphreeqcrun flags that take no value, as PhreeqcRun options
FLAGS = {
    "nodust": {"usedust": 0}, "usedust": {"usedust": 1},
    "nosorption": {"usesurface": 0, "useexchange": 0}, "noexchange": {"useexchange": 0},
    "nosurface": {"usesurface": 0}, "nosurfaces": {"usesurface": 0},
    "nouptake": {"includeuptake": 0}, "noplant": {"includeuptake": 0},
    "nodecomp": {"includehumus": 0, "includeresidues": 0, "includemanure": 0, "includerootresidues": 0},
    "cvode": {"usecvode": 1},
    "includedynexchange": {"includedynexchange": 1}, "dynexchange": {"includedynexchange": 1},
    "nodynamicsurfaces": {"includedynamicsurfaces": 0}, "includedynamicsurfaces": {"includedynamicsurfaces": 1},
    "noNkinetics": {"includeNkinetics": 0}, "notreatments": {"includetreatments": 0},
}
# getphreeqcrun arguments that take a value and map straight to an option
VALUED = {"runname", "runlab", "nyears", "ndays", "fixCO2gmol", "fixO2gmol", "noappend", "restart", "restartfile"}


def run_options(ds, runarg):
    """PhreeqcRun options for getphreeqcrun's arguments; 'exclude' filters ds['treatments']."""
    opts, i = {}, 0
    while i < len(runarg):
        a = runarg[i]
        if a in FLAGS:
            opts.update(FLAGS[a])
            i += 1
        elif a in VALUED:
            opts[a] = runarg[i + 1]
            i += 2
        elif a == "exclude":
            ex = runarg[i + 1]
            ex = [ex] if isinstance(ex, str) else ex
            if "treatments" in ds:
                ds["treatments"] = [t for t in ds["treatments"] if not any(e in t for e in ex)]
            i += 2
        else:
            raise NotImplementedError(f"getphreeqcrun argument {a!r}")
    # nyears is for the full run; day 1 only is written here
    opts.pop("nyears", None)
    opts.pop("ndays", None)
    return opts


def published(name, outdir):
    """The published .in_init files of a run (some runs were repeated), with their start times."""
    for f in sorted(Path(outdir).glob(f"{name}_*.in_init")):
        text = f.read_text()
        m = re.search(r"starting: (\S+ \S+)", text.splitlines()[0])
        yield f.name, text, datetime.datetime.strptime(m.group(1), "%d-%b-%Y %H:%M:%S")


def build(name, published_parms=True, fixes=None):
    """(ds, getphreeqcrun arguments, published runname) for a scenario, set up as its run
    script does, with the authors' paths; raises NotImplementedError for missing options.
    ``fixes`` is ``ds['fixes']`` (artemis_p.fixes)."""
    settings, runarg = ef2.scenario(name)
    if fixes:
        settings["fixes"] = tuple(fixes)
    if published_parms:
        settings["parm"].update(ef2.PUBLISHED_PARMS.get(name, {}))
    pub_name = (runarg[runarg.index("runname") + 1] if "runname" in runarg
                else ef2.POST_RUNNAME.get(name, settings.get("runname", name)))
    resolve, ef2dir = ef2_resolver()
    files = ef2.USEFILES.get(name, EF2_FILES)
    ds = complete_setup(setup_data_struc(settings, files, resolve, basedir=ef2dir), legacy_airpor=True)
    rs = rbcomp_run_settings()
    rs["runname"] = settings.get("runname", rs["runname"])
    rs["runlab"] = settings.get("runlab", rs["runlab"])
    rs.update(ef2.POST_DS.get(name, {}))
    ds.update(rs)
    return ds, runarg, pub_name


def check(name, outdir, published_parms=True):
    """(result, published runname) for one scenario; result starts with OK when the day-1
    input is byte-identical with one of the run's published inputs."""
    try:
        ds0, runarg, pub_name = build(name, published_parms)
    except NotImplementedError as e:
        return f"SETUP not ported: {e}", name
    cands = list(published(pub_name, outdir))
    if not cands:
        return "NO PUBLISHED FILE", pub_name
    best = None
    for fname, text, now in cands:
        ds = copy.deepcopy(ds0)
        try:
            opts = run_options(ds, runarg)
            ours = PhreeqcRun(ds, now=now, **opts).write_input(1)
        except NotImplementedError as e:
            return f"RUN not ported: {e}", pub_name
        if ours == text:
            return f"OK ({fname})", pub_name
        diff = list(difflib.unified_diff(text.splitlines(), ours.splitlines(), "published", "port", n=0, lineterm=""))
        diff = [l for l in diff[2:] if l[:1] in "+-"]
        if best is None or len(diff) < len(best[1]):
            best = (fname, diff)
    fname, diff = best
    return f"DIFF from {fname} ({len(diff)} lines)\n    " + "\n    ".join(l[:160] for l in diff[:12]), pub_name


def main(names):
    outdir = os.environ["ARTEMIS_PHREEQC"]
    for name in names or ef2.SCRIPTS:
        try:
            res, pub = check(name, outdir)
        except Exception as e:  # report and carry on with the next scenario
            res, pub = f"ERROR {type(e).__name__}: {e}\n    " + traceback.format_exc().splitlines()[-3].strip(), name
        print(f"{name:10s} ({pub}): {res}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
