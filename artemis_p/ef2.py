"""The EF2 example runs: EFrunsettings2.m and the r*script.m run scripts in examples/EF2.

Each run script sets a few variables (feedstock, exchange method, treatment
option, run category, apatite option, scalings), calls EFrunsettings2 to turn
them into arguments for setupdatastruc (``dsarg``) and getphreeqcrun
(``runarg``), and may add arguments of its own. :data:`SCRIPTS` records what
each script does; :func:`scenario` gives the settings for ``setup_data_struc``
and the options for ``PhreeqcRun``.

The settings start from defaultsettings.m (``reference.rbcomp_settings``
without the run's own parm values) and apply the arguments as setupdatastruc
does: a known field of ds or ds.parm takes the next argument, ``parm=X``
sets ds.parm.X, and ``fsname=``/``runname=``/``runlab=`` set themselves.
"""

import copy

from .reference import AUTHOR_DBDIR, EF2_FILES

BR_SCALINGS = {"RSAscale": 2, "impscale": 5}  # BRscalings.m


def ef_run_settings(*args):
    """EFrunsettings2.m: (dsarg, runarg) from a run script's arguments."""
    args = list(args)
    tuning, soflag, uaflag, upflag = 0, 1, 1, 1
    RSAscale, impscale = 1, 1
    runarg, dsarg = [], ["SAmethod", "traditional"]
    runname, runlab = "r", ""
    fixCO2gmol, fsrep = 10, 0
    fsname = apopt = EXCHmethod = None
    nyears = None
    for i, a in enumerate(args):
        if not isinstance(a, str):
            continue
        if "fsname=" in a:
            fsname = a.replace("fsname=", "")
            runname += fsname[0]
            if fsname[0].lower() == "b":
                runlab = "{\\bf BR}"
            elif fsname[0].lower() == "c":
                runlab = "{\\color{red}\\bf CTL}"
                runarg.append("nodust")
            else:
                runlab = "{\\bf " + fsname + "}"
            dsarg.append(a)
        if a.lower() in ("bucket", "richards"):
            runname += a[0].lower()
            dsarg += ["HYDROmethod", a]
        if "method=" in a:
            dsarg.append(a)
        if a.lower() in ("comp", "eq"):
            runname += a
            EXCHmethod = a
            dsarg += ["EXCHmethod", a]
        if a.lower() == "fixco2gmol":
            fixCO2gmol = args[i + 1]
        if "RSAscale" in a:
            RSAscale = args[i + 1]
        if "impscale" in a:
            impscale = args[i + 1]
        if "fractal" in a:
            dsarg += [a if "parm=" in a else "parm=" + a, args[i + 1]]
            v = args[i + 1]
            if "surface" in a:
                runname = (runname + "fDs" + _g(v)).replace(".", "p")
                runlab = f"{runlab} fDs={_g(v)}"
            elif "roughness" in a:
                runname = (runname + "fDr" + _g(v)).replace(".", "p")
                runlab = f"{runlab} fDr={_g(v)}"
        if "nyears" in a:
            runarg += [a, args[i + 1]]
            nyears = args[i + 1]
        if "ndays" in a:
            runarg += [a, args[i + 1]]
        if a.startswith("noso"):
            runarg.append("nosorption")
            if not tuning:
                soflag = 0
                runname = runname.replace("comp", "noso")
        if a.startswith("noup"):
            runarg.append("nouptake")
            if not tuning:
                upflag = 0
                runname = runname.replace("comp", "noup")
        if a.startswith("noua"):
            runarg += ["exclude", "UAN28"]
            if not tuning:
                uaflag = 0
                runname = runname.replace("comp", "noua")
        if a.lower() == "tuning" and fsname and fsname[0].lower() == "b":
            tuning = 1
        if a.lower() == "fsrep":
            fsrep = args[i + 1]
        if a.lower() in ("noap", "kinetic", "rkap", "cvap", "implicit"):
            apopt = a
        if "parm=" in a:
            if "fractal" not in a:  # MATLAB adds fractal parms twice; the second is the same value
                dsarg += [a, args[i + 1]]
            else:
                dsarg += [a, args[i + 1]]
    implicit = ["Titanite_blueridge", "Hydroxyapatite_blueridge"]
    isbr = bool(fsname) and fsname[0].lower() == "b"
    if tuning and isbr:
        ap = (apopt or "").lower()
        if ap == "noap":
            runlab = "Ignore Ap+Tit"
            runname = runname.replace("comp", "noap")
            dsarg += ["parm=ignore", list(implicit)]
        elif "rk" in ap:
            runlab = "Kinetic Ap (RK, ignore Tit)"
            runname = runname.replace("comp", "rkap")
            dsarg += ["parm=setkin", ["Hydroxyapatite_blueridge"], "parm=ignore", ["Titanite_blueridge"]]
        elif "cv" in ap:
            runarg.append("cvode")
            runlab = "Kinetic Ap (cvode, ignore Tit)"
            runname = runname.replace("comp", "cvap")
            dsarg += ["parm=setkin", ["Hydroxyapatite_blueridge"], "parm=ignore", ["Titanite_blueridge"]]
        elif "kin" in ap:
            runlab = "Kinetic Ap (ignore Tit)"
            runname = runname.replace("comp", "kap")
            dsarg += ["parm=setkin", ["Hydroxyapatite_blueridge"], "parm=ignore", ["Titanite_blueridge"]]
        elif impscale > 0:
            runname = runname.replace("comp", "") + "imp" + _g(impscale)
            runlab = f"Ap+Tit imp (imp X {_g(impscale)}, RSA X {_g(RSAscale)})"
            dsarg += ["parm=implicit", list(implicit)]
        runname = (runname + "s" + _g(RSAscale)).replace(".", "p")
    elif nyears is not None and nyears > 10:
        if impscale > 0 and isbr:
            dsarg += ["parm=implicit", list(implicit)]
        runname = f"rl{fsname[0]}rep{fsrep:d}"
        runlab = f"{runlab} ({fsrep:d} treatments)"
    else:
        if impscale > 0 and isbr:
            dsarg += ["parm=implicit", list(implicit)]
        runlab += "+UAN28" if uaflag else "{\\color{red}-UAN28}"
        runlab += "+uptake" if upflag else "{\\color{red}-uptake}"
        runlab += "+sorption" if soflag else "{\\color{red}-sorption}"
    if fixCO2gmol == 0:
        runname = runname.replace(EXCHmethod, "closed")
        runlab += "{\\color{red} closed system}"
        dsarg += ["parm=fixCO2gmol", fixCO2gmol]
    dsarg += ["runname=" + runname, "runlab=" + runlab, "parm=RSAscale", RSAscale, "parm=impscale", impscale]
    if nyears is not None and not any(isinstance(r, str) and r.lower() in ("nyears", "ndays") for r in runarg):
        runarg += ["nyears", nyears]
    return dsarg, runarg


def _g(v):
    """MATLAB %g."""
    return "%g" % v


def _script(fsname="Blueridge", EXCHmethod="comp", tropt="", runcat="main", apopt="implicit",
            RSAscale=None, impscale=None, extra=(), with_impscale=True, before=()):
    """The standard EFrunsettings2 call of a run script; RSAscale/impscale default to BRscalings."""
    RSAscale = BR_SCALINGS["RSAscale"] if RSAscale is None else RSAscale
    impscale = BR_SCALINGS["impscale"] if impscale is None else impscale
    args = ["fsname=" + fsname, *before, EXCHmethod, runcat, apopt, tropt, "RSAscale", RSAscale]
    if with_impscale:
        args += ["impscale", impscale]
    return args + list(extra)


# runname of the published input: (EFrunsettings2 arguments, extra setupdatastruc arguments,
# extra getphreeqcrun arguments), from the r*script.m files
SCRIPTS = {
    "rBcomp": (_script(), [], []),
    "rBeq": (_script(EXCHmethod="eq"), [], []),
    "rBnoso": (_script(tropt="noso"), [], []),
    "rBnoua": (_script(tropt="noua"), [], []),
    "rBnoup": (_script(tropt="noup"), [], []),
    "rBclosed": (_script(extra=["fixCO2gmol", 0]), [], []),
    "rBatm": (_script(), ["runname", "rBatm", "parm=soilco2opt", "atmospheric"], []),
    "rBbucket": (_script(extra=["bucket"]), [], ["runname", "rBbucket"]),
    "rBnodec": (_script(), ["runname", "rBnodec", "runlab", "No organic matter decomposition"], ["nodecomp"]),
    "rBnonit": (_script(extra=["parm=includeNkinetics", 0]), [], []),
    "rBcvap": (_script(runcat="tuning", apopt="cvap", with_impscale=False), ["runname", "rBcvap"], []),
    "rBnoap": (_script(runcat="tuning", apopt="noap", with_impscale=False), ["runname", "rBnoap"], []),
    "rBrkap": (_script(runcat="tuning", apopt="rkap", with_impscale=False), ["runname", "rBrkap"], []),
    "rBfDr2p33": (_script(runcat="tuning", before=["Richards"], extra=["parm=roughnessfractal", 2.33]),
                  ["runname", "rBfDr2p33", "runlab", "Roughness fractal fDr=2.33"], []),
    "rBfDs2p7": (_script(runcat="tuning", before=["Richards"], extra=["parm=surfacefractal", 2.7]),
                 ["runname", "rBfDs2p7", "runlab", "Surface fractal fDs=2.7"], []),
    "rBhumsx5": (_script(), ["parm=kNhumsday", 0.00005], ["runname", "rBhumsx5", "runlab", "kNhums=%g" % 0.00005]),
    "rBhumsx10": (_script(), ["parm=kNhumsday", 0.00010], ["runname", "rBhumsx10", "runlab", "kNhums=%g" % 0.00010]),
    "rBdynex": (_script(), ["runname", "rBdynex", "runlab", "Dynamic exchanger",
                            "INCLUDE", [AUTHOR_DBDIR + "Benson1982exchcoeffsXb"]],
                ["includedynexchange", "nodynamicsurfaces"]),
    "rCcomp": (_script(fsname="CTL", RSAscale=1, impscale=1), [], []),
    "rCnoso": (_script(fsname="CTL", tropt="noso", RSAscale=1, impscale=1), [], []),
    "rCclosed": (_script(fsname="CTL", extra=["fixCO2gmol", 0]), [], []),
    "rFcomp": (_script(fsname="Fo", RSAscale=1, impscale=1), [], []),
    "rFnoso": (_script(fsname="Fo", tropt="noso", RSAscale=1, impscale=1), [], []),
    "rFclosed": (_script(fsname="Fo", extra=["fixCO2gmol", 0]), [], []),
    "rWcomp": (_script(fsname="Wo", RSAscale=1, impscale=1), [], []),
    "rWnoso": (_script(fsname="Wo", tropt="noso", RSAscale=1, impscale=1), [], []),
    "rWclosed": (_script(fsname="Wo", extra=["fixCO2gmol", 0]), [], []),
    # long runs, 2015-2070, with repeated feedstock applications
    **{f"rlBrep{n}": (_script(before=["fsrep", n, "nyears", 56, "parm=ceasefeedstockyear", 2020]),
                      ["parm=t1", "2015-01-01", "parm=t2", "2070-12-30", "parm=repeattimetable", [2020.0, 2071.0]], [])
       for n in (4, 23, 55)},
    "rlC": (None, ["parm=t1", "2015-01-01", "parm=t2", "2070-12-30", "parm=repeattimetable", [2020.0, 2071.0]],
            ["runname", "rlC", "runlab", "No feedstock", "nodust"]),
    # spin-up runs, 2000-2015, from the spin-up forcing files
    "rspinmean": (None, ["parm=t1", "2000-Jan-14", "parm=t2", "2015-Dec-10", "parm=usemeanforcings", 1,
                         "parm=CO2atm", 395, "HYDROmethod", "Richards"],
                  ["runname", "rspinmean", "runlab", "Spinup mean forcings", "nodust", "notreatments"]),
    "rspinstd": (None, ["parm=t1", "2000-Jan-14", "parm=t2", "2015-Dec-10", "HYDROmethod", "Richards"],
                 ["runname", "rspinstd", "runlab", "Spinup std", "nodust", "notreatments"]),
}

# parm values that rbcomp_settings adds for the rBcomp run itself
_RUN_PARMS = ("fsname", "implicit", "RSAscale", "impscale")


def base_settings(artemis=None):
    """defaultsettings.m for EF2: the rBcomp settings without the run's own parm values."""
    from .reference import rbcomp_settings

    s = rbcomp_settings(artemis)
    for k in _RUN_PARMS:
        s["parm"].pop(k, None)
    return s


def apply_dsarg(settings, dsarg):
    """setupdatastruc.m's argument loop on the settings dict (fields and parm values)."""
    s = copy.deepcopy(settings)
    parm = s["parm"]
    i = 0
    while i < len(dsarg):
        a = dsarg[i]
        nxt = dsarg[i + 1] if i + 1 < len(dsarg) else None
        if isinstance(a, str):
            if a in s and a != "parm":
                s[a] = nxt
                i += 2
                continue
            if a in parm:
                parm[a] = nxt
                i += 2
                continue
            low = a.lower()
            if low in ("runname", "runlab"):
                s[low] = nxt
            elif low == "fsname":
                parm["fsname"] = nxt
            elif low in ("usefiles", "include"):
                s["usefiles" if low == "usefiles" else "INCLUDE"] = nxt
            elif a.startswith("fsname="):
                parm["fsname"] = a[len("fsname="):]
            elif a.startswith("parm="):
                parm[a[len("parm="):]] = nxt
            elif a.startswith(("scale=", "scheme=")):
                grp, key = a.split("=", 1)
                s.setdefault(grp, {})[key] = nxt
            elif "=" in a:
                k, v = a.split("=", 1)
                s[k] = v
        i += 1
    return s


def scenario(name, artemis=None):
    """(settings for setup_data_struc, getphreeqcrun arguments) of a published EF2 run."""
    efargs, dsextra, runextra = SCRIPTS[name]
    # scripts that call setupdatastruc without EFrunsettings2 get defaultsettings.m only
    dsarg, runarg = ef_run_settings(*efargs) if efargs is not None else ([], [])
    return apply_dsarg(base_settings(artemis), dsarg + list(dsextra)), runarg + list(runextra)

# ds fields the run scripts set between setupdatastruc and getphreeqcrun
POST_DS = {
    "rBnonit": {"runname": "rBnonit", "runlab": "{\\bf BR}{\\color{red}-N kinetics}"},
}
POST_RUNNAME = {k: v["runname"] for k, v in POST_DS.items() if "runname" in v}

# rBdynexscript.m reads secondary_phases3.csv instead of secondary_phases.csv
_SPINFILES = [f for f in EF2_FILES if "feedstock" not in f and "_obs_" not in f]
_SPINFILES = [f.replace("EF_maize_B0_timetable.csv", "EF_maize_B0_timetable_forspinup.csv")
              .replace("EF_CLM5_2015-2070", "EF_CLM5_2000-2014") for f in _SPINFILES]
USEFILES = {
    "rspinmean": _SPINFILES, "rspinstd": _SPINFILES,
    "rBdynex": [f.replace("/secondary_phases.csv", "/secondary_phases3.csv") for f in EF2_FILES],
}

# parm values that differ between the published inputs and the archived scripts: the
# 2 June 2026 runs predate defaultsettings.m's change of sdnco from 1.45 to 1.46 (3 June 2026)
PUBLISHED_PARMS = {"rBhumsx5": {"sdnco": 1.45}, "rBhumsx10": {"sdnco": 1.45}}
