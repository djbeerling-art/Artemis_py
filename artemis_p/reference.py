"""Rerun the published ARTEMIS v1.1 PHREEQC inputs and compare with the published outputs.

The Zenodo archive ARTEMIS_v1p1_PHREEQC-input-output-files_GMD_R1 holds
``<run>.in_init`` (the first daily iteration) for about 50 runs, and full
outputs (``.out_init``, ``.sel_init``) for rBcomp_Blueridge. The inputs carry
the authors' absolute paths, so they are rewritten before running.
"""

import os
import shutil
from pathlib import Path

from . import phreeqc
from .treatments import treatment_rates_text
from .selected_output import compare_sel

AUTHOR_OUTDIR = "/mnt/parscratch/users/bo1llt//EF2//"
AUTHOR_DBDIR = "/users/bo1llt/soilmodel//phreeqcdbincludefiles/"

# rBcomp includes this file with a relative path. It is missing from the v1.1
# archive; the v1.0 archive has it in examples/EF_CLM5.
FLANAGAN = "Flanagan_BlueRidge_PHASESRATES_withimplicitphases"


def artemis_root():
    """The v1.1 soilmodel directory ($ARTEMIS, as in the MATLAB code)."""
    root = os.environ.get("ARTEMIS")
    if not root:
        raise EnvironmentError("set $ARTEMIS to the ARTEMIS v1.1 soilmodel directory")
    return Path(root)


def stage_init_run(name, phreeqc_dir, workdir, artemis=None, flanagan=None, pellets=("UAN28",), text=None):
    """Copy ``<name>.in_init`` into ``workdir/EF2`` with paths rewritten.

    ``Pellets.RATES`` is rebuilt from ``Pellets.dat`` for the ``pellets``
    treatments, because the copy in the v1.1 archive (examples/EF2) was written
    for a different run (CaCO3_trmt).

    ``text`` replaces the published input (e.g. one written by the Python port).

    Returns the path of the staged .in file.
    """
    artemis = Path(artemis or artemis_root())
    workdir = Path(workdir)
    ef2 = workdir / "EF2"
    ef2.mkdir(parents=True, exist_ok=True)
    (workdir / "EF_CLM5").mkdir(exist_ok=True)
    (ef2 / "Pellets.RATES").write_text(
        treatment_rates_text(artemis / "phreeqcdbincludefiles", "Pellets", list(pellets)))
    flanagan = Path(flanagan) if flanagan else artemis / "examples" / "EF_CLM5" / FLANAGAN
    if flanagan.exists():
        shutil.copy(flanagan, workdir / "EF_CLM5" / FLANAGAN)
    if text is None:
        text = (Path(phreeqc_dir) / f"{name}.in_init").read_text()
    text = phreeqc.relocate_input(text, {
        AUTHOR_OUTDIR: f"{ef2}/",
        AUTHOR_DBDIR: f"{artemis / 'phreeqcdbincludefiles'}/",
    })
    infile = ef2 / f"{name}.in"
    infile.write_text(text)
    return infile


def rerun_init(name, phreeqc_dir, workdir, dbname="THAmmOrgAlP.dat", exe=None, **kw):
    """Stage and run one published first-iteration input. Returns (RunResult, comparison or None)."""
    artemis = Path(kw.get("artemis") or artemis_root())
    infile = stage_init_run(name, phreeqc_dir, workdir, **kw)
    res = phreeqc.run(infile, artemis / "phreeqcdbincludefiles" / dbname, exe=exe)
    ref = Path(phreeqc_dir) / f"{name}.sel_init"
    cmp = compare_sel(infile.with_suffix(".sel"), ref) if ref.exists() and res.ok else None
    return res, cmp


# Input files of the published EF2 runs, as in ds.usefiles (relative to examples/EF2)
EF2_FILES = [
    "../EF_data/Richland_loess_Batestown_till_phases_combined.csv",
    AUTHOR_DBDIR + "/secondary_phases.csv",
    AUTHOR_DBDIR + "/feedstock_phases_pure.csv",
    AUTHOR_DBDIR + "/feedstock_phases_Lewis_THERMODDEM.csv",
    "../EF_data/Lewis_feedstocks_PSD.csv",
    "../EF_data/initlayers_Flanagan_EF_0-183cm.csv",
    "../EF_data/rainwater_Champaign_NTN_2021_initialChemistry.csv",
    "../EF_data/EF_maize_B0_timetable.csv",
] + [f"../EF_data/EF_CLM5_2015-2070_timeseries_{f}.csv"
     for f in ("soiltemp", "airtemp", "soilresp", "CO2atm", "snow", "infil_ET")]


def ef2_resolver(artemis=None):
    """Map EF2 input names (relative to examples/EF2, or under the authors'
    include directory) to local paths."""
    artemis = Path(artemis or artemis_root())
    ef2 = artemis / "examples" / "EF2"

    def resolve(f):
        if f.startswith(AUTHOR_DBDIR):
            return artemis / "phreeqcdbincludefiles" / f[len(AUTHOR_DBDIR):].lstrip("/")
        return (ef2 / f).resolve() if not Path(f).is_absolute() else Path(f)

    return resolve, ef2


RBCOMP_PHREEQCBLOCK = ("SOLUTION_SPECIES; H2O + 0.01e- = H2O-0.01; log_k -9.0; KNOBS; -convergence_tolerance 1e-7; "
                       "-tolerance 1e-14; -step_size 10.; -diagonal_scale true; -iterations 200; -pe_step_size 5.;")


def rbcomp_run_settings(dbdir=AUTHOR_DBDIR, outdir=AUTHOR_OUTDIR[:-1]):
    """The ds fields from EF2 defaultsettings.m and EFrunsettings2.m that getphreeqcrun reads
    for rBcomp. The defaults are the authors' paths, as written into the published inputs
    (getphreeqcrun adds a slash to outdir)."""
    dbdir = str(dbdir)
    return {"runname": "rBcomp", "runlab": "{\\bf BR}+UAN28+uptake+sorption", "dbdir": dbdir,
            "dbfile": "THAmmOrgAlP.dat", "plantratesfile": dbdir + "plantRATES5",
            "schemeratesfile": dbdir + "SWATRATESPHASEStunable", "outdir": str(outdir),
            "phreeqcblock": RBCOMP_PHREEQCBLOCK}


def rbcomp_settings(artemis=None):
    """Initial ds for the rBcomp_Blueridge run: EF2 defaultsettings.m,
    BRscalings.m and the run's own arguments in rBcompscript.m."""
    artemis = Path(artemis or artemis_root())
    parm = {
        "writetrRATESfiles": 0, "t1": "2015-01-01", "t2": "2020-12-31",
        "Feactfrac": 0.01, "SOCactfrac": 0.629, "CECopt": "CECcarboxylic", "setexaciditysp": "Al",
        "nativeRSAscale": 1e-6, "tillagedepth": 18, "fixCO2gmol": 10, "releaseCO2": 0,
        "sdnco": 1.46, "maxsut": 3.0, "kNhumsday": 3.5e-4, "fsname": "Blueridge",
        "implicit": ["Titanite_blueridge", "Hydroxyapatite_blueridge"], "RSAscale": 2, "impscale": 5,
    }
    return {"parm": parm, "EXCHmethod": "comp", "WHAMSA": 12, "SAmethod": "traditional",
            "HYDROmethod": "Richards", "ETmethod": "MIXflownut", "MIXdir": "topdown",
            "dbdir": artemis / "phreeqcdbincludefiles",
            "phasesratesfile": "../EF_CLM5/" + FLANAGAN}
