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


def stage_init_run(name, phreeqc_dir, workdir, artemis=None, flanagan=None):
    """Copy ``<name>.in_init`` into ``workdir/EF2`` with paths rewritten.

    Returns the path of the staged .in file.
    """
    artemis = Path(artemis or artemis_root())
    workdir = Path(workdir)
    ef2 = workdir / "EF2"
    ef2.mkdir(parents=True, exist_ok=True)
    (workdir / "EF_CLM5").mkdir(exist_ok=True)
    shutil.copy(artemis / "examples" / "EF2" / "Pellets.RATES", ef2)
    flanagan = Path(flanagan) if flanagan else artemis / "examples" / "EF_CLM5" / FLANAGAN
    if flanagan.exists():
        shutil.copy(flanagan, workdir / "EF_CLM5" / FLANAGAN)
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
