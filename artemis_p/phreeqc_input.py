"""Write the daily PHREEQC input files, as getphreeqcrun.m does.

``PhreeqcRun`` holds the state that getphreeqcrun keeps in its workspace and
shares with its nested functions: the settings (``stg``), the decomposition
scheme (``scheme``), the kinetics flags and the ``out`` structure.  Its methods
keep the names of the nested MATLAB functions (writeselectedoutput,
writeinfluent, initsoilcells, writekinetics, ...), and each writes to the
current input text exactly as the MATLAB fprintf calls do, including MATLAB's
``%d``-of-a-non-integer behaviour (see :mod:`artemis_p.mformat`).

Day indices ``i`` are 1-based as in MATLAB; arrays are indexed with ``i - 1``.
``run_day`` does one pass of getphreeqcrun's daily loop: write the input, run
PHREEQC, read and reshape the selected output, and update moles and reactive
surface areas for the next day.
"""

import datetime as _dt
import math
from functools import partial
from pathlib import Path

import numpy as np

from . import phreeqc
from .fixes import active, enabled as fixed
from .mformat import sprintf
from .obs import assessrun2
from .rsa import update_rsa
from .selected_output import get_selected_output, reshape_selected_output
from .totals import gettotalfields

MONTHNAME = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

CARBSITE = ["H_a", "H_b", "H_c", "H_d"]
PHENSITE = ["H_e", "H_f", "H_g", "H_h"]
DIPSITE = ["H_ab", "H_ad", "H_af", "H_ah", "H_bc", "H_be", "H_bg", "H_cd", "H_cf", "H_ch", "H_de", "H_dg"]
PUNCHELT = ["Ca", "Mg", "K", "Na", "Fe", "Al", "P", "S", "Si", "Ti"]


def default_settings():
    """The ``stg`` defaults set at the top of getphreeqcrun.m."""
    return {
        "includetreatments": 1, "includeNfix": 0, "ureamoles": 0, "urease": 0.01, "usedust": 0,
        "usesurface": 1, "WHAMSA": 12, "includeWHAM": 1, "equilibratesurface": 1, "dipsp": ["H2"],
        "includedynamicsurfaces": 0, "useexchange": 1, "pitzer_exchange_gammas": 0, "includedynexchange": 0,
        "includekinetics": 1, "kinrelease": 1, "usecvode": 0, "tol": 1e-5, "showwarnings": 0,
        "skipSAinKINETICS": 0, "includenative": 1, "includeNkinetics": 1, "includehumus": 1,
        "includemanure": 0, "includeresidues": 1, "includerootresidues": 1, "Nfixprotons": 0.01,
        "releaseCO2": 0, "addresnow": 0, "addrootresnow": 0, "addmannow": 0, "includeuptake": 1,
        "includeAmmuptake": 1, "nocounterion": 0, "defineallsp": 1, "minconc": 1e-7, "fixCO2gmol": 10,
        "fixO2gmol": 10, "correctwater": 1, "bottomup": 1, "topdown": 0, "usewater": 1, "waxwane": 0,
        "ndailymix": 1, "oldway": 0, "testway": 0, "skipfreezing": 0, "savemat": 1, "keepall": 0,
        "saveallfiles": 0, "ndays": 0, "ndata": 0,
    }


def default_scheme():
    """The ``scheme`` defaults set at the top of getphreeqcrun.m."""
    return {
        "maxsut": 1, "sdnco": 1.1, "cdn": 1.4, "kresday": 0.05, "kmanday": 0.05, "kNhumaday": 0.0003,
        "kNhumsday": 1e-5, "humactfrac": 0.02, "SPmolrat": 0.6, "resCgfrac": 0.58, "resNgfrac": 0.0055,
        "resPgfrac": 0.001, "O2gprotocol": "reset", "CO2gprotocol": "reset",
        "Pellets": {"k1": 0.8869, "n1": 0.59, "lab": "Uncoated", "k2": 0, "n2": 0},
    }


# getphreeqcrun workspace variables that getsettings may set from ds.parm
# (``exist(parmfld,'var')``); only the scalar run controls are listed.
_PARM_VARS = ("nkinsteps", "stepdivide", "badstepmax", "cvode_steps", "SAmethod", "EXCHmethod",
              "vleachthresh", "saveinterval", "usefixedinfluent", "chargefld", "charge0fld")


class TableValues:
    """Scalar lookups ``tv(df, row, col)`` that skip DataFrame.loc, which deep-copies
    the table's attrs (units, labels) on every call. Columns are cached until a
    different table object is passed; call ``clear()`` after changing a table in place."""

    def __init__(self):
        self.df, self.cols = None, {}

    def clear(self):
        self.df, self.cols = None, {}

    def __call__(self, df, row, col):
        if df is not self.df:
            self.df, self.cols = df, {}
        c = self.cols.get(col)
        if c is None:
            c = self.cols[col] = dict(zip(df.index, df[col].to_numpy()))
        return c[row]


def gasmoles(ds, CO2gmol=0, O2gmol=0):
    """gasmoles.m: moles of CO2 and O2 gas per cell (PV=nRT, or fixed values)."""
    run, cell, parm = ds["run"], ds["cell"], ds["parm"]
    if CO2gmol > 0:
        parm["CO2gmol"] = CO2gmol
    elif "CO2gmol" in parm:
        CO2gmol = parm["CO2gmol"]
    if O2gmol > 0:
        parm["O2gmol"] = O2gmol
    elif "O2gmol" in parm:
        O2gmol = parm["O2gmol"]
    R = 0.08207
    if str(parm.get("CO2gmolvol", "")).lower() == "air":
        airvol = np.asarray(cell["SAT"])[:, None] - run["targH2O"]
    elif "airpor" in run:
        airvol = run["airpor"] * np.asarray(cell["vsoil"])[:, None]
    else:
        airvol = np.asarray(cell["SAT"])[:, None] + 0 * run["targH2O"]
    RTV = airvol / R / (run["soiltemp"] + 273.15)
    run["CO2gmol"] = run["pCO2g"] * RTV
    if CO2gmol > 0:
        run["CO2gmol"] = run["CO2gmol"] * 0 + CO2gmol
    run["O2gmol"] = run["pO2g"] * RTV
    if O2gmol > 0:
        run["O2gmol"] = run["O2gmol"] * 0 + O2gmol
    lab, units = run.setdefault("label", {}), run.setdefault("units", {})
    lab.update(CO2gmol="Soil CO_2 gas", O2gmol="Soil O_2 gas", si_CO2g="CO_2 gas", si_O2g="O_2 gas",
               dCO2gmol="{\\Delta}CO_2 gas", dO2gmol="{\\Delta}O_2 gas")
    units.update(CO2gmol="mol cell^{-1}", O2gmol="mol cell^{-1}", si_CO2g="log_{10}(atm)",
                 si_O2g="log_{10}(atm)", dCO2gmol="mol cell^{-1} day^{-1}", dO2gmol="mol cell^{-1} day^{-1}")
    with np.errstate(divide="ignore"):
        run["si_CO2g"] = np.log10(run["pCO2g"])
        run["si_O2g"] = np.log10(run["pO2g"])
    run["dCO2gmol"] = np.diff(run["CO2gmol"], axis=1)
    run["dO2gmol"] = np.diff(run["O2gmol"], axis=1)
    return ds


def _str(p):
    return str(p) if p is not None else ""


def _matlab_date(t):
    return f"{t.day:02d}-{MONTHNAME[t.month - 1]}-{t.year}"


def _matlab_datetime(t):
    return f"{_matlab_date(t)} {t.hour:02d}:{t.minute:02d}:{t.second:02d}"


class PhreeqcRun:
    """getphreeqcrun.m up to (and including) writing each day's PHREEQC input.

    ``ds`` is the structure returned by setupdatastruc (``complete_setup``);
    keyword arguments are getphreeqcrun's command-line options, e.g.
    ``topdown=True``, ``ndays=10``, ``runname='rBcomp'`` or any ``stg``/``scheme``
    field.  ``now`` fixes the timestamps written into the files.
    """

    def __init__(self, ds, now=None, **opts):
        self.ds = ds
        self.tv = TableValues()
        self.engine = None  # artemis_p.iphreeqc.Engine when engine="iphreeqc"
        self.now = now or _dt.datetime.now()
        self.stg = default_settings()
        self.scheme = default_scheme()
        self.out = {"metadata": {}, "run": {}, "units": {}, "label": {}}
        self.chargefld = "Cl"
        self.charge0fld = "Cl"
        self.nyear = 3
        self.nmonth = 12 * self.nyear
        self.saveinterval = 365
        self.usemix = 1
        self.usefixedinfluent = 0
        self.EXCHmethod = "species"
        self.nkinsteps = 3
        self.stepdivide = 0
        self.badstepmax = 1000
        self.cvode_steps = 100
        self.SAmethod = "inbuilt"
        self.schemeratesfile = "SWATRATESPHASEStunable"
        self.decompscheme = "SWAT"
        self.inpercfld = "inperc"
        self.runname = ""
        self.outdir = "$SCRATCH/"
        self.dbfile = "THAmmOrg.dat"
        self.dbdir = "$ARTEMIS/phreeqcdbincludefiles/"
        self.phrdir = "$HOME/bin/"
        self.phasesratesfile = "basaltRATESPHASES5"
        self.plantratesfile = "plantRATES5"
        self.restart = 0
        self.append = 0
        self.fname = ""
        self.kinphase = []
        self.nkinphase = 0
        self.nreversible = 0
        self.vleachthresh = 1e-6
        self.phreeqcblock = None
        self.includefile = None
        self.kinsteps = None
        self.carbsite, self.phensite, self.dipsite = list(CARBSITE), list(PHENSITE), list(DIPSITE)

        self.getsettings()
        self._apply_options(opts)
        self._setup()

    # ------------------------------------------------------------------ settings
    def getsettings(self):
        ds, stg, scheme = self.ds, self.stg, self.scheme
        for key, attr in (("dbdir", "dbdir"), ("phrdir", "phrdir"), ("outdir", "outdir"),
                          ("phasesratesfile", "phasesratesfile"), ("plantratesfile", "plantratesfile"),
                          ("schemeratesfile", "schemeratesfile"), ("dbfile", "dbfile"),
                          ("SAmethod", "SAmethod"), ("EXCHmethod", "EXCHmethod"), ("phreeqcblock", "phreeqcblock")):
            if key in ds:
                setattr(self, attr, ds[key])
        if "ndays" in ds:
            stg["ndays"] = ds["ndays"]
        if "runname" in ds:
            self.runname = ds["runname"]
            self.out["metadata"]["runname"] = ds["runname"]
        if "runlab" in ds:
            self.out["metadata"]["runlab"] = ds["runlab"]
        if "WHAMSA" in ds:
            stg["WHAMSA"] = ds["WHAMSA"]
        if ds.get("MIXdir") == "topdown":
            stg["topdown"], stg["bottomup"] = 1, 0
        elif ds.get("MIXdir") == "bottomup":
            stg["topdown"], stg["bottomup"] = 0, 1
        run = ds["run"]
        stg["usedust"] = np.nanmax(run["treated"])
        if "mmH2O" in run and "targH2O" not in run:
            run["targH2O"] = run["mmH2O"]
        if "rootresidue" not in run:
            run["rootresidue"] = run["residue"] * 0
        if "manure" not in run:
            run["manure"] = run["treated"] * 0
        if "exptday" not in run:
            run["exptday"] = np.arange(ds["nrun"])
        parm = ds.setdefault("parm", {"nativeRSAscale": 1})
        if "nocounterions" in parm:
            stg["nocounterions"] = 1
        parm.setdefault("impscale", 1)
        if "minconc" in parm:
            stg["minconc"] = parm["minconc"]
        if "releaseCO2" in parm:
            stg["releaseCO2"] = parm["releaseCO2"]
        for f, v in parm.items():
            if f in scheme:
                scheme[f] = v
            elif f in stg:
                stg[f] = v
            elif f in _PARM_VARS:
                setattr(self, f, v)
                self.out["metadata"][f] = v
        for f, v in ds.get("scheme", {}).items():
            scheme[f] = v

    def _apply_options(self, opts):
        stg = self.stg
        for k, v in opts.items():
            if k in self.scheme:
                self.scheme[k] = v
            elif k in ("topdown", "bottomup"):
                stg["topdown"], stg["bottomup"] = (1, 0) if (k == "topdown") == bool(v) else (0, 1)
            elif k in stg:
                stg[k] = v
            elif k == "nyears":
                stg["ndays"] = 365 * v
            elif k in ("printoutput", "cellworkers", "engine"):
                # port options: printoutput=0 keeps only warnings and errors in the .out file;
                # cellworkers>1 runs the cells' kinetics in parallel; engine="iphreeqc" keeps the
                # PHREEQC state in memory (artemis_p.iphreeqc) instead of running the executable
                stg[k] = v
            elif k == "runname":
                self.runname = v
                self.out["metadata"]["runname"] = v
            elif k == "runlab":
                self.out["metadata"]["runlab"] = v
            elif k in ("name",):
                self.fname = v
            elif k in ("noappend", "restart", "restartfile"):
                # restart from a dump without the previous results (spin-up runs)
                self.restart, self.append = 1, 0
                self.restartdump = self._restartdump(v)
            elif k == "append":
                raise NotImplementedError("getphreeqcrun 'append' (reloading the previous run's results)")
            elif hasattr(self, k):
                setattr(self, k, v)
            else:
                raise ValueError(f"unrecognised getphreeqcrun option {k!r}")

    @staticmethod
    def _restartdump(name):
        """The dump file of getphreeqcrun's restart argument: a .dmp, a .mat beside it, or no extension."""
        name = str(name)
        if name.endswith(".dmp"):
            return name
        if name.endswith(".mat"):
            return name[:-4] + ".dmp"
        if "." not in Path(name).name:
            return name + ".dmp"
        raise ValueError(f"Unrecognised restart filename extension: {name}")

    def _setup(self):
        ds, stg, scheme, out = self.ds, self.stg, self.scheme, self.out
        self.istart, self.ilast = self.setrunindices()
        self.dbfile = _str(self.dbdir) + self.dbfile
        if "/" not in _str(self.plantratesfile):
            self.plantratesfile = _str(self.dbdir) + self.plantratesfile
        if self.includefile is not None and "/" not in self.includefile:
            self.includefile = _str(self.dbdir) + self.includefile
        self.ptargH2O = ds["run"][self.inpercfld][0, :]
        out["run"]["ptargH2O"] = self.ptargH2O
        self.outpercfld = self.inpercfld.replace("in", "out")
        self.decompsetup(self.decompscheme)
        out["metadata"]["scheme"] = scheme
        self.setoptions()

        pt = ds.get("phasetable")
        self.eqsolid = []
        if pt is not None:
            self.eqsolid = list(pt.loc[pt["PHREEQCphasetype"] == "equilibrium", "MATLABphase"])
        self.neqsolids = len(self.eqsolid)
        if stg["usedust"]:
            tt = ds["timetable"]
            wfs = [j for j in range(len(tt["runindex"]))
                   if self.istart <= tt["runindex"][j] <= self.ilast and tt["dosekg"][j] > 0
                   and tt["feedstock"][j] != ""]
            if not wfs:
                stg["usedust"] = 0
            self.fsname = "-".join(sorted(set(tt["feedstock"][j] for j in wfs))) if wfs else ""
        out["run"]["RSA0"] = {k: v.copy() for k, v in ds["run"]["RSA"].items()}
        out["run"]["RSA"] = {k: v.copy() for k, v in ds["run"]["RSA"].items()}
        out["run"]["M0"] = {k: v.copy() for k, v in ds["run"]["moles"].items()}
        out["run"]["M"] = {k: v.copy() for k, v in ds["run"]["moles"].items()}
        if pt is not None and "affinity" not in pt.columns:
            pt["affinity"] = 1.0
        self.nnut = len(ds["nut"])
        self.ileach = ds["ncells"] + 1
        self.setmetadata()
        self.setfilename()
        out["metadata"]["istart0"] = self.istart
        self.ilast = min(self.ilast, stg["ndata"])
        gasmoles(ds, CO2gmol=stg["fixCO2gmol"], O2gmol=stg["fixO2gmol"])

        if not stg["equilibratesurface"]:
            self._uncharged_surface_sites()
        scheme["Nhumaformula"] = "Amm 1 CO2 16.3314" if stg["releaseCO2"] else "Amm 1"
        if stg["includehumus"] and not self.restart:
            self.kinflag[:, self.istart - 1] = np.fmax(1, self.kinflag[:, self.istart - 1])
        self.step = np.ones(ds["nrun"]) * ds["datastep"]
        run, cell = ds["run"], ds["cell"]
        if stg["usewater"]:
            self.kgH2O = run["targH2O"]
            self.targH2O = run["targH2O"] * 0 + 1
            self.cellH2O = np.asarray(cell["FC"]) * 0 + 1
        else:
            self.targH2O = run["targH2O"]
            self.cellH2O = np.asarray(cell["targH2O"])
            self.kgH2O = self.targH2O / self.cellH2O[:, None] if stg["waxwane"] else self.targH2O * 0 + 1
        out["run"]["mmrat"] = run["targH2O"] + np.nan
        if not self.restart:
            out["run"]["mmrat"][:, self.istart - 1] = run["targH2O"][:, self.istart - 1] / np.asarray(cell["targH2O"])
        out["metadata"]["istart"] = self.istart
        out["metadata"]["ilast"] = self.ilast
        out["metadata"]["H2Oscale"] = self.cellH2O * 0 + 1 if stg["usewater"] else self.cellH2O
        self.initw = 1
        self.iter = 0
        self.init = 0 if self.restart else 1
        if self.restart:
            out["metadata"]["restartdump"] = self.restartdump
        out["run"]["current"] = {}
        out["_nrun"], out["_ncells"] = ds["nrun"], ds["ncells"]
        if "grams" in run:
            out["run"]["grams"] = {k: np.array(v, dtype=float) for k, v in run["grams"].items()}
        for f in ("errdesc", "nwarnings", "stopwatch"):
            out["run"][f] = {}
        if self.restart and not self.append:  # the restart day starts the kinetics afresh
            col = self.kinflag[:, self.istart - 1]
            col[col > 1] = 1

    def _uncharged_surface_sites(self):
        stg, cell = self.stg, self.ds["cell"]
        self.carbsite = [s + "H" for s in self.carbsite]
        self.phensite = [s + "H" for s in self.phensite]
        ndip = len(self.dipsite)
        xmol, totx2 = [], 0.0
        for sp in stg["dipsp"]:
            x = np.asarray(cell["exchange"].get(sp + "X2", 0.0) * np.ones(self.ds["ncells"]))
            xmol.append(x)
            if sp + "X2" in cell["exchange"]:
                totx2 = totx2 + x
        j, k = 0, 0
        for sp, x in zip(stg["dipsp"], xmol):
            nrnd = int(np.round(x[0] / np.atleast_1d(totx2)[0] * ndip))
            if nrnd == 0:
                continue
            k += nrnd
            for ii in range(j, k):
                self.dipsite[ii] = self.dipsite[ii] + sp
            j = k
        self.out["metadata"]["dipsite"] = self.dipsite

    def setrunindices(self):
        ds, stg = self.ds, self.stg
        stg["ndata"] = ds["nrun"]
        istart = getattr(self, "startday", 1)
        if stg["ndays"] <= 0:
            stg["ndays"] = stg["ndata"]
        ilast = min(istart + stg["ndays"] - 1, ds["nrun"])
        return istart, ilast

    def decompsetup(self, decompscheme):
        ds, stg, scheme = self.ds, self.stg, self.scheme
        run, cell = ds["run"], ds["cell"]
        scheme["kin"] = decompscheme + "kin"
        if decompscheme.lower() != "swat":
            raise NotImplementedError("only the SWAT decomposition scheme is implemented in getphreeqcrun")
        scheme["Nresx"] = 0.8
        scheme["orgrates"] = ["Nhums", "Nhuma", "Phum", "Residue", "Manurefresh"]
        scheme["Nrates"] = ["Nfix", "Nitrification", "Denitrification"]
        # MATLAB builds the list of included rates and then overwrites it with all of them
        self.withNfix = stg["includeNfix"] or not fixed(ds, "nfix_flag")
        if not self.withNfix:
            scheme["Nrates"] = scheme["Nrates"][1:]
        scheme["rates"] = scheme["orgrates"] + scheme["Nrates"]
        scheme["resCNgrat"] = scheme["resCgfrac"] / scheme["resNgfrac"]
        scheme["resCPgrat"] = scheme["resCgfrac"] / scheme["resPgfrac"]
        scheme["NCrat"] = scheme["resNgfrac"] / scheme["resCgfrac"] * 12.011 / 14.007
        scheme["PCrat"] = scheme["resPgfrac"] / scheme["resCgfrac"] * 12.011 / 30.974
        scheme["etatmp"] = 0.41 * (run["soiltemp"] - 5) / 10
        WP, FC = np.asarray(cell["WP"]), np.asarray(cell["FC"])
        sw25 = WP + 0.25 * (FC - WP)
        targ = run["targH2O"]
        etaSW = targ * 0 + 1
        for kk in range(ds["ncells"]):
            w = targ[kk] < sw25[kk]
            etaSW[kk, w] = (targ[kk, w] - WP[kk]) / (sw25[kk] - WP[kk])
        scheme["etaSW"] = etaSW
        mm = np.asarray(cell["middepth"]) * 10
        scheme["etamidz"] = 1 - mm / (mm + np.exp(4.706 - 0.0305 * mm))
        scheme["etanit"] = scheme["etatmp"] * etaSW
        scheme["etacec"] = 0.15
        scheme["etavol"] = scheme["etatmp"] * scheme["etamidz"][:, None] * scheme["etacec"]
        scheme["frnit"] = 1 - np.exp(-scheme["etanit"])
        scheme["frvol"] = 1 - np.exp(-scheme["etavol"])
        scheme["frnitvol"] = 1 - np.exp(-scheme["etanit"] - scheme["etavol"])
        with np.errstate(invalid="ignore", divide="ignore"):
            Nnit = scheme["frnit"] / (scheme["frnit"] + scheme["frvol"])
            Nvol = scheme["frvol"] / (scheme["frnit"] + scheme["frvol"])
        # MATLAB max/min ignore NaN, so 0/0 becomes the bound
        scheme["Nvol"] = np.fmin(np.fmax(Nvol, 0), 1)
        scheme["Nnit"] = np.fmin(np.fmax(Nnit, 0), 1)
        sut = 0.1 + 0.9 * np.sqrt(np.fmax(1e-7, targ - WP[:, None]) / (FC - WP)[:, None])
        sut = np.fmin(scheme["maxsut"], sut)
        scheme["sut"] = np.fmax(0.05, sut)
        if self.schemeratesfile.lower() == "swatratesphases2":
            scheme["void"] = np.asarray(cell["por"])[:, None] * (1 - targ / np.asarray(cell["SAT"])[:, None])
        scheme["pctorgC"] = np.repeat((np.asarray(cell["gSOC"]) / np.asarray(cell["gsoil"]) * 100)[:, None],
                                      ds["nrun"], axis=1)
        if "/" not in _str(self.schemeratesfile):
            self.schemeratesfile = _str(self.dbdir) + "/" + self.schemeratesfile
        scheme["ratesfile"] = _str(self.schemeratesfile)
        scheme["plantratesfile"] = _str(self.plantratesfile)
        scheme["nrates"] = len(scheme["rates"])

    def setoptions(self):
        ds, stg, scheme = self.ds, self.stg, self.scheme
        run, cell = ds["run"], ds["cell"]
        nrun, ncells = ds["nrun"], ds["ncells"]
        self.nnut = len(ds["nut"])
        if ds.get("plantdensity") is not None and ds["plantdensity"] <= 0:
            stg["includeuptake"] = 0
        if "plant" not in run or np.sum(run["plant"]) == 0:
            stg["includeuptake"] = 0
        pt = ds.get("phasetable")
        if not stg["includenative"] and pt is not None:
            # MATLAB uses || on arrays here, which errors unless the table has one row
            wpnat = pt["category"].str.lower().isin(["native", "soil", "parent"])
            pt.loc[wpnat, "PHREEQCphasetype"] = "ignore"
        if "manure" not in run:
            run["manure"] = np.zeros((ncells, nrun))
            stg["includemanure"] = 0
        if "rootresidue" not in run:
            run["rootresidue"] = np.zeros((ncells, nrun))
            stg["includerootresidues"] = 0
        if np.nanmax(run["rootmass"]) <= 0:
            stg["includerootresidues"] = 0
        if "residue" not in run:
            run["residue"] = np.zeros((ncells, nrun))
            stg["includeresidues"] = 0
        stg["includemanure"] = stg["includemanure"] * (np.nanmax(run["manure"]) > 0)
        stg["includerootresidues"] = stg["includerootresidues"] * (np.nanmax(run["rootresidue"]) > 0)
        stg["includeresidues"] = stg["includeresidues"] * (np.nanmax(run["residue"]) > 0)
        if "gSOC" not in cell:
            stg["includehumus"] = 0
        if not stg["usedust"]:
            from .phasetable import check_phase_type
            # updates ds["phasetable"] in place
            check_phase_type(None, ds, ignorefeedstock="all")
            pt = ds["phasetable"]
        if pt is not None:
            kin = pt["PHREEQCphasetype"] == "kinetic"
            self.kinphase = list(pt.loc[kin, "MATLABphase"])
            self.nreversible = int(np.sum(kin & (pt["reversible"] > 0)))
        self.nkinphase = len(self.kinphase)
        stg["includekinetics"] = (stg["usedust"] + stg["includeuptake"] + stg["includeresidues"]
                                  + stg["includehumus"] + stg["includeNkinetics"] + stg["includenative"]
                                  + stg["includemanure"] + stg["includerootresidues"] + stg["includetreatments"])
        if not stg["includekinetics"]:
            for f in ("usedust", "includeuptake", "includeresidues", "includehumus", "includeNkinetics",
                      "includenative", "includemanure", "includerootresidues", "includetreatments"):
                stg[f] = 0
            self.nkinphase = 0
        if "N" not in ds["solute"] and "Amm" not in ds["solute"]:
            stg["includeNkinetics"] = 0
        if not stg["includemanure"]:
            stg["addmannow"] = 0
        ds.setdefault("ETmethod", "MIXflow")
        if not stg["useexchange"]:
            self.EXCHmethod = "none"
        if not stg["usedust"]:
            run["treated"] = run["treated"] * 0

        # kinetics flags: 0 none, 1 initialise, 2 USE existing, 3 KINETICS_MODIFY
        kinflag = np.zeros((ncells, nrun))
        xrun = np.arange(1, nrun + 1)
        inrun = (xrun >= self.istart) & (xrun <= self.ilast)
        kinflag[:, inrun] += 1
        treated = np.cumsum(run["treated"]) * stg["usedust"]
        self.planted = run["plant"] * stg["includeuptake"]
        if ds.get("trnutflag") is None:
            trn = np.zeros(len(ds["nut"]))
            if ds["ETmethod"].lower() == "mixflownut":
                trn = trn + 1
            ds["trnutflag"] = trn
        self.trnutfld = [f"p_{n}uptake" if ds["trnutflag"][k] else None for k, n in enumerate(ds["nut"])]
        self.trnutfld.append("p_Ammuptake")
        till = np.asarray(cell["tillagefrac"])
        if stg["includeresidues"] and np.sum(stg["addresnow"]) > 0:
            run["residue"][:, self.istart - 1] += stg["addresnow"] * till
        if stg["includerootresidues"] and np.sum(stg["addrootresnow"]) > 0:
            run["rootresidue"][:, self.istart - 1] += stg["addrootresnow"]
        resid = np.cumsum(run["residue"][0]) * stg["includeresidues"]
        rootresid = np.cumsum(run["rootresidue"], axis=1) * stg["includerootresidues"]
        if stg["addmannow"] > 0:
            run["manure"][:, self.istart - 1] = stg["addmannow"] * till
        manfert = np.zeros(nrun)
        if stg["includemanure"]:
            manfert = np.cumsum(run["manure"][0]) * stg["includemanure"]
        if stg["includeuptake"] or stg["includehumus"] or stg["includeNkinetics"] or self.nreversible > 0:
            wtill = np.arange(ncells)
        else:
            wtill = np.flatnonzero(till > 0)
        w = np.flatnonzero(((treated + self.planted + resid + manfert + stg["includehumus"]) != 0) & inrun)
        if stg["includekinetics"] and len(wtill) > 0 and len(w) > 0:
            kinflag[np.ix_(wtill, w)] = 2
            if not self.restart:
                kinflag[:, w[0]] = 1
            if (stg["includehumus"] or stg["includeresidues"] or stg["includerootresidues"] or stg["includeuptake"]
                    or stg["includemanure"] or stg["includetreatments"] or stg["skipSAinKINETICS"]
                    or self.nreversible > 0):
                w = np.flatnonzero(kinflag[0] > 1)
            else:
                w = np.flatnonzero(((resid + treated + manfert + rootresid[0]) > 0) & (kinflag[0] > 1))
            if len(w) > 0:
                kinflag[np.ix_(wtill, w)] = 3
            if stg["includeresidues"] or stg["includerootresidues"]:
                scheme["residueformula"] = self._residueformula(np.flatnonzero(run["residue"][0] > 0), "Residue")
            scheme["mnsolidfrac"] = 0
            scheme["mpsolidfrac"] = 0
            if stg["includemanure"]:
                scheme["manureformula"] = self._residueformula(np.flatnonzero(run["manure"][0] > 0), "Manure")
        self.kinflag = kinflag

        self.sp = list(ds["solute"])
        self.elt = [ds["solute"][s]["elt"] for s in self.sp]
        if self.usefixedinfluent:
            run.pop("influent", None)
        if "influent" not in run:
            run.setdefault("label", {})["influent"] = ds["influent"]["label"]["influent"]
            run.setdefault("units", {})["influent"] = ds["influent"]["units"]["msp"]
            run["influent"] = {s: np.zeros(nrun) + v for s, v in ds["influent"]["msp"].items()}
        self.spi = list(run["influent"])
        if self.charge0fld not in run["influent"] and self.chargefld:
            run["influent"][self.chargefld] = np.ones(nrun)
        if stg["includehumus"] and "gSON" not in cell and "gSOC" in cell:
            cell["gSON"] = np.asarray(cell["gSOC"]) * 0.1

    def _residueformula(self, wfirst, which):
        """Residue (or manure) KINETICS formula; MATLAB sets the loop index i to
        the first day with residue (or manure) to choose the year of nutrient ratios."""
        ds, stg, scheme = self.ds, self.stg, self.scheme
        run = ds["run"]
        i = wfirst[0] if len(wfirst) else self.istart - 1
        f = sprintf("Amm %g P %g", scheme["Nresx"] * scheme["NCrat"], scheme["Nresx"] * scheme["PCrat"])
        if stg["releaseCO2"]:
            f += sprintf(" CO3-2 %g", scheme["Nresx"])
        for nut in ds["nut"]:
            if nut in ("N", "P"):
                continue
            rat = np.atleast_2d(run["nutrat"][nut])[0]
            # MATLAB's manure branch uses && here, which errors on arrays; use & for both
            w = (run["year"] == run["year"][i]) & (rat > 0)
            if not w.any():
                continue
            f += sprintf(" %s %g", nut, scheme["NCrat"] * np.mean(rat[w]))
        return sprintf("%s Nhuma %g Phum %g", f, scheme["NCrat"] * (1 - scheme["Nresx"]),
                       scheme["PCrat"] * (1 - scheme["Nresx"]))

    def setmetadata(self):
        ds, out, stg = self.ds, self.out, self.stg
        md = out["metadata"]
        if not md.get("runname"):
            md["runname"] = ds.get("expID", "rtest")
            if stg["oldway"]:
                md["runname"] += "_old"
            if stg["testway"]:
                md["runname"] += "_test"
        md.setdefault("runlab", md["runname"])
        md["settings"] = dict(stg)
        md["perctype"] = "MIX" if self.usemix else "STAGNANT"
        md["chargefld"] = self.chargefld
        md["charge0fld"] = self.charge0fld
        md["dbfile"] = self.dbfile
        md["ileach"] = self.ileach
        if self.phreeqcblock is not None:
            md["phreeqcblock"] = self.phreeqcblock
        if "expname" not in ds:
            ds["expname"] = ds.get("expID", md["runname"])
        md["title"] = sprintf("Experiment %s", ds["expname"])
        if self.includefile is not None:
            md["includefile"] = self.includefile
        self.setruninfo()

    def setruninfo(self):
        """The parts of setruninfo that change settings, plus out.runinfo's including/excluding
        text (gettotalfields looks for ' dust' in it). The rest only builds text for stdout."""
        stg, md = self.stg, self.out["metadata"]
        stg["skipSAinKINETICS"] = 0 if self.SAmethod.lower() == "inbuilt" else 1
        md["SAmethod"] = self.SAmethod
        md["settings"]["skipSAinKINETICS"] = stg["skipSAinKINETICS"]
        md["EXCHmethod"] = self.EXCHmethod
        if "treatments" not in self.ds:
            stg["includetreatments"] = 0
        ri = self.out.setdefault("runinfo", {})
        ri["treatments"] = " ".join(np.atleast_1d(self.ds["treatments"])) if stg["includetreatments"] else "none"
        inc, exc = "Including ", "Excluding "
        if stg["useexchange"]:
            inc += "EXCHANGE " + ("(+Xb) " if stg["includedynexchange"] else "(static only) ")
        else:
            exc += "EXCHANGE "
        if stg["usesurface"]:
            inc += "SURFACE " + ("(+dynamic) " if stg["includedynamicsurfaces"] else "(static only) ")
        else:
            exc += "SURFACE "
        if stg["includeuptake"]:
            inc += "uptake "
            if stg["includeAmmuptake"] == 0:
                exc += "Ammuptake "
        else:
            exc += "uptake "
        if stg["nocounterion"]:
            exc += "counterions "
        else:
            inc += "counterions "
        if stg["includeNkinetics"]:
            inc += "N-kinetics "
        else:
            exc += "N-kinetics "
        if stg["includenative"]:
            inc += "native minerals"
        else:
            exc += "native minerals "
        if self.nkinphase > 0:
            treated = np.sum(np.ravel(self.ds["run"]["treated"])[self.istart - 1:self.ilast])
            if stg["usedust"] and treated:
                inc += " dust " + self.fsname
                md["phasesratesfile"] = self.phasesratesfile
                md["kinstoichmethod"] = "-formula" if stg["kinrelease"] else "PHASES"
                ri["kinstoichmethod"] = sprintf("stoichiometry from %s, using %s)", md["kinstoichmethod"],
                                                self.phasesratesfile)
            elif stg["usedust"]:
                exc += sprintf(" dust (%s, pre-treatment)", self.fsname)
            else:
                exc += " dust"
        if not stg["equilibratesurface"]:
            exc += " SURFACE equilibration"
        for flag, name in (("includehumus", "humus"), ("includeresidues", "residues"),
                           ("includerootresidues", "root residues"), ("includemanure", "manure")):
            if stg[flag]:
                inc += f" {name} "
            else:
                exc += f" {name} "
        ri["including"], ri["excluding"] = inc, exc

    def setfilename(self):
        ds, md, stg = self.ds, self.out["metadata"], self.stg
        outdir = _str(self.outdir) + "/"
        fname = self.fname + "_" if self.fname else ""
        fname = outdir + fname
        if md["runname"] and md["runname"] not in fname:
            fname += md["runname"].replace(" ", "") + "_"
        if stg["usedust"]:
            fname += self.fsname + "_"
        elif not md["runname"]:
            fname += "X_"
        else:
            fname += "nodust_"
        if stg["keepall"]:
            fname += "all_"
        run = ds["run"]
        if "month" in run:
            startstr = sprintf("%s%d_", MONTHNAME[int(run["month"][self.istart - 1]) - 1],
                               run["year"][self.istart - 1])
            laststr = sprintf("%s%d_", MONTHNAME[int(run["month"][self.ilast - 1]) - 1],
                              run["year"][self.ilast - 1])
        else:
            startstr = sprintf("Day%d", run["exptday"][self.istart - 1])
            laststr = sprintf("-%d_", run["exptday"][self.ilast - 1])
        md["settings"]["startstr"], md["settings"]["laststr"] = startstr, laststr
        fname += startstr
        if startstr != laststr:
            fname += laststr
        fname += _matlab_date(self.now)
        self.fname = fname
        md["fname"] = fname

    # ------------------------------------------------------------- daily input
    def _w(self, fmt, *args):
        self._buf.append(sprintf(fmt, *args))

    def write_input(self, i, now=None):
        """Text of the PHREEQC input for day ``i`` (1-based), as getphreeqcrun writes ``<fname>.in``."""
        ds, stg, out = self.ds, self.stg, self.out
        run = ds["run"]
        self.i = i
        self.iter += 1
        t = now or self.now
        if self.ilast - self.istart + 1 > 365:
            current = sprintf("iteration %d doy %d %d %s %d starting: %s", self.iter, run["doy"][i - 1],
                              run["day"][i - 1], MONTHNAME[int(run["month"][i - 1]) - 1], run["year"][i - 1],
                              _matlab_datetime(t))
        else:
            current = sprintf("iteration %d Day %d starting: %s", self.iter, run["exptday"][i - 1],
                              _matlab_datetime(t))
        out["run"]["current"][i] = current
        self._buf = []
        w = self._w
        w("#PHREEQC run for %s\n", current)
        w("#Experiment: %s %s runname: %s\n", ds["expname"], out["metadata"]["runlab"], out["metadata"]["runname"])
        w("PRINT\n")
        # printoutput=0 (not in MATLAB) leaves the .out file with just warnings and errors
        w("    -reset %s\n", "true" if stg.get("printoutput", 1) else "false")
        w("    -user_print true\n")
        self.writeselectedoutput()
        self.writepunch()
        if self.includefile is not None:
            w("INCLUDE$ %s\n", self.includefile)
        w("INCLUDE$ %s\n", _str(self.phasesratesfile))
        if "INCLUDE" in ds:
            for f in ([ds["INCLUDE"]] if isinstance(ds["INCLUDE"], str) else ds["INCLUDE"]):
                w("INCLUDE$ %s\n", _str(f))
        if (stg["includeuptake"] or stg["includehumus"] or stg["includeresidues"] or stg["includerootresidues"]
                or stg["includemanure"]):
            w("INCLUDE$ %s\n", _str(self.plantratesfile))
            if self.scheme["ratesfile"]:
                w("INCLUDE$ %s\n", self.scheme["ratesfile"])
        if self.restart:
            w("INCLUDE$ %s\n\n", self.restartdump)
            self.restart = 0
        elif not self.init:
            w("INCLUDE$ %s.dmp\n\n", self.fname)
        w("TITLE %s %s\n", out["metadata"]["title"], current)
        if self.phreeqcblock is not None:
            w("%s\n", self.phreeqcblock)
        self.writeinfluent(0, self.initw)
        if self.init:
            self.initsoilcells()
        else:
            self.updatesoilcells()
        if not run[self.inpercfld][0, i - 1] > 0:
            w("# No infiltration into topsoil during this timestep (%s)\n", current)
        ncells, ileach = ds["ncells"], self.ileach
        if self.usemix:
            mmrat = out["run"]["mmrat"]
            for k in range(1, ncells + 1):
                if "soil" in out and stg["correctwater"]:
                    mmrat[k - 1, i - 1] = self.kgH2O[k - 1, i - 1] / out["soil"]["mass_H2O"][k - 1, i - 2]
                    w("MIX %d; %d %f; SAVE solution %d; END # Correct amount of water to target %f kg\n ",
                      k, k, mmrat[k - 1, i - 1], k, self.kgH2O[k - 1, i - 1])
                elif np.isnan(mmrat[k - 1, i - 1]):
                    mmrat[k - 1, i - 1] = 1
            if stg["bottomup"]:
                self.percleachate()
                order = range(ncells, 0, -1)
            else:
                order = range(1, ncells + 1)
            mixfrac = run.get("mixfrac")
            for k in order:
                if mixfrac is not None:
                    if np.sum(mixfrac[k - 1, :, i - 1]) == 0:
                        w("# No %s\n", ds["MIXture"][k - 1])
                        continue
                    w("%s; ", ds["MIXture"][k - 1])
                    for kk in range(int(ds["nmixsoln"]), 0, -1):
                        if mixfrac[k - 1, kk - 1, i - 1] == 0:
                            continue
                        w("%d %f; ", ds["MIXsoln"][kk - 1], mixfrac[k - 1, kk - 1, i - 1])
                w("SAVE solution %d; END # Flow MIX\n", k)
            w("\n")
            if stg["topdown"]:
                self.percleachate()
        if "rhizon" in run:
            for k in range(1, ncells + 1):
                rh = run["rhizon"][k - 1, i - 1]
                if rh <= 0 or np.isnan(rh):
                    continue
                fc = run["cellFC"][k - 1, i - 1]
                w("MIX %d; %d %f; SAVE solution %d; END # Remove rhizon sample water (%g L)\n ",
                  k, k, (fc - rh) / fc, k, rh)
        for k in range(ileach - 1, 0, -1):
            if run["cellET"][k - 1, i - 1] == 0:
                w("# No evapotranspiration from cell %d timestep %d\n", k, i)
            self._use_cell(k)
            et = ds["ETmethod"]
            if et == "REACTION":
                w("REACTION %d # Evapotranspiration %g %s\n", k, run["cellET"][k - 1, i - 1], ds["units"]["mmscl"])
                w("     H2O -1.0\n")
                w("     %f moles\n", 55.506 * run["cellET"][k - 1, i - 1] / self.targH2O[k - 1, i - 1])
            elif "MIXflow" in et and "ETmix" in run:
                w("MIX %d; %d %f; # %s\n", k, k, run["ETmix"][k - 1, i - 1], ds["ETmethodlabel"])
            elif et == "MIXpure":
                w("SOLUTION %d # Pure water removed by ET\n", 100 + k)
                w("    %10s %f\n", "temp", run["soiltemp"][k - 1, i - 1])
                w("    %10s %f\n", "pH", 7)
                w("    %10s %g\n", "-water", run["cellET"][k - 1, i - 1] / self.targH2O[k - 1, i - 1])
                w("MIX %d; %d 1; %d -1\n", k, k, k + 100)
            if stg["useexchange"]:
                w("SAVE exchange %d; ", k)
            if stg["usesurface"]:
                w("SAVE surface %d; ", k)
            w("SAVE solution %d; SAVE equilibrium_phases %d\n", k, k)
            w("END\n")
        if stg["includekinetics"]:
            for k in range(ileach - 1, 0, -1):
                self._use_cell(k)
                self.k = k
                self.writekinetics()
                w("SAVE solution %d; SAVE equilibrium_phases %d\n", k, k)
                if stg["useexchange"]:
                    w("SAVE exchange %d; ", k)
                if stg["usesurface"]:
                    w("SAVE surface %d; ", k)
                w("END\n")
        w("DUMP\n")
        w("    -file %s.dmp\n", self.fname)
        w("    -all\n")
        w("END\n")
        return "".join(self._buf)

    # ------------------------------------------------------------- daily loop
    def run_day(self, i, now=None, exe=None, cwd=None):
        """One pass of getphreeqcrun's loop for day ``i``. Returns the PHREEQC RunResult.

        The input is written to ``<fname>.in`` (so ``outdir`` must be a local
        directory); ``cwd`` is where PHREEQC runs, which matters for relative
        INCLUDE$ paths. On the first day the files are also copied to ``*_init``.
        """
        import shutil
        import time

        ds, stg, out = self.ds, self.stg, self.out
        run = ds["run"]
        tic = time.time()
        self.tv.clear()
        text = self.write_input(i, now=now)
        fname = Path(self.fname)
        infile = fname.with_name(fname.name + ".in")
        infile.write_text(text)
        runner = phreeqc.run
        if stg.get("engine", "exe") == "iphreeqc":
            if self.engine is None:
                from .iphreeqc import Engine
                self.engine = Engine(self.dbfile, workers=stg.get("cellworkers", 1),
                                     printoutput=bool(stg.get("printoutput", 1)))
            runner = self.engine
        elif stg.get("cellworkers", 1) > 1:
            runner = partial(phreeqc.run_cells_parallel, workers=stg["cellworkers"])
        res = runner(infile, self.dbfile, outfile=fname.with_name(fname.name + ".out"),
                     screenout=fname.with_name(fname.name + ".screenout"), exe=exe, cwd=cwd)
        out["run"]["nwarnings"][i] = len(res.warnings)
        if res.errors:
            out["run"]["errdesc"][i] = "\n".join(res.errors)
            return res
        out["run"]["errdesc"][i] = "PHREEQC OK"
        m = get_selected_output(fname.with_name(fname.name + ".sel"))
        if "soln" not in m:
            return res
        m = self.addtoselectedoutput(m)
        reshape_selected_output(m, out, i)
        if self.init or self.iter == 1:
            for ext in ("in", "out", "dmp", "sel"):
                src = fname.with_name(f"{fname.name}.{ext}")
                if src.exists():
                    shutil.copy(src, fname.with_name(f"{fname.name}.{ext}_init"))
        soil = out["soil"]
        if stg["includeuptake"] and ds["ETmethod"].lower() == "mixflownut":
            k = 1  # MATLAB tests trnutflag(k) with the cell index left over from the loops above (1)
            for jnut in range(self.nnut):
                if ds["trnutflag"][jnut if fixed(ds, "passive_uptake_flag") else k - 1] == 0:
                    continue
                f = self.trnutfld[jnut]
                prev = soil[f][:, i - 2] * self.planted[i - 1] if i > self.istart else 0
                soil.setdefault(f, np.zeros((ds["ncells"], ds["nrun"])))[:, i - 1] = prev + soil["d" + f][:, i - 1]
            f = "p_Ammuptake"
            prev = soil[f][:, i - 2] * self.planted[i - 1] if i > self.istart else 0
            soil.setdefault(f, np.zeros((ds["ncells"], ds["nrun"])))[:, i - 1] = prev + soil["d" + f][:, i - 1]
            soil.setdefault("p_DINuptake", np.zeros((ds["ncells"], ds["nrun"])))[:, i - 1] = (
                soil["p_Nuptake"][:, i - 1] + soil["p_Ammuptake"][:, i - 1])
        for p in self.kinphase[: self.nkinphase]:
            out["run"]["M"][p][:, i - 1] = soil["k_" + p][:, i - 1] * self.targH2O[:, i - 1]
            if "Implicit_" in p:
                self.updateexplicitgrams(p)
        update_rsa(out, self.SAmethod, ds, i, self.kinphase[: self.nkinphase])
        out["run"]["stopwatch"][i] = time.time() - tic
        self.init = 0
        return res

    def finish(self):
        """getphreeqcrun's cleanup: end-of-run totals (gettotalfields), the comparison
        with observations (assessrun2) and run metadata."""
        ds, out, md = self.ds, self.out, self.out["metadata"]
        out["run"]["kinflag"] = self.kinflag
        md["scheme"] = self.scheme
        idone = len(out["run"]["stopwatch"])
        if idone < md["ilast"]:
            md["ilast"] = self.ilast = idone
        if idone < ds["nrun"]:  # MATLAB's output arrays grow day by day; ours are preallocated
            for grp in ("soil", "influent", "leachate"):
                for f, v in out[grp].items():
                    if isinstance(v, np.ndarray) and v.ndim == 2 and v.shape[1] == ds["nrun"]:
                        out[grp][f] = v[:, :idone]
        gettotalfields(out, ds)
        out["total"]["label"]["exptday"] = "Day of experiment"
        out["label"]["exptday"] = "Day of experiment"
        if "obs" in ds:
            assessrun2(out, ds)
        start = md.get("istart0", self.istart)
        mins = sum(v for d, v in out["run"]["stopwatch"].items() if start <= d <= idone) / 60
        errs = [d for d, e in out["run"]["errdesc"].items() if e != "PHREEQC OK"]
        what = f"ERROR at day {idone - self.istart + 1}" if errs else f"{idone - self.istart + 1} days"
        md["runtime"] = sprintf("Runtime for %s (%s, %d warnings): %g minutes", md["runname"], what,
                                sum(out["run"]["nwarnings"].values()), mins)
        md["dumpfile"] = self.fname + ".dmp"
        if self.engine is not None:  # the state was in memory: write the dump a restart reads
            self.engine.write_dump(md["dumpfile"])
            self.engine.close()
            self.engine = None
        md["fixes"] = active(ds)
        return out

    def addtoselectedoutput(self, m):
        ds, out, i = self.ds, self.out, self.i
        run = ds["run"]
        soln = np.asarray(m["soln"])
        mzer = soln * 0.0
        m["iter"] = mzer + self.iter
        m["runday"] = mzer + i
        for f in ("decyear", "year", "doy", "month", "day", "exptday"):
            m[f] = mzer + run[f][i - 1]
        if all(f in m for f in ("m_HCO3", "m_CO3", "m_CO2")):
            m["m_DIC"] = m["m_CO2"] + m["m_HCO3"] + m["m_CO3"]
        m["cateq"] = 2 * (m["Ca"] + m["Mg"]) + m["Na"] + m["K"]
        m["m_DIN"] = m["N"] + m["Amm"]
        for sp in ("N", "Amm"):
            if sp not in m:
                continue
            m["DIN"] = m[sp] if "DIN" not in m else m["DIN"] + m[sp]
            kf = f"k_{sp}uptake"
            if kf not in m:
                continue
            # only the first species found is used, so k_DINuptake is the nitrate uptake alone
            if "k_DINuptake" in m and fixed(self.ds, "din_uptake"):
                m["k_DINuptake"] = m["k_DINuptake"] + m[kf]
                m["dk_DINuptake"] = m["dk_DINuptake"] + m["d" + kf]
            if "k_DINuptake" not in m:
                m["k_DINuptake"] = m[kf]
                m["dk_DINuptake"] = m["d" + kf]
        m["targH2O"] = mzer.copy()
        m["targH2O"][soln == 0] = run["inperc"][0, i - 1]
        for f in self.trnutfld[: self.nnut]:
            m["d" + f] = np.zeros(len(soln))
        m["dp_Ammuptake"] = np.zeros(len(soln))
        m["dp_DINuptake"] = np.zeros(len(soln))
        mixflownut = ds["ETmethod"].lower() == "mixflownut"
        for jsol in range(1, ds["ncells"] + 1):
            w = soln == jsol
            if not w.any():
                continue
            m["targH2O"][w] = run["targH2O"][jsol - 1, i - 1]
            if mixflownut:
                tr = run["uptake"]["transpiration"][jsol - 1, i - 1]
                for inut, n in enumerate(ds["nut"]):
                    m["d" + self.trnutfld[inut]][w] = tr * m[n][w]
                m["dp_Ammuptake"][w] = tr * m["Amm"][w]
        if mixflownut:
            m["dp_DINuptake"] = m["dp_Nuptake"] + m["dp_Ammuptake"]
        w = np.flatnonzero(soln == self.ileach)
        operc = run[self.outpercfld][ds["ncells"] - 1]
        if len(w):
            m["targH2O"][w] = operc[i - 1]
        if np.sum(operc[:i]) == 0:
            m["mass_H2O"][w] = 0
        if operc[i - 1] < self.vleachthresh:
            m["pH"][w] = np.nan
        if run[self.inpercfld][0, i - 1] <= 0 and fixed(ds, "leachate_ph"):
            m["pH"][soln == 0] = np.nan
            m["mass_H2O"][w] = 0
        elif run[self.inpercfld][0, i - 1] <= 0:
            # MATLAB writes out.pH (not m.pH), and uses the leachate rows here, not the influent's
            if len(w):
                pH = out.get("pH", np.zeros(0))
                if pH.size < w.max() + 1:
                    pH = np.concatenate([pH, np.zeros(w.max() + 1 - pH.size)])
                pH[w] = np.nan
                out["pH"] = pH
            m["mass_H2O"][w] = 0
        return m

    def updateexplicitgrams(self, impname):
        ds, out, i = self.ds, self.out, self.i
        if i >= self.ilast:
            return
        pt = ds["phasetable"]
        expname = impname.replace("Implicit", "Explicit")
        fs = self.tv(pt, impname, "feedstock")
        soil = out["soil"]
        grams = out["run"]["grams"][expname][:, i].copy()
        same = (pt["feedstock"].str.lower() == str(fs).lower()) & (pt["PHREEQCphase"] != impname)
        for ph in pt.index[same & (pt["PHREEQCphasetype"].str.lower() == "kinetic")]:
            grams = grams + soil["k_" + self.tv(pt, ph, "MATLABphase")][:, i - 1] * self.tv(pt, ph, "molmass")
        for ph in pt.index[same & (pt["PHREEQCphasetype"].str.lower() == "equilibrium")]:
            grams = grams + soil[self.tv(pt, ph, "MATLABphase")][:, i - 1] * self.tv(pt, ph, "molmass")
        out["run"]["grams"][expname][:, i] = grams

    def _use_cell(self, k):
        self._w("USE solution %d; USE equilibrium_phases %d; ", k, k)
        if self.stg["useexchange"]:
            self._w("USE exchange %d; ", k)
        if self.stg["usesurface"]:
            self._w("USE surface %d; ", k)
        self._w("\n")

    def writeselectedoutput(self):
        ds, stg, out = self.ds, self.stg, self.out
        w = self._w
        cell = ds["cell"]
        w("\nSELECTED_OUTPUT \n")
        w("    -file %s.sel\n", self.fname)
        w("    -temperature true\n")
        w("    -molalities CO2 HCO3- CO3-2 AmmH+ NO3- Al+3")
        totfld = "Amm Si Urea"
        if stg["useexchange"]:
            exchfld = ["CaX2", "MgX2", "KX", "NaX", "AmmHX", "HX", "AlX3"]
            totfld = "X Amm Si Urea"
            # Only on the first call: later days keep the unsorted default list
            if "exchfld" not in cell and "exchange" in cell:
                cell["exchfld"] = list(cell["exchange"])
                exchfld = sorted(set(exchfld + cell["exchfld"]))
            for f in exchfld:
                w(" %s", f)
            pt = ds.get("phasetable")
            if stg["includedynexchange"] and pt is not None and "exchange" in pt.columns:
                if pt["exchange"].str.contains("Xb").any():
                    totfld += " Xb"
                    for f in exchfld:
                        w(" %s", f.replace("X", "Xb"))
        w("\n")
        w("    -totals %s", totfld)
        if "msp" in ds["influent"]:
            if self.chargefld not in ds["influent"]["msp"] and self.chargefld != "pH":
                w(" %s", self.chargefld)
        for e in self.elt:
            w(" %s", e)
            if e == "Fe":
                w(" Fe(2) Fe(3)")
        w(" Tracer\n")
        w("    -water\n")
        w("    -equilibrium_phases CO2(g) O2(g) ")
        pt = ds.get("phasetable")
        for e in self.eqsolid:
            w(" %s", self.tv(pt, e, "PHREEQCphase"))
        w("\n")
        w("    -saturation_indices CO2(g) O2(g)")
        for p in self.kinphase[: self.nkinphase]:
            w(" %s", p)
        for e in self.eqsolid:
            w(" %s", e)
        w("\n")
        rates = self.scheme["rates"]
        if (self.nkinphase > 0 or stg["includeuptake"] or stg["includehumus"] or stg["includeresidues"]
                or stg["includerootresidues"] or stg["includemanure"] or stg["includeNkinetics"]):
            w("    -kinetic_reactants ")
            for p in self.kinphase[: self.nkinphase]:
                w(" %s", p)
            if stg["includetreatments"]:
                for t in self._treatments():
                    w(" %s", t)
            if ds["ETmethod"] == "KINETICS":
                w(" Evapotranspiration")
            if stg["includeuptake"]:
                for n in ds["nut"]:
                    w(" %suptake", n)
                w(" Ammuptake")
            if stg["includetreatments"] and stg["includeNkinetics"]:
                w(" Ureolysis")
            if stg["includehumus"]:
                for r in rates:
                    if "hum" in r:
                        w(" %s", r)
            if stg["includeresidues"] or stg["includerootresidues"]:
                for r in rates:
                    if "esidue" in r:
                        w(" %s", r)
            if stg["includeNkinetics"]:
                for r in rates:
                    if "itrification" in r:
                        w(" %s", r)
                for r in rates:
                    if "Nfix" in r:
                        w(" %s", r)
            if stg["includemanure"]:
                for r in rates:
                    if "anure" in r:
                        w(" %s", r)
            w("\n")
        for line in ("-alkalinity true", "-high_precision true", "-charge_balance true", "-percent_error true",
                     "-time false", "-step false", "-distance false"):
            w("    %s\n", line)

    def _treatments(self):
        t = self.ds.get("treatments", [])
        return [t] if isinstance(t, str) else list(t)

    def writepunch(self):
        stg = self.stg
        w = self._w
        npunch = len(PUNCHELT)
        w("\nUSER_PUNCH\n -headings")
        if stg["useexchange"]:
            w(" CEC")
        for e in PUNCHELT:
            w(" %s", e + "sys")
        w(" Ammsys Nsys Clsys")
        if stg["usesurface"]:
            if stg["includeWHAM"]:
                w(" WHAM")
            w(" HFO EDL EDLH2O")
            for e in PUNCHELT:
                w(" %s", e + "surf")
            w(" AmmHsurf Hsurf")
        jj = 0
        w("\n -start\n")
        if stg["useexchange"]:
            jj += 10
            w('%4d PUNCH TOTMOLE("X") # CEC (mol/cell) \n', jj)
        lines = [f'PUNCH SYS("{e}")' for e in PUNCHELT] + ['PUNCH SYS("Amm")', 'PUNCH SYS("N")', 'PUNCH SYS("Cl")']
        if stg["usesurface"] and stg["includeWHAM"]:
            lines += ["WHAM=0",
                      'WHAM=WHAM+TOTMOLE("H_a")+TOTMOLE("H_b")+TOTMOLE("H_c")+TOTMOLE("H_d")',
                      'WHAM=WHAM+TOTMOLE("H_e")+TOTMOLE("H_f")+TOTMOLE("H_g")+TOTMOLE("H_h")',
                      'WHAM=WHAM+TOTMOLE("H_ab")+TOTMOLE("H_ad")+TOTMOLE("H_af")+TOTMOLE("H_ah")',
                      'WHAM=WHAM+TOTMOLE("H_bc")+TOTMOLE("H_be")+TOTMOLE("H_bg")',
                      'WHAM=WHAM+TOTMOLE("H_cd")+TOTMOLE("H_cf")+TOTMOLE("H_ch")',
                      'WHAM=WHAM+TOTMOLE("H_de")+TOTMOLE("H_dg")',
                      "PUNCH WHAM # Total WHAM moles/cell",
                      'totedl=EDL("Charge","H")+EDL("Charge","Hfo") # Diffuse layer charge',
                      'edlh2o=EDL("water","H")+EDL("water","Hfo") # Diffuse layer water']
        elif stg["usesurface"]:
            lines += ['totedl=EDL("Charge","Hfo") # Diffuse layer charge',
                      'edlh2o=EDL("water","Hfo") # Diffuse layer water']
        if stg["usesurface"]:
            lines += ['HFO=TOTMOLE("Hfo_w")+TOTMOLE("Hfo_s") # HFO', "PUNCH HFO # Total HFO moles/cell",
                      "PUNCH totedl # Total EDL surface charge/cell", "PUNCH edlh2o # Total EDL water kg/cell"]
        for line in lines:
            jj += 1
            w("%4d %s\n", jj, line)
        if stg["usesurface"]:
            w("# Total of each desired element associated with SURFACES and in the whole system\n")
        else:
            w("# Total of each desired element in the whole system\n")
        if stg["usesurface"]:
            lines = []
            for e in PUNCHELT:
                if stg["includeWHAM"]:
                    lines.append(f'PUNCH SURF("{e}","H")+SURF("{e}","Hfo")+EDL("{e}","H")+EDL("{e}","Hfo")')
                else:
                    lines.append(f'PUNCH SURF("{e}","Hfo")+EDL("{e}")')
            if stg["includeWHAM"]:
                lines += ['PUNCH SURF("AmmH","H")+SURF("AmmH","Hfo")+EDL("Amm","H")+EDL("Amm","Hfo")',
                          'PUNCH SURF("H","H")+SURF("H","Hfo")+EDL("H","H")+EDL("Amm","Hfo")']
            else:
                lines += ['PUNCH SURF("AmmH","Hfo")+EDL("Amm")', 'PUNCH SURF("H","Hfo")+EDL("H")']
            for line in lines:
                jj += 1
                w("%4d %s\n", jj, line)
        assert npunch == 10
        w("\n -end\n")

    def writeinfluent(self, insoln, includetracer):
        ds, stg = self.ds, self.stg
        run, i = ds["run"], self.i
        w = self._w
        if "influent" in run.get("label", {}):
            w("SOLUTION %d # %s\n", insoln, run["label"]["influent"])
        else:
            w("SOLUTION %d \n", insoln)
        if "influent" in run.get("units", {}):
            w("    %10s %s\n", "units", run["units"]["influent"])
        else:
            w("# No units available for influent!! CHECK!")
        iwater = 1
        chgfld = self.charge0fld if insoln == 0 else self.chargefld
        ptarg = self.ptargH2O[i - 1]
        if stg["usewater"]:
            if insoln == 0 and ptarg > 0:
                iwater = ptarg
            elif insoln > 0:
                iwater = ds["cell"]["targH2O"][insoln - 1]
            w("    %10s %d # kg\n", "water", iwater)
        spas = ds["influent"].get("spas")
        for s in self.spi:
            if s == chgfld or s.lower() in ("cond", "svol", "ppt"):
                continue
            v = run["influent"][s]
            if isinstance(v, dict):
                continue
            v = v[i - 1]
            if v <= 0 or np.isnan(v):
                continue
            elti = ds["solute"][s]["elt"] if s in ds["solute"] else s
            if spas is not None:
                w("    %10s %g %s\n", elti, v, spas[elti])
            elif elti in ds["influent"].get("as", {}):
                w("    %10s %g as %s\n", elti, v, ds["influent"]["as"][elti])
            else:
                w("    %10s %g \n", elti, v)
        if chgfld.lower() == "ph" or iwater <= 0:
            w("    %10s %f charge\n", "pH", 5.5)
        else:
            if "pH" not in run["influent"]:
                w("    %10s %d CO2(g) %f\n", "pH", 5.5, run["si_CO2atm"][i - 1])
            w("    %10s %g charge\n", self.charge0fld, 0.000001)
        w("    %19s %f\n", "C(4) 10 CO2(g)", run["si_CO2atm"][i - 1])
        # MATLAB compares the influent struct with 'pe'/'O2g' (always false), so these are fixed
        if fixed(ds, "influent_redox"):
            w("    %10s %g O2(g) %f\n", "pe", ds["influent"].get("pe", 12.5), ds["influent"].get("O2g", -0.677778))
        else:
            w("    %10s %d O2(g) %f\n", "pe", 12.5, -0.677778)
        w("    %10s %f\n", "temp", run["airtemp"][i - 1])
        if includetracer:
            w("    %10s %g %s # Same concentration as Cl\n", "Tracer", ds["influent"]["msp"]["Cl"], "")
            if ptarg > 0:
                self.initw = 0
        if stg["usewater"]:
            if insoln == 0 and ptarg > 0:
                w("    %10s %d # kg\n", "water", ptarg)
            elif insoln > 0:
                w("    %10s %d # kg\n", "water", ds["cell"]["targH2O"][insoln - 1])
        if insoln > 0:
            return
        w("SAVE solution %d; END\n", insoln)
        w("\n")

    def writesoilsoln(self, insoln):
        ds, stg = self.ds, self.stg
        run, cell, i = ds["run"], ds["cell"], self.i
        w = self._w
        if "soilsoln" not in cell:
            self.writeinfluent(insoln, 0)
            return
        ss = cell["soilsoln"]
        if "msp" in ss.get("label", {}):
            w("SOLUTION %d # %s\n", insoln, ss["label"]["msp"])
        else:
            w("SOLUTION %d \n", insoln)
        if "units" in ss:
            w("    %10s %s\n", "units", ss["units"]["msp"])
        else:
            w("# No units available for solution %d!! CHECK!", insoln)
        ss.setdefault("chargefld", self.chargefld)
        msp = ss["msp"]
        if ss["chargefld"] not in msp:
            msp[ss["chargefld"]] = np.full(ds["ncells"], 1e-3)
        if stg["defineallsp"]:
            ds["parm"].setdefault("defineallsp", 1)
            for s in ds["solute"]:
                if s not in msp:
                    msp[s] = np.zeros(ds["ncells"]) + stg["minconc"]
                    ss.setdefault("spas", {})[s] = ""
        for s, vals in msp.items():
            if isinstance(vals, dict):
                continue
            if s not in ds["pt2"] and s.lower() != "ph":
                continue
            v = np.atleast_1d(vals)[insoln - 1]
            if np.isnan(v) or v <= 0:
                continue
            w("    %10s %g ", s, v)
            if s.lower() == ss["chargefld"].lower():
                w(" charge ")
            elif "spas" in ss:
                w("%s ", ss["spas"][s])
            elif s in ss.get("as", {}):
                w("as %s ", ss["as"][s])
            w("\n")
        if "si_CO2g" in run:
            siCO2g = run["si_CO2g"][insoln - 1, i - 1]
        else:
            siCO2g = run["si_CO2atm"][0]
        siO2g = run["si_O2g"][insoln - 1, i - 1] if "si_O2g" in run else -0.677778
        w("    %19s %f\n", "C(4) 10 CO2(g)", siCO2g)
        w("    %19s %f\n", "O(0) 10 O2(g)", siO2g)
        if "pH" not in msp:
            w("    %10s %d CO2(g) %f\n", "pH", 5.5, siCO2g)
        w("    %10s %d O2(g) %f\n", "pe", 12.5, siO2g)
        if "soiltemp" in run:
            w("    %10s %f\n", "temp", run["soiltemp"][insoln - 1, i - 1])
        else:
            w("    %10s %f\n", "temp", run["airtemp"][i - 1])
        if stg["usewater"]:
            w("    %10s %d # kg\n", "water", cell["targH2O"][insoln - 1])

    def percleachate(self):
        ds, out = self.ds, self.out
        run, i = ds["run"], self.i
        k = self.ileach
        w = self._w
        w("# Cell %d Leachate \n", k)
        leachrat = out["run"].setdefault("leachrat", np.zeros(ds["nrun"]))
        if run[self.outpercfld][ds["ncells"] - 1, i - 1] < self.vleachthresh:
            w("SOLUTION %d # Placeholder with pure water\n", k)
            w("    %10s %f\n", "temp", run["soiltemp"][self.ileach - 2, i - 1])
            w("SAVE solution %d; END\n", k)
            leachrat[i - 1] = 0
            return
        mixfrac = run["mixfrac"]
        if np.sum(mixfrac[k - 1, :, i - 1]) > 0:
            w("%s; ", ds["MIXture"][k - 1])
            for kk in range(1, int(ds["nmixsoln"]) + 1):
                if mixfrac[k - 1, kk - 1, i - 1] > 0:
                    w("%d %f; ", ds["MIXsoln"][kk - 1], mixfrac[k - 1, kk - 1, i - 1])
                    leachrat[i - 1] = mixfrac[k - 1, kk - 1, i - 1]
        w("SAVE solution %d; END\n", k)
        if run["airtemp"][i - 1] != run["soiltemp"][k - 2, i - 1]:
            w("SOLUTION_MODIFY %d;\n", k)
            w("    %10s %f\n", "-temp", run["airtemp"][i - 1])
            w("SAVE solution %d; END\n", k)

    def initsurface(self):
        ds, stg = self.ds, self.stg
        if not stg["usesurface"]:
            return
        if not stg["includeWHAM"]:
            self.initFeSurface(1)
            return
        k = self.k
        cell, cellH2O = ds["cell"], self.cellH2O
        surf = cell["surface"]
        w = self._w
        h = cellH2O[k - 1]
        w("SURFACE %d After Example 19 in PHREEQC manual v. 3\n", k)
        w("# WHAM model: -7.1 meq/gOrganicCarbon in total\n")
        w("# WHAM carboxylic monoprotic sites (-2.84 meq/gOC)\n")
        w("#   total moles of sites, m2surface/gOC, total grams OC\n")
        w("#   SOC active fraction: %g\n", ds["parm"]["SOCactfrac"])
        carb = surf["CARBOXYLICSITES"][k - 1] / h
        w("%s %g %f %f\n", self.carbsite[0], carb, stg["WHAMSA"],
          cell["gSOC"][k - 1] * ds["parm"]["SOCactfrac"] / h)
        for s in self.carbsite[1:]:
            w("%s %g; ", s, carb)
        w("\n")
        w("# WHAM phenolic monoprotic sites (-1.42 meq/gHumicAcid)\n")
        for s in self.phensite:
            w("%s %g; ", s, surf["PHENOLICSITES"][k - 1] / h)
        w("\n")
        w("# WHAM diprotic sites (-2.84 meq/gHumicAcid)\n")
        dip = self.dipsite
        for kk, s in enumerate(dip):
            w("%s %g; ", s, surf["DIPROTICSITES"][k - 1] / h)
            if kk == len(dip) - 1:
                w("\n")
                continue
            if s[:3] != dip[kk + 1][:3]:
                w("\n")
        w("-Donnan 1e-8\n")
        w("\n")
        pt = ds.get("phasetable")
        if stg["includedynamicsurfaces"] and pt is not None and "surface" in pt.columns:
            for line in pt.loc[pt["surface"] != "", "surface"]:
                w("    %s\n", line)
            if not (~pt["surface"].str.contains("Hfo_")).any():
                self.initFeSurface(0)
        else:
            self.initFeSurface(0)
        if stg["equilibratesurface"]:
            w("  -equilibrate %d\n", k)

    def initFeSurface(self, only):
        ds, stg = self.ds, self.stg
        if not stg["usesurface"]:
            return
        k, w = self.k, self._w
        pt = ds.get("phasetable")
        if pt is not None and "Surfname" in pt.columns:
            if (pt["Surfname"] == "Hfo").any() and stg["includedynamicsurfaces"]:
                return
        surf = ds["cell"]["surface"]
        h = self.cellH2O[k - 1]
        if surf["gHFO"][k - 1] <= 0 or ds["parm"]["Feactfrac"] <= 0:
            return
        if only > 0:
            w("SURFACE %d Iron oxides only\n", k)
        w("# HFO: assume extractable soil Fe is 10%% HFO with mass of goethite FeO(OH)\n")
        w("# HFO active fraction: %g\n", ds["parm"]["Feactfrac"])
        w("# HFO: Total moles (sites/cell), specific surface area m2/g, total HFO mass (g/cell)\n")
        if "HFO" not in surf:
            surf["HFO"] = {"m2g": np.full(ds["ncells"], 600.0)}
        if "eqphase" in surf["HFO"]:
            eq = surf["HFO"]["eqphase"][k - 1]
            w("Hfo_w %s equilibrium_phase %g %g; Hfo_s %s equilibrium_phase %g\n", eq,
              ds["parm"]["Feactfrac"] / h, surf["HFO"]["m2g"][k - 1], eq, ds["parm"]["Feactfrac"] / h)
        else:
            w("Hfo_w %g 600 %g; Hfo_s %g\n", surf["Hfo_w"][k - 1] / h, surf["gHFO"][k - 1] / h,
              surf["Hfo_s"][k - 1] / h)
        if only > 0:
            w("  -equilibrate %d\n", k)

    def initeqphases(self):
        ds, k, i = self.ds, self.k, self.i
        run, w = ds["run"], self._w
        w("EQUILIBRIUM_PHASES %d; ", k)
        w("    %10s %g %g; ", "CO2(g)", run["si_CO2g"][k - 1, i - 1], run["CO2gmol"][k - 1, i - 1])
        w("    %10s %g %g\n", "O2(g)", run["si_O2g"][k - 1, i - 1], run["O2gmol"][k - 1, i - 1])
        pt = ds.get("phasetable")
        if pt is not None:
            for e in self.eqsolid:
                w("    %10s 0 %g; ", self.tv(pt, e, "PHREEQCphase"), run["moles"][e][k - 1, i - 1] / self.cellH2O[k - 1])
            w("\n")

    def initsoilcells(self):
        if not self.init:
            return
        ds, stg, i = self.ds, self.stg, self.i
        cell, w = ds["cell"], self._w
        for k in range(1, ds["ncells"] + 1):
            self.k = k
            if self.init and not self.append and self.kinflag[k - 1, i - 1] > 1:
                self.kinflag[k - 1, i - 1] = 1
            w("# Cell %d depth range %5.1f - %5.1f cm; \n", k, cell["topdepth"][k - 1], cell["bottomdepth"][k - 1])
            self.writesoilsoln(k)
            self.initeqphases()
            if stg["useexchange"] and "eq" in self.EXCHmethod:
                units = ds["parm"].get("units", {})
                if "CEC" in units:
                    w("# CEC = %g %s applied to %g grams of soil\n", ds["parm"]["CEC"], units["CEC"],
                      cell["gsoil"][k - 1])
                w("EXCHANGE %d; # %s %s\n", k, self.out["metadata"]["runlab"], cell["depthlab"][k - 1])
                w("    %10s %g\n", "X", cell["X"][k - 1] / self.targH2O[k - 1, i - 1])
                w(" -equilibrate with solution %d\n", k)
                if stg["pitzer_exchange_gammas"]:
                    w(" -pitzer_exchange_gammas true\n")
            elif stg["useexchange"]:
                w("EXCHANGE %d; # %s %s\n", k, self.out["metadata"]["runlab"], cell["depthlab"][k - 1])
                for f, v in cell["exchange"].items():
                    if v[k - 1] <= 0:
                        continue
                    w("    %10s %g\n", f, v[k - 1] / self.targH2O[k - 1, i - 1])
                pt = ds.get("phasetable")
                if stg["includedynexchange"] and pt is not None and "exchange" in pt.columns:
                    for line in pt.loc[pt["exchange"] != "", "exchange"]:
                        w("    %s\n", line)
                w("MIX %d; %d %f;\n", k, k, 1)
            w("SAVE equilibrium_phases %d; ", k)
            if stg["useexchange"]:
                w("SAVE exchange %d; ", k)
            w("SAVE solution %d; END\n", k)
            if stg["usesurface"]:
                w("USE solution %d; USE equilibrium_phases %d; ", k, k)
                if stg["useexchange"]:
                    w("USE exchange %d; ", k)
                w("\n")
                self.initsurface()
                w("SAVE surface %d; ", k)
                w("SAVE equilibrium_phases %d; ", k)
                if stg["useexchange"]:
                    w("SAVE exchange %d; ", k)
                w("SAVE solution %d; END\n", k)
        w("SELECTED_OUTPUT\n    -selected_out true\n")

    def updatesoilcells(self):
        if self.init:
            return
        ds, stg, scheme, out, i = self.ds, self.stg, self.scheme, self.out, self.i
        run, w = ds["run"], self._w
        pt = ds.get("phasetable")
        for k in range(1, ds["ncells"] + 1):
            K, I = k - 1, i - 1
            w("# Cell %d %s\n", k, out["run"]["current"][i])
            w("SOLUTION_MODIFY %d;\n", k)
            w("    %10s %f\n", "-temp", run["soiltemp"][K, I])
            w("EQUILIBRIUM_PHASES_MODIFY %d\n", k)
            resetCO2g = resetO2g = 0
            if i > 1:
                if scheme["CO2gprotocol"].lower() == "reset":
                    CO2gmol, resetCO2g = run["CO2gmol"][K, I], 1
                elif scheme["CO2gprotocol"].lower() == "topup" and out["soil"]["CO2g"][K, I - 1] < run["CO2gmol"][K, I]:
                    CO2gmol, resetCO2g = out["soil"]["CO2g"][K, I - 1], 1
                if scheme["O2gprotocol"].lower() == "reset":
                    O2gmol, resetO2g = run["O2gmol"][K, I], 1
                elif scheme["O2gprotocol"].lower() == "topup" and out["soil"]["O2g"][K, I - 1] < run["O2gmol"][K, I]:
                    O2gmol, resetO2g = out["soil"]["O2g"][K, I - 1], 1
            if resetCO2g:
                w("-component CO2(g)\n")
                w("    %10s %g\n", "-si", run["si_CO2g"][K, I])
                w("    %10s %g\n", "-moles", CO2gmol)
            if resetO2g:
                w("-component O2(g)\n")
                w("    %10s %g\n", "-si", run["si_O2g"][K, I])
                w("    %10s %g\n", "-moles", O2gmol)
            if stg["usedust"] and self.neqsolids > 0:
                for e in self.eqsolid:
                    fs = self.tv(pt, e, "feedstock")
                    if fs not in run.get("treatment", {}) or run["treatment"][fs][K, I] == 0:
                        continue
                    w("-component %s\n", self.tv(pt, e, "PHREEQCphase"))
                    mol = run["moles"][e][K, I] / self.targH2O[K, I] + out["soil"][e][K, I - 1]
                    w("    %10s %g\n", "-moles", mol)
            if stg["useexchange"]:
                w("USE exchange %d;\n", k)
            if stg["usesurface"]:
                w("USE surface %d;\n", k)
            w("SAVE equilibrium_phases %d; ", k)
            if stg["useexchange"]:
                w("SAVE exchange %d; ", k)
            if stg["usesurface"]:
                w("SAVE surface %d; ", k)
            w("SAVE solution %d; END\n", k)

    def writekinetics(self):
        ds, stg = self.ds, self.stg
        if not stg["includekinetics"]:
            return
        k, i = self.k, self.i
        flag = self.kinflag[k - 1, i - 1]
        if flag == 0:
            return
        if flag == 2:
            self._w("USE kinetics %d;\n", k)
            return
        if flag == 3:
            self._modify_kinetics()
            return
        self._init_kinetics()

    def _modify_kinetics(self):
        """The KINETICS_MODIFY branch of writekinetics (kinflag 3)."""
        ds, stg, scheme, out = self.ds, self.stg, self.scheme, self.out
        run, cell = ds["run"], ds["cell"]
        k, i = self.k, self.i
        K, I = k - 1, i - 1
        targ = self.targH2O[K, I]
        soil = out.get("soil")
        nout = k if soil is not None else 0
        pt = ds.get("phasetable")
        w = self._w
        w("KINETICS_MODIFY %d %s\n", k, out["run"]["current"][i])
        for p in self.kinphase[: self.nkinphase]:
            w("    -component %s\n", p)
            mol = out["run"]["M0"][p][K, I] / targ
            w("    %14s %g\n", "-m", mol)
            w("    %14s %g\n", "-m0", mol)
            if "Implicit" in p:
                w("    %14s %g \n", "-d_params", ds["parm"]["impscale"])
            else:
                w("    %14s %g %f %f %f %f %f %f\n", "-d_params", out["run"]["RSA0"][p][K, I], 1.0, 1.0, 1.0, targ,
                  stg["skipSAinKINETICS"], self.tv(pt, p, "precipdissrat") * run["targH2O"][K, I] * self.tv(pt, p, "reversible"))
            w("    %10s %g \n", "-tol", 1e-7 if self.usemix else 1e-8)
        if ds["ETmethod"] == "KINETICS":
            w("    -component %s\n", "Evapotranspiration")
            w("    %14s %s 1\n", "-namecoef", "H2O")
            w("    %14s %g\n", "-d_params", run["cellET"][K, I] / targ)
        if stg["includetreatments"] and cell["tillagefrac"][K] > 0:
            tt = ds["timetable"]
            for tr in self._treatments():
                wtr = [j for j, t in enumerate(tt["treatment"])
                       if tt["runindex"][j] == i and t.lower() == tr.lower() and tt["molpm2"][j] > 0]
                if not wtr:
                    continue
                mol = tt["molpm2"][wtr[0]] / targ * cell["tillagefrac"][K]
                if soil is not None and "k_" + tr in soil:
                    mol = soil["k_" + tr][K, I - 1] + mol
                w("    -component %s %% treatment\n", tr)
                w("    %14s %g\n", "-m", mol)
                w("    %14s %g\n", "-m0", mol)
        demand = out["run"].setdefault("demand", {})
        nrun, ncells = ds["nrun"], ds["ncells"]
        MM = ds.get("MMuptake")
        if stg["includeuptake"]:
            Ndemand = run["cellcuptake"]["N"][K, I] / targ
            if nout > 0:
                Ndemand = Ndemand - soil["k_Ammuptake"][K, I - 1] - soil["k_Nuptake"][K, I - 1]
            if nout > 0 and "k_Nfix" in soil:
                Ndemand = Ndemand - soil["k_Nfix"][K, I - 1]
            wnut = ds["nut"].index("N")
            if ds["trnutflag"][wnut] > 0 and i > 1:
                Ndemand = Ndemand - soil["p_Nuptake"][K, I - 1]
                Ndemand = Ndemand - soil["p_Ammuptake"][K, I - 1]
            Ndemand = max(0, Ndemand)
            demand.setdefault("N", np.zeros((ncells, nrun)))[K, I] = Ndemand
            w("    -component %suptake\n", "Amm")
            w("    %14s %s 1\n", "-namecoef", "Amm")
            if self.planted[I] == 0:
                w("    %14s 0\n", "-m")
                w("    %14s 0\n", "-m0")
            ind = int(max(run["Nindex"][I], 1))
            w("    %14s %g %g %g\n", "-d_params", MM["Km"][ind - 1], MM["cmin"][ind - 1],
              Ndemand * stg["includeAmmuptake"])
            for kk, n in enumerate(ds["nut"]):
                if n.lower() == "n":
                    dem = Ndemand
                else:
                    dem = run["cellcuptake"][n][K, I] / targ
                    if nout > 0:
                        dem = dem - soil[f"k_{n}uptake"][K, I - 1]
                    if ds["trnutflag"][kk] > 0 and i > 1:
                        dem = dem - soil[self.trnutfld[kk]][K, I - 1]
                dem = max(0, dem)
                demand.setdefault(n, np.zeros((ncells, nrun)))[K, I] = dem
                w("    -component %suptake\n", n)
                w("    %14s %s 1\n", "-namecoef", n)
                if self.planted[I] == 0:
                    w("    %14s 0\n", "-m")
                    w("    %14s 0\n", "-m0")
                ind = int(max(run[n + "index"][I], 1))
                w("    %14s %g %g %g\n", "-d_params", MM["Km"][ind - 1], MM["cmin"][ind - 1], dem)
        if self.decompscheme.lower() == "swat":
            if stg["includeNkinetics"]:
                # MATLAB tests isfield(out,'k_Nhuma') (not out.soil), so pctorgC keeps its initial value
                if nout > 0 and "k_Nhuma" in (soil if fixed(ds, "pctorgc_update") else out):
                    scheme["pctorgC"][K, I] = (ds["pt2"]["C"]["molmass"] * cell["gSOC"][K] / cell["gSON"][K]
                                               * (soil["k_Nhums"][K, I - 1] + soil["k_Nhuma"][K, I - 1])
                                               / cell["gsoil"][K] * 100)
                w("    -component %s\n", "Denitrification")
                if self.schemeratesfile.lower() == "swatratesphases2":
                    w("    %14s %f %f %f %f\n", "-d_params", scheme["void"][K, I], 1.4,
                      cell["gSOC"][K] / cell["gsoil"][K] * 100)
                else:
                    w("    %14s %f %f %f %f\n", "-d_params", scheme["sut"][K, I], scheme["sdnco"], scheme["cdn"],
                      scheme["pctorgC"][K, I])
                w("    -component %s\n", "Nitrification")
                w("    %14s %f %f\n", "-d_params", scheme["etaSW"][K, I], scheme["etamidz"][K])
                if stg["includeuptake"] and self.withNfix:
                    w("    -component %s\n", "Nfix")
                    if self.planted[I] == 0:
                        w("    %14s 0\n", "-m")
                        w("    %14s 0\n", "-m0")
                    w("    %14s %f %f %f\n", "-d_params", run["legstage"][I], Ndemand,
                      run["targH2O"][K, I] / cell["FC"][K] / 0.85)
            if stg["includehumus"]:
                w("    -component %s\n", "Nhuma")
                w("    %14s %f %f\n", "-d_params", scheme["sut"][K, I], scheme["kNhumaday"])
            if (stg["includeresidues"] and cell["tillagefrac"][K] > 0) or stg["includerootresidues"]:
                mol = scheme["resCgfrac"] / 12.01 * (run["residue"][K, I] + run["rootresidue"][K, I]) / targ
                w("    -component %s\n", "Residue")
                if mol > 0:
                    if nout > 0:
                        mol = mol + soil["k_Residue"][K, I - 1]
                    w("    %10s %g\n", "-m", mol)
                w("    %10s %f %f %f %f\n", "-d_params", scheme["sut"][K, I], scheme["resCNgrat"], scheme["resCPgrat"],
                  scheme["kresday"])
            if stg["includemanure"] and cell["tillagefrac"][K] > 0:
                mol = scheme["resCgfrac"] / 12.01 * run["manure"][K, I] / targ
                w("    -component %s\n", "Manurefresh")
                if mol > 0:
                    if nout > 0:
                        mol = mol + soil["k_Manurefresh"][K, I - 1]
                    w("    %10s %g\n", "-m", mol)
                w("    %10s %f %f %f %f\n", "-d_params", scheme["sut"][K, I], scheme["resCNgrat"], scheme["resCPgrat"],
                  scheme["kmanday"])
        elif stg["includeNkinetics"]:
            w("Nitrification\n")
            w("    %10s %g \n", "-d_params ", run["nitrifscl"][K, I])
            w("    %10s %g \n", "-tol", stg["tol"])
            w("FreeBioNfix\n")
            w("    %10s %g %g \n", "-d_params ", cell["CFannET"][K], targ)
            w("    %10s %g \n", "-tol", stg["tol"])
        w("INCREMENTAL_REACTIONS true\n")
        w("USE kinetics %d;\n", k)

    def _init_kinetics(self):
        ds, stg, scheme, out = self.ds, self.stg, self.scheme, self.out
        run, cell = ds["run"], ds["cell"]
        k, i = self.k, self.i
        K, I = k - 1, i - 1
        targ = self.targH2O[K, I]
        tol = stg["tol"]
        pt = ds.get("phasetable")
        w = self._w
        w("KINETICS %d\n", k)
        for p in self.kinphase[: self.nkinphase]:
            w("%s\n", p)
            if stg["kinrelease"]:
                w("    %10s %s \n", "-formula", self.tv(pt, p, "formula"))
            else:
                w("    %10s %s 1\n", "-formula", p)
            mol = out["run"]["M0"][p][K, I] / targ
            w("    %10s %g\n", "-m", mol)
            w("    %10s %g\n", "-m0", mol)
            if "Implicit" in p:
                w("    %14s %g \n", "-parms", ds["parm"]["impscale"])
            else:
                w("    %10s %f %f %f %f %f %f %f\n", "-parms", out["run"]["RSA0"][p][K, I], 1.0, 1.0, 1.0, targ,
                  stg["skipSAinKINETICS"], self.tv(pt, p, "precipdissrat") * run["targH2O"][K, I] * self.tv(pt, p, "reversible"))
            w("    %10s %g \n", "-tol", tol)
        if ds["ETmethod"] == "KINETICS":
            w("%s\n", "Evapotranspiration")
            w("    %10s %s 1\n", "-formula", "H2O")
            w("    %10s 0\n", "-m")
            w("    %10s 0\n", "-m0")
            w("    %10s %g\n", "-parms", run["cellET"][K, I] / targ)
            w("    %10s %g \n", "-tol", tol)
        MM = ds.get("MMuptake")
        if stg["includeuptake"]:
            ind = 1
            for n in ds["nut"]:
                w("%suptake\n", n)
                if stg["nocounterion"]:
                    w("    %10s %s 1\n", "-formula", n)
                else:
                    w("    %10s %s\n", "-formula", ds["uptakeformula"][n])
                w("    %10s 0\n", "-m")
                w("    %10s 0\n", "-m0")
                ind = int(max(run[n + "index"][I], 1))
                w("    %10s %g %g %g\n", "-parms", MM["Km"][ind - 1], MM["cmin"][ind - 1],
                  run["cellcuptake"][n][K, I] / targ)
                w("    %10s %g \n", "-tol", tol)
            # Ammuptake uses the Michaelis-Menten index left over from the last nutrient
            if fixed(ds, "amm_uptake_mm_row"):
                ind = int(max(run["Nindex"][I], 1))
            Ammdemand = run["cellcuptake"]["N"][K, I] / targ * stg["includeAmmuptake"]
            w("Ammuptake\n")
            if stg["nocounterion"]:
                w("    %10s %s 1\n", "-formula", "AmmH+")
            else:
                w("    %10s %s 1 H+ -1\n", "-formula", "AmmH+")
            w("    %10s 0\n", "-m")
            w("    %10s 0\n", "-m0")
            w("    %10s %g %g %g\n", "-parms", MM["Km"][ind - 1], MM["cmin"][ind - 1], Ammdemand)
            w("    %10s %g \n", "-tol", tol)
        if self.decompscheme.lower() == "swat":
            if stg["includeNkinetics"]:
                w("Denitrification\n")
                w("    %10s %s \n", "-formula", "NO3- 1 Nzero -1")
                w("    %10s 0\n", "-m")
                w("    %10s 0\n", "-m0")
                if self.schemeratesfile.lower() == "swatratesphases2":
                    w("    %10s %f %f %f %f\n", "-parms", scheme["void"][K, I], 1.4,
                      cell["gSOC"][K] / cell["gsoil"][K] * 100)
                else:
                    w("    %14s %f %f %f %f\n", "-parms", scheme["sut"][K, I], scheme["sdnco"], scheme["cdn"],
                      scheme["pctorgC"][K, 0])
                w("    %10s %g \n", "-tol", tol)
                w("Nitrification\n")
                w("    %10s %s %f\n", "-formula", "AmmH+ 1 NO3- ", -scheme["Nnit"][K, I])
                w("    %10s 0\n", "-m")
                w("    %10s 0\n", "-m0")
                w("    %10s %f %f\n", "-parms", scheme["etaSW"][K, I], scheme["etamidz"][K])
                w("    %10s %g \n", "-tol", tol)
            if stg["includeNkinetics"] and stg["includeuptake"] and self.withNfix:
                w("Nfix\n")
                w("    %10s %f %s %f\n", "-formula H", -stg["Nfixprotons"], "Br", -stg["Nfixprotons"])
                w("    %10s 0\n", "-m")
                w("    %10s 0\n", "-m0")
                w("    %14s %f %f %f\n", "-parms", run["legstage"][I], run["cellcuptake"]["N"][K, I] / targ,
                  run["targH2O"][K, I] / cell["FC"][K] / 0.85)
                w("    %10s %g \n", "-tol", tol)
            if stg["includehumus"]:
                mol = cell["gSON"][K] / targ / 14.01
                w("Nhums\n")
                w("    %10s Nhuma \n", "-formula")
                w("    %10s %g\n", "-m", mol * (1 - scheme["humactfrac"]))
                w("    %10s %g\n", "-m0", mol * (1 - scheme["humactfrac"]))
                w("    %10s %g \n", "-tol", tol)
                w("    %10s %f %f\n", "-parms", scheme["kNhumsday"], scheme["humactfrac"])
                w("Nhuma\n")
                w("    %10s %s \n", "-formula", scheme["Nhumaformula"])
                w("    %10s %g\n", "-m", mol * scheme["humactfrac"])
                w("    %10s %g\n", "-m0", mol * scheme["humactfrac"])
                w("    %10s %f %f\n", "-parms", scheme["sut"][K, I], scheme["kNhumaday"])
                w("    %10s %g \n", "-tol", tol)
                w("Phum\n")
                if "PSmolrat" in scheme:
                    w("    %10s PO4-2 1 SO4-2 %g\n", "-formula", scheme["PSmolrat"])
                else:
                    w("    %10s PO4-2 1 \n", "-formula")
                w("    %10s %g\n", "-m", mol * 0.125)
                w("    %10s %g\n", "-m0", mol * 0.125)
                w("    %10s %g \n", "-tol", tol)
            tillflag = cell["tillagefrac"][K] > 0
            if stg["includerootresidues"] or (stg["includeresidues"] and tillflag):
                w("Residue\n")
                mol = (run["rootresidue"][K, I] + run["residue"][K, I]) * scheme["resCgfrac"] / 12.01 / targ
                w("    %10s %g\n", "-m", mol)
                w("    %10s %g\n", "-m0", mol)
                w("    %10s %f %f %f %f\n", "-parms", scheme["sut"][K, I], scheme["resCNgrat"], scheme["resCPgrat"],
                  scheme["kresday"])
                w("    %10s  %s \n", "-formula", scheme["residueformula"])
                w("    %10s %g \n", "-tol", tol)
            if stg["includemanure"] and tillflag:
                w("Manurefresh\n")
                mol = run["manure"][K, I] * scheme["resCgfrac"] / 12.01 / targ
                w("    %10s %g\n", "-m", mol)
                w("    %10s %g\n", "-m0", mol)
                w("    %10s %f %f %f %f\n", "-parms", scheme["sut"][K, I], scheme["resCNgrat"], scheme["resCPgrat"],
                  scheme["kmanday"])
                w("    %10s  %s \n", "-formula", scheme["manureformula"])
                w("    %10s %g \n", "-tol", tol)
            if stg["includetreatments"] and cell["tillagefrac"][K] > 0:
                tt = ds["timetable"]
                for tr in self._treatments():
                    wtr = [j for j, t in enumerate(tt["treatment"]) if t.lower() == tr.lower() and tt["molpm2"][j] > 0]
                    if not wtr:
                        continue
                    j = wtr[0]
                    mol = tt["molpm2"][j] / targ * cell["tillagefrac"][K] if tt["runindex"][j] == i else 0
                    w("%s\n", tr)
                    w("    %10s %g\n", "-m", mol)
                    w("    %10s %g\n", "-m0", mol)
                    w("    %10s  %s \n", "-formula", tt["formula"][j])
                    if "pellet" in tt["doseform"][j].lower():
                        P = scheme["Pellets"]
                        w("    %10s %f %f %f %f\n", "-parms", P["k1"], P["n1"], P["k2"], P["n2"])
                    w("    %10s %g \n", "-tol", tol)
        else:
            if stg["includeNkinetics"]:
                w("Nitrification\n")
                w("    %s \n", "-formula AmmH+ -1.0 NO3- 0.9994 None2O 0.0006")
                w("    %10s %g \n", "-parms ", run["nitrifscl"][K, I])
                w("    %10s %g \n", "-tol", tol)
                w("FreeBioNfix\n")
                w("    %s \n", "-formula AmmH+ 1.0")
                w("    %10s %g %g\n", "-parms ", cell["CFannET"][K], targ)
                w("    %10s %g \n", "-tol", tol)
        if stg["includetreatments"]:
            w("Ureolysis\n")
            w("    %s \n", "-formula Urea 1.0 H2O 2 AmmH+ -2.0 CO3-2 -1")
            w("    %14s %g\n", "-m", 0)
            w("    %14s %g\n", "-m0", 0)
            w("    %10s %g \n", "-parms ", stg["urease"])
            w("    %10s %g \n", "-tol", tol)
        md = self.out["metadata"].setdefault("PHREEQC", {})
        if stg["usecvode"]:
            w("-cvode true\n-cvode_order 5\n-cvode_steps %d\n", self.cvode_steps)
            md.update(solver="cvode", cvode_steps=self.cvode_steps)
        else:
            w("-runge_kutta 6\n-bad_step_max %d\n", self.badstepmax)
            md.update(solver="runge_kutta 6", bad_step_max=self.badstepmax)
        if self.usemix:
            if self.kinsteps is not None:
                ks = self.kinsteps if isinstance(self.kinsteps, str) else sprintf(" %g", self.kinsteps)
                w("-steps %s\n", ks)
                md["kinsteps"] = ks
            else:
                w("-steps %13.4g in %d steps\n", self.step[I] * 86400, self.nkinsteps)
                md.update(steps=self.step[I] * 86400, nkinsteps=self.nkinsteps)
            if self.stepdivide != 0:
                w("-step_divide %g\n", self.stepdivide)
                md["step_divide"] = self.stepdivide
        w("INCREMENTAL_REACTIONS true\n")
