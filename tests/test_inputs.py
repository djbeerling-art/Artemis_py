"""Input-reading functions against the authors' saved rBcomp ``ds`` struct."""

import os
from functools import lru_cache
from pathlib import Path

import pytest

from artemis_p import matref
from artemis_p.periodictable import get_periodic_table

RBCOMP = "rBcomp_Blueridge_Jan2015_Dec2020_03-Jun-2026"

pytestmark = pytest.mark.skipif(
    not (os.environ.get("ARTEMIS") and os.environ.get("ARTEMIS_PHREEQC")),
    reason="reference data not configured",
)


@lru_cache(maxsize=1)
def ref_ds():
    return matref.load_run(Path(os.environ["ARTEMIS_PHREEQC"]) / f"{RBCOMP}.mat")["ds"]


def test_periodic_table():
    pt, pt2 = get_periodic_table()
    assert matref.differences(pt, ref_ds()["pt"], "pt") == []
    assert pt2["Ca"]["molmass"] == 40.08  # value in periodictable.csv


def test_timetable():
    from artemis_p.timetable import read_timetable

    fname = Path(os.environ["ARTEMIS"]) / "examples" / "EF_data" / "EF_maize_B0_timetable.csv"
    ds = read_timetable({"parm": {}, "nrun": 2192}, fname)
    ref = ref_ds()["timetable"]
    # Later preprocessing steps add more columns to ds.timetable; check those readtimetable sets.
    # setupdatastruc/processtimetable later renames feedstock treatments to ds.parm.fsname.
    ours = {k: v for k, v in ds["timetable"].items() if k != "treatment"}
    assert matref.differences(ours, {k: ref[k] for k in ours}, "timetable") == []


def test_model_species():
    from artemis_p.species import set_model_species

    assert matref.differences(set_model_species(), ref_ds()["solute"], "solute") == []


def test_init_layer_data():
    from artemis_p.layers import read_init_layer_data

    fname = Path(os.environ["ARTEMIS"]) / "examples" / "EF_data" / "initlayers_Flanagan_EF_0-183cm.csv"
    ds = read_init_layer_data({"parm": {}}, fname)
    ref = ref_ds()["layer"]
    assert ds["nlayers"] == ref_ds()["nlayers"]
    # soilsoln is completed later from the initialChemistry file and wcell by
    # getcellsfromlayers; source holds the authors' relative paths.
    skip = {"soilsoln", "wcell", "source"}
    ours = {k: v for k, v in ds["layer"].items() if k not in skip}
    assert matref.differences(ours, {k: v for k, v in ref.items() if k not in skip}, "layer") == []
    assert set(ds["layer"]["source"]) == set(ref["source"]) - {"soilsoln"}


def test_phase_table():
    from artemis_p.phasetable import get_phase_table

    root = Path(os.environ["ARTEMIS"])
    files = [root / "examples" / "EF_data" / "Richland_loess_Batestown_till_phases_combined.csv"] + [
        root / "phreeqcdbincludefiles" / f
        for f in ("secondary_phases.csv", "feedstock_phases_pure.csv", "feedstock_phases_Lewis_THERMODDEM.csv")
    ]
    t = None
    for f in files:
        t = get_phase_table(f, t)
    ref = ref_ds()["phasetable"]
    rows = ref["_rownames"]
    # setupphases later keeps only the phases in use, adds columns and the
    # implicit phase, and marks Titanite/Hydroxyapatite implicit (parm=implicit).
    later = {"Implicit_blueridge"}
    implicit = {"Titanite_blueridge", "Hydroxyapatite_blueridge"}
    assert set(rows) - set(t.index) == later
    for c in t.columns:
        if c == "csvfile":  # authors' paths
            continue
        keep = [i for i, r in enumerate(rows) if r not in later and not (c == "PHREEQCphasetype" and r in implicit)]
        ours = [t.loc[rows[i], c] for i in keep]
        assert matref.differences(ours, [ref[c][i] for i in keep], c) == []
    assert (t.loc[list(implicit), "PHREEQCphasetype"] == "ignore").all()


def test_initial_chemistry():
    from artemis_p.chemistry import read_initial_chemistry
    from artemis_p.layers import read_init_layer_data

    data = Path(os.environ["ARTEMIS"]) / "examples" / "EF_data"
    ds = read_init_layer_data({"parm": {}}, data / "initlayers_Flanagan_EF_0-183cm.csv")
    ds = read_initial_chemistry(ds, data / "rainwater_Champaign_NTN_2021_initialChemistry.csv")
    ref = ref_ds()
    assert matref.differences(ds["init"], ref["init"], "init") == []
    assert matref.differences(ds["influent"], ref["influent"], "influent") == []
    assert matref.differences(ds["layer"]["soilsoln"], ref["layer"]["soilsoln"], "soilsoln") == []


def test_climate_forcing():
    from artemis_p.forcing import read_timeseries, set_run_time_arrays
    from artemis_p.layers import read_init_layer_data

    data = Path(os.environ["ARTEMIS"]) / "examples" / "EF_data"
    ds = read_init_layer_data({"parm": {}}, data / "initlayers_Flanagan_EF_0-183cm.csv")
    ds = set_run_time_arrays(ds, "2015-01-01", "2020-12-31")
    ref = ref_ds()
    run, rrun = ds["run"], ref["run"]
    assert ds["nrun"] == ref["nrun"]
    for k in ("year", "month", "day", "decyear", "doy", "exptday"):
        assert matref.differences(run[k], rrun[k], k) == []
    assert (run["tt"] == rrun["tt"]).all()
    for f in ("soiltemp", "airtemp", "soilresp", "CO2atm", "snow", "infil_ET"):
        ds = read_timeseries(ds, data / f"EF_CLM5_2015-2070_timeseries_{f}.csv")
    forced = ["soiltemp", "airtemp", "soilresp", "hetresp", "rootresp", "CO2atm",
              "snowfall", "snowmelt", "snowpack", "sublimation", "ET", "infil"]
    for k in forced:
        assert matref.differences(run[k], rrun[k], k, rtol=1e-12, atol=1e-12) == [], k
        assert run["units"][k] == rrun["units"][k]
        assert Path(ds["provided"][k]).name == Path(ref["provided"][k]).name


def _ef_phase_table():
    from artemis_p.phasetable import get_phase_table

    root = Path(os.environ["ARTEMIS"])
    files = [root / "examples" / "EF_data" / "Richland_loess_Batestown_till_phases_combined.csv"] + [
        root / "phreeqcdbincludefiles" / f
        for f in ("secondary_phases.csv", "feedstock_phases_pure.csv", "feedstock_phases_Lewis_THERMODDEM.csv")
    ]
    t = None
    for f in files:
        t = get_phase_table(f, t)
    return t


def test_cells_from_layers():
    import numpy as np

    from artemis_p.cells import get_cells_from_layers
    from artemis_p.chemistry import read_initial_chemistry
    from artemis_p.layers import read_init_layer_data

    data = Path(os.environ["ARTEMIS"]) / "examples" / "EF_data"
    ds = {"parm": {"nativeRSAscale": 1e-6, "implicit": ["Titanite_blueridge", "Hydroxyapatite_blueridge"]},
          "phasetable": _ef_phase_table()}
    ds = read_init_layer_data(ds, data / "initlayers_Flanagan_EF_0-183cm.csv")
    ds = read_initial_chemistry(ds, data / "rainwater_Champaign_NTN_2021_initialChemistry.csv")
    ds = get_cells_from_layers(ds)
    ref = ref_ds()
    cell, rcell = ds["cell"], ref["cell"]
    # Changed later: soilsoln (setupphases adds species), CEC (CECsorb converts
    # units), moles/RSA (moved into ds.run and zeroed), labels/units/source.
    later = {"soilsoln", "moles", "RSA", "CEC", "label", "units", "source"}
    for k in cell:
        if k not in later:
            assert matref.differences(cell[k], rcell[k], k, rtol=1e-10) == []
    for ph, v in cell["moles"].items():
        assert np.allclose(v, np.asarray(ref["run"]["moles"][ph])[:, 0], rtol=1e-10), ph
        assert np.allclose(cell["RSA"][ph], np.asarray(ref["run"]["RSA"][ph])[:, 0], rtol=1e-10), ph
    pt, rpt = ds["phasetable"], ref["phasetable"]
    for i, r in enumerate(rpt["_rownames"]):
        if r in pt.index:
            assert pt.loc[r, "PHREEQCphasetype"] == rpt["PHREEQCphasetype"][i], r


@lru_cache(maxsize=1)
def _ef_cells():
    """EF ds after getcellsfromlayers, built the way test_cells_from_layers does."""
    from artemis_p.cells import get_cells_from_layers
    from artemis_p.chemistry import read_initial_chemistry
    from artemis_p.layers import read_init_layer_data

    data = Path(os.environ["ARTEMIS"]) / "examples" / "EF_data"
    ds = {"parm": {"nativeRSAscale": 1e-6, "implicit": ["Titanite_blueridge", "Hydroxyapatite_blueridge"],
                   "tillagedepth": 18, "Feactfrac": 0.01, "SOCactfrac": 0.629, "CECopt": "CECcarboxylic",
                   "setexaciditysp": "Al"},
          "phasetable": _ef_phase_table(), "EXCHmethod": "comp", "scale": {"pCO2g": 1.0}}
    ds = read_init_layer_data(ds, data / "initlayers_Flanagan_EF_0-183cm.csv")
    ds = read_initial_chemistry(ds, data / "rainwater_Champaign_NTN_2021_initialChemistry.csv")
    return get_cells_from_layers(ds)


def test_cec_sorb():
    import copy

    from artemis_p.sorption import cec_sorb, set_soil_texture

    ds = set_soil_texture(copy.deepcopy(_ef_cells()))
    ds = cec_sorb(ds)
    cell, rcell = ds["cell"], ref_ds()["cell"]
    for k in ("tillagefrac", "X", "basesatsum", "cecsum", "CECcmolckgsoil", "eqCEC", "SOCeq", "XCECfrac",
              "exchange", "CEC"):
        assert matref.differences(cell[k], rcell[k], k, rtol=1e-12) == [], k
    for k, v in cell["surface"].items():
        assert matref.differences(v, rcell["surface"][k], k, rtol=1e-12) == [], k
    for k in ("X", "CEC", "exchange", "basesatsum", "cecsum", "eqCEC", "SOCeq"):
        assert cell["units"][k] == rcell["units"][k], k
    assert cell["units"]["surface"] == {k: v for k, v in rcell["units"]["surface"].items() if k in cell["units"]["surface"]}
    wham = ref_ds()["parm"]["WHAM"]
    assert ds["parm"]["WHAM"] == {k: (int(v) if isinstance(v, float) else v) for k, v in wham.items()}


def test_process_timetable(tmp_path):
    from artemis_p.forcing import set_run_time_arrays
    from artemis_p.timetable import read_timetable
    from artemis_p.treatments import process_timetable, write_treatment_rates

    root = Path(os.environ["ARTEMIS"])
    ds = {"parm": {"fsname": "Blueridge", "writetrRATESfiles": 0}, "dbdir": root / "phreeqcdbincludefiles"}
    ds = set_run_time_arrays(ds, "2015-01-01", "2020-12-31")
    ds = read_timetable(ds, root / "examples" / "EF_data" / "EF_maize_B0_timetable.csv")
    ds = process_timetable(ds)
    ref = ref_ds()
    tt, rtt = ds["timetable"], ref["timetable"]
    for k in ("runindex", "molpm2", "molmass", "formula", "treatment", "dosage"):
        assert matref.differences(tt[k], rtt[k], k, rtol=1e-12) == [], k
    assert ds["treatments"] == [ref["treatments"]]
    assert ds["INCLUDE"] == [ref["INCLUDE"]]
    rft = ref["ferttable"]
    for k in ("treatment", "PHREEQCformula", "molmass", "ttindex"):
        assert matref.differences(ds["ferttable"][k], rft[k], k, rtol=1e-12) == [], k
    # The archived examples/EF2/Pellets.RATES was written for another run
    # (CaCO3_trmt), so check against the copy PHREEQC echoed in the rBcomp .out.
    (path,) = write_treatment_rates(ds, tmp_path)

    def norm(lines):
        return [x for x in (ln.split("#")[0].strip().lstrip("-") for ln in lines) if x and x != "RATES"]

    mine = norm(path.read_text().splitlines())
    out = (Path(os.environ["ARTEMIS_PHREEQC"]) / f"{RBCOMP}.out").read_text(encoding="latin-1").splitlines()
    i = next(k for k, ln in enumerate(out) if "UAN28 # Slow-release" in ln)
    assert mine == norm(out[i:i + 200])[:len(mine)]


def test_setup_phases():
    import copy

    import numpy as np
    import pandas as pd

    from artemis_p.forcing import set_run_time_arrays
    from artemis_p.phases import setup_phases
    from artemis_p.sorption import cec_sorb, set_soil_texture
    from artemis_p.timetable import read_timetable
    from artemis_p.treatments import process_timetable

    root = Path(os.environ["ARTEMIS"])
    ds = copy.deepcopy(_ef_cells())
    ds["parm"].update(fsname="Blueridge", writetrRATESfiles=0, RSAscale=2)
    ds["dbdir"] = root / "phreeqcdbincludefiles"
    ds["feedstockPSDfile"] = "../EF_data/Lewis_feedstocks_PSD.csv"
    ds = cec_sorb(set_soil_texture(ds))
    ds = set_run_time_arrays(ds, "2015-01-01", "2020-12-31")
    ds = read_timetable(ds, root / "examples" / "EF_data" / "EF_maize_B0_timetable.csv")
    ds = setup_phases(process_timetable(ds), basedir=root / "examples" / "EF2")
    ref = ref_ds()

    def col(df, c):
        return df[c].to_numpy() if pd.api.types.is_numeric_dtype(df[c]) else list(df[c])

    t, rt = ds["phasetable"], ref["phasetable"]
    assert list(t.index) == rt["_rownames"]
    assert [c for c in t] == [c for c in rt if c != "_rownames"]
    for c in t:
        if c == "csvfile":  # MATLAB keeps the relative path it was given
            assert [Path(p).name for p in t[c]] == [Path(p).name for p in rt[c]]
            continue
        assert matref.differences(col(t, c), rt[c], c, rtol=1e-10) == [], c
    ft, rft = ds["feedstocktable"], ref["feedstocktable"]
    assert [c for c in ft] == [c for c in rft if c != "_rownames"]
    for c in ft:
        assert matref.differences(col(ft, c), rft[c], c, rtol=1e-10) == [], c
    run, rrun = ds["run"], ref["run"]
    for k in ("moles", "RSA", "treatment", "grams"):
        assert sorted(run[k]) == sorted(rrun[k]), k
        for p, v in rrun[k].items():
            assert np.allclose(run[k][p], v, rtol=1e-10, atol=0), (k, p)
    for k in ("nucSA", "nucN", "treatment"):
        assert matref.differences(ds["cell"][k], ref["cell"][k], k, rtol=1e-12) == [], k
    assert matref.differences(ds["stoich"], ref["stoich"], "stoich", rtol=1e-12) == []
    assert matref.differences(ds["timetable"]["dosekg"], ref["timetable"]["dosekg"], "dosekg") == []


@lru_cache(maxsize=1)
def _rbcomp_preprocessed():
    from artemis_p.preprocess import setup_data_struc
    from artemis_p.reference import EF2_FILES, ef2_resolver, rbcomp_settings

    resolve, ef2 = ef2_resolver()
    return setup_data_struc(rbcomp_settings(), EF2_FILES, resolve, basedir=ef2)


def test_preprocess_rbcomp():
    """setupdatastruc up to setrunarrays, end to end from the EF2 inputs."""
    ds, ref = _rbcomp_preprocessed(), ref_ds()
    assert ds["inputfiles"] == [f for f in ref["inputfiles"] if "obs" not in f]
    assert matref.differences(ds["provided"], ref["provided"], "provided") == []
    run, rrun = ds["run"], ref["run"]
    # residue, rootresidue and mixfrac are filled in later by setuprun/setupnutuptake
    for k in ("crop", "plant", "wtt", "sowdepth", "harvest", "peaktot", "peakroot", "treated",
              "tillagefrac", "si_CO2atm", "precip", "airtemp", "soiltemp", "ET", "infil",
              "moles", "RSA", "treatment", "grams", "gdust"):
        assert matref.differences(run[k], rrun[k], k, rtol=1e-12) == [], k
    for k in ("nmix", "nmixsoln", "MIXture", "MIXsoln", "MIXsolnstr", "plantdensity", "mark"):
        assert matref.differences(ds[k], ref[k], k) == [], k
    # 1x1 MATLAB cell arrays come back from the .mat file as plain strings
    assert ds["treatments"] == [ref["treatments"]]
    assert ds["INCLUDE"] == [ref["INCLUDE"]]
    tt, rtt = ds["timetable"], ref["timetable"]
    assert sorted(tt) == sorted(rtt)
    for k in tt:
        assert matref.differences(tt[k], rtt[k], k) == [], k
    # Set later: evapfrac, pO2g, si_O2g (setuprun) and exchfld (getphreeqcrun)
    later = {"evapfrac", "pO2g", "si_O2g", "exchfld"}
    assert sorted(ds["cell"]) == sorted(k for k in ref["cell"] if k not in later)
