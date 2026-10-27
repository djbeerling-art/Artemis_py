# Artemis_p

Python port of **ARTEMIS v1.1**, the reactive transport enhanced-weathering
model for soils by Lyla Taylor (Zenodo
[10.5281/zenodo.21371426](https://doi.org/10.5281/zenodo.21371426), CC-BY-4.0;
GMD discussion paper [egusphere-2025-5823](https://doi.org/10.5194/egusphere-2025-5823)).

ARTEMIS is MATLAB code that drives the USGS PHREEQC geochemical model, one
PHREEQC run per simulated day. This port reproduces it in Python, checking each
piece against the published reference files. See
[docs/PORTING_PLAN.md](docs/PORTING_PLAN.md) for the stages and status.

## Setup

```bash
pip install -e .[test]
scripts/build_phreeqc.sh v3.7.3          # the PHREEQC version the published outputs used
export PHREEQC=$PWD/phreeqc-v3.7.3/build/phreeqc
```

Reference data come from the two Zenodo archives:

```bash
export ARTEMIS=/path/to/soilmodel                       # ARTEMIS_V1p1_code tarball
export ARTEMIS_PHREEQC=/path/to/PHREEQCfiles_GMD_R1     # PHREEQC input/output tarball
# missing from v1.1; take it from the v1.0 archive (Zenodo 10.5281/zenodo.17682252)
export ARTEMIS_FLANAGAN=/path/to/examples/EF_CLM5/Flanagan_BlueRidge_PHASESRATES_withimplicitphases
pytest
```

## Bug fixes

By default the port keeps the MATLAB code's bugs, so its output matches the
published runs. Set `ds["fixes"] = "all"` in the settings to fix the clear
coding errors; see the Fixes section of
[docs/PORTING_PLAN.md](docs/PORTING_PLAN.md) and `artemis_p/fixes.py`.

## Speed

A day takes about 1.2 s: about 0.8 s of PHREEQC kinetics (one simulation
per soil cell), 0.1 s of mixing and 0.15 s of Python. Loading PHREEQC and the
database takes only 0.02 s. Options for `PhreeqcRun`:

- `cellworkers=N` runs the cells' kinetics in N parallel PHREEQC processes
  after the day's mixing (`phreeqc.run_cells_parallel`). It agrees with a
  serial run to about 1e-5, except in cells with almost no water (rBcomp
  cell 2, frozen, 0.01 mm), where it is about 1e-4.
- `printoutput=0` writes only warnings and errors to the `.out` file (about
  0.1 s a day).
- `engine="iphreeqc"` runs PHREEQC in memory through the IPhreeqc library
  (`artemis_p.iphreeqc`) instead of starting the executable each day. The
  state stays in memory rather than going through the daily dump file, and
  the PHASES/RATES includes and the SELECTED_OUTPUT/USER_PUNCH definitions
  are read again only when they change. With `cellworkers=N` the cells'
  kinetics run in N threads, each cell in its own IPhreeqc instance, so the
  results are the same for any N. The cell with the least water runs where
  the state already is: started from a dump, it can take twice as long to
  converge. Results agree with the executable to about 3e-4 of each field's
  range over 60 days (the executable rounds the state to text every day).
  The whole state is dumped at the end of the run, for restarts.

On 4 cores, 60 days of rBcomp take 1.14 s a day with the executable, 0.94
with `engine="iphreeqc"`, and 0.56 with `engine="iphreeqc", cellworkers=4,
printoutput=0`.

Building IPhreeqc: there is no 3.7.3 release of IPhreeqc, so use v3.7.4 of
github.com/usgs-coupled/IPhreeqc with `src/phreeqcpp/pitzer.cpp` and
`transport.cpp` copied from PHREEQC v3.7.3 (the only engine files that
differ), then `cmake -DBUILD_SHARED_LIBS=ON -DBUILD_TESTING=OFF` and
`make IPhreeqc`. Point `$IPHREEQC_LIB` at the `libIPhreeqc.so` it builds.

Separate sites or scenarios can also run in parallel.

## Status

- [x] PHREEQC runner, selected-output reader, MATLAB `.mat` reference loader
- [x] rBcomp day 1 reproduced from the published input (PHREEQC 3.7.3)
- [x] Input reading, soil, phases and chemistry setup (rBcomp path)
- [x] Climate and Richards hydrology (rBcomp path)
- [x] PHREEQC input writer and daily driver loop: the full 2,192-day rBcomp
  run matches the published one (last-day `.in` byte-identical with the
  executable; end-of-run totals within 1e-7 with `engine="iphreeqc"`)
- [x] End-of-run totals and comparison with observations
- [x] Bug fixes as switches (`artemis_p/fixes.py`)
- [x] Day-1 input of the other published EF2 scenarios: 36 of the 37 run
  scripts that have a published `.in_init` give it byte for byte
  (`artemis_p/ef2.py`, `scripts/check_day1.py`, `tests/test_scenarios.py`),
  including the controls, forsterite and wollastonite runs, the 2015-2070 long
  runs, the spin-ups and the runs restarted from them. Only rBcomp is checked
  end to end.
- [x] Restarts from a spin-up (`noappend`: rspinstdf, rwms, rwss, rwssf). The
  authors' spin-up states were not published, so `scripts/run_scenario.py`
  takes `noappend=PATH` to restart from the port's own spin-up dump.
- [x] Layer fluxes and net CDR (`artemis_p/cdr.py`)
- [ ] Not ported: bucket hydrology (rBbucket; not needed), bottom-up mixing,
  SWAT snow, potential ET, soil respiration, drains, plotting
