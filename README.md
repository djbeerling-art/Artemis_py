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

## Status

- [x] PHREEQC runner, selected-output reader, MATLAB `.mat` reference loader
- [x] rBcomp day 1 reproduced from the published input (PHREEQC 3.7.3)
- [ ] Input reading
- [ ] Soil, phases and chemistry setup
- [ ] Climate and hydrology
- [ ] PHREEQC input writer (match published `.in_init` files)
- [ ] Daily driver loop
