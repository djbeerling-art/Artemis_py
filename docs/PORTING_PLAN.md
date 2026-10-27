# Porting plan

The goal is a Python version of the ARTEMIS v1.1 MATLAB code (Taylor, 2026,
Zenodo 10.5281/zenodo.21371426, CC-BY-4.0). It should produce the same PHREEQC
inputs and the same daily results as the MATLAB code, so that it can be
debugged and extended without MATLAB.

## How ARTEMIS runs

Each example script (e.g. `examples/EF2/rBcompscript.m`) does:

1. `EFrunsettings2(...)`: turns run options into argument lists.
2. `ds = setupdatastruc(...)`: the preprocessor. It reads the input CSVs, sets
   up soil layers and cells, mineral phases, sorption, fertiliser, plant
   uptake, climate forcing and hydrology (water mixing fractions between cells).
3. `r = getphreeqcrun(ds, ...)`: the driver. For each day it writes a PHREEQC
   `.in` file, runs `phreeqc`, reads the selected output and updates the state
   (reactive surface areas, gas moles, nutrient uptake) for the next day.
   rBcomp runs 2192 daily iterations.

## What we can check against

| Reference | Source | Checks |
|---|---|---|
| `<run>.in_init` (~50 runs) | PHREEQCfiles_GMD_R1 | Day-1 PHREEQC input written by the preprocessor and driver, across many option combinations |
| `rBcomp...mat` | PHREEQCfiles_GMD_R1 | Full `ds` struct (preprocessor output) and daily results for all 2192 days |
| `rBcomp...in`, `.sel`, `.out` | PHREEQCfiles_GMD_R1 | Last-day input and output, after the whole driver loop |
| `rBcomp...sel_init` | PHREEQCfiles_GMD_R1 | Day-1 PHREEQC output. Reproduced with PHREEQC 3.7.3 (pH within 2e-6) |

## Stages

Line counts are from the MATLAB call graph (`tools/matlab_callgraph.py`).
Plotting code (`getfig`, `show*`, figure scripts) is not ported.

1. **Infrastructure (done).** PHREEQC runner, `.sel` reader and comparison,
   `.mat` loader, rBcomp day-1 regression test.
2. **Input reading (~2,000 lines; in progress: `readcsvdata`, `getperiodictable`, `readtimetable`, `getdecyear`, `getscreenedtable`, `setmodelspecies`, `readinitlayerdata`, `getphasetable`, `checkphasetype` (table defaults), initial chemistry, run calendar and climate time series done; `readobsdata` (validation data only) still to do).** `readinputfiles`, `readcsvdata`,
   `readtimetable`, `readinitlayerdata`, `readobsdata`, `getscreenedtable`,
   `getphasetable`, `getperiodictable`, `setmodelspecies`, `getdecyear`.
   Checked against `ds.timetable`, `ds.layer`, `ds.phasetable` in the `.mat`.
3. **Soil, phases and chemistry setup (~4,000 lines; in progress: `getcellsfromlayers`, `setupsoilphases`, `pedotransfer`, full `checkphasetype`, `setsoiltexture`, `CECsorb`, `processtimetable`, `getfert`, `setupphases` (vintage PSD path), `setnucSA`, full `implicitcalc`, `parse_formula`, `setfeedstock`, `setrunarrays`, `setcropfromtimetable` done; `preprocess.setup_data_struc` runs the whole chain from the EF2 input files and matches the saved ds up to this point).** `setupphases`,
   `setupsoilphases`, `checkphasetype`, `implicitcalc`, `getcellsfromlayers`,
   `pedotransfer`, `markPSDinit`, `getRosinRammlerPSD`, `setnucSA`,
   `setstoichiometry`, `CECsorb`, `Cerling1991`, `getfert`, `setupnutuptake`,
   `writePHASESRATES`. Checked against the matching `ds` fields.
4. **Climate and hydrology (~6,000 lines; done for the rBcomp path:
   `daySWAT7`, `setuprun`, `getrichquick`/`soilwatermovement` (CLM5
   zengdecker_2009 Richards solver, phase change, water table, drainage),
   `setmixfracs2` (top-down), `getsoilCO2g` (Cerling 1991 time series) and
   `setupnutuptake`; all 2,192 days match the saved ds to ~1e-11. Not ported:
   `bucket`, bottom-up mixing, `SWATsnow`, `potET`, `getsoilresp`, drains).** `setuprun`, `daySWAT7`,
   `SWATsnow`, `potET`, `getsoilCO2g`, `getsoilresp`, `setmixfracs2`,
   `bucket`, `getrichquick`, `soilwatermovement`, `drainrate`. Checked against
   `ds.run` (mixing fractions, temperatures, CO2).
5. **PHREEQC input writer (part of `getphreeqcrun`, 3,846 lines).** Day-1
   inputs must match every published `.in_init`, apart from paths and
   timestamps.
6. **Driver loop (~2,500 lines).** The rest of `getphreeqcrun`, plus
   `getselectedoutput`, `reshapeselectedoutput`, `updateRSA`, `gasmoles`,
   `assessrun2`, `gettotalfields`, `assessuptake`, `etaRH`, `getbasesatpct`,
   `checkuselog`. Checked against the daily series in the `.mat` and the
   final-day `.in`.

## Known issues in the published v1.1 archive

- `examples/EF_CLM5/Flanagan_BlueRidge_PHASESRATES_withimplicitphases` is
  included by the EF2 runs but missing from v1.1. The v1.0 copy gives the
  same echoed input.
- The published outputs were made with PHREEQC 3.7.3. PHREEQC 3.8.6 and
  3.9.0 differ by up to 0.008 pH from the first kinetic step onward.
- `examples/EF2/Pellets.RATES` defines the rate law for `CaCO3_trmt`, not
  `UAN28`, so it is left over from another run. The published rBcomp `.out`
  echoes a UAN28 version, which the port rebuilds from `Pellets.dat`
  (`treatments.treatment_rates_text`) and matches.
- **The published runs predate the v1.1 hydrology code.** The rBcomp `.mat`
  is dated 3 June 2026, but `getrichquick.m` and `setuprun.m` were changed on
  20 and 21 June 2026. The saved ds has no `run.sat` or
  `run.qflx_rootsoi_day`, has `run.trueairpor`, and keeps
  `parm.draindepthm = 1.3` (v1.1 sets it to -1). Most importantly, v1.1
  getrichquick sets `run.airpor = por - Sws`, so getsoilCO2g no longer uses
  its own `max(por - fc, por - Sws)`. With the v1.1 code, soil pCO2 for
  rBcomp is a median 1.8 to 5.4 times higher per cell than in the published
  run (up to 12 times; mean 0.14 atm instead of 0.025 to 0.034 atm in the
  30 to 124 cm cells). `setup_run(..., legacy_airpor=True)` reproduces the
  published values; the default follows v1.1.

## Bugs found in the MATLAB code

### Likely to affect published results

- **Native mineral surface area ignores mineral abundance (`setupsoilphases.m`).**
  `RSAD = PMD./psdthresh./1000./G(m).*6` uses the mass of *all* soil
  minerals in each size class (`PMD`) and is never multiplied by the
  mineral's proportion `prop.(phase)`. As a result every native mineral gets
  the reactive surface area of the whole mineral fraction, differing only by
  density. In the published rBcomp day-1 input, cell 8 has 0.29 mol of
  MagnetiteRLBT with RSA 0.072 m2 and 4055 mol of QuartzRLBT with RSA
  0.142 m2. The saved `ds.run.RSA` shows identical `RSA x G` for every
  native phase, including CalciteBT in layers with zero calcite. Trace
  minerals (magnetite, ilmenite, garnet, epidote, hornblende, pyroxenes) are
  therefore given far too much surface relative to their mass. The global
  `nativeRSAscale = 1e-6` may have been tuned around this. The fix is
  presumably `mass`-weighting, `RSAD.*prop.(phase{m})`.
- **Hornblende formula typo (`examples/EF_data/Richland_loess_Batestown_till_phases_combined.csv`).**
  HornblendeRLBT's `Formula` is `Ca2Mg52Si8O22(OH)2`; it should be
  `Ca2Mg5Si8O22(OH)2` (tremolite, as in its reaction and molar mass). The
  formula goes into the KINETICS `-formula` line, so each mole of native
  hornblende that dissolves releases 52 Mg instead of 5.

### Other issues

These are ported as written, so outputs keep matching the reference. They are
flagged here to fix once the port is verified.

- **Leap years (`getdecyear.m`).** It tests "divisible by 4 and not divisible
  by 400", which makes 2000 a non-leap year (and 1900 and 2100 leap years).
  This doesn't affect 2001-2099, but the spin-up runs (`rspin*`) start in
  January 2000. Python: `timetable.is_leap_matlab`.

- **NaN filling (`readinitlayerdata.m` `fillNaN`).** `for inan=wnan` loops
  over a column vector, so it runs once with the whole vector. Gaps fill
  correctly only when a parameter has exactly one missing layer. The EF input
  has no gaps, so the published runs are unaffected. Python fills each gap
  (`layers.fill_nan`).

- **Density units (`setupphases.m` line 341).** When the feedstock density is
  in g/cm3, it sets `rho_kgm3 = G./1000`. It should be `G.*1000`, so the value
  is a factor of 10^6 too small (e.g. 0.0052 instead of 5206 for magnetite in
  the saved rBcomp phase table). `rho_kgm3` is only read in `markPSDinit.m`
  line 175, on the alternative ("Mark") surface-area path. That line also
  indexes the table with an unquoted `rho_kgm3`. The default path
  (`SAmethod='traditional'`, used by the published runs) computes
  `G.*1000` correctly at line 225, so the published runs look unaffected.
  To be confirmed when `markPSDinit` is ported.

- **Measured soil pH discarded (`readinputfiles.m` `readinitialchemistry`).**
  `readinitlayerdata` stores the initlayers pH profile (6.27 to 7.9 for EF) in
  `ds.layer.soilsoln.msp.pH`. `readinitialchemistry` then replaces the whole
  `ds.layer.soilsoln` struct with the rainwater composition. With no
  porewater column it also deletes the rainwater pH. The soil cells therefore
  start from rainwater equilibrated with soil CO2 (`pH 5.5 CO2(g) ...` in the
  published `.in_init`), and the measured profile is never used. This may be
  intended, but the code's comments suggest not. Worth asking the authors.
- **Depth units in time series (`checktsdepths`).** The test
  `numel(contains(ucell,'cm'))` is always non-zero, so depths given in metres
  would never be converted to cm. All EF files use cm.
- **Stale index (`settimeseries`).** For depth-resolved files with per-month
  units it reads `ds.run.units.(var{i})` with an `i` left over from an
  enclosing loop. No published input takes this path.

- **Vector `if` (`pedotransfer.m`).** `if (om_frac>pcalpha)` and
  `if (om_frac<1)` test vectors, so they are true only when *every* cell
  qualifies. A profile mixing organic and mineral layers would get no
  percolation term at all. The EF soils are all below the threshold, so they
  are unaffected.
- **Label slip (`getcellsfromlayers.m`).** It sets `units.thickness='fraction'`
  where `units.frac` was meant. Cosmetic only.
- **`eqCEC` is 100 times too large (`CECsorb.m` line 391).**
  `eqCEC=cecsum.*gsoil/1000` uses `cecsum` while it is still in meq/100 g, so
  the result is 100x the eq/cell in its label (EF top cell: 1863.6 instead of
  18.6). It feeds only the `ccateqCECrat` diagnostic in `gettotalfields.m`,
  so the PHREEQC run is unaffected but that reported ratio is 100x too small.
- **HFO iron from FeOOH phases (`CECsorb.m`).** Iron for the HFO surface is
  taken as `moles/100*55.85/molmass`, which treats a mole amount as a soil mass
  percent. EF has no native FeOOH phase, so it falls back to 2.79 mg Fe/g soil
  and is unaffected.
- **Stale treatment RATES files (`setupdatastruc.m` `processtimetable`).**
  With the default `writetrRATESfiles=0`, an existing `Pellets.RATES` in the
  run directory is reused without checking which treatments it defines. A
  file left by an earlier run with other fertilisers is silently used, and
  PHREEQC then has no rate law for the new treatment (the archived EF2 copy
  is an example of such a leftover).
- **Geometric area from lambda (`setupphases.m` `setRSA`).** With a BET area
  and lambda but no geometric area it sets `geo_m2g=lambda.*BET_m2g`. Since
  lambda is BET/geo, this should be `BET_m2g./lambda`. EF gives all three, so
  it is unaffected.
- **Per-feedstock RSAscale never used (`setupphases.m` `setRSA`).** The test
  `numel(strcmpi(fld,'RSAscale'))==0` is never true, and `ds.parm.RSAscale`
  is always set (default 1) at the top of the function, so an `RSAscale`
  column in a feedstock table is ignored.
- **Timetabled phase treatments (`setupphases.m` `checktimetablephase`).**
  `ds.timetable.RSA=zeros(numel(runindex))` makes an n-by-n matrix, and an
  assumed RSA is written to every phase entry (`RSA(wph)`) rather than the
  current one (`RSA(jph)`). Not used by EF (no `phase` dose forms).
- **Hydrated formulas (`setupphases.m` `setstoichiometry`).** `f{ws}=[...]`
  only works when exactly one formula has a `:nH2O` part; two or more would
  error. EF has one (nesquehonite).
- **Feedstock name case (`markPSDinit.m` `getgoodandready`).** It matches
  timetable treatments to feedstocks with case-sensitive `strcmp`, so EF's
  `Blueridge` never matches `blueridge` and no particle-size tracking is set
  up. Harmless for the published runs (`SAmethod='traditional'`), but the
  `discrete` method would silently track nothing.
- **Crop start day (`setupdatastruc.m` `setcropfromtimetable`).** Timetable
  decimal years use the middle of the day (`doy-0.5`) and run decimal years the
  start, so `run.decyear>=timetable.decyear` puts each crop change on the
  following day, while peak biomass, harvest and sowing depth go on the
  timetable day itself.
- **Maize silicon uptake skipped (`daySWAT7.m` `getnutuptake`).** Crop-specific
  rows in cropstoichiometry.csv are matched with case-sensitive `contains`, so
  `Maize` (Si 2.1% of dry mass, B) never matches the timetable crop `maize` and
  maize takes up no Si. This may matter for the Si budget of silicate weathering.
- **P uptake rescaling (`daySWAT7.m` `getcellup`).** `celluptake.P(1,w)` rescales
  only the top cell (N uses `(:,w)`). In rBcomp the P arrays are then replaced
  by the N-ratio values from cropstoichiometry.csv, so the result is unaffected.
- **Upward-flow mixing fractions (`setmixfracs2.m`, top-down).** Column
  `i+2` of MIX `i` is cell `i+1`'s solution, so upward flow into cell `i` should
  be `upflow(i)/targH2O(i+1)`. The code first sets `inflow(i+1)/targH2O(i+1)`
  (and `upflow(2)/targH2O(2)` for MIX 1), then overwrites MIXes 1 to n-2 with
  `upflow(i)/targH2O(i)`, the wrong cell's water. MIX n-1 keeps
  `inflow(n)/targH2O(n)`, which is downward inflow, not upflow. Nothing is
  taken out of cell `i+1` for the water that moves up.
- **Bottom soil layer misses its inflow (`soilwatermovement.m`,
  `soilwater_zengdecker2009`).** When the water table is below the soil column
  (as in rBcomp), the flux into the bottom layer is stored in a local `qin`
  instead of `c.qin`, so `rmx(nlevsoi)` uses 0 for the inflow and water
  leaving layer `nlevsoi-1` never arrives. In rBcomp this is the deepest added
  CLM5 layer (about 8 m), so the top 8 cells are not directly affected.
- **Added CLM5 layers (`soilwatermovement.m` `getds`).** The overlap fix for
  the first layer below the soil cells is applied at index `nlevsoi+1` (after
  `nlevsoi` already counts the added layers), so the layer at 1.83 to 2.08 m
  overlaps the bottom soil cell and a deeper layer gets a 0.57 m thickness.
  Initial water in the added layers uses `clm5.dz(nwadd)` (one thickness for
  all of them) rather than `clm5.dz(wadd)`. CLM5 `dz` is also built as the
  distance between node depths, not interface depths.
- **Frost table (`soilwatermovement.m` `WaterTable`).** The search sets `kfrz`
  instead of `k_frz`, so the frost table is always the top or bottom layer.
  (`Drainage` has it right.)
- **Day-1 soil moisture is zero (`getrichquick.m`).** `run.Sws` and
  `run.Swstot` are only filled from day 2, so day 1 has zero moisture,
  `airpor = por` and `sat = 0`.
- **Frozen cells have almost no water.** Frozen water is held as ice inside
  `soilwatermovement`, and `run.targH2O` is the liquid part only. In rBcomp
  cells 2 and 3 (10 to 38 cm) start frozen and are refilled to 0.01 mm
  (`watmin`) from below, so PHREEQC sees 0.01 mm of water there. This may be
  intended.
- **Last-day terms in `SL` (`getsoilCO2g.m`).** Inside the day loop
  `p.SL=p.phi0.*...` uses the whole `phi0` row with day `i`'s `cdepth` and
  `Dsoil`, so only the last day's values survive. `SL` is only used when
  `correctforfinitesoildepth` is set.
- **CEC scaling ignored (`CECsorb.m`).** If `ds.scale` exists without a `CEC`
  field (EF has `ds.scale.pCO2g`), `parm.CECscl` is ignored.

## MATLAB call graph (rBcomp entry points)

```
  257 examples/EF2/EFrunsettings2.m  -> bucket
  947 setupdatastruc.m  -> CECsorb,Cerling1991,defaultsettings,getcellsfromlayers,getfert,getperiodictable,markPSDinit,potET,readinputfiles,setmodelspecies,setruncellindices,setupnutuptake,setupphases,setuprun,writePHASESRATES
 3846 getphreeqcrun.m  -> CECsorb,assessrun2,checkphasetype,gasmoles,getperiodictable,getselectedoutput,gettotalfields,reshapeselectedoutput,updateRSA
    4 examples/EF2/BRscalings.m  -> 
 1395 bucket.m  -> arrayrange,showmixfracs
  743 CECsorb.m  -> 
  542 Cerling1991.m  -> arraystats,getfig
  166 examples/EF2/defaultsettings.m  -> 
  473 getcellsfromlayers.m  -> pedotransfer,setupsoilphases
  187 getfert.m  -> 
   57 getperiodictable.m  -> readcsvdata
  259 markPSDinit.m  -> 
  361 potET.m  -> 
  715 readinputfiles.m  -> getphasetable,readinitlayerdata,readobsdata,readtimetable
   97 setmodelspecies.m  -> getperiodictable
   93 setruncellindices.m  -> 
   66 setupnutuptake.m  -> readcsvdata
 1180 setupphases.m  -> checkphasetype,getRosinRammlerPSD,readcsvdata,setnucSA,setstoichiometry
  488 setuprun.m  -> daySWAT7,getsoilCO2g,getsoilresp,setmixfracs2
  303 writePHASESRATES.m  -> 
  489 assessrun2.m  -> checkuselog,getperiodictable
  392 checkphasetype.m  -> implicitcalc
  149 gasmoles.m  -> 
  211 getselectedoutput.m  -> getfig
 1008 gettotalfields.m  -> assessuptake,etaRH,getbasesatpct
  123 reshapeselectedoutput.m  -> getselectedoutput
  315 updateRSA.m  -> 
   45 arrayrange.m  -> 
  427 showmixfracs.m  -> arrayrange,bucket,getfig,showhydrology
   27 arraystats.m  -> 
   22 getfig.m  -> 
  266 pedotransfer.m  -> 
  253 setupsoilphases.m  -> checkphasetype,setupphases
   73 readcsvdata.m  -> 
  365 getphasetable.m  -> checkphasetype
  431 readinitlayerdata.m  -> getperiodictable,setmodelspecies,setupsoilphases
  347 readobsdata.m  -> getperiodictable,getscreenedtable
  148 readtimetable.m  -> getdecyear,getscreenedtable
   64 getRosinRammlerPSD.m  -> 
  120 setnucSA.m  -> 
   54 setstoichiometry.m  -> 
  975 daySWAT7.m  -> SWATsnow,potET,readcsvdata
  352 getsoilCO2g.m  -> 
  342 getsoilresp.m  -> arraystats,getfig
  488 setmixfracs2.m  -> bucket,getrichquick
   53 checkuselog.m  -> 
  322 implicitcalc.m  -> checkphasetype,setupphases
  139 assessuptake.m  -> 
   39 etaRH.m  -> 
   32 getbasesatpct.m  -> 
  188 showhydrology.m  -> arrayrange,getfig
  121 getscreenedtable.m  -> 
   67 getdecyear.m  -> 
   68 SWATsnow.m  -> getfig
  168 getrichquick.m  -> drainrate,soilwatermovement
   40 drainrate.m  -> 
 2169 soilwatermovement.m  -> bucket,gettablefromstruct
   84 gettablefromstruct.m  -> 
total 23155 files 58
```
