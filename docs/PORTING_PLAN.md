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
3. **Soil, phases and chemistry setup (~4,000 lines).** `setupphases`,
   `setupsoilphases`, `checkphasetype`, `implicitcalc`, `getcellsfromlayers`,
   `pedotransfer`, `markPSDinit`, `getRosinRammlerPSD`, `setnucSA`,
   `setstoichiometry`, `CECsorb`, `Cerling1991`, `getfert`, `setupnutuptake`,
   `writePHASESRATES`. Checked against the matching `ds` fields.
4. **Climate and hydrology (~6,000 lines).** `setuprun`, `daySWAT7`,
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

## Bugs found in the MATLAB code

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
