# Rice paddy chemistry

The flooded-soil reactions of David Beerling's SCEPTER paddy module, rewritten as
PHREEQC `RATES` from the process spec (`paddy/docs/riceland_process_spec.md` on the
SCEPTER branch `paddy/fe3-reduction`, v15, commit 771da05). The spec is the source,
not the Fortran. Step 1 of the design note is one cell, compared against a SCEPTER run. Step 2 puts the
rates in every cell of an ARTEMIS column (`artemis_p/paddy.py`).

| File | What it is |
|---|---|
| `paddy.dat` | INCLUDE$ file: decoupled CH4 (`Mtg`), O2 (`Oxg`) and sulfide (`Hsg`), and the RATES |
| `kinetics_block.py` | the KINETICS block: each reaction's net change to the solution, and the solid pools |
| `batch.py` | a PHREEQC batch cell set up like cell 1 of a SCEPTER reference run |
| `read_scepter.py` | reads a SCEPTER run's profile files |
| `figure.py` | the batch comparison figure |
| `column_figure.py` | the step 2 figure: three depths of an ARTEMIS column through a season |

## How the reactions are written

Every reaction is a KINETICS entry whose `-formula` is its net change to the
solution, balanced in mass and oxidation state. PHREEQC therefore never has to
move electrons between pathways, and its pe cannot override the kinetic redox
ladder. CH4, O2 and sulfide are decoupled master species. N was already decoupled
in `THAmmOrgAlP.dat`. Fe stays coupled: every formula adds Fe at the right
valence. The solids (OM, ferrihydrite, green rust, goethite, FeS and urea) are
KINETICS pools that do not change the solution. `Paddy_rates` evaluates every
rate law once per call, in SCEPTER units, and PUTs the results. The reaction
entries and the pools then GET their own rates. `tests/test_paddy.py` checks that
each formula gives exactly its products.

## Result (113-day flood, cell at 1 cm)

`python paddy/batch.py OUT scx` runs the cell. The `scx` option uses SCEPTER's
exchange selectivities and leaves out NH4+ exchange, as SCEPTER does.
`python paddy/figure.py SCEPTER_RUNDIR OUT.sel fig.png` draws the comparison.

These quantities agree closely:

- the acceptor sequence: O2 is gone by day 7, the nitrate peak is about 80 µM on day 6, sulfate is gone by day 30;
- CH4 onset;
- green rust;
- FeS formation;
- N2O.

CH4 at day 112 is 2.8 atm, against 2.5 in SCEPTER. SCEPTER's cell loses CH4 to
the cells below by diffusion, and the batch cell has no neighbours.

pH matches to day 30, but by day 112 it is 7.11, against 6.94 in SCEPTER.
Dissolved Fe(II) is a third of SCEPTER's. Both differences come from the
exchanger: SCEPTER's CEC depends on pH, and PHREEQC's exchanger here has a fixed
CEC. In ARTEMIS the exchanger is ARTEMIS's own, so a match on these two is not
expected.

With the database's own exchanger, NH4+ goes onto the exchange sites, so
nitrification and the loss of O2 are slower than in SCEPTER.

## Step 2: the paddy chemistry in an ARTEMIS column

`ds['paddy']` turns it on (options in the `artemis_p/paddy.py` docstring), and
`python scripts/run_paddy.py rCcomp WORKDIR NDAYS floodday=152,flooddays=113` runs the
control column as a paddy. In paddy mode:

- **Hydrology.** On flooded days every cell is saturated and `q` (0.5 m/yr) percolates out
  of the bottom. The Richards solution stands on drained days. Irrigation makes up what
  the soil needs, and rain the soil cannot take runs off.
- **Gases.** Flooded cells have no gas phase. CH4 leaves by ebullition above 0.4 atm, and
  the top cell exchanges O2, CO2 and CH4 with the atmosphere through the floodwater
  (spec 6.1 and 6.2). Both are kinetic: a gas phase that switched on and off within a
  day's integration stalled the solver. Drained cells keep ARTEMIS's soil CO2 and take O2
  as Oxg(g), so Fe(II) and CH4 oxidation stay kinetic. New solutions carry Oxg, not
  coupled O(0), so nothing oxidises Fe(II) at equilibrium.
- **Pools** are carried from day to day in moles per cell, and the rates see them per
  m3 of soil whatever water the cell holds.
- SWAT's OM and N kinetics are off; crop uptake stays. Paddy runs use CVODE.
- **Not in SCEPTER:** reactions ramp to zero below 15% saturation and stop below 5%.
  ARTEMIS's frozen cells keep 0.01 mm of liquid, and the full rates in that water stalled
  the solver.

Runs with the executable (the default engine) and no cell workers. The in-memory engine
with cellworkers=4 failed on the surface model's convergence in some cells.

`paddy_column_winter_test.png` is a test of the flood and drain transitions: 9 drained
days, a 113-day flood from 10 January, then drained (rCcomp, 2015). Soil temperature
averaged 5-6 C over the flood, so the chemistry is slow. The acceptor sequence is right:
O2 goes first, then nitrate, then Fe(III) is reduced, and CH4 builds slowly. The 5 cm cell
(10 cm thick) stays oxic because floodwater O2 reaches it all. The nitrate jump on day
~118 is the season's UAN28 fertiliser. Over the flood, 0.09 mol CH4 m-2 (14 kg ha-1) was
produced. The SCEPTER line is a closed flood at 25 C in a 1 cm cell, so it is a guide
only.

`paddy_column_summer.png` is the summer season: a 113-day flood from 1 June 2015 (soil
22-24 C), drained before and after. Over the flood the column made 0.86 mol CH4 m-2
(138 kg ha-1). The spec calibrates SCEPTER's emission, with plants, to 132-135 kg ha-1.
Ebullition vented 0.49 mol m-2 (78 kg ha-1), and almost none was oxidised. The rest
stays dissolved until the drain, when it degasses. Without plants (step 3) there is no
plant transport. Nitrate (fertiliser) goes in about three weeks, then Fe(II) rises to
about 0.5 mM and CH4 reaches the 0.4 atm ebullition threshold. The soil is more acid
than the SCEPTER reference site (pH 6.0-6.5, against 7-7.6), so green rust does not form
and the pH stays low. The rain carries little sulfate, so little FeS forms. These are
differences between the sites. The top cell stays oxic, as in the winter test.
