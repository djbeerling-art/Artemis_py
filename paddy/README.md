# Rice paddy chemistry (step 1: batch cell)

The flooded-soil reactions of David Beerling's SCEPTER paddy module, rewritten as
PHREEQC `RATES` from the process spec (`paddy/docs/riceland_process_spec.md` on the
SCEPTER branch `paddy/fe3-reduction`, v15, commit 771da05). The spec is the source,
not the Fortran. This is step 1 of the design note: one cell, compared against a
SCEPTER run. The rates are not yet wired into `PhreeqcRun`.

| File | What it is |
|---|---|
| `paddy.dat` | INCLUDE$ file: decoupled CH4 (`Mtg`), O2 (`Oxg`) and sulfide (`Hsg`), and the RATES |
| `kinetics_block.py` | the KINETICS block: each reaction's net change to the solution, and the solid pools |
| `batch.py` | a PHREEQC batch cell set up like cell 1 of a SCEPTER reference run |
| `read_scepter.py` | reads a SCEPTER run's profile files |
| `figure.py` | the comparison figure |

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
