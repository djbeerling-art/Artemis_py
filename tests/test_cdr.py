"""getlayerfluxes and getnetCDRfluxes (artemis_p.cdr) on the saved rBcomp run."""

import copy

import numpy as np

from test_inputs import saved_totals_inputs
from test_parallel import needs_data


@needs_data
def test_layer_and_net_fluxes():
    from artemis_p.cdr import get_layer_fluxes, get_net_cdr_fluxes

    run, _, ds = saved_totals_inputs(None)
    out = get_layer_fluxes(copy.deepcopy(run), ds)
    soil = out["soil"]
    # gross exports at the bottom of layer 5 (76 cm), as in the manuscript
    assert abs(soil["cCO2cons"][4, -1] - 1.716) < 1e-3
    assert abs(soil["cCO2cateq"][4, -1] - 3.817) < 1e-3
    # the per-layer solute fluxes use the top cell's concentration unless fixed
    assert np.allclose(soil["Caflux"][4], np.asarray(run["soil"]["Ca"])[0] * np.asarray(ds["run"]["outperc"])[4, :2192],
                       equal_nan=True)
    net = get_net_cdr_fluxes(out, copy.deepcopy(out), ds, ds)
    assert np.all(net["soil"]["netcCO2cons"] == 0) and np.all(net["total"]["netcCO2cateq"] == 0)
