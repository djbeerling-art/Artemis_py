"""Soil gas diffusion step (artemis_p.gasdiff)."""

import numpy as np

from artemis_p import gasdiff


def test_step_steady_state_matches_analytic():
    # uniform column, source spread evenly: C(z) = Ctop + S/(D L) (L z - z^2/2) at the cell centres
    n, dz, D, cap = 40, 5.0, 0.01, 0.2
    L = n * dz
    q = 1e-9  # mol cm-3 s-1
    S = np.full(n, q * dz)
    C = np.zeros(n)
    for _ in range(200):
        C = gasdiff.step(C, np.full(n, D), np.full(n, cap), np.full(n, dz), S, 0.0, dt=1e7)
    z = (np.arange(n) + 0.5) * dz
    exact = q / D * (L * z - z ** 2 / 2)
    assert np.allclose(C, exact, rtol=0, atol=1e-3 * exact.max())  # the half-cell top boundary is first order


def test_step_conserves_mass_when_closed():
    rng = np.random.default_rng(1)
    n = 8
    dz, D, cap = rng.uniform(5, 50, n), rng.uniform(1e-4, 1e-2, n), rng.uniform(0.05, 0.4, n)
    C0 = rng.uniform(1e-6, 1e-5, n)
    S = rng.uniform(0, 1e-9, n)
    C = gasdiff.step(C0, D, cap, dz, S, 0.0, closed=True)
    assert np.isclose(np.sum(cap * dz * (C - C0)), S.sum() * gasdiff.SEC, rtol=1e-9)


def test_wet_soil_slows_diffusion():
    TK = np.full(2, 293.15)
    por = np.full(2, 0.5)
    D, cap = gasdiff.properties("CO2", TK, np.array([0.3, 0.05]), np.array([0.2, 0.45]), por)
    assert D[1] < D[0] / 50  # Millington-Quirk: strongly non-linear in air-filled porosity
    Dm, _ = gasdiff.properties("CO2", TK, np.array([0.3, 0.05]), np.array([0.2, 0.45]), por, "Moldrup")
    assert Dm[1] > D[1]
    assert np.all(cap > np.array([0.3, 0.05]))  # dissolved CO2 adds storage
