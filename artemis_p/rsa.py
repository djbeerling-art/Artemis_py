"""updateRSA.m: end-of-day reactive surface area of the kinetic phases."""

import numpy as np


def update_rsa(out, meth, ds, irun, phases):
    """Update ``out['run']['RSA']`` for day ``irun`` (1-based) from the day's
    mole loss, and set the next days' starting moles and RSA (with topups).

    Only the 'traditional' and 'inbuilt' methods (shrinking spheres,
    (M/M0)^0.67) and reversible-phase nucleation are ported; the 'continuum'
    and 'discrete' particle-size methods are not.
    """
    pt, parm, run = ds["phasetable"], ds["parm"], ds["run"]
    I = irun - 1
    for ph in phases:
        oldRSA = out["run"]["RSA0"][ph][:, I].copy()
        newRSA = oldRSA.copy()
        M = out["run"]["M"][ph][:, I]
        M0 = out["run"]["M0"][ph][:, I]
        wgood = M0 > 1e-15
        Mrat = np.zeros_like(M0)
        Mrat[wgood] = M[wgood] / M0[wgood]
        Mratexp = 0.67
        if "surfacefractal" in parm:
            Mratexp = (parm["surfacefractal"] - 3) / 3
        elif "roughnessfractal" in parm:
            Mratexp = (parm["roughnessfractal"] - 2) / 3
        # MATLAB's per-phase fractal columns are tested with ismember on all
        # the table's variable names, which is never all true, so are never used
        lam = Mrat ** Mratexp
        if pt.loc[ph, "reversible"]:
            scl = pt.loc[ph, "molmass"] / pt.loc[ph, "G"] / 1e6
            for j in range(ds["ncells"]):
                M0j, Mj = M0[j], M[j]
                if M0j == 0 and Mj > 0:  # nucleation
                    Nj = ds["cell"]["nucN"][j]
                    rj = (scl * Mj / Nj * 3 / 4 / np.pi) ** (1 / 3)
                    newRSA[j] = Nj * 4 * np.pi * rj ** 2
                elif M0j > 0:
                    newRSA[j] = oldRSA[j] * (Mj / M0j) ** 0.67
        elif meth.lower() in ("traditional", "inbuilt"):
            newRSA[wgood] = oldRSA[wgood] * lam[wgood]
        elif meth.lower() in ("continuum", "discrete"):
            raise NotImplementedError(f"updateRSA method {meth!r} is not ported")
        out["run"]["RSA"][ph][:, I] = newRSA
        if irun < ds["nrun"]:
            out["run"]["M0"][ph][:, irun:] = M[:, None] + run["moles"][ph][:, irun:]
            out["run"]["RSA0"][ph][:, irun:] = newRSA[:, None] + run["RSA"][ph][:, irun:]
    return out
