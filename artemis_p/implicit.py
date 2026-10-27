"""Implicitly dissolving phases. Python counterpart of implicitcalc.m.

Phases flagged ``implicit`` (e.g. apatite and titanite in the Blue Ridge
feedstock) are lumped into one pseudo-phase per feedstock,
``Implicit_<feedstock>``, that dissolves in step with the feedstock's
explicitly dissolving phases.
"""


def implicit_calc(ds, quiet=False):
    """implicitcalc.m. Only the early exits are ported so far.

    Existing ``Implicit_*`` rows (and the ``impmolrat`` column) are removed
    first. The pseudo-phase is only built once ``ds["stoich"]`` exists, i.e.
    after setupphases; before that implicitcalc returns.
    """
    t = ds["phasetable"]
    imp_rows = [str(p).startswith("Implicit_") or "Implicit_" in str(p) for p in t["PHREEQCphase"]]
    if any(imp_rows):
        t = t[[not r for r in imp_rows]]
        t = t.drop(columns=["impmolrat"], errors="ignore")
        ds["phasetable"] = t
    if not any(str(p).lower() == "implicit" for p in t["PHREEQCphasetype"]):
        return ds
    if "stoich" not in ds:
        if not quiet:
            print("(implicitcalc) Unable to proceed without ds.stoich and ds.run.treatment")
        return ds
    raise NotImplementedError("implicitcalc pseudo-phase construction is not ported yet")
