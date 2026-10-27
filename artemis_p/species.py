"""Model solute species. Python counterpart of setmodelspecies.m (default options)."""

from .periodictable import get_periodic_table

GOODELT = ["Ca", "Mg", "Na", "K", "Al", "Fe", "Mn", "Amm", "N", "P", "S", "Cl", "Si", "Ti"]


def set_model_species(pt2=None, include_amm=True):
    """Return the solute dict: one entry per modelled element plus ammonium (``Amm``).

    Only the default call (no ``setgoodelt``/``major``/``all`` options) is ported.
    """
    if pt2 is None:
        _, pt2 = get_periodic_table()
    s = {e: dict(pt2[e]) for e in GOODELT if e in pt2}
    if include_amm and "Amm" not in s:
        s["Amm"] = {
            "solutionspecies": "AmmH+",
            "charge": 1.0,
            "elt": "Amm",
            "molmass": pt2["N"]["molmass"] + 3 * pt2["H"]["molmass"],
            "exchangespecies": "AmmHX",
            "name": "Ammonia",
        }
    s.pop("H", None)
    return s
