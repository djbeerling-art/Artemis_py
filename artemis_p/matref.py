"""Load the authors' saved MATLAB run structures as reference data for the port.

getphreeqcrun.m saves the run struct (named after the run, e.g. ``rBcomp``)
to ``<fname>.mat``. Its ``ds`` field is the full output of setupdatastruc.m,
so it gives intermediate values to test the Python preprocessor against.
"""

import numpy as np
import scipy.io as sio
from scipy.io.matlab import mat_struct


def _convert(obj):
    if isinstance(obj, mat_struct):
        return {f: _convert(getattr(obj, f)) for f in obj._fieldnames}
    if isinstance(obj, np.ndarray) and obj.dtype == object:
        return [_convert(o) for o in obj.ravel()]
    return obj


def load_run(path, name=None):
    """Return the run struct in ``path`` as nested dicts.

    ``name`` is the MATLAB variable name; by default the first non-internal one.
    """
    m = sio.loadmat(path, squeeze_me=True, struct_as_record=False)
    if name is None:
        name = next(k for k in m if not k.startswith("__"))
    return _convert(m[name])
