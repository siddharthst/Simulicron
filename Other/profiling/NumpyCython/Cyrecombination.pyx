import numpy as np
cimport numpy as np
import cython
@cython.boundscheck(False)

cpdef recombination(np.ndarray v1, np.ndarray v2, np.ndarray rates):
    assert len(v1) == len(v2) and len(v1) == (
        len(rates) + 1), "Length mismatch"
    rec = np.random.uniform(size=len(rates)) < rates
    start = [0 if (np.random.uniform() < 0.5) else 1]
    whichcol = 1 + np.cumsum(np.concatenate((start, rec))) % 2
    return (np.where(whichcol == 1, v1, v2))
