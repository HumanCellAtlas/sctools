import numpy as np


def base4_entropy(x, axis=1):
    """return entropy of x in base 4, calulated across axis.

    Useful for measuring DNA entropy (with 4 nucleotides) as the output is restricted to [0, 1]

    :param np.array x: array of dimension one or more containing numeric types
    :param axis: (default 1) calculate entropy across nucleotides, assuming a 4-column matrix
      where each row consists of base frequencies
    :return np.array: array of input dimension - 1 containin entropy values bounded in [0, 1]
    """

    # convert to probabilities
    if axis == 1:
        x = np.divide(x, np.sum(x, axis=axis)[:, None])
    else:
        x = np.divide(x, np.sum(x, axis=axis))

    with np.errstate(divide='ignore'):
        r = np.log(x) / np.log(4)

    # convention: 0 * log(0) = 0, != -INF.
    r[np.isinf(r)] = 0

    return np.abs(-1 * np.sum(x * r, axis=axis))
