from typing import Tuple
import numpy as np


# todo see if this can be removed to eliminate the numpy dependency
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


class OnlineGaussianSufficientStatistic:
    """
    Implementation of Welford's online mean and variance algorithm
    """

    __slots__ = ['_count', '_mean', '_mean_squared_error']

    def __init__(self):
        self._mean_squared_error: float = 0.
        self._mean: float = 0.
        self._count: int = 0

    def update(self, new_value: float) -> None:
        self._count += 1
        delta = new_value - self._mean
        self._mean += delta / self._count
        delta2 = new_value - self._mean
        self._mean_squared_error += delta * delta2

    @property
    def mean(self) -> float:
        return self._mean

    def calculate_variance(self):
        if self._count < 2:
            return float("nan")
        else:
            return self._mean_squared_error / (self._count - 1)

    def mean_and_variance(self) -> Tuple[float, float]:
        return self.mean, self.calculate_variance()
