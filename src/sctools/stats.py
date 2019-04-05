"""
Statistics Functions for Sequence Data Analysis
===============================================

.. currentmodule:: sctools

This module implements statistical modules for sequence analysis

Methods
-------
base4_entropy(x: np.array, axis: int=1)
    calculate the entropy of a 4 x sequence length base frequency matrix

Classes
-------
OnlineGaussianSuficientStatistic        Empirical (online) calculation of mean and variance

"""

from typing import Tuple
import numpy as np


def base4_entropy(x, axis=1):
    """Calculate entropy in base four of a data matrix x

    Useful for measuring DNA entropy (with 4 nucleotides) as the output is restricted to [0, 1]

    Parameters
    ----------
    x : np.ndarray
        array of dimension one or more containing numeric types
    axis : int, optional
        axis to calculate entropy across. Values in this axis are treated as observation frequencies

    Returns
    -------
    entropy : np.ndarray
        array of input dimension - 1 containin entropy values bounded in [0, 1]

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

    Methods
    -------
    update(new_value: float)
        incorporate new_value into the online estimate of mean and variance
    mean()
        return the mean value
    calculate_variance()
        calculate and return the variance
    mean_and_variance()
        return both mean and variance

    """

    __slots__ = ['_count', '_mean', '_mean_squared_error']

    def __init__(self):
        self._mean_squared_error: float = 0.0
        self._mean: float = 0.0
        self._count: int = 0

    def update(self, new_value: float) -> None:
        self._count += 1
        delta = new_value - self._mean
        self._mean += delta / self._count
        delta2 = new_value - self._mean
        self._mean_squared_error += delta * delta2

    @property
    def mean(self) -> float:
        """return the mean value"""
        return self._mean

    def calculate_variance(self):
        """calculate and return the variance"""
        if self._count < 2:
            return float("nan")
        else:
            return self._mean_squared_error / (self._count - 1)

    def mean_and_variance(self) -> Tuple[float, float]:
        """calculate and return the mean and variance"""
        return self.mean, self.calculate_variance()
