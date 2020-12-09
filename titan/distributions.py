from numpy import log  # type: ignore

"""
This file contains distributions that don't exist in numpy.
"""


def set_value(np_random, value):
    """
    A distribution that always returns the value passed

    args:
        np_random: random number generator (to conform to distribution interface)
        value: value to return
    """
    return value


def pert(np_random, low, peak, high, temperature):
    """
    A pert distribution, inspired by [tensorflow](https://github.com/tensorflow/probability/blob/c833ee5cd9f60f3257366b25447b9e50210b0590/tensorflow_probability/python/distributions/pert.py#L137)

    arguments must be so that:

    * low < peak < high
    * temperature > 0

    The support is `[low, high]`.  The `peak` must fit in that interval:
      `low < peak < high`.  The `temperature` is a positive parameter that
      controls the shape of the distribution. Higher values yield a sharper peak.

    args:
        np_random: random number generator (used to get beta)
        low: distribution low value
        peak: modal point in distribution
        high: distribution high value
        temperature: scaling factor
    """
    assert low < peak < high
    assert temperature > 0

    scale = high - low
    alpha = 1.0 + temperature * (peak - low) / scale
    beta = 1.0 + temperature * (high - peak) / scale
    return low + scale * np_random.beta(alpha, beta)


def weibull_modified(np_random, shape, scale):
    """
    Modified version of numpy's (single parameter) weibull distribution to use the 2-parameter weibull.

    args:
        np_random: random number generator
        shape: weibull shape parameter
        scale: weibull scale parameter
    """
    random_number = np_random.random()
    return scale * (-log(1 - random_number)) ** (1 / shape)
