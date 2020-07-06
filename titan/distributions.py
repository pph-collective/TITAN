"""
This file contains distributions that don't exist in numpy.
"""


def set_value(np_random, value):
    """
    A distribution that always returns the value passed
    """
    return value


def pert(np_random, low, peak, high, temperature):
    """
    A pert distribution, inspired by https://github.com/tensorflow/probability/blob/c833ee5cd9f60f3257366b25447b9e50210b0590/tensorflow_probability/python/distributions/pert.py#L137

    arguments must be so that:
        low < peak < high
        temperature > 0

    The support is `[low, high]`.  The `peak` must fit in that interval:
      `low < peak < high`.  The `temperature` is a positive parameter that
      controls the shape of the distribution. Higher values yield a sharper peak.
    """
    assert low < peak < high
    assert temperature > 0

    scale = high - low
    alpha = 1.0 + temperature * (peak - low) / scale
    beta = 1.0 + temperature * (high - peak) / scale
    return low + scale * np_random.beta(alpha, beta)
