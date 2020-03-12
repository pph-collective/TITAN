import random
from typing import Sequence, TypeVar, Optional

import numpy as np #type: ignore


def get_check_rand_int(seed):
    """
    Check the value passed of a seed, make sure it's an int, if 0, get a random seed
    """
    if type(seed) is not int:
        raise ValueError("Random seed must be integer")
    elif seed == 0:
        return random.randint(1, 1000000)
    else:
        return seed


def safe_divide(numerator: int, denominator: int):
    """
    Default 0 if denominator is 0, otherwise divide as normal
    """
    if denominator == 0:
        return 0.0
    else:
        return 1.0 * numerator / denominator


# Requirement for safe_random_choice function
T = TypeVar("T")


def safe_random_choice(seq: Sequence[T], rand_gen) -> Optional[T]:
    """
    Return None or a random choice
    """
    if seq:
        return rand_gen.choice(seq)
    else:
        return None

def binom(k, n, p, np_rand):
    """
        mirrors scipy binom.pmf function but from random distribution from numpy
    """
    samples = 10000
    return sum(np_rand.binomial(n, p, samples) == k)/samples

def poisson(mu, np_rand, size=1):
    """
        mirrors scipy poisson.rvs function as used in code
    """
    return round(np_rand.poisson(mu, size)[0])
