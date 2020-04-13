import random
from typing import TypeVar, Optional, Collection
from functools import wraps
from math import factorial

import numpy as np  # type: ignore


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


def safe_random_choice(seq: Collection[T], rand_gen) -> Optional[T]:
    """
    Return None or a random choice
    """
    if seq:
        if isinstance(seq, set):
            return rand_gen.choice(tuple(seq))
        else:
            return rand_gen.choice(seq)
    else:
        return None


def safe_shuffle(seq: Collection[T], rand_gen) -> Optional[T]:
    """
    Return None or a shuffled sequence
    """
    if seq:
        if isinstance(seq, set):
            return rand_gen.shuffle(tuple(seq))
        else:
            return rand_gen.shuffle(seq)
    else:
        return None


def binom_0(n: int, p: float):
    """
        mirrors scipy binom.pmf as used in code
    """
    return (1 - p) ** n


def poisson(mu: float, np_rand):
    """
        mirrors scipy poisson.rvs function as used in code
    """
    return np_rand.poisson(mu)


def memo(f):
    """
    decorator to memoize a function (caches results given args, only use if deterministic)
    """
    cache = {}

    @wraps(f)
    def wrap(*arg):
        if arg not in cache:
            cache[arg] = f(*arg)
        return cache[arg]

    return wrap
