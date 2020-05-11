import random
from functools import wraps
from typing import TypeVar, Optional, Collection
from math import ceil


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


def safe_shuffle(seq: Collection[T], rand_gen):
    """
    Return None or a shuffled sequence
    """
    if seq:
        if isinstance(seq, set):
            rand_gen.shuffle(list(seq))
            return seq
        else:
            rand_gen.shuffle(seq)
            return seq
    else:
        return None


def safe_dist(dist_info, rand_gen, dist_type=None):
    if not dist_type:
        dist_type = dist_info.distribution
    print(dist_info)

    print("\n\n", dist_type)

    if dist_type == "set_value":
        return int(dist_info.var_1)
    else:
        dist = getattr(rand_gen, dist_type)

    try:
        value = dist(dist_info.var_1, dist_info.var_2)
    except TypeError:  # If second param of function is shape, must be int
        value = dist(dist_info.var_1, int(dist_info.var_2))

    if hasattr(value, "__iter__"):  # check if value is any type of sequence
        return value[0]
    else:
        return value


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
    decorator to memoize a function
    caches results given args, only use if deterministic)
    """
    cache = {}

    @wraps(f)
    def wrap(*arg):
        if arg not in cache:
            cache[arg] = f(*arg)
        return cache[arg]

    return wrap
