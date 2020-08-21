import random
from functools import wraps
from typing import TypeVar, Optional, Collection

from . import distributions


def get_check_rand_int(seed):
    """
    Check the value passed of a seed, make sure it's an int, if 0, get a random seed
    """
    if type(seed) is not int or seed < 0:
        raise ValueError("Random seed must be positive integer")
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


def safe_random_choice(seq: Collection[T], rand_gen, weights=None) -> Optional[T]:
    """
    Return None or a random choice
    """
    if not seq:
        return None

    if isinstance(seq, set):
        seq = tuple(seq)

    choices = rand_gen.choices(seq, weights=weights)
    return choices[0]


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


def safe_dist(dist_info, rand_gen):
    # gather arguments
    args = []
    for i in range(1, len(dist_info.vars) + 1):
        val = dist_info.vars[i].value
        type_caster = eval(dist_info.vars[i].value_type)
        val = type_caster(val)
        args.append(val)

    dist_type = dist_info.dist_type

    try:  # does dist exist in numpy?
        dist = getattr(rand_gen, dist_type)
        value = dist(*args)
    except AttributeError:
        try:  # does dist exist in distributions.py
            dist = getattr(distributions, dist_type)
            value = dist(rand_gen, *args)
        except AttributeError:
            raise AttributeError(f"Distribution type {dist_type} not found!")

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


def get_param_from_path(params, param_path, delimiter):
    """
    Given a params object and a delimited path, get the leaf of the params tree
    and the last key to access it
    """
    path = param_path.split(delimiter)
    path_params = params
    for p in path[:-1]:
        path_params = path_params[p]

    return path_params, path[-1]


def scale_param(params, param_path, scalar, delimiter="|"):
    """
    Given the params and a parameter path in the format prep|target, scale the
    current value by the scalar
    """
    scaling_item, last_key = get_param_from_path(params, param_path, delimiter)

    old_val = scaling_item[last_key]
    print(f"scaling - {param_path}: {old_val} => {old_val * scalar}")
    scaling_item[last_key] = old_val * scalar


def override_param(params, param_path, value, delimiter="|"):
    """
    Given the params and a parameter path in the format prep|target, scale the
    current value by the scalar
    """
    override_item, last_key = get_param_from_path(params, param_path, delimiter)

    old_val = override_item[last_key]
    print(f"overriding - {param_path}: {old_val} => {value}")
    override_item[last_key] = value
