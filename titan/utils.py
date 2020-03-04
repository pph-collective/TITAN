import random
from typing import Sequence, TypeVar, Optional


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
