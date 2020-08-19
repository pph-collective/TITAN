import pytest
import os
import random
import numpy as np

import titan.utils as utils
from titan.parse_params import ObjMap

from conftest import FakeRandom


@pytest.mark.unit
def test_check_rand_int():
    assert utils.get_check_rand_int(1) == 1
    assert utils.get_check_rand_int(0) > 0

    with pytest.raises(ValueError):
        utils.get_check_rand_int(3.3)


@pytest.mark.unit
def test_safe_divide():
    assert utils.safe_divide(1, 0) == 0.0
    assert utils.safe_divide(1, 2) == 0.5


@pytest.mark.unit
def test_safe_random_choice():
    rand_gen = random.Random(123)
    assert utils.safe_random_choice([], rand_gen) is None

    assert utils.safe_random_choice({1, 2, 3}, rand_gen) in (1, 2, 3)
    assert utils.safe_random_choice([1, 2, 3], rand_gen) == 1

    assert utils.safe_random_choice([1, 2, 3], rand_gen, [1000, 1, 1]) == 1


@pytest.mark.unit
def test_safe_shuffle():
    rand_gen = random.Random(123)
    assert utils.safe_shuffle([], rand_gen) is None

    assert utils.safe_shuffle({1, 2, 3}, rand_gen) != [1, 2, 3]
    assert utils.safe_shuffle([1, 2, 3], rand_gen) != [1, 2, 3]


@pytest.mark.unit
def test_safe_dist():
    rand_gen = np.random.RandomState(123)

    dist_info = ObjMap(
        {"dist_type": "poisson", "vars": {1: {"value": 20, "value_type": "int"}}}
    )

    # numpy dist
    assert utils.safe_dist(dist_info, rand_gen) > 0

    # custom dist
    dist_info["dist_type"] = "set_value"
    assert utils.safe_dist(dist_info, rand_gen) == 20

    # pert
    low = 2
    peak = 5
    high = 100
    temp = 4
    pert_info = ObjMap(
        {
            "dist_type": "pert",
            "vars": {
                1: {"value": low, "value_type": "int"},
                2: {"value": peak, "value_type": "int"},
                3: {"value": high, "value_type": "int"},
                4: {"value": temp, "value_type": "int"},
            },
        }
    )

    val = utils.safe_dist(pert_info, rand_gen)
    assert val > low
    assert val < high
