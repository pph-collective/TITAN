import pytest
import numpy as np

import titan.probabilities as probs


class FakeRandom:
    def __init__(self, num: float):
        self.num = num

    def random(self):
        return self.num

    def randrange(self, start, stop, step=1):
        return start

    def randint(self, start, stop):
        return start


def test_unsafe_sex():
    # initiate result dict with 2 time steps
    assert probs.unsafe_sex(0) == 0.443
    assert probs.unsafe_sex(1) == 0.481
    assert probs.unsafe_sex(5) == 0.514
    assert probs.unsafe_sex(100) == 0.759


def test_adherence_prob():
    # initiate result dict with 2 time steps
    assert probs.adherence_prob(1) == 0.0051
    assert probs.adherence_prob(2) == 0.0039
    assert probs.adherence_prob(3) == 0.0032
    assert probs.adherence_prob(4) == 0.0025
    assert probs.adherence_prob(5) == 0.0008
    assert probs.adherence_prob(6) == 0.0051


def test_get_death_rate():
    for hiv in [True, False]:
        for aids in [True, False]:
            for race in ["WHITE", "BLACK"]:
                for adh in [0, 1]:
                    assert probs.get_death_rate(hiv, aids, race, adh) > 0


def test_get_mean_num_partners():
    for du in ["IDU", "NIDU"]:
        for i in np.arange(0, 1, 0.01):
            assert probs.get_mean_num_partners(du, FakeRandom(i)) > 0
