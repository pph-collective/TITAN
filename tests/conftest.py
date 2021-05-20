import pytest

import os
import shutil

from titan.parse_params import create_params
from titan.population import Population
from titan.model import TITAN
from titan.agent import Agent, Relationship
from titan.location import Location

# helper method to generate a fake number deterministically
class FakeRandom:
    def __init__(self, num: float, fake_choice: int = 0):
        self.num = num
        self.fake_choice = fake_choice

    def random(self):
        return self.num

    def randrange(self, start, stop, step=1):
        return start

    def sample(self, seq, rate):
        return seq

    def choice(self, seq):
        return seq[0]

    def choices(self, seq, weights=None, k=1):
        if weights is None:
            return [seq[self.fake_choice]]
        else:
            selection = weights.index(max(weights))
            return [seq[selection]]

    def randint(self, start, stop, size=1):
        return start

    def shuffle(self, seq):
        return seq

    def poisson(self, var: float, size: int = 1):
        return int(round(var))


# test fixtures used throughout unit tests
@pytest.fixture
def params(tmpdir):
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "basic.yml"
    )
    return create_params(None, param_file, tmpdir)


@pytest.fixture
def world_location(params):
    return Location("world", params.classes.locations.world, params)


@pytest.fixture
def make_agent(params, world_location):
    def _make_agent(
        SO="MSM", age=30, race="black", DU="None", location=None, init_bond_fields=True
    ):
        if location is None:
            location = world_location
        agent = Agent(SO, age, race, DU, location)
        if init_bond_fields:
            for bond in params.classes.bond_types:
                agent.target_partners[bond] = 0
                agent.mean_num_partners[bond] = 0
                agent.partners[bond] = set()

        return agent

    return _make_agent


@pytest.fixture
def make_relationship():
    def _make_relationship(id1, id2, bond_type="Sex", duration=2):
        return Relationship(id1, id2, duration, bond_type)

    return _make_relationship


@pytest.fixture
def make_model(params):
    def _make_model(p=params):
        return TITAN(p)

    return _make_model


@pytest.fixture
def make_population(params):
    def _make_population(n=0, p=params):
        p.model.num_pop = n
        return Population(p)

    return _make_population


@pytest.fixture
def setup_results_dir():
    outfile_dir = os.path.join(os.getcwd(), "results", "network")
    if os.path.isdir(outfile_dir):
        shutil.rmtree(outfile_dir)
    os.makedirs(outfile_dir)
    yield outfile_dir
    shutil.rmtree(outfile_dir)
