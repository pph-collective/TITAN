import pytest
import os

from titan.partnering import *
from titan.agent import Agent
from titan.population import Population
from titan.parse_params import create_params


@pytest.fixture
def params(tmpdir, paramname="basic.yml"):
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", paramname
    )
    return create_params(None, param_file, tmpdir)


@pytest.fixture
def make_agent():
    def _make_agent(SO="MSM", age=30, race="WHITE", DU="None"):
        return Agent(SO, age, race, DU)

    return _make_agent


@pytest.fixture
def make_population(params):
    def _make_population(n=0):
        params.model.num_pop = n
        return Population(params)

    return _make_population


# helper method to generate a fake number deterministically
class FakeRandom:
    def __init__(self, num: float, fake_choice: int=0):
        self.num = num
        self.fake_choice = fake_choice

    def random(self):
        return self.num

    def randrange(self, start, stop, step=1):
        return start

    def randint(self, start, stop):
        return start

    def choice(self, seq):
        return seq[-1]

    def choices(self, seq, weights=None, k=1):
        return list(seq)[self.fake_choice]

def test_get_random_pwid_partner_no_PWID(make_population, make_agent, params):
    empty_pop = make_population()
    idu_agent = make_agent(DU="Inj")
    assert(idu_agent)
    nidu_agent = make_agent()
    empty_pop.add_agent(idu_agent)
    empty_pop.add_agent(nidu_agent)
    print("++++", select_partner(idu_agent, empty_pop.all_agents, params, FakeRandom(1.0)))

    assert (
        select_partner(idu_agent, empty_pop.all_agents, params, FakeRandom(1.0))[0]
        is None
    )

@pytest.mark.skip
def test_get_random_pwid_partner_w_PWID(make_population, make_agent, params):
    empty_pop = make_population()
    idu_agent = make_agent(DU="Inj")
    idu_partner = make_agent(DU="Inj")
    empty_pop.add_agent(idu_agent)
    empty_pop.add_agent(idu_partner)
    assert (
        get_random_pwid_partner(idu_agent, empty_pop.all_agents, empty_pop.pop_random)
        == idu_partner
    )

@pytest.mark.skip
def test_get_random_sex_partner_valid(make_population, make_agent, params):
    empty_pop = make_population()
    hm_agent = make_agent(SO="HM")
    hf_partner = make_agent(SO="HF")
    empty_pop.add_agent(hm_agent)
    empty_pop.add_agent(hf_partner)
    assert (
        get_random_sex_partner(
            hm_agent, empty_pop.all_agents, params, empty_pop.pop_random
        )
        == hf_partner
    )

@pytest.mark.skip
def test_get_random_sex_partner_bad(make_population, make_agent, params):
    empty_pop = make_population()
    hm_agent = make_agent(SO="HM")
    hf_partner = make_agent(SO="MSM")
    empty_pop.add_agent(hm_agent)
    empty_pop.add_agent(hf_partner)
    assert (
        get_random_sex_partner(
            hm_agent, empty_pop.all_agents, params, empty_pop.pop_random
        )
        is None
    )


def test_sex_possible(params):
    # agent sex types are ["HM", "MSM", "WSW", "HF", "MTF"]
    assert sex_possible("HM", "HM", params) == False
    assert sex_possible("HM", "MSM", params) == False
    assert sex_possible("HM", "HF", params) == True
    assert sex_possible("HM", "WSW", params) == False
    assert sex_possible("HM", "MTF", params) == True

    assert sex_possible("MSM", "HM", params) == False
    assert sex_possible("MSM", "MSM", params) == True
    assert sex_possible("MSM", "HF", params) == False
    assert sex_possible("MSM", "WSW", params) == False
    assert sex_possible("MSM", "MTF", params) == True

    assert sex_possible("WSW", "HM", params) == False
    assert sex_possible("WSW", "MSM", params) == False
    assert sex_possible("WSW", "HF", params) == False
    assert sex_possible("WSW", "WSW", params) == True
    assert sex_possible("WSW", "MTF", params) == True

    assert sex_possible("HF", "HM", params) == True
    assert sex_possible("HF", "MSM", params) == False
    assert sex_possible("HF", "HF", params) == False
    assert sex_possible("HF", "WSW", params) == False
    assert sex_possible("HF", "MTF", params) == False

    assert sex_possible("MTF", "HM", params) == True
    assert sex_possible("MTF", "MSM", params) == True
    assert sex_possible("MTF", "HF", params) == False
    assert sex_possible("MTF", "WSW", params) == True
    assert sex_possible("MTF", "MTF", params) == False

    with pytest.raises(ValueError, match=r"Invalid .*_sex_type.*"):
        sex_possible("HM", "XYZ", params)
        sex_possible("XYZ", "HM", params)
