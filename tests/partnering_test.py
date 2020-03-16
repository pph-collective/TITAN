import pytest
import os

from titan.partnering import *
from titan.agent import Agent, Relationship
from titan.population import Population
from titan.parse_params import create_params


@pytest.fixture
def params(tmpdir):
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "basic.yml"
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
    def __init__(self, num: float, fake_choice: int = 0):
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
        if weights is None:
            return [seq[self.fake_choice]]
        else:
            selection = weights.index(max(weights))
            return [seq[selection]]


def test_get_random_pwid_partner_no_PWID(make_population, make_agent, params):
    empty_pop = make_population()
    idu_agent = make_agent(DU="Inj")
    assert idu_agent
    nidu_agent = make_agent()
    empty_pop.add_agent(idu_agent)
    empty_pop.add_agent(nidu_agent)

    partner, bond_type = select_partner(
        idu_agent,
        empty_pop.all_agents,
        empty_pop.sex_partners,
        empty_pop.pwid_agents,
        params,
        FakeRandom(1.0),
    )
    assert partner is None


def test_get_random_pwid_partner_w_PWID(make_population, make_agent, params):
    empty_pop = make_population()
    idu_agent = make_agent(DU="Inj")
    idu_partner = make_agent(DU="Inj")
    idu_agent.target_partners = 10
    idu_partner.target_partners = 10
    empty_pop.add_agent(idu_agent)
    empty_pop.add_agent(idu_partner)

    empty_pop.params.partnership.bonds.PWID.Inj.prob = (
        10  # make this weight really high so it will get chosen
    )

    partner, bond_type = select_partner(
        idu_agent,
        empty_pop.all_agents,
        empty_pop.sex_partners,
        empty_pop.pwid_agents,
        params,
        FakeRandom(1.0),
    )
    assert partner == idu_partner


def test_get_random_sex_partner_valid(make_population, make_agent, params):
    empty_pop = make_population()
    hm_agent = make_agent(SO="HM")
    hf_partner = make_agent(SO="HF")
    empty_pop.add_agent(hm_agent)
    empty_pop.add_agent(hf_partner)

    hm_agent.target_partners = 10
    hf_partner.target_partners = 10

    partner, bond_type = select_partner(
        hm_agent,
        empty_pop.all_agents,
        empty_pop.sex_partners,
        empty_pop.pwid_agents,
        params,
        FakeRandom(1.0),
    )
    assert partner == hf_partner

    rel = Relationship(partner, hm_agent, 10, "Sex")
    empty_pop.add_relationship(rel)

    # no match after bonded
    partner, bond_type = select_partner(
        hm_agent,
        empty_pop.all_agents,
        empty_pop.sex_partners,
        empty_pop.pwid_agents,
        params,
        FakeRandom(1.0),
    )
    assert partner is None


def test_get_random_sex_partner_bad(make_population, make_agent, params):
    empty_pop = make_population()
    hm_agent = make_agent(SO="HM")
    msm_partner = make_agent(SO="MSM")
    empty_pop.add_agent(hm_agent)
    empty_pop.add_agent(msm_partner)

    partner, bond_type = select_partner(
        hm_agent,
        empty_pop.all_agents,
        empty_pop.sex_partners,
        empty_pop.pwid_agents,
        params,
        FakeRandom(1.0),
    )
    assert partner is None


def test_sex_possible(params):
    sex_types = params.classes.sex_types
    # agent sex types are ["HM", "MSM", "WSW", "HF", "MTF"]
    assert sex_possible("HM", "HM", sex_types) == False
    assert sex_possible("HM", "MSM", sex_types) == False
    assert sex_possible("HM", "HF", sex_types) == True
    assert sex_possible("HM", "WSW", sex_types) == False
    assert sex_possible("HM", "MTF", sex_types) == True

    assert sex_possible("MSM", "HM", sex_types) == False
    assert sex_possible("MSM", "MSM", sex_types) == True
    assert sex_possible("MSM", "HF", sex_types) == False
    assert sex_possible("MSM", "WSW", sex_types) == False
    assert sex_possible("MSM", "MTF", sex_types) == True

    assert sex_possible("WSW", "HM", sex_types) == False
    assert sex_possible("WSW", "MSM", sex_types) == False
    assert sex_possible("WSW", "HF", sex_types) == False
    assert sex_possible("WSW", "WSW", sex_types) == True
    assert sex_possible("WSW", "MTF", sex_types) == True

    assert sex_possible("HF", "HM", sex_types) == True
    assert sex_possible("HF", "MSM", sex_types) == False
    assert sex_possible("HF", "HF", sex_types) == False
    assert sex_possible("HF", "WSW", sex_types) == False
    assert sex_possible("HF", "MTF", sex_types) == False

    assert sex_possible("MTF", "HM", sex_types) == True
    assert sex_possible("MTF", "MSM", sex_types) == True
    assert sex_possible("MTF", "HF", sex_types) == False
    assert sex_possible("MTF", "WSW", sex_types) == True
    assert sex_possible("MTF", "MTF", sex_types) == False

    with pytest.raises(ValueError, match=r"Invalid .*_sex_type.*"):
        sex_possible("HM", "XYZ", sex_types)
        sex_possible("XYZ", "HM", sex_types)
