import pytest
import os

from titan.ABM_partnering import *
from titan.agent import Agent
from titan.HIVABM_Population import PopulationClass
from titan.params_parse import create_params


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
        return PopulationClass(0, params)

    return _make_population


def test_get_random_PWID_partner_no_PWID(make_population, make_agent):
    empty_pop = make_population()
    idu_agent = make_agent(DU="Inj")
    nidu_agent = make_agent()
    empty_pop.add_agent_to_pop(idu_agent)
    empty_pop.add_agent_to_pop(nidu_agent)
    assert get_random_PWID_partner(idu_agent, empty_pop.All_agentSet) is None


def test_get_random_PWID_partner_w_PWID(make_population, make_agent):
    empty_pop = make_population()
    idu_agent = make_agent(DU="Inj")
    idu_partner = make_agent(DU="Inj")
    empty_pop.add_agent_to_pop(idu_agent)
    empty_pop.add_agent_to_pop(idu_partner)
    assert get_random_PWID_partner(idu_agent, empty_pop.All_agentSet) == idu_partner


def test_get_random_sex_partner_valid(make_population, make_agent, params):
    empty_pop = make_population()
    hm_agent = make_agent(SO="HM")
    hf_partner = make_agent(SO="HF")
    empty_pop.add_agent_to_pop(hm_agent)
    empty_pop.add_agent_to_pop(hf_partner)
    assert (
        get_random_sex_partner(hm_agent, empty_pop.All_agentSet, params) == hf_partner
    )


def test_get_random_sex_partner_bad(make_population, make_agent, params):
    empty_pop = make_population()
    hm_agent = make_agent(SO="HM")
    hf_partner = make_agent(SO="MSM")
    empty_pop.add_agent_to_pop(hm_agent)
    empty_pop.add_agent_to_pop(hf_partner)
    assert get_random_sex_partner(hm_agent, empty_pop.All_agentSet, params) is None


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
