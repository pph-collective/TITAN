import pytest

from titan.ABM_partnering import *
from titan.agent import Agent
from titan.HIVABM_Population import PopulationClass


@pytest.fixture
def make_agent():
    def _make_agent(SO="MSM", age=30, race="WHITE", DU="NDU"):
        return Agent(SO, age, race, DU)

    return _make_agent


@pytest.fixture
def make_population():
    def _make_population(n=0):
        return PopulationClass(n, 0, "PrEP")

    return _make_population


def test_get_random_IDU_partner_no_IDU(make_population, make_agent):
    empty_pop = make_population()
    idu_agent = make_agent(DU="IDU")
    nidu_agent = make_agent()
    empty_pop.add_agent_to_pop(idu_agent)
    empty_pop.add_agent_to_pop(nidu_agent)
    assert get_random_IDU_partner(idu_agent, empty_pop.All_agentSet) is None


def test_get_random_IDU_partner_w_IDU(make_population, make_agent):
    empty_pop = make_population()
    idu_agent = make_agent(DU="IDU")
    idu_partner = make_agent(DU="IDU")
    empty_pop.add_agent_to_pop(idu_agent)
    empty_pop.add_agent_to_pop(idu_partner)
    assert get_random_IDU_partner(idu_agent, empty_pop.All_agentSet) == idu_partner


def test_get_random_sex_partner_valid(make_population, make_agent):
    empty_pop = make_population()
    hm_agent = make_agent(SO="HM")
    hf_partner = make_agent(SO="HF")
    empty_pop.add_agent_to_pop(hm_agent)
    empty_pop.add_agent_to_pop(hf_partner)
    assert get_random_sex_partner(hm_agent, empty_pop.All_agentSet) == hf_partner


def test_get_random_sex_partner_bad(make_population, make_agent):
    empty_pop = make_population()
    hm_agent = make_agent(SO="HM")
    hf_partner = make_agent(SO="MSM")
    empty_pop.add_agent_to_pop(hm_agent)
    empty_pop.add_agent_to_pop(hf_partner)
    assert get_random_sex_partner(hm_agent, empty_pop.All_agentSet) is None


def test_sex_possible():
    # agent sex types are ["HM", "MSM", "WSW", "HF", "MTF"]
    assert sex_possible("HM", "HM") == False
    assert sex_possible("HM", "MSM") == False
    assert sex_possible("HM", "HF") == True
    assert sex_possible("HM", "WSW") == True
    assert sex_possible("HM", "MTF") == True

    assert sex_possible("MSM", "HM") == False
    assert sex_possible("MSM", "MSM") == True
    assert sex_possible("MSM", "HF") == True
    assert sex_possible("MSM", "WSW") == True
    assert sex_possible("MSM", "MTF") == True

    assert sex_possible("WSW", "HM") == True
    assert sex_possible("WSW", "MSM") == True
    assert sex_possible("WSW", "HF") == False
    assert sex_possible("WSW", "WSW") == True
    assert sex_possible("WSW", "MTF") == False

    assert sex_possible("HF", "HM") == True
    assert sex_possible("HF", "MSM") == True
    assert sex_possible("HF", "HF") == False
    assert sex_possible("HF", "WSW") == False
    assert sex_possible("HF", "MTF") == False

    assert sex_possible("MTF", "HM") == True
    assert sex_possible("MTF", "MSM") == True
    assert sex_possible("MTF", "HF") == False
    assert sex_possible("MTF", "WSW") == False
    assert sex_possible("MTF", "MTF") == False

    with pytest.raises(ValueError, match=r"Invalid .*_sex_type.*"):
        sex_possible("HM", "XYZ")
        sex_possible("XYZ", "HM")
