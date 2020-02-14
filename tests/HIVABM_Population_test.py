import pytest
import os

from titan.HIVABM_Population import *
from titan.agent import Agent
from titan.params_parse import create_params
from titan.network_graph_tools import NetworkClass


@pytest.fixture
def params(tmpdir):
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "basic.yml"
    )
    return create_params({}, param_file, tmpdir)


@pytest.fixture
def make_agent():
    def _make_agent(SO="MSM", age=30, race="WHITE", DU="None"):
        return Agent(SO, age, race, DU)

    return _make_agent


@pytest.fixture
def make_population(params):
    def _make_population(n=100):
        params.model.num_pop = n
        return PopulationClass(0, params)

    return _make_population


# helper method to generate a fake number deterministically
class FakeRandom:
    def __init__(self, num: float):
        self.num = num

    def random(self):
        return self.num

    def randrange(self, start, stop, step=1):
        return start

    def randint(self, start, stop):
        return start


def test_pop_init(make_population):
    n_pop = 100
    pop = make_population(n=n_pop)

    assert pop.All_agentSet.num_members() == n_pop

    # test umbrella sets are all consistent
    parent_sets = ["drugUse_agentSet", "SO_agentSet", "racial_agentSet"]
    for set_name in parent_sets:
        set = getattr(pop, set_name)

        assert set.num_members() == n_pop

        child_pops = 0
        for child_set in set.iter_subset():
            child_pops += child_set.num_members()

        assert child_pops == n_pop


def test_create_agent(make_population):
    pop = make_population()

    a1 = pop.create_agent("WHITE")
    assert a1.race == "WHITE"
    assert a1.opinion in range(
        5
    ), f"Agents opinion of injectible PrEP is out of bounds {a1.opinion}"

    a2 = pop.create_agent("BLACK")
    assert a2.race == "BLACK"

    a3 = pop.create_agent("WHITE", "HM")
    assert a3.so == "HM"
    assert a3.race == "WHITE"

    # check PWID and HIV and high risk
    pop.pop_random = FakeRandom(-0.1)
    a4 = pop.create_agent("WHITE")
    assert a4.drug_use == "Inj"
    assert a4.hiv
    assert a4.aids
    assert a4.hiv_dx
    assert a4.haart
    assert a4.haart_adherence == 5
    assert a4.haart_time == 0
    assert a4.intervention_ever
    assert a4.high_risk
    assert a4.high_risk_ever

    # check not PWID and HIV
    pop.pop_random = FakeRandom(0.999)
    a4 = pop.create_agent("WHITE")
    assert a4.drug_use == "None"
    assert a4.hiv is False
    assert a4.prep is False
    assert a4.intervention_ever is False


def test_create_agent_proportions(make_population, params):
    pop = make_population()

    n = 1000
    race = "WHITE"
    # check proportions
    pop.pop_weights[race] = {"values": ["HM", "HF"], "weights": [0.1, 0.9]}
    prop_idu = round(params.demographics[race]["PWID"].ppl * n)
    num_HM = 0
    num_HF = 0
    num_PWID = 0
    for i in range(n):
        a = pop.create_agent(race)
        if a.drug_use == "Inj":
            num_PWID += 1

        if a.so == "HF":
            num_HF += 1
        elif a.so == "HM":
            num_HM += 1
        else:
            assert False

    assert num_HM > 80 and num_HM < 120
    assert num_HF > 880 and num_HF < 920
    assert num_PWID > prop_idu - 20 and num_PWID < prop_idu + 20


def test_add_agent_to_pop(make_population):
    pop = make_population()
    agent = pop.create_agent("WHITE", "HM")
    agent.drug_use = "Inj"
    agent.hiv = True
    agent.aids = True
    agent.intervention_ever = True
    agent.haart = True
    agent.prep = True
    agent.hiv_dx = True
    agent.incar = True
    agent.high_risk = True

    pop.add_agent_to_pop(agent)

    assert agent in pop.All_agentSet.members
    assert agent in pop.racial_agentSet.members
    assert agent in pop.Race_WHITE_agentSet.members
    assert agent in pop.SO_agentSet.members
    assert agent in pop.SO_HM_agentSet.members
    assert agent in pop.drugUse_agentSet.members
    assert agent in pop.DU_Inj_agentSet.members
    assert agent in pop.HIV_agentSet.members
    assert agent in pop.HIV_AIDS_agentSet.members
    assert agent in pop.treatment_agentSet.members
    assert agent in pop.Trt_ART_agentSet.members
    assert agent in pop.Trt_ART_agentSet.members
    assert agent in pop.Trt_PrEP_agentSet.members
    assert agent in pop.Trt_Tstd_agentSet.members
    assert agent in pop.incarcerated_agentSet.members
    assert agent in pop.highrisk_agentsSet.members

    # check not in all agent sets
    assert agent not in pop.Race_BLACK_agentSet.members


def test_get_age(make_population, params):
    pop = make_population()

    race = "WHITE"

    expected_ages = [15, 25, 35, 45, 55]
    for i in range(1, 6):
        # make sure rand is less than the setting
        pop.pop_random = FakeRandom(params.demographics[race].age[i].prob - 0.001)
        age, ageBin = pop.get_age(race)
        assert age == expected_ages[i - 1]
        assert ageBin == i


def test_update_agent_partners_no_match(make_population, params):
    pop = make_population(n=1)
    params.model.num_pop = 0
    net = NetworkClass(params)

    agent = pop.All_agentSet.members[0]  # the only agent in the pop

    pop.update_agent_partners(net.G, agent)  # noMatch == True
    assert agent in net.G.nodes()
    assert len(net.G.edges()) == 0


def test_update_agent_partners_match(make_population, params):
    pop = make_population(n=0)
    a = pop.create_agent("WHITE", "MSM")
    p = pop.create_agent("WHITE", "MSM")
    pop.add_agent_to_pop(a)
    pop.add_agent_to_pop(p)

    params.model.num_pop = 0
    net = NetworkClass(params)

    pop.update_agent_partners(net.G, a)
    assert a in net.G.nodes()
    assert p in net.G.nodes()
    assert len(net.G.edges()) == 1


def test_update_partner_assignments_match(make_population, params):
    pop = make_population(n=0)
    a = pop.create_agent("WHITE", "MSM")
    p = pop.create_agent("WHITE", "MSM")
    pop.add_agent_to_pop(a)
    pop.add_agent_to_pop(p)
    a.mean_num_partners = 100
    p.mean_num_partners = 100

    params.model.num_pop = 0
    net = NetworkClass(params)

    pop.update_partner_assignments(net.G)
    assert a in net.G.nodes()
    assert p in net.G.nodes()
    assert len(net.G.edges()) == 1


def test_update_partner_assignments_no_match(make_population, params):
    pop = make_population(n=0)
    a = pop.create_agent("WHITE", "MSM")
    p = pop.create_agent("WHITE", "MSM")
    pop.add_agent_to_pop(a)
    pop.add_agent_to_pop(p)

    params.model.num_pop = 0
    net = NetworkClass(params)

    pop.update_partner_assignments(net.G)
    assert a in net.G.nodes()
    assert p in net.G.nodes()
    assert len(net.G.edges()) == 0
