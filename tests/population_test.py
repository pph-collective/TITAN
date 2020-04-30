import pytest
import os

from titan.population import *
from titan.agent import Agent
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
    def _make_population(n=100):
        params.model.num_pop = n
        return Population(params)

    return _make_population


# helper method to generate a fake number deterministically
class FakeRandom:
    def __init__(self, num: float, fake_choice: int = 0):
        self.num = num

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


def test_create_agent(make_population):
    pop = make_population()

    a1 = pop.create_agent("WHITE")
    assert a1.race == "WHITE"
    assert a1.prep_opinion in range(
        5
    ), f"Agents opinion of injectible PrEP is out of bounds {a1.prep_opinion}"

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

    assert num_HM > 70 and num_HM < 130
    assert num_HF > 830 and num_HF < 930
    assert num_PWID > prop_idu - 50 and num_PWID < prop_idu + 50


def test_add_remove_agent_to_pop(make_population):
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

    pop.add_agent(agent)

    assert agent in pop.all_agents.members
    assert agent in pop.hiv_agents.members
    assert agent in pop.high_risk_agents.members

    assert pop.graph.has_node(agent)

    pop.remove_agent(agent)

    assert agent not in pop.all_agents.members
    assert agent not in pop.hiv_agents.members
    assert agent not in pop.high_risk_agents.members

    assert not pop.graph.has_node(agent)


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


def test_update_agent_partners_one_agent(make_population, params):
    pop = make_population(n=1)
    params.model.num_pop = 0

    agent = next(iter(pop.all_agents))  # the only agent in the pop

    for bond in params.classes.bond_types:
        pop.update_agent_partners(agent, bond)  # noMatch == True
    assert agent in pop.graph.nodes()
    assert len(pop.graph.edges()) == 0


def test_update_agent_partners_PWID_no_match(make_population, params):
    pop = make_population(n=0)
    a = pop.create_agent("WHITE", "MSM")
    p = pop.create_agent("WHITE", "HF")
    pop.pop_random = FakeRandom(1.1)
    a.drug_use = "Inj"
    p.drug_use = "Inj"
    pop.add_agent(a)
    pop.add_agent(p)

    p.so = "HF"
    for bond in params.classes.bond_types.keys():
        assert pop.update_agent_partners(a, bond)
        assert a in pop.graph.nodes()
        assert p in pop.graph.nodes()
        assert not a.partners[bond]
        assert len(pop.graph.edges()) == 0

    p.so = "MSM"
    p.drug_use = "None"
    for bond in params.classes.bond_types.keys():
        assert pop.update_agent_partners(a, bond)
        assert a in pop.graph.nodes()
        assert p in pop.graph.nodes()
        assert not a.partners[bond]
        assert len(pop.graph.edges()) == 0


def test_update_agent_partners_MSM_no_match(make_population, params):
    pop = make_population(n=0)
    a = pop.create_agent("WHITE", "MSM")
    p = pop.create_agent("WHITE", "HF")
    pop.pop_random = FakeRandom(1.1)
    a.drug_use = "None"
    p.drug_use = "None"
    pop.add_agent(a)
    pop.add_agent(p)

    assert pop.update_agent_partners(a, "Sex")
    assert a in pop.graph.nodes()
    assert p in pop.graph.nodes()
    assert not a.partners["Sex"]
    assert len(pop.graph.edges()) == 0


def test_update_agent_partners_PWID_match(make_population, params):
    pop = make_population(n=0)
    a = pop.create_agent("WHITE", "MSM")
    p = pop.create_agent("WHITE", "MSM")
    # ensure random sex partner no assorting
    pop.pop_random = FakeRandom(1.1)
    a.drug_use = "Inj"
    p.drug_use = "Inj"
    pop.add_agent(a)
    p.target_partners["Inj"] = 10
    pop.add_agent(p)
    assert pop.partnerable_agents["Inj"]

    no_match = pop.update_agent_partners(a, "Inj")
    assert no_match is False
    assert a in pop.graph.nodes()
    assert p in pop.graph.nodes()
    assert a.partners["Inj"]
    assert len(pop.graph.edges()) == 1


def test_update_agent_partners_MSM_match(make_population, params):
    pop = make_population(n=0)
    a = pop.create_agent("WHITE", "MSM")
    p = pop.create_agent("WHITE", "MSM")
    # ensure random sex partner no assorting
    pop.pop_random = FakeRandom(1.1)
    a.drug_use = "None"
    p.drug_use = "None"
    pop.add_agent(a)
    pop.add_agent(p)

    no_match = pop.update_agent_partners(a, "Sex")
    assert no_match is False
    assert a in pop.graph.nodes()
    assert p in pop.graph.nodes()
    assert a.partners
    assert len(pop.graph.edges()) == 1


def test_update_agent_partners_NDU_PWID_match(make_population, params):
    pop = make_population(n=0)
    a = pop.create_agent("WHITE", "MSM")
    p = pop.create_agent("WHITE", "MSM")
    # ensure random sex partner no assorting
    pop.pop_random = FakeRandom(1.1)
    a.drug_use = "None"
    p.drug_use = "Inj"
    pop.add_agent(a)
    pop.add_agent(p)

    p.target_partners["Sex"] = 10

    no_match = pop.update_agent_partners(a, "Sex")
    assert no_match is False
    assert a in pop.graph.nodes()
    assert p in pop.graph.nodes()
    assert a.partners
    assert len(pop.graph.edges()) == 1


def test_update_partner_assignments_MSM_match(make_population, params):
    pop = make_population(n=0)
    a = pop.create_agent("WHITE", "MSM")
    p = pop.create_agent("WHITE", "MSM")
    # ensure random sex partner no assorting
    pop.pop_random = FakeRandom(1.1)
    a.drug_use = "None"
    p.drug_use = "None"
    pop.add_agent(a)
    pop.add_agent(p)
    a.mean_num_partners["Sex"] = 100
    p.mean_num_partners["Sex"] = 100
    assert params.model.network.enable == True
    assert pop.enable_graph

    pop.update_partner_assignments()
    assert a in pop.graph.nodes()
    assert p in pop.graph.nodes()
    assert a.partners
    assert len(pop.graph.edges()) == 1


def test_update_partner_assignments_PWID_match(make_population, params):
    pop = make_population(n=0)
    a = pop.create_agent("WHITE", "MSM")
    p = pop.create_agent("WHITE", "MSM")
    # ensure random sex partner no assorting
    pop.pop_random = FakeRandom(1.1)
    a.drug_use = "Inj"
    p.drug_use = "Inj"
    pop.add_agent(a)
    pop.add_agent(p)
    a.mean_num_partners["Inj"] = 100
    p.mean_num_partners["Inj"] = 100
    assert params.model.network.enable == True
    assert pop.enable_graph

    pop.update_partner_assignments()
    assert a in pop.graph.nodes()
    assert p in pop.graph.nodes()
    assert a.partners
    assert pop.relationships
    assert len(pop.graph.edges()) == 1


def test_update_partner_assignments_NDU_PWID_match(make_population, params):
    pop = make_population(n=0)
    a = pop.create_agent("WHITE", "MSM")
    p = pop.create_agent("WHITE", "MSM")
    # ensure random sex partner no assorting
    pop.pop_random = FakeRandom(1.1)
    a.drug_use = "None"
    p.drug_use = "Inj"
    pop.add_agent(a)
    pop.add_agent(p)
    a.mean_num_partners["Sex"] = 100
    p.mean_num_partners["Sex"] = 100
    assert params.model.network.enable == True
    assert pop.enable_graph

    pop.update_partner_assignments()
    assert a in pop.graph.nodes()
    assert p in pop.graph.nodes()
    assert a.partners
    assert pop.relationships
    assert len(pop.graph.edges()) == 1


def test_update_partner_assignments_no_match(make_population, params):
    pop = make_population(n=0)
    a = pop.create_agent("WHITE", "MSM")
    a.id = 1
    p = pop.create_agent("WHITE", "HM")
    p.id = 2
    # ensure random sex partner no assorting
    pop.pop_random = FakeRandom(1.1)
    a.drug_use = "None"
    p.drug_use = "None"
    pop.add_agent(a)
    pop.add_agent(p)
    for bond in params.classes.bond_types:
        assert bond in ["Sex", "Inj", "SexInj", "Social"]
        if bond == "Sex":
            a.target_partners[bond] = 50
            p.target_partners[bond] = 50
        else:
            a.target_partners[bond] = 0
            p.target_partners[bond] = 0

    params.model.num_pop = 0
    assert p not in pop.sex_partners[a.so]
    assert a not in pop.sex_partners[p.so]

    pop.update_partner_assignments(1)
    assert not a.partners["Social"]
    assert not p.partners["Social"]
    assert not p.partners["Sex"]
    assert not a.partners["Sex"]
    assert a in pop.graph.nodes()
    assert p in pop.graph.nodes()
    assert pop.all_agents.num_members() == 2
    assert len(pop.graph.edges()) == 0


def test_network_init_scale_free(params):
    """Test if all Inj,NonInj,None drug use agents are in the population"""
    net = Population(params)
    assert params.model.num_pop == net.all_agents.num_members()

    for agent in net.all_agents:
        assert agent in net.graph.nodes()

    for agent in net.all_agents:
        assert agent.drug_use in ["Inj", "NonInj", "None"]
        assert agent.so in params.classes.sex_types


def test_network_init_max_k(params):
    """Test if all Inj,NonInj,None drug use agents are in the population"""
    params.model.network.type = "max_k_comp_size"
    net = Population(params)
    assert params.model.num_pop == net.all_agents.num_members()

    for agent in net.all_agents:
        assert agent in net.graph.nodes()

    for agent in net.all_agents:
        assert agent.drug_use in ["Inj", "NonInj", "None"]
        assert agent.so in params.classes.sex_types


def test_population_consistency(params):
    """Test if Drug users add up"""
    net = Population(params)
    assert net.graph.number_of_edges() == len(net.relationships)
    assert net.graph.number_of_nodes() == net.all_agents.num_members()
    assert params.model.num_pop == net.all_agents.num_members()


def test_population_consistency_HIV(params):
    """Test HIV consistency"""
    net = Population(params)
    for agent in net.all_agents:
        if agent.hiv:
            assert agent in net.hiv_agents

    for agent in net.hiv_agents:
        assert agent.hiv
