import pytest
import os

from titan.partnering import *
from titan.agent import Agent, Relationship
from titan.population import Population
from titan.parse_params import create_params, ObjMap


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

    def poisson(self, var: float, size: int):
        return var


def test_partnership_duration(params):
    # test duration with randint
    assert get_partnership_duration(params, FakeRandom(1.0), "Inj") == 1
    # test duration with bins
    assert get_partnership_duration(params, FakeRandom(0.1), "Sex") == 1


def test_get_random_pwid_partner_no_PWID(make_population, make_agent, params):
    empty_pop = make_population()
    idu_agent = make_agent(DU="Inj")
    assert idu_agent
    nidu_agent = make_agent()
    for bond in params.classes.bond_types:
        idu_agent.target_partners[bond] = 0
        nidu_agent.target_partners[bond] = 0
        idu_agent.partners[bond] = set()
        nidu_agent.partners[bond] = set()

    empty_pop.add_agent(idu_agent)
    empty_pop.add_agent(nidu_agent)

    partner = select_partner(
        idu_agent,
        empty_pop.all_agents.members,
        empty_pop.sex_partners,
        empty_pop.pwid_agents,
        params,
        FakeRandom(1.0),
        "Inj",
    )
    assert partner is None


def test_get_random_pwid_partner_w_PWID(make_population, make_agent, params):
    empty_pop = make_population()
    idu_agent = make_agent(DU="Inj")
    idu_partner = make_agent(DU="Inj")
    for bond in params.classes.bond_types:
        idu_agent.target_partners[bond] = 0
        idu_partner.target_partners[bond] = 0
        idu_agent.partners[bond] = set()
        idu_partner.partners[bond] = set()
    idu_agent.target_partners["Inj"] = 10
    idu_partner.target_partners["Inj"] = 10
    empty_pop.add_agent(idu_agent)
    empty_pop.add_agent(idu_partner)

    for race in params.classes.races:
        for bond in copy(params.classes.bond_types):
            if bond != "Sex":
                params.demographics[race]["PWID"]["num_partners"][bond]["var_1"] = 0.0

    partner = select_partner(
        idu_agent,
        empty_pop.all_agents.members,
        empty_pop.sex_partners,
        empty_pop.pwid_agents,
        params,
        FakeRandom(1.0),
        "Inj",
    )
    assert partner == idu_partner


def test_get_random_sex_partner_valid(make_population, make_agent, params):
    empty_pop = make_population()
    hm_agent = make_agent(SO="HM")
    hf_partner = make_agent(SO="HF")
    for bond in params.classes.bond_types:
        hm_agent.target_partners[bond] = 0
        hf_partner.target_partners[bond] = 0
        hm_agent.partners[bond] = set()
        hf_partner.partners[bond] = set()
    empty_pop.add_agent(hm_agent)
    empty_pop.add_agent(hf_partner)

    hm_agent.target_partners = 10
    hf_partner.target_partners = 10

    partner = select_partner(
        hm_agent,
        empty_pop.all_agents.members,
        empty_pop.sex_partners,
        empty_pop.pwid_agents,
        params,
        FakeRandom(1.0),
        "Sex",
    )
    assert partner == hf_partner

    rel = Relationship(partner, hm_agent, 10, "Sex")
    empty_pop.add_relationship(rel)

    # no match after bonded
    partner = select_partner(
        hm_agent,
        empty_pop.all_agents.members,
        empty_pop.sex_partners,
        empty_pop.pwid_agents,
        params,
        FakeRandom(1.0),
        "Sex",
    )
    assert partner is None


def test_get_random_sex_partner_bad(make_population, make_agent, params):
    empty_pop = make_population()
    hm_agent = make_agent(SO="HM")
    msm_partner = make_agent(SO="MSM")
    for bond in params.classes.bond_types:
        hm_agent.target_partners[bond] = 0
        msm_partner.target_partners[bond] = 0
        hm_agent.partners[bond] = set()
        msm_partner.partners[bond] = set()
    empty_pop.add_agent(hm_agent)
    empty_pop.add_agent(msm_partner)

    partner = select_partner(
        hm_agent,
        empty_pop.all_agents.members,
        empty_pop.sex_partners,
        empty_pop.pwid_agents,
        params,
        FakeRandom(1.0),
        "Sex",
    )
    assert partner is None


def test_get_assort_partner_race(make_population, make_agent, params):
    pop = make_population()
    a = make_agent(SO="MSM", race="WHITE")
    p1 = make_agent(SO="MSM", race="WHITE")
    p2 = make_agent(SO="MSM", race="BLACK")
    for bond in params.classes.bond_types:
        a.target_partners[bond] = 1
        p1.target_partners[bond] = 1
        p2.target_partners[bond] = 1
        a.partners[bond] = set()
        p1.partners[bond] = set()
        p2.partners[bond] = set()
    pop.add_agent(a)
    pop.add_agent(p1)
    pop.add_agent(p2)

    params.features.assort_mix = True

    # assrot with WHITE
    test_rule = ObjMap(
        {
            "attribute": "race",
            "agent_value": "WHITE",
            "partner_values": {"WHITE": 0.9, "__other__": 0.1},
        }
    )
    params.assort_mix["test_rule"] = test_rule

    partner = select_partner(
        a,
        pop.all_agents.members,
        pop.sex_partners,
        pop.pwid_agents,
        params,
        FakeRandom(0.5),
        "Sex",
    )

    assert partner == p1

    # get __other__
    test_rule = ObjMap(
        {
            "attribute": "race",
            "agent_value": "WHITE",
            "partner_values": {"WHITE": 0.9, "__other__": 10},
        }
    )
    params.assort_mix["test_rule"] = test_rule

    partner = select_partner(
        a,
        pop.all_agents.members,
        pop.sex_partners,
        pop.pwid_agents,
        params,
        FakeRandom(0.5),
        "Sex",
    )

    assert partner == p2


def test_get_assort_partner_high_risk(make_population, make_agent, params):
    pop = make_population()
    a = make_agent(SO="MSM", race="WHITE")
    p1 = make_agent(SO="MSM", race="WHITE")
    p2 = make_agent(SO="MSM", race="BLACK")

    a.high_risk = True
    p1.high_risk = True
    p2.high_risk = False

    for bond in params.classes.bond_types:
        a.target_partners[bond] = 0
        p1.target_partners[bond] = 0
        p2.target_partners[bond] = 0
        a.partners[bond] = set()
        p1.partners[bond] = set()
        p2.partners[bond] = set()
    a.target_partners["Sex"] = 1
    p1.target_partners["Sex"] = 1
    p2.target_partners["Sex"] = 1
    pop.add_agent(a)
    pop.add_agent(p1)
    pop.add_agent(p2)

    params.features.assort_mix = True

    # assrot with high_risk
    test_rule = ObjMap(
        {
            "attribute": "high_risk",
            "agent_value": True,
            "partner_values": {"True": 0.9, "__other__": 0.1},
        }
    )
    params.assort_mix["test_rule"] = test_rule

    partner = select_partner(
        a,
        pop.all_agents.members,
        pop.sex_partners,
        pop.pwid_agents,
        params,
        FakeRandom(0.5),
        "Sex",
    )

    assert partner == p1

    # get __other__
    test_rule = ObjMap(
        {
            "attribute": "high_risk",
            "agent_value": True,
            "partner_values": {"True": 0.9, "__other__": 10},
        }
    )
    params.assort_mix["test_rule"] = test_rule

    partner = select_partner(
        a,
        pop.all_agents.members,
        pop.sex_partners,
        pop.pwid_agents,
        params,
        FakeRandom(0.5),
        "Sex",
    )

    assert partner == p2


def test_get_assort_partner_drug_use(make_population, make_agent, params):
    pop = make_population()
    a = make_agent(SO="MSM", race="WHITE", DU="Inj")
    p1 = make_agent(SO="MSM", race="WHITE", DU="Inj")
    p2 = make_agent(SO="MSM", race="BLACK", DU="None")
    p3 = make_agent(SO="MSM", race="BLACK", DU="NonInj")

    for bond in params.classes.bond_types:
        a.target_partners[bond] = 0
        p1.target_partners[bond] = 0
        p2.target_partners[bond] = 0
        a.partners[bond] = set()
        p1.partners[bond] = set()
        p2.partners[bond] = set()

    pop.add_agent(a)
    pop.add_agent(p1)
    pop.add_agent(p2)

    params.features.assort_mix = True

    # make sure partnering on sex_type
    # params.partnership.bonds["PWID"]["Sex"]["prob"] = 10
    for race in params.classes.races:
        for bond in params.classes.bond_types:
            if bond != "Sex":
                params.demographics[race]["PWID"]["num_partners"][bond]["var_1"] = 0.0

    # assort with Inj
    test_rule = ObjMap(
        {
            "attribute": "drug_use",
            "agent_value": "Inj",
            "partner_values": {"Inj": 0.8, "NonInj": 0.1, "__other__": 0.1},
        }
    )
    params.assort_mix["test_rule"] = test_rule

    partner = select_partner(
        a,
        pop.all_agents.members,
        pop.sex_partners,
        pop.pwid_agents,
        params,
        FakeRandom(0.5),
        "Inj",
    )

    assert partner == p1

    # get __other__
    test_rule = ObjMap(
        {
            "attribute": "drug_use",
            "agent_value": "Inj",
            "partner_values": {"Inj": 0.8, "NonInj": 0.1, "__other__": 10},
        }
    )
    params.assort_mix["test_rule"] = test_rule

    partner = select_partner(
        a,
        pop.all_agents.members,
        pop.sex_partners,
        pop.pwid_agents,
        params,
        FakeRandom(0.5),
        "Sex",
    )
    assert partner == p2


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
