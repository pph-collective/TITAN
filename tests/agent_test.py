import pytest
import os

from titan.agent import *
from titan.parse_params import create_params

import random


@pytest.fixture
def params(tmpdir):
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "basic.yml"
    )
    return create_params(None, param_file, tmpdir)


@pytest.fixture
def make_agent():
    def _make_agent(SO="MSM", age=30, race="BLACK", DU="None"):
        return Agent(SO, age, race, DU)

    return _make_agent


@pytest.fixture
def make_relationship():
    def _make_relationship(id1, id2, bond_type="#REVIEW", duration=2):
        return Relationship(id1, id2, duration, bond_type)

    return _make_relationship


# helper method to generate a fake number deterministically
class FakeRandom:
    def __init__(self, num: float):
        assert num >= 0 and num <= 1
        self.num = num

    def random(self):
        return self.num

    def randrange(self, start, stop, step):
        return start


# ============================= AGENT TESTS ============================


def test_agent_init(make_agent):
    a = make_agent()
    b = make_agent()
    assert b.id == a.id + 1

    # demographics
    assert a.so == "MSM"
    assert a.age == 30
    assert a.race == "BLACK"
    assert a.drug_use == "None"
    assert a.age_bin == 0
    assert a.msmw is False

    # partner params
    assert a.relationships == []
    assert a.partners == []
    assert a.mean_num_partners == 0

    # STI params
    assert a.hiv is False
    assert a.hiv_time == 0
    assert a.aids is False

    # treatment params
    assert a.haart is False
    assert a.haart_time == 0
    assert a.haart_adherence == 0
    assert a.sne is False
    assert a.intervention_ever is False
    assert a.prep_reason == []
    assert a.vaccine_time == 0
    assert a.vaccine_type == ""
    assert a.partner_traced is False
    assert a.prep_awareness is False
    assert a.prep_opinion == 0.0
    assert a.prep_type == ""
    assert a.pca is False
    assert a.pca_suitable is False

    # prevention parameters
    assert a.hiv_dx is False
    assert a.prep is False
    assert a.prep_adherence == 0

    # prep pharmacokinetics
    assert a.prep_load == 0.0
    assert a.prep_last_dose == 0

    # high risk params
    assert a.high_risk is False
    assert a.high_risk_time == 0
    assert a.high_risk_ever is False

    # incarceration
    assert a.incar is False
    assert a.incar_ever is False
    assert a.incar_time == 0


def test_partner_list(make_agent, make_relationship):
    a = make_agent()

    assert a.partner_list() == []

    p = make_agent()
    r = make_relationship(a, p)

    assert a.partner_list() == [p.id]
    assert p.partner_list() == [a.id]


def test_get_acute_status(make_agent):
    a = make_agent()  # no HIV on init
    assert a.get_acute_status() == False
    a.hiv_time = 1  # manually force this to test logic
    assert a.get_acute_status() == True


def test_update_prep_load(make_agent, params):
    a = make_agent()
    assert a.prep_last_dose == 0
    assert a.prep_load == 0
    a.update_prep_load(params)
    assert a.prep_last_dose == 1
    assert a.prep_load > 0

    # make time pass
    for i in range(12):
        a.update_prep_load(params)

    assert a.prep_last_dose == 13
    assert a.prep_load == 0.0


def test_get_transmission_probability(make_agent, params):
    a = make_agent(race="WHITE", SO="MSM")
    a.haart_adherence = 1  # set this explicitly

    p_needle = params.partnership.needle.transmission[1].prob
    p_sex = params.partnership.sex.transmission["MSM"][1].prob
    scale = params.calibration.transmission

    # test base case (not tested, not HAART, "WHITE")
    assert a.get_transmission_probability("NEEDLE", params) == p_needle * scale
    assert a.get_transmission_probability("SEX", params) == p_sex * scale

    # test acute
    a.hiv_time = 1
    assert (
        a.get_transmission_probability("SEX", params)
        == p_sex * scale * params.hiv.acute.infectivity
    )
    a.hiv_time = 0

    # test tested status
    a.hiv_dx = True
    assert a.get_transmission_probability("SEX", params) == p_sex * scale * (
        1 - params.hiv.dx.risk_reduction
    )
    a.hiv_dx = False

    # test HAART
    a.haart = True
    assert (
        a.get_transmission_probability("SEX", params)
        == p_sex * scale * params.haart.transmission.prob
    )
    a.haart = False

    # test Black
    a.race = "BLACK"
    assert (
        a.get_transmission_probability("SEX", params)
        == p_sex * scale * params.demographics[a.race].hiv.transmission
    )


def test_get_number_of_sex_acts(make_agent, params):
    a = make_agent()

    rand_gen_low = FakeRandom(0.0)
    min_val_low = params.partnership.sex.frequency[1].min

    rand_gen_high = FakeRandom(1.0)

    assert a.get_number_of_sex_acts(rand_gen_low, params) == min_val_low

    # test fallthrough
    assert a.get_number_of_sex_acts(rand_gen_high, params) == 37


# ============== RELATIONSHIP TESTS ===================


def test_relationship(make_agent, make_relationship):
    a = make_agent()
    p1 = make_agent()
    p2 = make_agent()
    r1 = make_relationship(a, p1)
    r2 = make_relationship(a, p2)

    assert r1.agent1 == a
    assert r1.agent2 == p1

    # properties
    assert r1.duration == 2
    assert r1.total_sex_acts == 0

    assert r2.duration == 2
    assert r2.total_sex_acts == 0

    assert p1.id in a.partner_list()
    assert p2.id in a.partner_list()
    assert a.id in p1.partner_list()
    assert a.id in p2.partner_list()

    assert r1 in a.relationships
    assert r1 in p1.relationships
    assert r2 in a.relationships
    assert r2 in p2.relationships

    # move forward one time step in the relationship, duration 2 -> 1
    ended = r1.progress()
    assert ended == False
    assert r1.duration == 1
    assert p1.id in a.partner_list()
    assert p2.id in a.partner_list()
    assert a.id in p1.partner_list()
    assert a.id in p2.partner_list()

    assert r1 in a.relationships
    assert r1 in p1.relationships
    assert r2 in a.relationships
    assert r2 in p2.relationships

    # move forward one more timestep, duration 1 -> 0, rel over on next progress
    ended = r1.progress()
    assert r1.duration == 0
    ended = r1.progress()
    assert ended == True
    assert r1.duration == 0
    assert p1.id not in a.partner_list()
    assert p2.id in a.partner_list()
    assert a.id not in p1.partner_list()
    assert a.id in p2.partner_list()

    assert r1 not in a.relationships
    assert r1 not in p1.relationships
    assert r2 in a.relationships
    assert r2 in p2.relationships


def test_get_partner(make_agent, make_relationship):
    a = make_agent()
    p = make_agent()
    rel = make_relationship(a, p)

    assert rel.get_partner(a) == p
    assert rel.get_partner(p) == a


# ============================== AGENT SET TESTS ===============================


def test_AgentSet_init(make_agent):
    s = AgentSet("test")

    assert s.id == "test"
    assert s.members == set()
    assert s.subset == {}

    assert s.parent_set is None
    assert s.numerator == s

    # add another agent set as the child of s
    c = AgentSet("child", s, s)

    assert c.id == "child"
    assert c.parent_set == s
    assert s.subset["child"] == c
    assert c.numerator == s


def test_add_remove_agent(make_agent):
    a = make_agent()
    s = AgentSet("test")
    c = AgentSet("child", s)

    assert s.id == "test"

    c.add_agent(a)
    s.add_agent(a)

    assert s.members == {a}
    assert s.is_member(a)
    assert s.num_members() == 1

    assert c.members == {a}
    assert c.is_member(a)
    assert c.num_members() == 1

    s.remove_agent(a)

    assert s.members == set()
    assert s.is_member(a) is False
    assert s.num_members() == 0

    assert c.members == set()
    assert c.is_member(a) is False
    assert c.num_members() == 0


def test_clear_set(make_agent):
    a = make_agent()
    s = AgentSet("test")
    s.add_agent(a)

    assert s.members == {a}
    assert s.is_member(a)
    assert s.num_members() == 1

    s.clear_set()

    assert s.members == set()
    assert s.is_member(a) == False
    assert s.num_members() == 0
