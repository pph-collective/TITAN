import pytest
import os

from titan.agent import *
from titan.parse_params import create_params

from conftest import FakeRandom

# ============================= AGENT TESTS ============================


@pytest.mark.unit
def test_agent_init(make_agent):
    a = make_agent(init_bond_fields=False)
    b = make_agent(init_bond_fields=False)
    assert b.id == a.id + 1

    # demographics
    assert a.sex_type == "MSM"
    assert a.age == 30
    assert a.race == "black"
    assert a.drug_type == "None"
    assert a.age_bin == 0
    assert a.msmw is False
    assert a.sex_role is "versatile"

    # partner params
    assert a.relationships == set()
    assert a.partners == {}
    assert a.mean_num_partners == {}

    # STI params
    assert a.hiv is False
    assert a.hiv_time == 0
    assert a.aids is False

    # treatment params
    assert a.haart is False
    assert a.haart_time == 0
    assert a.haart_adherence == 0
    assert a.ssp is False
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
    assert a.incar_time == 0


@pytest.mark.unit
def test_get_acute_status(make_agent, params):
    a = make_agent()  # no HIV on init
    assert a.get_acute_status(params.partnership.ongoing_duration) == False
    a.hiv_time = 1  # manually force this to test logic
    assert a.get_acute_status(params.partnership.ongoing_duration) == True


@pytest.mark.unit
def test_iter_partners(make_agent):
    a = make_agent()
    total_partners = 0
    for i in range(3):
        for bond in a.partners:
            total_partners += 1
            a.partners[bond].add(make_agent())

    itered_partners = 0
    for p in a.iter_partners():
        itered_partners += 1

    assert itered_partners == total_partners


@pytest.mark.unit
def test_enroll_prep_choice(make_agent, params):
    params.prep.type = ["Oral", "Inj"]
    params.prep.peak_load = 0.3
    rand_gen = FakeRandom(-0.1)
    a = make_agent()
    a.prep_load = 10

    a.enroll_prep(params, rand_gen)

    assert a.prep
    assert a.intervention_ever
    assert a.prep_last_dose == 0
    assert a.prep_load == 0.3
    assert a.prep_adherence == 1
    assert a.prep_type == "Inj"


@pytest.mark.unit
def test_enroll_prep_one(make_agent, params):
    params.prep.type = ["Oral"]
    params.prep.peak_load = 0.3
    rand_gen = FakeRandom(1.1)
    a = make_agent()
    a.prep_load = 10

    a.enroll_prep(params, rand_gen)

    assert a.prep
    assert a.intervention_ever
    assert a.prep_last_dose == 0
    assert a.prep_load == 0.3
    assert a.prep_adherence == 0
    assert a.prep_type == "Oral"


@pytest.mark.unit
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

    assert a.prep_last_dose == 0
    assert a.prep_load == 0.0


@pytest.mark.unit
def test_get_number_of_sex_acts(make_agent, params):
    a = make_agent()

    rand_gen_low = FakeRandom(0.0)
    min_val_low = params.partnership.sex.frequency[1].min

    rand_gen_high = FakeRandom(1.0)

    assert a.get_number_of_sex_acts(rand_gen_low, params) == min_val_low

    # test fallthrough
    assert a.get_number_of_sex_acts(rand_gen_high, params) == 37


# ============== RELATIONSHIP TESTS ===================


@pytest.mark.unit
def test_relationship(make_agent, make_relationship):
    a = make_agent()
    a.partners["Sex"] = set()
    p1 = make_agent()
    p1.partners["Sex"] = set()
    p2 = make_agent()
    p2.partners["Sex"] = set()
    r1 = make_relationship(a, p1)
    r2 = make_relationship(a, p2)

    assert r1.agent1 == a
    assert r1.agent2 == p1

    # properties
    assert r1.duration == 2
    assert r1.total_sex_acts == 0

    assert r2.duration == 2
    assert r2.total_sex_acts == 0

    assert p1 in a.partners["Sex"]
    assert p2 in a.partners["Sex"]
    assert a in p1.partners["Sex"]
    assert a in p2.partners["Sex"]

    assert r1 in a.relationships
    assert r1 in p1.relationships
    assert r2 in a.relationships
    assert r2 in p2.relationships

    # move forward one time step in the relationship, duration 2 -> 1
    ended = r1.progress()
    assert ended == False
    assert r1.duration == 1
    assert p1 in a.partners["Sex"]
    assert p2 in a.partners["Sex"]
    assert a in p1.partners["Sex"]
    assert a in p2.partners["Sex"]

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
    assert p1 not in a.partners["Sex"]
    assert p2 in a.partners["Sex"]
    assert a not in p1.partners["Sex"]
    assert a in p2.partners["Sex"]

    assert r1 not in a.relationships
    assert r1 not in p1.relationships
    assert r2 in a.relationships
    assert r2 in p2.relationships


@pytest.mark.unit
def test_get_partner(make_agent, make_relationship):
    a = make_agent()
    p = make_agent()
    a.partners["Sex"] = set()
    p.partners["Sex"] = set()
    rel = make_relationship(a, p)

    assert rel.get_partner(a) == p
    assert rel.get_partner(p) == a


# ============================== AGENT SET TESTS ===============================


@pytest.mark.unit
def test_AgentSet_init(make_agent):
    s = AgentSet("test")

    assert s.id == "test"
    assert s.members == set()
    assert s.subset == {}

    assert s.parent_set is None

    # add another agent set as the child of s
    c = AgentSet("child", s)

    assert c.id == "child"
    assert c.parent_set == s
    assert s.subset["child"] == c


@pytest.mark.unit
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


@pytest.mark.unit
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
