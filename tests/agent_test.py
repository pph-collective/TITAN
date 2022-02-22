import pytest
import os

from titan.agent import *

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
    assert a.sex_role is "versatile"

    # partner params
    assert a.relationships == set()
    assert a.partners == {}
    assert a.mean_num_partners == {}

    # STI params
    assert a.hiv.active is False
    assert a.hiv.dx is False
    assert a.hiv.time is None
    assert a.hiv.aids is False


@pytest.mark.unit
def test_get_partners(make_agent):
    a = make_agent()
    p1 = make_agent()
    p2 = make_agent()

    assert "Inj" in a.partners.keys()
    assert "Sex" in a.partners.keys()
    a.partners["Sex"].update({p1})
    a.partners["Inj"].update({p2})

    assert a.get_partners() == {p1, p2}
    assert a.get_partners(["Sex"]) == {p1}
    assert a.get_partners(["Inj"]) == {p2}


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
def test_is_msm(make_agent):
    a = make_agent()
    assert a.is_msm()
    a.sex_type = "MTF"
    assert not a.is_msm()


@pytest.mark.unit
def test_has_partners(make_agent, make_relationship):
    a = make_agent()

    assert a.has_partners() is False

    p = make_agent()
    r = make_relationship(a, p)

    assert a.has_partners() is True


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

    assert r2.duration == 2

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
    a2 = make_agent()
    a.partners["Sex"] = set()
    p.partners["Sex"] = set()
    rel = make_relationship(a, p)

    assert rel.get_partner(a) == p
    assert rel.get_partner(p) == a
    with pytest.raises(ValueError):
        rel.get_partner(a2)


@pytest.mark.unit
def test_get_number_of_sex_acts(make_agent, make_relationship, params):
    a = make_agent()
    p = make_agent()
    rel = make_relationship(a, p)

    rand_gen_low = FakeRandom(0.0)
    min_val_low = params.partnership.sex.frequency.Sex.bins[1].min
    max_val_high = params.partnership.sex.frequency.Sex.bins[2].max

    rand_gen_high = FakeRandom(1.0)

    assert rel.get_number_of_sex_acts(rand_gen_low) == min_val_low

    # test fallthrough
    assert rel.get_number_of_sex_acts(rand_gen_high) == max_val_high

    # test with distribution; should be independent of random
    a.location.params.partnership.sex.frequency.Sex.type = "distribution"

    assert rel.get_number_of_sex_acts(rand_gen_low) == 0
    assert rel.get_number_of_sex_acts(rand_gen_high) == 0

    a.location.params.partnership.sex.frequency.Sex.type = "not a thing"
    with pytest.raises(Exception):
        rel.get_number_of_sex_acts(rand_gen_low)


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
    assert a in s
    assert s.num_members() == 1

    assert c.members == {a}
    assert a in c
    assert c.num_members() == 1

    s.remove_agent(a)

    assert s.members == set()
    assert a not in s
    assert s.num_members() == 0

    assert c.members == set()
    assert a not in c
    assert c.num_members() == 0


@pytest.mark.unit
def test_clear_set(make_agent):
    a = make_agent()
    s = AgentSet("test")
    s.add_agent(a)

    assert s.members == {a}
    assert a in s
    assert s.num_members() == 1

    s.clear_set()

    assert s.members == set()
    assert a not in s
    assert s.num_members() == 0


@pytest.mark.unit
def test_agebin(make_agent):
    a = make_agent()
    a.age = 30
    assert a.agebin == 1
