import pytest

from titan.agent import *


@pytest.fixture
def make_agent():
    def _make_agent(SO="MSM", age=30, race="BLACK", DU="NDU"):
        return Agent(SO, age, race, DU)

    return _make_agent


@pytest.fixture
def make_relationship():
    def _make_relationship(id1, id2, SO="MSM", rel_type="#REVIEW", duration=2):
        return Relationship(id1, id2, SO, rel_type, duration)

    return _make_relationship


def test_agent_init(make_agent):
    a = make_agent()
    b = make_agent()
    assert b._ID == a._ID + 1
    assert a._timeAlive == 0

    # demographics
    assert a._SO == "MSM"
    assert a._age == 30
    assert a._race == "BLACK"
    assert a._DU == "NDU"
    assert a._ageBin == 0
    assert a._MSMW is False

    # partner params
    assert a._relationships == []
    assert a._partners == []
    assert a._num_sex_partners == 0
    assert a._num_NE_partners == 0
    assert a._mean_num_partners == 0
    assert a._sexualRole == "Vers"

    # STI params
    assert a._HIV_bool is False
    assert a._HIV_time == 0
    assert a._AIDS_bool is False
    assert a._AIDS_time == 0
    assert a._PrEPresistance == 0

    # treatment params
    assert a._tested is False
    assert a._HAART_bool is False
    assert a._HAART_time == 0
    assert a._HAART_adh == 0
    assert a._SNE_bool is False
    assert a._PrEP_bool is False
    assert a._PrEP_time == 0
    assert a._PrEP_adh == 0
    assert a._treatment_bool is False
    assert a._treatment_time == 0
    assert a._PrEP_reason == []

    # prep pharmacokinetics
    assert a._PrEP_load == 0.0
    assert a._PrEP_lastDose == 0

    # high risk params
    assert a._highrisk_bool is False
    assert a._highrisk_time == 0
    assert a._everhighrisk_bool is False

    # incarceration
    assert a._incar_bool is False
    assert a._ever_incar_bool is False
    assert a._incar_time == 0
    assert a._incar_treatment_time == 0


def test_get_id(make_agent):
    a = make_agent()
    assert a.get_ID() == a._ID


def test_partner_list(make_agent, make_relationship):
    a = make_agent()

    assert a.partner_list() == []

    p = make_agent()
    r = make_relationship(a, p)

    assert a.partner_list() == [p._ID]
    assert p.partner_list() == [a._ID]


def test_relationship(make_agent, make_relationship):
    a = make_agent()
    p1 = make_agent()
    p2 = make_agent()
    r1 = make_relationship(a, p1)
    r2 = make_relationship(a, p2)

    assert r1._ID1 == a
    assert r1._ID2 == p1
    assert r2.get_ID() == r1.get_ID() + 1

    # properties
    assert r1._duration == 2
    assert r1._total_sex_acts == 0

    assert r2._duration == 2
    assert r2._total_sex_acts == 0

    assert a._num_sex_partners == 2
    assert p1._num_sex_partners == 1
    assert p2._num_sex_partners == 1

    assert p1._ID in a.partner_list()
    assert p2._ID in a.partner_list()
    assert a._ID in p1.partner_list()
    assert a._ID in p2.partner_list()

    assert r1 in a._relationships
    assert r1 in p1._relationships
    assert r2 in a._relationships
    assert r2 in p2._relationships

    # move forward one time step in the relationship, duration 2 -> 1
    ended = r1.progress()
    assert ended == False
    assert r1._duration == 1
    assert p1._ID in a.partner_list()
    assert p2._ID in a.partner_list()
    assert a._ID in p1.partner_list()
    assert a._ID in p2.partner_list()

    assert r1 in a._relationships
    assert r1 in p1._relationships
    assert r2 in a._relationships
    assert r2 in p2._relationships

    # move forward one more timestep, duration 1 -> 0, rel over on next progress
    ended = r1.progress()
    assert r1._duration == 0
    ended = r1.progress()
    assert ended == True
    assert r1._duration == 0
    assert p1._ID not in a.partner_list()
    assert p2._ID in a.partner_list()
    assert a._ID not in p1.partner_list()
    assert a._ID in p2.partner_list()

    assert r1 not in a._relationships
    assert r1 not in p1._relationships
    assert r2 in a._relationships
    assert r2 in p2._relationships


def test_Agent_set_init(make_agent):
    s = Agent_set("test")

    assert s._ID == "test"
    assert s._members == []
    assert s._subset == {}

    assert s._parent_set is None
    assert s._numerator == s


def test_add_remove_agent(make_agent):
    a = make_agent()  # 9
    s = Agent_set("test")

    s.add_agent(a)

    assert s._members == [a]
    assert s.is_member(a)
    assert s.num_members() == 1

    s.remove_agent(a)

    assert s._members == []
    assert s.is_member(a) is False
    assert s.num_members() == 0
