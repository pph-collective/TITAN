import pytest

from titan.agent import *
from titan import params

import random


@pytest.fixture
def make_agent():
    def _make_agent(SO="MSM", age=30, race="BLACK", DU="NDU"):
        return Agent(SO, age, race, DU)

    return _make_agent


@pytest.fixture
def make_relationship():
    def _make_relationship(id1, id2, rel_type="#REVIEW", duration=2):
        return Relationship(id1, id2, duration, rel_type)

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
    assert a._mean_num_partners == 0
    assert a._sexualRole == "Vers"

    # STI params
    assert a._HIV_bool is False
    assert a._HIV_time == 0
    assert a._AIDS_bool is False
    assert a._AIDS_time == 0
    assert a._PrEPresistance == 0

    # treatment params
    assert a._HAART_bool is False
    assert a._HAART_time == 0
    assert a._HAART_adh == 0
    assert a._SNE_bool is False
    assert a._treatment_bool is False
    assert a._treatment_time == 0
    assert a._PrEP_reason == []
    assert a.vaccine_time == 0
    assert a.vaccine_type == ""
    assert a.partnerTraced is False
    assert a.awareness is False
    assert a.opinion == 0.0
    assert a.PrEP_type == ""
    assert a._pca is False

    # prevention parameters
    assert a._tested is False
    assert a._PrEP_bool is False
    assert a._PrEP_time == 0
    assert a._PrEP_adh == 0

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


def test_get_acute_status(make_agent):
    a = make_agent()  # no HIV on init
    assert a.get_acute_status() == False
    a._HIV_time = 1  # manually force this to test logic
    assert a.get_acute_status() == True


def test_update_PrEP_load(make_agent):
    a = make_agent()
    assert a._PrEP_lastDose == 0
    assert a._PrEP_load == 0
    a.update_PrEP_load()
    assert a._PrEP_lastDose == 1
    assert a._PrEP_load > 0

    # make time pass
    for i in range(12):
        a.update_PrEP_load()

    assert a._PrEP_lastDose == 13
    assert a._PrEP_load == 0.0


def test_get_transmission_probability(make_agent):
    a = make_agent(race="WHITE", SO="HM")
    a._HAART_adh = 1  # set this explicitly

    p_needle = params.TransmissionProbabilities["NEEDLE"]["1"]
    p_sex = params.TransmissionProbabilities["SEX"]["HM"]["1"]
    scale = params.cal_pXmissionScaling

    # test base case (not tested, not HAART, "WHITE")
    assert a.get_transmission_probability("NEEDLE") == p_needle * scale
    assert a.get_transmission_probability("SEX") == p_sex * scale

    # test acute
    a._HIV_time = 1
    assert (
        a.get_transmission_probability("SEX") == p_sex * scale * params.cal_AcuteScaling
    )
    a._HIV_time = 0

    # test tested status
    a._tested = True
    assert a.get_transmission_probability("SEX") == p_sex * scale * (
        1 - params.cal_RR_Dx
    )
    a._tested = False

    # test HAART
    a._HAART_bool = True
    assert a.get_transmission_probability("SEX") == p_sex * scale * params.cal_RR_HAART
    a._HAART_bool = False

    # test Black
    a._race = "BLACK"
    assert (
        a.get_transmission_probability("SEX") == p_sex * scale * params.cal_raceXmission
    )
    a._race = "WHITE"


def test_get_number_of_sex_acts(make_agent):
    a = make_agent()

    rand_gen_low = FakeRandom(0.0)
    min_val_low = params.sexualFrequency[1]["min"]

    rand_gen_high = FakeRandom(1.0)

    assert a.get_number_of_sexActs(rand_gen_low) == min_val_low

    # test fallthrough
    assert a.get_number_of_sexActs(rand_gen_high) == 37


# ============== RELATIONSHIP TESTS ===================


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


def test_get_partner(make_agent, make_relationship):
    a = make_agent()
    p = make_agent()
    rel = make_relationship(a, p)

    assert rel.get_partner(a) == p
    assert rel.get_partner(p) == a


# ============================== AGENT SET TESTS ===============================


def test_Agent_set_init(make_agent):
    s = Agent_set("test")

    assert s._ID == "test"
    assert s._members == []
    assert s._subset == {}

    assert s._parent_set is None
    assert s._numerator == s

    # add another agent set as the child of s
    c = Agent_set("child", s, s)

    assert c._ID == "child"
    assert c._parent_set == s
    assert s._subset["child"] == c
    assert c._numerator == s


def test_add_remove_agent(make_agent):
    a = make_agent()
    s = Agent_set("test")
    c = Agent_set("child", s)

    assert s.get_ID() == "test"

    c.add_agent(a)
    s.add_agent(a)

    assert s._members == [a]
    assert s.is_member(a)
    assert s.num_members() == 1

    assert c._members == [a]
    assert c.is_member(a)
    assert c.num_members() == 1

    s.remove_agent(a)

    assert s._members == []
    assert s.is_member(a) is False
    assert s.num_members() == 0

    assert c._members == []
    assert c.is_member(a) is False
    assert c.num_members() == 0


def test_clear_set(make_agent):
    a = make_agent()
    s = Agent_set("test")
    s.add_agent(a)

    assert s._members == [a]
    assert s.is_member(a)
    assert s.num_members() == 1

    s.clear_set()

    assert s._members == []
    assert s.is_member(a) == False
    assert s.num_members() == 0
