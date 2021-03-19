import pytest
import os

from titan.partnering import *
from titan.agent import Agent, Relationship
from titan.population import Population
from titan.parse_params import create_params, ObjMap

from conftest import FakeRandom


def test_partnership_duration(params):
    # test duration with randint
    assert get_partnership_duration(params, FakeRandom(1.0), "Inj", "white") == 3
    # test duration with bins
    assert get_partnership_duration(params, FakeRandom(0.1), "Sex", "white") == 1
    # test duration with second race
    assert get_partnership_duration(params, FakeRandom(0.1), "Sex", "black") == 3


@pytest.mark.unit
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


@pytest.mark.unit
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
                params.demographics[race].sex_type.MSM.drug_type.Inj.num_partners[
                    bond
                ].vars[1].value = 0.0

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


@pytest.mark.unit
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


@pytest.mark.unit
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


@pytest.mark.unit
def test_get_assort_partner_race(make_population, make_agent, params):
    pop = make_population()
    a = make_agent(SO="MSM", race="white")
    p1 = make_agent(SO="MSM", race="white")
    p2 = make_agent(SO="MSM", race="black")
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

    # assort with white
    test_rule = ObjMap(
        {
            "attribute": "race",
            "partner_attribute": "__agent__",
            "bond_types": [],
            "agent_value": "white",
            "partner_values": {"white": 0.9, "__other__": 0.1},
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
            "partner_attribute": "__agent__",
            "bond_types": [],
            "agent_value": "white",
            "partner_values": {"white": 0.9, "__other__": 10},
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


@pytest.mark.unit
def test_get_assort_partner_high_risk(make_population, make_agent, params):
    pop = make_population()
    a = make_agent(SO="MSM", race="white")
    p1 = make_agent(SO="MSM", race="white")
    p2 = make_agent(SO="MSM", race="black")

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
            "partner_attribute": "__agent__",
            "bond_types": [],
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
            "partner_attribute": "__agent__",
            "bond_types": [],
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


@pytest.mark.unit
def test_get_assort_partner_drug_type(make_population, make_agent, params):
    pop = make_population()
    a = make_agent(SO="MSM", race="white", DU="Inj")
    p1 = make_agent(SO="MSM", race="white", DU="Inj")
    p2 = make_agent(SO="MSM", race="black", DU="None")
    p3 = make_agent(SO="MSM", race="black", DU="NonInj")

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
                params.demographics[race].sex_type.MSM.drug_type.Inj.num_partners[
                    bond
                ].vars[1].value = 0.0

    # assort with Inj
    test_rule = ObjMap(
        {
            "attribute": "drug_type",
            "partner_attribute": "__agent__",
            "bond_types": [],
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
            "attribute": "drug_type",
            "partner_attribute": "__agent__",
            "bond_types": [],
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


@pytest.mark.unit
def test_get_assort_same_race(make_population, make_agent, params):
    pop = make_population()
    a = make_agent(SO="MSM", race="white")
    p1 = make_agent(SO="MSM", race="white")
    p2 = make_agent(SO="MSM", race="black")
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

    # assort with white
    test_rule = ObjMap(
        {
            "attribute": "race",
            "partner_attribute": "__agent__",
            "bond_types": [],
            "agent_value": "__any__",
            "partner_values": {"__same__": 0.9, "__other__": 0.1},
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
    params.assort_mix["test_rule"]["partner_values"]["__other__"] = 10

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

    # throw error on key that's not __same__ or __other__
    params.assort_mix["test_rule"]["partner_values"]["a_key"] = 100

    with pytest.raises(ValueError):
        select_partner(
            a,
            pop.all_agents.members,
            pop.sex_partners,
            pop.pwid_agents,
            params,
            FakeRandom(0.5),
            "Sex",
        )


@pytest.mark.unit
def test_get_assort_hiv(make_population, make_agent, params):
    pop = make_population()
    a = make_agent(SO="MSM", race="white")
    p1 = make_agent(SO="MSM", race="white")
    p2 = make_agent(SO="MSM", race="black")
    a.hiv.active = True
    p1.hiv.active = True
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

    # assort with white
    test_rule = ObjMap(
        {
            "attribute": "hiv.active",
            "partner_attribute": "__agent__",
            "bond_types": [],
            "agent_value": True,
            "partner_values": {"True": 0.9, "False": 0.1},
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
            "attribute": "hiv.active",
            "partner_attribute": "__agent__",
            "bond_types": [],
            "agent_value": True,
            "partner_values": {"True": 0.9, "False": 10},
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


@pytest.mark.unit
def test_get_assort_hiv_prep(make_population, make_agent, params):
    pop = make_population()
    a = make_agent(SO="MSM", race="white")
    p1 = make_agent(SO="MSM", race="white")
    p2 = make_agent(SO="MSM", race="black")
    a.hiv.active = True
    p1.prep.active = True
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

    # assort with white
    test_rule = ObjMap(
        {
            "attribute": "hiv.active",
            "partner_attribute": "prep.active",
            "bond_types": [],
            "agent_value": True,
            "partner_values": {"True": 0.9, "False": 0.1},
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
            "attribute": "hiv.active",
            "partner_attribute": "prep.active",
            "bond_types": [],
            "agent_value": True,
            "partner_values": {"True": 0.9, "False": 10},
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


@pytest.mark.unit
def test_get_assort_bond_type(make_population, make_agent, params):
    pop = make_population()
    a = make_agent(SO="MSM", race="white")
    p = make_agent(SO="MSM", race="white")
    a.hiv.active = True
    p.hiv.active = True
    for bond in params.classes.bond_types:
        a.target_partners[bond] = 1
        p.target_partners[bond] = 1
        a.partners[bond] = set()
        p.partners[bond] = set()
    pop.add_agent(a)
    pop.add_agent(p)

    params.features.assort_mix = True

    # assort with white
    test_rule = ObjMap(
        {
            "attribute": "hiv.active",
            "partner_attribute": "__agent__",
            "bond_types": ["Social"],
            "agent_value": True,
            "partner_values": {"True": 0.1, "False": 0.9},
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
        "Social",
    )

    assert partner is None

    # get __other__
    test_rule = ObjMap(
        {
            "attribute": "hiv.active",
            "partner_attribute": "__agent__",
            "bond_types": ["Sex"],
            "agent_value": True,
            "partner_values": {"True": 0.9, "False": 0.1},
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

    assert partner == p


@pytest.mark.unit
def test_get_str_attr(make_agent):
    a = make_agent(race="white", age=42)
    assert get_str_attr(a, "race") == "white"
    assert get_str_attr(a, "age") == "42"
    assert get_str_attr(a, "location.name") == "world"
    assert get_str_attr(a, "hiv.active") == "False"


@pytest.mark.unit
def test_get_partner_attr():
    test_rule = ObjMap(
        {
            "attribute": "drug_type",
            "partner_attribute": "__agent__",
            "agent_value": "Inj",
            "partner_values": {"Inj": 0.8, "NonInj": 0.1, "__other__": 10},
        }
    )
    assert get_partner_attr(test_rule) == "drug_type"

    test_rule.partner_attribute = "race"
    assert get_partner_attr(test_rule) == "race"


@pytest.mark.unit
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
